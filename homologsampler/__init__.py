import os, gzip, warnings, sys
from collections import Counter

import click

from cogent import LoadSeqs, DNA
from cogent.db.ensembl import Compara, Genome, HostAccount, Species

from scitrack import CachingLogger

from homologsampler.util import (species_names_from_csv, missing_species_names,
        get_chrom_names, load_coord_names, display_available_dbs)

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2014, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "0.11"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

LOGGER = CachingLogger(create_dir=True)

def get_one2one_orthologs(compara, ref_genes, outpath, not_strict, force_overwrite, test):
    """writes one-to-one orthologs of protein coding genes to outpath"""
    
    species = Counter(compara.Species)
    written = 0
    with click.progressbar(ref_genes,
        label="Finding 1to1 orthologs") as ids:
        for gene in ids:
            outfile_name = os.path.join(outpath, "%s.fa.gz" % gene)
            if os.path.exists(outfile_name) and not force_overwrite:
                written += 1
                continue
            
            syntenic = compara.getRelatedGenes(StableId=gene,
                            Relationship='ortholog_one2one')
        
            if not not_strict and (syntenic is None or Counter(syntenic.getSpeciesSet()) != species):
                # skipping, not all species had a 1to1 ortholog for this gene
                continue
            
            seqs = []
            for m in syntenic.Members:
                name = Species.getCommonName(m.genome.Species)
                cds = m.CanonicalTranscript.Cds.withoutTerminalStopCodon(allow_partial=True)
                cds.Name = name
                seqs.append([name, cds])
            
            seqs = LoadSeqs(data=seqs, aligned=False)
            if test:
                print()
                print(gene)
                print(seqs.toFasta())
            else:
                with gzip.open(outfile_name, 'w') as outfile:
                    outfile.write(seqs.toFasta() + '\n')
                LOGGER.output_file(outfile_name)
            
            written += 1
    if test:
        msg = "Would have written %d files to %s" % (written, outpath)
    else:
        msg = "Wrote %d files to %s" % (written, outpath)
    
    click.echo(msg)
    
    return

def get_latin_from_label(label):
    """returns latin name from the sequence label"""
    return label.split(':')[0]

def renamed_seqs(aln):
    """renames sequences to be just species common name"""
    new = []
    names = Counter()
    for seq in aln.Seqs:
        latin = get_latin_from_label(seq.Name)
        common = Species.getCommonName(latin)
        names[common] += 1
        seq.Name = common
        new.append((seq.Name, seq))
    
    if max(names.values()) > 1:
        return None
    
    return LoadSeqs(data=new, moltype=DNA)

def with_masked_features(aln, reverse=False):
    """returns an alignment with the tandem repeat sequences masked"""
    for name in aln.Names:
        seq = aln.getSeq(name)
        remove = []
        for a in seq.annotations:
            if a.Name not in ('trf', 'cpg', 'splice'):
                remove.append(a)
                break
        
        if remove:
            seq.detachAnnotations(remove)
    
    if reverse:
        aln = aln.rc()
    
    aln = aln.withMaskedAnnotations(['repeat', 'cpg'], mask_char='?')
    aln = aln.withMaskedAnnotations(['exon'], mask_char='?')
    return aln

def get_syntenic_alignments_introns(compara, ref_genes, outpath, method_clade_id,
                mask_features, outdir, force_overwrite, test):
    """writes Ensembl `method` syntenic alignments to ref_genes"""
    species = Counter(compara.Species)
    common_names = list(map(Species.getCommonName, compara.Species))
    filler = LoadSeqs(data=[(n, 'N') for n in common_names], moltype=DNA)
    
    written = 0
    with click.progressbar(ref_genes,
        label="Finding 1to1 intron orthologs") as ids:
        for gene_id in ids:
            gene = _get_gene_from_compara(compara, gene_id)
            if not gene:
                continue
            
            if gene.CanonicalTranscript.Introns is None:
                continue
            
            outfile_name = os.path.join(outpath, "%s.fa.gz" % gene.StableId)
            if os.path.exists(outfile_name) and not force_overwrite:
                written += 1
                continue
            
            syntenic_regions = []
            regions = list(compara.getSyntenicRegions(region=gene.CanonicalTranscript,
                                                method_clade_id=str(method_clade_id)))
            alignments = []
            for index, region in enumerate(regions):
                if region is None:
                    continue
                
                try:
                    got = Counter(region.getSpeciesSet())
                except (AttributeError, AssertionError):
                    # this is a PyCogent bug
                    error = sys.exc_info()
                    err_type = str(error[0]).split('.')[-1][:-2]
                    err_msg = str(error[1])
                    msg = 'gene_stable_id=%s; err_type=%s; msg=%s' % (gene.StableId, err_type, err_msg)
                    click.echo(click.style("ERROR:" + msg, fg="red"))
                    LOGGER.log_message(msg, label="ERROR")
                    continue
                
                if got != species:
                    continue
                
                if mask_features:
                    aln = region.getAlignment(feature_types=['gene', 'repeat', 'cpg'])
                    aln = with_masked_features(aln, reverse=gene.Location.Strand == -1)
                else:
                    aln = region.getAlignment()
                
                if aln is None:
                    continue
                
                aln = renamed_seqs(aln)
                if aln is not None:
                    alignments.append(aln)
            
            if not alignments:
                continue
            
            # we put a column of Ns between syntenic regions so that subsequent
            # sampling for tuple aligned columns does not construct artificial
            # motifs
            align = None
            for aln in alignments:
                if align is None:
                    align = aln
                    continue
                
                align += (filler + aln)
            
            if test:
                print(align)
            else:
                with gzip.open(outfile_name, 'w') as outfile:
                    outfile.write(align.toFasta())
                LOGGER.output_file(outfile_name)
            
            written += 1
        
    print("Wrote %d files to %s" % (written, outpath))
    return


def display_ensembl_alignment_table(compara):
    """prints the method_species_link table and then exits"""
    compara.method_species_links.Legend = \
    "Assign the desired value from method_link_species_set_id to the"\
    " method_clade_id argument"
    print(compara.method_species_links)
    exit(0)

def _get_account(ensembl_account):
    """returns HostAccount or None"""
    try:
        acc = HostAccount(*ensembl_account.split())
    except (KeyError, AttributeError):
        warnings.warn("ENSEMBL_ACCOUNT environment variable not set, defaulting to UK sever. Slow!!")
        acc = None
    
    return acc

def _get_gene_from_compara(compara, stable_id):
    """returns gene instance from a compara db"""
    for sp in list(compara._genomes.values()):
        gene = sp.getGeneByStableId(stable_id)
        if gene:
            break
    return gene

def _get_ref_genes(ref_genome, chroms, limit, biotype='protein_coding'):
    """returns stable ID's for genes from reference genome"""
    print("Sampling %s genes" % ref_genome)
    all_genes = ref_genome.getGenesMatching(BioType=biotype)
    
    ref_genes = []
    with click.progressbar(all_genes,
        label="Finding genes") as genes:
        for index, g in enumerate(genes):
            if limit is not None and index >= limit:
                break
            if chroms and g.Location.CoordName not in chroms:
                continue
            
            ref_genes.append(g.StableId)
    return ref_genes

class Config(object):
    def __init__(self):
        super(Config, self).__init__()
        self.force_overwrite = False
        self.test = False

pass_config = click.make_pass_decorator(Config, ensure=True)

@click.group()
@click.option('--ensembl_account', envvar='ENSEMBL_ACCOUNT',
    help="shell variable with MySQL account details, e.g. export ENSEMBL_ACCOUNT='myhost.com jill jills_pass'")
@click.option('-F', '--force_overwrite', is_flag=True, help="Overwrite existing files.")
@click.option('--test', is_flag=True,
    help="sets limit # queries to 2, does not write files, prints seqs and exits.")
@click.version_option(version=__version__)
@pass_config
def cli(ctx, ensembl_account, force_overwrite, test):
    ctx.ensembl_account = _get_account(ensembl_account)
    ctx.force_overwrite = force_overwrite
    ctx.test = test

@cli.command()
@click.option('--release', help='Ensembl release.')
@pass_config
def show_available_species(ctx, release):
    """shows available species and Ensembl release at ENSEMBL_ACCOUNT"""
    available = display_available_dbs(ctx.ensembl_account, release)
    available.Title = "Species available at: %s" % str(ctx.ensembl_account)
    print(available)
    sys.exit(0)
    
@cli.command()
@click.option('--species', required=True, help='Comma separated list of species names.')
@click.option('--release', required=True, help='Ensembl release.')
@pass_config
def show_align_methods(ctx, species, release):
    """Shows the align methods in release ENSEMBL_ACCOUNT and exits."""
    species = species_names_from_csv(species)
    missing_species = missing_species_names(species)
    if missing_species:
        msg = ["The following species names don't match an Ensembl record. Check spelling!",
              str(missing_species),
              "\nAvailable species are at this server are:",
              str(display_available_dbs(account))]
        
        click.echo(click.style("\n".join(msg), fg="red"))
        sys.exit(-1)
    
    compara = Compara(species, Release=release, account=ctx.ensembl_account)
    
    if show_align_methods:
        display_ensembl_alignment_table(compara)
    
    

@cli.command()
@click.option('--species', required=True, help='Comma separated list of species names.')
@click.option('--release', required=True, help='Ensembl release.')
@click.option('--outdir', required=True, type=click.Path(resolve_path=True), help='Path to write files.') ##
@click.option('--ref', default=None, help='Reference species.')
@click.option('--ref_genes_file', default=None, type=click.File('rb'),
    help='File containing Ensembl stable identifiers for genes of interest. '
         'One identifier per line.')
@click.option('--coord_names', default=None, type=click.Path(resolve_path=True),
                help='File containing chrom/coord names to restrict sampling to, one per line.')
@click.option('--not_strict', is_flag=True,
    help="Genes with an ortholog in any species are exported. "
         "Default is all species must have a ortholog.")
@click.option('--introns', is_flag=True, help="Sample syntenic alignments of introns, requires --method_clade_id.")
@click.option('--method_clade_id',
    help="The value of method_link_species_set_id to use (see ) "
    "required if sampling introns.")
@click.option('--mask_features', is_flag=True, help="Intron masks repeats, exons, CpG islands.")
@click.option('--limit', type=int, default=0, help="Limit to this number of genes.")
@click.option('--logfile_name', default="one2one.log", help="Name for log file, written to outdir.")
@pass_config
def one2one(ctx, species, release, outdir, ref, ref_genes_file, coord_names, not_strict, introns, method_clade_id, mask_features, logfile_name, limit):
    """Command line tool for sampling homologous sequences from Ensembl."""
    if not any([ref, ref_genes_file]):
        # just the command name, indicate they need to display help
        msg = "%s\n\n--help to see all options\n" % ctx.get_usage()
        click.echo(msg)
        exit(-1)
        
    
    args = locals()
    args.pop("ctx")
    args.update(vars(ctx))
    args['ensembl_account'] = str(args['ensembl_account'])
    LOGGER.log_message(str(args), label="params")
    
    if ctx.test:
        limit = 2
    else:
        limit = limit or None
    
    if (introns and not method_clade_id) or (mask_features and not introns):
        msg = ["Must specify the introns and method_clade_id in order to export introns.",
               "Use show_align_methods to see the options"]
        click.echo(click.style("\n".join(msg), fg="red"))
        exit(-1)
    
    species = species_names_from_csv(species)
    species_missing = missing_species_names(species)
    if species_missing:
        msg = ["The following species names don't match an Ensembl record. Check spelling!",
              str(species_missing),
              "\nAvailable species are at this server are:",
              str(display_available_dbs(ctx.ensembl_account))]
        
        click.echo(click.style("\n".join(msg), fg="red"))
        exit(-1)
    
    if ref and ref not in species:
        print("The reference species not in species names")
        exit(-1)
    
    compara = Compara(species, Release=release, account=ctx.ensembl_account)
    runlog_path = os.path.join(outdir, logfile_name)
    
    if os.path.exists(runlog_path) and not ctx.force_overwrite:
        msg = ["Log file (%s) already exists!" % runlog_path,
               "Use force_overwrite or provide logfile_name"]
        click.echo(click.style("\n".join(msg), fg="red"))
        exit(-1)
    
    if not ctx.test:
        LOGGER.log_file_path = runlog_path
    
    chroms = None
    if coord_names:
        chroms = load_coord_names(coord_names)
        LOGGER.input_file(coord_names)
    elif coord_names and ref:
        chroms = get_chrom_names(ref, compara)
    
    if not os.path.exists(outdir) and not ctx.test:
        os.makedirs(outdir)
        print("Created", outdir)
    
    if ref and not ref_genes_file:
        ref_genome = Genome(ref, Release=release, account=ctx.ensembl_account)
        ref_genes = _get_ref_genes(ref_genome, chroms, limit)
    else:
        ref_genes = [l.strip() for l in ref_genes_file if l.strip()]
    
    if not introns:
        print("Getting orthologs")
        get_one2one_orthologs(compara, ref_genes, outdir, not_strict, ctx.force_overwrite, ctx.test)
    else:
        print("Getting orthologous introns")
        get_syntenic_alignments_introns(compara, ref_genes, outdir, method_clade_id,
                mask_features, outdir, ctx.force_overwrite, ctx.test)
    

if __name__ == "__main__":
    cli()
