import os, gzip, warnings, sys
from collections import Counter

import click

from cogent import LoadSeqs, DNA
from cogent.db.ensembl import Compara, Genome, HostAccount, Species

from scitrack import CachingLogger

from homologsampler.util import (missing_species_names, get_chrom_names,
                        load_coord_names, display_available_dbs)

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2014, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

LOGGER = CachingLogger(create_dir=True)

def get_one2one_orthologs(compara, ref_genes, outpath, force_overwrite):
    """writes one-to-one orthologs of protein coding genes to outpath"""
    
    species = Counter(compara.Species)
    written = 0
    with click.progressbar(ref_genes,
        label="Finding 1to1 orthologs") as ids:
        for gene in ids:
            outfile_name = os.path.join(outpath, "%s.fa.gz" % gene.StableId)
            if os.path.exists(outfile_name) and not force_overwrite:
                written += 1
                continue
            
            syntenic = compara.getRelatedGenes(StableId=gene.StableId,
                            Relationship='ortholog_one2one')
        
            if syntenic is None or Counter(syntenic.getSpeciesSet()) != species:
                # skipping, not all species had a 1to1 ortholog for this gene
                continue
            
            seqs = []
            for m in syntenic.Members:
                name = Species.getCommonName(m.genome.Species)
                cds = m.CanonicalTranscript.Cds.withoutTerminalStopCodon(allow_partial=True)
                cds.Name = name
                seqs.append([name, cds])
            seqs = LoadSeqs(data=seqs, aligned=False)
            with gzip.open(outfile_name, 'w') as outfile:
                outfile.write(seqs.toFasta() + '\n')
            
            written += 1
            LOGGER.output_file(outfile_name)
        
    print "Wrote %d files to %s" % (written, outpath)
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
                mask_features, outdir, force_overwrite):
    """writes Ensembl `method` syntenic alignments to ref_genes"""
    species = Counter(compara.Species)
    common_names = map(Species.getCommonName, compara.Species)
    filler = LoadSeqs(data=[(n, 'N') for n in common_names], moltype=DNA)
    
    written = 0
    with click.progressbar(ref_genes,
        label="Finding 1to1 intron orthologs") as ids:
        for gene in ids:
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
            
            with gzip.open(outfile_name, 'w') as outfile:
                outfile.write(align.toFasta())
            
            written += 1
            LOGGER.output_file(outfile_name)
        
    print "Wrote %d files to %s" % (written, outpath)
    return


def display_ensembl_alignment_table(compara):
    """prints the method_species_link table and then exits"""
    compara.method_species_links.Legend = \
    "Assign the desired value from method_link_species_set_id to the"\
    " method_clade_id argument"
    print compara.method_species_links
    exit(0)

@click.command()
@click.option('--ref', help='Reference species.')
@click.option('--species', help='Comma separated list of species names.')
@click.option('--release', help='Ensembl release.')
@click.option('--outdir', type=click.Path(resolve_path=True), help='Path to write files.') ##
@click.option('--ensembl_account', envvar='ENSEMBL_ACCOUNT',
    help="shell variable with MySQL account details, e.g. export ENSEMBL_ACCOUNT='myhost.com jill jills_pass'")
@click.option('--coord_names', default=None, type=click.Path(resolve_path=True),
                help='File containing chrom/coord names to restrict sampling to, one per line.')
@click.option('--introns', is_flag=True, help="Sample syntenic alignments of introns, requires --method_clade_id.")
@click.option('--method_clade_id',
    help="The align method ID to use, required if sampling introns.")
@click.option('--mask_features', is_flag=True, help="Intron masks repeats, exons, CpG islands.")
@click.option('--force_overwrite', is_flag=True, help="Overwrite existing files.")
@click.option('--show_align_methods', is_flag=True, help="Shows the align methods and exits.")
@click.option('--show_available_species', is_flag=True, help="Shows the available db's at ENSEMBL_ACCOUNT.")
@click.option('--limit', type=int, default=0, help="Limit to this number of genes.")
@click.option('--logfile_name', default="one2one.log", help="Name for log file, written to outdir.")
@click.option('--test', is_flag=True)
@click.pass_context
def main(ctx, ref, species, release, outdir, ensembl_account, coord_names, introns, method_clade_id, mask_features, force_overwrite, show_align_methods, logfile_name, limit, show_available_species, test):
    """Command line tool for sampling homologous sequences from Ensembl."""
    # There are XX possible uses
    # 1 - just the command name, in which case indicate they need to display help
    # 2 - list available databases at a server
    #   - they need the --show_available_species arg only
    # 3 - list --show_align_methods
    #   - they also need the species and release
    # 4 - query for orthologs
    #   - introns
    #       - need method_clade_id, ref species, species, release, outdir
    #   - else
    #       - need ref species, species, release, outdir
    if not any([show_align_methods, show_available_species, ref, species]):
        # just the command name, indicate they need to display help
        msg = "%s\n\n--help to see all options\n" % ctx.get_usage()
        click.echo(msg)
        exit(-1)
    elif show_align_methods and not all([species, release]):
        msg = ["The following arguments are required for show_align_methods:",
               "--species", "--release"]
        click.echo(click.style("\n".join(msg), fg="red"))
        exit(-1)
        
    try:
        acc = HostAccount(*ensembl_account.split())
    except (KeyError, AttributeError):
        warnings.warn("ENSEMBL_ACCOUNT environment variable not set, defaulting to UK sever. Slow!!")
        acc = None
        pass
    
    if show_available_species:
        available = display_available_dbs(acc)
        print
        print available
        exit(0)
        
    
    args = locals()
    args.pop("ctx")
    args.update(ctx.params)
    LOGGER.log_message(str(args), label="params")
    
    limit = limit or None
    if not show_align_methods and not all([species, release, outdir]):
        msg = ["You are missing an argument. Either",
               "--outdir (for exporting orthologs)",
               "OR",
               "--show_align_methods (if you intend exporting introns)"]
        
        msg.append("")
        msg.append("Use --help to see options")
        
        click.echo(click.style("\n".join(msg), fg="red"))
        exit(-1)
    
    if introns and not method_clade_id:
        msg = ["Must specify the method_clade_id in order to export introns.",
               "Use the --show_align_methods argument to see the options"]
        click.echo(click.style("\n".join(msg), fg="red"))
        exit(-1)
    
    species = species.split(',')
    species_missing = missing_species_names(species)
    if species_missing:
        msg = ["The following species names don't match an Ensembl record. Check spelling!",
              str(species_missing),
              "\nAvailable species are at this server are:",
              str(display_available_dbs(acc))]
        
        click.echo(click.style("\n".join(msg), fg="red"))
        exit(-1)
    
    compara = Compara(species, Release=release, account=acc)
    
    if show_align_methods:
        display_ensembl_alignment_table(compara)
    
    if not ref in species:
        print "The reference species not in species names"
        exit(-1)
    
    ref_genome = Genome(ref, Release=release, account=acc)
    runlog_path = os.path.join(outdir, logfile_name)
    
    if os.path.exists(runlog_path) and not force_overwrite:
        msg = ["Log file (%s) already exists!" % runlog_path,
               "Use force_overwrite or provide logfile_name"]
        click.echo(click.style("\n".join(msg), fg="red"))
        exit(-1)
    
    LOGGER.log_file_path = runlog_path
    
    if coord_names:
        chroms = load_coord_names(coord_names)
        LOGGER.input_file(coord_names)
    else:
        chroms = get_chrom_names(ref, compara)
    
    chroms = chroms or None
    
    if not os.path.exists(outdir) and not test:
        os.makedirs(outdir)
        print "Created", outdir
    
    print "Sampling %s genes" % ref_genome
    all_genes = ref_genome.getGenesMatching(BioType='protein_coding')
    
    ref_genes = []
    with click.progressbar(all_genes,
        label="Finding %s genes" % ref) as genes:
        for index, g in enumerate(genes):
            if limit is not None and index >= limit:
                break
            if chroms and g.Location.CoordName not in chroms:
                continue
            
            ref_genes.append(g)
    
    if not introns:
        print "Getting orthologs"
        get_one2one_orthologs(compara, ref_genes, outdir, force_overwrite)
    else:
        print "Getting orthologous introns"
        get_syntenic_alignments_introns(compara, ref_genes, outdir, method_clade_id,
                mask_features, outdir, force_overwrite)
    

if __name__ == "__main__":
    main()
