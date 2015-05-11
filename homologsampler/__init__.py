import os, gzip, warnings
from collections import Counter

import click

from cogent import LoadSeqs, DNA
from cogent.db.ensembl import Compara, Genome, HostAccount, Species

from homologsampler.util import missing_species_names, load_coord_names

def get_one2one_orthologs(compara, ref_genes, outpath, force_overwrite):
    """writes one-to-one orthologs of protein coding genes to outpath"""
    
    species = set(compara.Species)
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
        
            if syntenic is None or syntenic.getSpeciesSet() != species:
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
        
    print "Wrote %d files to %s" % (written, outpath)
    return

def get_latin_from_label(label):
    """returns latin name from the sequence label"""
    return label.split(':')[0]

def renamed_seqs(aln):
    """renames sequences to be just species common name"""
    new = []
    for seq in aln.Seqs:
        latin = get_latin_from_label(seq.Name)
        common = Species.getCommonName(latin)
        seq.Name = seq.Name.replace(latin, common)
        new.append((seq.Name, seq))
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
    species = set(compara.Species)
    common_names = map(Species.getCommonName, compara.Species)
    filler = LoadSeqs(data=[(n, 'N') for n in common_names], moltype=DNA)
    
    written = 0
    with click.progressbar(ref_genes,
        label="Finding 1to1 intron orthologs") as ids:
        for gene in ids:
            print "gene", gene
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
                if region.getSpeciesSet() != species:
                    continue
                
                if mask_features:
                    aln = region.getAlignment(feature_types=['gene', 'repeat', 'cpg'])
                    aln = with_masked_features(aln, reverse=gene.Location.Strand == -1)
                else:
                    aln = region.getAlignment()
                
                if aln is None:
                    continue
                
                aln = renamed_seqs(aln)
                alignments.append(aln)
                print aln.Names, common_names
            
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
            align.writeToFile(outfile_name)


def display_ensembl_alignment_table(compara):
    """prints the method_species_link table and then exits"""
    compara.method_species_links.Legend = \
    "Assign the desired value from method_link_species_set_id to the"\
    " method_clade_id argument"
    print compara.method_species_links
    exit(0)

@click.command()
@click.option('--ref', required=True, help='Reference species.')
@click.option('--species', required=True, help='Comma separated list of species names.')
@click.option('--release', required=True, help='Ensembl release.')
@click.option('--outdir', required=True, type=click.Path(resolve_path=True), help='Path to write files.')
@click.option('--coord_names', required=True, type=click.Path(resolve_path=True),
                    help='File containing chrom/coord names, one per line.')
@click.option('--introns', is_flag=True, help="Sample syntenic alignments of introns.")
@click.option('--method_clade_id', help="The align method ID to use.")
@click.option('--mask_features', is_flag=True, help="Intron masks repeats, exons, CpG islands.")
@click.option('--force_overwrite', is_flag=True, help="Overwrite existing files.")
@click.option('--show_align_methods', is_flag=True, help="Hows the align methods and exists.")
@click.option('--limit', type=int, default=0, help="Limit to this number of genes.")
@click.option('--test', is_flag=True)
def main(ref, species, release, outdir, coord_names, introns, method_clade_id, mask_features, force_overwrite, show_align_methods, limit, test):
    limit = limit or None
    try:
        acc = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
    except KeyError:
        warnings.warn("ENSEMBL_ACCOUNT environment variable not set, defaulting to UK sever. Slow!!")
        acc = None
    
    if introns and not method_clade_id:
        print "Must specify the method_clade_id in order to export introns"
        exit(-1)
    
    species = species.split(',')
    species_missing = missing_species_names(species)
    if species_missing:
        print "The following species names don't match an Ensembl record. Check spelling!"
        print species_missing
        exit(-1)
    
    if not ref in species:
        print "The reference species not in species names"
        exit(-1)
    
    compara = Compara(species, Release=release, account=acc)
    ref_genome = Genome(ref, Release=release, account=acc)
    
    if show_align_methods:
        display_ensembl_alignment_table(compara)
    
    chroms = load_coord_names(coord_names)
    chroms = chroms or None
    
    if not os.path.exists(outdir) and not test:
        os.makedirs(outdir)
        print "Created", outdir
    
    print "Sampling %s genes" % ref_genome
    # all_genes = ref_genome.getGenesMatching(BioType='protein_coding')
    all_genes = ref_genome.getGenesMatching(Symbol="BRCA1")
    ref_genes = []
    with click.progressbar(all_genes,
        label="Finding genes") as genes:
        for index, g in enumerate(genes):
            print g
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
