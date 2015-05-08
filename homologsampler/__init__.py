import os, gzip, warnings
from collections import Counter

import click

from cogent import LoadSeqs
from cogent.db.ensembl import Compara, Genome, HostAccount

from homologsampler.util import missing_species_names, load_coord_names

def get_one2one_orthologs(compara, ref_genes, outpath):
    """returns Table of one-to-one orthologs of protein coding genes"""
    
    species = set(compara.Species)
    written = 0
    with click.progressbar(ref_genes,
        label="Finding 1to1 orthologs") as ids:
        for gene in ids:
            outfile_name = os.path.join(outpath, "%s.fa.gz" % gene.StableId)
            if os.path.exists(outfile_name):
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

@click.command()
@click.option('--ref', required=True, help='Reference species.')
@click.option('--species', required=True, help='Comma separated list of species names.')
@click.option('--release', required=True, help='Ensembl release.')
@click.option('--outdir', required=True, type=click.Path(resolve_path=True), help='Path to write files.')
@click.option('--coord_names', required=True, type=click.Path(resolve_path=True),
    help='File containing chrom/coord names, one per line.')
@click.option('--test', is_flag=True)
def main(ref, species, release, outdir, coord_names, test):
    try:
        acc = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
    except KeyError:
        warnings.warn("ENSEMBL_ACCOUNT environment variable not set, defaulting to UK sever. Slow!!")
        acc = None
    
    species = species.split(',')
    species_missing = missing_species_names(species)
    if species_missing:
        print "The following species names don't match an Ensembl record. Check spelling!"
        print species_missing
        exit(-1)
    
    if not ref in species:
        print "The reference species not in species names"
        exit(-1)
    
    chroms = load_coord_names(coord_names)
    chroms = chroms or None
    
    if not os.path.exists(outdir) and not test:
        os.makedirs(outdir)
        print "Created", outdir
    
    compara = Compara(species, Release=release, account=acc)
    ref_genome = Genome(ref, Release=release, account=acc)
    
    print "Sampling %s genes" % ref_genome
    all_genes = ref_genome.getGenesMatching(BioType='protein_coding')
    
    if chroms:
        ref_genes = [g for g in all_genes if g.Location.CoordName in chroms]
    else:
        ref_genes = [g for g in all_genes]
    
    print "Getting orthologs"
    get_one2one_orthologs(compara, ref_genes, outpath)


if __name__ == "__main__":
    main()
