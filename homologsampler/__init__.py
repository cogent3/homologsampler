import os, gzip
from collections import Counter

import click

from cogent import LoadSeqs
from cogent.db.ensembl import Compara, HostAccount, Species

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
@click.option('--outdir', required=True, type=click.Path(resolve_path=True), help='Path to write.')
@click.option('--test', is_flag=True)
def main(outdir, test):
    acc = HostAccount(*os.environ['ENSEMBL_ACCOUNT'].split())
    if not os.path.exists(outdir) and not test:
        os.makedirs(outdir)
        print "Created", outpath
    
    compara = Compara(['human', 'dog', 'mouse', 'opossum', 'platypus'], Release=76, account=acc)
    
    print "Sampling genes"
    all_genes = compara.Human.getGenesMatching(BioType='protein_coding')
    chroms = map(str, range(1, 23)) + ['X', 'Y']
    ref_genes = [g for g in all_genes if g.Location.CoordName in chroms]
    print "Getting orthologs"
    get_one2one_orthologs(compara, ref_genes, outpath)


if __name__ == "__main__":
    main()
