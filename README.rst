##############
HomologSampler
##############

This is a command line tool for sampling related sequences from Ensembl. It requires PyCogent, Numpy, and a currently unpublished tool SciTrack. The latter basically logs commands, file inputs and outputs, to assist with reproducible research.

HomologSampler has one command line tool ``one2one``. At present, this tool is limited to sampling one-to-one orthologs from protein coding genes. It can either write out the protein coding sequences from the canonical CDS, or it can write out the Ensembl multiple sequence alignment of the entire gene with exons masked.

Increasing the number of biotypes and the homology relationships will be done if requested. But for now, it serves my purpose.

Basic usage is::

    $ one2one --help
    Usage: one2one [OPTIONS]

      Command line tool for sampling homologous sequences from Ensembl.

    Options:
      --ref TEXT              Reference species.  [required]
      --species TEXT          Comma separated list of species names.  [required]
      --release TEXT          Ensembl release.  [required]
      --outdir PATH           Path to write files.  [required]
      --coord_names PATH      File containing chrom/coord names, one per line.
                              [required]
      --introns               Sample syntenic alignments of introns.
      --method_clade_id TEXT  The align method ID to use.
      --mask_features         Intron masks repeats, exons, CpG islands.
      --force_overwrite       Overwrite existing files.
      --show_align_methods    Shows the align methods and exits.
      --limit INTEGER         Limit to this number of genes.
      --logfile_name TEXT     Name for log file, written to outdir.
      --test
      --help                  Show this message and exit.

Here's a sample usage::

    $ one2one --ref=Mouse --species="Mouse,Ground Squirrel,Rat" --release=81 --outdir=~/Desktop/Outbox/tmp/ --force_overwrite
