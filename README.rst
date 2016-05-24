##############
HomologSampler
##############

This is a command line tool for sampling related sequences from Ensembl. It requires PyCogent, Numpy, and a currently unpublished tool SciTrack. The latter basically logs commands, file inputs and outputs, to assist with reproducible research.

HomologSampler has one command line tool ``one2one``. At present, this tool is limited to sampling one-to-one orthologs from protein coding genes. It can either write out the protein coding sequences from the canonical CDS, or it can write out the Ensembl multiple sequence alignment of the entire gene with exons masked.

Increasing the number of biotypes and the homology relationships will be done if requested. But for now, it serves my purpose.

************
Installation
************

Because we rely on PyCogent, whose install depends on numpy in a way that standard package installers don't cope with, the following is the required order of statements. (Note, I assume you already have pip_ installed.)

Install numpy

::
    
    $ pip install numpy --upgrade

Then install HomologSampler directly from the bitbucket, specifying to follow dependency links

::
    
    $ pip install --process-dependency-links hg+ssh://hg@bitbucket.org/gavin.huttley/homologsampler

.. _pip: https://pip.pypa.io/en/stable/installing/


**************
Basic usage is
**************

::

    $ one2one --help
    Usage: one2one [OPTIONS]

      Command line tool for sampling homologous sequences from Ensembl.

    Options:
      --ref TEXT                Reference species.
      --species TEXT            Comma separated list of species names.
      --release TEXT            Ensembl release.
      --outdir PATH             Path to write files.
      --ensembl_account TEXT    shell variable with MySQL account details, e.g.
                                export ENSEMBL_ACCOUNT='myhost.com jill
                                jills_pass'
      --coord_names PATH        File containing chrom/coord names to restrict
                                sampling to, one per line.
      --introns                 Sample syntenic alignments of introns.
      --method_clade_id TEXT    The align method ID to use, required if sampling
                                introns.
      --mask_features           Intron masks repeats, exons, CpG islands.
      --force_overwrite         Overwrite existing files.
      --show_align_methods      Shows the align methods and exits.
      --show_available_species  Shows the available db's at ENSEMBL_ACCOUNT.
      --limit INTEGER           Limit to this number of genes.
      --logfile_name TEXT       Name for log file, written to outdir.
      --test
      --help                    Show this message and exit.
