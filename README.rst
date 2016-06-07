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

    $ homolog_sampler
    Usage: homolog_sampler [OPTIONS] COMMAND [ARGS]...

    Options:
      --ensembl_account TEXT  shell variable with MySQL account details, e.g.
                              export ENSEMBL_ACCOUNT='myhost.com jill jills_pass'
      -F, --force_overwrite   Overwrite existing files.
      --test
      --help                  Show this message and exit.

    Commands:
      one2one                 Command line tool for sampling homologous...
      show_align_methods      Shows the align methods in release...
      show_available_species  shows available species and Ensembl release...

To show available species
=========================

::

    $ homolog_sampler show_available_species
    Species available at: username:passwd@your.mysql.server
    =================================================================================================
    Release                                   Db Name                     Species         Common Name
    -------------------------------------------------------------------------------------------------
         20              ashbya_gossypii_core_20_73_1                        None                None
         20         aspergillus_clavatus_core_20_73_1        Aspergillus clavatus          A.clavatus
         20           aspergillus_flavus_core_20_73_1          Aspergillus flavus            A.flavus
         20        aspergillus_fumigatus_core_20_73_2       Aspergillus fumigatus         A.fumigatus...

To sample orthologs CDS
=======================

::

    $ homolog_sampler one2one --release 81 --ref human --species human,mouse,opossum --outdir sampled_cds

This command will write CDSs sequences to the directory ``sampled_cds`` as gzip compressed fasta files with the file prefix as the human Ensembl gene stable identifier (e.g. ``ENSG00000012048.fa.gz``).

To sample syntenic introns
==========================

We first need to know the Ensembl alignment method identifier. We get this as follows ::

    $ homolog_sampler show_align_methods --species human,mouse,opossum --release=81
    Align Methods/Clades
    ======================================================================================================
    method_link_species_set_id  method_link_id  species_set_id  align_method                   align_clade
    ------------------------------------------------------------------------------------------------------
                           788              10           36176         PECAN  23 amniota vertebrates Pecan
    ------------------------------------------------------------------------------------------------------
    Assign the desired value from method_link_species_set_id to the method_clade_id argument

We then use the ``one2one`` subcommand ::

    $ homolog_sampler one2one --release 81 --ref human --species human,mouse,opossum --outdir sampled_intron --introns --method_clade_id 788

