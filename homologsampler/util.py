from cogent import LoadTable
from cogent.db.ensembl import Species

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2014, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

def missing_species_names(names):
    '''returns a Table of missing species names, or None'''
    missing = []
    for name in names:
        n = Species.getSpeciesName(name)
        if n == 'None':
            missing.append([name])
    
    if missing:
        result = LoadTable(header=["MISSING SPECIES"], rows=missing)
    else:
        result = None
    return result

def load_coord_names(infile_path):
    """loads chrom names, assumes separate name per file"""
    with open(infile_path) as infile:
        coord_names = [l.strip() for l in infile]
    return coord_names