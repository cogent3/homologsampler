import sqlalchemy as sql
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

def get_chrom_names(ref_species, compara):
    """returns the list of chromosome names"""
    genome_db = compara.ComparaDb.getTable("genome_db")
    dnafrag = compara.ComparaDb.getTable("dnafrag")
    joined = genome_db.join(dnafrag, onclause=genome_db.c.genome_db_id==dnafrag.c.genome_db_id)
    condition = sql.and_(dnafrag.c.coord_system_name=="chromosome",
                    genome_db.c.name==Species.getEnsemblDbPrefix(ref_species),
                    dnafrag.c.is_reference==1)
    query = sql.select([dnafrag.c.name], condition).select_from(joined)
    chroms = [r[0] for r in query.execute()]
    return chroms

def load_coord_names(infile_path):
    """loads chrom names, assumes separate name per file"""
    with open(infile_path) as infile:
        coord_names = [l.strip() for l in infile]
    return coord_names
