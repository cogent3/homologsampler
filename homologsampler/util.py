import sqlalchemy as sql
from cogent import LoadTable
from cogent.db.ensembl import Species
from cogent.db.ensembl.host import get_db_name

__author__ = "Gavin Huttley"
__copyright__ = "Copyright 2014, Gavin Huttley"
__credits__ = ["Gavin Huttley"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Gavin Huttley"
__email__ = "Gavin.Huttley@anu.edu.au"
__status__ = "Development"

def display_available_dbs(account):
    """displays the available Ensembl databases at the nominated host"""
    db_list = get_db_name(account=account, db_type='core')
    db_list += get_db_name(account=account, db_type='compara')
    rows = []
    for db_name in db_list:
        species_name = db_name.Species
        if species_name:
            common_name = Species.getCommonName(db_name.Species, level='ignore')
    
        if 'compara' in db_name.Name:
            species_name = common_name = '-'
        rows.append([db_name.Release, db_name.Name, species_name, common_name])

    table = LoadTable(header=["Release", "Db Name", "Species", "Common Name"], rows=rows, space=2)
    table = table.sorted(["Release", "Db Name"])
    table.Legend = "Values of 'None' indicate cogent does not have a value for that database name."
    return table

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
