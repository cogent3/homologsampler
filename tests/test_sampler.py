import os, shutil

import click
from click.testing import CliRunner

from cogent import LoadTable

from homologsampler import cli

def _parse_db_display(output, columns):
    """finds the table display and accumulates the content"""
    result = output.splitlines()
    has_table = False
    records = []
    header = []
    for index, line in enumerate(result):
        if not header and columns[0] in line:
            header = columns
            break
    
    if header:
        rows = []
        for i in range(index+2, len(result)):
            line = result[i].strip()
            if line.startswith('----------'):
                break
            
            line = line.split()
            rows.append(line[:len(columns)])
        table = LoadTable(header=header, rows=rows)
    else:
        table = None
    
    return table

def test_show_databases():
    """exercises showing all databases"""
    # any ensembl release
    runner = CliRunner()
    r = runner.invoke(cli, ["show_available_species"])
    assert r.exit_code == 0
    result = _parse_db_display(r.output, ["Release", "Db Name"])
    assert result.Shape[0] > 0

def test_show_release_databases():
    """exercises showing all databases"""
    # any ensembl release
    runner = CliRunner()
    r = runner.invoke(cli, ["show_available_species", "--release=81"])
    result = _parse_db_display(r.output, ["Release", "Db Name"])
    db_val = result.getDistinctValues('Release')
    assert db_val == set(['81'])

def test_show_align_methods():
    """show align methods works"""
    runner = CliRunner()
    r = runner.invoke(cli, ["show_align_methods", "--species=human,chimp", "--release=81"])
    assert r.exit_code == 0
    result = _parse_db_display(r.output, ["method_link_species_set_id"])
    assert result.Shape[0] > 0

def test_one2one_cds():
    """samples a CDS sequence"""
    runner = CliRunner()
    r = runner.invoke(cli, ["--test=1", "one2one", "--species=human,mouse", "--release=81",
                "--ref=human", "--outdir=delme"])
    
    assert ">Human" in r.output and ">Mouse" in r.output
    
def test_one2one_intron():
    """samples intron sequence"""
    runner = CliRunner()
    r = runner.invoke(cli, ["--test=14", "one2one", "--species=human,chimp", "--release=81",
                "--ref=human", "--outdir=delme", "--introns", "--method_clade_id=688"])
    assert ">Human" in r.output and ">Chimp" in r.output

