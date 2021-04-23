import argparse
import sqlite3
import sys, os
parser = argparse.ArgumentParser(description='merge dbs to one db or add new db to old db')
parser.add_argument('-i', '--in_db', nargs='+', help='input new dbs e.g chr1.db chr2.db (sqlite3 db file)')
parser.add_argument('-d', '--db', help='sqlite3 db file')
parser.add_argument('-m', '--mode', choices=['add', 'new'],
                    help='create: rm exist table and build a new table. add: only add data to table',
                    default='add')
args = parser.parse_args()

in_new_db_addr_l = args.in_db
in_db_addr = args.db


conn = sqlite3.connect(in_db_addr)
c = conn.cursor()
if args.mode == 'new':
    c.execute('DROP TABLE IF EXISTS anno')
    c.execute(
        'CREATE TABLE anno (variant text, chr text, pos integer, rsid text, af real, info real, '
        'enrichment_nfsee_genomes real, enrichment_nfsee_exomes real, gene_most_severe text, most_severe text, '
        'consequence_gnomad text, in_data )')
    c.execute('DROP TABLE IF EXISTS chip')
    c.execute('CREATE TABLE chip (variant text, chip text)')

for each_new_db_addr in in_new_db_addr_l:
    c.execute('ATTACH DATABASE ? AS db2', (each_new_db_addr,))
    c.execute('INSERT INTO anno SELECT * FROM db2.anno')
    c.execute('INSERT INTO chip SELECT * FROM db2.chip')
    conn.commit()
    c.execute('DETACH DATABASE db2')
    print('Done', each_new_db_addr)
conn.commit()
c.close()
conn.close()
