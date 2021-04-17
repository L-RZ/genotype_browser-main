import sqlite3
import sys, os
import argparse
import re
# make fake anno db from only VCF file
in_vcf_addr = sys.argv[1]
in_db_addr = sys.argv[2]

parser = argparse.ArgumentParser(description='create db or add variant to db')
parser.add_argument('-i', '--in_vcf', nargs='+', help='input vcf file')
parser.add_argument('-e', '--end_pos', help='output a file with: chr file_name end_pos')
parser.add_argument('-d', '--db', help='sqlite3 db file')
parser.add_argument('-m', '--mode', choices=['add', 'new'],
                    help='create: rm exist table and build a new table. add: only add data to table',
                    default='add')
args = parser.parse_args()

in_vcf_addr_l = args.in_vcf
in_db_addr  = args.db


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
chr_end_l = []
for in_vcf_addr in in_vcf_addr_l:
    f_in_vcf = open(in_vcf_addr)
    for line in f_in_vcf:
        if not line.startswith('#'):
            line_l = line.split('\t')
            snp_chr = line_l[0].replace('chr', '')
            snp_pos = int(line_l[1])
            snp_id = line_l[2]
            # snp_ref = line_l[3]
            # snp_alt = line_l[4]
            snp_info = line_l[7]
            snp_af = float(snp_info.split(';')[1].split('=')[1])

            snp_anno = (snp_id, snp_chr, snp_pos, 'NA', snp_af, 1, 0, 0, 'NA', 'NA', 'NA', 3)
            c.execute('INSERT INTO anno VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', snp_anno)
            c.execute('INSERT INTO chip VALUES (?,?)', (snp_id, 'exon'))
    vcf_file_name = os.path.basename(in_vcf_addr)
    chr_end_l.append('{}\t{}/{}\t{}\n'.format(snp_chr, 'data', vcf_file_name, snp_pos))
    print('Done', vcf_file_name)


conn.commit()
c.close()
conn.close()

f_chr_end = open(args.end_pos, 'w')
f_chr_end.write(''.join(chr_end_l))
f_chr_end.close()