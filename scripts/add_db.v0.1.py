import sqlite3
import sys
import re
# make fake anno db from only VCF file
in_vcf_addr = sys.argv[1]
in_db_addr = sys.argv[2]

conn = sqlite3.connect(in_db_addr)
c = conn.cursor()
c.execute('DROP TABLE IF EXISTS anno')
c.execute(
    'CREATE TABLE anno (variant text, chr text, pos integer, rsid text, af real, info real, enrichment_nfsee_genomes real, enrichment_nfsee_exomes real, gene_most_severe text, most_severe text, consequence_gnomad text, in_data )')

c.execute('DROP TABLE IF EXISTS chip')

c.execute('CREATE TABLE chip (variant text, chip text)')
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

conn.commit()
c.close()
conn.close()
