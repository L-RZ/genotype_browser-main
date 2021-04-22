import sys
# make a dump data for each sample in VCF
in_addr = sys.argv[1]
out_addr = sys.argv[2]
header = ['ID', 'DEATH', 'SEX', 'AGE_AT_DEATH_OR_NOW', 'regionofbirthname', 'cohort', 'BATCH', 'CHIP', 'ARRAY']

f_in = open(in_addr)
f_out = open(out_addr, 'w')
f_out.write('\t'.join(header) + '\n')
for line in f_in:
    sample_id = line.strip()
    out_line_l = [sample_id, '-9', '-9',  'NA', 'NA', 'exome', 'exome', 'exome', '1']
    f_out.write('\t'.join(out_line_l) + '\n')
f_out.close()