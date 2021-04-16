import sys

in_vcf = sys.argv[1]
out_vcf = sys.argv[2]

f_in_vcf = open(in_vcf)
f_out_vcf = open(out_vcf, 'w')
for line in f_in_vcf:
    if line.startswith('#'):
        f_out_vcf.write(line)
    else:
        line_l = line.split()
        snp_chr = line_l[0].replace('chr', '')
        snp_pos = line_l[1]
        snp_ref = line_l[3]
        snp_alt = line_l[4]

        snp_new_id = ':'.join([snp_chr, snp_pos, snp_ref, snp_alt])
        line_l[2] = snp_new_id
        f_out_vcf.write('\t'.join(line_l))
f_out_vcf.close()

