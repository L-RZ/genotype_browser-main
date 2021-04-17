vcf_name=''


bcftools annotate -x ^INFO/AC,INFO/AF,INFO/AN,^FORMAT/GT,FORMAT/DP -O v ${vcf_name}.vcf.gz  | \
bcftools norm -m - -O v -o ${vcf_name}.s_m.vcf

python ~/fill_snpID_vcf.v0.1.py ${vcf_name}.s_m.vcf ${vcf_name}.s_m_id.vcf

python3 ~/add_db.v0.2.py -i ${vcf_name}.s_m_id.vcf \
 -e ${vcf_name}_end.txt -d exomes_test.db -m new



