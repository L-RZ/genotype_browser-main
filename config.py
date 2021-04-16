import logging

sqlite_db = 'data/toy.db'

basic_info_file = 'data/genotype_browser_dummy_data_toy.txt'

vcf_files = {
    "imputed_data": [
        'data/toy_chr10.vcf.gz'

    ],
    "chip_data": [
        'data/toy_chr1.vcf.gz',
        'data/toy_chr2.vcf.gz',
        'data/toy_chr3.vcf.gz',
        'data/toy_chr4.vcf.gz',
        'data/toy_chr5.vcf.gz',
        'data/toy_chr6.vcf.gz',
        'data/toy_chr7.vcf.gz',
        'data/toy_chr8.vcf.gz',
        'data/toy_chr9.vcf.gz',
        'data/toy_chr10.vcf.gz'
    ]
}

write_dir = 'out'

log_level = logging.INFO
