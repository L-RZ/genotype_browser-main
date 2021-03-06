from flask import make_response
import gzip, pysam, threading, logging, timeit, os, sqlite3
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import bisect
import utils
import re

class Datafetch(object):

    def _init_tabix(self):
        self.tabix_files_imputed = defaultdict(lambda: [pysam.TabixFile(file, parser=None) for file in self.conf['vcf_files']['imputed_data']])
        self.tabix_files_chip = defaultdict(lambda: [pysam.TabixFile(file, parser=None) for file in self.conf['vcf_files']['chip_data']])


    def _init_chr_end(self):
        chr_end_dict = {}
        f_in_chr_end = open(self.conf['chr_end_pos_file'])
        init_type = 0
        for each_line in f_in_chr_end:
            each_line_l = each_line.strip().split()
            chr_id, chr_sub_addr, end_pos = each_line_l[:3]
            dat_type = each_line_l[3].replace('_data', '')  # data type or db type
            chr_id.replace('chr', '')
            if init_type == 0:
                i_ = 0
                init_type = dat_type
            elif init_type == dat_type:
                i_ += 1
            else:
                i_ = 0
                init_type = dat_type

            if dat_type in chr_end_dict:
                chr_end_dict[dat_type]['file_addr'].append(chr_sub_addr)
                if chr_id in chr_end_dict[dat_type]:
                    chr_end_dict[dat_type][chr_id]['vcf_index'].append(i_)
                    chr_end_dict[dat_type][chr_id]['end_pos'].append(int(end_pos))

                else:
                    chr_end_dict[dat_type][chr_id] = {'vcf_index': [i_],
                                                       'end_pos': [int(end_pos)]}


            else:
                chr_end_dict = {dat_type: {chr_id: {'vcf_index': [i_],
                                                    'end_pos': [int(end_pos)]},
                                           'file_addr': [chr_sub_addr]
                                           }
                                }
        for each_type in chr_end_dict:
            chr_end_dict[each_type]['tabix'] = defaultdict(
                lambda: [pysam.TabixFile(each_vcf, parser=None)
                         for each_vcf in chr_end_dict[each_type]['file_addr']]
            )


            # if dat_type in chr_end_dict:
            #     if chr_id in chr_end_dict[dat_type]:
            #         chr_end_dict[dat_type][chr_id]['vcf_index'].append(i_)
            #         chr_end_dict[dat_type][chr_id]['end_pos'].append(int(end_pos))
            #     else:
            #         chr_end_dict[dat_type] = {chr_id: {'vcf_index': [i_],
            #                                  'end_pos': [int(end_pos)]}
            #                                   }
            # else:
            #     chr_end_dict = {dat_type: {chr_id: {'vcf_index': [i_],
            #                                        'end_pos': [int(end_pos)]}}
            #                     }
            #
            # if chr_id in chr_end_dict:
            #     chr_end_dict[chr_id]['vcf_index'].append(i_)
            #     chr_end_dict[chr_id]['end_pos'].append(int(end_pos))
            # else:
            #     chr_end_dict = {chr_id: {'vcf_index': [i_],
            #                              'end_pos': [int(end_pos)]}
            #                     }
        self.chr_end_dict = chr_end_dict

    def _init_db(self):
        self.conn = defaultdict(lambda: sqlite3.connect(self.conf['sqlite_db']))

    def _init_info(self, select_chip):
        info = pd.read_csv(self.conf['basic_info_file'], sep='\t').fillna('NA')
        if select_chip:
            info.index = info['ID']
            # first_vcf = self.chr_end_dict['chip']['file_addr'][0]
            # tabix_iter = defaultdict(lambda: [pysam.TabixFile(first_vcf, parser=None)])[threading.get_ident()][0]

            tabix_iter = self.chr_end_dict['chip']['tabix'][threading.get_ident()][0]
            # tabix_iter = self.tabix_files_chip[threading.get_ident()][0]
            h = tabix_iter.header[len(tabix_iter.header) - 1].split('\t')
            chip_samples = h[9:]
            info = info.loc[chip_samples, :]

        # make a copy for writing data out with male/female texts
        info_orig = info.copy()
        # info_orig['SEX'] = np.where(info_orig['SEX'] == 1, 'female', 'male')
        info_orig['DEATH'] = np.where(info_orig['DEATH']=="NA", np.nan, info_orig['DEATH'])
        info_orig['DEATH'] = info_orig['DEATH'].astype('Int64')
        
        # change cohort/region names to index
        cohort_list = sorted(list(set(info['cohort'].tolist())))
        cohort_idx = {cohort:i for i,cohort in enumerate(cohort_list)}
        info['cohort'] = info['cohort'].map(cohort_idx)
        
        region_list = sorted(list(set(info['regionofbirthname'].tolist())))
        region_idx = {region:i for i,region in enumerate(region_list)}
        info['regionofbirth'] = info['regionofbirthname'].map(region_idx)
        info = info.drop(columns=['ID', 'regionofbirthname'])
        info.reset_index(drop=True, inplace=True)
        info_columns = info.columns.tolist()

        # remove index
        info.reset_index(drop=True, inplace=True)
        info_orig.reset_index(drop=True, inplace=True)
        return info, info_orig, cohort_list, region_list, info_columns

    def __init__(self, conf):
        self.conf = conf
        # self._init_tabix()
        self._init_chr_end()
        self._init_db()
        self.info, self.info_orig, self.cohort_list, self.region_list, self.info_columns = self._init_info(select_chip=False)
        self.info_chip, self.info_orig_chip, self.cohort_list_chip, self.region_list_chip, self.info_columns_chip = self._init_info(select_chip=True)


    def _get_annotation(self, chr, pos, ref, alt, data_type):
        in_data = 1 if data_type == 'imputed' else 2
        if self.conn[threading.get_ident()].row_factory is None:
            self.conn[threading.get_ident()].row_factory = sqlite3.Row
        c = self.conn[threading.get_ident()].cursor()         
        query = 'SELECT * FROM anno WHERE variant = "%s" AND (in_data=%s OR in_data=3);' % (str(chr) + ':' + str(pos) + ':' + ref + ':' + alt, in_data)
        c.execute(query)
        res = c.fetchone()
        if len(res) > 0:
            return dict(res)
        else:
            raise utils.NotFoundException(str(chr) + '-' + str(pos) + '-' + ref + '-' + alt)

    # each chr will include multi vcf files, get each vcf file end pos and build a index bin.
    # using bisect.bisect to locate which bin(vcf file) include the SNP and fetch on it.
    def _get_genotype_data(self, chr, pos, ref, alt, data_type):
        chr_var = chr if chr != 23 else 'X'
        vcf_index_l = self.chr_end_dict[data_type][str(chr_var)]['vcf_index']
        end_pos_l = self.chr_end_dict[data_type][str(chr_var)]['end_pos']
        chr_bin_index = bisect.bisect_left(end_pos_l, pos)
        vcf_index = vcf_index_l[chr_bin_index]

        # if data_type == 'imputed':
        #     tabix_iter = self.tabix_files_imputed[threading.get_ident()][vcf_index].fetch('chr'+str(chr_var), pos-1, pos)
        # else:
        #     tabix_iter = self.tabix_files_chip[threading.get_ident()][vcf_index].fetch('chr'+str(chr_var), pos-1, pos)
        tabix_iter = self.chr_end_dict[data_type]['tabix'][threading.get_ident()][vcf_index].fetch('chr'+str(chr_var), pos-1, pos)
        # vcf_addr = self.chr_end_dict[data_type]['file_addr'][vcf_index]
        # tabix_iter = defaultdict(lambda: [pysam.TabixFile(vcf_addr, parser=None)])[threading.get_ident()][0].fetch('chr'+str(chr_var), pos-1, pos)
        var_data = None
        for row in tabix_iter:
            data = row.split('\t')            
            if data[3] == ref and data[4] == alt:
                var_data = data
                break
        return var_data[9:] if var_data is not None else None

    def _get_chips(self, chr, pos, ref, alt):
        if self.conn[threading.get_ident()].row_factory is None:
            self.conn[threading.get_ident()].row_factory = sqlite3.Row
        c = self.conn[threading.get_ident()].cursor()
        res = c.execute('SELECT * FROM chip WHERE variant = ?', [str(chr) + ':' + str(pos) + ':' + ref + ':' + alt])
        chips = set([r['chip'] for r in res])
        if len(chips) > 0:
            return chips
        else: # variant doesn't exist in the table if it's not on any chip
            return set()
    
    def _filter(self, df, filters, chips):
        _df = df.copy()
        if 'alive' in filters:
            if filters['alive'] == 'alive':
                _df = _df.loc[_df['DEATH'] == 0]
            elif filters['alive'] == 'dead':
                _df = _df.loc[_df['DEATH'] == 1]
        if 'sex' in filters:
            if filters['sex'] == 'female':
                _df = _df.loc[_df['SEX'] == 1]
            elif filters['sex'] == 'male':
                _df = _df.loc[_df['SEX'] == 0]
        if 'array' in filters:
            if filters['array'] == 'finngen':
                _df = _df.loc[_df['ARRAY'] == 1]
            elif filters['array'] == 'legacy':
                _df = _df.loc[_df['ARRAY'] == 0]
        if 'impchip' in filters:
            if filters['impchip'] == 'chip':
                _df = _df.loc[_df['BATCH'].isin(chips)]
            elif filters['impchip'] == 'imp':
                _df = _df.loc[~_df['BATCH'].isin(chips)]
        return _df

    def _is_homozygous_alt(self, gt):
        return gt == '1|1' or gt == '1/1'

    def _is_heterozygous(self, gt):
        return gt == '1|0' or gt == '1/0' or gt == '0|1' or gt == '0/1'

    def _is_homozygous_wt(self, gt):
        return gt == '0|0' or gt =='0/0'

    def _gt_passes_threshold(self, gt_probability, gt_probability_thres):
        return gt_probability >= float(gt_probability_thres)

    def _get_het_hom_index(self, data, index, use_gt, gp_thres, calc_info):
        het_i = []
        hom_alt_i = []
        sum_eij = 0
        sum_fij_minus_eij2 = 0
        wt_hom_i = []
        missing_i = []
        for i in index:
            #GT:DS:GP
            #0|0:0:1,0,0
            s = data[i].split(':')
            gt = s[0]
            if not use_gt or calc_info:
                gp = [float(p) for p in s[2].split(',')]
            if calc_info:
                dosage = float(s[1])
                sum_eij = sum_eij + dosage
                fij_minus_eij2 = 4*gp[2] + gp[1] - dosage*dosage
                sum_fij_minus_eij2 = sum_fij_minus_eij2 + fij_minus_eij2
            if use_gt:
                if self._is_homozygous_alt(gt):
                    hom_alt_i.append(i)
                elif self._is_heterozygous(gt):
                    het_i.append(i)
                elif self._is_homozygous_wt(gt):
                    wt_hom_i.append(i)
                else:
                    missing_i.append(i)
            else:
                if self._is_homozygous_alt(gt) and self._gt_passes_threshold(gp[2], gp_thres):
                    hom_alt_i.append(i)
                elif self._is_heterozygous(gt) and self._gt_passes_threshold(gp[1], gp_thres):
                    het_i.append(i)
                elif self._is_homozygous_wt(gt) and self._gt_passes_threshold(gp[0], gp_thres):
                    wt_hom_i.append(i)
                else:
                    missing_i.append(i)

        if calc_info and len(index)>0:
            theta_hat = sum_eij / (2*len(index))
            info = 1 if theta_hat == 0 or theta_hat == 1 else 1 - sum_fij_minus_eij2 / (2*len(index)*theta_hat*(1-theta_hat))
        else:
            info = -1

        return (het_i, hom_alt_i, wt_hom_i, missing_i, info)
    
    def _aggregate_het_hom(self, het, hom, full, data_type):
        agg = {'regions': {}, 'cohorts': {}}
        if data_type == 'imputed':
            cohort_list = self.cohort_list
            region_list = self.region_list
        else:
            cohort_list = self.cohort_list_chip
            region_list = self.region_list_chip

        for type in [('regions', 'regionofbirth', region_list), ('cohorts', 'cohort', cohort_list)]:
            # gt_count: list of length 2 (het/hom), each is a dict from region/cohort index to het/hom count
            gt_count = [het.groupby(type[1]).size().to_dict(), hom.groupby(type[1]).size().to_dict()]
            # add zeros
            gt_count = [
                [gt_count[0][i] if i in gt_count[0] else 0 for i in range(len(type[2]))],
                [gt_count[1][i] if i in gt_count[1] else 0 for i in range(len(type[2]))]
            ]
            num_indiv = full.groupby(type[1]).size().to_dict()
            num_indiv = [num_indiv[i] if i in num_indiv else 0 for i in range(len(type[2]))]
            # calculate allele frequency
            af = [(cnt+2*gt_count[1][i])/num_indiv[i]/2 if num_indiv[i] > 0 else 0 for i,cnt in enumerate(gt_count[0])]
            agg[type[0]] = {'names': type[2], 'gt_counts': gt_count, 'num_indiv': num_indiv, 'af': af}
        return agg
    
    def _count_gt(self, data, filters, chips, data_type):
        if data_type == 'imputed':
            filtered_basic_info = self._filter(self.info, filters, chips)
            info_df = self.info
            info_columns = self.info_columns
        else:
            filtered_basic_info = self._filter(self.info_chip, filters, chips)
            info_df = self.info_chip
            info_columns = self.info_columns_chip

        id_index = list(filtered_basic_info.index)
        het_i = []
        hom_i = []
        wt_hom_i = []
        missing_i = []
        for d in data:
            if data_type == 'imputed':
                het_i_d, hom_i_d, wt_hom_i_d, missing_i_d, info = self._get_het_hom_index(d, id_index, filters['gtgp'] == 'gt', filters['gpThres'], len(data) == 1)
            else:
                het_i_d, hom_i_d, wt_hom_i_d, missing_i_d, info = self._get_het_hom_index(d, id_index, True, None, False)
            het_i.extend(het_i_d)
            hom_i.extend(hom_i_d)
            wt_hom_i.extend(wt_hom_i_d)
            missing_i.extend(missing_i_d)
        het_cnt = Counter(het_i)        
        het_i = set(het_i)
        hom_i = set(hom_i)
        wt_hom_i = set(wt_hom_i)
        missing_i = set(missing_i)
        # maybe treat multiheterozygotes as homozygotes
        if 'hethom' in filters and filters['hethom']:
            multihet = [i for i in het_cnt if het_cnt[i] > 1]
            hom_i = set().union(hom_i, multihet)
        # if an individual is homozygous for a variant, don't count as heterozygous for other variants
        het = info_df.iloc[[i for i in het_i if i not in hom_i]]
        hom = info_df.iloc[list(hom_i)]
        wt_hom = info_df.iloc[list(wt_hom_i)]
        missing = info_df.iloc[list(missing_i)]
        agg = self._aggregate_het_hom(het, hom, filtered_basic_info, data_type)
        total_af = (len(het) + 2*len(hom))/len(filtered_basic_info)/2 if len(filtered_basic_info) > 0 else -1
        het = het.to_numpy().T.tolist()
        hom = hom.to_numpy().T.tolist()

        # add wt hom and missing data
        wt_hom = wt_hom.to_numpy().T.tolist()
        missing = missing.to_numpy().T.tolist()        
        return {
            'het': het if len(het) > 0 else [[] for i in info_columns],
            'hom_alt': hom if len(hom) > 0 else [[] for i in info_columns],
            'wt_hom': wt_hom if len(wt_hom) > 0 else [[] for i in info_columns],
            'missing': missing if len(missing) > 0 else [[] for i in info_columns],
            'columns': info_columns,
            'agg': agg,
            'total_af': total_af,
            'info': info,
            'total_indiv': len(filtered_basic_info),
            'filters': filters
        }
    
    def _count_gt_for_write(self, variants, data, filters, chips, data_type):
        start_time = timeit.default_timer()
        df_list = []
        if data_type == 'imputed':
            filtered_basic_info = self._filter(self.info, filters, chips)
            info_orig = self.info_orig
        else:
            filtered_basic_info = self._filter(self.info_chip, filters, chips)
            info_orig = self.info_orig_chip
        id_index = list(filtered_basic_info.index) 
        for i, d in enumerate(data):
            # get indices of het/hom individuals in genotype data
            if data_type == 'imputed':
                het_i, hom_alt_i, wt_hom_i, missing_i, info = self._get_het_hom_index(d, id_index, filters['gtgp'] == 'gt', filters['gpThres'], len(data) == 1)
            else:
                het_i, hom_alt_i, wt_hom_i, missing_i, info = self._get_het_hom_index(d, id_index, True, None, False)
            gt = [element.split(':')[0] for element in d ]
            gt_arr = np.array(gt)
            gt_het = list(gt_arr[het_i])
            gt_hom = list(gt_arr[hom_alt_i])
            gt_wt_hom = list(gt_arr[wt_hom_i])
            gt_missing = list(gt_arr[missing_i])

            # subset dataframe for het/hom individuals, copy needed as this will be mutated
            het = info_orig.iloc[het_i].copy()
            hom_alt = info_orig.iloc[hom_alt_i].copy()
            wt_hom = info_orig.iloc[wt_hom_i].copy()
            missing = info_orig.iloc[missing_i].copy()

            # extract gt probs and gts
            if data_type == 'imputed':
                gt_probs = [element.split(':')[2] for element in d ]
                gt_probs_arr = np.array(gt_probs)
                gt_probs_het = list(gt_probs_arr[het_i])
                gt_probs_hom = list(gt_probs_arr[hom_alt_i])
                gt_probs_wt_hom = list(gt_probs_arr[wt_hom_i])
                gt_probs_missing = list(gt_probs_arr[missing_i])
                het['three_gt_probs'] = gt_probs_het
                hom_alt['three_gt_probs'] = gt_probs_hom
                wt_hom['three_gt_probs'] = gt_probs_wt_hom
                missing['three_gt_probs'] = gt_probs_missing                

            # add main gt and probs for three genotypes (and missing data for raw chip)
            het['gt'] = gt_het
            hom_alt['gt'] = gt_hom
            wt_hom['gt'] = gt_wt_hom
            missing['gt'] = gt_missing

            # if specified the type of variants to be saved
            if 'hethom' in filters:
                if filters['hethom'] == 'hom':
                    df = hom_alt
                elif filters['hethom'] == 'het':
                    df = het
                elif filters['hethom'] == 'wt_hom':
                    df = wt_hom
                else:
                    # append all data frames: wt, het, hom, wt_hom
                    df = hom_alt.append(het, ignore_index=True).append(wt_hom, ignore_index=True).append(missing, ignore_index=True)

            # append data frames
            df = self._filter(df, filters, chips)
            df['variant'] = variants[i].replace('-', ':')
            df_list.append(df)
        
        elapsed = timeit.default_timer() - start_time
        return pd.concat(df_list)

    def get_variants(self, variants, filters, data_type):
        start_time = timeit.default_timer()
        filters = {k:(True if v=="true" else v) for k,v in filters.items()}
        filters = {k:(False if v=="false" else v) for k,v in filters.items()}
        if 'gpThres' in filters:
            filters['gpThres'] = float(filters['gpThres'])
        vars_data = []
        anno = []
        vars = []
        for variant in variants.split(','):
            chr, pos, ref, alt = utils.parse_variant(variant)
            var_data = self._get_genotype_data(chr, pos, ref, alt, data_type)
            if var_data is not None:
                vars_data.append(var_data)
                vars.append('-'.join([str(s) for s in [chr, pos, ref, alt]]))
                anno.append(self._get_annotation(chr, pos, ref, alt, data_type))
        if len(vars_data) == 0:
            raise utils.NotFoundException()
        fetch_time = timeit.default_timer() - start_time
        start_time = timeit.default_timer()
        chips = self._get_chips(chr, pos, ref, alt) if len(vars_data) == 1 else set()
        if data_type == 'chip':
            if 'impchip' in filters:
                del filters['impchip']
            if 'gtgp' in filters:
                del filters['gtgp']
        data = self._count_gt(vars_data, filters, chips, data_type)
        munge_time = timeit.default_timer() - start_time
        return {
            'variants': vars,
            'annotation': anno,
            'data': data,
            'time': {
                'fetch': fetch_time,
                'munge': munge_time
            },
            'data_type': data_type
        }

    def check_var_in_chip(self, variant):
        chr, pos, ref, alt = utils.parse_variant(variant)
        var_data = self._get_genotype_data(chr, pos, ref, alt, 'chip')
        if var_data is not None:
            return True
        else:
            return False

    def write_variants(self, variants, filters, data_type):
        vars_data = []
        for variant in variants.split(','):
            chr, pos, ref, alt = utils.parse_variant(variant)
            var_data = self._get_genotype_data(chr, pos, ref, alt, data_type)
            if var_data is not None:
                vars_data.append(var_data)
        chips = self._get_chips(chr, pos, ref, alt) if len(vars_data) == 1 else set()
        data = self._count_gt_for_write(variants.split(','), vars_data, filters, chips, data_type)
        if data_type == 'imputed':
            del filters['data_type']
            filename = variants.replace(',', '_') + '__imputed_data__' + '_'.join([k+'_'+v for k,v in filters.items()]) + '.tsv'
        else:
            for key in ['array', 'impchip', 'data_type']:
                del filters[key]
            # filename = variants.replace(',', '_') + '__rawchip_data__' + '_'.join([k+'_'+v for k,v in filters.items()]) + '.tsv'
            filename = variants.replace(',', '_') + '__raw_data__' + '_'.join(
                [k + '_' + v for k, v in filters.items() if k not in ('alive', 'sex')]) + '.tsv'
        # data = data.drop(columns=['AGE_AT_DEATH_OR_NOW'])
        # data['SEX'] = np.where(data['SEX'] == 1, 'female', 'male')
        data = data.drop(columns=['AGE_AT_DEATH_OR_NOW', 'DEATH', 'regionofbirthname',
                                   'BATCH', 'CHIP', 'ARRAY'])
        # output Genotype
        var_id_dict = {}
        for each_var_id in set(data['variant'].to_list()):
            snp_chr, snp_pos, snp_ref, snp_alt = each_var_id.split(':')
            var_id_dict[each_var_id] = {'./.': './.',
                                        '1/1': '{}/{}'.format(snp_alt, snp_alt),
                                        '0/1': '{}/{}'.format(snp_ref, snp_alt),
                                        '1/0': '{}/{}'.format(snp_ref, snp_alt),
                                        '0/0': '{}/{}'.format(snp_ref, snp_ref),
                                        }
        gt_real = [var_id_dict[gt_var[1]][gt_var[0]]
            for gt_var in zip(data['gt'].to_list(), data['variant'].to_list())]
        data['genotype'] = gt_real
        try:
            data.to_csv(sep='\t', index=False, na_rep='NA')
            output = make_response(data.to_csv(sep='\t', index=False, na_rep='NA'))

            output.headers["Content-Disposition"] = "attachment; filename=" + filename
            output.headers["Content-type"] = "text/tab-separated-values"
            return output

        except Exception as e:
            # TODO should return non-200 maybes
            return {'status': 'failed', 'message': str(e)}

    def get_gene_variants(self, gene, data_type):
        if self.conn[threading.get_ident()].row_factory is None:
            self.conn[threading.get_ident()].row_factory = sqlite3.Row
        c = self.conn[threading.get_ident()].cursor()
        c.execute('SELECT * FROM genes WHERE gene_name = ? OR gene_name = ?', [gene, gene.upper()])
        gene_db = [dict(row) for row in c.fetchall()]
        if len(gene_db) == 0:
            raise utils.NotFoundException()
        else:
            vars_db = self.get_genomic_range_variants(gene_db[0]['chr'], gene_db[0]['start'], gene_db[0]['end'], data_type)
            res_vars = vars_db['data']
        # drop columns we don't show
        exclude_cols = ['gene_most_severe', 'consequence_gnomad', 'chr',  'info', 'in_data',  'enrichment_nfsee_genomes', 'enrichment_nfsee_exomes']
        cols = [col for col in res_vars[0].keys() if col not in exclude_cols]
        return {
            'gene': gene,
            'columns': cols,
            'data': res_vars,
            'data_type': data_type
        }

    def get_genomic_range_variants(self, chr, start, end, data_type):
        in_data = 1 if data_type == 'imputed' else 2
        if self.conn[threading.get_ident()].row_factory is None:
            self.conn[threading.get_ident()].row_factory = sqlite3.Row
        c = self.conn[threading.get_ident()].cursor()
        query = 'SELECT * FROM anno WHERE chr="%s" AND pos>=%s AND pos<=%s AND (in_data=%s OR in_data=3);' % (chr, start, end, in_data)
        c.execute(query)
        res = [dict(row) for row in c.fetchall()]
        if len(res) == 0:
            raise utils.NotFoundException()
        data = []
        for item in res:
            item['variant'] = '-'.join(item['variant'].split(':'))
            data.append(item)
        genomic_range = "%s:%s-%s" % (chr, start, end)                     
        exclude_cols = ['gene_most_severe', 'consequence_gnomad', 'chr', 'info', 'in_data', 'enrichment_nfsee_genomes', 'enrichment_nfsee_exomes']
        cols = [col for col in data[0].keys() if col not in exclude_cols]
        return {
            'range': genomic_range,
            'columns': cols,
            'data': data,
            'data_type': data_type
        }