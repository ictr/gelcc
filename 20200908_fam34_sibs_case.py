#!/usr/bin/env python

import vcf
import pandas as pd
import csv
import argparse

def calculate_core_count(vcf_path, output_path, phen_path):
    phen_df = pd.read_csv(phen_path, sep = '\t')
    # create a dictionary that maps CID_CIND to samples names in vcf
    cind2sample_name = dict(list(zip(phen_df['CID_CIND'], phen_df['Musolf_sampleID'])))

    vcf_reader = vcf.Reader(open(vcf_path))

    out_file = open(output_path, 'w')
    carrier_counts_H = ['MIC_REF_COUNT','CORE_CARRIER_COUNT','NON_CORE_CARRIER_COUNT']
    identifiers_H = ['CHROM', 'POS', 'ID', 'REF','ALT']
##module1  fam_gts_H is list of sequenced fam102 samples
    fam_gts_H = ['102_101','102_112','102_121','102_124','102_132','102_136','102_137','102_180','102_301','102_311','102_313','102_346']
    headers =  identifiers_H + carrier_counts_H + fam_gts_H
    tsv_writer = csv.writer(out_file, delimiter='\t', lineterminator='\n')
    tsv_writer.writerow(headers)
    for variant in vcf_reader:
        are_carriers = lambda x: list(map(is_carrier, get_sample_gts(variant, x)))
        are_refs = lambda x: list(map(is_ref, get_sample_gts(variant, x)))
        c2s = lambda cinds: [cind2sample_name[s] for s in cinds]
        
        # observed number of married-in controls (unrelated to sibs)
##module2  mic_ids is c2s(list of sequenced fam102 MCIs)
        mic_ids = c2s(['102_346', '102_124'])
        mic_ref_count = sum(are_refs(mic_ids))
##module3 seq_core_ids is c2s(list of sequenced fam102 LCS (affected sibs))
        seq_core_ids  = c2s(['102_112'])
        seq_core_carrier_count = sum(are_carriers(seq_core_ids))
##module4: defines desc_ids for the 2 inferred core LCs in fam102
# number of inferred core carriers based on sequence of descendants is:  inf-core-LC = (313, 311, 301) +  (121, 132, 136, 137):
        desc111_ids = c2s(['102_313', '102_311', '102_301'])
        desc145_ids = c2s(['102_121', '102_132', '102_136', '102_137'])
##module5 inf_core_carrier_count is count for inferred core LCs
        inf_core_carrier_count = any(are_carriers(desc111_ids)) + any(are_carriers(desc145_ids))
        core_carrier_count = seq_core_carrier_count + inf_core_carrier_count

        """
        all_inf_desc = [desc111_ids, desc145_ids]
        test2 = sum(map(lambda x: any(are_carriers(x)), all_inf_desc))
        print test2 == inf_core_carrier_count
        """
        
        #observed number of sequenced non-core carriers:seq-non-core-LC = 101 + 180 
#module6 seq_noncore_ids is c2s for sequenced non-core carriers in fam102:
        seq_noncore_ids = c2s(['102_101', '102_180'])
        seq_non_core_carrier_count = sum(are_carriers(seq_noncore_ids))
        
        #Number of inferred non-core carriers:inf-non-core-LC = (136,137)
#module7 defines desc_ids for the single inferred non-core LC in fam102: 102_123
        desc123_ids = c2s(['102_136', '102_137'])
#module8 inf-non-core_carrier_count for inferred non-core carriers in fam102
        inf_non_core_carrier_count = any(are_carriers(desc123_ids))

        non_core_carrier_count = seq_non_core_carrier_count + inf_non_core_carrier_count

        identifiers = [variant.CHROM, variant.POS, variant.ID, variant.REF, variant.ALT]
        carrier_counts = [mic_ref_count, core_carrier_count, non_core_carrier_count]
        fam_gts = get_sample_gts(variant, c2s(fam_gts_H))
        row = identifiers + carrier_counts + fam_gts
        tsv_writer.writerow(row)
    
    out_file.close()


def is_carrier(gt):##Boolean function, returns True if gt = "*/0" or gt = "*/1" or gt = "*/2"
    return '/1' in gt or '/2' in gt or '/3' in gt
    
def is_ref(gt):  ##Boolean function, returns True only if gt = "0/0"
    return gt == '0/0'


def get_sample_gts(variant, ids):
    return [variant.genotype(s)['GT'] for s in ids]

def test():
    phen_path = '20190717_Musolf_phenped.txt'
    vcf_path = '190829_for_test_data_gelcc_test1.txt'#'190826_for_test_data_gelcc.vcf'
    output_path = 'core_counts_fam102.txt'
    calculate_core_count(vcf_path, output_path, phen_path)

if __name__ == "__main__":
    """
    examples:
    test:
        python 20191107_calculate_core_counts_cmd_fam102_1.py test_head10.vcf 20190717_Musolf_phenped.txt 20191108_test_counts_out.txt
    real:
        python 20191107_calculate_core_counts_cmd_fam102_1.py /home/d1381v4/20191010_AGM_Musolf_fam102/20191009_Musolf_allchr_name_AF_CADD_PHRED_fixed.vcf 20190717_Musolf_phenped.txt 20191108_fam102_out.txt
    """
    parser = argparse.ArgumentParser(description='Run core carrier analysis on fam102.')
    parser.add_argument('vcf_path', action="store",help='Path to VCF file algotithm is to be run on')
    parser.add_argument('phen_path', action="store", help='Path to phenotype file.')
    parser.add_argument('output_path', action="store", help='Path that results are written to.')
    args = parser.parse_args()

    vcf_path = args.vcf_path
    phen_path = args.phen_path
    output_path = args.output_path
    calculate_core_count(vcf_path, output_path, phen_path)
