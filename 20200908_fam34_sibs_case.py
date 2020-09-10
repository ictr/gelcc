#!/usr/bin/env python
#20200908 fam34 siblings 3-34: case

import vcf
import pandas as pd
import csv
import argparse

def calculate_sib_count(vcf_path, output_path, phen_path):
    phen_df = pd.read_csv(phen_path, sep = '\t')
    # create a dictionary that maps CID_CIND to samples names in vcf
    cind2sample_name = dict(list(zip(phen_df['CID_CIND'], phen_df['Musolf_sampleID'])))

    vcf_reader = vcf.Reader(open(vcf_path))

    out_file = open(output_path, 'w')
    carrier_counts_H = ['MIC_REF_COUNT','total_sib_CARRIER_COUNT','NON_sib_CARRIER_COUNT']
    identifiers_H = ['CHROM', 'POS', 'ID', 'REF','ALT']
##module1  fam_gts_H is list of sequenced fam34 samples
    fam_gts_H = ['34_12','34_14','34_15','34_21','34_22','34_29','34_3','34_30','34_32','34_36','34_8']
    headers =  identifiers_H + carrier_counts_H + fam_gts_H
    tsv_writer = csv.writer(out_file, delimiter='\t', lineterminator='\n')
    tsv_writer.writerow(headers)
    for variant in vcf_reader:
        are_carriers = lambda x: list(map(is_carrier, get_sample_gts(variant, x)))
        are_refs = lambda x: list(map(is_ref, get_sample_gts(variant, x)))
        c2s = lambda cinds: [cind2sample_name[s] for s in cinds]
        
        # observed number of married-in controls (unrelated to sibs)
##module2  mic_ids is c2s(list of sequenced fam34 MICs)
        mic_ids = c2s(['34_12', '34_14','34_36'])
        mic_ref_count = sum(are_refs(mic_ids))
##module3 seq_sib_ids is c2s(list of sequenced fam34 LCS (affected sibs))
        seq_sib_ids  = c2s(['34_3','34_8','34_15'])
        seq_sib_carrier_count = sum(are_carriers(seq_sib_ids))
##module4: defines desc_ids for the 3 inferred sib LCs in fam34
        desc13_ids = c2s(['34_32','34_30'])
        desc5_ids = c2s(['34_21','34_22'])
        desc10_ids = c2s(['34_29'])
##module5 inf_sib_carrier_count is count for inferred sib LCs
        inf_sib_carrier_count = any(are_carriers(desc13_ids)) + any(are_carriers(desc5_ids)) + any(are_carriers(desc10_ids))
        sib_carrier_count = seq_sib_carrier_count + inf_sib_carrier_count

        """
        all_inf_desc = [desc111_ids, desc145_ids]
        test2 = sum(map(lambda x: any(are_carriers(x)), all_inf_desc))
        print test2 == inf_sib_carrier_count
        """
        
        #observed number of sequenced non-sib carriers:seq-non-sib-LC = 101 + 180 
#module6 seq_nonsib_ids is c2s for sequenced non-sib carriers in fam34:
        #seq_nonsib_ids = c2s(['102_101', '102_180'])
        #seq_non_sib_carrier_count = sum(are_carriers(seq_nonsib_ids))
        
        #Number of inferred non-sib carriers:inf-non-sib-LC = 0 #(136,137)
#module7 defines desc_ids for the single inferred non-sib LC in fam34: 102_123
        #desc123_ids = c2s(['102_136', '102_137'])
#module8 inf-non-sib_carrier_count for inferred non-sib carriers in fam34
        #inf_non_sib_carrier_count = any(are_carriers(desc123_ids))

        non_sib_carrier_count = 0 # seq_non_sib_carrier_count + inf_non_sib_carrier_count

        identifiers = [variant.CHROM, variant.POS, variant.ID, variant.REF, variant.ALT]
        carrier_counts = [mic_ref_count, sib_carrier_count, non_sib_carrier_count]
        fam_gts = get_sample_gts(variant, c2s(fam_gts_H))
        #row is line that gets written using 3 lists above
        row = identifiers + carrier_counts + fam_gts
        tsv_writer.writerow(row)
    
    out_file.close()


def is_carrier(gt):##Boolean function, returns True if gt = "*/0" or gt = "*/1" or gt = "*/2"
    return '/1' in gt or '/2' in gt or '/3' in gt or '/4' in gt or '/5' in gt or '/6' in gt
    
def is_ref(gt):  ##Boolean function, returns True only if gt = "0/0"
    return gt == '0/0'


def get_sample_gts(variant, ids):
    return [variant.genotype(s)['GT'] for s in ids]

def test():
    phen_path = '20190717_Musolf_phenped.txt'
    vcf_path = '190829_for_test_data_gelcc_test1.txt'#'190826_for_test_data_gelcc.vcf'
    output_path = 'sib_counts_fam34.txt'
    calculate_sib_count(vcf_path, output_path, phen_path)

if __name__ == "__main__":
    """
    examples:
    test:
        python 20191107_calculate_sib_counts_cmd_fam34_1.py test_head10.vcf 20190717_Musolf_phenped.txt 20191108_test_counts_out.txt
    real:
        python 20191107_calculate_sib_counts_cmd_fam34_1.py /home/d1381v4/20191010_AGM_Musolf_fam34/20191009_Musolf_allchr_name_AF_CADD_PHRED_fixed.vcf 20190717_Musolf_phenped.txt 20191108_fam34_out.txt
    """
    parser = argparse.ArgumentParser(description='Run sib carrier analysis on fam34.')
    parser.add_argument('vcf_path', action="store",help='Path to VCF file algotithm is to be run on')
    parser.add_argument('phen_path', action="store", help='Path to phenotype file.')
    parser.add_argument('output_path', action="store", help='Path that results are written to.')
    args = parser.parse_args()

    vcf_path = args.vcf_path
    phen_path = args.phen_path
    output_path = args.output_path
    calculate_sib_count(vcf_path, output_path, phen_path)
