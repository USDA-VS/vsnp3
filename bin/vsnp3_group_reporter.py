#!/usr/bin/env python

__version__ = "3.07"

import os
import sys
import pandas as pd
import vcf
import argparse
import textwrap

from vsnp3_reference_options import Ref_Options

class GroupReporter:

    def __init__(self, vcf, reference_type=None):
        if reference_type and vcf:
            reference_options = Ref_Options(reference_type)
            self.vcf = vcf
        else:
            print(f'VCF file and reference option must be provided')
            sys.exit(0)
                
        excel_path = reference_options.defining_snps
        xl = pd.ExcelFile(excel_path)
        sheet_names = xl.sheet_names
        ws = pd.read_excel(excel_path, sheet_name=sheet_names[0])
        defining_snps = ws.iloc[0]
        defsnp_iterator = iter(defining_snps.iteritems())
        next(defsnp_iterator)
        defining_snps={}
        inverted_defining_snps={}
        for abs_pos, group in defsnp_iterator:
            if '!' in abs_pos:
                inverted_defining_snps[abs_pos.replace('!', '')] = group
            else:
                defining_snps[abs_pos.replace('###', '')] = group #capture groups blocked out by ###
        self.defining_snps = defining_snps
        self.inverted_defining_snps = inverted_defining_snps

    def find_initial_positions(self, filename):
        found_positions = {}
        found_positions_mix = {}
        AC = 2
        qual_threshold = 150
        MQ = 56
        try:
            vcf_reader = vcf.Reader(open(filename, 'r'))
            try:
                for record in vcf_reader:
                    try:
                        record_qual = int(record.QUAL)
                    except TypeError:
                        record_qual = 0 
                    try:
                        # Freebayes VCFs place MQ values are placed into a list.  GATK as a float
                        record.INFO['MQ'] = record.INFO['MQ'][0]
                    except TypeError:
                        pass
                    except KeyError:
                        pass
                    chrom = record.CHROM
                    position = record.POS
                    absolute_positon = str(chrom) + ":" + str(position)
                    try:
                        if str(record.ALT[0]) != "None" and record.INFO['AC'][0] == AC and len(record.REF) == 1 and record_qual > qual_threshold and record.INFO['MQ'] > MQ:
                            found_positions.update({absolute_positon: record.REF})
                        if str(record.ALT[0]) != "None" and record.INFO['AC'][0] == 1 and len(record.REF) == 1 and record_qual > qual_threshold and record.INFO['MQ'] > MQ:
                            found_positions_mix.update({absolute_positon: record.REF})
                    except KeyError as e:
                        pass
                return filename, found_positions, found_positions_mix
            except (ZeroDivisionError, ValueError, UnboundLocalError, TypeError) as e:
                return filename, f'see error', {'': ''}, {'': ''}
        except (SyntaxError, AttributeError) as e:
            # print(type(e)(str(e) + f'\n### VCF SyntaxError {filename} File Removed'))
            os.remove(filename)
            return filename, f'see error', {'': ''}, {'': ''}
    
    def bin_and_html_table(self, filename, found_positions, found_positions_mix):
            sample_groups_list = []
            tablename = os.path.basename(filename)
            defining_snps = self.defining_snps
            inverted_defining_snps = self.inverted_defining_snps
            try:
                defining_snp = False
                for abs_position in list(defining_snps.keys() & (found_positions.keys() | found_positions_mix.keys())): #absolute positions in set union of two list
                    group = defining_snps[abs_position]
                    sample_groups_list.append(group)
                    if len(list(defining_snps.keys() & found_positions_mix.keys())) > 0:
                        tablename = f'{os.path.basename(filename)} <font color="red">[[MIXED]]</font>'
                    defining_snp = True
                if not set(inverted_defining_snps.keys()).intersection(found_positions.keys() | found_positions_mix.keys()):
                    for abs_position in list(inverted_defining_snps.keys()):
                        group = inverted_defining_snps[abs_position]
                        sample_groups_list.append(group)
                        defining_snp = True
                if defining_snp:
                    sample_groups_list = sorted(sample_groups_list)
                else:
                    sample_groups_list = ['No defining SNPs']
            except TypeError:
                message = f'File TypeError'
                print(f'{message}: {filename}')
                sample_groups_list = [f'{message}: {filename}']
                pass
            return sample_groups_list

    def get_groups(self):
        filename, found_positions, found_positions_mix = self.find_initial_positions(self.vcf)
        sample_groups_list = self.bin_and_html_table(filename, found_positions, found_positions_mix)
        return sample_groups_list

if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------
    
    Used by vsnp3 step 1 to report on group

    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-vcf', '--vcf', action='store', dest='vcf', required=True, help='Required: vcf file')
    parser.add_argument('-r', '--reference_type', action='store', dest='reference_type', required=True, default=None, help='Provide reference option.  See vsnp3_path_adder.py -s for options')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
 
    group_reporter = GroupReporter(args.vcf, args.reference_type)
    sample_groups_list = group_reporter.get_groups()
    print(sample_groups_list)
    
    print("Done")