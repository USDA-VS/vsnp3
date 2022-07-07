#!/usr/bin/env python

__version__ = "3.06"

import os
import gzip
import re
import glob
import shutil
import regex
import argparse
import textwrap
from collections import OrderedDict
import multiprocessing
multiprocessing.set_start_method('spawn', True)
from dask import delayed
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import subprocess

from vsnp3_file_setup import Setup
from vsnp3_file_setup import bcolors
from vsnp3_file_setup import Banner
from vsnp3_file_setup import Latex_Report
from vsnp3_file_setup import Excel_Stats


class Spoligo(Setup):

    def __init__(self, FASTQ_R1=None, FASTQ_R2=None, debug=False):
        Setup.__init__(self, FASTQ_R1=FASTQ_R1, FASTQ_R2=FASTQ_R2, debug=debug)
        self.print_run_time('Spoligotype')
        self.cpu_count_half = int(self.cpus / 2)
        real_path = os.path.dirname(os.path.realpath(__file__))
        self.spoligo_db = real_path + "/../dependencies/spoligotype_db.txt"
        self.spoligo_fasta = real_path + "/../dependencies/spoligo_spacers.fasta"
    
    def finding_sp(self, spacer_sequence):
        # spacer_id, spacer_sequence = spacer_id_and_spacer_sequence
        total_count = 0
        total_finds = 0

        total_finds = [len(regex.findall("(" + spacer + "){s<=1}", self.seq_string)) for spacer in spacer_sequence]

        for number in total_finds:
            total_count += number

        return (total_count)

    def count_spacer_occurence(self):
        otal_count = 0
        total_finds = 0

        if len(self.FASTQ_list) == 1:
            cmd = ['bbduk.sh',
                  'in={}'.format(self.FASTQ_list[0]),
                  'ref={}'.format(self.spoligo_fasta),
                  'k=21', 'rcomp=t',
                  'edist=1',  # up to 1 missmatch
                  'maskmiddle=f', # Do not treat the middle base of a kmer as a wildcard
                  'stats=stats.txt',
                  'ow=t']
        else:
            cmd = ['bbduk.sh',
                  'in={}'.format(self.FASTQ_list[0]),
                  'in2={}'.format(self.FASTQ_list[1]),
                  'ref={}'.format(self.spoligo_fasta),
                  'k=21', 'rcomp=t',
                  'edist=1',  # up to 1 missmatch
                  'maskmiddle=f', # Do not treat the middle base of a kmer as a wildcard
                  'stats=stats.txt',
                  'ow=t']

        subprocess.run(cmd)

        # Parse stats file
        """
        #File /home/bioinfo/analyses/vsnp3_test_spoligo/SRR16058435_R1.fastq.gz
        #Total  822714
        #Matched    799 0.09712%
        #Name   Reads   ReadsPct
        spacer25    62  0.00754%
        spacer02    48  0.00583%
        """

        count_summary = dict()
        with open('stats.txt', 'r') as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue

                if line.startswith('#'):
                    continue
                else:
                    field_list = line.split('\t')
                    count_summary[field_list[0]] = int(field_list[1])

        # Fill any spacers with zero counts
        # Parse spacer fasta file
        spacer_list = list()
        with open(self.spoligo_fasta, 'r') as f:
            for line in f:
                line = line.rstrip()
                if not line:
                    continue
                if line.startswith('>'):
                    spacer_list.append(line[1:])  # Drop the leading ">"

        for spacer_id in spacer_list:
            if spacer_id not in count_summary:
                count_summary[spacer_id] = 0

        return count_summary

    def binary_to_octal(self, binary):
        #binary_len = len(binary)
        i = 0
        ie = 1
        octal = ""
        while ie < 43:
            ie = i + 3
            # print(binary[i:ie])
            region = binary[i:ie]
            region_len = len(region)
            i += 3
            if int(region[0]) == 1:
                if region_len < 2: # for the lone spacer 43.  When present needs to be 1 not 4.
                    oct = 1
                else:
                    oct = 4
            else:
                oct = 0
            try:
                if int(region[1]) == 1:
                    oct += 2
                if int(region[2]) == 1:
                    oct += 1
            except IndexError:
                pass
            octal = octal + str(oct)
        return(octal)

    def spoligo(self):

        octal = None
        sbcode = None
        db_binarycode = None
        sample_binary = None

        count_summary = self.count_spacer_occurence()
        count_summary = OrderedDict(sorted(count_summary.items()))

        spoligo_binary_dictionary = {}
        self.call_cut_off = 4
        for k, v in count_summary.items():
            if v > self.call_cut_off:
                spoligo_binary_dictionary.update({k: 1})
            else:
                spoligo_binary_dictionary.update({k: 0})
        spoligo_binary_dictionary = OrderedDict(sorted(spoligo_binary_dictionary.items()))
        spoligo_binary_list = []
        for v in spoligo_binary_dictionary.values():
            spoligo_binary_list.append(v)
        sample_binary = ''.join(str(e) for e in spoligo_binary_list)  #sample_binary correct
        self.sample_binary = sample_binary
        self.octal = self.binary_to_octal(sample_binary)
        found = False
        with open(self.spoligo_db) as spoligo_db_file: # put into dictionary or list
            for line in spoligo_db_file:
                line = line.rstrip()
                sbcode = line.split()[1]
                db_binarycode = line.split()[2]
                if sample_binary == db_binarycode:
                    found = True
                    self.sbcode = sbcode
        if not found:
            if sample_binary == '0000000000000000000000000000000000000000000':
                self.sbcode = "spoligo not found, binary all zeros, see spoligo file"
            else:
                self.sbcode = "Not Found"
        self.sample_binary = sample_binary
        self.count_summary_list=[]
        for spacer, count in count_summary.items():
            self.count_summary_list.append(count)

    def latex(self, tex):
        blast_banner = Banner("Spoligotype")
        print(r'\begin{table}[ht!]', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{center}', file=tex)
        print('\includegraphics[scale=1]{' + blast_banner.banner + '}', file=tex)
        print(r'\end{center}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\begin{adjustbox}{width=1\textwidth}', file=tex)
        print(r'\begin{tabular}{ l | l | l }', file=tex)
        print(r'\multicolumn{3}{l}{Spacer Counts} \\', file=tex)
        print(r'\hline', file=tex) 
        count_summary = ":".join(map(str, self.count_summary_list))
        print(r'\multicolumn{3}{l}{' + f'{count_summary}' + r' } \\', file=tex)
        print(r'\hline', file=tex)
        print(f'Binary Code, threshold greater than {str(self.call_cut_off)} spacer counts & Octal Code & SB Number \\\\', file=tex)
        print(r'\hline', file=tex)
        print(f'{self.sample_binary} & {self.octal} & {self.sbcode} \\\\', file=tex)
        print(r'\hline', file=tex)
        print(r'\end{tabular}', file=tex)
        print(r'\end{adjustbox}', file=tex)
        print(r'\end{table}', file=tex)
    
    def excel(self, excel_dict):
        excel_dict['Spoligotype Spacer Counts'] = f'{":".join(map(str, self.count_summary_list))}'
        excel_dict['Spoligotype Binary Code'] = f'binary-{self.sample_binary}'
        excel_dict['Spoligotype Octal Code'] = f'octal-{self.octal}'
        excel_dict['Spoligotype SB Number'] = f'{self.sbcode}'


if __name__ == "__main__": # execute if directly access by the interpreter

    parser = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent('''\

    ---------------------------------------------------------

    Mycobacterium bovis spoligotype from WGS
   
    '''), epilog='''---------------------------------------------------------''')
    
    parser.add_argument('-r1', '--FASTQ_R1', action='store', dest='FASTQ_R1', required=True, help='Required: single read')
    parser.add_argument('-r2', '--FASTQ_R2', action='store', dest='FASTQ_R2', required=False, default=None, help='Optional: paired read')
    parser.add_argument('-d', '--debug', action='store_true', dest='debug', default=False, help='turn off map.pooling of samples')
    parser.add_argument('-v', '--version', action='version', version=f'{os.path.abspath(__file__)}: version {__version__}')

    args = parser.parse_args()
    spoligo = Spoligo(FASTQ_R1=args.FASTQ_R1, FASTQ_R2=args.FASTQ_R2, debug=args.debug)
    spoligo.spoligo()

    #Latex report
    latex_report = Latex_Report(spoligo.sample_name)
    spoligo.latex(latex_report.tex)
    latex_report.latex_ending()

    #Excel Stats
    excel_stats = Excel_Stats(spoligo.sample_name)
    spoligo.excel(excel_stats.excel_dict)
    excel_stats.post_excel()

    temp_dir = './temp'
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)
    files_grab = []
    for files in ('*.aux', '*.log', '*tex', '*png', '*out'):
        files_grab.extend(glob.glob(files))
    for each in files_grab:
        shutil.move(each, temp_dir)

    if args.debug is False:
        shutil.rmtree(temp_dir)

# Updated July 2022 by Marc-Olivier Duceppe
