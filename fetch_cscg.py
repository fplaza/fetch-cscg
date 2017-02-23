#!/usr/bin/python
# -*- coding: utf-8 -*-

""""""

from __future__ import print_function
import argparse
import os
import subprocess
import multiprocessing
import tempfile

__author__ = "Florian Plaza Oñate"
__copyright__ = "Copyright 2017, Enterome"
__version__ = "1.0.0"
__maintainer__ = "Florian Plaza Oñate"
__email__ = "fplaza-onate@enterome.com"
__status__ = "Development"

def is_dir(path):
    """Check if path is an existing file.
    """

    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file.".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)

    return path

def is_file(path):
    """Check if path is an existing file.
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_parameters():
    """Parse command line parameters.
    """
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('-i', '--fasta-file', dest='fasta_file', type=is_file, required=True,
            help='Multi-FASTA file with protein sequences from which CSCG should be extracted')

    parser.add_argument('-d', '--cscg-db', dest='cscg_db', type=is_dir, default=os.path.join(os.getcwd(), 'cscg_db'),
            help='Path to directory that contains hmm models.')

    parser.add_argument('-m', '--hmm-model', dest='hmm_model', choices=['bacteria', 'archaea'], default='bacteria',
            help='HMM model which will be used to searh for CSCG')

    parser.add_argument('-o', '--output-file', dest='output_file', required=True,
            help='Output file listing CSCG found in input multi-FASTA file.')

    return parser.parse_args()

def check_hmm_model(cscg_db, hmm_model):
    extensions=['.h3f','.h3i','.h3m','.h3p','.tab']

    for f in (os.path.join(cscg_db, hmm_model + extension) for extension in extensions):
        if not os.path.exists(f):
            raise Exception('Error: {} not found'.format(f))

def run_hmmer(cscg_db, hmm_model, fasta_file):
    hmmsearch_results_file=tempfile.NamedTemporaryFile(delete=False).name

    with open(os.devnull, 'w') as devnull:
        subprocess.check_call(['hmmsearch',\
                '--tblout', hmmsearch_results_file,\
                '--cpu', str(multiprocessing.cpu_count()),\
                os.path.join(cscg_db, hmm_model), fasta_file],
                stdout=devnull)

    return hmmsearch_results_file

def load_pfam_cutoffs(cscg_db, hmm_model):
    pfam_cutoffs=dict()
    cutoffs_file=os.path.join(cscg_db, hmm_model + '.tab')
    with open(cutoffs_file, 'r') as istream:
        istream.next()
        for line in istream:
            pfam_name, pfam_cutoff = line.rstrip().split()
            pfam_cutoffs[pfam_name]=float(pfam_cutoff)

    return pfam_cutoffs

def filter_hmmsearch_results(pfam_cutoffs, hmmsearch_results_file):
    hmmsearch_filtered_results=list()

    with open(hmmsearch_results_file, 'r') as istream:

        for line in istream:
            if line.startswith('#'):
                continue

            line_items=line.split()
            gene_name=line_items[0]
            pfam_name=line_items[3].split('.')[0]
            score=float(line_items[5])

            if score >= pfam_cutoffs[pfam_name]:
                hmmsearch_filtered_results.append((gene_name, pfam_name, score))

    return hmmsearch_filtered_results

def write_hmmsearch_results(hmmsearch_results, output_file):
    with open(output_file, 'w') as ostream:
        print('\t'.join(['gene_name', 'pfam_name', 'score']), file=ostream)

        for gene_name, pfam_name, score in hmmsearch_results:
            print('\t'.join([gene_name, pfam_name, str(score)]), file=ostream)

def main():
    parameters = get_parameters()
    parameters.hmm_model = 'cscg_' + parameters.hmm_model

    print('Loading {} hmm model...'.format(parameters.hmm_model))
    check_hmm_model(parameters.cscg_db, parameters.hmm_model)
    pfam_cutoffs=load_pfam_cutoffs(parameters.cscg_db, parameters.hmm_model)
    print('Done. Model contains {} protein families.\n'.format(len(pfam_cutoffs)))

    print('Running hmmsearch...')
    hmmsearch_results_file=run_hmmer(parameters.cscg_db, parameters.hmm_model, parameters.fasta_file)
    print('Done.\n')

    print('Filtering results...')
    hmmsearch_filtered_results=filter_hmmsearch_results(pfam_cutoffs, hmmsearch_results_file)
    print('Done. {} cscg found\n'.format(len(hmmsearch_filtered_results)))

    print('Writing results...')
    write_hmmsearch_results(hmmsearch_filtered_results, parameters.output_file)
    print('Done.\n')

    os.remove(hmmsearch_results_file)

if __name__ == '__main__':
    main()
