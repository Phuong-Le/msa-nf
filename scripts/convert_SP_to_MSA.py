""" convert_SP_to_MSA.py
Module to convert mutation and signature tables from SigProfiler output to comma-separated, multi-indexed tables.
Reindex according to existing signature tables.
"""

import os
import glob
import copy
import pandas as pd
from argparse import ArgumentParser

def compare_index(first, second):
    if not first.index.equals(second.index):
        print('Converted index:', first.index.to_list())
        print('Target index:', second.index.to_list())
        raise ValueError("Index mismatch, check your input data.")
    return

def convert_index(input_dataframe, context=96):
    input_table = copy.deepcopy(input_dataframe)
    if context==96:
        input_table = input_table.sort_index(level=0)
        for element in input_table.index:
            sub = element.split('[', 1)[1].split(']')[0]
            replaced_element = element.replace('['+sub+']', sub[0])
            replaced_element = sub + ':' + replaced_element
            input_table.rename(index={element:replaced_element}, inplace=True)
        input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
    elif context==192 or context==384:
        if context==192:
            input_table = input_table[~input_table.index.str.contains("B:")]
            input_table = input_table[~input_table.index.str.contains("N:")]
        for element in input_table.index:
            sub = element.split('[', 1)[1].split(']')[0]
            replaced_element = element.replace('['+sub+']', sub[0])
            replaced_element = replaced_element.replace(':',':' + sub + ':')
            input_table.rename(index={element:replaced_element}, inplace=True)
        input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
    elif context==288:
        for element in input_table.index:
            sub = element.split('[', 1)[1].split(']')[0]
            replaced_element = element.replace('['+sub+']', sub[0])
            replaced_element = replaced_element.replace(':',':' + sub + ':')
            input_table.rename(index={element:replaced_element}, inplace=True)
        input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
    return input_table

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_folder", dest="input_path", default='./',
                      help="set path to input mutation tables")
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='SBS',
                      help="set mutation type (SBS, DBS, ID)")
    parser.add_argument("-e", "--exome", dest="exome", action="store_true",
                      help="Treat input matrices of the exome regions of the genome, otherwise assumme WGS")
    parser.add_argument("-s", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                      help="set path to signature tables")
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                      help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_argument("-c", "--context", dest="context", type=int, default=288,
                      help="set context (default: 288)")
    parser.add_argument("-S", "--reindex_signatures", dest="reindex_signatures", default=None,
                      help="reindex signature tables instead of mutation tables (provide full path to file)")

    options = parser.parse_args()
    mutation_type = options.mutation_type
    context = options.context
    input_path = options.input_path
    signature_tables_path = options.signature_tables_path
    signatures_prefix = options.signatures_prefix

    if mutation_type not in ['SBS', 'DBS', 'ID']:
        raise ValueError("Unsupported mutation type: %s. Supported types: SBS, DBS, ID" % mutation_type)

    signatures = None
    if mutation_type=='SBS':
        if context==96:
            index_col = [0,1]
            signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=index_col)
        elif context in [192, 288]:
            index_col = [0,1,2]
            signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type, context), sep=None, index_col=index_col)
    else:
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=0)

    if options.exome:
        suffix = '.exome'
    else:
        suffix = '.all'

    # reindex signatures
    if options.reindex_signatures:
        signature_table = pd.read_csv(options.reindex_signatures, sep=None, index_col=0)
        if mutation_type=='SBS':
            reindexed_signatures = convert_index(signature_table, context=context)
            reindexed_signatures = reindexed_signatures.reindex(signatures.index)
        else:
            # simply overwrite index for other mutation types (equality assumption)
            reindexed_signatures = signature_table
            reindexed_signatures.index = signatures.index
        reindexed_signatures.to_csv(options.reindex_signatures + '_reindexed', sep = ',')
    else:
        input_files = glob.glob(input_path + '/*' + suffix)
        for file in input_files:
            if not mutation_type in file:
                continue
            if context==192 and not '384' in file:
                continue
            if context!=192 and not str(context) in file:
                continue
            input_table = pd.read_csv(file, sep=None, index_col=0)

            if mutation_type=='SBS':
                input_table = convert_index(input_table, context=context)
                if signatures is not None:
                    input_table = input_table.reindex(signatures.index)
                    compare_index(input_table, signatures)
            else:
                # simply overwrite index for other mutation types (equality assumption)
                input_table.index = signatures.index

            new_filename = file.replace(suffix, '.csv')
            if context==192:
                new_filename = new_filename.replace('384','192')

            input_table.to_csv(new_filename, sep = ',')
