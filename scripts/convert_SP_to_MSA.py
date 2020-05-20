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
    parser.add_argument("-c", "--context", dest="context", type=int, default=192,
                      help="set context for SBS mutation type (default: 192)")

    options = parser.parse_args()
    mutation_type = options.mutation_type
    context = options.context
    input_path = options.input_path
    signature_tables_path = options.signature_tables_path
    signatures_prefix = options.signatures_prefix

    if options.mutation_type not in ['SBS']:
        raise ValueError("Unsupported mutation type: %s. Supported type: SBS" % mutation_type)

    if options.mutation_type=='SBS':
        if context==96:
            index_col = [0,1]
            signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=index_col)
        elif context==192 or context==384:
            index_col = [0,1,2]
            signatures = pd.read_csv('%s/%s_%s_192_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=index_col)
        elif context==288:
            index_col = 0
            signatures = pd.read_csv('%s/%s_%s_288_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=index_col)
        else:
            raise ValueError("Context %i is not supported." % context)
    else:
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep=None, index_col=0)

    if options.exome:
        suffix = 'exome'
    else:
        suffix = 'all'

    # reindex 288 signatures
    if context==288:
        for element in signatures.index:
            sub = element.split('[', 1)[1].split(']')[0]
            replaced_element = element.replace('['+sub+']', sub[0])
            replaced_element = replaced_element.replace(':',':' + sub + ':')
            signatures.rename(index={element:replaced_element}, inplace=True)
        signatures.index = pd.MultiIndex.from_tuples(signatures.index.str.split(':').tolist())
        signatures.to_csv('%s/%s_%s_288_new_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), sep = ',')

    input_files = glob.glob(input_path + '/*.' + suffix)

    for file in input_files:
        input_table = pd.read_csv(file, sep=None, index_col=0)
        if context==96:
            if not '96' in file:
                continue
            input_table = input_table.sort_index(level=0)
            for element in input_table.index:
                sub = element.split('[', 1)[1].split(']')[0]
                replaced_element = element.replace('['+sub+']', sub[0])
                replaced_element = sub + ':' + replaced_element
                input_table.rename(index={element:replaced_element}, inplace=True)
            input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
            input_table = input_table.reindex(signatures.index)
        if context==192 or context==384:
            if not '384' in file:
                continue
            input_table = input_table[~input_table.index.str.contains("B:")]
            input_table = input_table[~input_table.index.str.contains("N:")]
            for element in input_table.index:
                sub = element.split('[', 1)[1].split(']')[0]
                replaced_element = element.replace('['+sub+']', sub[0])
                replaced_element = replaced_element.replace(':',':' + sub + ':')
                input_table.rename(index={element:replaced_element}, inplace=True)
            input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
            input_table = input_table.reindex(signatures.index)
        if context==288:
            if not '384' in file:
                continue
            # move bidirectional mutations to non-trancribed ones
            bidirectional_mutations = copy.deepcopy(input_table)
            bidirectional_mutations[~bidirectional_mutations.index.str.contains("B:")] = 0
            non_transcribed_mutations = copy.deepcopy(bidirectional_mutations)
            non_transcribed_mutations.index = non_transcribed_mutations.index.str.replace("N:","X:")
            non_transcribed_mutations.index = non_transcribed_mutations.index.str.replace("B:","N:")
            non_transcribed_mutations.index = non_transcribed_mutations.index.str.replace("X:","B:")
            input_table = input_table + non_transcribed_mutations
            # bidirectional_mutations_halved = round(bidirectional_mutations/2,0)
            # bidirectional_mutations_halved = bidirectional_mutations_halved.astype(int)
            input_table = input_table[~input_table.index.str.contains("B:")]
            for element in input_table.index:
                sub = element.split('[', 1)[1].split(']')[0]
                replaced_element = element.replace('['+sub+']', sub[0])
                replaced_element = replaced_element.replace(':',':' + sub + ':')
                input_table.rename(index={element:replaced_element}, inplace=True)
            input_table.index = pd.MultiIndex.from_tuples(input_table.index.str.split(':').tolist())
            input_table = input_table.reindex(signatures.index)

        compare_index(input_table, signatures)
        input_table.to_csv(file.replace(suffix, 'csv'), sep = ',')
