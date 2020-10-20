""" convert_SP_to_MSA.py
Module to convert mutation and signature tables from SigProfiler output to
comma-separated, multi-indexed tables. Output files are saved in MSA directories
specifiable by -o flag for mutation tables and -s for signatures.
Reindexing according to existing signature tables specified using -s and -p flags.
"""

import os
import glob
import copy
import pandas as pd
from argparse import ArgumentParser
from common_methods import make_folder_if_not_exists

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
    parser = ArgumentParser(description='Convert SigProfilerMatrixGenerator output to MSA format')
    parser.add_argument("-i", "--input_folder", dest="input_path",
                      help="set path to SigProfiler output with mutation tables")
    parser.add_argument("-d", "--dataset_name", dest="dataset_name", default='',
                        help="set dataset name to use in converted filenames")
    parser.add_argument("-t", "--mutation_types", nargs='+', dest="mutation_types", default=['SBS','DBS','ID'],
                      help="set mutation types, e.g. -t SBS DBS ID (default)")
    parser.add_argument("-c", "--contexts", nargs='+', dest="contexts", type=int, default=[96, 288],
                      help="set SBS contexts e.g. -c 96 288 (default). Supported contexts: 96, 192, 288")
    parser.add_argument("-o", "--output_path", dest="output_path", default='input_mutation_tables',
                        help="set output path for converted mutation tables (default: input_mutation_tables)")
    # parser.add_argument("-e", "--exome", dest="exome", action="store_true",
    #                   help="Treat input matrices of the exome regions of the genome, otherwise assumme WGS")
    parser.add_argument("-s", "--signature_path", dest="signature_tables_path", default='signature_tables',
                      help="set path to signature tables to extract indexes (default: signature_tables)")
    parser.add_argument("-p", "--input_signatures_prefix", dest="input_signatures_prefix", default='sigProfiler',
                      help="set prefix in signature filenames to extract indexes (sigProfiler by default)")
    parser.add_argument("-n", "--reindexed_signatures_prefix", dest="reindexed_signatures_prefix", default='sigProfilerNew',
                      help="set prefix for reindexed signatures filenames (sigProfilerNew by default)")
    parser.add_argument("-S", "--reindex_signatures", dest="reindex_signatures", default=None,
                      help="reindex signature tables instead of mutation tables (provide full path to file)")

    options = parser.parse_args()
    input_path = options.input_path
    dataset_name = options.dataset_name
    mutation_types = options.mutation_types
    contexts = options.contexts
    output_path = options.output_path
    signature_tables_path = options.signature_tables_path
    input_signatures_prefix = options.input_signatures_prefix
    reindexed_signatures_prefix = options.reindexed_signatures_prefix

    if not options.reindex_signatures:
        if not input_path:
            raise ValueError("Please specify the input path to SigProfiler tables using -i flag. Alternatively, use -S flag to reindex a signature table (see -h for help)")
        if not dataset_name:
            raise ValueError("Please specify an arbitrary dataset name for use in MSA pipeline execution (see -h for help)")
        print('Converting mutation tables for dataset', dataset_name, ' mutation types', mutation_types, ', considering SBS contexts', contexts)
        print('SigProfilerMatrixGenerator output to parse:', input_path)
    else:
        print('Converting signature table', options.reindex_signatures)

    for mutation_type in mutation_types:
        if mutation_type not in ['SBS', 'DBS', 'ID']:
            raise ValueError("Unsupported mutation type: %s. Supported types: SBS, DBS, ID" % mutation_type)

    signatures = {}
    for mutation_type in mutation_types:
        if mutation_type=='SBS':
            for context in contexts:
                if context==96:
                    index_col = [0,1]
                    signatures[mutation_type + str(context)] = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, input_signatures_prefix, mutation_type), sep=',', index_col=index_col)
                elif context in [192, 288]:
                    index_col = [0,1,2]
                    signatures[mutation_type + str(context)] = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, input_signatures_prefix, mutation_type, context), sep=',', index_col=index_col)
        else:
            signatures[mutation_type] = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, input_signatures_prefix, mutation_type), sep=',', index_col=0)

    if options.reindex_signatures:
        # reindex signatures
        signature_table_to_reindex = pd.read_csv(options.reindex_signatures, sep='\t', index_col=0)
        mutation_type = mutation_types[0]
        context = contexts[0]
        if mutation_type=='DBS':
            context = 78
        elif mutation_type=='ID':
            context = 83
        print('Assuming %s mutation type' % mutation_type, '(%i context)' % context)
        if mutation_type=='SBS':
            reindexed_signatures = convert_index(signature_table_to_reindex, context=context)
            template_signatures = signatures[mutation_type + str(context)]
            reindexed_signatures = reindexed_signatures.reindex(template_signatures.index)
            compare_index(reindexed_signatures, template_signatures)
        else:
            # simply overwrite index for other mutation types (equality assumption)
            reindexed_signatures = signature_table_to_reindex
            template_signatures = signatures[mutation_type]
            reindexed_signatures.index = template_signatures.index
            compare_index(reindexed_signatures, template_signatures)
        reindexed_signatures_filename = '%s/%s_%s_%i_signatures.csv' % (signature_tables_path, reindexed_signatures_prefix, mutation_type, context)
        if mutation_type!='SBS' or (mutation_type=='SBS' and context==96):
            reindexed_signatures_filename = reindexed_signatures_filename.replace('_%i' % context,'')
        reindexed_signatures.to_csv(reindexed_signatures_filename, sep = ',')
        print('Done. Check the output signature table:', reindexed_signatures_filename)
    else:
        # reindex mutation tables
        make_folder_if_not_exists(output_path + '/' + dataset_name)
        for mutation_type in mutation_types:
            input_files = glob.glob(input_path + '/%s/*%s*' % (mutation_type, mutation_type))
            for file in input_files:
                if mutation_type=='SBS':
                    for context in contexts:
                        if context==192 and not '384' in file:
                            continue
                        if context!=192 and not str(context) in file:
                            continue
                        input_table = pd.read_csv(file, sep='\t', index_col=0)
                        print('Converting:', mutation_type, context, file)
                        input_table = convert_index(input_table, context=context)
                        signature_table = signatures[mutation_type + str(context)]
                        if signature_table is not None:
                            input_table = input_table.reindex(signature_table.index)
                            compare_index(input_table, signature_table)
                        new_filename = output_path + '/%s/WGS_%s.%i.csv' % (dataset_name, dataset_name, context)
                        if context==192:
                            new_filename = new_filename.replace('384','192')
                        input_table.to_csv(new_filename, sep = ',')
                else:
                    if mutation_type=='DBS' and not '78' in file:
                        continue
                    if mutation_type=='ID' and not '83' in file:
                        continue
                    # simply overwrite index for other mutation types (equality assumption)
                    print('Converting:', mutation_type, file)
                    input_table = pd.read_csv(file, sep='\t', index_col=0)
                    input_table.index = signatures[mutation_type].index
                    if mutation_type == 'DBS':
                        new_filename = output_path + '/%s/WGS_%s.dinucs.csv' % (dataset_name, dataset_name)
                    elif mutation_type == 'ID':
                        new_filename = output_path + '/%s/WGS_%s.indels.csv' % (dataset_name, dataset_name)
                    input_table.to_csv(new_filename, sep = ',')
        print('Done, please check the outputs:', output_path + '/' + dataset_name)
