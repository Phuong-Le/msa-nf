from argparse import ArgumentParser
import copy
import numpy as np
import pandas as pd
from common_methods import make_folder_if_not_exists, write_data_to_JSON, calculate_confidence_interval

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_attributions_folder", dest="input_attributions_folder", default='output_tables/',
                        help="set path to NNLS output data")
    parser.add_argument("-S", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                        help="set path to signature tables")
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                        help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='SIM',
                        help="set the dataset name (e.g. SIM)")
    parser.add_argument("-o", "--output_folder", dest="output_folder", default='output_tables/',
                        help="set path to save plots")
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='',
                        help="set mutation type (SBS, DBS, ID)")
    parser.add_argument("-c", "--context", dest="context", default=192, type=int,
                        help="set SBS context (96, 192)")
    parser.add_argument("-a", "--plot_absolute_numbers", dest="abs_numbers", action="store_true",
                        help="show absolute numbers of mutations")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="print additional information for debugging")
    parser.add_argument("-n", "--number_of_b_samples", dest="number_of_b_samples", default=1000, type=int,
                        help="Number of bootstrapped samples (1000 by default)")
    parser.add_argument("-N", "--number_of_samples", dest="number_of_samples", default=-1, type=int,
                        help="limit the number of samples to analyse (all by default)")

    args = parser.parse_args()

    dataset_name = args.dataset_name
    mutation_type = args.mutation_type
    context = args.context
    number_of_b_samples = args.number_of_b_samples
    signature_tables_path = args.signature_tables_path
    signatures_prefix = args.signatures_prefix
    input_attributions_folder = args.input_attributions_folder + '/' + dataset_name + '/'
    output_folder = args.output_folder + '/' + dataset_name
    make_folder_if_not_exists(output_folder)

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (DBS, ID).")
    elif mutation_type not in ['SBS', 'DBS', 'ID']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID" % mutation_type)

    print("*"*50)
    print("Making plots for %s mutation type, %s dataset" % (mutation_type, dataset_name))

    central_attribution_table_abs = pd.read_csv(input_attributions_folder + '/output_%s_%s_mutations_table.csv' % (dataset_name, mutation_type), index_col=0)
    central_attribution_table_weights = pd.read_csv(input_attributions_folder + '/output_%s_%s_weights_table.csv' % (dataset_name, mutation_type), index_col=0)
    central_stat_table = pd.read_csv(input_attributions_folder + '/output_%s_%s_stat_info.csv' % (dataset_name, mutation_type), index_col=0)
    bootstrap_attribution_table_abs_filename = input_attributions_folder + '/bootstrap_output/output_%s_%s_i_mutations_table.csv' % (dataset_name, mutation_type)
    bootstrap_attribution_table_weights_filename = input_attributions_folder + '/bootstrap_output/output_%s_%s_i_weights_table.csv' % (dataset_name, mutation_type)
    bootstrap_stat_table_filename = input_attributions_folder + '/bootstrap_output/output_%s_%s_i_stat_info.csv' % (dataset_name, mutation_type)

    if mutation_type == 'SBS':
        if context == 96:
            signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=[0, 1])
        elif context in [192, 288]:
            signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type, context), index_col=[0, 1, 2])
        else:
            raise ValueError("Context %i is not supported." % context)
    elif mutation_type == 'DBS':
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
    elif mutation_type == 'ID':
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)

    signatures_to_consider = list(central_attribution_table_abs.columns)
    signatures = signatures[signatures_to_consider]

    if args.abs_numbers:
        central_attribution_table = central_attribution_table_abs
        bootstrap_attribution_table_filename = bootstrap_attribution_table_abs_filename
        filename = mutation_type + '_bootstrap_output_abs_mutations'
    else:
        central_attribution_table = central_attribution_table_weights
        bootstrap_attribution_table_filename = bootstrap_attribution_table_weights_filename
        filename = mutation_type + '_bootstrap_output_weights'

    # limit the number of samples to analyse (if specified by -N option)
    if args.number_of_samples != -1:
        central_attribution_table = central_attribution_table.head(args.number_of_samples)

    main_title = dataset_name.replace('_', '/') + ' data, ' + mutation_type + ' mutation type'

    # samples and metrics to consider
    samples = central_attribution_table.index.to_list()
    stat_metrics = central_stat_table.columns.to_list()

    # initialise statistical metrics to fill with bootstrap values
    stat_metrics_dict = {}
    for metric in stat_metrics:
        stat_metrics_dict[metric] = pd.DataFrame(index=range(number_of_b_samples), columns=samples, dtype=float)

    # initialise per-sample attribution dictionary
    attributions_per_sample_dict = {}
    for sample in samples:
        attributions_per_sample_dict[sample] = pd.DataFrame(index=range(number_of_b_samples), columns=signatures_to_consider, dtype=float)

    # initialise per-signature attribution dictionary
    attributions_per_signature_dict = {}
    for signature in signatures_to_consider:
        attributions_per_signature_dict[signature] = pd.DataFrame(index=range(number_of_b_samples), columns=samples, dtype=float)

    # mutation categories from signatures table
    categories = signatures.index.to_list()
    # initialise mutation spectra dictionary
    mutation_spectra_dict = {}
    for sample in samples:
        mutation_spectra_dict[sample] = pd.DataFrame(index=range(number_of_b_samples), columns=categories, dtype=float)

    # fill bootstrap dataframes
    for i in range(number_of_b_samples):
        bootstrap_attribution_table = pd.read_csv(bootstrap_attribution_table_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        bootstrap_attribution_table_abs = pd.read_csv(bootstrap_attribution_table_abs_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        stat_table = pd.read_csv(bootstrap_stat_table_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        # back-calculate mutation spectra from bootstrap attributions and signatures
        mutation_spectra = bootstrap_attribution_table_abs.dot(signatures.T)
        for sample in samples:
            for signature in signatures_to_consider:
                attributions_per_sample_dict[sample].loc[i, signature] = bootstrap_attribution_table.loc[sample, signature]
                attributions_per_signature_dict[signature].loc[i, sample] = bootstrap_attribution_table.loc[sample, signature]
            for category in categories:
                mutation_spectra_dict[sample].loc[i, category] = mutation_spectra.loc[sample, category]
            for metric in stat_metrics:
                stat_metrics_dict[metric].loc[i, sample] = stat_table.loc[sample, metric]
                stat_metrics_dict[metric].loc[i, sample] = stat_table.loc[sample, metric]

    confidence_intervals = pd.DataFrame(index=samples, columns=signatures_to_consider, dtype=object)
    for sample in samples:
        for signature in signatures_to_consider:
            confidence_interval = calculate_confidence_interval(attributions_per_sample_dict[sample][signature])
            central_value = central_attribution_table.loc[sample, signature]
            confidence_intervals.loc[sample, signature] = [central_value, confidence_interval]

    confidence_intervals.to_csv(output_folder + '/CIs_' + filename + '.csv')

    write_data_to_JSON(attributions_per_sample_dict, output_folder + '/attributions_per_sample_' + filename + '.json')
    write_data_to_JSON(attributions_per_signature_dict, output_folder + '/attributions_per_signature_' + filename + '.json')
    write_data_to_JSON(stat_metrics_dict, output_folder + '/stat_metrics_' + filename + '.json')
    # # pd.to_json function does not fully support MultiIndex (present in signature tables), hence commented out for now
    # write_data_to_JSON(mutation_spectra_dict, output_folder + '/mutation_spectra_' + filename + '.json')
