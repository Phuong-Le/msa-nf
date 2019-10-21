from optparse import OptionParser
import os, copy
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
from common_methods import make_folder_if_not_exists, calculate_similarity, calculate_stat_scores

def get_comma_separated_args(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

def make_boxplot(input_table, title, xlabel, ylabel, show_mean = False, savepath = "./boxplot.pdf"):
    output_folder = savepath.rsplit('/',1)[0]
    make_folder_if_not_exists(output_folder)

    table = copy.deepcopy(input_table)
    columns = input_table.columns.to_list()

    fig = plt.figure(figsize=(0.2*len(columns) if len(columns)>30 else 6, 4))
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel(xlabel, size = 12)
    ax.set_ylabel(ylabel, size = 12)
    ax.tick_params(axis='x', rotation=90)
    ax.set_ylim(0, 1.05)
    plt.axhline(y=1, color='r', linestyle='--')
    table.boxplot(column = columns, ax=ax, grid=False, return_type='axes', sym='.',
                 meanline = show_mean, showmeans = show_mean, boxprops = dict(linewidth=3))
    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

def make_efficiency_comparison_plot(input_data, title, xlabel, ylabel, savepath = "./efficiencies.pdf"):
    output_folder = savepath.rsplit('/',1)[0]
    make_folder_if_not_exists(output_folder)

    if isinstance(input_data, pd.DataFrame):
        methods = input_data.index.to_list()
        metrics = input_data.columns.to_list()
    else:
        methods = input_data.keys()
        metrics = next(iter(input_data.values())).columns.to_list()

    f, axes = plt.subplots(1, len(metrics), sharey=True, figsize=(len(metrics)*1.3, 4))
    f.suptitle(title, fontsize=12, y=0.999)
    axes[0].set_xlabel(xlabel, size = 12)
    axes[0].set_ylabel(ylabel, size = 12)

    for metric, axis in zip(metrics, axes):
        i = 0
        for method in methods:
            i += 1
            if isinstance(input_data, pd.DataFrame):
                data_to_plot = input_data.loc[method, metric]
            else:
                data_to_plot = input_data[method][metric]
            # print(data_to_plot)
            # print(data_to_plot.mean(), data_to_plot.std())
            if isinstance(data_to_plot, float):
                mean = data_to_plot
                err = 0
            else:
                mean = data_to_plot.mean()
                err = data_to_plot.std()
            if 'unoptimised' in method:
                marker = "v"
                linestyle = '--'
            else:
                marker = "o"
                linestyle = '-'
            errorbar = axis.errorbar(i, mean, yerr=err, fmt=marker, label=method)
            errorbar[-1][0].set_linestyle(linestyle)

        axis.set_ylim(0, 1.05)
        axis.axhline(y=1, color='r', linestyle='--')
        axis.set_xlabel(metric, fontsize=10)
        axis.get_xaxis().set_ticks([])


    handles, labels = axes[0].get_legend_handles_labels()
    f.legend(handles, methods, loc='lower center', framealpha=0.0, borderaxespad=0.1, fancybox=False, shadow=False, ncol=2)

    plt.tight_layout()
    f.subplots_adjust(bottom=0.06 + 0.03*len(methods))
    plt.savefig(savepath, transparent=True)
    plt.close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", "--input_attributions_folder", dest="input_attributions_folder", default='output_tables/',
                      help="set path to NNLS output data")
    parser.add_option("-I", "--input_mutations_folder", dest="input_mutations_folder", default='input_mutation_tables/',
                      help="set path to datasets with input mutation tables")
    parser.add_option("-S", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                      help="set path to signature tables")
    parser.add_option("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                      help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_option("-d", "--dataset", dest="dataset_name", default='ESCC',
                      help="set the dataset name (e.g. ESCC, TCE)")
    parser.add_option("-o", "--output_folder", dest="output_folder", default='plots/',
                      help="set path to save plots")
    parser.add_option("-t", "--mutation_type", dest="mutation_type", default='',
                      help="set mutation type (SBS, DBS, ID)")
    parser.add_option("-c", "--context", dest="context", default=192, type='int',
                      help="set SBS context (96, 192, 1536)")
    parser.add_option("-a", "--plot_absolute_numbers", dest="abs_numbers", action="store_true",
                      help="show absolute numbers of mutations")
    parser.add_option("-s", "--signatures_to_consider", type='string', action='callback', callback=get_comma_separated_args,
                      dest = "signatures_to_consider", default = [],
                      help="set signatures for box plots as a comma-separated list, e.g. -s SBS5,SBS40")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      help="print additional information for debugging")
    parser.add_option("-n", "--number_of_b_samples", dest="number_of_b_samples", default=1000, type='int',
                      help="Number of bootstrapped samples (1000 by default)")
    parser.add_option("-N", "--number_of_samples", dest="number_of_samples", default=-1, type='int',
                      help="limit the number of samples to analyse (all by default)")


    (options, args) = parser.parse_args()

    dataset_name = options.dataset_name
    mutation_type = options.mutation_type
    context = options.context
    number_of_samples = options.number_of_samples
    number_of_b_samples = options.number_of_b_samples
    signatures_to_consider = options.signatures_to_consider
    signature_tables_path = options.signature_tables_path
    signatures_prefix = options.signatures_prefix
    input_mutations_folder = options.input_mutations_folder
    input_attributions_folder = options.input_attributions_folder + '/' + dataset_name + '/'
    output_folder = options.output_folder + '/' + dataset_name + '/' + mutation_type + '/bootstrap_plots/'
    make_folder_if_not_exists(output_folder)

    metrics = ['Cosine', 'Correlation', 'L1', 'L2', 'Chebyshev', 'Jensen-Shannon']
    scores = ['Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'F1', 'MCC']

    if mutation_type not in ['SBS','DBS','ID']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID" % mutation_type)

    if not input_attributions_folder:
        parser.error("Please specify the input path for reconstructed weights tables using -i option.")


    central_attribution_table_abs = pd.read_csv(input_attributions_folder + '/output_%s_mutations_table.csv' % mutation_type, index_col=0)
    central_attribution_table_weights = pd.read_csv(input_attributions_folder + '/output_%s_weights_table.csv' % mutation_type, index_col=0)
    bootstrap_attribution_table_abs_filename = input_attributions_folder + '/output_%s_%s_i_mutations_table.csv' % (dataset_name, mutation_type)
    bootstrap_attribution_table_weights_filename = input_attributions_folder + '/output_%s_%s_i_weights_table.csv' % (dataset_name, mutation_type)
    if 'SIM' in dataset_name and 'PCAWG' in dataset_name:
        bootstrap_attribution_table_weights_filename = input_attributions_folder + '/output_%s_i_weights_table.csv' % (mutation_type)

    if mutation_type=='SBS':
        if context==96:
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_mutations_folder, dataset_name, dataset_name, context), index_col=[0,1])
            signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=[0,1])
        elif context==192:
            # input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_mutations_folder, dataset_name, dataset_name, context), index_col=[0,1,2])
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_mutations_folder, dataset_name, dataset_name, context), index_col=0)
            signatures = pd.read_csv('%s/%s_%s_192_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=[0,1,2])
        else:
            raise ValueError("Context %i is not supported." % context)
    elif mutation_type=='DBS':
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
        input_mutations = pd.read_csv('%s/%s/WGS_%s.dinucs.csv' % (input_mutations_folder, dataset_name, dataset_name), index_col=0)
    elif mutation_type=='ID':
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
        input_mutations = pd.read_csv('%s/%s/WGS_%s.indels.csv' % (input_mutations_folder, dataset_name, dataset_name), index_col=0)

    if 'SIM' in dataset_name and not 'PCAWG' in dataset_name:
        if mutation_type=='SBS':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.%i.weights.csv' % (dataset_name, context), index_col=0)
        elif mutation_type=='DBS':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.dinucs.weights.csv' % dataset_name, index_col=0)
        elif mutation_type=='ID':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.indels.weights.csv' % dataset_name, index_col=0)
    elif 'SIM' in dataset_name:
        truth_attribution_table = pd.DataFrame(index=central_attribution_table_abs.index, columns=central_attribution_table_abs.columns, dtype=float)
        truth_attribution_table['SBS1'] = 0.065
        truth_attribution_table['SBS5'] = 0.055
        truth_attribution_table['SBS22'] = 0.230
        truth_attribution_table['SBS40'] = 0.650
        truth_attribution_table.fillna(0, inplace=True)

    # use all available signatures for boxplots unless specified
    if not signatures_to_consider:
        signatures_to_consider = list(central_attribution_table_abs.columns)

    central_attribution_table = central_attribution_table_weights
    bootstrap_attribution_table_filename = bootstrap_attribution_table_weights_filename
    colormap_label = 'Relative contribution'
    filename = 'bootstrap_plot_weights'

    # limit the number of samples to analyse (if specified by -N option)
    if options.number_of_samples!=-1:
        central_attribution_table = central_attribution_table.head(options.number_of_samples)
        if 'SIM' in dataset_name:
            truth_attribution_table = truth_attribution_table.head(options.number_of_samples)

    main_title = dataset_name.replace('_','/') + ' data, ' + mutation_type + ' mutation type'

    # samples to consider
    samples = central_attribution_table.index.to_list()

    # initialise attribution dictionary
    signature_attributions_dict = {}
    for sample in samples:
        signature_attributions_dict[sample] = pd.DataFrame(index=range(number_of_b_samples), columns=signatures_to_consider, dtype=float)

    # mutation categories from signatures table
    categories = signatures.index.to_list()
    # fill bootstrap dataframes
    for i in range(number_of_b_samples):
        bootstrap_attribution_table = pd.read_csv(bootstrap_attribution_table_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        if options.number_of_samples!=-1:
            bootstrap_attribution_table = bootstrap_attribution_table.head(options.number_of_samples)
        # bootstrap_attribution_table_abs = pd.read_csv(bootstrap_attribution_table_abs_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        for sample in samples:
            for signature in signatures_to_consider:
                signature_attributions_dict[sample].loc[i, signature] = bootstrap_attribution_table.loc[sample, signature]


    if 'SIM' in dataset_name:
        acting_signatures = []
        for signature in list(truth_attribution_table.columns):
            if truth_attribution_table[signature].max()>0:
                acting_signatures.append(signature)

    stat_scores_tables = pd.DataFrame(index = central_attribution_table.index, columns = scores, dtype=float)
    signatures_prevalences = pd.DataFrame(index = central_attribution_table.index, columns = signatures_to_consider, dtype=float)
    signatures_CPs = pd.DataFrame(index = central_attribution_table.index, columns = acting_signatures, dtype=float)
    for sample in samples:
        for signature in signatures_to_consider:
            array = signature_attributions_dict[sample][signature]
            confidence_interval = [np.percentile(array, 2.5), np.percentile(array, 97.5)]
            # TODO: create an uncertainty dataframe and output as JSON

            if 'SIM' in dataset_name:
                if signature in acting_signatures:
                    if confidence_interval[0]<=truth_attribution_table.loc[sample, signature]<=confidence_interval[1]:
                        signatures_CPs.loc[sample,signature] = 1
                    else:
                        signatures_CPs.loc[sample,signature] = 0
                else:
                    if confidence_interval[0]<=truth_attribution_table.loc[sample, signature]<=confidence_interval[1]:
                        signatures_CPs.loc[sample,'Others'] = 1
                    else:
                        signatures_CPs.loc[sample,'Others'] = 0

            # disregard attributions consistent with 0
            if confidence_interval[0]==0:
                signature_attributions_dict[sample][signature] = 0

    for sample in samples:
        signatures_prevalences.loc[sample, :] = signature_attributions_dict[sample].astype(bool).sum(axis=0)
        if 'SIM' in dataset_name:
            truth_one_sample_table = truth_attribution_table.loc[sample, :].to_frame().T
            truth_table = truth_one_sample_table.loc[truth_one_sample_table.index.repeat(number_of_b_samples)]
            print(signatures)
            print(signature_attributions_dict[sample])
            print(truth_table)
            stat_scores_tables.loc[sample, :] = calculate_stat_scores(signatures, signature_attributions_dict[sample], truth_table)

    signatures_prevalences = signatures_prevalences/number_of_b_samples
    make_boxplot(signatures_prevalences, 'Signature prevalence (%s)' % mutation_type, 'Signatures', 'Reconstructed (%)', savepath = output_folder + '/signature_prevalences.pdf')

    if 'SIM' in dataset_name:
        signatures_CPs = signatures_CPs.sum(axis=0).to_frame().T.div(len(samples))
        make_boxplot(signatures_CPs, 'Confidence probability (%s)' % mutation_type, 'Signatures', 'CP (%)', savepath = output_folder + '/signature_confidence_probabilities.pdf')
        make_boxplot(stat_scores_tables, 'Metrics (%s)' % mutation_type, 'Metrics', 'Values', savepath = output_folder + '/sensitivity_specificity_metrics.pdf')
