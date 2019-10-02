from optparse import OptionParser
import os, copy
import warnings
import numpy as np
import pandas as pd
from common_methods import get_comma_separated_args, make_folder_if_not_exists
import matplotlib.pyplot as plt
plt.rc('text', usetex=False)
import matplotlib.patches as mpatches
import seaborn as sns

def plot_bootstrap_attribution_histograms(input_table, centrals, title, savepath='./hist.pdf'):
    table = copy.deepcopy(input_table)
    output_folder = savepath.rsplit('/',1)[0]
    make_folder_if_not_exists(output_folder)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    if max(centrals.values())<=1:
        bins = [0.01*k for k in range(101)]
    else:
        bins = None

    for column in table.columns:
        color = next(ax._get_lines.prop_cycler)['color']
        plt.hist(table[column].to_list(), label = column, bins = bins, histtype='step', color=color, density=True)
        plt.axvline(x=centrals[column], color=color, ls = '--')

    legend = plt.legend()
    plt.gca().add_artist(legend)

    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel('Attribution', fontsize=12)
    ax.set_ylabel('Number of samples (normalised)', fontsize=12)
    if max(centrals.values())<=1:
        ax.set_xticks([0.1*k for k in range(11)])

    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

def make_bootstrap_attribution_boxplot(input_table, centrals = None, truth = None, prefix = 'Bootstrap', method = 'NNLS', title='', x_label = '', y_label = '', savepath='./boxplot.pdf'):
    table = copy.deepcopy(input_table)
    output_folder = savepath.rsplit('/',1)[0]
    make_folder_if_not_exists(output_folder)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_xlabel(x_label, size = 12)
    ax.set_ylabel(y_label, size = 12)
    ax.set_title(title, fontsize=14, pad=10)

    sns.set(style="whitegrid")
    ax = sns.boxplot(data=table, palette="colorblind", whis=[2.5, 97.5], showfliers=False)

    # add counts to labels
    # sig_labels = [t.get_text() for t in ax.get_xticklabels()]
    # sig_labels = [sig + ' (N = %i)' % len(table[sig].dropna()) for sig in sig_labels]
    # ax.set_xticklabels(sig_labels)
    ax.tick_params(axis='x', rotation=90, labelsize=8)
    column_range = range(len(table.columns))

    if centrals or truth:
        # for id, name in zip(column_range, table.columns.to_list()):
        for id in column_range:
            if centrals:
                # plt.plot([column_range[id]], [centrals[name]], marker="o", linestyle="None", markersize=3, color='red', label = method)
                plt.plot([column_range[id]], [list(centrals.values())[id]], marker="o", linestyle="None", markersize=3, color='red', label = method)
            if truth:
                # plt.plot([column_range[id]], [truth[name]], marker="o", linestyle="None", markersize=3, color='green', label = "Truth")
                plt.plot([column_range[id]], [list(truth.values())[id]], marker="o", linestyle="None", markersize=3, color='green', label = "Truth")

    # Add a boxplot patch for legend
    boxplot_patch = mpatches.Patch(color='gray', edgecolor='black', lw=2, label=prefix + ' ' + method)
    handles, labels = ax.get_legend_handles_labels()
    if truth and centrals:
        handles = handles[:2]
    else:
        handles = handles[:1]

    handles.append(boxplot_patch)
    plt.legend(handles=handles)

    # adjust figure size for long data
    if len(column_range)>100:
        fig = ax.get_figure()
        fig.set_size_inches(20, 5)

    # if relative:
        # ax.set_ylim(0, 1)

    # Add a horizontal grid to the plot
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    plt.tight_layout()
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
    parser.add_option("-d", "--dataset", dest="dataset_name", default='SIM',
                      help="set the dataset name (e.g. SIM)")
    parser.add_option("-o", "--output_folder", dest="output_folder", default='plots/',
                      help="set path to save plots")
    parser.add_option("-t", "--mutation_type", dest="mutation_type", default='',
                      help="set mutation type (SBS, DBS, ID)")
    parser.add_option("-c", "--context", dest="context", default=192, type='int',
                      help="set SBS context (96, 192)")
    parser.add_option("-a", "--plot_absolute_numbers", dest="abs_numbers", action="store_true",
                      help="show absolute numbers of mutations")
    parser.add_option("-s", "--signatures_for_boxplots", type='string', action='callback', callback=get_comma_separated_args,
                      dest = "signatures_for_boxplots", default = [],
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
    number_of_b_samples = options.number_of_b_samples
    signatures_for_boxplots = options.signatures_for_boxplots
    signature_tables_path = options.signature_tables_path
    signatures_prefix = options.signatures_prefix
    input_mutations_folder = options.input_mutations_folder
    input_attributions_folder = options.input_attributions_folder + '/' + dataset_name + '/'
    output_folder = options.output_folder + '/' + dataset_name + '/' + mutation_type + '/bootstrap_plots/'
    make_folder_if_not_exists(output_folder)

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (DBS, ID).")
    elif mutation_type not in ['SBS','DBS','ID']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID" % mutation_type)

    print("*"*50)
    print("Making plots for %s mutation type, %s dataset" % (mutation_type, dataset_name))

    central_attribution_table_abs = pd.read_csv(input_attributions_folder + '/output_%s_mutations_table.csv' % mutation_type, index_col=0)
    central_attribution_table_weights = pd.read_csv(input_attributions_folder + '/output_%s_weights_table.csv' % mutation_type, index_col=0)
    central_stat_table = pd.read_csv(input_attributions_folder + '/output_%s_stat_info.csv' % mutation_type, index_col=0)
    bootstrap_attribution_table_abs_filename = input_attributions_folder + '/output_%s_%s_i_mutations_table.csv' % (dataset_name, mutation_type)
    bootstrap_attribution_table_weights_filename = input_attributions_folder + '/output_%s_%s_i_weights_table.csv' % (dataset_name, mutation_type)
    bootstrap_stat_table_filename = input_attributions_folder + '/output_%s_%s_i_stat_info.csv' % (dataset_name, mutation_type)

    if mutation_type=='SBS':
        if context==96:
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_mutations_folder, dataset_name, dataset_name, context), index_col=[0,1])
            signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=[0,1])
        elif context==192:
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_mutations_folder, dataset_name, dataset_name, context), index_col=[0,1,2])
            signatures = pd.read_csv('%s/%s_%s_192_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=[0,1,2])
        else:
            raise ValueError("Context %i is not supported." % context)
    elif mutation_type=='DBS':
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
        input_mutations = pd.read_csv('%s/%s/WGS_%s.dinucs.csv' % (input_mutations_folder, dataset_name, dataset_name), index_col=0)
    elif mutation_type=='ID':
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
        input_mutations = pd.read_csv('%s/%s/WGS_%s.indels.csv' % (input_mutations_folder, dataset_name, dataset_name), index_col=0)

    if 'SIM' in dataset_name:
        if mutation_type=='SBS':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.%i.weights.csv' % (dataset_name, context), index_col=0)
        elif mutation_type=='DBS':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.dinucs.weights.csv' % dataset_name, index_col=0)
        elif mutation_type=='ID':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.indels.weights.csv' % dataset_name, index_col=0)

    # use all available signatures for boxplots unless specified
    if not signatures_for_boxplots:
        signatures_for_boxplots = list(central_attribution_table_abs.columns)

    if options.abs_numbers:
        central_attribution_table = central_attribution_table_abs
        bootstrap_attribution_table_filename = bootstrap_attribution_table_abs_filename
        colormap_label = 'Absolute mutations number'
        filename = 'bootstrap_plot_abs_mutations'
    else:
        central_attribution_table = central_attribution_table_weights
        bootstrap_attribution_table_filename = bootstrap_attribution_table_weights_filename
        colormap_label = 'Relative contribution'
        filename = 'bootstrap_plot_weights'

    # limit the number of samples to analyse (if specified by -N option)
    if options.number_of_samples!=-1:
        central_attribution_table = central_attribution_table.head(options.number_of_samples)

    main_title = dataset_name.replace('_','/') + ' data, ' + mutation_type + ' mutation type'

    # samples and metrics to consider
    samples = central_attribution_table.index.to_list()
    stat_metrics = central_stat_table.columns.to_list()

    # initialise statistical metrics to fill with bootstrap values
    stat_metrics_dict = {}
    for metric in stat_metrics:
        stat_metrics_dict[metric] = pd.DataFrame(index=range(number_of_b_samples), columns=samples, dtype=float)

    # initialise attribution dictionary
    signature_attributions_dict = {}
    for sample in samples:
        signature_attributions_dict[sample] = pd.DataFrame(index=range(number_of_b_samples), columns=signatures_for_boxplots, dtype=float)

    # mutation categories from signatures table
    categories = signatures.index.to_list()
    # initialise mutation spectra dictionary
    mutation_spectra_dict = {}
    for sample in samples:
        mutation_spectra_dict[sample] = pd.DataFrame(index=range(number_of_b_samples), columns=categories, dtype=float)

    # fill bootstrap dataframes for boxplots
    for i in range(number_of_b_samples):
        bootstrap_attribution_table = pd.read_csv(bootstrap_attribution_table_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        bootstrap_attribution_table_abs = pd.read_csv(bootstrap_attribution_table_abs_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        stat_table = pd.read_csv(bootstrap_stat_table_filename.replace('_i_', '_%i_' % (i+1)), index_col=0)
        # back-calculate mutation spectra from bootstrap attributions and signatures
        mutation_spectra = bootstrap_attribution_table_abs.dot(signatures.T)
        for sample in samples:
            for signature in signatures_for_boxplots:
                signature_attributions_dict[sample].loc[i, signature] = bootstrap_attribution_table.loc[sample, signature]
            for category in categories:
                mutation_spectra_dict[sample].loc[i, category] = mutation_spectra.loc[sample, category]
            for metric in stat_metrics:
                stat_metrics_dict[metric].loc[i, sample] = stat_table.loc[sample, metric]
                stat_metrics_dict[metric].loc[i, sample] = stat_table.loc[sample, metric]


    make_bootstrap_attribution_boxplot(central_attribution_table,
        title = main_title + ': central attributions',
        prefix = 'Optimised',
        y_label = 'Signature attribution',
        savepath = output_folder + '/' + filename + '_average.pdf')

    print(truth_attribution_table)
    # make bootstrap attribution plots for each sample
    for sample in samples:
        make_bootstrap_attribution_boxplot(signature_attributions_dict[sample],
            centrals = central_attribution_table.loc[sample].to_dict(),
            title = main_title + ': sample ' + str(sample).replace("_", "-"),
            truth = truth_attribution_table.loc[sample].to_dict() if 'SIM' in dataset_name else None,
            y_label = 'Signature attribution',
            savepath = output_folder + '/' + filename + '_' + str(sample) + '.pdf')
        make_bootstrap_attribution_boxplot(mutation_spectra_dict[sample],
            truth = input_mutations[str(sample)].to_dict(),
            title = main_title + ': sample ' + str(sample).replace("_", "-"),
            method = "NNLS (recalculated)",
            y_label = 'Mutation count',
            savepath = output_folder + '/bootstrap_plot_' + str(sample) + '_spectrum.pdf')
        for signature in signatures_for_boxplots:
            # do not plot signatures with both central (or truth) and mean bootstrap attributions of zero:
            if 'SIM' in dataset_name:
                if truth_attribution_table.loc[sample, signature] == 0 and np.mean(signature_attributions_dict[sample][signature] == 0):
                    signature_attributions_dict[sample].drop(signature, axis=1, inplace=True)
            elif central_attribution_table.loc[sample, signature] == 0 and np.mean(signature_attributions_dict[sample][signature] == 0):
                signature_attributions_dict[sample].drop(signature, axis=1, inplace=True)
        plot_bootstrap_attribution_histograms(signature_attributions_dict[sample],
            centrals = central_attribution_table.loc[sample].to_dict(),
            title = main_title + ': sample ' + str(sample).replace("_", "-"),
            savepath = output_folder + '/' + filename + '_' + str(sample) + '_dist.pdf')

    # make bootstrap plots for each statistical metric available
    for metric in stat_metrics:
        make_bootstrap_attribution_boxplot(stat_metrics_dict[metric],
            centrals = central_stat_table[metric].to_dict(),
            title = main_title + ': %s' % (metric),
            x_label = 'Samples',
            y_label = metric,
            method = metric, savepath = output_folder + '/' + metric + '_bootstrap_plot.pdf')
