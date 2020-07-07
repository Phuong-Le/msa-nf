""" plot_metrics.py
Module to plot various metrics such as sensitivity, specificity (if simulated data) and other metrics
Only signature prevalences are plotted on non-simulated data.
"""

from argparse import ArgumentParser
import copy
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from common_methods import make_folder_if_not_exists, read_data_from_JSON
matplotlib.use("agg")

def make_lineplot(x_values, y_values, xlabel, ylabel, title='', savepath="./lineplot.pdf"):
    make_folder_if_not_exists(savepath.rsplit('/', 1)[0])

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel(xlabel, size=12)
    ax.set_ylabel(ylabel, size=12)
    ax.tick_params(axis='x', rotation=90)

    plt.scatter(x_values, y_values)
    plt.plot(x_values, y_values, color='r', linestyle='--')
    # Add a grid to the plot
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    # Add a red line at 0.95
    plt.axhline(y=0.95, color='r', linestyle='-')
    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

def make_boxplot(input_table, title, xlabel, ylabel, show_mean=False, ylim_zero_to_one=True, add_jitter = False, savepath="./boxplot.pdf"):
    make_folder_if_not_exists(savepath.rsplit('/', 1)[0])

    table = copy.deepcopy(input_table)
    columns = input_table.columns.to_list()

    fig = plt.figure(figsize=(0.2*len(columns) if len(columns) > 30 else 6, 4))
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel(xlabel, size=12)
    ax.set_ylabel(ylabel, size=12)
    ax.tick_params(axis='y', which='minor', bottom=False)
    if ylim_zero_to_one:
        ax.set_ylim(0, 1.05)
        plt.axhline(y=1, color='r', linestyle='--')
        plt.axhline(y=0.95, color='g', linestyle='--')
    table.boxplot(column=columns, ax=ax, grid=False, return_type='axes', sym='.',
                  meanline=show_mean, showmeans=show_mean, boxprops=dict(linewidth=3))
    if add_jitter:
        for i, d in enumerate(table):
            y = table[d]
            x = np.random.normal(i+1, 0.04, len(y))
            plt.plot(x, y, mec='k', ms=3, marker="o", linestyle="None")

    plt.xticks(rotation=90)
    # Add a horizontal grid to the plot
    ax.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
    plt.minorticks_on()
    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

def make_efficiency_comparison_plot(input_data, title, xlabel, ylabel, savepath="./efficiencies.pdf"):
    make_folder_if_not_exists(savepath.rsplit('/', 1)[0])

    if isinstance(input_data, pd.DataFrame):
        methods = input_data.index.to_list()
        metrics = input_data.columns.to_list()
    else:
        methods = input_data.keys()
        metrics = next(iter(input_data.values())).columns.to_list()

    f, axes = plt.subplots(1, len(metrics), sharey=True, figsize=(len(metrics)*1.3, 4))
    f.suptitle(title, fontsize=12, y=0.999)
    axes[0].set_xlabel(xlabel, size=12)
    axes[0].set_ylabel(ylabel, size=12)

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
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_attributions_folder", dest="input_attributions_folder", default='output_tables/',
                        help="set path to NNLS output data")
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='ESCC',
                        help="set the dataset name (e.g. ESCC, TCE)")
    parser.add_argument("-o", "--output_folder", dest="output_folder", default='plots/',
                        help="set path to save plots")
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='',
                        help="set mutation type (SBS, DBS, ID)")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="print additional information for debugging")

    args = parser.parse_args()

    dataset_name = args.dataset_name
    mutation_type = args.mutation_type
    input_attributions_folder = args.input_attributions_folder + '/' + dataset_name + '/'
    output_folder = args.output_folder + '/' + dataset_name + '/' + mutation_type + '/bootstrap_plots/'
    make_folder_if_not_exists(output_folder)

    scores = ['Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'F1', 'MCC']

    if mutation_type not in ['SBS', 'DBS', 'ID']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID" % mutation_type)

    if not input_attributions_folder:
        parser.error("Please specify the input path for reconstructed weights tables using -i option.")

    central_attribution_table_abs = pd.read_csv(input_attributions_folder + '/output_%s_%s_mutations_table.csv' % (dataset_name, mutation_type), index_col=0)

    main_title = dataset_name.replace('_', '/') + ' data, ' + mutation_type + ' mutation type'

    # signatures to consider
    signatures_to_consider = list(central_attribution_table_abs.columns)

    # read truth-related dataframes and dictionaries of calculated variables
    if 'SIM' in dataset_name:
        truth_and_measured_difference = pd.read_csv(input_attributions_folder + '/truth_studies/truth_and_measured_difference_' + mutation_type + '.csv', index_col=0)
        signatures_CPs = pd.read_csv(input_attributions_folder + '/truth_studies/signatures_CPs_' + mutation_type + '.csv', index_col=0)
        sensitivity_thresholds = pd.read_csv(input_attributions_folder + '/truth_studies/sensitivity_thresholds_' + mutation_type + '.csv', index_col=0)
        stat_scores_from_CI_tables = pd.read_csv(input_attributions_folder + '/truth_studies/stat_scores_from_CI_tables_' + mutation_type + '.csv', index_col=0)
        signatures_CPs_dict = read_data_from_JSON(input_attributions_folder + '/truth_studies/signatures_CPs_dict_' + mutation_type + '.json')
        signatures_scores = read_data_from_JSON(input_attributions_folder + '/truth_studies/signatures_scores_' + mutation_type + '.json')
        stat_scores_from_CI_per_sig = read_data_from_JSON(input_attributions_folder + '/truth_studies/stat_scores_from_CI_per_sig_' + mutation_type + '.json')

    # plot sensitivity threshold vs CP plot for each sig
    if 'SIM' in dataset_name:
        for signature in signatures_to_consider:
            signature_attribution_thresholds = sorted(signatures_CPs_dict.keys())
            CPs = [signatures_CPs_dict[threshold][signature][0] for threshold in signature_attribution_thresholds]
            make_lineplot(signature_attribution_thresholds, CPs, 'Sensitivity threshold', 'Confidence probability', 'Signature %s' % signature, output_folder + '/truth_studies/CP_curves/signature_%s_curve.pdf' % signature)

    signatures_prevalences = pd.read_csv(input_attributions_folder + '/signatures_prevalences_' + mutation_type + '.csv', index_col=0)
    make_boxplot(signatures_prevalences, 'Signature prevalence (%s)' % mutation_type, 'Signatures', 'Reconstructed (%)', savepath=output_folder + '/signature_prevalences.pdf')

    if 'SIM' in dataset_name:
        make_boxplot(signatures_CPs, 'Confidence probability (%s)' % mutation_type, 'Signatures', 'CP (%)', savepath=output_folder + '/truth_studies/signature_confidence_probabilities.pdf')
        make_boxplot(sensitivity_thresholds, 'Sensitivity thresholds (%s)' % mutation_type, 'Signatures', 'Threshold (%)', ylim_zero_to_one = False, savepath=output_folder + '/truth_studies/sensitivity_thresholds.pdf')
        # make_boxplot(stat_scores_tables, 'Metrics (%s)' % mutation_type, 'Metrics', 'Values', savepath=output_folder + '/sensitivity_specificity_metrics.pdf')
        make_boxplot(stat_scores_from_CI_tables, 'Metrics (%s)' % mutation_type, 'Metrics', 'Values', savepath=output_folder + '/truth_studies/sensitivity_specificity_from_CI_metrics.pdf')
        for signature in signatures_to_consider:
            # make_boxplot(stat_scores_per_sig[signature], 'Signature %s metrics (%s)' % (signature, mutation_type), 'Metrics', 'Values', savepath=output_folder + '/per_signature_metrics/%s_metrics.pdf' % signature)
            make_boxplot(stat_scores_from_CI_per_sig[signature], 'Signature %s metrics (%s)' % (signature, mutation_type), 'Metrics', 'Values', savepath=output_folder + '/truth_studies/per_signature_from_CI_metrics/%s_metrics.pdf' % signature)
        for score in scores:
            make_boxplot(signatures_scores[score], '%s (%s)' % (score, mutation_type), 'Signatures', 'Values', savepath=output_folder + '/truth_studies/signature_%s.pdf' % score)
        make_boxplot(truth_and_measured_difference, 'Truth and measured difference (%s)' % mutation_type, '', 'Difference', ylim_zero_to_one=False, add_jitter = True, savepath=output_folder + '/truth_studies/truth_and_measured_difference.pdf')
