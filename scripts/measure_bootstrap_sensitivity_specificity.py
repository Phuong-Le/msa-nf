""" measure_bootstrap_sensitivity_specificity.py
Module to measure and plot bootstrap sensitivity, specificity and other metrics
"""

from argparse import ArgumentParser
import copy
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from common_methods import make_folder_if_not_exists, read_data_from_JSON, calculate_stat_scores, calculate_confidence_interval
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
    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

def make_boxplot(input_table, title, xlabel, ylabel, show_mean=False, ylim_zero_to_one=True, savepath="./boxplot.pdf"):
    make_folder_if_not_exists(savepath.rsplit('/', 1)[0])

    table = copy.deepcopy(input_table)
    columns = input_table.columns.to_list()

    fig = plt.figure(figsize=(0.2*len(columns) if len(columns) > 30 else 6, 4))
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel(xlabel, size=12)
    ax.set_ylabel(ylabel, size=12)
    ax.tick_params(axis='x', rotation=90)
    if ylim_zero_to_one:
        ax.set_ylim(0, 1.05)
        plt.axhline(y=1, color='r', linestyle='--')
    table.boxplot(column=columns, ax=ax, grid=False, return_type='axes', sym='.',
                  meanline=show_mean, showmeans=show_mean, boxprops=dict(linewidth=3))
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
    parser.add_argument("-S", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                        help="set path to signature tables")
    parser.add_argument("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                        help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_argument("-d", "--dataset", dest="dataset_name", default='ESCC',
                        help="set the dataset name (e.g. ESCC, TCE)")
    parser.add_argument("-o", "--output_folder", dest="output_folder", default='plots/',
                        help="set path to save plots")
    parser.add_argument("-t", "--mutation_type", dest="mutation_type", default='',
                        help="set mutation type (SBS, DBS, ID)")
    parser.add_argument("-c", "--context", dest="context", default=192, type=int,
                        help="set SBS context (96, 192, 1536)")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true",
                        help="print additional information for debugging")
    parser.add_argument("-N", "--number_of_samples", dest="number_of_samples", default=-1, type=int,
                        help="limit the number of samples to analyse (all by default)")

    args = parser.parse_args()

    dataset_name = args.dataset_name
    mutation_type = args.mutation_type
    context = args.context
    number_of_samples = args.number_of_samples
    signature_tables_path = args.signature_tables_path
    signatures_prefix = args.signatures_prefix
    input_attributions_folder = args.input_attributions_folder + '/' + dataset_name + '/'
    output_folder = args.output_folder + '/' + dataset_name + '/' + mutation_type + '/bootstrap_plots/'
    make_folder_if_not_exists(output_folder)

    metrics = ['Cosine', 'Correlation', 'L1', 'L2', 'Chebyshev', 'Jensen-Shannon']
    scores = ['Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'F1', 'MCC']

    if mutation_type not in ['SBS', 'DBS', 'ID']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID" % mutation_type)

    if not input_attributions_folder:
        parser.error("Please specify the input path for reconstructed weights tables using -i option.")

    central_attribution_table_abs = pd.read_csv(input_attributions_folder + '/output_%s_%s_mutations_table.csv' % (dataset_name, mutation_type), index_col=0)
    central_attribution_table_weights = pd.read_csv(input_attributions_folder + '/output_%s_%s_weights_table.csv' % (dataset_name, mutation_type), index_col=0)

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

    if 'SIM' in dataset_name:
        if mutation_type == 'SBS':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.%i.weights.csv' % (dataset_name, context), index_col=0)
        elif mutation_type == 'DBS':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.dinucs.weights.csv' % dataset_name, index_col=0)
        elif mutation_type == 'ID':
            truth_attribution_table = pd.read_csv(input_attributions_folder + '/WGS_%s.indels.weights.csv' % dataset_name, index_col=0)

    signatures_to_consider = list(central_attribution_table_abs.columns)

    central_attribution_table = central_attribution_table_weights
    colormap_label = 'Relative contribution'
    filename = 'bootstrap_plot_weights'

    # limit the number of samples to analyse (if specified by -N option)
    if args.number_of_samples != -1:
        central_attribution_table = central_attribution_table.head(args.number_of_samples)
        if 'SIM' in dataset_name:
            truth_attribution_table = truth_attribution_table.head(args.number_of_samples)

    # convert columns and index to str
    central_attribution_table.columns = central_attribution_table.columns.astype(str)
    central_attribution_table.index = central_attribution_table.index.astype(str)
    if 'SIM' in dataset_name:
        truth_attribution_table.columns = truth_attribution_table.columns.astype(str)
        truth_attribution_table.index = truth_attribution_table.index.astype(str)

    main_title = dataset_name.replace('_', '/') + ' data, ' + mutation_type + ' mutation type'

    # samples and signatures to consider
    samples = central_attribution_table.index.to_list()
    signatures_to_consider = list(central_attribution_table_abs.columns)

    # read dictionaries from JSON files produced by make_bootstrap_tables.py script
    common_filename_suffix = "_%s_bootstrap_output_weights.json" % mutation_type
    attributions_per_sample_dict = read_data_from_JSON(input_attributions_folder + "attributions_per_sample" + common_filename_suffix)
    # attributions_per_signature_dict = read_data_from_JSON(input_attributions_folder + "attributions_per_signature" + common_filename_suffix)
    number_of_b_samples = len(next(iter(attributions_per_sample_dict.values())).index)

    if 'SIM' in dataset_name:
        acting_signatures = []
        for signature in list(truth_attribution_table.columns):
            if truth_attribution_table[signature].max() > 0:
                acting_signatures.append(signature)

    stat_scores_tables = pd.DataFrame(index=central_attribution_table.index, columns=scores, dtype=float)
    signatures_prevalences = pd.DataFrame(index=central_attribution_table.index, columns=signatures_to_consider, dtype=float)
    if 'SIM' in dataset_name:
        signatures_CPs = pd.DataFrame(index=central_attribution_table.index, columns=acting_signatures, dtype=float)
        truth_and_measured_difference = pd.DataFrame(0, index=central_attribution_table.index, columns=['Truth - mean', 'Truth - median', 'Truth - NNLS'], dtype=float)

    # calculate signature CPs
    if 'SIM' in dataset_name:
        for sample in samples:
            for signature in signatures_to_consider:
                array = attributions_per_sample_dict[sample][signature]
                confidence_interval = calculate_confidence_interval(array)
                if signature in acting_signatures:
                    if confidence_interval[0] <= truth_attribution_table.loc[sample, signature] <= confidence_interval[1]:
                        signatures_CPs.loc[sample, signature] = 1
                    else:
                        signatures_CPs.loc[sample, signature] = 0
                else:
                    if confidence_interval[0] <= truth_attribution_table.loc[sample, signature] <= confidence_interval[1]:
                        signatures_CPs.loc[sample, 'Others'] = 1
                    else:
                        signatures_CPs.loc[sample, 'Others'] = 0

                # fill residuals dataframe
                truth_and_measured_difference.loc[sample, 'Truth - mean'] += np.abs(truth_attribution_table.loc[sample, signature] - np.mean(array))
                truth_and_measured_difference.loc[sample, 'Truth - median'] += np.abs(truth_attribution_table.loc[sample, signature] - np.median(array))
                truth_and_measured_difference.loc[sample, 'Truth - NNLS'] += np.abs(truth_attribution_table.loc[sample, signature] - central_attribution_table.loc[sample, signature])


    # calculate CPs for a range of sensitivity thresholds for each sig
    if 'SIM' in dataset_name:
        signature_attribution_thresholds = [0.01 * i for i in range(20)]
        signatures_CPs_dict = {}
        for threshold in signature_attribution_thresholds:
            signatures_CPs_dict[threshold] = pd.DataFrame(index=central_attribution_table.index, columns=signatures_to_consider, dtype=float)
            for sample in samples:
                for signature in signatures_to_consider:
                    array = attributions_per_sample_dict[sample][signature]
                    confidence_interval = calculate_confidence_interval(array)
                    # increase upper bound of CI if less than threshold
                    if confidence_interval[0] == 0 and confidence_interval[1] < threshold:
                        confidence_interval[1] = threshold
                    if confidence_interval[0] <= truth_attribution_table.loc[sample, signature] <= confidence_interval[1]:
                        signatures_CPs_dict[threshold].loc[sample, signature] = 1
                    # elif confidence_interval[0] > 0 and confidence_interval[1] > 0 and truth_attribution_table.loc[sample, signature] > 0:
                    #     # this case includes cases when CI is not consistent with 0, signature is simulated but CI does not contain the true value
                    #     signatures_CPs_dict[threshold].loc[sample, signature] = 1
                    else:
                        signatures_CPs_dict[threshold].loc[sample, signature] = 0
            # normalise to calculate CPs
            signatures_CPs_dict[threshold] = signatures_CPs_dict[threshold].sum(axis=0).to_frame().T.div(len(samples))

    # plot sensitivity threshold vs CP plot for each sig
    if 'SIM' in dataset_name:
        for signature in signatures_to_consider:
            CPs = [signatures_CPs_dict[threshold][signature][0] for threshold in signature_attribution_thresholds]
            make_lineplot(signature_attribution_thresholds, CPs, 'Sensitivity threshold', 'Confidence probability', 'Signature %s' % signature, output_folder + '/CP_curves/signature_%s_curve.pdf' % signature)

    for sample in samples:
        for signature in signatures_to_consider:
            array = attributions_per_sample_dict[sample][signature]
            # disregard attributions consistent with 0
            # if confidence_interval[0]==0:
            # disregard attribution with 0 CI median
            if np.median(array) == 0:
                attributions_per_sample_dict[sample][signature] = 0

    for sample in samples:
        signatures_prevalences.loc[sample, :] = attributions_per_sample_dict[sample].astype(bool).sum(axis=0)
        if 'SIM' in dataset_name:
            truth_one_sample_table = truth_attribution_table.loc[sample, :].to_frame().T
            truth_table = truth_one_sample_table.loc[truth_one_sample_table.index.repeat(number_of_b_samples)]
            stat_scores_tables.loc[sample, :] = calculate_stat_scores(signatures, attributions_per_sample_dict[sample], truth_table)


    signatures_prevalences = signatures_prevalences/number_of_b_samples
    make_boxplot(signatures_prevalences, 'Signature prevalence (%s)' % mutation_type, 'Signatures', 'Reconstructed (%)', savepath=output_folder + '/signature_prevalences.pdf')

    if 'SIM' in dataset_name:
        signatures_CPs = signatures_CPs.sum(axis=0).to_frame().T.div(len(samples))
        make_boxplot(signatures_CPs, 'Confidence probability (%s)' % mutation_type, 'Signatures', 'CP (%)', savepath=output_folder + '/signature_confidence_probabilities.pdf')
        make_boxplot(stat_scores_tables, 'Metrics (%s)' % mutation_type, 'Metrics', 'Values', savepath=output_folder + '/sensitivity_specificity_metrics.pdf')
        make_boxplot(truth_and_measured_difference, 'Truth and measured difference (%s)' % mutation_type, '', 'Difference abs. sum', ylim_zero_to_one=False, savepath=output_folder + '/truth_and_measured_difference.pdf')
