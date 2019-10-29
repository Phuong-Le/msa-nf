from optparse import OptionParser
import os
import warnings
import pandas as pd
import numpy as np
import copy
import time
from scipy.optimize import nnls
from scipy.spatial import distance
from common_methods import make_folder_if_not_exists

def bootstrap_mutation_table(input_dataframe, method="classic"):
    input_mutations = copy.deepcopy(input_dataframe)
    if method=="classic":
        # classic bootstrap (resampling with replacement, i.i.d. assumption)
        original_index = input_mutations.index
        # bootstrap using pd.sample
        bootstrap_mutations = input_mutations.sample(n=len(input_mutations.index), replace=True)
        # keep the old index (mutation categories)
        bootstrap_mutations.index = original_index
    elif method=="poisson":
        # Poisson bootstrap (n>100 and i.i.d. assumption)
        bootstrap_mutations = input_mutations.applymap(lambda x: x*np.random.poisson(1))
    elif method=="binomial" or "multinomial" in method:
        # mutational burdens for each sample (column) for normalisation
        sums = input_mutations.sum(axis=0)
        normalised_input_mutations = input_mutations.div(sums, axis=1)
        bootstrap_mutations = pd.DataFrame(index = input_mutations.index, columns = input_mutations.columns)
        for i in range(len(bootstrap_mutations.columns)):
            if method=="binomial":
                # simple independent binomials for each category
                bootstrap_mutations.iloc[:,i] = np.random.binomial(sums.iloc[i], normalised_input_mutations.iloc[:,i])
            elif method=="multinomial":
                # Weighted multinomial with fixed number of draws equalling mutational burden
                bootstrap_mutations.iloc[:,i] = np.random.multinomial(sums.iloc[i], normalised_input_mutations.iloc[:,i])
            elif method=="multinomial_weight":
                # Mulinomial counts with equal probabilities, equivalent to classic bootstrap
                n = len(normalised_input_mutations.iloc[:,i])
                counts = np.random.multinomial(n, [1/n for i in range(n)])
                input_data = input_mutations.iloc[:,i]
                bootstrap = np.concatenate([[input_data[j]]*counts[j] for j in range(len(input_data))]).ravel().astype(int)
                np.random.shuffle(bootstrap)
                bootstrap_mutations.iloc[:,i] = bootstrap
            else:
                raise ValueError("Unknown bootstrap method: %s" % method)
    else:
        raise ValueError("Unknown bootstrap method: %s" % method)
    return bootstrap_mutations

def optimise_signatures(selected_mutations, initial_signatures, weak_threshold=0.001, strong_threshold=0.001, verbose=False):
    """
    Perform signature optimisation for NNLS attribution method.
    The method is outlined in the PCAWG paper, with the main idea as follows:
    The NNLS is intially applied on a given input set of signatures and mutations.
    The similarity of the reconstructed profile (linear combination of input signatures)
    to the initial sample is calculated, called base similarity (last column in stat_info dataframe).
    Afterwards, two loops are executed: in the first one, all signatures that decrease the
    similarity by less than a given threshold (weak_threshold) are removed; in the second one,
    all excluded so far signatures that increase the similarity by more than a given
    threshold (strong_threshold) are included again.
    The output is the final list of remaining signatures.

    Parameters:
    ----------
    selected_mutations: list of ints
        List of mutation counts for a single sample. Example: list of 96 integers,
        for a 96-context SBS mutations.

    initial_signatures: pandas dataframe
        Dataframe of input signatures with the index of mutation categories, with the order
        that has to be the same as for the input mutations.

    weak_threshold: float
        Similarity decrease threshold to exclude weakest signatures (default: 0.01)

    strong_threshold: float
        Similarity increase threshold to include strongest signatures (default: 0.05)

    verbose: boolean
        Verbosity flag: lots of output for debugging if set to True

    Returns:
    -------
    final_signatures: list of strings
        The final list of remaining signatures passing the optimisation.
    -------
    """
    # Signature optimisation routine starts here
    # calculate the base similarity with an initial set of signatures
    _, _, stat_info = perform_signature_attribution(selected_mutations, initial_signatures)
    base_similarity = stat_info[-1]
    significant_signatures = copy.deepcopy(initial_signatures)

    # a loop to remove all weak signatures
    all_weak_signatures_removed = False
    while all_weak_signatures_removed == False:
        signatures_contribution = {}
        if verbose:
            print('Current significant signatures:')
            print(significant_signatures.columns.tolist())
        for signature in significant_signatures.columns.tolist():
            signatures_to_run = copy.deepcopy(significant_signatures)
            signatures_to_run = signatures_to_run.drop(signature, axis=1)
            if verbose:
                print('Running without signature:', signature)
            _, _, stat_info = perform_signature_attribution(selected_mutations, signatures_to_run)
            new_similarity = stat_info[-1]
            contribution = base_similarity-new_similarity
            signatures_contribution[signature] = contribution

        if verbose:
            print('Signatures contribution for sample', sample)
            print(signatures_contribution)

        least_contributing_signature = min(signatures_contribution, key=signatures_contribution.get)
        if signatures_contribution[least_contributing_signature]<weak_threshold and len(significant_signatures.columns.tolist())>1:
            significant_signatures = significant_signatures.drop(least_contributing_signature, axis=1)
            if verbose:
                print('Dropping signature %s from sample %s' % (least_contributing_signature, sample))
                print('Number of signatures left:',len(significant_signatures.columns.tolist()))
        else:
            if verbose:
                print('All cleaned up!')
            all_weak_signatures_removed = True

        # check if only one significant signature is left
        if len(significant_signatures.columns.tolist())==1:
            if verbose:
                print('Only one signature left:',significant_signatures.columns.tolist())
            all_weak_signatures_removed = True

    # calculate the base similarity again with a given set of significant signatures
    _, _, stat_info = perform_signature_attribution(selected_mutations, significant_signatures)
    base_similarity = stat_info[-1]

    # final signatures to update
    final_signatures = copy.deepcopy(significant_signatures)
    # remaining signatures to loop through in search for most contributing ones
    remaining_signatures = initial_signatures.drop(significant_signatures.columns, axis=1)

    # a loop to add strong signatures if there is anything to add
    if remaining_signatures.empty:
        all_strong_signatures_added = True
    else:
        all_strong_signatures_added = False

    while all_strong_signatures_added == False:
        signatures_contribution = {}
        if verbose:
            print('Current significant signatures:')
            print(final_signatures.columns.tolist())

        # loop through the remaining signatures to add contributing ones
        for signature in remaining_signatures.columns.tolist():
            signatures_to_run = copy.deepcopy(final_signatures)
            signatures_to_run[signature] = pd.Series(remaining_signatures[signature], index = signatures_to_run.index)
            if verbose:
                print('Running with signature:', signature)
            _, _, stat_info = perform_signature_attribution(selected_mutations, signatures_to_run)
            new_similarity = stat_info[-1]
            contribution = new_similarity-base_similarity
            signatures_contribution[signature] = contribution

        if verbose:
            print('Signatures contribution for sample', sample)
            print(signatures_contribution)

        most_contributing_signature = max(signatures_contribution, key = signatures_contribution.get)
        if signatures_contribution[most_contributing_signature]>strong_threshold:
            if verbose:
                print('Adding signature %s to sample %s' % (most_contributing_signature, sample))
            final_signatures[most_contributing_signature] = pd.Series(remaining_signatures[most_contributing_signature], index = final_signatures.index)
            remaining_signatures = remaining_signatures.drop(most_contributing_signature, axis=1)
            if remaining_signatures.empty:
                if verbose:
                    print('All strong signatures added!')
                all_strong_signatures_added = True
        else:
            if verbose:
                print('All strong signatures added!')
            all_strong_signatures_added = True

    if verbose:
        print('Final significant signatures:')
        print(final_signatures.columns.tolist())

    return final_signatures

def perform_signature_attribution(selected_mutations, signatures, verbose=False):
    """
    Apply the NNLS attribution method using NNLS function from scipy.optimize:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html

    Parameters:
    ----------
    selected_mutations: list of ints
        List of mutation counts for a single sample. Example: list of 96 integers,
        for a 96-context SBS mutations. The order has to be identical to input signatures

    signatures: pandas dataframe
        Dataframe of input signatures with the index of mutation categories, with the order
        that has to be the same as for the input mutations

    verbose: boolean
        Verbosity flag: lots of output for debugging if set to True

    Returns:
    -------
    normalised_weights, mutation_numbers, stat_info: tuple
    normalised_weights: list of floats
        The list weights normalised by the sum of weights (i.e. all add up to one).
        Corresponds to a probability of all signatures to contribute to the given sample
    mutation_numbers: list of floats
        The list of normalised weights mentioned above, multiplied by total mutational burden
    stat_info: a list of floats
        A list of statistical info variables, as follows:
        rss: residual sum of squares (RSS),
        chi2: Chi-2 of the fit,
        r2: determination coefficient (R2)
        cosine_similarity: cosine similarity of the reconstructed profile with the input sample
        correlation: correlation of the reconstructed profile with the input sample
        chebyshev_similarity: Chebyshev similarity for mentioned profiles
        L1/L2/L3 similarity: 1 minus L1/L2/L3 Minkowski distances between mentioned profiles
        jensenshannon_similarity: Jensen-Shannon similarity of mentioned profiles
    -------
    """
    reg = nnls(signatures, selected_mutations)
    weights = reg[0]

    normalised_weights = weights/sum(weights)
    mutation_numbers = normalised_weights*sum(selected_mutations)

    # calculate RSS, chi2, r2, similarity metrics
    observed = selected_mutations
    fitted = np.matmul(signatures.values,weights)
    residuals = observed - fitted
    chi2 = sum(residuals*residuals/fitted)

    mean = np.mean(observed)
    tot_squares = (observed-mean)*(observed-mean)
    tss = sum(tot_squares)
    rss = sum(residuals*residuals)
    r2 = 1-rss/tss

    cosine_similarity = 1 - distance.cosine(observed, fitted)
    correlation = 1 - distance.correlation(observed, fitted)
    chebyshev_similarity = 1 - distance.chebyshev(observed, fitted)
    L1_similarity = 1 - distance.minkowski(observed, fitted, p=1)
    L2_similarity = 1 - distance.euclidean(observed, fitted)
    L3_similarity = 1 - distance.minkowski(observed, fitted, p=3)
    jensenshannon_similarity = 1 - distance.jensenshannon(observed, fitted)

    stat_info = [rss, chi2, r2, cosine_similarity, chebyshev_similarity, L1_similarity, L2_similarity, L3_similarity, jensenshannon_similarity]

    if verbose:
        if chi2>10e10:
            print('************* High chi2 sample *************')
        print('Signatures:')
        print(signatures.columns.tolist())
        print('[rss, chi2, r2, cosine_similarity, chebyshev_similarity, L1_similarity, L2_similarity, L3_similarity, jensenshannon_similarity]:',stat_info)
        print('Observed:', observed)
        print('Sum observed:', np.sum(observed))
        print('Fitted:', fitted)
        print('Sum Fitted:', np.sum(fitted))
        print('Residuals:', residuals)

    return normalised_weights, mutation_numbers, stat_info

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-t", "--mutation_type", dest="mutation_type", default='',
                      help="set mutation type (SBS, DBS, ID)")
    parser.add_option("-c", "--context", dest="context", default=96, type='int',
                      help="set SBS context (96, 192)")
    parser.add_option("-i", "--input_path", dest="input_path", default='input_mutation_tables/',
                      help="set path to input mutation tables")
    parser.add_option("-s", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                      help="set path to signature tables")
    parser.add_option("-p", "--signature_prefix", dest="signatures_prefix", default='sigProfiler',
                      help="set prefix in signature filenames (sigProfiler by default)")
    parser.add_option("-o", "--output_path", dest="output_path", default='output_tables/',
                      help="set path to save output tables")
    parser.add_option("-x", "--optimise_signatures", dest="optimise_signatures", action="store_true",
                      help="perform signature optimisation (remove weak signature, add strong ones)")
    parser.add_option("-W", "--weak_threshold", dest="weak_threshold", default=0.002, type='float',
                      help="Similarity decrease threshold to exclude weakest signatures in optimisation (default: 0.001)")
    parser.add_option("-S", "--strong_threshold", dest="strong_threshold", default=0.002, type='float',
                      help="Similarity increase threshold to include strongest signatures in optimisation (default: 0.001)")
    parser.add_option("-d", "--dataset", dest="dataset_name", default='SIM',
                      help="set the dataset name (SIM by default)")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      help="verbosity flag for debugging (lots of output)")
    parser.add_option("-n", "--number", dest="number_of_samples", default=-1, type='int',
                      help="limit the number of samples to analyse (all by default)")
    parser.add_option("-B", "--bootstrap", dest="bootstrap", action="store_true",
                      help="Reshuffle input mutation tables for bootstrapping")
    parser.add_option("--bootstrap_method", dest="bootstrap_method", default='binomial',
                      help="choose a method for bootstrapping (perturbing) samples (classic, binomial, multinomial)")
    parser.add_option("--add_suffix", dest="add_suffix", action="store_true",
                      help="add suffixes to output folders (useful during optimisation scan analysis)")

    (options, args) = parser.parse_args()

    dataset_name = options.dataset_name
    mutation_type = options.mutation_type
    context = options.context
    input_path = options.input_path
    signature_tables_path = options.signature_tables_path
    signatures_prefix = options.signatures_prefix
    output_path = options.output_path + '/' + dataset_name

    if options.add_suffix:
        output_path += '_%i_NNLS' % context
        if options.optimise_signatures:
            output_path += '_%.4f_%.4f' % (options.weak_threshold, options.strong_threshold)
        else:
            output_path += '_unoptimised'

    make_folder_if_not_exists(output_path)

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (or '-t DBS', '-t ID').")
    elif mutation_type not in ['SBS','DBS','ID']:
        raise ValueError("Error: Unknown mutation type: %s. Known types: SBS, DBS, ID" % mutation_type)

    if mutation_type=='SBS':
        if context==96:
            signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=[0,1])
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_path, dataset_name, dataset_name, context), index_col=[0,1])
        elif context==192:
            signatures = pd.read_csv('%s/%s_%s_192_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=[0,1,2])
            input_mutations = pd.read_csv('%s/%s/WGS_%s.%i.csv' % (input_path, dataset_name, dataset_name, context), index_col=[0,1,2])
        else:
            raise ValueError("Context %i is not supported." % context)
    elif mutation_type=='DBS':
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
        input_mutations = pd.read_csv('%s/%s/WGS_%s.dinucs.csv' % (input_path, dataset_name, dataset_name), index_col=0 if 'SIM' in dataset_name else [0,1])
    elif mutation_type=='ID':
        signatures = pd.read_csv('%s/%s_%s_signatures.csv' % (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
        input_mutations = pd.read_csv('%s/%s/WGS_%s.indels.csv' % (input_path, dataset_name, dataset_name), index_col=0 if 'SIM' in dataset_name else [0,1,2,3])

    print("Performing NNLS for %s dataset, %s mutation type." % (dataset_name, mutation_type))
    if mutation_type=='SBS':
        print("SBS context: %i" % context)
    if options.optimise_signatures:
        print("Optimised NNLS method, weak/strong thresholds: %f/%f" % (options.weak_threshold, options.strong_threshold))
    if options.bootstrap:
        print("Perturbing the input mutation sample, bootstrap (simulation) method: %s" % options.bootstrap_method)

    if options.bootstrap:
        input_mutations = bootstrap_mutation_table(input_mutations, method=options.bootstrap_method)

    # limit the number of samples to analyse (if specified by -n option)
    if options.number_of_samples!=-1:
        input_mutations = input_mutations.iloc[:,0:options.number_of_samples]

    samples = input_mutations.columns

    n_types = signatures.shape[0]
    num_ref_sigs = signatures.shape[1]

    sel_sig_nums=list(range(0,num_ref_sigs))
    # sel_sig_nums=[0,4,26,44]
    # sel_sig_nums=range(0,44)+range(46,65)  # all 65 signatures except SBS40

    signature_columns_list = signatures.columns[sel_sig_nums].tolist()
    print("Analysing signatures:", signature_columns_list)
    # print("Analysing samples:", input_mutations.columns.tolist())

    output_weights = pd.DataFrame(index=samples, columns=signature_columns_list)
    output_mutations = pd.DataFrame(index=samples, columns=signature_columns_list)
    output_stat_info = pd.DataFrame(index=samples, columns=['RSS', 'Chi2', 'R2', 'Cosine similarity', 'Chebyshev similarity', 'L1 similarity', 'L2 similarity', 'L3 similarity', 'Jensen-Shannon similarity'])

    output_weights.index.name = output_mutations.index.name = output_stat_info.index.name = 'Sample'

    start_time = time.process_time()

    # sort indexes to prevent signatures/samples mismatch
    input_mutations.sort_index(inplace = True)
    signatures.sort_index(inplace = True)

    # print(input_mutations.index.tolist())
    # print(signatures.index.tolist())
    for sample in input_mutations.columns:
        print("Processing sample", sample)
        selected_mutations = input_mutations[sample].tolist()
        if not sum(selected_mutations)>0:
            warnings.warn("Zero mutations: skipping")
            continue

        initial_signatures = signatures.iloc[:,sel_sig_nums]

        if options.optimise_signatures:
            final_signatures = optimise_signatures(selected_mutations, initial_signatures, options.weak_threshold, options.strong_threshold, verbose=options.verbose)
        else:
            # skip signature optimisation
            final_signatures = initial_signatures

        if not final_signatures.empty:
            normalised_weights, mutation_numbers, stat_info = perform_signature_attribution(selected_mutations, final_signatures, verbose=options.verbose)
            signatures_to_fill = final_signatures.columns
            output_weights.loc[sample,signatures_to_fill] = normalised_weights
            output_mutations.loc[sample,signatures_to_fill] = mutation_numbers
            output_stat_info.loc[sample] = stat_info

    end_time = time.process_time()
    process_elapsed_time = end_time - start_time
    print("Signature attribution took %.2f s (%.2f s per sample)" % (process_elapsed_time, process_elapsed_time/len(samples.tolist())))

    # replace NANs by zeros
    output_weights = output_weights.fillna(0.0)
    output_mutations = output_mutations.fillna(0.0)
    output_stat_info = output_stat_info.fillna(0.0)

    # write to files
    output_weights.to_csv(output_path + '/output_%s_weights_table.csv' % mutation_type)
    output_mutations.to_csv(output_path + '/output_%s_mutations_table.csv' % mutation_type)
    output_stat_info.to_csv(output_path + '/output_%s_stat_info.csv' % mutation_type)
