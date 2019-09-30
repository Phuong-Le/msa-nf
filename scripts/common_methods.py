import os, math, copy
import numpy as np
from scipy.spatial import distance
import matplotlib.pyplot as plt
import scipy.stats as stats

def make_folder_if_not_exists(folder):
    if not os.path.exists(folder):
        try:
            os.makedirs(folder)
        except:
            warnings.warn("Could not create a folder ", folder)

def get_comma_separated_args(option, opt, value, parser):
    setattr(parser.values, option.dest, value.split(','))

def plot_array_as_histogram(arrays, labels, title, savepath='./hist.pdf'):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # bins = [0.025*k for k in range(42)]
    bins = [0.01*k for k in range(101)]
    for array, label in zip(arrays, labels):
        plt.hist(array, label = label, histtype='step', bins = bins)

    legend = plt.legend()
    plt.gca().add_artist(legend)

    ax.set_title(title, fontsize=14, pad=10)
    ax.set_xlabel('Attribution ($\%$)', fontsize=12)
    ax.set_ylabel('Number of attributions', fontsize=12)
    ax.set_xticks([0.1*k for k in range(11)])

    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

def calculate_similarity(first_sample, second_sample, metric='Cosine'):
    """
    Calculate various similarity metrics as 1 minus distance functions provided
    by scipy.spatial package (see https://docs.scipy.org/doc/scipy/reference/spatial.distance.html).
    If input samples are zeroes, the similarity can be undefined (nan) for some metrics.

    Parameters:
    ----------
    first_sample: array_like
        First input array (e.g. reconstructed sample)
    sedond_sample: array_like
        Second input array (e.g. true sample)
    metric: string
        One of the supported metrics: Jensen-Shannon, Cosine, Correlation, Chebyshev, L1, L2, L3
        (see https://docs.scipy.org/doc/scipy/reference/spatial.distance.html)

    Returns:
    -------
    similarity: double.
        Smaller distance gives similarity closer to one.
        Note: distance is not always normalised to 1 (e.g. JS vs minkowski metrics).
    -------
    """
    if 'Cosine' in metric or 'cosine' in metric:
        similarity = 1 - distance.cosine(first_sample, second_sample)
    elif 'Correlation' in metric or 'correlation' in metric:
        similarity = 1 - distance.correlation(first_sample, second_sample)
    elif 'Chebyshev' in metric or 'chebyshev' in metric:
        similarity = 1 - distance.chebyshev(first_sample, second_sample)
    elif 'L1' in metric or 'Manhattan' in metric or 'manhattan' in metric:
        similarity = 1 - distance.minkowski(first_sample, second_sample, p=1)
    elif 'L2' in metric or 'Euclidean' in metric or 'euclidean' in metric:
        similarity = 1 - distance.euclidean(first_sample, second_sample)
    elif 'L3' in metric:
        similarity = 1 - distance.minkowski(first_sample, second_sample, p=3)
    elif 'Jensen-Shannon' in metric or 'jensen-shannon' in metric:
        similarity = 1 - distance.jensenshannon(first_sample, second_sample)
    elif 'Canberra' in metric or 'canberra' in metric:
        similarity = 1 - distance.canberra(first_sample, second_sample)
    else:
        raise ValueError("Unknown metric: ", metric)

    return similarity

def calculate_stat_scores(signatures, reco_table, truth_table, number_of_samples = None):
    """
    Calculate various statistical scores.

    Parameters:
    ----------
    signatures: list
        List of strings with signature names.
        Should be available in input reco and truth tables as columns.
    reco_table: pandas dataframe
        Dataframe containing the signature attribution table for a set of samples
    truth_table: pandas dataframe
        Dataframe containing the true signature activities (known from simulations)
        for the same set of samples (order has to be the same)
    number_of_samples: int or None
        If provided as an int, cap the number of samples by this number.
        Otherwise, all provided samples are used

    Returns:
    -------
    sensitivity, specificity, precision, accuracy, F1, MCC: list of doubles
    List of various scores (see e.g. https://en.wikipedia.org/wiki/Confusion_matrix for info)
    MCC stands for Matthews correlation coefficient.
    -------
    """
    true_positives = true_negatives = false_positives = false_negatives = 0
    assert len(reco_table) == len(truth_table)

    true_positives_signatures = []
    false_negatives_signatures = []
    false_positives_signatures = []

    if not number_of_samples:
        number_of_samples = len(reco_table)

    for i in range(number_of_samples):
        reco_sample = reco_table.iloc[i]
        truth_sample = truth_table.iloc[i]
        for sig in signatures:
            if truth_sample[sig] > 0 and reco_sample[sig] > 0:
                true_positives += 1
                true_positives_signatures.append(truth_sample[sig])
            elif truth_sample[sig] == 0 and reco_sample[sig] == 0:
                true_negatives += 1
            elif truth_sample[sig] > 0 and reco_sample[sig] == 0:
                false_negatives += 1
                false_negatives_signatures.append(truth_sample[sig])
            elif truth_sample[sig] == 0 and reco_sample[sig] > 0:
                false_positives += 1
                false_positives_signatures.append(reco_sample[sig])
            else:
                raise ValueError(
                    "Invalid signature attribution values -- please check the tables")
    sensitivity = true_positives / (true_positives + false_negatives)
    specificity = true_negatives / (true_negatives + false_positives)
    precision = true_positives / (true_positives + false_positives)

    # print("Missed sigs mean/median/stdev/max:",np.mean(false_negatives_signatures),np.median(false_negatives_signatures),np.std(false_negatives_signatures),np.max(false_negatives_signatures))
    # plot_array_as_histogram([true_positives_signatures, false_positives_signatures, false_negatives_signatures], ['True positives', 'False positives', 'False negatives'], title = 'Signature attribution distributions', savepath="distributions.pdf")

    accuracy = (true_positives + true_negatives) / (true_positives + true_negatives + false_positives + false_negatives)
    F1 = 2 * true_positives / (2 * true_positives + false_positives + false_negatives)
    MCC = (true_positives * true_negatives - false_positives * false_negatives) / math.sqrt((true_positives + false_positives)
            * (true_positives + false_negatives) * (true_negatives + false_positives) * (true_negatives + false_negatives))

    return sensitivity, specificity, precision, accuracy, F1, MCC
