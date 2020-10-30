from optparse import OptionParser
import os
import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from common_methods import make_folder_if_not_exists

signatures_to_generate = {
    # Dictionary with signature burdens to be generated.
    # Each signature is assigned with a list of type [mu, sigma], where mu is the mean
    # and sigma is the standard deviation of the normal distribution used for simulating
    # the mutational burden attributed to the signature.
    # Signature names (strings) have to be present in the input signature tables.
    'SBS1':[150, 100],
    'SBS2':[150, 130],
    'SBS5':[0, 50],
    'SBS13':[150, 130],
    'SBS18':[100, 100],
    'SBS21':[0, 50],
    'SBS40':[200, 300],
    # 'ID1':[140,270],
    # 'ID2':[2000,1000],
}

def plot_mutational_burden(mutational_burden, mu=None, sigma=None, title='Total', savepath = './burden.pdf'):
    """
    Plot a histogram of the input mutational burden, plot the Gaussian on top if
    the mean (mu) and standard deviation (sigma) are provided, otherwise fit a
    Gaussian to the histogram.

    Parameters:
    ----------
    mutational_burden: array_like
        list of numbers (ints or floats)

    mu: float (Optional)
        Mean of the normal distribution to plot

    sigma: float (Optional)
        Standard deviation of the normal distribution to plot

    Returns:
    -------
    Nothing. Creates and saves a plot in specified path.
    -------
    """
    make_folder_if_not_exists(savepath.rsplit('/',1)[0])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    n, bins, patches = plt.hist(mutational_burden, density=True, bins=20, color='blue', histtype='step')
    if mu and sigma:
        # create a Gaussian with given mu and sigma
        y = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bins - mu)**2 / (2 * sigma**2) )
        text = '$\\mu$=%.2f, $\\sigma=$%.2f' %  (mu, sigma)
    else:
        # fit a Gaussian to input distribution
        (mu, sigma) = stats.norm.fit(mutational_burden)
        y = stats.norm.pdf( bins, mu, sigma)
        text = 'Fitted Gaussian:\n $\\mu$=%.2f, $\\sigma=$%.2f' %  (mu, sigma)

    number_of_positive_samples = sum(x > 0 for x in mutational_burden)
    text += '\n %i/%i samples' % (number_of_positive_samples, len(mutational_burden))

    ax.plot(bins, y, 'r--', linewidth=2)
    # remove indentation from latex-based text
    ax.text(0.65, 0.9, text, fontsize=9, transform=ax.axes.transAxes)

    ax.set_ylabel("Probability", fontsize=12)
    ax.set_xlabel("Mutation count", fontsize=12)
    ax.set_title(title, fontsize=14, pad=10)
    plt.tight_layout()
    plt.savefig(savepath, transparent=True)
    plt.close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-t", "--mutation_type", dest="mutation_type", default='',
                      help="set mutation type (SBS, DBS, ID)")
    parser.add_option("-c", "--context", dest="context", default=96, type='int',
                      help="set SBS context (96, 192)")
    parser.add_option("-s", "--signature_path", dest="signature_tables_path", default='signature_tables/',
                      help="set path to signature tables")
    parser.add_option("-o", "--output_path", dest="output_path", default='input_mutation_tables/',
                      help="set path to save output simulatied mutation tables")
    parser.add_option("-d", "--dataset", dest="dataset_name", default='SIM',
                      help="set the dataset name ('SIM' by default)")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      help="verbosity flag for debugging (lots of output)")
    parser.add_option("-n", "--number_of_samples", dest="number_of_samples", default=100, type='int',
                      help="set the number of samples to generate (100 by default)")
    parser.add_option("-z", "--noise", dest="add_noise", action="store_true",
                      help="Add noise of type specified by noise_type parameter")
    parser.add_option("--noise_type", dest="noise_type", default='poisson',
                      help="Choose the type of noise: Poisson, Gaussian or negative_binomial (Poisson variation by default)")
    parser.add_option("-Z", "--noise_sigma", dest="noise_sigma", default=2, type='float',
                      help="Set standard deviation of Gaussian noise if used (2 by default)")
    parser.add_option("-r", "--random", dest="random_signatures", action="store_true",
                      help="Use randomly generated signatures instead of PCAWG reference ones")
    parser.add_option("-N", "--number_of_random_sigs", dest="number_of_random_sigs", default=5, type='int',
                      help="Number of random signatures to consider")
    parser.add_option("-B", "--bootstrap_input_mutation_table", dest="bootstrap_input_mutation_table", action="store_true",
                      help="Bootstrap (reshuffle with replacement) the input mutation table instead if using signatures")
    parser.add_option("-i", "--input_table", dest="input_table", default='',
                      help="set path to input mutation table to bootstrap")
    (options, args) = parser.parse_args()

    mutation_type = options.mutation_type
    dataset_name = options.dataset_name
    context = options.context
    signature_tables_path = options.signature_tables_path
    output_path = options.output_path + '/' + dataset_name
    number_of_samples = options.number_of_samples
    random_signatures = options.random_signatures
    number_of_random_sigs = options.number_of_random_sigs

    make_folder_if_not_exists(output_path)

    if not mutation_type:
        parser.error("Please specify the mutation type using -t option, e.g. add '-t SBS' to the command (or '-t DBS', '-t ID').")
    elif mutation_type not in ['SBS','DBS','ID']:
        raise ValueError("Unknown mutation type: %s. Known types: SBS, DBS, ID" % mutation_type)

    if random_signatures:
        signatures_prefix = 'sigRandom'
    else:
        signatures_prefix = 'sigProfiler'

    if mutation_type=='SBS':
        output_filename = '%s/WGS_%s.%i.csv' % (output_path, dataset_name, context)
        if context==96:
            reference_signatures = pd.read_csv('%s/%s_%s_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type), index_col=[0,1])
        elif context in [192, 288]:
            reference_signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type, context), index_col=[0,1,2])
        elif context==1536:
            reference_signatures = pd.read_csv('%s/%s_%s_%i_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type, context), index_col=0)
        else:
            raise ValueError("Context %i is not supported." % context)
    elif mutation_type=='DBS':
        output_filename = '%s/WGS_%s.dinucs.csv' % (output_path, dataset_name)
        reference_signatures = pd.read_csv('%s/%s_%s_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type), index_col=0)
    elif mutation_type=='ID':
        output_filename = '%s/WGS_%s.indels.csv' % (output_path, dataset_name)
        reference_signatures = pd.read_csv('%s/%s_%s_signatures.csv' %
                                            (signature_tables_path, signatures_prefix, mutation_type), index_col=0)

    for signature in reference_signatures.columns:
        if not np.isclose(reference_signatures.sum()[signature], 1, rtol=1e-2):
            raise ValueError("Probabilities for signature %s do not add up to 1: %.3f" % (signature, reference_signatures.sum()[signature]))

    generated_signature_burdens = {}
    if options.bootstrap_input_mutation_table:
        if not options.input_table:
            parser.error("Please provide the input mutation table for bootstrap with -i option.")
        input_table = pd.read_csv(options.input_table, index_col=0, sep=None)
        bootstrapped_table = input_table.sample(n=number_of_samples, replace=True)
        for signature in input_table.columns:
            generated_signature_burdens[signature] = bootstrapped_table[signature].to_numpy()
            plot_mutational_burden(generated_signature_burdens[signature], title=signature,
            savepath='%s/%s_plots/generated_burden_%s.pdf' % (output_path, dataset_name, signature))
    else:
        if random_signatures:
            random_signature_indexes = random.sample(reference_signatures.columns.to_list(), number_of_random_sigs)
            for random_sig in random_signature_indexes:
                mu = np.random.normal(2000, 1000)
                sigma = np.random.normal(500, 100)
                mu = mu if mu>=0 else 0
                sigma = sigma if sigma>=0 else 0
                # draw burdens from normal distribution
                generated_burdens = np.random.normal(mu, sigma, number_of_samples)
                # replace negative values of the Gaussian with zeros
                generated_burdens[generated_burdens<0] = 0
                # plot the generated burden
                plot_mutational_burden(generated_burdens, mu, sigma, random_sig,
                                    savepath='%s/%s_plots/generated_burden_%s.pdf' % (output_path, dataset_name, random_sig))
                generated_signature_burdens[random_sig] = generated_burdens
        else:
            # generate burdens for each signature from the 'signatures_to_generate' dictionary
            for signature in signatures_to_generate.keys():
                mu = signatures_to_generate[signature][0]
                sigma = signatures_to_generate[signature][1]
                # draw burdens from normal distribution
                generated_burdens = np.random.normal(mu, sigma, number_of_samples)
                # replace negative values of the Gaussian with zeros
                generated_burdens[generated_burdens<0] = 0
                # plot the generated burden
                plot_mutational_burden(generated_burdens, mu, sigma, signature,
                                    savepath='%s/%s_plots/generated_burden_%s.pdf' % (output_path, dataset_name, signature))
                generated_signature_burdens[signature] = generated_burdens

    samples_range = range(0,number_of_samples)

    generated_weights = pd.DataFrame(0, index=samples_range, columns=reference_signatures.columns)
    generated_mutations = pd.DataFrame(0, index=reference_signatures.index, columns=samples_range)

    # loop to fill mutation tables:
    generated_mutational_burdens = []
    for i in samples_range:
        mutational_burden = sum(signature_burden[i] for signature_burden in list(generated_signature_burdens.values()))
        generated_mutational_burdens.append(mutational_burden)
        weights = {}
        for signature in generated_signature_burdens.keys():
            if mutational_burden!=0:
                weights[signature] = generated_signature_burdens[signature][i]/mutational_burden
            else:
                weights[signature] = 0
            generated_weights.loc[i,signature] = weights[signature]
            generated_mutations[i] += reference_signatures[signature]*weights[signature]

        # multiply by mutational burden
        generated_mutations[i] = mutational_burden*generated_mutations[i]

        if options.add_noise:
            if options.noise_type=="gaussian" or options.noise_type=="Gaussian" or options.noise_type=="normal" or options.noise_type=="Normal":
                noise_term = np.random.normal(0, options.noise_sigma, len(generated_mutations[i]))
                generated_mutations[i] += noise_term
            elif options.noise_type=="poisson" or options.noise_type=="Poisson":
                for category_index in range(len(generated_mutations[i])):
                    generated_mutations[i][category_index] = np.random.poisson(generated_mutations[i][category_index])
            elif options.noise_type=="negative_binomial" or options.noise_type=="Negative_binomial":
                noise_term = np.random.negative_binomial(2, 0.5, len(reference_signatures[signature].values))
                generated_mutations[i] += noise_term
            # make sure there are no negative mutation counts
            generated_mutations.loc[generated_mutations[i]<0, i] = 0

        # rounding and converting to to integer counts
        generated_mutations[i] = round(generated_mutations[i],0)
        generated_mutations[i] = generated_mutations[i].astype(int)

    # plot total mutational burden distribution:
    plot_mutational_burden(generated_mutational_burdens, savepath='%s/%s_plots/generated_burden_%s_total.pdf' % (output_path, dataset_name, mutation_type) )

    # save dataframes
    generated_mutations.to_csv(output_filename)
    generated_weights.to_csv(output_filename.replace('.csv','.weights.csv'))
