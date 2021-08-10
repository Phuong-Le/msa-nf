# Mutational Signature Attribution pipeline

![logo](MSA.png)

Mutational signature attribution analysis, featuring automised optimisation study with simulated data.

## Running with Nextflow
The best way to run the code is by using [Nextflow](https://www.nextflow.io/).
Once you have installed Nextflow, run the test job locally or on your favourite cluster:

```
nextflow run https://gitlab.com/s.senkin/MSA -profile docker
```

If you don't have [docker](https://www.docker.com/) installed, you can also use [conda](https://conda.io) or [singularity](https://sylabs.io/singularity/) profiles.
The pipeline should run everything and produce all the results automatically. You can retrieve the code ([see below](https://gitlab.com/s.senkin/MSA#getting-started)) in order to adjust the inputs and parameters. In the [run_auto_optimised_analysis.nf](run_auto_optimised_analysis.nf) file various parameters can be specified.

## Running on SigProfiler output

MSA natively supports [SigProfilerExtractor](https://github.com/AlexandrovLab/SigProfilerExtractor) and [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator) outputs.

```
nextflow run https://gitlab.com/s.senkin/MSA -profile docker --dataset SP_test --SP_matrix_generator_output_path path/to/SP_ME/ --SP_extractor_output_path path/to/SP/
```

## Options

### General parameters

| Parameters  | Default value | Description |
|-----------|-------------|-------------|
| --help | null | print usage and optional parameters |
| --SP_matrix_generator_output_path | null | optionally use SigProfilerMatrixGenerator output from specified path |
| --SP_extractor_output_path | null | optionally use SigProfilerExtractor output from specified path to attribute signatures extracted by SigProfiler |
| --dataset | SIM_test | set the name of the dataset. If no SigProfiler output is provided, the matrices must exist in params.input_tables folder (see example) |
| --input_tables | $baseDir/input_mutation_tables | if not using SigProfiler outputs, location of input mutation tables, the repository one is used by default |
| --signature_tables | $baseDir/signature_tables | if not using SigProfiler outputs, location of input signature tables, the repository one is used by default |
| --signature_prefix | sigProfiler |  if not using SigProfiler outputs, prefix of signature files to use, must be located in signature_tables folder (e.g. sigProfiler, sigRandom) |
| --output_path | . | output path for plots and tables |
| --mutation_types | ['SBS', 'DBS', 'ID'] | mutation types to analyse. Only one can be specified from command line, or a list in the run_auto_optimised_analysis.nf file |
| --number_of_samples | -1 | number of samples to analyse (-1 means all available) |
| --SBS_context | 96 | SBS context to use (96, 192 or 288) |
| --COSMIC_signatures | false | if set to true, COSMIC signatures are used form SigProfiler output, otherwise de-novo ones are used |

### Attribution and optimisation options

| Parameters  | Default value | Description |
|-----------|-------------|-------------|
| --run_only_optimisation | false | set to true if only optimisation results are required, without final attributions |
| --optimisation_plots_output_path | params.plots_output_path + "/optimisation_plots" | directory to store optimisation plots |
| --number_of_bootstrapped_samples | 100 | number of bootstrap samples variations (at least 100 is recommended) |
| --optimised | true | perform signature optimisation for NNLS attribution method; if set to false, optimisation will run but not be used in final attributions |
| --optimisation_strategy | removal | parameter defining the strategy: 'removal' (default), 'addition' or 'add-remove', determining the method of executing signature addition and/or removal loops |
| --weak_thresholds | ['0.0000', '0.0100', '0.0200', '0.0300', '0.0400', '0.0500'] | range of L2 similarity decrease thresholds to be scanned, excluding weakest signatures |
| --strong_thresholds | ['0.0000'] | range of L2 similarity increase thresholds to be scanned, including strongest signatures: only one is sufficient in default removal strategy |
| --bootstrap_method | binomial | method of parametric bootstrap (binomial, multinomial, residuals, classic, bootstrap_residuals) |
| --number_of_bootstrapped_samples_in_optimisation | 100 | number of bootstrap variations used in optimisation steps |
| --number_of_bootstrapped_samples | 1000 | number of bootstrap variations in final attribution, at least 100 is recommended |
| --confidence_level | 95 | confidence level for CI calculation, in percent |
| --use_absolute_attributions | false | use absolute mutation counts in final bootstrap outputs (relative by default) |
| --metric_to_prioritise | "specificity" | metric to prioritise (default: specificity), requiring at least the specified threshold or closest alternative |
| --metric_threshold | 0.95 | minimum threshold of the prioritised metric |
| --signatures_to_prioritise | [] | list of signatures to prioritise (empty list means all, by default) |
| --no_CI_for_penalties | false | do not use confidence intervals for optimal penalties calculation, rather do standard count experiment |
| --calculate_penalty_on_average | false | apply criteria based on signatures overall (on average, less conservative), rather than maximising prioritised metric for every signature (default, more conservative) |

### Simulation options
| Parameters  | Default value | Description |
|-----------|-------------|-------------|
| --run_only_simulations | false | set to true if only simulations are needed, which will be produced in $baseDir/output_tables folder |
| --number_of_simulated_samples | 1000 | number of simulated samples in data-driven model, at least 1000 is recommended |
| --add_noise | true | add noise in simulations (recommended) |
| --noise_type | "gaussian" | type of noise in simulations: gaussian, poisson or negative_binomial (Gaussian by default) |
| --noise_stdev | 10 | standard deviation of gaussian noise, in percentage of sample mutation burden (10 percent by default) |

### Plotting options

| Parameters  | Default value | Description |
|-----------|-------------|-------------|
| --plot_optimisation_plots | true | make plots showing various optimisation metrics, showing justification for optimal penalties |
| --plot_signatures | true | plot provided signatures (SBS, DBS or ID supported) |
| --plot_input_spectra | true | plot mutation spectra for provided samples |
| --plot_fitted_spectra | false | plot fitted mutation spectra using NNLS output |
| --plot_residuals | false | plot residuals spectra (fitted-input) |
| --show_poisson_errors | false | show Poisson errors in spectra plots |
| --show_strands | false | show different strands, only works for SBS with higher contexts (192, 288) |
| --show_nontranscribed_region | false | show non-transcibed region, only works for 288 context in SBS |

## Running manually

### Getting started

Retrieve the code:
```
git clone https://gitlab.com/s.senkin/MSA.git
cd MSA
```

Run the fully automised pipeline with optimisation on the test dataset (or adjust parameters/input matrices accordingly)
```
nextflow run run_auto_optimised_analysis.nf --dataset SIM_test --input_tables ./input_tables --output_path ./test
```

Alternatively, to run the pipeline without optimisation, using fixed penalties on the test dataset:
```
nextflow run run_analysis.nf --dataset SIM_test --weak_threshold 0.02 --input_tables ./input_tables --output_path ./test
```

### Setting up dependencies

If you can not run *nextflow*, you can still run some basic analysis manually (scripts in the *./bin* folder).
Dependencies so far are: *pandas*, *numpy*, *scipy*, *matplotlib* and *seaborn*. If you don't have them, the easiest way is to set up the virtual environment using [conda](https://conda.io) package manager:

```
conda env create -f environment.yml
```

This only needs to be done once. Afterwards, just activate the environment whenever needed:

```
source activate msa
```

Alternatively, you can use *docker* yourself with the *Dockerfile* provided, or use ready-made images ([docker](https://hub.docker.com/r/ssenkin/msa/tags) or [singularity](https://cloud.sylabs.io/library/ssenkin/default/msa)).


### Simulating data (if not using automised pipeline)

* [input_mutation_tables/SIM](input_mutation_tables/SIM) folder contains a set produced with existing (PCAWG or COSMIC) signatures. [np.random.normal](https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.normal.html) function was used to generate normal distributions of mutational burdens corresponding to each PCAWG signature mentioned in the [signatures_to_generate](bin/simulate_data.py#L9) dictionary in the script, containing Gaussian means and standard deviations for each signature.
* Note that the distributions are not strictly Gaussian since negative numbers of burdens are replaced by zeros
* To reproduce the simulated set of samples with reshuffled *SBS1/5/22/40* PCAWG signatures, one can run the following script (without the *-r* option):
```
python bin/simulate_data.py -t SBS -c 96 -n 100 -s signature_tables
```

* [input_mutation_tables/SIMrand](input_mutation_tables/SIMrand) folder contains a set of 100 simulated samples for 96/192 contexts SBS, as well as dinucs and indels, where each sample contains contributions from **5** randomly selected signatures out of **100** Poisson-generated signatures. To reproduce (e.g. for 96-context SBS, 100 signatures and samples), run:
```
python bin/generate_random_signatures.py -t SBS -c 96 -n 100
python bin/simulate_data.py -r -t SBS -c 96 -n 100 -d SIMrand
```
In both scripts, a normal distribution can be used to generate white noise using *-z* option, with a Gaussian centred around **0** for each category of mutations, with standard deviation set by *-Z* option (**2** by default). Additional flags can be viewed in the code or using *-h* option in each script.

The file format produced is the same as that of the existing PCAWG dataset.

### Running NNLS manually (deprecated, use with care)

To run NNLS on simulated data:
```
./run_all_NNLS.sh
```

To run NNLS on your own $dataset, create your own $dataset_folder with mutation tables:
```
./run_all_NNLS.sh -d $dataset -i $dataset_folder -s signature_tables
```

To limit the running to the first **n** samples, use the *-n* option.

To manually run NNLS (*-t* to set mutation type, *-x* flag for optimised method):

```
python bin/run_NNLS.py -t SBS -x
```

To change the context, use the *-c* flag (only relevant for SBS):

```
python bin/run_NNLS.py -t SBS -c 192 -x
```

Bootstrap option *-B* allows to run this script for a perturbed mutation table, using method specified with *--bootstrap_method* option (see more in the script).
