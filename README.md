# Mutational Signature Analysis

Mutational signature attribution analysis, including code used for optimisation study with simulated data.

## Getting started

Retrieve the code:
```
git clone git@gitlab.com:s.senkin/MSA.git
cd MSA
```


## Running with Nextflow
The easiest way to run the code is by using [Nextflow](https://www.nextflow.io/).
Once you have installed Nextflow, run it locally or on your favourite cluster:

```
nextflow run run_analysis.nf
```

The pipeline should run everything and produce all the results automatically.

## Running manually

### Setting up dependencies

Dependencies so far are: *pandas*, *numpy*, *scipy* and *matplotlib*. If you don't have them, the easiest way is to set up the virtual environment using [conda](https://conda.io) package manager:

```
conda create -n MSA
source activate MSA
conda install pandas numpy scipy matplotlib
```

This only needs to be done once. Afterwards, just activate the environment whenever needed:

```
source activate MSA
```


### Simulating data

* [simulations](simulations) folder contains a set of 100 simulated samples for 96/192 contexts SBS, as well as dinucs and indels, where each sample contains contributions from **5** randomly selected signatures out of **100** Poisson-generated signatures. To reproduce (e.g. for 96-context SBS, 100 signatures and samples), run:
```
python simulation_code/generate_random_signatures.py -t SBS -c 96 -n 100
python simulation_code/simulate_data.py -r -t SBS -c 96 -n 100
```
In both scripts, a normal distribution can be used to generate white noise using *-z* option, with a Gaussian centred around **0** for each category of mutations, with standard deviation set by *-Z* option (**2** by default). Additional flags can be viewed in the code or using *-h* option in each script.

* Alternatively, [np.random.normal](https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.normal.html) function can be used to generate normal distributions of mutational burdens corresponding to each PCAWG signature mentioned in the [signatures_to_generate](https://gitlab.com/s.senkin/signature_attribution_algorithms/blob/simulations/simulate_PCAWG_signatures.py#L8) dictionary in the script, containing Gaussian means and standard deviations for each signature.
* Note that the distributions are not strictly Gaussian since negative numbers of burdens are replaced by zeros
* To reproduce the simulated set of samples with reshuffled *SBS1/5/22/40* PCAWG signatures, one can run the above-mentioned script without the *-r* option:
```
python simulation_code/simulate_data.py -t SBS -c 96 -n 100 -s signatures_PCAWG
```

The file format produced is the same as that of the existing PCAWG dataset.

### Running NNLS

To run NNLS on simulated data:
```
./run_all_NNLS.sh
```

To run NNLS on your own $dataset, create your own $dataset_folder with mutation tables:
```
./run_all_NNLS.sh -d $dataset -i $dataset_folder -s signatures_PCAWG
```

To limit the running to the first **n** samples, use the *-n* option.

To manually run NNLS (*-t* to set mutation type, *-x* flag for optimised method):

```
python NNLS/run_NNLS.py -t SBS -x
```

To change the context, use the *-c* flag (only relevant for SBS):

```
python NNLS/run_NNLS.py -t SBS -c 192 -x
```

## Similarity and attribution efficiency measurements

### Optimisation thresholds parameter space scan

[Nextflow](https://www.nextflow.io/) scripts have been implemented to run the parameter space scan of the optimisation thresholds for NNLS routine. This can be run locally or on your favourite cluster.

```
nextflow run run_NNLS_optimisation.nf
```

The output will be produced in **output_opt_check** folder, upon which one may run the plotting script producing efficiency heatmaps, e.g.:

```
python scripts/make_efficiency_heatmaps.py -m NNLS -c 96
python scripts/make_efficiency_heatmaps.py -m NNLS -c 192
```

### Fixed attribution results comparison

Upon running signature attribution, the script to measure and plot attribution efficiencies for various methods can be run as follows (*-t* flag to choose other mutations types: *DBS* or *ID*, help on more flags with *-h*):
```
python scripts/measure_attribution_efficiency.py -t SBS
```

All the efficiency plots are be produced in **efficiency_plots** folder.
