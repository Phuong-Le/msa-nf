#!/usr/bin/env nextflow

// Copyright (C) 2018 Sergey Senkin

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null

log.info ''
log.info '--------------------------------------------------------'
log.info '              __   __  _____                             '
log.info '             |  \\/  |/ ____|  /\\                      '
log.info '             | \\  / | (___   /  \\                     '
log.info '             | |\\/| |\\___ \\ / /\\ \\                 '
log.info '             | |  | |____) / ____ \\                    '
log.info '             |_|  |_|_____/_/    \\_\\                  '
log.info '                                                        '
log.info '          MUTATIONAL SIGNATURE ANALYSIS v1.1            '
log.info '--------------------------------------------------------'
log.info 'Copyright (C) Sergey Senkin'
log.info 'This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE'
log.info 'This is free software, and you are welcome to redistribute it'
log.info 'under certain conditions; see LICENSE for details.'
log.info '--------------------------------------------------------'
log.info ''

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run run_analysis.nf"
    log.info ""
    log.info "Nextflow currently does not support list parameters,"
    log.info "so please specify the parameters directly in the script."
    log.info ""
    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
}

// input data
params.datasets = ['SIM_test']
params.mutation_types = ['SBS', 'DBS', 'ID'] // add or remove SBS/DBS/ID if needed
params.input_tables = "$baseDir/input_mutation_tables"
params.SBS_context = 96 // 96, 192, and 288 context matrices can be provided (SBS only)
params.number_of_samples = -1 // number of samples to analyse (-1 means all available)

// output paths
params.NNLS_output_path = "$baseDir/output_tables" //_" + params.weak_threshold + "_" + params.strong_threshold
params.plots_output_path = "$baseDir/plots" //_" + params.weak_threshold + "_" + params.strong_threshold

// signatures to use
params.signature_tables = "$baseDir/signature_tables"
params.signature_prefix = "sigRandom" // prefix of signature files to use (e.g. sigProfiler, sigRandom)

// optimisation flag and parameters
params.optimised = false
params.optimisation_strategy = "removal" // optimisation strategy (removal, addition or add-remove)
params.weak_threshold = 0.02
params.strong_threshold = 0.02

// plotting flags
params.plot_signatures = true
params.plot_input_spectra = true
params.plot_fitted_spectra = true
params.plot_residuals = true
params.show_poisson_errors = true
params.show_strands = false // only works with higher contexts (192, 288)
params.show_nontranscribed_region = false // only wortks with higher contexts (288)

// bootstrap flag and method (binomial, multinomial, residuals, classic, bootstrap_residuals)
params.perform_bootstrapping = true
params.bootstrap_method = "binomial"
params.number_of_bootstrapped_samples = 10 // at least 100 is recommended
params.use_absolute_attributions = false // use absolute mutation counts in bootstrap analysis (relative by default)

// if SIM in dataset name (synthetic data), use the following percentage range for measuring signature attirbution sensitivities
params.signature_attribution_thresholds = 0..20

// helper flags for scripts (automatic based on parameters)
optimised_flag = (params.optimised) ? "-x" : ''
abs_flag = (params.use_absolute_attributions) ? "-a" : ''
suffix = (params.use_absolute_attributions) ? "abs_mutations" : 'weights'
error_flag = (params.show_poisson_errors) ? "-e" : ''
strands_flag = (params.show_strands) ? "-b" : ''
nontranscribed_flag = (params.show_nontranscribed_region) ? "-n" : ''

// add parameter values to log output (.nextflow.log)
log.info params.collect { k,v -> "${k.padRight(34)}: $v" }.join("\n")

process plot_input_spectra {
  publishDir "${params.plots_output_path}"

  input:
  each mutation_type from params.mutation_types
  each dataset from params.datasets

  output:
  file '*/*/*.pdf' optional true
  file '*/*/*/*.pdf' optional true
  file '*/*/*/*/*.pdf' optional true
  file '*/*/*/*/*/*.pdf' optional true

  when:
  params.plot_input_spectra

  script:
  """
  python $baseDir/scripts/plot_mutation_spectra.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -i ${params.input_tables} ${strands_flag} ${nontranscribed_flag} -o "./"
  python $baseDir/scripts/plot_mutation_spectra.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -r -i ${params.input_tables} ${strands_flag} ${nontranscribed_flag} -o "./"
  """
}

process plot_signatures {
  publishDir "${params.plots_output_path}"

  input:
  each mutation_type from params.mutation_types
  each dataset from params.datasets

  output:
  file '*/*/*.pdf' optional true
  file '*/*/*/*.pdf' optional true
  file '*/*/*/*/*.pdf' optional true
  file '*/*/*/*/*/*.pdf' optional true

  when:
  params.plot_signatures

  script:
  """
  python $baseDir/scripts/plot_mutation_spectra.py -S -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} \
                                               -p ${params.signature_prefix} -s ${params.signature_tables} \
                                               -r ${strands_flag} ${nontranscribed_flag} -o "./"
  """
}

process run_NNLS_normal {
  publishDir "${params.NNLS_output_path}"

  input:
  each mutation_type from params.mutation_types
  each dataset from params.datasets

  output:
  file("./${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv")
  file("./${dataset}/output_${dataset}_${mutation_type}_weights_table.csv")
  file("./${dataset}/output_${dataset}_${mutation_type}_stat_info.csv")
  file("./${dataset}/output_${dataset}_${mutation_type}_fitted_values.csv")
  file("./${dataset}/output_${dataset}_${mutation_type}_residuals.csv") into central_NNLS_residuals
  set dataset, mutation_type into attribution_for_bootstrap_plots
  set dataset, mutation_type into attribution_for_spectra_plots
  set dataset, mutation_type into attribution_for_residuals
  set dataset, mutation_type into attribution_for_metrics
  set dataset, mutation_type into attribution_for_tables

  script:
  """
  mkdir -p ${params.NNLS_output_path}/${dataset}
  [[ ${dataset} == *"SIM"* ]] && [[ ${mutation_type} == "SBS" ]] && \
    cp ${params.input_tables}/${dataset}/WGS_${dataset}.${params.SBS_context}.weights.csv ${params.NNLS_output_path}/${dataset}/
  [[ ${dataset} == *"SIM"* ]] && [[ ${mutation_type} == "DBS" ]] && \
    cp ${params.input_tables}/${dataset}/WGS_${dataset}.dinucs.weights.csv ${params.NNLS_output_path}/${dataset}/
  [[ ${dataset} == *"SIM"* ]] && [[ ${mutation_type} == "ID" ]] && \
    cp ${params.input_tables}/${dataset}/WGS_${dataset}.indels.weights.csv ${params.NNLS_output_path}/${dataset}/
  python $baseDir/scripts/run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \
                                  --optimisation_strategy ${params.optimisation_strategy} \
                                  -W ${params.weak_threshold} -S ${params.strong_threshold} -n ${params.number_of_samples} \
                                  -p ${params.signature_prefix} -i ${params.input_tables} -s ${params.signature_tables} -o "./"
  cp ${dataset}/output_${dataset}_${mutation_type}_residuals.csv ${params.input_tables}/${dataset}/
  cp ${dataset}/output_${dataset}_${mutation_type}_fitted_values.csv ${params.input_tables}/${dataset}/
  """
}


process plot_fitted_spectra {
  publishDir "${params.plots_output_path}"

  input:
  set dataset, mutation_type from attribution_for_spectra_plots

  output:
  file '*/*/*.pdf' optional true
  file '*/*/*/*.pdf' optional true
  file '*/*/*/*/*.pdf' optional true
  file '*/*/*/*/*/*.pdf' optional true

  when:
  params.plot_fitted_spectra

  script:
  """
  python $baseDir/scripts/plot_mutation_spectra.py -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -i ${params.NNLS_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  python $baseDir/scripts/plot_mutation_spectra.py -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -r -i ${params.NNLS_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  python $baseDir/scripts/plot_mutation_spectra.py -C -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -i ${params.NNLS_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  python $baseDir/scripts/plot_mutation_spectra.py -C -f -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -r -i ${params.NNLS_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  """
}

process plot_residuals {
  publishDir "${params.plots_output_path}"

  input:
  set dataset, mutation_type from attribution_for_residuals

  output:
  file '*/*/*.pdf' optional true
  file '*/*/*/*.pdf' optional true
  file '*/*/*/*/*.pdf' optional true
  file '*/*/*/*/*/*.pdf' optional true

  when:
  params.plot_residuals

  script:
  """
  python $baseDir/scripts/plot_mutation_spectra.py -H -R -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -i ${params.NNLS_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  python $baseDir/scripts/plot_mutation_spectra.py -C -R -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} --number ${params.number_of_samples} \
                                                -i ${params.NNLS_output_path} ${error_flag} ${strands_flag} ${nontranscribed_flag} -o "./"
  """
}

// plot errors from bootstrap, add stat.info on plots


process run_NNLS_bootstrapping {
  publishDir "${params.NNLS_output_path}"

  input:
  each i from 1..params.number_of_bootstrapped_samples
  each mutation_type from params.mutation_types
  each dataset from params.datasets
  // file residuals from central_NNLS_residuals // uncomment in using residuals bootstrapping

  output:
  file("./${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_mutations_table.csv")
  file("./${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_stat_info.csv")
  file("./${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_weights_table.csv") into bootstrap_output_tables

  when:
  params.perform_bootstrapping

  script:
  """
  python $baseDir/scripts/run_NNLS.py -B -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \
                                  --optimisation_strategy ${params.optimisation_strategy} --bootstrap_method ${params.bootstrap_method} \
                                  -W ${params.weak_threshold} -S ${params.strong_threshold} -n ${params.number_of_samples} \
                                  -p ${params.signature_prefix} -i ${params.input_tables} -s ${params.signature_tables} -o "./"
  mkdir -p ${dataset}/bootstrap_output
  mv ${dataset}/output_${dataset}_${mutation_type}_mutations_table.csv ${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_mutations_table.csv
  mv ${dataset}/output_${dataset}_${mutation_type}_weights_table.csv ${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_weights_table.csv
  mv ${dataset}/output_${dataset}_${mutation_type}_stat_info.csv ${dataset}/bootstrap_output/output_${dataset}_${mutation_type}_${i}_stat_info.csv
  """
}

process make_bootstrap_tables {
  publishDir "${params.NNLS_output_path}"

  input:
  set dataset, mutation_type from attribution_for_tables
  file bootstrap_weights from bootstrap_output_tables.collect()

  output:
  file("./${dataset}/signatures_prevalences_${mutation_type}.csv") into signature_prevalences
  file("./${dataset}/attributions_per_sample_${mutation_type}_bootstrap_output_${suffix}.json") into attributions_per_sample
  file '*/*.csv'
  file '*/*.json'
  file '*/truth_studies/*.csv' optional true
  file '*/truth_studies/*.json' optional true

  when:
  params.perform_bootstrapping

  script:
  """
  python $baseDir/scripts/make_bootstrap_tables.py -d ${dataset} -t ${mutation_type} -p ${params.signature_prefix} ${abs_flag} \
                                               -c ${params.SBS_context} -S ${params.signature_tables} \
                                               -T ${params.signature_attribution_thresholds.join(' ')} \
                                               -i ${params.NNLS_output_path} -o "./" -n ${params.number_of_bootstrapped_samples}
  """
}

process plot_bootstrap_attributions {
  publishDir "${params.plots_output_path}"

  input:
  set dataset, mutation_type from attribution_for_bootstrap_plots
  file bootstrap_attributions from attributions_per_sample.collect()

  output:
  file '*/*/bootstrap_plots/*.pdf' optional true
  file '*/*/bootstrap_plots/*/*.pdf' optional true
  file '*/*/bootstrap_plots/*/*/*.pdf' optional true

  when:
  params.perform_bootstrapping

  script:
  """
  python $baseDir/scripts/plot_bootstrap_attributions.py -d ${dataset} -t ${mutation_type} -p ${params.signature_prefix} ${abs_flag} \
                                                     -c ${params.SBS_context} -S ${params.signature_tables} -I ${params.input_tables} \
                                                     -i ${params.NNLS_output_path} -o "./" -n ${params.number_of_bootstrapped_samples}
  """
}


process plot_metrics {
  publishDir "${params.plots_output_path}"

  input:
  set dataset, mutation_type from attribution_for_metrics
  file prevalences from signature_prevalences.collect()

  output:
  file '*/*/bootstrap_plots/*.pdf' optional true
  file '*/*/bootstrap_plots/*/*.pdf' optional true
  file '*/*/bootstrap_plots/*/*/*.pdf' optional true

  when:
  params.perform_bootstrapping

  script:
  """
  python $baseDir/scripts/plot_metrics.py -d ${dataset} -t ${mutation_type} -i ${params.NNLS_output_path} -o "./"
  """
}
