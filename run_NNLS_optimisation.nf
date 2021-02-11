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

log.info ""
log.info "--------------------------------------------------------"
log.info "      NEXTFLOW NNLS OPTIMISATION ANALYSIS v1.1          "
log.info "--------------------------------------------------------"
log.info "Copyright (C) Sergey Senkin"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run run_NNLS_optimisation.nf"
    log.info ""
    log.info "Nextflow currently does not support list parameters,"
    log.info "so please specify the parameters directly in the script."
    log.info ""
    exit 0
} else {
/* Software information */
log.info "help:                               ${params.help}"
}

params.mutation_type = 'SBS'
params.dataset = 'SIM_ESCC'
params.input_tables = "$baseDir/input_mutation_tables"
params.signature_tables = "$baseDir/signature_tables"
params.NNLS_output_path = "$baseDir/output_opt_check"
params.plots_output_path = "$baseDir/plots_opt_check"
params.signature_prefix = "sigProfilerESCC"
params.SBS_context = 96

// bootstrap flag and method (binomial, multinomial, residuals, classic, bootstrap_residuals)
params.perform_bootstrapping = true
params.bootstrap_method = "binomial"
params.number_of_bootstrapped_samples = 100

params.weak_thresholds = ['0.0000', '0.0100', '0.0200', '0.0300', '0.0400', '0.0500', '0.0600', '0.0700', '0.0800', '0.0900']
params.strong_threshold = '0.0000'
params.optimisation_strategy = "removal" // optimisation strategy (removal, addition or add-remove)

params.signature_attribution_thresholds = 0..20

all_thresholds = params.weak_thresholds

process run_NNLS {
  publishDir "${params.NNLS_output_path}"

  input:
  val weak_threshold from all_thresholds

  output:
  file("${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/*.csv") into all_outputs

  script:
  """
  mkdir -p ${params.NNLS_output_path}/${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}
  [[ ${params.dataset} == *"SIM"* ]] && [[ ${params.mutation_type} == "SBS" ]] && \
    cp ${params.input_tables}/${params.dataset}/WGS_${params.dataset}.${params.SBS_context}.weights.csv ${params.NNLS_output_path}/${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/
  [[ ${params.dataset} == *"SIM"* ]] && [[ ${params.mutation_type} == "DBS" ]] && \
    cp ${params.input_tables}/${params.dataset}/WGS_${params.dataset}.dinucs.weights.csv ${params.NNLS_output_path}/${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/
  [[ ${params.dataset} == *"SIM"* ]] && [[ ${params.mutation_type} == "ID" ]] && \
    cp ${params.input_tables}/${params.dataset}/WGS_${params.dataset}.indels.weights.csv ${params.NNLS_output_path}/${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/
  python $baseDir/bin/run_NNLS.py -d ${params.dataset} -t ${params.mutation_type} -p ${params.signature_prefix} \
                                  --optimisation_strategy ${params.optimisation_strategy} \
                                  -W ${weak_threshold} -S ${params.strong_threshold} \
                                  -i ${params.input_tables} -s ${params.signature_tables} \
                                  -o "./" -x -c ${params.SBS_context} --add_suffix
  """
}

process run_NNLS_bootstrapping {
  publishDir "${params.NNLS_output_path}"

  input:
  each i from 1..params.number_of_bootstrapped_samples
  val weak_threshold from all_thresholds

  output:
  file("./${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/bootstrap_output/output_${params.dataset}_${params.mutation_type}_${weak_threshold}_${params.strong_threshold}_${i}_mutations_table.csv")
  file("./${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/bootstrap_output/output_${params.dataset}_${params.mutation_type}_${weak_threshold}_${params.strong_threshold}_${i}_stat_info.csv")
  file("./${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/bootstrap_output/output_${params.dataset}_${params.mutation_type}_${weak_threshold}_${params.strong_threshold}_${i}_weights_table.csv") into bootstrap_output_tables

  when:
  params.perform_bootstrapping

  script:
  """
  python $baseDir/bin/run_NNLS.py -B -d ${params.dataset} -t ${params.mutation_type} -c ${params.SBS_context} -x \
                                  --optimisation_strategy ${params.optimisation_strategy} \
                                  --bootstrap_method ${params.bootstrap_method} \
                                  -W ${weak_threshold} -S ${params.strong_threshold} --add_suffix \
                                  -p ${params.signature_prefix} -i ${params.input_tables} -s ${params.signature_tables} -o "./"
  mkdir -p ${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/bootstrap_output
  mv ${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/output_${params.dataset}_${params.mutation_type}_mutations_table.csv ${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/bootstrap_output/output_${params.dataset}_${params.mutation_type}_${weak_threshold}_${params.strong_threshold}_${i}_mutations_table.csv
  mv ${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/output_${params.dataset}_${params.mutation_type}_weights_table.csv ${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/bootstrap_output/output_${params.dataset}_${params.mutation_type}_${weak_threshold}_${params.strong_threshold}_${i}_weights_table.csv
  mv ${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/output_${params.dataset}_${params.mutation_type}_stat_info.csv ${params.dataset}_${params.SBS_context}_NNLS_${weak_threshold}_${params.strong_threshold}/bootstrap_output/output_${params.dataset}_${params.mutation_type}_${weak_threshold}_${params.strong_threshold}_${i}_stat_info.csv
  """
}

process make_bootstrap_tables {
  publishDir "${params.NNLS_output_path}"

  input:
  file bootstrap_weights from bootstrap_output_tables.collect()
  val weak_threshold from all_thresholds

  output:
  file("*/*.csv") into all_bootstrap_outputs
  file '*/truth_studies/*.csv' optional true
  file '*/truth_studies/*.json' optional true

  when:
  params.perform_bootstrapping

  script:
  """
  echo 1
  python $baseDir/bin/make_bootstrap_tables.py -d ${params.dataset} -t ${params.mutation_type} -p ${params.signature_prefix} \
          --suffix ${weak_threshold}_${params.strong_threshold} \
          -c ${params.SBS_context} -S ${params.signature_tables} \
          -T ${params.signature_attribution_thresholds.join(' ')} \
          -i ${params.NNLS_output_path} -o "./" -n ${params.number_of_bootstrapped_samples}
  """
}

process plot_heatmaps {
  publishDir "${params.plots_output_path}"

  input:
  file ('*.csv') from all_bootstrap_outputs.collect()

  output:
  file '*/*.pdf' optional true
  file '*/*/*.pdf' optional true
  file '*/*/*/*.pdf' optional true

  script:
  """
  python $baseDir/bin/plot_metric_heatmaps.py -d ${params.dataset} -t ${params.mutation_type} \
                                                  -I ${params.input_tables}/${params.dataset} \
                                                  -i ${params.NNLS_output_path} -o "./" \
                                                  -c ${params.SBS_context} \
                                                  -W ${params.weak_thresholds.join(' ')} \
                                                  -S ${params.strong_threshold} \
                                                  -T ${params.signature_attribution_thresholds.join(' ')} \
                                                  --signature_path ${params.signature_tables} \
                                                  -p ${params.signature_prefix}
  """
}
