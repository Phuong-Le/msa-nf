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
log.info "      NEXTFLOW MUTATIONAL SIGNATURE ANALYSIS v1.0       "
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

params.datasets = ['SIM']
params.mutation_types = ['SBS']
params.input_tables = "$PWD/input_mutation_tables"
params.signature_tables = "$PWD/signature_tables"
params.SBS_context = 192
params.optimised = true
params.perform_bootstrapping = true
params.number_of_bootstrapped_samples = 1000
params.NNLS_output_path = "$PWD/output_tables"
params.plots_output_path = "$PWD/plots"

optimised_flag = (params.optimised) ? "-x" : ''

process run_NNLS_normal {
  publishDir "${params.NNLS_output_path}"

  input:
  each mutation_type from params.mutation_types
  each dataset from params.datasets

  output:
  file("./${dataset}/output_${mutation_type}_mutations_table.csv")
  file("./${dataset}/output_${mutation_type}_weights_table.csv")
  file("./${dataset}/output_${mutation_type}_stat_info.csv")
  set dataset, mutation_type into attribution_for_bootstrap_plots

  script:
  """
  python $PWD/scripts/run_NNLS.py -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \
                                  -i ${params.input_tables} -s ${params.signature_tables} -o "./"
  """
}

process run_NNLS_bootstrapping {
  publishDir "${params.NNLS_output_path}"

  input:
  each i from 1..params.number_of_bootstrapped_samples
  each mutation_type from params.mutation_types
  each dataset from params.datasets

  output:
  file("./${dataset}/output_${dataset}_${mutation_type}_${i}_mutations_table.csv")
  file("./${dataset}/output_${dataset}_${mutation_type}_${i}_stat_info.csv")
  file("./${dataset}/output_${dataset}_${mutation_type}_${i}_weights_table.csv") into bootstrap_output

  when:
  params.perform_bootstrapping

  script:
  """
  python $PWD/scripts/run_NNLS.py -B -d ${dataset} -t ${mutation_type} -c ${params.SBS_context} ${optimised_flag} \
                                  -i ${params.input_tables} -s ${params.signature_tables} -o "./"
  mv ${dataset}/output_${mutation_type}_mutations_table.csv ${dataset}/output_${dataset}_${mutation_type}_${i}_mutations_table.csv
  mv ${dataset}/output_${mutation_type}_weights_table.csv ${dataset}/output_${dataset}_${mutation_type}_${i}_weights_table.csv
  mv ${dataset}/output_${mutation_type}_stat_info.csv ${dataset}/output_${dataset}_${mutation_type}_${i}_stat_info.csv
  """
}

process plot_bootstrap_attributions {
  publishDir "${params.plots_output_path}"

  input:
  set dataset, mutation_type from attribution_for_bootstrap_plots
  file bootstrap_weights from bootstrap_output.collect()

  output:
  file '*/*/bootstrap_plots/*.pdf' optional true

  when:
  params.perform_bootstrapping

  script:
  """
  mkdir -p ${params.NNLS_output_path}/${dataset}
  [[ ${dataset} == *"SIM"* ]] && [[ ${mutation_type} == "SBS" ]] && \
    cp ${params.input_tables}/${dataset}/WGS_${dataset}.${params.SBS_context}.weights.csv ${params.NNLS_output_path}/${dataset}/
  [[ ${dataset} == *"SIM"* ]] && [[ ${mutation_type} == "DBS" ]] && \
    cp ${params.input_tables}/${dataset}/WGS_${dataset}.dinucs.weights.csv ${params.NNLS_output_path}/${dataset}/
  [[ ${dataset} == *"SIM"* ]] && [[ ${mutation_type} == "ID" ]] && \
    cp ${params.input_tables}/${dataset}/WGS_${dataset}.indels.weights.csv ${params.NNLS_output_path}/${dataset}/
  python $PWD/scripts/plot_bootstrap_attributions.py -d ${dataset} -t ${mutation_type} \
          -c ${params.SBS_context} -S ${params.signature_tables} -I ${params.input_tables} \
          -i ${params.NNLS_output_path} -o "./" -n ${params.number_of_bootstrapped_samples}
  """
}
