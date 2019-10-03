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
log.info "      NEXTFLOW NNLS OPTIMISATION ANALYSIS v1.0       "
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

params.mutation_types = ['SBS']
params.input_tables = "$PWD/input_mutation_tables"
params.input_signatures = "$PWD/signature_tables"
params.NNLS_output_path = "$PWD/output_opt_check"
params.signature_prefix = "sigProfiler"

params.weak_thresholds = ['0.0010', '0.0020', '0.0030', '0.0040', '0.0050', '0.0060', '0.0070', '0.0080', '0.0090', '0.0100']
params.strong_thresholds = ['0.0010', '0.0020', '0.0030', '0.0040', '0.0050', '0.0060', '0.0070', '0.0080', '0.0090', '0.0100']

Channel.from (params.weak_thresholds).combine( params.strong_thresholds ).set { all_thresholds }
all_thresholds.into { simple_context_thresholds ; TSB_context_thresholds }

process run_NNLS {
  publishDir "${params.NNLS_output_path}"

  input:
  set weak_threshold, strong_threshold from simple_context_thresholds

  output:
  file("SIM_96_NNLS_${weak_threshold}_${strong_threshold}/*.csv")

  script:
  """
  python $PWD/scripts/run_NNLS.py -d SIM -t SBS -W ${weak_threshold} -p ${params.signature_prefix} -S ${strong_threshold} -i ${params.input_tables} -s ${params.input_signatures} -o "./" -x --add_suffix
  """
}

process run_NNLS_TSB {
  publishDir "${params.NNLS_output_path}"

  input:
  set weak_threshold, strong_threshold from TSB_context_thresholds

  output:
  file("SIM_192_NNLS_${weak_threshold}_${strong_threshold}/*.csv")

  script:
  """
  python $PWD/scripts/run_NNLS.py -d SIM -t SBS -W ${weak_threshold} -p ${params.signature_prefix} -S ${strong_threshold} -i ${params.input_tables} -s ${params.input_signatures} -o "./" -x -c 192 --add_suffix
  """
}
