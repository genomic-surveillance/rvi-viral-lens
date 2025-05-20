#!/usr/bin/env nextflow
// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include {check_generate_consensus_params} from './workflows/GENERATE_CONSENSUS.nf'
include {check_sort_reads_params} from './workflows/SORT_READS_BY_REF.nf'
include {check_classification_report_params} from './workflows/GENERATE_CLASSIFICATION_REPORT.nf'
include { validateParameters; paramsSummaryLog} from 'plugin/nf-schema'

include {SORT_READS_BY_REF} from './workflows/SORT_READS_BY_REF.nf'
include {GENERATE_CONSENSUS} from './workflows/GENERATE_CONSENSUS.nf'
include {SCOV2_SUBTYPING} from './workflows/SCOV2_SUBTYPING.nf'
include {COMPUTE_QC_METRICS} from './workflows/COMPUTE_QC_METRICS.nf'
include {GENERATE_CLASSIFICATION_REPORT} from './workflows/GENERATE_CLASSIFICATION_REPORT.nf'


/*
* ANSI escape codes to color output messages
*/
ANSI_GREEN = "\033[1;32m"
ANSI_RED = "\033[1;31m"
ANSI_RESET = "\033[0m"
ANSI_BOLD = "\033[1m"

log.info """${ANSI_RESET}
  ===========================================
  Viral Lens [v1.0]
  Used parameters:
  -------------------------------------------
  --> general pipeline parameters:
    --use_local_containers     : ${params.use_local_containers}
    --use_registry_containers  : ${params.use_registry_containers}
    --entry_point              : ${params.entry_point}
    --containers_dir           : ${params.containers_dir}
    --outdir                   : ${params.outdir}
    --run_preprocessing         : ${params.run_preprocessing}

  --> PREPROCESSING workflow parameters:
    --run_trimmomatic          : ${params.run_trimmomatic}
    --run_trf                  : ${params.run_trf}
    --run_hrr                  : ${params.run_hrr}

  --> SORT_READS_BY_REF workflow parameters:
    --manifest                   : ${params.manifest}
    --db_path                    : ${params.db_path}
    --db_library_fa_path         : ${params.db_library_fa_path}
    --min_reads_for_taxid        : ${params.min_reads_for_taxid}
    --k2r_max_total_reads_per_fq : ${params.max_total_reads_per_fq}
    --k2r_dump_fq_mem            : ${params.k2r_dump_fq_mem}

  --> GENERATE_CONSENSUS workflow parameters:
    --consensus_mnf            : ${params.consensus_mnf}
    --ivar_min_depth           : ${params.ivar_min_depth}
    --ivar_freq_threshold      : ${params.ivar_freq_threshold}

  --> GENERATE_CLASSIFICATION_REPORT workflow parameters:
    --min_coverage_percent  : ${params.min_coverage_percent}

  --> viral subtyping branching parameters:
    --scv2_keyword             : ${params.scv2_keyword}

  --> resource management:
    --default_error_strategy   : ${params.default_error_strategy}
    --mem_k2r_b0_offset        : ${params.mem_k2r_b0_offset}
    --mem_k2r_b0               : ${params.mem_k2r_b0}
    --mem_k2r_b0_final         : ${params.mem_k2r_b0_final}
    --mem_k2r_b1               : ${params.mem_k2r_b1}
    --mem_k2r_f1               : ${params.mem_k2r_f1}
    --mem_k2r_a2               : ${params.mem_k2r_a2}
    --max_attempts             : ${params.max_attempts}
  ------------------------------------------
  Runtime data:
  -------------------------------------------
  Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
  Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
  Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
  Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
  ------------------------------------------
""".stripIndent()

// Validate input parameters
validateParameters()
// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// Main entry-point workflow
workflow {

    // === 1 - Process input ===
    check_main_params()
    // ==========================
    reads_ch = parse_mnf(params.manifest)

    // ==========================
    // === 2 - Map reads to taxid
    SORT_READS_BY_REF(reads_ch)
    sample_taxid_ch = SORT_READS_BY_REF.out.sample_taxid_ch // tuple (meta, reads)
    sample_pre_report_ch = SORT_READS_BY_REF.out.sample_pre_report_ch

    // === 3 - Generate consensus ==
    GENERATE_CONSENSUS(sample_taxid_ch)

    // === 4 - Compute QC Metrics
    COMPUTE_QC_METRICS(GENERATE_CONSENSUS.out)

    // === 5 - branching output from QC for viral specific subtyping

    // 5.1 - process pre_report files
    // NOTE: if the consensus_gen entry point is removed,
    //           this processing should be moved back to 
    //           SORT_READS_BY_REF workflow
    sample_pre_report_ch
      .filter{it -> (it.size() > 1)} // remove empty pre_reports
      .splitCsv(header: true, sep:"\t")
      .map{it -> 
        id="${it.sample_id}.${it.selected_taxid}"
        tuple(id, it)
      }
      .set{sample_report_ch}

    // raise warning for sample taxids which had empty pre_reports
    sample_pre_report_ch
      .filter{it -> (it.size() <= 1)}
      .view(it -> log.warn("Excluding ${it} as input due to small size ( < 1 byte)"))

    // 5.2 - add report info to out qc metric chanel and branch for SCOV2 subtyping
    COMPUTE_QC_METRICS.out // tuple (meta, bam)
      .map { meta, bam -> tuple(meta.id, meta, bam)}
      .join(sample_report_ch)//, by: 0) // tuple (id, meta, bam, report)
      .map {id, meta, bam, report ->
        meta.putAll(report)
        tuple(meta, bam)
      }
      .branch{ it ->
        scv2_subtyping_workflow_in_ch: it[0].ref_selected.contains("${params.scv2_keyword}")
        no_subtyping_ch: true
      }
      .set {qc_metrics_out_ch}

    // 5.3 - do SCOV2 subtyping
    if (params.do_scov2_subtyping == true){
      qc_metrics_out_ch.scv2_subtyping_workflow_in_ch
        .map {it -> tuple(it[0], it[0].consensus_fa)}
        .set {scov_2_subt_In_ch}
      SCOV2_SUBTYPING(scov_2_subt_In_ch)
      SCOV2_SUBTYPING.out.set{scov2_subtyped_ch}
    }

  // === 6 - write final classification report

  if (!params.do_scov2_subtyping == true){
    scov2_subtyped_ch = Channel.empty()
  }
  qc_metrics_out_ch.no_subtyping_ch.concat(scov2_subtyped_ch)
    .set{report_in_ch}

  GENERATE_CLASSIFICATION_REPORT(report_in_ch)
}

def __check_if_params_file_exist(param_name, param_value){
  def error = 0

  if (!(param_value==null)){
    param_file = file(param_value)
    if (!param_file.exists()){
      log.error("${param_file} does not exist")
      error +=1
    }
  }

  if (param_value==null){
    log.error("${param_name} must be provided")
    error +=1
  }
  return error
}

def check_main_params(){

    def errors = 0
    def valid_entry_points = ["sort_reads", "consensus_gen"]
    
    // check if execution mode is valid
    if (!valid_entry_points.contains(params.entry_point)){
        log.error("The execution mode provided (${params.entry_point}) is not valid. valid modes = ${valid_entry_points}")
        errors += 1
    }

    if (params.entry_point == "sort_reads"){
        errors += check_sort_reads_params()
    }

    if (params.entry_point=="consensus_gen"){
        // check if manifest was provided
        errors += __check_if_params_file_exist("consensus_mnf", params.consensus_mnf)
    }

    errors += check_classification_report_params()
    if (errors > 0) {
        log.error("Parameter errors were found, the pipeline will not run.")
        exit 1
    }
}
/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html
 */
workflow.onComplete {
  // Log colors ANSI codes
  
  println """
  Pipeline execution summary
  ---------------------------
  Completed at : ${ANSI_GREEN}${workflow.complete}${ANSI_RESET}
  Duration     : ${ANSI_GREEN}${workflow.duration}${ANSI_RESET}
  Success      : ${workflow.success ? ANSI_GREEN : ANSI_RED}${workflow.success}${ANSI_RESET}
  Results Dir  : ${ANSI_GREEN}${file(params.outdir)}${ANSI_RESET}
  Work Dir     : ${ANSI_GREEN}${workflow.workDir}${ANSI_RESET}
  Exit status  : ${ANSI_GREEN}${workflow.exitStatus}${ANSI_RESET}
  Error report : ${ANSI_GREEN}${workflow.errorReport ?: '-'}${ANSI_RESET}
  """.stripIndent()
}

def parse_mnf(consensus_mnf) {
    /*
    -----------------------------------------------------------------
    Parses the manifest file to create a channel of metadata and 
    FASTQ file pairs.

    Also, checks if there are sample_id duplicated and/or containing
    non alphanumeric characters. Only exception accepted is "_", as
    long as it is not two consecutives "__".

    -----------------------------------------------------------------

    - **Input**:
        consensus_mnf (path to the manifest file)

    - **Output**: 
        Channel with tuples of metadata and FASTQ file pairs.

    -----------------------------------------------------------------
    */
    // Read manifest file into a list of rows
    def mnf_rows = Channel.fromPath(consensus_mnf)
                          | splitCsv(header: true, sep: ',')

    // Collect sample IDs and validate
    def sample_ids = []
    def errors = 0
    
    def errors_ch = mnf_rows.map { row ->
        def sample_id = row.sample_id

        // Check if sample_id is empty
        if (!sample_id) {
            log.error("Empty sample_id detected.")
            errors += 1
        } else {
            // Check for unique sample IDs
            if (sample_ids.contains(sample_id)) {
                log.error("${sample_id} is duplicated")
                errors += 1
            } else {
                sample_ids << sample_id
            }
        
            // Check if sample_id is alphanumeric, allows underscores but not consecutive
            if (!sample_id.matches(/^(?!.*__)[A-Za-z0-9_]+$/)) {
                log.error("Non alphanumeric sample id ${sample_id} ['_' is permitted]")
                errors += 1
            }
            return errors
        }
        }
        // be sure that the number of errors is evaluated after all rows are processed
        .collect() 
        // kill the pipeline if errors are found
        .subscribe{ v ->
        if (errors > 0) {
            log.error("${errors} critical errors in the manifest were detected. Please check README for more details.")
            exit 1
        }
    }

    // If validation passed, create the channel as before
    def mnf_ch = mnf_rows.map { row -> 
                    // set meta
                    def meta = [id: row.sample_id]
                    // set files
                    def reads = [row.reads_1, row.reads_2]
                    // declare channel shape
                    tuple(meta, reads)
                 }

    return mnf_ch // tuple(meta, [fastq_pairs])
}