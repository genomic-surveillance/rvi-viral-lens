
//
// Since files can (currently) only be published via processes, these 
// dummy process acts to publsh any file we like. Two processes are 
// provided: one for publishing files with an associated meta
// ( for per-consensus publishing), and one without (for run-level)
//

process publish_consensus_files {
    label "consensus_output"

    input:
        tuple val(meta), path(files)

    output:
        tuple val(meta), path(files)

    script:
    """
    echo "Published: ${files}"
    """
}

process publish_run_files {
   label "run_output"

    input:
        path(files)

    output:
        path(files)

    script:
    """
    echo "Published: ${files}"
    """
}
