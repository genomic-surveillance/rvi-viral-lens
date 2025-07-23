

process write_single_properties_json {
    label "sample_output"

    input:
        val(meta)

    output:
        tuple val(meta), path("${meta.id}.properties.json")

    script:
        def metaJson = groovy.json.JsonOutput.toJson(meta)
        def formatted = groovy.json.JsonOutput.prettyPrint(metaJson)

        """
        echo '${formatted}' > ${meta.id}.properties.json
        """
}

process write_collated_properties_json {
    label "run_output"

    input:
        val(meta)
        val(output_filename_base)

    output:
        tuple val(meta), path(output_filename)

    script:
        output_filename = "${output_filename_base}.json"
        def metaJson = groovy.json.JsonOutput.toJson(meta)
        def formatted = groovy.json.JsonOutput.prettyPrint(metaJson)
    
        """
        echo '${formatted}' > ${output_filename}
        """
}