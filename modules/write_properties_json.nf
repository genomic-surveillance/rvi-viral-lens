

process write_single_properties_json {
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
    input:
        val(metalist)
        val(output_filename_base)

    output:
        path(output_filename)

    script:
        output_filename = "${output_filename_base}.json"
        def metaJson = groovy.json.JsonOutput.toJson(metalist)
        def formatted = groovy.json.JsonOutput.prettyPrint(metaJson)
    
        """
        echo '${formatted}' > ${output_filename}
        """
}