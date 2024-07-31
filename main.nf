process run_stenrich {

    tag ""
    publishDir "${params.output_dir}/STEnrich.pathway_results", mode:"copy", pattern:"${output_filename}"
    container "ghcr.io/web-mev/mev-spatialge-stenrich"
    cpus 4
    memory '16 GB'

    input:
        path raw_counts
        path coords_metadata

    output:
        path "${output_filename}"

    script:
        output_filename = "stenrich_pathways.json"
        """
        /usr/local/bin/run.sh ${raw_counts} ${coords_metadata} "${params.sample_name}" ${params.normalization_method} "${params.gene_set_database}" ${params.organism} ${params.gene_id_choice} ${params.xpos_col} ${params.ypos_col} ${output_filename}
        """
}

workflow {
    run_stenrich(params.raw_counts, params.coords_metadata)
}