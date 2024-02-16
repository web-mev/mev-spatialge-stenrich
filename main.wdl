workflow STEnrich {

    File raw_counts
    File coords_metadata
    String sample_name
    String normalization_method
    String gene_set_database
    String organism
    String gene_id_choice

    call runSTEnrich {
        input:
            raw_counts = raw_counts,
            coords_metadata = coords_metadata,
            sample_name = sample_name,
            normalization_method = normalization_method,
            gene_set_database = gene_set_database,
            organism = organism
            gene_id_choice = gene_id_choice
    }

    output {
        File pathway_results = runSTEnrich.enriched_pathways
    }

}

task runSTEnrich {

    File raw_counts
    File coords_metadata
    String sample_name
    String normalization_method
    String gene_set_database
    String organism
    String gene_id_choice
    String output_filename = "stenrich_pathways.json"
    Int disk_size = 50

    command <<<
        /usr/local/bin/run.sh \
            ${raw_counts} \ 
            ${coords_metadata} \
            "${sample_name}" \
            ${normalization_method} \
            "${gene_set_database}" \
            ${organism} \
            ${gene_id_choice} \
            ${output_filename}
    >>>

    output {
        File enriched_pathways = "${output_filename}"
    }

    runtime {
        docker: "ghcr.io/web-mev/mev-spatialge-stenrich"
        cpu: 4
        memory: "30 G"
        disks: "local-disk " + disk_size + " HDD"
    }
}


