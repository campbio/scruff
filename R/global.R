# Declare global vars to get rid of NOTE
vars <- c(".", "gene_id", "umi", "inferred_umi", "median", "experiment", "lane",
          "filename", "project", "cell_index", "reads", "percent_assigned",
          "read1_path", "read2_path", "min.phred1", "length1", "bc_correct",
          "barcode", "read2", "qtring2", "rname2", "fastq_path", "read1", 
          "qtring1", "rname1", "type", "gene_name", "gene_biotype", 
          "avg_reads_per_corrected_umi", "number_of_cells", ".x",
          "avg_reads_per_umi", "total_counts", "protein_coding_genes", 
          "number_of_cells", "reads_mapped_to_genes", "reads_mapped_to_genome",
          "median_reads_per_corrected_umi", "median_reads_per_umi", "mt_counts",
          "CB", "readname", "NH", "GX", "MM", "cell_barcode", "cbtop10000",
          "v2chemistry", "v1chemistry", "v3chemistry", "geneid", ".getLevel",
          ".getTxdt", "transcript_id", ".getTxdt", "transcript_name",
          ".getTxdt", ".transRect", "exon_number", ".transArrow", ".transText",
          "i", "seqid", "x1", "x2", "y1", "y2", "x", "y", "cb",
          "protein_coding_counts", "level", "bamExample")
utils::globalVariables(vars)

