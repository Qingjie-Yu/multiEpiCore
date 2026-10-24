# Peak Profiling Pipeline
#
# End-to-end peak calling and annotation pipeline: calls peaks from bedGraph
# files, then optionally runs downstream annotation.
#
# Parameters:
#   bedgraph_path:
#     Path(s) to input bedGraph files
#   out_dir:
#     Output directory for results. Default: "./"
#   ref_genome:
#     Reference genome version. Default: "hg38". Supported: "hg38", "mm10"
#   min_cov:
#     Minimum mean coverage (auc / length) pre-filter applied to candidate
#     peak blocks before statistical testing. Default: 2
#   auc_top_pct:
#     Retain only the top X% of blocks by AUC after statistical filtering.
#     Default: 0.1
#   qvalue_cutoff:
#     BH-adjusted q-value cutoff for peak significance. Default: 0.05
#   fc_cutoff:
#     Minimum fold change (signal vs. background) cutoff. Default: 2
#   apply_annotation:
#     Logical. Controls whether to annotate the called peaks. When TRUE, runs
#     all three annotation types on the peak set:
#       - genomic distribution annotation (genic / ccre / chromhmm / repeat)
#       - pathway annotation (MSigDB gene set enrichment)
#       - TFBS enrichment (summit-anchored, matched-control based)
#     Default: TRUE
#   plot:
#     Whether to generate annotation plots. Default: TRUE
peak_profiling <- function(bedgraph_path, out_dir = "./", ref_genome = "hg38", min_cov = 2, auc_top_pct = 0.1, qvalue_cutoff = 0.05, fc_cutoff = 2, apply_annotation = TRUE, plot = TRUE) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Step1: peak calling
  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("  Peak Calling")
  cat("\n", strrep("=", 40), "\n", sep = "")
  peak_result <- peak_calling(bedgraph_path = bedgraph_path, out_dir = out_dir, ref_genome = ref_genome, min_cov = min_cov, auc_top_pct = auc_top_pct, qvalue_cutoff = qvalue_cutoff, fc_cutoff = fc_cutoff)

  # Step2: annotation
  if (apply_annotation) {
    cat("\n", strrep("=", 40), "\n", sep = "")
    cat("  Annotation")
    cat("\n", strrep("=", 40), "\n", sep = "")
    peak_path <- unlist(peak_result)

    distribution_dir <- file.path(out_dir, "genomic_distribution")
    peak_genomic_distribution(peak_path = peak_path, out_dir = distribution_dir, ref_genome = ref_genome, plot = plot)

    gene_dir <- file.path(out_dir, "gene_annotation")
    peak_gene_annotation(peak_path = peak_path, out_dir = gene_dir, ref_genome = ref_genome)

    pathway_dir <- file.path(out_dir, "pathway_annotation")
    peak_pathway_annotation(peak_path = peak_path, out_dir = pathway_dir, ref_genome = ref_genome, plot = plot)

    tfbs_dir <- file.path(out_dir, "TFBS_enrichment")
    peak_TFBS_enrichment(peak_path = peak_path, out_dir = tfbs_dir, ref_genome = ref_genome, plot = plot)
  }
}