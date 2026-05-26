#!/usr/bin/env Rscript

# ============================================================
# ADMIXTOOLS2 qpGraph analysis for Jewish population dataset
# ============================================================
#
# Input:
#   EIGENSTRAT dataset:
#     jews_plus_ref.geno
#     jews_plus_ref.snp
#     jews_plus_ref.ind
#
# Output:
#   f2 statistics
#   fitted admixture graphs
#   model comparison tables
#   PDF  plots
suppressPackageStartupMessages({
  library(admixtools)
  library(igraph)
})

# -----------------------------
# Configuration
# -----------------------------

prefix <- "jews_plus_ref"

f2_dir <- "results/qpgraph/f2_target_like"
outdir <- "results/qpgraph/graphs_target_like"

dir.create(f2_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Population set selected for the target-like qpGraph analysis.
# These labels must exactly match the third column of jews_plus_ref.ind.
pops <- c(
  "Mbuti.DG",
  "Natufian",
  "Ust_Ishim",
  "Levant_BA",
  "CHG",
  "Iran_N",
  "WHG",
  "Ashkenazi",
  "Bukharan",
  "Georgian_Jews",
  "Kurdistani",
  "Mountain"
)

outpop <- "Mbuti.DG"

numadmix_values <- c(3, 4, 5)
seeds <- c(101, 202, 303, 404, 505)

stop_gen <- 150

# -----------------------------
# Helper functions
# -----------------------------

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path)
  }
}

write_tsv <- function(x, file) {
  write.table(
    x,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

read_edges <- function(file) {
  read.table(
    file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

safe_qpgraph <- function(f2_blocks, graph) {
  tryCatch(
    {
      qpgraph(
        f2_blocks,
        graph,
        return_fstats = TRUE
      )
    },
    error = function(e) {
      message("qpgraph failed: ", conditionMessage(e))
      return(NULL)
    }
  )
}

plot_graph_pdf <- function(edges, file, width = 18, height = 12) {
  pdf(file, width = width, height = height)
  plot_graph(edges, textsize = 2)
  dev.off()
}

# -----------------------------
# Input checks
# -----------------------------

geno_file <- paste0(prefix, ".geno")
snp_file <- paste0(prefix, ".snp")
ind_file <- paste0(prefix, ".ind")

check_file_exists(geno_file)
check_file_exists(snp_file)
check_file_exists(ind_file)

ind <- read.table(
  ind_file,
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(ind) <- c("IID", "SEX", "POP")

available_pops <- sort(unique(ind$POP))
missing_pops <- setdiff(pops, available_pops)

cat("Available selected populations:\n")
print(table(ind$POP[ind$POP %in% pops]))

if (length(missing_pops) > 0) {
  cat("\nMissing populations:\n")
  print(missing_pops)
  stop("Some requested populations are absent from the .ind file.")
}

cat("\nSNP count in merged dataset:\n")
print(system(paste("wc -l", snp_file), intern = TRUE))

cat("\nIndividual count in merged dataset:\n")
print(system(paste("wc -l", ind_file), intern = TRUE))

# -----------------------------
# Compute f2 statistics
# -----------------------------

cat("\nExtracting f2 statistics...\n")

extract_f2(
  prefix,
  f2_dir,
  pops = pops,
  maxmiss = 1,
  overwrite = TRUE
)

f2_blocks <- f2_from_precomp(f2_dir)

cat("\nf2 statistics loaded successfully.\n")

# -----------------------------
# Run one graph search
# -----------------------------

run_one_graph_search <- function(numadmix, seed) {

  run_name <- paste0("targetlike_numadmix", numadmix, "_seed", seed)
  run_prefix <- file.path(outdir, run_name)

  score_file <- paste0(run_prefix, "_score.tsv")

  if (file.exists(score_file)) {
    cat("\nSkipping existing run:", run_name, "\n")

    old_score <- read.table(
      score_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )

    return(old_score)
  }

  cat("\n============================================================\n")
  cat("Running:", run_name, "\n")
  cat("numadmix:", numadmix, "\n")
  cat("seed:", seed, "\n")
  cat("stop_gen:", stop_gen, "\n")
  cat("============================================================\n")

  set.seed(seed)

  opt_results <- tryCatch(
    {
      find_graphs(
        f2_blocks,
        numadmix = numadmix,
        stop_gen = stop_gen,
        outpop = outpop
      )
    },
    error = function(e) {
      message("find_graphs failed for ", run_name, ": ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(opt_results)) {
    return(NULL)
  }

  saveRDS(
    opt_results,
    paste0(run_prefix, "_all_results.rds")
  )

  winner_idx <- which.min(opt_results$score)
  winner_score <- opt_results$score[[winner_idx]]
  winner_edges <- opt_results$edges[[winner_idx]]

  write_tsv(
    winner_edges,
    paste0(run_prefix, "_edges.tsv")
  )

  plot_graph_pdf(
    winner_edges,
    paste0(run_prefix, "_graph.pdf")
  )

  graph <- edges_to_igraph(winner_edges)

  out <- safe_qpgraph(f2_blocks, graph)

  if (is.null(out)) {
    score_table <- data.frame(
      run = run_name,
      numadmix = numadmix,
      seed = seed,
      score = winner_score,
      worst_residual = NA,
      abs_worst_residual = NA
    )

    write_tsv(score_table, score_file)
    return(score_table)
  }

  write_tsv(
    out$edges,
    paste0(run_prefix, "_qpgraph_refitted_edges.tsv")
  )

  score_table <- data.frame(
    run = run_name,
    numadmix = numadmix,
    seed = seed,
    score = out$score,
    worst_residual = out$worst_residual,
    abs_worst_residual = abs(out$worst_residual)
  )

  write_tsv(score_table, score_file)

  saveRDS(
    out,
    paste0(run_prefix, "_qpgraph_out.rds")
  )

  # Save additional data.frame elements if present in the qpgraph output.
  for (nm in names(out)) {
    obj <- out[[nm]]

    if (is.data.frame(obj)) {
      write_tsv(
        obj,
        paste0(run_prefix, "_", nm, ".tsv")
      )
    }
  }

  cat("\nFinished:", run_name, "\n")
  print(score_table)

  return(score_table)
}

# -----------------------------
# Run batch of models
# -----------------------------

all_scores <- list()

for (k in numadmix_values) {
  for (s in seeds) {
    res <- run_one_graph_search(
      numadmix = k,
      seed = s
    )

    if (!is.null(res)) {
      all_scores[[length(all_scores) + 1]] <- res
    }
  }
}

scores <- do.call(rbind, all_scores)

scores_by_score <- scores[order(scores$score), ]
scores_by_residual <- scores[order(scores$abs_worst_residual, scores$score), ]

write_tsv(
  scores_by_score,
  file.path(outdir, "model_comparison_sorted_by_score.tsv")
)

write_tsv(
  scores_by_residual,
  file.path(outdir, "model_comparison_sorted_by_abs_worst_residual.tsv")
)

cat("\n============================================================\n")
cat("Best models by score\n")
cat("============================================================\n")
print(head(scores_by_score, 10))

cat("\n============================================================\n")
cat("Best models by absolute worst residual\n")
cat("============================================================\n")
print(head(scores_by_residual, 10))

# -----------------------------
# Save final best model separately
# -----------------------------

best_run <- scores_by_score$run[1]

cat("\nBest run by score:\n")
print(best_run)

best_edges_file <- file.path(outdir, paste0(best_run, "_qpgraph_refitted_edges.tsv"))

if (!file.exists(best_edges_file)) {
  best_edges_file <- file.path(outdir, paste0(best_run, "_edges.tsv"))
}

best_edges <- read_edges(best_edges_file)

write_tsv(
  best_edges,
  file.path(outdir, "best_model_edges.tsv")
)

plot_graph_pdf(
  best_edges,
  file.path(outdir, "best_model_graph.pdf")
)

cat("\nFinal outputs written to:\n")
cat(outdir, "\n")

cat("\nMain output files:\n")
cat(file.path(outdir, "model_comparison_sorted_by_score.tsv"), "\n")
cat(file.path(outdir, "model_comparison_sorted_by_abs_worst_residual.tsv"), "\n")
cat(file.path(outdir, "best_model_edges.tsv"), "\n")
cat(file.path(outdir, "best_model_graph.pdf"), "\n")
