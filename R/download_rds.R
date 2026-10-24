# Download and prepare reference RDS file from GitHub Release
download_rds <- function(rds_name, release_tag = "data-v1", force = TRUE) {
  sha256_list <- c(                                                                     
    "ChromHMM_hg38.rds" = "04a71a1fc8c3575acd3b34b312fe5e47d975d9089e7839606b7f9aa9f60b6771",
    "ChromHMM_mm10.rds" = "2aa86b1b3f05065ea20ce74a92ddbd0044064bab6566c05367124e8977b045e8",
    "ENCODE_cCRE_hg38.rds" = "e19b38bdc2eff9233488a603cd3bfcfee400ab07b5179d71cc45dced4ffd2023",
    "ENCODE_cCRE_mm10.rds" = "673b790f2ed501c15776be7e6a740c3d8c5ed7a9b4f5a573d0a940fddb73dd83",
    "GENCODE_gene_hg38.rds" = "0933580985d8caf5cfdf403c05bb72fb6d19e570ac3a61fb12cd96147f5b358f",
    "GENCODE_gene_mm10.rds" = "6a25c2db7f334a84f3db604710315c6583bb94a29531dbfc46fa00c28501cc2e",
    "GENCODE_genic_hg38.rds" = "f240cf0f9d8f42f8f5060fa67d063672963d25b11ee8b06cd849142dfb64f171",
    "GENCODE_genic_mm10.rds" = "3623621d97ed97a6de1b11cb9e9171caf60250fdb7b98255eef62404400dc183",
    "GENCODE_rep_hg38.rds" = "a054def0660c0852d411bfe633db929a1f75a9f715cc361db62b7b4118c001d4",
    "GENCODE_rep_mm10.rds" = "71b3a642a68e58f5f21bfb402f06fe83c092184877f054d0113864d48cb3e1c4",
    "knownGene_genic_hg38.rds" = "050b5da897e80f09e1a86acdc356b6354452157de3833671a1d9023aa0cfaf2a",
    "knownGene_genic_mm10.rds" = "3700fc5c6c6962530d4379c20a1879cff02cf9b3cabb7beaf8585b7538637506",
    "RepeatMasker_hg38.rds" = "217560420b545bee8977a50b7fe5bc4abca859fe1c588e738faca34292ed3a94",
    "RepeatMasker_mm10.rds" = "cbfa9f5d4d05f87500634676a72428cc5a1faf00553163cfcc3da99bfe6ff1e6",
    "TFBS_lib_hg38.rds" = "0f7f1d2407ad429a9617e11cf71654c7e2af4b2ca43a92ca8537f3cb99540d69",
    "TFBS_lib_mm10.rds" = "b36315592000f2bd3d7a102007dea1942d3b175b1d8fbdad2e607be132aeb203"
  )

  # Validate rds_name
  if (!rds_name %in% names(sha256_list)) {
    stop("Unknown file: ", rds_name)
  }

  cache_dir <- rappdirs::user_cache_dir("multiEpiCore")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

  rds_path <- file.path(cache_dir, rds_name)
  tmp_path <- paste0(rds_path, ".part")

  rds_url <- sprintf(
    "https://github.com/Qingjie-Yu/multiEpiCore/releases/download/%s/%s",
    release_tag, rds_name
  )

  needs_download <- force || !file.exists(rds_path)
  if (!needs_download) {
    sha256_local <- digest::digest(file = rds_path, algo = "sha256")
    if (!identical(tolower(sha256_local), tolower(sha256_list[[rds_name]]))) {
      message("Checksum mismatch for ", rds_name, ". Re-downloading...")
      needs_download <- TRUE
      unlink(rds_path, force = TRUE)
    }
  }

  if (needs_download) {
    max_try <- 3
    success <- FALSE
    for (i in seq_len(max_try)) {
      tryCatch(
        {
          message("Downloading: ", rds_name, " (attempt ", i, "/", max_try, ")")
          utils::download.file(rds_url, tmp_path, mode = "wb", quiet = TRUE)
          sha256_local <- digest::digest(file = tmp_path, algo = "sha256")
          if (!identical(tolower(sha256_local), tolower(sha256_list[[rds_name]]))) {
            stop("sha256 mismatch for ", rds_name,
                 "\nExpected: ", sha256_list[[rds_name]], "\nObserved: ", sha256_local)
          }
          file.rename(tmp_path, rds_path)
          success <- TRUE
          break
        },
        error = function(e) {
          message("Attempt ", i, " failed: ", e$message)
          unlink(tmp_path, force = TRUE)
        }
      )
      if (!success && i < max_try) Sys.sleep(2^i)
    }
    if (!success) stop("Failed to download ", rds_name, " after ", max_try, " attempts.")
  }

  normalizePath(rds_path)
}
