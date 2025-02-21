papalexi_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/papalexi/eccite_screen/")

##########################
# STEP 0: LOAD ONDISC DATA
##########################
# gene and grna odm files
gene_odm_fp <- paste0(papalexi_dir, "gene/matrix.odm")
grna_odm_fp <- paste0(papalexi_dir, "grna_expression/matrix.odm")

# mm odm metadata fp
mm_metadata_fp <- paste0(papalexi_dir, "multimodal_metadata.rds")

# construct mm odm
mm_odm <- ondisc::read_multimodal_odm(odm_fps = c(gene_odm_fp, grna_odm_fp),
                                      multimodal_metadata_fp = mm_metadata_fp)

#######################################
# STEP 1: CREATE GENE EXPRESSION MATRIX
#######################################
set.seed(9)
feat_ids <- mm_odm |>
  ondisc::get_modality("gene") |>
  ondisc::get_feature_ids()

mt_feats <- grep(pattern = "^MT-", x = feat_ids, value = TRUE)
my_feats <- c(sample(feat_ids, 285), sample(mt_feats, 5))

exp_mat <- mm_odm@modalities$gene[my_feats,]
response_matrix_lowmoi <- exp_mat[[my_feats,]]
rownames(response_matrix_lowmoi) <- my_feats

#######################################
# STEP 2: CREATE GRNA EXPRESSION MATRIX
#######################################
grna_odm <- mm_odm |>
  ondisc::get_modality("grna_expression")
grna_matrix_lowmoi <- grna_odm[[seq(1, nrow(grna_odm)),]]
rownames(grna_matrix_lowmoi) <- grna_odm |> ondisc::get_feature_ids()

####################################
# STEP 3: THE GLOBAL CELL COVARIATES
####################################
global_cell_covariates <- mm_odm |> ondisc::get_cell_covariates()
rownames(global_cell_covariates) <- NULL
covariate_data_frame_lowmoi <- global_cell_covariates |>
  dplyr::select(gene_n_umis, gene_n_nonzero, bio_rep, p_mito) |>
  dplyr::rename("response_n_umis" = "gene_n_umis", "response_n_nonzero" = "gene_n_nonzero")

##########################
# STEP 4: GRNA GROUP TABLE
##########################
grna_features <- mm_odm |>
   ondisc::get_modality("grna_expression") |>
   ondisc::get_feature_covariates()
grna_group_data_frame_lowmoi <- data.frame(grna_id = rownames(grna_features),
                                           grna_group = factor(grna_features$target)) |>
  dplyr::arrange(grna_group)

################################################
# STEP 5: SORT ACCORDING TO BIOLOGICAL REPLICATE 
################################################
cell_order <- order(global_cell_covariates$bio_rep)
response_matrix_lowmoi <- response_matrix_lowmoi[,cell_order]
grna_matrix_lowmoi <- grna_matrix_lowmoi[,cell_order]
covariate_data_frame_lowmoi <- covariate_data_frame_lowmoi[cell_order,]

######################################
# STEP 6: SAVE THE DATA IN THE PACKAGE
######################################
usethis::use_data(response_matrix_lowmoi, grna_matrix_lowmoi, covariate_data_frame_lowmoi, grna_group_data_frame_lowmoi, overwrite = TRUE)
