# replace the values of mat_template by the one of mat_values (keeping dimnames)
change_matrix_values <- function(mat_template, mat_values){
  dimnames(mat_values) <- dimnames(mat_template)
  mat_values
}

# takes a template mofa_obj and replaces values of the W and Z matrices
# with the ones of the MCMC simulation stored in the `draws` object
update_Z_and_W <- function(mofa_obj, draws){
  
  # make Z and W matrices from MCMC simulations
  posterior_mean <- summary(draws)$statistics[,"Mean"] %>%
    enframe(name = "element", value="mean") %>% 
    separate(element, into=c("matrix", "row_idx", "col_idx","end"), sep="\\[|\\,|\\]") %>% 
    mutate(row_idx = as.integer(row_idx)) %>% 
    mutate(col_idx = as.integer(col_idx)) %>% 
    group_by(matrix) %>% 
    group_split()
  
  names(posterior_mean) <- map_chr(posterior_mean, ~.x$matrix[1])
  
  Z_mcmc <- matrix(nrow = max(posterior_mean[["Z"]]$row_idx), ncol = max(posterior_mean[["Z"]]$col_idx))
  walk(1:nrow(posterior_mean[["Z"]]),function(i){
    Z_mcmc[posterior_mean[["Z"]]$row_idx[i], posterior_mean[["Z"]]$col_idx[i]] <<- posterior_mean[["Z"]]$mean[i]
  })
  
  W_drug_mcmc <- matrix(nrow = max(posterior_mean[["W_drug"]]$row_idx), ncol = max(posterior_mean[["W_drug"]]$col_idx))
  walk(1:nrow(posterior_mean[["W_drug"]]),function(i){
    W_drug_mcmc[posterior_mean[["W_drug"]]$row_idx[i], posterior_mean[["W_drug"]]$col_idx[i]] <<- posterior_mean[["W_drug"]]$mean[i]
  })
  
  W_met_mcmc <- matrix(nrow = max(posterior_mean[["W_met"]]$row_idx), ncol = max(posterior_mean[["W_met"]]$col_idx))
  walk(1:nrow(posterior_mean[["W_met"]]),function(i){
    W_met_mcmc[posterior_mean[["W_met"]]$row_idx[i], posterior_mean[["W_met"]]$col_idx[i]] <<- posterior_mean[["W_met"]]$mean[i]
  })
  
  W_rna_mcmc <- matrix(nrow = max(posterior_mean[["W_rna"]]$row_idx), ncol = max(posterior_mean[["W_rna"]]$col_idx))
  walk(1:nrow(posterior_mean[["W_rna"]]),function(i){
    W_rna_mcmc[posterior_mean[["W_rna"]]$row_idx[i], posterior_mean[["W_rna"]]$col_idx[i]] <<- posterior_mean[["W_rna"]]$mean[i]
  })
  
  # make new MOFA object and overwrite W and Z matrices with new ones
  mofa_obj_mcmc <- mofa_obj
  mofa_obj_mcmc@expectations$Z$group1 <- change_matrix_values(mofa_obj_mcmc@expectations$Z$group1, Z_mcmc)
  mofa_obj_mcmc@expectations$W_drug$group1 <- change_matrix_values(mofa_obj_mcmc@expectations$W_drug$group1, W_drug_mcmc)
  mofa_obj_mcmc@expectations$W_met$group1 <- change_matrix_values(mofa_obj_mcmc@expectations$W_met$group1, W_met_mcmc)
  mofa_obj_mcmc@expectations$W_rna$group1 <- change_matrix_values(mofa_obj_mcmc@expectations$W_rna$group1, W_rna_mcmc)
  
  # re-arrange factors by variance explained
  mofa_obj_mcmc <- calculate_variance_explained(mofa_obj_mcmc)$r2_per_factor$group1 %>% 
    rowSums() %>% 
    order(decreasing = T) %>% 
    subset_factors(mofa_obj_mcmc, .)
  
  mofa_obj_mcmc
}

# like update Z and W, but with a result of opt() instead of mcmc()
update_Z_and_W_opt <- function(mofa_obj, opt_res) {
  # make new MOFA object and overwrite W and Z matrices with new computed ones
  mofa_obj_opt <- mofa_obj
  mofa_obj_opt@expectations$Z$group1 <- change_matrix_values(mofa_obj_opt@expectations$Z$group1, opt_res$par$Z)
  mofa_obj_opt@expectations$W_drug$group1 <- change_matrix_values(mofa_obj_opt@expectations$W_drug$group1, opt_res$par$W_drug)
  mofa_obj_opt@expectations$W_met$group1 <- change_matrix_values(mofa_obj_opt@expectations$W_met$group1, opt_res$par$W_met)
  mofa_obj_opt@expectations$W_rna$group1 <- change_matrix_values(mofa_obj_opt@expectations$W_rna$group1, opt_res$par$W_rna)
  
  # re-arrange factors by variance explained
  mofa_obj_opt <- calculate_variance_explained(mofa_obj_opt)$r2_per_factor$group1 %>% 
    rowSums() %>% 
    order(decreasing = T) %>% 
    subset_factors(mofa_obj_opt, .)
  
  mofa_obj_opt
}