glm_analy <- function(data_analysis_glm_design, 
                      trait_name, 
                      outcome,
                      outcome_abbr, 
                      covar_model1, covar_model2, covar_model3 ,
                      cata_x ,
                      cata_y 
                      
                      ) {
  options(warn = 1)
  print("run svy_tableone")


  covar_vector <- unlist(strsplit(gsub("\\s*\\+\\s*", ",", covar_model3), ","))
  covar_vector <- trimws(covar_vector)

  covar_vector <- c(covar_vector)
  covar_vector <- covar_vector[covar_vector != "ageyr"]
  
  
  file_path_svy_tableone <- paste0("./regression/", outcome_abbr, "/table_one_", outcome_abbr, "_", trait_name, ".xlsx")

  cv <- c(if(cata_x == "c") trait_name else NULL, if(cata_y == "c") outcome else NULL , "ageyr")
  gv <- c(if(cata_x == "g") trait_name else NULL, if(cata_y == "g") outcome else NULL, covar_vector)
  
  svy_tableone(design = data_analysis_glm_design,
               cv = cv[!is.null(cv)],
               gv = gv[!is.null(gv)],
               view = FALSE, xlsx = file_path_svy_tableone)

  covar_model1 <- gsub("\n", "", covar_model1)
  covar_model2 <- gsub("\n", "", covar_model2)
  covar_model3 <- gsub("\n", "", covar_model3)

  model1 <- as.formula(paste(outcome," ~", trait_name, "+", covar_model1))
  model2 <- as.formula(paste(outcome," ~", trait_name, "+", covar_model2))
  model3 <- as.formula(paste(outcome," ~", trait_name, "+", covar_model3))
  

  print("run model1")
  glm1 <- svyglm(model1, design = data_analysis_glm_design)
  print("run model2")
  glm2 <- svyglm(model2, design = data_analysis_glm_design)
  print("run model3")
  glm3 <- svyglm(model3, design = data_analysis_glm_design)

  source("F:/Research_Skills_and_Resources/lc_code_package/reg_table_svyglm_lc.R")
  glm1_table <- reg_table_svyglm_lc(glm1)
  glm2_table <- reg_table_svyglm_lc(glm2)
  glm3_table <- reg_table_svyglm_lc(glm3)
  


  fwrite(glm1_table, sep = "\t", file = paste0("./regression/",outcome_abbr,"/",outcome_abbr,"1_table_", trait_name, ".txt"))
  fwrite(glm2_table, sep = "\t", file = paste0("./regression/",outcome_abbr,"/",outcome_abbr,"2_table_", trait_name, ".txt"))
  fwrite(glm3_table, sep = "\t", file = paste0("./regression/",outcome_abbr,"/",outcome_abbr,"3_table_", trait_name, ".txt"))
  

  broom::tidy(glm1, exponentiate = TRUE, conf.int = TRUE) |> 
    write_xlsx(paste0("./regression/",outcome_abbr,"/",outcome_abbr,"1_table_", trait_name, ".xlsx"))
  broom::tidy(glm2, exponentiate = TRUE, conf.int = TRUE) |> 
    write_xlsx(paste0("./regression/",outcome_abbr,"/",outcome_abbr,"2_table_", trait_name, ".xlsx"))
  broom::tidy(glm3, exponentiate = TRUE, conf.int = TRUE) |> 
    write_xlsx(paste0("./regression/",outcome_abbr,"/",outcome_abbr,"3_table_", trait_name, ".xlsx"))
  

  
  assign("glm1_table", glm1_table, envir = .GlobalEnv)
  assign("glm2_table", glm2_table, envir = .GlobalEnv)
  assign("glm3_table", glm3_table, envir = .GlobalEnv)
  gc()
}

