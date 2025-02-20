cox_analy <- function(data_analysis_cox_design, 
                      trait_name, covar_model1, covar_model2, covar_model3,
                      cata_x
                      
                      
                      ) {
  options(warn = 1)
  print("run svy_tableone")

  covar_vector <- unlist(strsplit(gsub("\\s*\\+\\s*", ",", covar_model3), ","))
  covar_vector <- trimws(covar_vector)

  covar_vector <- c(covar_vector)
  covar_vector <- covar_vector[covar_vector != "ageyr"]
  
  
  file_path_svy_tableone <- paste0("./regression/cox/table_one_cox_", trait_name, ".xlsx")

  cv <- c(if(cata_x == "c") trait_name else NULL, "ageyr")
  gv <- c(if(cata_x == "g") trait_name else NULL, covar_vector)
  
  svy_tableone(design = data_analysis_cox_design,
               cv = cv[!is.null(cv)],
               gv = gv[!is.null(gv)],
               view = FALSE, xlsx = file_path_svy_tableone)
  
  covar_model1 <- gsub("\n", "", covar_model1)
  covar_model2 <- gsub("\n", "", covar_model2)
  covar_model3 <- gsub("\n", "", covar_model3)

  model1 <- as.formula(paste("Surv(permth_exm, mortstat) ~", trait_name, "+", covar_model1))
  model2 <- as.formula(paste("Surv(permth_exm, mortstat) ~", trait_name, "+", covar_model2))
  model3 <- as.formula(paste("Surv(permth_exm, mortstat) ~", trait_name, "+", covar_model3))
  

  print("run model1")
  cox1 <- svycoxph(model1, design = data_analysis_cox_design)
  print("run model2")
  cox2 <- svycoxph(model2, design = data_analysis_cox_design)
  print("run model3")
  cox3 <- svycoxph(model3, design = data_analysis_cox_design)
  
  source("F:/Research_Skills_and_Resources/lc_code_package/reg_table_svycoxph_lc.R")
  cox1_table <- reg_table_svycoxph_lc(cox1)
  cox2_table <- reg_table_svycoxph_lc(cox2)
  cox3_table <- reg_table_svycoxph_lc(cox3)
  

  fwrite(cox1_table, sep = "\t", file = paste0("./regression/cox/cox1_table_", trait_name, ".txt"))
  fwrite(cox2_table, sep = "\t", file = paste0("./regression/cox/cox2_table_", trait_name, ".txt"))
  fwrite(cox3_table, sep = "\t", file = paste0("./regression/cox/cox3_table_", trait_name, ".txt"))
  

  broom::tidy(cox1, exponentiate = TRUE, conf.int = TRUE) |> 
    write_xlsx(paste0("./regression/cox/cox1_table_", trait_name, ".xlsx"))
  broom::tidy(cox2, exponentiate = TRUE, conf.int = TRUE) |> 
    write_xlsx(paste0("./regression/cox/cox2_table_", trait_name, ".xlsx"))
  broom::tidy(cox3, exponentiate = TRUE, conf.int = TRUE) |> 
    write_xlsx(paste0("./regression/cox/cox3_table_", trait_name, ".xlsx"))
  

  assign("cox1_table", cox1_table, envir = .GlobalEnv)
  assign("cox2_table", cox2_table, envir = .GlobalEnv)
  assign("cox3_table", cox3_table, envir = .GlobalEnv)
  gc()
}


