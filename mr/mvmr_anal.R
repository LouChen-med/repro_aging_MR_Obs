
args = commandArgs(TRUE)
mergefold = as.character(args[1])

library(writexl)
library(MVMR)
library(MendelianRandomization)
library(openxlsx)
library(TwoSampleMR)
library(data.table)
library(writexl)
library(dplyr)
library(coloc)
library(openxlsx)
library(magrittr)
library(rlang)
library(ieugwasr)
library(stringr)
library(dplyr)
library(vroom)
library(data.table)


file_list <- list.files(mergefold, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)


for(filename in file_list) {
  
  datamvivw <- data.frame()
  dataMREgger <- data.frame()
  datamlasso <- data.frame()
  sres <- data.frame()
  my_df <- data.frame()
  
  
  SNPfinally<-read.csv(filename)
  
  nrowdata<-nrow(SNPfinally)
  print(nrowdata)
  exp<-SNPfinally
  
  filenamebase<-basename(filename)
  print(filenamebase)
  


  
  
  matches <- strcapture("^O_(.*?)_X_(.*?)_Y_(.*?).csv$", filenamebase, proto = list(str1 = character(), str2 = character(), str3 = character()))
  o<-matches$str1
  x<-matches$str2
  y<-matches$str3
  print(o)
  print(x)
  print(y)

  
  r_input <- mr_mvinput(
    bx = cbind(exp$Beta.x,exp$Beta.y),
    by= exp$Beta.o,
    bxse =cbind(exp$SE.x,exp$SE.y),
    byse = exp$SE.o,
    snp = exp$SNP,
    exposure = c(x,y),
    outcome = o)
  

  F.data <- format_mvmr(BXGs = cbind(exp$Beta.x,exp$Beta.y),
                        BYG =exp$Beta.o,
                        seBXGs = cbind(exp$SE.x,exp$SE.y),
                        seBYG = exp$SE.o,
                        RSID = exp$SNP)

  if (nrowdata > 1) {
    tryCatch({
      sres <- strength_mvmr(r_input = F.data, gencov = 0)
      sres <- rename(sres, `F-statistic_X` = exposure1, `F-statistic_Y` = exposure2)
      
      sres$X<-x
      sres$Y<-y
      sres$Outcome<-o
      
      sres$Group<-filenamebase11

      
      print(sres)
    }, error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    })
  }
  

  if (nrowdata > 1) {
    tryCatch({
      pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)
      pres

      my_df <- data.frame(matrix(0, nrow = 1, ncol = 2))

      row.names(my_df) <- c("")
      colnames(my_df) <- c("Q-Statistic", "P-value")

      my_df[1,1] <- pres$Qstat
      my_df[1,2] <- pres$Qpval
      my_df$X<-x
      my_df$Y<-y
      my_df$Outcome<-o
      my_df$Group<-filenamebase11


      print(my_df)
    }, error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    })
  }

  

  if (nrowdata > 2) {
    tryCatch({
      mvivw<-mr_mvivw(r_input,model="default",robust = FALSE,correl = FALSE,distribution = "normal",alpha = 0.05)
      print(mvivw)
      datamvivw<-data.frame(Exposure = c(x,y),
                            Method = c ("IVW","IVW"),
                            Outcome = c(o,o),
                            Beta = mvivw@Estimate,
                            SE= mvivw@StdError,
                            OR= exp(mvivw@Estimate),
                            ORlower = exp(mvivw@CILower),
                            ORupper = exp(mvivw@CIUpper),
                            P = mvivw@Pvalue)
      datamvivw$`Number_of_Variants` <- mvivw@SNPs
      datamvivw$`Residual_standard_error` <- mvivw@RSE
      datamvivw$`Heterogeneity_test_statistic` <- mvivw@Heter.Stat[1]
      datamvivw$`P_of_Heterogeneity` <- mvivw@Heter.Stat[2]
      datamvivw$Group<-filenamebase11

      print(datamvivw)
    }, error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    })
  }
  
  

  if (nrowdata > 2) {
    tryCatch({
      lasso<-mr_mvlasso(
        r_input,
        orientate = 1,
        distribution = "normal",
        alpha = 0.05,
        lambda = numeric(0))
      print(lasso)
      datamlasso<-data.frame(Exposure = c(x,y),
                             Method = c ("Lasso","Lasso"),
                             Outcome = c(o,o),
                             Beta = lasso@Estimate,
                             SE=lasso@StdError,
                             OR= exp(lasso@Estimate),
                             ORlower = exp(lasso@CILower),
                             ORupper = exp(lasso@CIUpper),
                             P = lasso@Pvalue)
      datamlasso$`Orientated_to_exposure`<-lasso@Orientate
      datamlasso$`Number_of_variants`<-lasso@SNPs
      datamlasso$`Number_of_valid_instruments `<-lasso@Valid
      datamlasso$`Tuning_parameter `<-lasso@Lambda
      datamlasso$Group<-filenamebase11

      print(datamlasso)
    }, error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    })
  }
  
  
  
  

  if (nrowdata > 3) {
    tryCatch({
      egger<-mr_mvegger(
        r_input,
        orientate = 1,
        correl = FALSE,
        distribution = "normal",
        alpha = 0.05)
      print(egger)
      dataMREgger<-data.frame(Exposure = c(x,y),
                              Method = c ("MR-Egger","MR-Egger"),
                              Outcome = c(o,o),
                              Beta = egger@Estimate,
                              SE= egger@StdError.Est,
                              OR= exp(egger@Estimate),
                              ORlower = exp(egger@CILower.Est),
                              ORupper = exp(egger@CIUpper.Est),
                              P = egger@Pvalue.Est)
      
      dataMREgger$`Orientated_to_exposure`<-egger@Orientate
      dataMREgger$`Number_of_Variants`<-egger@SNPs
      dataMREgger$`Intercept_Estimate`<-egger@Intercept
      dataMREgger$`Intercept_SE`<-egger@StdError.Int
      dataMREgger$`Intercept_CILower`<-egger@CILower.Int
      dataMREgger$`Intercept_CIUpper`<-egger@CIUpper.Int
      dataMREgger$`Intercept_Pvalue`<-egger@Pvalue.Int
      dataMREgger$`Residual_standard_error`<-egger@RSE
      dataMREgger$`Heterogeneity_test_statistic` <- egger@Heter.Stat[1]
      dataMREgger$`P_of_Heterogeneity` <- egger@Heter.Stat[2]
      dataMREgger$Group<-filenamebase11

      print(dataMREgger)
    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
      dataMREgger <- data.frame()
      NULL
    })
  }
  
  
  
  
  

  wb <- createWorkbook()

  addWorksheet(wb, "IVW")
  addWorksheet(wb, "MR-Egger")
  addWorksheet(wb, "LASSO")
  addWorksheet(wb, "Conditional_F-statistics")
  addWorksheet(wb, "Pleiotropy")
  
  data_frames <- c("datamvivw", "dataMREgger", "datamlasso", "sres", "my_df")
  

  for (df_name in data_frames) {

    if (!exists(df_name)) {

      assign(df_name, data.frame())
    }
  }
  

  writeData(wb, "IVW", datamvivw)
  writeData(wb, "MR-Egger", dataMREgger)
  writeData(wb, "LASSO", datamlasso)
  writeData(wb, "Conditional_F-statistics", sres)
  writeData(wb, "Pleiotropy", my_df)
  
  
  folderpathXYO<-paste0("//MVMR/MVMR_result","/","EXP_",x,"/","OUT_",o)

  dir.create(folderpathXYO, recursive = TRUE)
  
  filenamebase1 <- sub("\\.csv$", "", filenamebase)
  
  
  filenamexlsx <- paste0(folderpathXYO,"/Result_",filenamebase1,".xlsx")
  
  saveWorkbook(wb, filenamexlsx, overwrite = TRUE)
}

