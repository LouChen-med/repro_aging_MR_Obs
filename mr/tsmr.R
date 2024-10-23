
args = commandArgs(TRUE)
expfold = as.character(args[1])
outcomepath = as.character(args[2])
bigpath = as.character(args[3])



if (!file.exists(bigpath)) {
  dir.create(bigpath, recursive = TRUE)
}
setwd(bigpath)
library(TwoSampleMR)
library(data.table)
library(writexl)
library(MVMR)
library(MendelianRandomization)
library(openxlsx)
library(dplyr)
library(coloc)
library(tidyverse)
library(openxlsx)
library(LDlinkR) 

find_proxy <- function(snp = "rs1", r2_threshold = 0.8, build = "37", pop = "EUR") {
   print(paste0("find proxy for ",snp))
  proxies <- NULL

  tokens <- c(" "
  )
  

  random_token <- tokens[1]
  #random_token <- tokens[1]
  print(paste0("token is ",random_token))
  proxies<-NULL

  success <- FALSE
  while(!success) {
    tryCatch({
      proxies <- LDlinkR::LDproxy(snp = snp, 
                                  pop = pop, 
                                  r2d = "r2", 
                                  token = random_token,  
                                  genome_build = "grch37", 
                                  file = FALSE)
      if (!is.null(proxies)) {
        success <- TRUE  
      }
    }, error = function(e) {
      if(grepl("Timeout was reached", e$message) || grepl("Error in the HTTP2 framing layer", e$message)) {
        message("link wifi problem ")
        message(e$message)
        message("Error in the HTTP2 framing layer sleep and retry...")
        Sys.sleep(60)  
      } else {
        stop("Error other", e$message)  
      }
    })
  }
  
  
  
  if (nrow(proxies) == 1 & length(proxies) != 10) {
    if (substr(proxies$`data_out[1, 1]`, 3, 7) == "error") {
      print(proxies)
      print("there is error, proxy=null")  
      proxies<-NULL
    }
  }
  
  if(is.null(proxies)) {
    print(paste0("no proxy for ",snp))
    return(NULL)
  } else {
    proxy_df<-as.data.frame(proxies)
    head(proxy_df)
    proxy_df<- proxy_df[proxy_df$R2 > r2_threshold, ]
    
    if(nrow(proxy_df)==0 ){ 
      print(paste0("no proxy for ",snp))
      return(NULL)}
    if(nrow(proxy_df)>0 ){
      print(paste0("there has proxy for ",snp))
      
      proxy_df <- proxy_df %>%
        filter(grepl("^[ATCG]+=[ATCG]+,[ATCG]+=[ATCG]+$", Correlated_Alleles)) %>%# Match one or more ATCG characters, followed by an equals sign (=), then another group of one or more ATCG characters, separated by a comma (,).
        separate(Correlated_Alleles, into = c("target_a1.outcome", "proxy_a1.outcome", "target_a2.outcome", "proxy_a2.outcome"), sep = "[=,]") %>%
        mutate(
          target_snp.outcome = snp, 
          proxy_snp.outcome = RS_Number,
        )
      
      head(proxy_df)
      
      filtered_df<-subset(proxy_df, select = c(R2, target_a1.outcome, proxy_a1.outcome ,
                                               target_a2.outcome ,proxy_a2.outcome,
                                               target_snp.outcome, proxy_snp.outcome
      ))
      result_dt <- filtered_df
      result_dt <- result_dt[order(-result_dt$R2), ]
      result_dt$row_number <- seq_len(nrow(result_dt))
      print(head(result_dt,n=2))
      return(result_dt)
    }
  }
}


out_data<-fread(outcomepath)
names(out_data) <- gsub("exposure", "outcome", names(out_data))
out_data$effect_allele.outcome <- toupper(out_data$effect_allele.outcome)
out_data$other_allele.outcome <- toupper(out_data$other_allele.outcome)
outcome_name<-sub(".*/f_(.*).txt(.*$)?", "\\1", outcomepath)
out_data<-as.data.frame(out_data)
print(head(out_data))
print(str(out_data))


folder_path <- expfold

last_field <- basename(folder_path)
print(last_field)  
file_list <- list.files(folder_path, pattern = "\\.csv$")




output_directory1  <- paste0("结局因素_", outcome_name)
output_directory <- file.path(paste0(bigpath,"/result/harmon_data/", last_field,"/",output_directory1))
if (!file.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}
wb <- createWorkbook()

OR_results <- data.frame()
het_ex_results <- data.frame()
pleio_ex_results <- data.frame()
harmon_ex_results <- data.frame()
mr_presso_results <- data.frame()

addWorksheet(wb, "OR_results")
addWorksheet(wb, "Het_ex_results")
addWorksheet(wb, "Pleio_ex_results")
addWorksheet(wb, "MR_presso_ex_results")
addWorksheet(wb, "Harmon_ex")




column_names <- c("SNP", "effect_allele.exposure", "other_allele.exposure",  
                  "effect_allele.outcome", "other_allele.outcome", "beta.exposure", 
                  "beta.outcome", "eaf.exposure", "eaf.outcome", "remove", "palindromic", 
                  "ambiguous", "id.outcome", "chr.outcome", "pos.outcome", "se.outcome", 
                  "samplesize.outcome", "pval.outcome", "outcome", "originalname.outcome", 
                  "outcome.deprecated", "mr_keep.outcome", "data_source.outcome", "proxy.outcome", 
                  "target_snp.outcome", "proxy_snp.outcome", "target_a1.outcome", "target_a2.outcome", 
                  "proxy_a1.outcome", "proxy_a2.outcome", "samplesize.exposure", "chr.exposure", 
                  "pval.exposure", "se.exposure", "pos.exposure", "id.exposure", "pval_origin.outcome"  ,
                  "exposure", "mr_keep.exposure", "pval_origin.exposure", "data_source.exposure", 
                  "action", "mr_keep","ncase.outcome" ,"ncontrol.outcome")


harmon_ex_all <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(harmon_ex_all) <- column_names
str(harmon_ex_all)

for (file_name in file_list) {
 
  print(file_name)
  file_path <- file.path(folder_path, file_name)
  file_name_no_ext <- tools::file_path_sans_ext(file_name)
  
  if (file_name_no_ext==outcome_name){
    print("exp_out are the same,next!")
    next
    
  }
  
  out_dataSNP_proxy_final<-data.frame()

  exposure_name <-file_name_no_ext
  print(exposure_name)

  exp_dat <- read_exposure_data(
    filename = file_path,
    sep= ",",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col ="effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    pval_col = "pval.exposure",
    pos_col = "pos.exposure",
    id_col = "id.exposure",
    phenotype_col ="exposure",
    chr_col="chr.exposure",
    min_pval = 1e-200, 
    log_pval = FALSE, 
    samplesize_col = "samplesize.exposure",
    ncontrol_col = "ncontrol.exposure",
    ncase_col = "ncase.exposure"
  )
  
  if (!("eaf.exposure" %in% names(exp_dat))) {
    exp_dat$eaf.exposure <- NA
  }
  if (!("ncontrol.exposure" %in% names(exp_dat))) {
    exp_dat$ncontrol.exposure<- NA
  }
  if (!("ncase.exposure" %in% names(exp_dat))) {
    exp_dat$ncase.exposure <- NA
  }
  if (!("samplesize.exposure" %in% names(exp_dat))) {
    exp_dat$samplesize.exposure <- NA
  }
  
  exp_dat$exposure<-exposure_name
  exp_dat$id.exposure<-exposure_name
  head(exp_dat)

  desired_orderexp <- c("SNP", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", 
                        "beta.exposure", "se.exposure", "pval.exposure", "id.exposure", "chr.exposure", 
                        "samplesize.exposure", "exposure", "mr_keep.exposure", "pval_origin.exposure", 
                        "data_source.exposure", "ncontrol.exposure", "ncase.exposure")
  

  exp_dat <- exp_dat[, desired_orderexp]
  
  
  
  out_dataSNP_try <- try(
    format_data(
      out_data,
      type = "outcome",
      phenotype_col = "outcome",
      snps = exp_dat$SNP,
      header = TRUE,
      snp_col = "SNP",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      eaf_col = "eaf.outcome",
      effect_allele_col = "effect_allele.outcome",
      other_allele_col = "other_allele.outcome",
      pval_col = "pval.outcome",
      units_col = "units",
      ncase_col = "ncase.outcome",
      ncontrol_col = "ncontrol.outcome",
      samplesize_col = "samplesize.outcome",
      id_col = "id.outcome",
      min_pval = 1e-200,
      chr_col = "chr.outcome",
      pos_col = "pos.outcome",
      log_pval = FALSE), silent = F)
  

  if (inherits(out_dataSNP_try, "try-error")) {
    cat("An error occurred while processing", file_name, "\n")
    column_names <- c("SNP", "chr.outcome", "pos.outcome", "beta.outcome", "se.outcome", "samplesize.outcome", "pval.outcome", 
                      "eaf.outcome", "effect_allele.outcome", "other_allele.outcome", 
                      "outcome", "id.outcome", "originalname.outcome", "outcome.deprecated", "mr_keep.outcome", 
                      "data_source.outcome", "proxy.outcome", "target_snp.outcome", 
                      "proxy_snp.outcome", "target_a1.outcome", "target_a2.outcome", "proxy_a1.outcome", "proxy_a2.outcome")
    

    out_dataSNP <- data.frame(matrix(ncol = length(column_names), nrow = 0))
    colnames(out_dataSNP) <- column_names        
  } else {
    out_dataSNP <- out_dataSNP_try
  }
  
  snp_proxy <- as.data.table(unique(exp_dat$SNP)[!unique(exp_dat$SNP) %in% unique(out_dataSNP$SNP)])
  head(snp_proxy)
  
  
  
  if (nrow(snp_proxy) > 0) {
    print(paste0("there are ",nrow(snp_proxy), "SNPs to proxy..."))

    all_proxy <- data.table()
    

    for (snp1 in snp_proxy$V1) {

      a <- tryCatch({
        find_proxy(snp = snp1, r2_threshold = 0.8, build = "37", pop = "EUR")
      }, error = function(e) {

        print(paste("Error with SNP:", snp1, "; Error message:", e$message))
        NULL
      })
      
      if (!is.null(a)) {

        all_proxy <- rbind(all_proxy, a)
      }
    }
    head(all_proxy)
    
    if (nrow(all_proxy)>0){
      print("SNPs to proxy format in outcome")
      head(out_data)
      
      snp_proxy1<-unique(all_proxy$proxy_snp.outcome)
      out_dataSNP_proxy_try <- try(
        format_data(
          out_data,
          type = "outcome",
          phenotype_col = "outcome",
          snps = snp_proxy1,
          header = TRUE,
          snp_col = "SNP",
          beta_col = "beta.outcome",
          se_col = "se.outcome",
          eaf_col = "eaf.outcome",
          effect_allele_col = "effect_allele.outcome",
          other_allele_col = "other_allele.outcome",
          pval_col = "pval.outcome",
          units_col = "units",
          ncase_col = "ncase.outcome",
          ncontrol_col = "ncontrol.outcome",
          samplesize_col = "samplesize.outcome",
          id_col = "id.outcome",
          chr_col = "chr.outcome",
          pos_col = "pos.outcome",
          log_pval = FALSE), silent = F)
      

      if (inherits(out_dataSNP_proxy_try, "try-error")) {
        print(paste0("An error occurred while processing ", file_name))
        out_dataSNP_proxy <- data.frame()
      } else {
        out_dataSNP_proxy <- out_dataSNP_proxy_try
      }
      
      if (nrow(out_dataSNP_proxy)>0){
        print("SNPs to proxy format with exposure")
        

        out_proxy<-left_join(out_dataSNP_proxy,all_proxy,by=c("SNP"="proxy_snp.outcome"))
        
        head(out_proxy)
        
        out_proxy$proxy_snp.outcome<-out_proxy$SNP
        out_proxy$SNP<-out_proxy$target_snp.outcome
        out_proxy <- subset(out_proxy, (  
          (out_proxy$effect_allele.outcome == out_proxy$proxy_a1.outcome &  
             out_proxy$other_allele.outcome == out_proxy$proxy_a2.outcome) |  
            (out_proxy$effect_allele.outcome == out_proxy$proxy_a2.outcome &  
               out_proxy$other_allele.outcome == out_proxy$proxy_a1.outcome)    
        ))
        head(out_proxy)
        index1 <- which(out_proxy$effect_allele.outcome == out_proxy$proxy_a1.outcome & 
                          out_proxy$other_allele.outcome == out_proxy$proxy_a2.outcome)
        index2 <- which(out_proxy$effect_allele.outcome == out_proxy$proxy_a2.outcome & 
                          out_proxy$other_allele.outcome == out_proxy$proxy_a1.outcome)
        
        out_proxy$effect_allele.outcome[index1] <- out_proxy$target_a1.outcome[index1]
        out_proxy$other_allele.outcome[index1] <- out_proxy$target_a2.outcome[index1]
        
        out_proxy$effect_allele.outcome[index2] <- out_proxy$target_a2.outcome[index2]
        out_proxy$other_allele.outcome[index2] <- out_proxy$target_a1.outcome[index2]
        
        
        head(out_proxy)
        
        if (nrow(out_proxy)>0){
          
          out_dataSNP_proxy_origin_min <- as.data.frame(out_proxy%>%
                                                          group_by(target_snp.outcome) %>%
                                                          filter(row_number == min(row_number)))
          
          exp_dat_chr_pos<-dplyr::select(exp_dat,c(SNP,chr.exposure,pos.exposure))
          out_dataSNP_proxy_origin_min <- merge(out_dataSNP_proxy_origin_min, exp_dat_chr_pos, by.x="target_snp.outcome",by.y = "SNP", all.x = TRUE)  
          out_dataSNP_proxy_origin_min$chr.outcome <- NULL 
          out_dataSNP_proxy_origin_min$pos.outcome <- NULL
          names(out_dataSNP_proxy_origin_min)[names(out_dataSNP_proxy_origin_min) =="chr.exposure"] <-  "chr.outcome" 
          names(out_dataSNP_proxy_origin_min)[names(out_dataSNP_proxy_origin_min) == "pos.exposure" ] <- "pos.outcome"

          out_dataSNP_proxy_final<-out_dataSNP_proxy_origin_min
          out_dataSNP_proxy_final$originalname.outcome<-out_dataSNP_proxy_final$id.outcome
          out_dataSNP_proxy_final$outcome.deprecated<-paste0(unique(out_dataSNP_proxy_final$id.outcome),"||")
          out_dataSNP_proxy_final$data_source.outcome<-out_dataSNP_proxy_final$id.outcome
          out_dataSNP_proxy_final$proxy.outcome<-TRUE
          out_dataSNP_proxy_final<-subset(out_dataSNP_proxy_final, select = -c(R2,row_number))
        }
      }
    }
  }
  
  
  
  if(nrow(out_dataSNP)==0 & nrow(out_dataSNP_proxy_final)>0) {

    out_dataSNP <- setNames(data.frame(matrix(ncol = length(names(out_dataSNP_proxy_final)), nrow = 0)), names(out_dataSNP_proxy_final))
  } else if  (nrow(out_dataSNP)==0 & nrow(out_dataSNP_proxy_final) == 0)  { 
    print(paste0(file_name," no SNP in outcome"))
    next
  } else {
    out_dataSNP$originalname.outcome<-out_dataSNP$id.outcome
    out_dataSNP$outcome.deprecated<-paste0(unique(out_dataSNP$id.outcome),"||")
    out_dataSNP$data_source.outcome<-out_dataSNP$id.outcome
    out_dataSNP$proxy.outcome<-NA
    out_dataSNP$proxy_snp.outcome<-NA
    out_dataSNP$target_snp.outcome<-NA
    out_dataSNP$target_a1.outcome<-NA
    out_dataSNP$target_a2.outcome<-NA
    out_dataSNP$proxy_a1.outcome<-NA
    out_dataSNP$proxy_a2.outcome<-NA
  }
  
  if (nrow(out_dataSNP_proxy_final)>0) {
    
    out_data_all<-rbind(out_dataSNP,out_dataSNP_proxy_final)
    
  } else {
    
    out_data_all<-out_dataSNP
  }
  
  
  harmon_ex <- harmonise_data(
    exposure_dat=exp_dat,
    outcome_dat=out_data_all,
    action=3
  )
  
  
  if (nrow(harmon_ex) == 0 ) { 
    print(paste0("after harmonise ",file_name," no SNP"))
    next
  } else {
    missing_columns <- setdiff(names(harmon_ex_all), names(harmon_ex))
    for (col in missing_columns) {
      harmon_ex[[col]] <- NA
    }
    harmon_ex_all <- rbind(harmon_ex_all, harmon_ex)
  }
}


harmon_ex<-harmon_ex_all

################# run MR_presso############

# MR_presso directly uses harmon_ex to run. It cannot run with one SNP or even three SNPs. It requires at least 4 SNPs to avoid errors, so we need to filter the data first.

filtered_data_harmon_ex <- harmon_ex %>%
  group_by(exposure) %>%                    # Group by the exposure column
  filter(sum(mr_keep) >= 4) %>%             # Keep groups where the number of mr_keep rows that are True is greater than or equal to 4. Less than 4 won't work for running presso; Boolean values are automatically calculated.
  ungroup()

filtered_data_harmon_ex <- as.data.frame(filtered_data_harmon_ex)
# MR_presso directly runs, but it cannot handle one or three SNPs. It requires at least 4 SNPs to avoid errors.
presso <- run_mr_presso(filtered_data_harmon_ex, NbDistribution = 1000, SignifThreshold = 0.05)
pressoall <- presso
presso_list_num <- length(presso)
presso1 <- NULL
# First summarize presso[[i]]$`MR-PRESSO results`$`Global Test`$RSSobs, and presso[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue

for (i in 1:presso_list_num) {  
  
  outliers_indices <- presso[[i]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
  
  if (is.null(outliers_indices)) { outliers_indices=NA}
  
  distortion_coefficient <- presso[[i]]$`MR-PRESSO results`$`Distortion Test`$`Distortion Coefficient`  
  if (is.null(distortion_coefficient)) { distortion_coefficient=NA}
  pvalue <- presso[[i]]$`MR-PRESSO results`$`Distortion Test`$Pvalue  
  if (is.null(pvalue )) { pvalue =NA}
  rssobs <- presso[[i]]$`MR-PRESSO results`$`Global Test`$RSSobs  
  if (is.null(rssobs )) { rssobs =NA}
  pvalueglobal<-presso[[i]]$`MR-PRESSO results`$`Global Test`$Pvalue
  if (is.null(pvalueglobal)) { pvalueglobal=NA}
  

  if (is.character(outliers_indices) & any(outliers_indices == "No significant outliers")) {  

    outliers_length <- 0  
  } else if (is.na(outliers_indices[1])) {
    outliers_length <- 0  
  }else{  

    outliers_length <- length(outliers_indices)  
  }  
  
  

  presso2 <- c(rssobs,  
               pvalueglobal,   
               paste(as.character(outliers_indices), collapse = "_"),   
               outliers_length,   
               as.numeric(distortion_coefficient),   
               pvalue)  
  
  presso1<-rbind(presso1,presso2)
}



presso1<-as.data.frame(presso1)
presso1<-cbind(presso1,attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])


presso3<-NULL
for (i in 1:presso_list_num) {  
  presso4<-presso[[i]][["Main MR results"]]
  presso3<-rbind(presso3,presso4)
}
presso2<-cbind(attributes(presso)[["exposure"]],attributes(presso)[["outcome"]])
presso2<-as.data.frame(presso2)
presso2 <- presso2[rep(1:nrow(presso2), each = 2), ]

presso5<-cbind(presso2,presso3)
OR = exp(presso5$`Causal Estimate`)
OR_upper = exp(presso5$`Causal Estimate`+1.96*presso5$Sd)
OR_lower = exp(presso5$`Causal Estimate`-1.96*presso5$Sd)

presso6<-cbind(presso5,OR,OR_lower,OR_upper)

presso6<-dplyr::rename(presso6,
                       id.exposure=V1,
                       id.outcome=V2,
                       MR_analysis=`MR Analysis`,
                       `Causal_estimate`=`Causal Estimate`)
presso1<-dplyr::rename(presso1,
                       Global_Test_RSSobs=V1,
                       Global_Test_Pvalue=V2,
                       Outliers_indices=V3,
                       N_outliers=V4,
                       Distortion_Test_Pvalue  =V6,
                       Distortion_Test_Distortion_Coefficient=V5,
                       id.exposure=`attributes(presso)[["exposure"]]`,
                       id.outcome=`attributes(presso)[["outcome"]]`
                       
)


result_presso <- left_join(presso6, presso1, by = c("id.exposure", "id.outcome"))

result_presso$Exposure<-NULL



result_presso<-dplyr::rename(result_presso,
                             Exposure=id.exposure,
                             Outcome=id.outcome
                             
)


head(result_presso,n=2)
mr_presso_ex_results<-result_presso
writeData(wb, "MR_presso_ex_results", mr_presso_ex_results)




####################run 2SMR######################

split_data <- split(harmon_ex, harmon_ex$exposure)



lapply(names(split_data), function(x) {

  valid_name <- make.names(x, unique = TRUE)

  file_name <- paste0("exposure_", valid_name, ".csv")
  file_path <- file.path(output_directory, file_name)
  write.csv(split_data[[x]], file_path, row.names = FALSE)
})


folder_path <- output_directory
file_list <- list.files(folder_path, pattern = "\\.csv$")

res_ex <- mr(harmon_ex, method_list=c("mr_wald_ratio",
                                      'mr_ivw_mre', 
                                      'mr_ivw_fe', 
                                      "mr_ivw",
                                      'mr_egger_regression', 
                                      'mr_weighted_median', 
                                      'mr_simple_mode', 
                                      'mr_weighted_mode'))

OR <- generate_odds_ratios(res_ex)
het_ex <- mr_heterogeneity(harmon_ex)
pleio_ex <- mr_pleiotropy_test(harmon_ex)



cols <- c( "exposure" , "outcome" , "method" ,
           "nsnp" , "b"  ,"se", "pval" , "lo_ci"  ,
           "up_ci"  , "or"   , "or_lci95"  ,  "or_uci95" )
new_cols <- c( "Exposure" ,"Outcome" ,  "Method" ,
               "Nsnp"   ,  "Beta"  ,  "SE", "P" , "LowerCI"  ,
               "UpperCI"  ,"OR", "ORlower", "ORupper" )
OR_1 <-OR[, cols]
names(OR_1) <- new_cols
OR_results <- OR_1


cols <- c(  "exposure" ,  "outcome" ,"method", "Q" ,  "Q_df" ,  "Q_pval" )
new_cols <-c(  "Exposure" ,  "Outcome" , "Method", "Q_statistic", "Degrees_of_freedom_(df)_for_Q", "P-value_for_Q" )
het_ex_1 <-het_ex[, cols]
names(het_ex_1) <- new_cols
het_ex_results <- het_ex_1

cols <- c( "exposure","outcome" ,"egger_intercept", "se","pval" )
new_cols <-c( "Exposure", "Outcome", "Egger_intercept", "SE", "P-value")
pleio_ex_1 <-pleio_ex[, cols]
names(pleio_ex_1) <- new_cols
pleio_ex_results <- pleio_ex_1

harmon_ex_results <- harmon_ex



writeData(wb, "OR_results", OR_results)
writeData(wb, "Het_ex_results", het_ex_results)
writeData(wb, "Pleio_ex_results", pleio_ex_results)
writeData(wb, "Harmon_ex", harmon_ex_results)


resultfolderpath1<-paste0(bigpath,"/result/tsmr_result/", last_field)
dir.create(resultfolderpath1, recursive = TRUE)
filename1ast <- paste0(resultfolderpath1,"/",last_field,"_results_", outcome_name, ".xlsx")
saveWorkbook(wb, filename1ast, overwrite = TRUE)








