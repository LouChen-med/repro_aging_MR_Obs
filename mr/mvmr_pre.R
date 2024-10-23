
setwd("")
args = commandArgs(TRUE)
exp1 = as.character(args[1])
exp2 = as.character(args[2])
out = as.character(args[3])

library(tidyverse)
library(data.table)
library(meta)
library(devtools)
library(pacman)
library(TwoSampleMR)
library(ieugwasr)
library(phenoscanner)
library(LDlinkR)
library(mr.raps)
library(MRPRESSO)
library(extrafont)
library(writexl)
library(MVMR)
library(MendelianRandomization)
library(openxlsx)
library(dplyr)
library(coloc)
library(magrittr)
library(rlang)
library(stringr)
library(vroom)
extract_part <- function(path) {
  match <- sub(".*f_(.*)\\.txt\\.gz.*", "\\1", path)
  return(match)
}

name_exp1 <- extract_part(exp1)
name_exp2 <- extract_part(exp2)
name_out <- extract_part(out)


print(name_exp1)
print(name_exp2)
print(name_out)


clump_local<-function (data_filter_pval,clump_r2 = 0.001, clump_kb = 10000,clump_p1 =5e-8){
  print(paste0(nrow(data_filter_pval)," SNP after pval.exposure < 5e-8,independent SNPs might >2"))

  data_filter_pval_to_clump <- dplyr::rename(
    data_filter_pval,
    rsid = SNP,
    pval = pval.exposure,
    id = id.exposure
  )
  

  tryCatch({
    data_p_clump <- ld_clump(
      data_filter_pval_to_clump, 
      clump_kb = clump_kb,
      clump_r2 = clump_r2,
      clump_p =clump_p1,
      pop = "EUR",
      plink_bin = "/home/louchen/plink_software/plink",
      bfile = "/home/louchen/Genome_reference_files/ieu_1000g_maf0_01/EUR"
    )
  }, error = function(e) {
    message("Clumping failed or no significant SNPs: ", e$message)

    data_p_clump <<- setNames(data.frame(matrix(ncol = length(names(data_filter_pval_to_clump)), nrow = 0)), names(data_filter_pval_to_clump))
  })
  
  data_final <- dplyr::rename(
    data_p_clump,
    SNP = rsid,
    pval.exposure = pval,
    id.exposure = id
  )
  print(paste0("After clump ",nrow(data_final)," SNPs,independent SNPs on 5e-8"))
  return(data_final)
}


mv_extract_exposures_local_lc<-function (filenames_exposure, sep = " ", phenotype_col = "Phenotype", 
                                         snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf", 
                                         effect_allele_col = "effect_allele", other_allele_col = "other_allele", 
                                         pval_col = "pval", units_col = "units", ncase_col = "ncase", 
                                         ncontrol_col = "ncontrol", samplesize_col = "samplesize", 
                                         gene_col = "gene", id_col = "id", min_pval = 1e-200, log_pval = FALSE, 
                                         pval_threshold = 5e-08, clump_r2 = 0.001, clump_kb = 10000, 
                                         harmonise_strictness = 2) 
{
  print("SNP: MVMR_extract")
  l_full <- list()
  l_inst <- list()

  for (i in 1:length(filenames_exposure)) {
    l_full[[i]] <- fread(filenames_exposure[i])
    l_full[[i]]<-as.data.frame(l_full[[i]])
    name<-unique(l_full[[i]]$exposure)
    print(paste0("now_extract_for:",name))
    print(nrow(l_full[[i]]))
    l_full[[i]] <- subset(l_full[[i]], !is.na(SNP) & SNP != "")
    print(nrow(l_full[[i]]))
    l_full[[i]] <- subset(l_full[[i]], !is.na(pval.exposure) & pval.exposure != "")
    print(nrow(l_full[[i]]))
    l_full[[i]] <- subset(l_full[[i]], pval.exposure >= 0 & pval.exposure <= 1)

    l_inst[[i]] <- subset(l_full[[i]], pval.exposure < pval_threshold)

    names(l_full[[i]]) <- gsub("exposure","outcome",names(l_full[[i]]))
    
    names(l_inst[[i]]) <- gsub("outcome","exposure",names(l_inst[[i]]))
    
    l_inst[[i]] <- subset(l_inst[[i]], pval.exposure < pval_threshold)
    if (nrow(l_inst[[i]])==0) {
      print(paste0("there is no significant SNP for : ",name))
    }

    if (nrow(l_inst[[i]])>0) {
      l_inst[[i]] <- clump_local(l_inst[[i]], clump_p1 = pval_threshold, clump_r2 = clump_r2, clump_kb = clump_kb)
    }
    print(head(l_inst[[i]]))
    nameall<-c("other_allele.exposure" , "effect_allele.exposure", "SNP" ,"pval.exposure","beta.exposure" ,"se.exposure",           
               "eaf.exposure", "samplesize.exposure" , "exposure" ,"mr_keep.exposure", "pval_origin.exposure", "id.exposure" ,  
               "ncase.exposure","ncontrol.exposure","chr.exposure","pos.exposure" )
    
    missing_columns <- setdiff(nameall, names(l_inst[[i]]))
    print(names(l_inst[[i]]))
    print(missing_columns)
    for (col in missing_columns) {
      l_inst[[i]][[col]] <- NA
    }
    print(head(l_inst[[i]]))
  }

  exposure_dat <-  do.call(rbind, l_inst)
  
  if(nrow(exposure_dat) == 0) {
    print("after bind ,there is no SNP")
    return(NULL)
  }

  temp <- exposure_dat
  temp$id.exposure <- 1
  temp <- temp[order(temp$pval.exposure, decreasing = FALSE), ]
  temp <- subset(temp, !duplicated(SNP))
  if(nrow(temp) == 0) {
    print("after !duplicated(SNP) for temp,there is no SNP")
    return(NULL)
  }
  
  temp <-  clump_local(temp, clump_p1 = pval_threshold, clump_r2 = clump_r2, 
                       clump_kb = clump_kb)
  if(nrow(temp) == 0) {
    print("after clump for temp,there is no SNP")
    return(NULL)
  }
  
  

  exposure_dat <- subset(exposure_dat, SNP %in% temp$SNP)

  d1 <- lapply(l_full, function(x) {
    subset(x, SNP %in% exposure_dat$SNP)
  }) %>% dplyr::bind_rows()
  
  id_exposure <- unique(d1$id.outcome)

  #stopifnot(length(unique(d1$id)) == length(unique(id_exposure)))
  
  if(!length(unique(d1$id.outcome)) == length(unique(filenames_exposure))) {
    print("after extract,the snps are not in all gwas")
    return(NULL)
  }

  d2 <- subset(d1, id.outcome != id_exposure[1])
  d1 <- convert_outcome_to_exposure(subset(d1, id.outcome == id_exposure[1]))
  d <- harmonise_data(d1, d2, action = harmonise_strictness)
 
  tab <- table(d$SNP)
  keepsnps <- names(tab)[tab == length(id_exposure) - 1]
  
  if(length(keepsnps) == 0) {
    print("after harmonise,the snps are not in all gwas")
    return(NULL)
  }
  
  d <- subset(d, SNP %in% keepsnps)
  
  dh1 <- subset(d, id.outcome == id.outcome[1], select = c(SNP, exposure, id.exposure, 
                                                           effect_allele.exposure, other_allele.exposure, 
                                                           eaf.exposure, beta.exposure, se.exposure, 
                                                           pval.exposure))
  dh2 <- subset(d, select = c(SNP, outcome, id.outcome, effect_allele.outcome, 
                              other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, 
                              pval.outcome))
  names(dh2) <- gsub("outcome", "exposure", names(dh2))
  dh <- rbind(dh1, dh2)
  return(dh)
}


exposure_dat <- mv_extract_exposures_local_lc(
  filenames_exposure = c(exp1, exp2), 
  sep = "\t",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col ="effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",
  pval_col = "pval.exposure",
  id_col = "id.exposure",
  phenotype_col ="exposure",
  samplesize_col = "samplesize.exposure",
  ncontrol_col = "ncontrol.exposure",
  ncase_col = "ncase.exposure",
  min_pval = 1e-200, 
  log_pval = FALSE, 
  pval_threshold = 5e-8, 
  clump_r2 = 0.001,
  clump_kb = 10000,
  harmonise_strictness = 3
)

if (is.null(exposure_dat)) {
  stop("no SNP for exp1 and exp2")
}


print(paste0("Number of IVs for exps: ", as.character(length(exposure_dat$SNP)/2)))

out_data<-fread(out)
names(out_data) <- gsub("exposure", "outcome", names(out_data))
out_data$effect_allele.outcome <- toupper(out_data$effect_allele.outcome)
out_data$other_allele.outcome <- toupper(out_data$other_allele.outcome)
outcome_name<-sub(".*/f_(.*).txt(.*$)?", "\\1", out)
out_data<-as.data.frame(out_data)
print(head(out_data))
print(str(out_data))



out_dataSNP_try <- try(
  format_data(
    out_data,
    type = "outcome",
    phenotype_col = "outcome",
    snps = exposure_dat$SNP,
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

outcome_dat<-out_dataSNP

print(paste0("Number of IVs for out: ", nrow(outcome_dat)))

if (nrow(outcome_dat) == 0) {
  stop("If there are no matching SNPs in the outcome, it will cause an error in MV-harmonise, terminating the execution.")
}


# Harmonising exposure and outcome datasets
dat <- tryCatch({

  mv_harmonise_data(
    exposure_dat = exposure_dat, 
    outcome_dat = outcome_dat, 
    harmonise_strictness = 3
  )
}, error = function(e) {

  message("An error occurred during the execution of harmonise_data: insufficient valid SNPs.", e$message)
  stop("stop")
})


if (length(dat[["exposure_beta"]][, 1]) == 0) {
  stop("结局中没有匹配的SNP，终止执行。")
}


# Get the names of the first and second columns. This is not based on the input order of exp1 and exp2, but rather sorted alphabetically by the first letter. For example, if Menarche and CAD are input, the order becomes CAD, Menarche.

names1 <- colnames(dat[["exposure_beta"]])[1]  
names2 <- colnames(dat[["exposure_beta"]])[2]  
outname1 <-dat["outname"]$outname$outcome

str(dat)

mvmr_dat <- data.frame(
  rsid = rownames(dat[["exposure_beta"]]),
  name1_beta = dat[["exposure_beta"]][, 1],
  name1_se = dat[["exposure_se"]][, 1], 
  name1_p = dat[["exposure_pval"]][, 1], 
  name2_beta = dat[["exposure_beta"]][, 2],
  name2_se = dat[["exposure_se"]][, 2],
  name2_p = dat[["exposure_pval"]][, 2],
  outname_beta = dat[["outcome_beta"]],
  outname_se = dat[["outcome_se"]],
  outname_p = dat[["outcome_pval"]]
)
mvmr_dat$SNP<-mvmr_dat$rsid
if (names1==name_exp1) {  
  mvmr_dat$Beta.x<-mvmr_dat$name1_beta
  mvmr_dat$SE.x<-mvmr_dat$name1_se
  mvmr_dat$Beta.y<-mvmr_dat$name2_beta
  mvmr_dat$SE.y<-mvmr_dat$name2_se
}
if (names2==name_exp1) {   
  mvmr_dat$Beta.x<-mvmr_dat$name2_beta
  mvmr_dat$SE.x<-mvmr_dat$name2_se
  mvmr_dat$Beta.y<-mvmr_dat$name1_beta
  mvmr_dat$SE.y<-mvmr_dat$name1_se
}

mvmr_dat$Beta.o<-mvmr_dat$outname_beta
mvmr_dat$SE.o<-mvmr_dat$outname_se

if (names1==name_exp1) {
  colnames(mvmr_dat)[2:4] <- c(paste0("exp1_",names1, "_beta"), paste0("exp1_",names1, "_se"), paste0("exp1_",names1, "_p"))
  colnames(mvmr_dat)[5:7] <- c(paste0("exp2_",names2, "_beta"), paste0("exp2_",names2, "_se"), paste0("exp2_",names2, "_p"))
}
if (names2==name_exp1) { 
  colnames(mvmr_dat)[2:4] <- c(paste0("exp2_",names1, "_beta"), paste0("exp2_",names1, "_se"), paste0("exp2_",names1, "_p"))
  colnames(mvmr_dat)[5:7] <- c(paste0("exp1_",names2, "_beta"), paste0("exp1_",names2, "_se"), paste0("exp1_",names2, "_p"))
}
colnames(mvmr_dat)[8:10] <- c(paste0("out_",outname1, "_beta"), paste0("out_",outname1, "_se"), paste0("out_",outname1, "_p"))

rownames(mvmr_dat) <- c()
head(mvmr_dat)


folderpathXYO<-paste0("/MVMR/X_Y_O_merge_pre","/","EXP_",name_exp1,"/","OUT_",outname1)


dir.create(folderpathXYO, recursive = TRUE)

filenameout <- paste0(folderpathXYO,"/O_",outname1,"_",paste0("X_", name_exp1, "_Y_", name_exp2, ".csv"))
filenameout

write.csv(mvmr_dat, file = filenameout,row.names = FALSE)

