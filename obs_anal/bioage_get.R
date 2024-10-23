library(BioAge)
library(dplyr)
library(readr)
library(foreign)
library(data.table)
library(stringr)
library(dplyr)
library(sas7bdat)
library(data.table)


######DNAm data######
data <- read.sas7bdat("/dnmepi.sas7bdat")
gc()

data_epiage<-dplyr::select(data,c(
  SEQN,HorvathAge,HannumAge,SkinBloodAge,PhenoAge,GrimAgeMort,GrimAge2Mort,WTDN4YR,
  CD8TPP,CD4TPP,Nkcell, Bcell,MonoPP,NeuPP
))
data_epiage<-dplyr::rename(data_epiage,
                           seqn=SEQN,
                           DNAmHorvathAge=HorvathAge,
                           DNAmHannumAge=HannumAge,
                           DNAmSkinBloodAge=SkinBloodAge,
                           DNAmPhenoAge=PhenoAge,
                           DNAmGrimAge=GrimAgeMort,
                           DNAmGrimAge2=GrimAge2Mort
)
head(data_epiage)

data_epiage[data_epiage == "NaN"] <- NA
head(data_epiage)

fwrite(data_epiage,"",sep="\t",quote=F)

#####DNAmageAccel#####


library(rlang)
library(dplyr)

DNAmAge_list <- c("DNAmHorvathAge", "DNAmHannumAge",  "DNAmPhenoAge",
                  "DNAmSkinBloodAge", "DNAmGrimAge", "DNAmGrimAge2","IEAA")


for (DNAmAge in DNAmAge_list) {
  
  
  print(paste0("now run from DNAmAge ",trait_name))
  
  
  DNAmAgeAccel <- ifelse(DNAmAge == "IEAA", "IEAA", paste0(DNAmAge, "Accel"))
  
  print(paste0("now run from ",DNAmAgeAccel))
  

  if (DNAmAge == "IEAA") {
    data_analysis_DNAm <- filter(data_analysis, !is.na(DNAmHorvathAge))
  } else {
    data_analysis_DNAm <- filter(data_analysis, !is.na(!!sym(DNAmAge)))
  }
  

  if (DNAmAge == "IEAA") {

    DNAm_model <- lm(DNAmHorvathAge ~ ageyr +  CD8TPP +  CD4TPP+ Nkcell + Bcell + MonoPP +NeuPP
                     , data = data_analysis_DNAm)
    

    data_analysis_DNAm[[DNAmAgeAccel]] <- residuals(DNAm_model)
  } else {

    DNAm_model <- lm(as.formula(paste(DNAmAge, "~ ageyr")), data = data_analysis_DNAm)
    

    data_analysis_DNAm[[DNAmAgeAccel]] <- residuals(DNAm_model)
  }
  
  
 

}



#####bioage data######

#specify bioamarkers included in the algorithms 
biomarkers_kdm = c("sbp","totchol","hba1c","albumin","creat","lncrp","alp","bun")

biomarkers_phenoage = c("albumin_gL","lymph","mcv","glucose_mmol",
                        "rdw","creat_umol","lncrp","alp","wbc")

#KDM bioage using NHANES (separate training for men and women)
kdm = kdm_nhanes(biomarkers_kdm)
#phenoage uinsg NHANES
phenoage = phenoage_nhanes(biomarkers_phenoage)

kdm_data = as.data.frame(kdm[["data"]])
kdm_data<-subset(kdm_data,select = c(sampleID, year,kdm))
kdm_data<-dplyr::rename(kdm_data,
                        seqn=sampleID,
                        Year=year,
                        BioKDMAge=kdm
)
kdm_data <- kdm_data %>%
  mutate(seqn = sub(".*_", "", seqn))
kdm_data1<-subset(kdm_data,!is.na(kdm_data$BioKDMAge))
unique(kdm_data1$Year)

kdm_data1_filtered <- kdm_data1 %>%
  filter(!Year %in% c(2011, 2013))

unique(kdm_data1_filtered$Year)
head(kdm_data1_filtered)

phenoage_data = as.data.frame(phenoage[["data"]])
head(phenoage_data)
phenoage_data<-subset(phenoage_data,select = c(sampleID, year,phenoage))
phenoage_data<-dplyr::rename(phenoage_data,
                             seqn=sampleID,
                             Year=year,
                             BioPhenoAge=phenoage
)
phenoage_data <- phenoage_data %>%
  mutate(seqn = sub(".*_", "", seqn))

phenoage_data1<-subset(phenoage_data,!is.na(phenoage_data$BioPhenoAge))
unique(phenoage_data1$Year)


all_data <- full_join(kdm_data1_filtered, phenoage_data1, by = "seqn")



fwrite(all_data,".txt",sep="\t",quote=F)

















