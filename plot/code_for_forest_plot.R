source("F:/Research_Skills_and_Resources/lc_code_package/change_name.R")

library(ggview)
plot_forest_lc_beta<-function(data_o,
                              out_plot_name)
{
  data_o<-data_o %>% 
    dplyr::select(-c(Method,Nsnp,SE , OR, ORlower,ORupper))  %>% 
    dplyr::rename(x= Exposure   ,
                  pvalue=P,
                  med=Beta,
                  min = LowerCI     ,
                  max =UpperCI
    )
  data_o$x
  head(data_o)
  
  data_o<-replace_strings(data_frame = data_o, 
                          rules_txt = "", 
                          columns = c("x"), 
                          whole_cell = T)
  
  group1_vars <- c("Endometriosis", "GH",  "Gestational diabetes", "Infertility", "Preterm delivery", 
                   "Post-term delivery",
                   "Pre-eclampsia", "Ectopic pregnancy", "POF", "Ovarian dysfunction", 
                   "UL", "EVP","Spontaneous abortion", "Habitual aborter", 
                   "Eclampsia", "Adenomyosis", "Menorrhagia", "Irregular menses", 
                   "Oligomenorrhoea", "PCOS")
  
  group2_vars <- c("AMH", "Estrogen", "BAT", "SHBG", "Total testosterone", 
                   "FSH", "LH", "Prolactin", "Progesterone")
  
  
  group3_vars <- c("AFB", "AFS", "ALB", "Years ovulating", "NLB", "Gestational duration",  "AAM", "ANM")
  
  

  data_o$group <- ifelse(data_o$x %in% group1_vars, "group1",
                         ifelse(data_o$x %in% group2_vars, "group2",
                                ifelse(data_o$x %in% group3_vars, "group3", NA)))
  
  #data_o$x <- factor(data_o$x, levels = c(group1_vars, group2_vars, group3_vars))

  head(data_o)
  color_mapping <- c("group1" = "#3D85BC", "group2" = "#F98F34", "group3" = "#CD3A39")
  

  data_o$group_col <- color_mapping[data_o$group]
  
  
  data_o$p <- ifelse(data_o$pvalue > 0.05, "",
                     ifelse(data_o$pvalue > ( 0.05/37) & data_o$pvalue <= 0.05, "*",
                            ifelse(data_o$pvalue > ( 0.05/(37*6)) & data_o$pvalue <= ( 0.05/37), "**",
                                   ifelse(data_o$pvalue <= ( 0.05/(37*6)), "***", ""))))
  

  head(data_o)
  unique(data_o$p)
  
  data_o$p_col <- NA
  
  data_o$p_col[which(data_o$p != "" & data_o$med > 0)] <- "Postive effect(P<0.05)"
  data_o$p_col[which(data_o$p == "" & data_o$med > 0)] <- "Postive effect(P>=0.05)"
  data_o$p_col[which(data_o$p != "" & data_o$med <= 0)] <- "Negtive effect(P<0.05)"
  data_o$p_col[which(data_o$p == "" & data_o$med <= 0)] <- "Negtive effect(P>=0.05)"
  
  
  group1_vars_sorted <- sort(group1_vars)
  group2_vars_sorted <- sort(group2_vars)
  group3_vars_sorted <- sort(group3_vars)
  
  

  group1_vars_filtered <- intersect(group1_vars_sorted, data_o$x)
  group2_vars_filtered <- intersect(group2_vars_sorted, data_o$x)
  group3_vars_filtered <- intersect(group3_vars_sorted, data_o$x)
  
  

  
  data_o <- data_o %>%
    dplyr::mutate(x = factor(x, levels = c(group1_vars_filtered, group2_vars_filtered,group3_vars_filtered))) %>%
    dplyr::arrange(x)
  
  
  levels(data_o$x)
  
  
  
  head(data_o)
  # Create a new data frame to hold arrow information
  data_arrows1 <- data_o %>%
    mutate(arrow_y1 = ifelse(min < -1,-0.92, NA),
           
           arrow_x = ifelse(min < -1 | max > 1, x, NA)) %>%
    filter(!is.na(arrow_y1))
  
  
  
  data_arrows2 <- data_o %>%
    mutate(arrow_y2 = ifelse(max > 1, 0.92, NA),
           arrow_x = ifelse(min < -1 | max > 1, x, NA)) %>%
    filter(!is.na(arrow_y2))
  
  head(data_o)
  
  length1 <- length(group1_vars_filtered)
  length2 <- length(group2_vars_filtered)
  length3 <- length(group3_vars_filtered)
  
  
  fill_color1 <- ifelse(length1==0 & length2==0, "#FFF1E9", 
                        ifelse(length1==0, "#FDFDF5", "#E5EBF0"))
  
  fill_color2 <- ifelse(length3==0 & length2==0 , "#E5EBF0", 
                        ifelse(length3==0 , "#FDFDF5", "#FFF1E9"))
  
  
  p <- ggplot(data_o) +
    annotate("rect", xmin = 0.4, xmax = 0.5, ymin = -1, ymax = 1, alpha = 0.8, fill = fill_color1) +
    annotate("rect", xmin = 0.5  + length1 + length2 + length3, xmax =  0.6 + length1 + length2 + length3, ymin = -1, ymax = 1, alpha = 0.8, fill = fill_color2)+
    annotate("rect", xmin = c(0.5, 
                              0.5 + length1, 
                              0.5 + length1 + length2), 
             xmax = c(0.5 + length1, 
                      0.5 + length1 + length2, 
                      0.5 + length1 + length2 + length3), 
             ymin = -1, ymax = 1, alpha = 0.8, fill = c("#E5EBF0", "#FDFDF5", "#FFF1E9")) +
    
    geom_hline(yintercept = 0, linewidth = 0.3) +

    geom_linerange(aes(x, ymin = ifelse(min < -1, -1, min), ymax = ifelse(max > 1, 1, max), color = p_col, 
                       size = factor(p_col, levels = c("Postive effect(P<0.05)", "Postive effect(P>=0.05)", "Negtive effect(P<0.05)", "Negtive effect(P>=0.05)")),
                       alpha = factor(p_col, levels = c("Postive effect(P<0.05)", "Postive effect(P>=0.05)", "Negtive effect(P<0.05)", "Negtive effect(P>=0.05)"))), show.legend = F) +
    # Scatter points
    geom_point(aes(x, med, 
                   alpha = factor(p_col, levels = c("Postive effect(P<0.05)", "Postive effect(P>=0.05)", "Negtive effect(P<0.05)", "Negtive effect(P>=0.05)")),
                   color = p_col), size = 1) +
    # Significance
    geom_text(aes(x, y = ifelse(max > 0.8, min - 0.15, max + 0.15), label = p, color = p_col), size = 10,size.unit = "pt", show.legend = F) +
    # Arrows
    geom_segment(data = data_arrows1, aes(x = arrow_x, xend = arrow_x, y = arrow_y1, yend = arrow_y1 + sign(arrow_y1)*0.07, color =p_col,
                                          size = factor(p_col, levels = c("Postive effect(P<0.05)", "Postive effect(P>=0.05)", "Negtive effect(P<0.05)", "Negtive effect(P>=0.05)")),
                                          alpha = factor(p_col, levels = c("Postive effect(P<0.05)", "Postive effect(P>=0.05)", "Negtive effect(P<0.05)", "Negtive effect(P>=0.05)"))), 
                 arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
    
    geom_segment(data = data_arrows2, aes(x = arrow_x, xend = arrow_x, y = arrow_y2, yend = arrow_y2 + sign(arrow_y2)*0.07, color =p_col,
                                          size = factor(p_col, levels = c("Postive effect(P<0.05)", "Postive effect(P>=0.05)", "Negtive effect(P<0.05)", "Negtive effect(P>=0.05)")),
                                          alpha = factor(p_col, levels = c("Postive effect(P<0.05)", "Postive effect(P>=0.05)", "Negtive effect(P<0.05)", "Negtive effect(P>=0.05)"))), 
                 arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
    
    

    scale_color_manual(name = "",
                       values = c("Postive effect(P<0.05)" = "#FF4500",
                                  "Postive effect(P>=0.05)" = "#B2B2B2",
                                  "Negtive effect(P<0.05)" = "#1E90FF",
                                  "Negtive effect(P>=0.05)" = "#B2B2B2"), guide = "none") +

    scale_size_manual(values = c("Postive effect(P<0.05)" = 0.4,
                                 "Postive effect(P>=0.05)" = 0.4,
                                 "Negtive effect(P<0.05)" = 0.4,
                                 "Negtive effect(P>=0.05)" = 0.4), guide = "none") +

    scale_alpha_manual(values = c("Postive effect(P<0.05)" = 1,
                                  "Postive effect(P>=0.05)" = 1,
                                  "Negtive effect(P<0.05)" = 1,
                                  "Negtive effect(P>=0.05)" = 1), guide = "none") +

    scale_y_continuous(expand = c(0, 0)) +
    xlab("") +
    ylab(paste0("Beta for ", out_plot_name)) +
    theme_bw() 
  
  p_noytext <- p +   theme(axis.text.y =  element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.text.x = element_text(size =6),#	text size in pts.
                           axis.title.x = element_text(size = 6),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank()) +
    coord_flip(ylim = c(-1, 1)) +
    theme(aspect.ratio = 1.9)
  
  print(p_noytext)
  
  
  ggsave(paste0("tsmr_", Outcome_name, "_forest_no_ytext.pdf"), plot = p_noytext, width = 60, height = 100, units = "mm")
  
  
  p_ytext <- p +  theme(axis.text.y = element_text(color = data_o$group_col, size = 6),
                        axis.text.x = element_text(size = 6),
                        axis.title.x = element_text(size = 6),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank()) +
    coord_flip(ylim = c(-1, 1)) +
    theme(aspect.ratio = 1.9)
  
  print(p_ytext)
  ggview(p_noytext, width = 60, height = 100, units = "mm")
  ggview(p_ytext, width = 120, height = 100, units = "mm")
  # 
  ggsave(paste0("tsmr_",Outcome_name,"_forest.pdf"),  plot = p_ytext,width = 120, height = 100, units = "mm")
  
  
}

