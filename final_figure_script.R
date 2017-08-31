library(tidyverse)
library(ggplot2)
library(cowplot)

######## Figure 1: Unique QTL plot ########
load("~/Dropbox/AndersenLab/QTLpaper/data/uniqueQTL_588corrcut.Rda")

textsize = 10
titlesize = 12

threetrait <- uniqueQTL %>%
  dplyr::ungroup() %>%
  dplyr::group_by(condition, chr) %>%
  arrange(condition, chr, ci_l_pos, ci_r_pos) %>%
  na.omit() %>%
  distinct(chr, trait, condition, .keep_all = T) %>%
  mutate(class = NA) %>%
  select(-unique, -chrdrugtrait, -drugtrait)

#Add drug class
for(i in 1:nrow(threetrait)) {
  if(threetrait$condition[i] %in% c("cadmium", "silver", "cisplatin", "copper")) {
    threetrait$class[i] <- "Heavy metal"
  } else if(threetrait$condition[i] %in% c("chlorpyrifos", "deiquat", "paraquat")) {
    threetrait$class[i] <- "Pesticide"
  } else if(threetrait$condition[i] == "fluoxetine") {
    threetrait$class[i] <- "Neuroactive"
  } else if(threetrait$condition[i] %in% c("carmustine", "chlorothanilil", "irinotecan", "mechlorethamine", "topotecan", "tunicamycin", "vincristine", "FUdR")) {
    threetrait$class[i] <- "Chemotherapeutic"
  }
}

#Factor drug class
threetrait$class <- factor(threetrait$class)
threetrait <- threetrait %>%
  arrange(class, condition)
threetrait$condition <- factor(threetrait$condition, levels = c("carmustine", "chlorothanilil", "FUdR", "irinotecan", "mechlorethamine", "topotecan", "tunicamycin", "vincristine", 
                                                                "cadmium", "cisplatin", "copper", "silver", "chlorpyrifos", "deiquat", "paraquat", "fluoxetine"),
                               labels = c("Carmustine", "Chlorothanilil", "FUdR", "Irinotecan", "Mechlorethamine", "Topotecan", "Tunicamycin", "Vincristine", 
                                          "Cadmium", "Cisplatin", "Copper", "Silver", "Chlorpyrifos", "Deiquat", "Paraquat", "Fluoxetine"))

#Set chromosome boundaries
newrows <- threetrait[1,] 
newrows[1,] = c(NA,"I",5000000,"q25.EXT",0,NA,NA,NA,NA,NA,1,NA,14972282,"Cadmium", "Heavy metal")
newrows[2,] = c(NA,"II",5000000,"q25.EXT",0,NA,NA,NA,NA,NA,1,NA,15173999,"Cadmium","Heavy metal")
newrows[3,] = c(NA,"III",5000000,"q25.EXT",0,NA,NA,NA,NA,NA,1,NA,13829314,"Cadmium","Heavy metal")
newrows[4,] = c(NA,"IV",5000000,"q25.EXT",0,NA,NA,NA,NA,NA,1,NA,17450860,"Cadmium","Heavy metal")
newrows[5,] = c(NA,"V",5000000,"q25.EXT",0,NA,NA,NA,NA,NA,1,NA,20914693,"Cadmium","Heavy metal")
newrows[6,] = c(NA,"X",5000000,"q25.EXT",0,NA,NA,NA,NA,NA,1,NA,17748731,"Cadmium","Heavy metal")
newrows$ci_l_pos <- as.numeric(newrows$ci_l_pos)
newrows$ci_r_pos <- as.numeric(newrows$ci_r_pos)
newrows$pos <- as.numeric(newrows$pos)
newrows$lod <- as.numeric(newrows$lod)

#Plot
unique_all <- ggplot(threetrait)+
  aes(x=pos/1E6, y=trait)+
  theme_bw() +
  viridis::scale_fill_viridis(name = "LOD") + viridis::scale_color_viridis(name = "LOD") +
  geom_segment(aes(x = ci_l_pos/1e6, y = trait, xend = ci_r_pos/1e6, yend = trait, color = lod), size = 2, alpha = 1) +
  geom_segment(data=newrows,aes(x = 0, y = trait, xend = ci_r_pos/1e6, yend = trait), size = 2.5, alpha = 0) +
  geom_point(aes(fill=lod),colour = "black",size = 2, alpha = 1, shape = 21)+
  xlab("Genomic position (Mb)") + ylab("") +
  theme(axis.text.x = element_text(size=textsize, face="bold", color="black"),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = titlesize, face = "bold"), legend.text = element_text(size = textsize),
        legend.key.size = unit(.75, "cm"),
        panel.grid.major.x = element_line(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=titlesize, face="bold", color= "black"),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size=titlesize, face="bold", color="black"),
        strip.text.y = element_text(size=titlesize, face="bold", color="black", angle = 0),
        strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
        plot.title = element_text(size=titlesize, face="bold")) +
  facet_grid(condition ~ chr, scales = "free_x", space = "free")

#make into a Grob to find the coordinates of where we want to add colored boxes
elements <- ggplotGrob(unique_all)
panels = subset(elements$layout, grepl("strip-r", elements$layout$name), t:r)
chemos = rectGrob(gp = gpar(col = "#1B9E77", lwd = 5, fill = NA))
gt <- gtable::gtable_add_grob(elements, chemos ,
                      t = 7, l = 15, b = 21, r = 17)
metals = rectGrob(gp = gpar(col = "#D95F02", lwd = 5, fill = NA))
gt <- gtable::gtable_add_grob(gt, metals,
                      t = 23, l = 15, b = 29, r = 17)
pesticides = rectGrob(gp = gpar(col = "#7570B3", lwd = 5, fill = NA))
gt <- gtable::gtable_add_grob(gt, pesticides,
                      t = 31, l = 15, b = 35, r = 17)
neuro = rectGrob(gp = gpar(col = "#E7298A", lwd = 5, fill = NA))
gt <- gtable::gtable_add_grob(gt, neuro,
                      t = 37, l = 15, b = 37, r = 17)
pdf("~/Dropbox/AndersenLab/QTLpaper/figures/final/figure1_uniqueQTL.pdf", width = 7.5, height = 7.5)
grid.newpage()
grid.draw(gt)
dev.off()

######## Figure 2: VE v. H2 and H2 v. h2 ########
load("~/Dropbox/AndersenLab/QTLpaper/data/papermappings.Rda") 
load("~/Dropbox/AndersenLab/QTLpaper/data/20170626_additive_interactive_QTL.RData")
load("~/Dropbox/AndersenLab/QTLpaper/data/uniqueQTL_588corrcut.Rda")
drugclasses <- read.csv("~/Dropbox/AndersenLab/QTLpaper/data/drugclasses.csv") %>%
  dplyr::mutate(class = sub("Heavy Metal", "Heavy metal", class))
load("~/Dropbox/AndersenLab/QTLpaper/data/HTA Dose Data/pruned_doseresponses.Rda")

textsize = 8
titlesize = 10

alltraits <- papermaps %>%
  na.omit() %>%
  dplyr::filter(trait %in% uniqueQTL$drugtrait) %>%
  dplyr::arrange(trait)

additive <- VE_SE_QTL %>% 
  dplyr::ungroup() %>%
  dplyr::mutate(trait = paste0(condition, ".", trait)) %>%
  dplyr::filter(type == "Additive", trait %in% alltraits$trait) %>%
  dplyr::arrange(trait)

#Add VE for all QTL within a trait, if applicable
addQTL <- alltraits %>%
  dplyr::mutate(condition = stringr::str_split_fixed(trait, "\\.", 2)[,1])%>%
  dplyr::distinct(condition, trait, chr, .keep_all = T) %>%
  dplyr::select(trait, var_exp) 
addQTL <- stats::aggregate(var_exp~trait,data=addQTL,FUN=sum)

allQTLtraits <- dplyr::left_join(addQTL, additive, by = 'trait') %>%
  dplyr::mutate(condition = stringr::str_split_fixed(trait, "\\.", 2)[,1]) %>%
  dplyr::mutate(trait = stringr::str_split_fixed(trait, "\\.", 2)[,2]) %>%
  dplyr::left_join(drugclasses, by = c('condition' = 'drug')) %>%
  dplyr::mutate(traitclass = ifelse(grepl("median|mean|q75|q90", trait), "Large size", ifelse(trait %in% c("n", "norm.n"), "Brood size", 
                                                                                              ifelse(grepl("cv|var", trait), "Variance", "Small size"))))

#plot narrow sense heritability on x versus the ve per all QTL for a trait on y
narrow_VE <- ggplot2::ggplot(allQTLtraits) +
  ggplot2::geom_point(ggplot2::aes(x = VE, y = var_exp, color = factor(class), shape = factor(traitclass)), alpha = .8) +
  # ggplot2::scale_colour_manual("Drug Class") + #Change these colors
  ggplot2::scale_color_brewer("Drug class", palette="Dark2") +
  ggplot2::scale_shape_manual("Trait class", values = c(16, 17, 15, 4))+
  ggplot2::geom_abline(slope = 1) +
  ggplot2::ylim(0, 0.61) +
  ggplot2::xlim(0, 0.61) +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.background = ggplot2::element_rect(colour = "black", size = 1),
                 axis.text = ggplot2::element_text(size = textsize, colour = "black", face = "bold"),
                 axis.title = ggplot2::element_text(size = titlesize, color = "black", face = "bold"), 
                 legend.text = element_text(size = (textsize - 2)),
                 legend.key.size = unit(0.4, "cm")) +
  ggplot2::labs(x = "Narrow-sense heritability", y = "Variance explained by QTL")
narrow_VE

#H2 v. h2

load("~/Dropbox/AndersenLab/QTLpaper/data/20170629_uniqueQTLH2fromdose.Rda")
H2vsh2 <- left_join(H2, additive) %>%
  dplyr::mutate(condition = stringr::str_split_fixed(trait, "\\.", 2)[,1]) %>%
  dplyr::mutate(trait = stringr::str_split_fixed(trait, "\\.", 2)[,2]) %>%
  dplyr::left_join(drugclasses, by = c('condition' = 'drug')) %>%
  dplyr::mutate(traitclass = ifelse(grepl("median|mean|q75|q90", trait), "Large Size", ifelse(trait %in% c("n", "norm.n"), "Brood Size", 
                                                                                              ifelse(grepl("cv|var", trait), "Variance", "Small Size"))))

narrow_broad <- ggplot2::ggplot(H2vsh2) +
  ggplot2::geom_point(ggplot2::aes(x = H2, y = VE, color = factor(class), shape = factor(traitclass)), alpha = .8) +
  # ggplot2::scale_colour_manual("Drug Class", values = c("purple", "red", "green", "blue")) + #Change these colors
  ggplot2::scale_color_brewer("Drug class", palette="Dark2") +
  ggplot2::scale_shape_manual("Trait class", values = c(16, 17, 15, 4))+
  ggplot2::geom_abline(slope = 1) +
  theme_bw() +
  labs(y = "Narrow-sense heritability", x = "Broad-sense heritability") + 
  ylim(0, 1) +
  xlim(0, 1) +
  #geom_errorbar(aes(x = trait, ymin = lower, ymax = higher, group = n, linetype = type), width = 0.4, position = position_dodge(0.4)) +
  #geom_errorbarh(aes(x = H2, xmin = ci_l, xmax = ci_r, y = VE)) +
  ggplot2::theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "solid", size = 1),
                 strip.text = element_text(size = titlesize, face = "bold"),
                 axis.text = element_text(size = textsize, face = "bold", color = "black"),
                 axis.title = element_text(size = titlesize, face = "bold", color = "black"),
                 panel.background = ggplot2::element_rect(colour = "black", size = 1),
                 legend.position = "none")

fig4 <- plot_grid(narrow_broad, narrow_VE, labels = c("A", "B"), nrow = 1, ncol =2, rel_widths = c(1, 1.4))
fig4
ggsave(fig4, filename = "~/Dropbox/AndersenLab/QTLpaper/figures/final/fig4_narrow_broad_hert.pdf", width = 7.5, height = 3)


######## Figure 3: Hotspots #########
library(linkagemapping)
data("N2xCB4856cross")

load("~/Dropbox/AndersenLab/RCode/GWAS/Ancillary/marker_pos_conversion.Rda")

map <- pull.map(N2xCB4856cross, as.table = TRUE) %>%
  dplyr::mutate(id = rownames(.))

mappos <- merge(map, conv, by= "id")
colnames(mappos)[3:4] <- c("cM", "pos")

#load in dataframe of unique QTL (annotated_lod df for each trait you are using for the plot, distinct will be called later in case of multiple peaks per QTL)
load("~/Dropbox/AndersenLab/QTLpaper/data/uniqueQTL_588corrcut.Rda")

QTLperbin <- function(mappos, uniqueQTL, numcm = 50, pval = .05, bonferroni = "FALSE") {
  
  #round to the nearest bin (first bin is bin "0")
  rounded <- mappos %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(bin_cM = numcm*(floor(cM/numcm)))
  
  #find cutoffs for each bin
  grouped <- rounded %>%
    dplyr::group_by(chr, bin_cM) %>%
    dplyr::summarise(start = min(pos), end = max(pos))
  
  
  #Count how many "unique" QTL we found
  allunique <- uniqueQTL %>%
    na.omit() %>%
    dplyr::group_by(chr, drugtrait) %>%
    distinct()
  
  #count number of bins in each chr (keeping in mind that bin order starts at "0" not "1")
  binsperchr <- grouped %>%
    dplyr::group_by(chr) %>%
    dplyr::summarize(totalbin = max(bin_cM)/numcm + 1)
  
  totalbins <- sum(binsperchr$totalbin)
  
  #number of QTL per bin based on poisson distribution with mean of...:
  lambda <- nrow(allunique)/totalbins
  
  pval <- ifelse(bonferroni == "TRUE", 1-(pval/totalbins), 1-pval)
  
  
  #plot by bin
  mappedunique <- uniqueQTL %>%
    na.omit()
  
  #add bin information to map
  mappedwithbin <- mappedunique %>%
    dplyr::mutate(bin = NA)
  
  
  for (i in 1:nrow(mappedwithbin)) {
    for(j in 1:nrow(grouped)) {
      mappedwithbin$bin[i] <- 
        ifelse(mappedwithbin$chr[i] == grouped$chr[j] &&
                 mappedwithbin$pos[i] <= grouped$end[j] &&
                 mappedwithbin$pos[i] >= grouped$start[j],
               grouped$bin_cM[j],
               mappedwithbin$bin[i])
    }
  }
  
  mappedwithbin <- mappedwithbin %>%
    ungroup() %>%
    dplyr::group_by(chr, trait) %>%
    distinct(.keep_all = TRUE)
  
  peaksperbin <- rounded %>%
    dplyr::mutate(peaks = 0)
  
  for(i in 1:nrow(peaksperbin)) {
    for(j in 1:nrow(mappedwithbin)) {
      peaksperbin$peaks[i] <- ifelse(peaksperbin$chr[i] == mappedwithbin$chr[j] &&
                                       peaksperbin$bin_cM[i] == mappedwithbin$bin[j],
                                     peaksperbin$peaks[i]+1, peaksperbin$peaks[i])
    }
  }
  
  peaksperbin <- peaksperbin %>%
    ungroup() %>%
    dplyr::group_by(chr) %>%
    arrange(chr,pos)
  
  plot <- ggplot2::ggplot(peaksperbin) +
    ggplot2::aes(x = pos/1e6, y = peaks) +
    ggplot2::theme_bw() +
    ggplot2::geom_line()+
    ggplot2::geom_hline(yintercept = qpois(pval, lambda), color = "red") +
    ggplot2::facet_grid(.~chr, scales = "free", space = "free") +
    ggplot2::labs(x = "Genomic position (Mb)", y = "Number of peaks in bin") +
    ggplot2::theme(axis.text.x = element_text(size=10, face="bold", color="black"),
                   axis.text.y = element_text(size=10, face="bold", color="black"),
                   axis.title.x = element_text(size=12, face="bold", color="black"),
                   axis.title.y = element_text(size=12, face="bold", color="black"),
                   strip.text.x = element_text(size=12, face="bold", color="black"),
                   strip.text.y = element_text(size=12, color="black", face="bold", angle = 0),
                   strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
                   plot.title = element_blank())
  
  summary <- list(plot, c(data.frame(c(nrow(allunique)), totalbins, lambda, bonferroni, pval)), peaksperbin)
  assign("summary", summary, envir = globalenv())
  
}

QTLperbin(mappos, uniqueQTL)

ggsave(summary[[1]], filename = "~/Dropbox/AndersenLab/QTLpaper/figures/final/fig2_50cMbins.pdf", width = 7.5, height = 2.5, units = "in")



####### Figure 4: NIL and CSS overview figure #####
load("~/Dropbox/AndersenLab/QTLpaper/data/HTA NIL CSS Data/allNILCSSregressed.Rdata")
load("~/Dropbox/AndersenLab/QTLpaper/data/HTA NIL CSS Data/20170716_allstats.Rda")
load("~/Dropbox/AndersenLab/QTLpaper/data/papermappings.Rda")
load("~/Dropbox/AndersenLab/QTLpaper/data/20170626_additive_interactive_QTL.RData")

#Values for plot
x1 <- .13
y1 <- .07
y2 <- 0.02
blck <- .078
spce <- .025
x2 <- x1+blck+spce
x3 <- x2+blck+spce
x4 <- x3+blck+spce
x5 <- x4+blck+spce*2
x6 <- x5+blck+spce
x7 <- x6+blck+spce
x8 <- x7+blck+spce
width <- .03
chr <- 0.015
textsize = 10
titlesize = 12
pointsize = 0.2

# Plot the CSS and NILS together #
plot_cssIV_IVR <- function(cond, pltrt, signif, ylab = pltrt){
  phen_gen <- allNILCSSregressed %>%
    dplyr::filter(experiment %in% c("CSSIV", "IVR"), trait == pltrt, condition == cond) %>%
    dplyr::mutate(experiment = ifelse(experiment == "CSSIV", "CSS Assay", "NIL Assay"))
  
  if(is.null(signif)) {
    test <- phen_gen %>%
      ggplot(.)+
      aes(x = factor(strain, levels = c("N2", "CB4856", "ECA575", "ECA574", "ECA241", "ECA240"), labels = c("N2", "CB4856", "ECA575", "ECA574", "ECA241", "ECA240"), ordered = T), y = phenotype,fill=strain)+
      scale_fill_manual(values = c("N2" = "orange", "CB4856"= "blue", "ECA574" = "gray", "ECA575" = "gray", "ECA240" = "gray", "ECA241" = "gray"))+
      geom_jitter(alpha = 1, size = pointsize)+
      # geom_text(data = signif, aes(x = strain, y = pheno*1.2), label = "*", size = 12, color = "black") +
      geom_boxplot(outlier.colour = NA, aes(alpha = 0.8))+
      theme_bw()+
      facet_grid(~experiment, scales = "free_x", space = "free") +
      theme(axis.text.x = element_text(size=textsize, face="bold", color="black"),
            axis.text.y = element_text(size=textsize,  face = "bold", color="black"),
            axis.title.x = element_text(size=titlesize, color="black", vjust=-.3, face = "bold"),
            axis.title.y = element_text(size=titlesize,  color="black", face = "bold"),
            strip.text = element_text(size=textsize,  color="black", face = "bold"),
            legend.position="none",
            plot.title = element_text(size=titlesize, color = "black", face = "bold", hjust = 0.5),
            panel.background = element_rect(color="black",size=1),
            strip.background = element_rect(fill = "white", color = "black", size = 1)) + 
      labs(y= ylab, x = "", title = paste(toupper(substr(cond, 1, 1)), substr(cond, 2, nchar(cond)), sep="")) 
  } else {
    signif$experiment <- factor(signif$experiment)
    signif <- signif %>%
      dplyr::mutate(experiment = ifelse(experiment == "CSSIV", "CSS Assay", "NIL Assay"))
    
    test <- phen_gen %>%
      ggplot(.)+
      aes(x = factor(strain, levels = c("N2", "CB4856", "ECA575", "ECA574", "ECA241", "ECA240"), labels = c("N2", "CB4856", "ECA575", "ECA574", "ECA241", "ECA240"), ordered = T), y = phenotype,fill=strain)+
      scale_fill_manual(values = c("N2" = "orange", "CB4856"= "blue", "ECA574" = "gray", "ECA575" = "gray", "ECA240" = "gray", "ECA241" = "gray"))+
      geom_jitter(alpha = 0.5, size = pointsize)+
      geom_text(data = signif, aes(x = strain, y = pheno*1.2), label = "*", size = 8, color = "black") +
      geom_boxplot(outlier.colour = NA, aes(alpha = 0.8))+
      theme_bw()+
      facet_grid(~experiment, scales = "free", space = "free", drop = TRUE) +
      theme(axis.text.x = element_text(size=textsize, face="bold", color="black"),
            axis.text.y = element_text(size=textsize,  face = "bold", color="black"),
            axis.title.x = element_text(size=titlesize, color="black", vjust=-.3, face = "bold"),
            axis.title.y = element_text(size=titlesize,  color="black", face = "bold"),
            strip.text = element_text(size=textsize,  color="black", face = "bold"),
            legend.position="none",
            plot.title = element_text(size=titlesize, color = "black", face = "bold", hjust = 0.5),
            panel.background = element_rect(color="black",size=1),
            strip.background = element_rect(fill = "white", color = "black", size = 1, margin = margin(0.2, 0.2, 0.2, 0.2))) + 
      ggplot2::ylim(NA, max(phen_gen$phenotype)*1.5)+
      labs(y= ylab, x = "", title = paste(toupper(substr(cond, 1, 1)), substr(cond, 2, nchar(cond)), sep=""))  
  }
  
  ggdraw() +
    draw_plot(test, 0, 0, 1, 1) +
    draw_label("Chr IV", x = 0.01, y = .10, size = textsize, hjust = -.16)+
    geom_rect(aes(xmin = x1+chr, xmax = x1+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #N2 chr
    geom_rect(aes(xmin = x2+chr, xmax = x2+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #CB chr
    geom_rect(aes(xmin = x3+chr, xmax = x3+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #First NIL chr
    geom_rect(aes(xmin = x4+chr, xmax = x4+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #Last NIL chr
    draw_label("Genome", x = 0.01, y = .04, hjust = -.1, size = textsize) +
    geom_rect(aes(xmin = x1, xmax = x1+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange")+ #N2 genome
    geom_rect(aes(xmin = x2, xmax = x2+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #CB genome
    geom_rect(aes(xmin = x3, xmax = x3+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #First NIL gehome
    geom_rect(aes(xmin = x4, xmax = x4+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange") + #Last NIL genome
    ### ADD CHRIVR ###
    geom_rect(aes(xmin = x5+chr, xmax = x5+blck-chr, ymin = y1, ymax = y1+ width),
              colour = "black", fill = "orange")+ #N2 chr
    geom_rect(aes(xmin = x6+chr, xmax = x6+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #CB chr
    geom_rect(aes(xmin = x7+chr, xmax = x7+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #First NIL chr
    geom_rect(aes(xmin = x7+blck*.6, xmax = x7+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #first NIL introgressed
    geom_rect(aes(xmin = x8+chr, xmax = x8+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #Last NIL chr
    geom_rect(aes(xmin = x8+blck*.6, xmax = x8+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #Last NIL introgressed
    geom_rect(aes(xmin = x5, xmax = x5+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange")+ #N2 genome
    geom_rect(aes(xmin = x6, xmax = x6+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #CB genome
    geom_rect(aes(xmin = x7, xmax = x7+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #First NIL genome
    geom_rect(aes(xmin = x8, xmax = x8+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange") #last NIL genome
}
plot_cssIV_IVL <- function(cond, pltrt, signif, ylab = pltrt){
  phen_gen <- allNILCSSregressed %>%
    dplyr::filter(experiment %in% c("CSSIV", "IVL"), trait == pltrt, condition == cond) %>%
    dplyr::mutate(experiment = ifelse(experiment == "CSSIV", "CSS Assay", "NIL Assay"))
  
  if(is.null(signif)) {
    test <- phen_gen %>%
      ggplot(.)+
      aes(x = factor(strain, levels = c("N2", "CB4856", "ECA575", "ECA574", "ECA229", "ECA231"), labels = c("N2", "CB4856", "ECA575", "ECA574", "ECA229", "ECA231"), ordered = T), y = phenotype,fill=strain)+
      scale_fill_manual(values = c("N2" = "orange", "CB4856"= "blue", "ECA574" = "gray", "ECA575" = "gray", "ECA229" = "gray", "ECA231" = "gray"))+
      geom_jitter(alpha = 1, size = pointsize)+
      # geom_text(data = signif, aes(x = strain, y = pheno*1.2), label = "*", size = 12, color = "black") +
      geom_boxplot(outlier.colour = NA, aes(alpha = 0.8))+
      theme_bw()+
      facet_grid(~experiment, scales = "free", space = "free") +
      # facet_wrap(~experiment, scales = "free", nrow = 1, ncol = 2) +
      theme(axis.text.x = element_text(size=textsize, face="bold", color="black"),
            axis.text.y = element_text(size=textsize,  face = "bold", color="black"),
            axis.title.x = element_text(size=titlesize, color="black", vjust=-.3, face = "bold"),
            axis.title.y = element_text(size=titlesize,  color="black", face = "bold"),
            strip.text = element_text(size=textsize,  color="black", face = "bold"),
            legend.position="none",
            plot.title = element_text(size=titlesize, color = "black", face = "bold", hjust = 0.5),
            panel.background = element_rect(color="black",size=1),
            strip.background = element_rect(fill = "white", color = "black", size = 1)) + 
      labs(y= ylab, x = "", title = paste(toupper(substr(cond, 1, 1)), substr(cond, 2, nchar(cond)), sep="")) 
  } else {
    signif$experiment <- factor(signif$experiment)
    signif <- signif %>%
      dplyr::mutate(experiment = ifelse(experiment == "CSSIV", "CSS Assay", "NIL Assay"))
    
    test <- phen_gen %>%
      ggplot(.)+
      aes(x = factor(strain, levels = c("N2", "CB4856", "ECA575", "ECA574", "ECA229", "ECA231"), labels = c("N2", "CB4856", "ECA575", "ECA574", "ECA229", "ECA231"), ordered = T), y = phenotype,fill=strain)+
      scale_fill_manual(values = c("N2" = "orange", "CB4856"= "blue", "ECA574" = "gray", "ECA575" = "gray", "ECA229" = "gray", "ECA231" = "gray"))+
      geom_jitter(alpha = 0.5, size = pointsize)+
      geom_text(data = signif, aes(x = strain, y = pheno*1.2), label = "*", size = 8, color = "black") +
      geom_boxplot(outlier.colour = NA, aes(alpha = 0.8))+
      theme_bw()+
      facet_grid(~experiment, scales = "free", space = "free") +
      # facet_wrap(~experiment, scales = "free", nrow = 1, ncol = 2) +
      theme(axis.text.x = element_text(size=textsize, face="bold", color="black"),
            axis.text.y = element_text(size=textsize,  face = "bold", color="black"),
            axis.title.x = element_text(size=titlesize, color="black", vjust=-.3, face = "bold"),
            axis.title.y = element_text(size=titlesize,  color="black", face = "bold"),
            strip.text = element_text(size=textsize,  color="black", face = "bold"),
            legend.position="none",
            plot.title = element_text(size=titlesize, color = "black", face = "bold", hjust = 0.5),
            panel.background = element_rect(color="black",size=1),
            strip.background = element_rect(fill = "white", color = "black", size = 1)) + 
      # ggplot2::ylim(NA, max(phen_gen$phenotype)*1.5)+
      labs(y= ylab, x = "", title = paste(toupper(substr(cond, 1, 1)), substr(cond, 2, nchar(cond)), sep=""))  
  }
  
  ggdraw() +
    draw_plot(test, 0, 0, 1, 1) +
    draw_label("Chr IV", x = 0.01, y = .10, size = textsize, hjust = -.16)+
    geom_rect(aes(xmin = x1+chr, xmax = x1+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #N2 chr
    geom_rect(aes(xmin = x2+chr, xmax = x2+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #CB chr
    geom_rect(aes(xmin = x3+chr, xmax = x3+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #First NIL chr
    geom_rect(aes(xmin = x4+chr, xmax = x4+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #Last NIL chr
    draw_label("Genome", x = 0.01, y = .04, hjust = -.1, size = textsize) +
    geom_rect(aes(xmin = x1, xmax = x1+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange")+ #N2 genome
    geom_rect(aes(xmin = x2, xmax = x2+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #CB genome
    geom_rect(aes(xmin = x3, xmax = x3+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #First NIL gehome
    geom_rect(aes(xmin = x4, xmax = x4+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange") + #Last NIL genome
    ### ADD CHRIVR ###
    geom_rect(aes(xmin = x5+chr, xmax = x5+blck-chr, ymin = y1, ymax = y1+ width),
              colour = "black", fill = "orange")+ #N2 chr
    geom_rect(aes(xmin = x6+chr, xmax = x6+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #CB chr
    geom_rect(aes(xmin = x7+chr, xmax = x7+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #First NIL chr
    geom_rect(aes(xmin = x7+blck*.4, xmax = x7+blck*.6, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #first NIL introgressed
    geom_rect(aes(xmin = x8+chr, xmax = x8+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #Last NIL chr
    geom_rect(aes(xmin = x8+blck*.4, xmax = x8+blck*.6, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #Last NIL introgressed
    geom_rect(aes(xmin = x5, xmax = x5+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange")+ #N2 genome
    geom_rect(aes(xmin = x6, xmax = x6+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #CB genome
    geom_rect(aes(xmin = x7, xmax = x7+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #First NIL genome
    geom_rect(aes(xmin = x8, xmax = x8+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange") #last NIL genome
}
plot_cssV_V <- function(cond, pltrt, signif, ylab = pltrt){
  phen_gen <- allNILCSSregressed %>%
    dplyr::filter(experiment %in% c("CSSV", "V"), trait == pltrt, condition == cond) %>%
    dplyr::mutate(experiment = ifelse(experiment == "CSSV", "CSS Assay", "NIL Assay"))
  
  if(is.null(signif)) {
    test <- phen_gen %>%
      ggplot(.)+
      aes(x = factor(strain, levels = c("N2", "CB4856", "ECA554", "ECA573", "ECA230", "ECA232"), labels = c("N2", "CB4856", "ECA554", "ECA573", "ECA230", "ECA232"), ordered = T), y = phenotype,fill=strain)+
      scale_fill_manual(values = c("N2" = "orange", "CB4856"= "blue", "ECA554" = "gray", "ECA573" = "gray", "ECA230" = "gray", "ECA232" = "gray"))+
      geom_jitter(alpha = 1, size = pointsize)+
      # geom_text(data = signif, aes(x = strain, y = pheno*1.2), label = "*", size = 12, color = "black") +
      geom_boxplot(outlier.colour = NA, aes(alpha = 0.8))+
      theme_bw()+
      facet_grid(~experiment, scales = "free", space = "free") +
      theme(axis.text.x = element_text(size=textsize, face="bold", color="black"),
            axis.text.y = element_text(size=textsize,  face = "bold", color="black"),
            axis.title.x = element_text(size=titlesize, color="black", vjust=-.3, face = "bold"),
            axis.title.y = element_text(size=titlesize,  color="black", face = "bold"),
            strip.text = element_text(size=textsize,  color="black", face = "bold"),
            legend.position="none",
            plot.title = element_text(size=titlesize, color = "black", face = "bold", hjust = 0.5),
            panel.background = element_rect(color="black",size=1),
            strip.background = element_rect(fill = "white", color = "black", size = 1)) + 
      labs(y= ylab, x = "", title = paste(toupper(substr(cond, 1, 1)), substr(cond, 2, nchar(cond)), sep="")) 
  } else {
    signif$experiment <- factor(signif$experiment)
    signif <- signif %>%
      dplyr::mutate(experiment = ifelse(experiment == "CSSV", "CSS Assay", "NIL Assay"))
    
    test <- phen_gen %>%
      ggplot(.)+
      aes(x = factor(strain, levels = c("N2", "CB4856", "ECA554", "ECA573", "ECA230", "ECA232"), labels = c("N2", "CB4856", "ECA554", "ECA573", "ECA230", "ECA232"), ordered = T), y = phenotype,fill=strain)+
      scale_fill_manual(values = c("N2" = "orange", "CB4856"= "blue", "ECA554" = "gray", "ECA573" = "gray", "ECA230" = "gray", "ECA232" = "gray"))+
      geom_jitter(alpha = 0.5, size = pointsize)+
      geom_text(data = signif, aes(x = strain, y = pheno*1.2), label = "*", size = 8, color = "black") +
      geom_boxplot(outlier.colour = NA, aes(alpha = 0.8))+
      theme_bw()+
      facet_grid(~experiment, scales = "free", space = "free") +
      theme(axis.text.x = element_text(size=textsize, face="bold", color="black"),
            axis.text.y = element_text(size=textsize,  face = "bold", color="black"),
            axis.title.x = element_text(size=titlesize, color="black", vjust=-.3, face = "bold"),
            axis.title.y = element_text(size=titlesize,  color="black", face = "bold"),
            strip.text = element_text(size=textsize,  color="black", face = "bold"),
            legend.position="none",
            plot.title = element_text(size=titlesize, color = "black", face = "bold", hjust = 0.5),
            panel.background = element_rect(color="black",size=1),
            strip.background = element_rect(fill = "white", color = "black", size = 1)) + 
      ggplot2::ylim(NA, max(phen_gen$phenotype)*1.5)+
      labs(y= ylab, x = "", title = paste(toupper(substr(cond, 1, 1)), substr(cond, 2, nchar(cond)), sep=""))  
  }
  
  ggdraw() +
    draw_plot(test, 0, 0, 1, 1) +
    draw_label("Chr V", x = 0.01, y = .10, size = textsize, hjust = -.16)+
    geom_rect(aes(xmin = x1+chr, xmax = x1+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #N2 chr
    geom_rect(aes(xmin = x2+chr, xmax = x2+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #CB chr
    geom_rect(aes(xmin = x3+chr, xmax = x3+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #First NIL chr
    geom_rect(aes(xmin = x4+chr, xmax = x4+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #Last NIL chr
    draw_label("Genome", x = 0.01, y = .04, hjust = -.1, size = textsize) +
    geom_rect(aes(xmin = x1, xmax = x1+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange")+ #N2 genome
    geom_rect(aes(xmin = x2, xmax = x2+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #CB genome
    geom_rect(aes(xmin = x3, xmax = x3+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #First NIL gehome
    geom_rect(aes(xmin = x4, xmax = x4+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange") + #Last NIL genome
    ### ADD CHRIVR ###
    geom_rect(aes(xmin = x5+chr, xmax = x5+blck-chr, ymin = y1, ymax = y1+ width),
              colour = "black", fill = "orange")+ #N2 chr
    geom_rect(aes(xmin = x6+chr, xmax = x6+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #CB chr
    geom_rect(aes(xmin = x7+chr, xmax = x7+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #First NIL chr
    geom_rect(aes(xmin = x7+blck*.4, xmax = x7+blck*.6, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #first NIL introgressed
    geom_rect(aes(xmin = x8+chr, xmax = x8+blck-chr, ymin = y1, ymax = y1+width),
              colour = "black", fill = "orange")+ #Last NIL chr
    geom_rect(aes(xmin = x8+blck*.4, xmax = x8+blck*.6, ymin = y1, ymax = y1+width),
              colour = "black", fill = "blue")+ #Last NIL introgressed
    geom_rect(aes(xmin = x5, xmax = x5+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange")+ #N2 genome
    geom_rect(aes(xmin = x6, xmax = x6+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #CB genome
    geom_rect(aes(xmin = x7, xmax = x7+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "blue")+ #First NIL genome
    geom_rect(aes(xmin = x8, xmax = x8+blck, ymin = y2, ymax = y2+width),
              colour = "black", fill = "orange") #last NIL genome
}

cidefiner <- function(cis, map) {
  ci_lod <- vapply(1:nrow(map), function(marker) {
    pos <- map$pos[marker]
    chr <- map$chr[marker]
    lod <- map$lod[marker]
    s <- sum(chr == cis$chr & (pos >= cis$ci_l_pos & pos <= cis$ci_r_pos))
    inci <- as.logical(s)
    cilodscore <- ifelse(inci, lod, 0)
    return(cilodscore)
  }, numeric(1))
  map$ci_lod <- ci_lod
  return(map)
}
maxlodplot_kt <- function (map, t) {
  map1 <- map %>% dplyr::group_by(marker) %>% dplyr::filter(lod ==max(lod))
  cis <- map %>% dplyr::group_by(marker) %>% dplyr::mutate(maxlod = max(lod)) %>%
    dplyr::group_by(iteration) %>% dplyr::filter(!is.na(var_exp)) %>%
    dplyr::do(head(., n = 1))
  if (nrow(cis) == 0) {
    plot <- ggplot2::ggplot(map1, ggplot2::aes(x = genotype,y = pheno)) + ggplot2::geom_blank()
    return(plot)
  }
  map1 <- cidefiner(cis, map1)
  plot <- ggplot2::ggplot(map1) + ggplot2::aes(x = pos/1e+06,y = lod)
  if (nrow(cis) != 0) {
    plot <- plot + ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e+06,ymin = 0, ymax = ci_lod), fill = "blue", alpha = 0.5) +
      ggplot2::geom_point(data = cis, ggplot2::aes(x = pos/1e+06,y = (1.05 * maxlod)), fill = "red", shape = 25,
                          size = 2, show.legend = FALSE) + 
      ggplot2::geom_text(data = cis, ggplot2::aes(x = pos/1e+06, y = (1.2 * maxlod), label = paste0(100 *round(var_exp, digits = 4), "%")), 
                         colour = "black",size = 3, hjust = "inward")
  }
  split <- cis %>%
    mutate(condition = stringr::str_split_fixed(trait, "\\.", 2)[[1]],
           trait = stringr::str_split_fixed(trait, "\\.", 2)[[2]]) 
  plot <- plot + ggplot2::geom_line(size = 0.75, alpha = 0.85) +
    ggplot2::facet_grid(. ~ chr, scales = "free", space = "free") +
    ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD", title =  paste0(toupper(substr(t, 1, 1)), substr(t, 2, nchar(t)))) +
    ggplot2::scale_colour_discrete(name = "Mapping\nIteration") +
    ggplot2::theme_bw() +
    ggplot2::ylim(0, max(cis$maxlod)*1.3)+
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
                   axis.text.y = ggplot2::element_text(size = textsize, face = "bold", color = "black"),
                   axis.title.x = ggplot2::element_text(size = titlesize, face = "bold", color="black"),
                   axis.title.y = ggplot2::element_text(size = titlesize, face = "bold",color="black"),
                   strip.text.x = ggplot2::element_text(size = titlesize, face = "bold", color = "black", margin = margin(2, 2, 2, 2)),
                   strip.text.y = ggplot2::element_text(size = titlesize, face = "bold", color = "black"),
                   strip.background = element_rect(colour = "black", fill = "white", size = 0.75, linetype = "solid"),
                   plot.title = element_text(color = "black", face = "bold", size = titlesize),
                   # plot.margin = unit(c(0, 0, 0, 0), "cm"),
                   panel.background = ggplot2::element_rect(color = "black", size = 1.2))
  return(plot)
}
plot_add_int <- function(df) {
  ggplot2::ggplot(df) +
    ggplot2::geom_bar(data = df, ggplot2::aes(x = trait, 
                                              y = VE, fill = factor(type, levels = c("Interactive", "Additive"))), stat = "identity") +
    ggplot2::scale_fill_manual("Key", values = c("slateblue2", "cadetblue3"))+
    ggplot2::labs(y = "VE", x = "") + 
    ggplot2::geom_errorbar(ggplot2::aes(x = trait, ymin = lower, ymax = higher, group = n, linetype = type), width = 0.7, position = position_dodge(0.4)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = textsize, face = "bold"),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_text(size = titlesize, face = "bold"),
                   axis.title.y = ggplot2::element_text(size = titlesize, face = "bold"),
                   legend.position = "none")
}

#Stats calcs for each hotspot
getNIL <- function(sp, exp, cond, trt) {
  if(exp == "IVL" && sp == "N2_nil_sig") { nil <- "ECA229" } else
    if(exp == "IVR" && sp == "N2_nil_sig") { nil <- "ECA241" } else
      if(exp == "V" && sp == "N2_nil_sig") { nil <- "ECA230" } else
        if(exp == "CSSIV" && sp == "N2_nil_sig") { nil <- "ECA575" } else
          if(exp == "CSSV" && sp == "N2_nil_sig") { nil <- "ECA554" } else
            if(exp == "IVL" && sp == "CB_nil_sig") { nil <- "ECA231" } else
              if(exp == "IVR" && sp == "CB_nil_sig") { nil <- "ECA240" } else
                if(exp == "V" && sp == "CB_nil_sig") { nil <- "ECA232" } else
                  if(exp == "CSSIV" && sp == "CB_nil_sig") { nil <- "ECA574" } else
                    if(exp == "CSSV" && sp == "CB_nil_sig") { nil <- "ECA573" }
  regressed <- allNILCSSregressed %>%
    dplyr::filter(condition == cond, trait == trt, experiment == exp, strain == nil)
  pheno <- max(regressed$phenotype)
  ans <- c(nil, pheno)
  return(ans)
}
plotstats <- function(drug, trt, expr) {
  fullstats <- NULL
  for(j in 1:length(expr)) {
    e <- expr[j]
    stats <- allstats %>%
      dplyr::filter(exp == e, condition == drug, trait == trt) %>%
      dplyr::select(condition, trait, N2_nil_sig, CB_nil_sig, experiment = exp) %>%
      tidyr::gather(strain_pair, significant, N2_nil_sig:CB_nil_sig) %>%
      na.omit() %>%
      dplyr::mutate(strain = c(NA), pheno = c(NA)) 
    if(nrow(stats) != 0) {
      for(i in 1:nrow(stats)) {
        test <- getNIL(stats$strain_pair[i], e, drug, trt)
        stats$strain[i] <- test[1]
        stats$pheno[i] <- as.numeric(test[2])
      } 
      fullstats <- rbind(fullstats, stats)
    } else {
      stats <- NULL
      fullstats <- rbind(fullstats, stats)}
  }
  return(fullstats)
}

#perfect category
drug <- "deiquat"
trt <- "q75.EXT"

VE <- plot_grid(NULL, plot_add_int(VE_SE_QTL %>% dplyr::filter(condition == drug, trait == trt)), nrow = 2, rel_heights = c(.2,1.5))
nil_css1 <- plot_cssIV_IVR(drug, trt, plotstats(drug, trt, c("IVR", "CSSIV")))
a <- plot_grid(nil_css1, VE, align = "h", nrow = 1, rel_widths = c(1,0.2)) + 
  draw_plot_label("(i)", x = .03, y = 0.97, size = 14, fontface = "plain") +
  draw_plot_label("(ii)", x = .83, y = 0.97, size = 14, fontface = "plain")

##### interchromosomal external bidirectional
drug <- "silver"
trt <- "median.TOF"

VE <- plot_grid(NULL, plot_add_int(VE_SE_QTL %>% dplyr::filter(condition == drug, trait == trt)), nrow = 2, rel_heights = c(.2,1.5))
nil_css2 <- plot_cssV_V(drug, trt, plotstats(drug, trt, c("V", "CSSV")))
b <- plot_grid(nil_css2, VE, align = "h", nrow = 1, rel_widths = c(1,0.2))+ 
  draw_plot_label("(i)", x = .03, y = 0.97, size = 14, fontface = "plain") +
  draw_plot_label("(ii)", x = .83, y = 0.97, size = 14, fontface = "plain")

##### interchromosomal internal unidirectional
drug <- "silver"
trt <- "q10.EXT"

VE <- plot_grid(NULL, plot_add_int(VE_SE_QTL %>% dplyr::filter(condition == drug, trait == trt)), nrow = 2, rel_heights = c(.2,1.5))
nil_css3 <- plot_cssIV_IVL(drug, trt, plotstats(drug, trt, c("IVL", "CSSIV")))
c <- plot_grid(nil_css3, VE, align = "h", nrow = 1, rel_widths = c(1,0.2))+ 
  draw_plot_label("(i)", x = .03, y = 0.97, size = 14, fontface = "plain") +
  draw_plot_label("(ii)", x = .83, y = 0.97, size = 14, fontface = "plain")

##### interchromosomal external unidirectional
drug <- "mechlorethamine"
trt <- "norm.n"

VE <- plot_grid(NULL, plot_add_int(VE_SE_QTL %>% dplyr::filter(condition == drug, trait == trt)), nrow = 2, rel_heights = c(.2,1.5))
nil_css4 <- plot_cssIV_IVR(drug, trt, plotstats(drug, trt, c("IVR", "CSSIV")))
d <- plot_grid(nil_css4, VE, align = "h", nrow = 1, rel_widths = c(1,0.2))+ 
  draw_plot_label("(i)", x = .03, y = 0.97, size = 14, fontface = "plain") +
  draw_plot_label("(ii)", x = .83, y = 0.97, size = 14, fontface = "plain")


##### intrachromosomal unidirectional
drug <- "tunicamycin"
trt <- "q90.TOF"
# Could also do cisplatin.q90.EXT, looks very extreme! (and would be on V)

VE <- plot_grid(NULL, plot_add_int(VE_SE_QTL %>% dplyr::filter(condition == drug, trait == trt)), nrow = 2, rel_heights = c(.2,1.5))
nil_css5 <- plot_cssIV_IVL(drug, trt, plotstats(drug, trt, c("IVL", "CSSIV")))
e <- plot_grid(nil_css5, VE, align = "h", nrow = 1, rel_widths = c(1,0.2))+ 
  draw_plot_label("(i)", x = .03, y = 0.97, size = 14, fontface = "plain") +
  draw_plot_label("(ii)", x = .83, y = 0.97, size = 14, fontface = "plain")

all <- plot_grid(a, b, c, d, e, nrow = 5, ncol = 1, labels = c("A", "B", "C", "D", "E"))

ggsave(all, file = "~/Dropbox/AndersenLab/QTLpaper/figures/final/fig3_cssnil.pdf", height = 10, width = 7.5, units = "in")



######## Supplement figure for all interactions/additive stacked bar chart ########

load("~/Dropbox/AndersenLab/QTLpaper/data/20170626_additive_interactive_QTL.RData")
load("~/Dropbox/AndersenLab/QTLpaper/data/uniqueQTL_588corrcut.Rda")

textsize <- 8
titlesize <- 10

VE_SE_QTL <- VE_SE_QTL %>%
  dplyr::mutate(condtrt = paste0(condition, ".", trait)) %>%
  dplyr::filter(condtrt %in% unique(uniqueQTL$drugtrait)) %>%
  # group_by(condition, trait) %>%
  dplyr::arrange(condition, trait)

allinteractions <- ggplot(VE_SE_QTL) +
  geom_bar(data = VE_SE_QTL, aes(x = trait, 
                                 y = VE, fill = factor(type, levels = c("Interactive", "Additive"))), stat = "identity") +
  # scale_fill_manual("Key", values = c("slateblue4", "cadetblue3"))+
  scale_fill_manual("Key", values = c("slateblue2", "cadetblue3"))+
  coord_flip() +
  facet_grid(~condition, scales = "free_y", switch = "both") +
  labs(y = "Variance Explained", x = "") + 
  geom_errorbar(aes(x = trait, ymin = lower, ymax = higher, group = n, linetype = type), width = 0.7, position = position_dodge(0.4)) +
  theme(legend.position = "none",
        strip.background = element_rect(colour = "black", fill = "white", linetype = "solid", size = 1),
        strip.text = element_text(size = titlesize, face = "bold"),
        axis.text.x = element_text(size = textsize),
        axis.text.y = element_text(size = textsize),
        axis.title.x = element_blank()) +
  facet_wrap(~condition, scales = "free", ncol = 3, nrow = 6)
allinteractions 

legstartx <- 0.5 #position where the top left corner of the legend begins
legstarty <- 0.1 #position where the top left corner of the legend begins
width <- 0.02 #width of the boxes/lines of legend
sp <- 0.01 #spacing within legend

#Make my own legend because I don't know how to combine line and color like this automatically
# plot <- ggdraw() +
#   draw_plot(allinteractions, 0, 0, 1, 1) +
#   geom_rect(aes(xmin = legstartx - sp, xmax = legstartx + width*12+sp*2, ymin = legstarty-width-sp*2, ymax = legstarty + sp), colour = "black", fill = "white") +
#   geom_rect(aes(xmin = legstartx, xmax = legstartx + width, ymax = legstarty, ymin = legstarty - width*0.8), colour = "black", fill = "cadetblue3") +
#   draw_label("Additive", x = legstartx + width*3 + sp, y = legstarty - width*0.8, size = titlesize, fontface = "bold") +
#   geom_segment(aes(x = legstartx, xend = legstartx + width, y = legstarty - width*1.2, yend = legstarty - width*1.2), colour = "black", linetype = "solid")+ 
#   geom_rect(aes(xmin = legstartx+width*6, xmax = legstartx + width*7, ymax = legstarty, ymin = legstarty - width*0.8), colour = "black", fill = "slateblue2") +
#   draw_label("Interactive", x = legstartx + width*9.5 + sp, y = legstarty - width*0.8, size = titlesize, fontface = "bold") +
#   geom_segment(aes(x = legstartx+width*6, xend = legstartx + width*7, y = legstarty - width*1.2, yend = legstarty - width*1.2), colour = "black", linetype = "dashed") 
# 
# ggsave(plot, filename = "~/Dropbox/AndersenLab/QTLpaper/figures/final/supplement_allinteractions.pdf", width = 7.5, height = 10.5)

plot <- ggdraw() +
  draw_plot(allinteractions, 0, 0, 1, 1) +
  geom_rect(aes(xmin = legstartx - sp, xmax = legstartx + width*12+sp*2, ymin = legstarty-width-sp, ymax = legstarty + sp), colour = "black", fill = "white") +
  geom_rect(aes(xmin = legstartx, xmax = legstartx + width, ymax = legstarty, ymin = legstarty - width*0.8), colour = "black", fill = "cadetblue3") +
  draw_label("Additive", x = legstartx + width*3 + sp, y = legstarty - width*0.35, size = titlesize, fontface = "bold") +
  geom_segment(aes(x = legstartx - width*.25, xend = legstartx + width*1.35, y = legstarty - width*.35, yend = legstarty - width*.35), colour = "black", linetype = "solid")+ 
  geom_rect(aes(xmin = legstartx+width*6, xmax = legstartx + width*7, ymax = legstarty, ymin = legstarty - width*0.8), colour = "black", fill = "slateblue2") +
  draw_label("Interactive", x = legstartx + width*9.5 + sp, y = legstarty - width*0.35, size = titlesize, fontface = "bold") +
  geom_segment(aes(x = legstartx+width*5.75, xend = legstartx + width*7.5, y = legstarty - width*.35, yend = legstarty - width*.35), colour = "black", linetype = "dashed") 

ggsave(plot, filename = "~/Dropbox/AndersenLab/QTLpaper/figures/final/supplement_allinteractions.pdf", width = 7.5, height = 10.5)


########## Flowchart R ############
library(DiagrammeR)

png(file = "~/Dropbox/AndersenLab/QTLpaper/figures/final/flowchartR.png")
grViz("
      digraph flowchart {
      
      # a 'graph' statement
      graph [overlap = false, fontsize = 10, fontname = Helvetica]
      
      # several 'node' statements
      node [shape = box, fontname = Helvetica, fillcolor = CornflowerBlue, style = filled, fontcolor = White]
      a [label = 'All traits']
      b [label = 'Parents are significantly \\n different in CSS and NIL']
      c [label = 'CSS and NIL work']
      d [label = 'CSS and NIL \\n do not both work']
      e [label = 'CSS does not work']
      f [label = 'CSS works']
      g [label = 'NIL does not work']
      h [label = 'NIL works']
      i [label = 'Parents are not significantly \\n different for CSS and NIL']
      j [label = 'CSS parents \\n are significantly different']
      k [label = 'CSS works']
      l [label = 'CSS does not work']
      m [label = 'CSS parents are not \\n significantly different']
      n [label = 'NIL parents are \\n significantly different']
      o [label = 'NILs work']
      p [label = 'NILs do not work']
      q [label = 'NIL parents are not \\n significantly different']
      r [label = 'NILs and CSS are \\n significantly different \\n and in same direction']
      s [label = 'NILs and CSS are \\n not significantly different \\n and/or not in the \\n same direction']
      
      node [shape = box, fontname = Helvetica, fillcolor = DarkBlue, style = filled, fontcolor = White] 
      1 [label = 'Additive']
      2 [label = 'Power/Complexity']
      3 [label = 'Interchromosomal \\n Epistasis']
      4 [label = 'Intrachromosomal \\n Epistasis']
      5 [label = 'Transgressive']
      6 [label = 'Power/Complexity']
      7 [label = 'Power/Complexity']
      8 [label = 'Power/Complexity']
      # 9 [label = 'Interchromosomal \\n Epistasis']
      # 10 [label = 'Intrachromosomal \\n Epistasis']
      
      
      # several 'edge' statements
      a->b b->c b->d c->1 d->e d->f f->4 e->g e->h g->2 h->3 
      a->i i->j i->m j->k j->l k->4 l->6 m->n m->q n->o n->p o->3 p->7 q->r q->s r->5 s->8
      }
      
      
      ")
dev.off()



