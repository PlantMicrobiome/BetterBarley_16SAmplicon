library(phyloseq); packageVersion("phyloseq")
library("DESeq2");packageVersion("DESeq2") 
library ("ape")
library(ggplot2); packageVersion("ggplot2")
library (tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr") 
library ("vegan"); packageVersion("vegan") 
library("MicrobiotaProcess")
library("microbiome")

colorblind_Palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
myshapes=c(0,1,2,3,5,6,8,15,16,17,18,21,22,23,24,25)

#otu table
otu <- read.csv(file = "feature-table.tsv",row.names = 1, sep="\t",header = TRUE)
head(otu)

#taxonomy table
taxonomy <- read.table(file = "taxonomy.tsv", sep = "\t", header = T ,row.names = 1)

#clean taxonomy table
tax <- taxonomy %>%
  select(Taxon) %>% 
  separate(Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "; ")
colnames(tax)
head (taxonomy)
tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "d__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

#metadata
mapfile<-import_qiime_sample_data("metadata.tsv")

# sequences
seq=Biostrings::readDNAStringSet(filepath = "sequences.fasta")

#Physeq obj
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(as.matrix(tax.clean))
META <- sample_data(mapfile)
sequneces=refseq(seq)

#Combine OTU table, taxonomy table and metadata into a phyloseq object. 
physeq <- phyloseq(OTU, TAX, META, sequneces)
physeq
rank_names(physeq)

# Pre-processing: remove Chloroplast and Mitochondria
no_chlor <-subset_taxa(physeq, (Order!="Chloroplast") | is.na(Order))
no_chlor
no_mito <-subset_taxa(no_chlor, (Family!="Mitochondria") | is.na(Family))
no_mito

# Remove ASVs assigned to NA at Genus level
clean_1 <- subset_taxa(no_mito, Genus!= "NA")

#########################
# Remove contamination
#########################
blanks=subset_samples(clean_1, sample_or_control == "control")

# remove ASVs with zero reads
blanks_with_reads=prune_taxa(taxa_sums(blanks)>0,blanks)

# remove ASVs from all that are present in blanks
blanks_taxa=taxa_names(blanks_with_reads)
non_contam_taxa=prune_taxa(blanks_taxa, clean_1)
contam_taxa=blanks_taxa
all_taxa=taxa_names(clean_1)
keepTaxa=all_taxa[!(all_taxa %in% contam_taxa)]
clean_2=prune_taxa(keepTaxa,clean_1)

# remove all zero read/ASVs samples,i.e. controls/blanks (10)
clean_3=prune_samples(sample_sums(clean_2)>0, clean_2)

# rarefy
rare_final_clean_ps=rarefy_even_depth(final_clean_ps,rngseed = TRUE, replace=FALSE, trimOTUs = TRUE, sample.size=11000)

# Aggregate samples at genus level 
ps_genus <- tax_glom(rare_final_clean_ps, taxrank= "Genus", NArm=FALSE, bad_empty=c(NA, "", " ", "\t"))
rare_ps_genus=ps_genus

#alpha diversity
#Index calculations
alpha_rare <-  estimate_richness(rare_ps_genus, measures = c("Observed", "Shannon", "Chao1"))
alpha_rare

#design file
design <- rare_ps_genus@sam_data
design

#treatment
design_treatment <- as.data.frame(design[, 4])
rownames(design_treatment) <- rownames(design)
colnames(design_treatment) <- c("fertilizer")
design_treatment

#old_or_new
design_ovn <- as.data.frame(design[, 6])
rownames(design_ovn) <- rownames(design)
colnames(design_ovn) <- c("old_or_new")
design_ovn

#description
design_description <- as.data.frame(design[, 3])
rownames(design_description) <- rownames(design)
colnames(design_description) <- c("type")
design_description 

#data frame Genotype_Description
design_TD <- cbind(design_treatment, design_description, design_ovn)

#Observed Genera
alpha_rare_Observed <- as.data.frame(alpha_rare[ ,1])
rownames(alpha_rare_Observed) <- rownames(alpha_rare)
colnames(alpha_rare_Observed) <- c("Observed")

#Combine the dataset sample description and Observed Genera
alpha_rare_Observed_TD <- cbind(design_TD, alpha_rare_Observed)
alpha_rare_Observed_TD <- as.data.frame(alpha_rare_Observed_TD)
alpha_rare_Observed_TD$fertilizer
alpha_rare_Observed_TD$type

# Order the levels according to a defined order
alpha_rare_Observed_TD$fertilizer <- ordered(alpha_rare_Observed_TD$fertilizer, levels=c("no_fertilizer", "organic_2", "organic_4", "organic_6","mineral")) 
alpha_rare_Observed_TD$type <- ordered(alpha_rare_Observed_TD$type, levels=c("babuska","feedway", "flair", "irina", "langeland", "rgt","salka" )) 

description_names <- c(
  `babuska` = "Babushka",
  `langeland` = "Langeland",
  `salka` = "Salka",
  `feedway` = "Feedway",
  `flair` = "Flair",
  `irina` = "Irina",
  `rgt` = "RGT")

fert_names <- c(
  `no_fertilizer` = "No fertilizer",
  `mineral` = "Mineral fertilizer",
  `organic_2` = "Organic 2",
  `organic_4` = "Organic 4",
  `organic_6` = "Organic 6")

old_or_new_names<-c('old'="Old",
                    'new'="Modern")

ggplot(alpha_rare_Observed_TD, aes(x=old_or_new, y=Observed, fill=old_or_new))+
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 4.5)+
  geom_point(size=3, shape = 16)+
  scale_fill_manual(name = "Cultivar", labels = old_or_new_names, values = colorblind_Palette)+
  facet_wrap(~ fertilizer, labeller = as_labeller(fert_names))+
  theme_bw()+
  ylab("Observed Genera")+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_blank(),
        strip.text = element_text(size=18),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_blank(),
        legend.position="bottom",
        legend.title = element_text(size = 16))


# CHAO1
alpha_rare_Chao1 <- as.data.frame(alpha_rare[ ,2])
rownames(alpha_rare_Chao1) <- rownames(alpha_rare)
colnames(alpha_rare_Chao1) <- c("Chao1")

#Fig 5A
alpha_rare_Chao1_TD$fertilizer <- ordered(alpha_rare_Chao1_TD$fertilizer, levels=c("no_fertilizer", "organic_2", "organic_4", "organic_6","mineral")) 

(Chao1_box_treatment <- ggplot(data = alpha_rare_Chao1_TD, aes(x=fertilizer,y=Chao1))+
    geom_boxplot(aes(fill=fertilizer),outlier.shape = NA, lwd=1)+
    scale_fill_manual(values = c("no_fertilizer" = "#E69F00", "organic_2" = "#56B4E9", "organic_4" = "#009E73", "organic_6" = "#D55E00", "mineral" = "#999999"))+
    theme_classic()+
    theme(axis.title=element_text(size=24))+
    theme(axis.title.x = element_blank())+
    theme(axis.text.x = element_blank())+
    theme(axis.text.y = element_text(size = 23, color = "black", margin=margin(r=10)))+
    theme(legend.text = element_text(size = 10))+
    theme(legend.title = element_text(size = 10))+
    theme(axis.line = element_line(linewidth=1))+
    theme(axis.ticks = element_line(linewidth=1))+
    theme(legend.position="none")+
    guides(fill = guide_legend(nrow = 1))+
    guides(color = guide_legend(nrow = 1)))

#Combine the dataset sample description and Chao1 Genera
alpha_rare_Chao1_TD <- cbind(design_TD, alpha_rare_Chao1)
alpha_rare_Chao1_TD <- as.data.frame(alpha_rare_Chao1_TD)

#Order the levels according to a defined order
alpha_rare_Chao1_TD$fertilizer <- ordered(alpha_rare_Chao1_TD$fertilizer, levels=c("no_fertilizer", "mineral", "organic_2", "organic_4", "organic_6")) 
alpha_rare_Chao1_TD$type <- ordered(alpha_rare_Chao1_TD$type, levels=c("babuska","feedway", "flair", "irina", "langeland", "rgt","salka" )) 

ggplot(alpha_rare_Chao1_TD, aes(x=old_or_new, y=Chao1, fill=old_or_new))+
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 4.5)+
  #geom_jitter(size = 3, shape = 16)+
  geom_point(size=3, shape = 16)+
  scale_fill_manual(name = "Cultivar", labels = old_or_new_names, values = colorblind_Palette)+
  facet_wrap(~ fertilizer, labeller = as_labeller(fert_names))+
  theme_bw()+
  ylab("Chao1")+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_blank(),
        strip.text = element_text(size=18),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.position="bottom",
        legend.title = element_text(size = 16))

#Shannon Genera
alpha_rare_Shannon <- as.data.frame(alpha_rare[ ,4])
rownames(alpha_rare_Shannon) <- rownames(alpha_rare)
colnames(alpha_rare_Shannon) <- c("Shannon")

#Combine the dataset sample description and Shannon Genera
alpha_rare_Shannon_TD <- cbind(design_TD, alpha_rare_Shannon)
alpha_rare_Shannon_TD <- as.data.frame(alpha_rare_Shannon_TD)
alpha_rare_Shannon_TD$fertilizer
alpha_rare_Shannon_TD$type
alpha_rare_Shannon_TD$fertilizer <- ordered(alpha_rare_Shannon_TD$fertilizer, levels=c("no_fertilizer", "organic_2", "organic_4", "organic_6","mineral")) 

# Fig 5A
Shannon_box_treatment <- ggplot(data = alpha_rare_Shannon_TD, aes(x=fertilizer,y=Shannon))+
  geom_boxplot(aes(fill=fertilizer),outlier.shape = NA, lwd=1)+
  scale_fill_manual(values = c("no_fertilizer" = "#E69F00", "organic_2" = "#56B4E9", "organic_4" = "#009E73", "organic_6" = "#D55E00", "mineral" = "#999999"))+
  theme_classic()+
  theme(axis.title=element_text(size=23,))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size = 23, color = "black", margin=margin(r=10)))+
  theme(legend.text = element_text(size = 10))+
  theme(legend.title = element_text(size = 10))+
  theme(axis.line = element_line(linewidth=1))+
  theme(axis.ticks = element_line(linewidth=1))+
  theme(legend.position="none")+
  guides(fill = guide_legend(nrow = 1))+
  guides(color = guide_legend(nrow = 1))

Shannon_box_treatment

alpha_rare_Shannon_TD$fertilizer <- ordered(alpha_rare_Shannon_TD$fertilizer, levels=c("no_fertilizer", "mineral", "organic_2", "organic_4", "organic_6")) 
alpha_rare_Shannon_TD$type <- ordered(alpha_rare_Shannon_TD$type, levels=c("babuska","feedway", "flair", "irina", "langeland", "rgt","salka" )) 

ggplot(alpha_rare_Shannon_TD, aes(x=old_or_new, y=Shannon, fill=old_or_new))+
  geom_boxplot(outlier.colour = "red", outlier.shape = 1, outlier.size = 4.5)+
  #geom_jitter(size = 3, shape = 16)+
  geom_point(size=3, shape = 16)+geom_point(size=3, shape = 16)+
  scale_color_manual(name="Cultivar", labels=description_names, values=colorblind_Palette)+
  scale_fill_manual(name = "Age", labels = old_or_new_names, values = colorblind_Palette)+
  facet_wrap(~ fertilizer, labeller = as_labeller(fert_names))+
  theme_bw()+
  ylab("Shannon")+
  theme(axis.title.y = element_text(color="Black", size=16),
        axis.title.x = element_blank(),
        strip.text = element_text(size=18),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 14),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_blank(),
        legend.position="left",
        legend.title = element_text(size = 16))

# kruskal wallis test
data_alphadiv=cbind(design_TD,alpha_rare_Observed, alpha_rare_Chao1, alpha_rare_Shannon)
data_alphadiv_df=as.data.frame(data_alphadiv)
kruskal.test(Shannon ~ fertilizer , data = data_alphadiv)
kruskal.test(Shannon ~ old_or_new, data = data_alphadiv)
kruskal.test(Shannon ~ type, data = data_alphadiv)

library ("FSA")
library ("rcompanion")
# Dunn test
data_alphadiv_df$fertilizer=as.factor(data_alphadiv_df$fertilizer)
PT=dunnTest(Shannon ~data_alphadiv_df$fertilizer, data=data_alphadiv_df, method="bh")#"bh" the FDR is controlled using the Benjamini-Hochberg adjustment (1995), a step-down procedure appropriate to independent tests or tests that are positively dependent.
PT2=PT$res
cldList(comparison = PT2$Comparison, p.value = PT2$P.adj, threshold  = 0.05)

# beta diversity
# Fig 5b
# CAP (Canonical analysis of principal coordinates)
colorblind_Palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#D55E00", "#CC79A7", "#000000")
non_agglo_CAP=ordinate(beta_div_non_agglo, "CAP", "bray", ~ old_or_new*fertilizer)
plot_ordination(beta_div_non_agglo, non_agglo_CAP, color = "fertilizer", shape="old_or_new")

#assign color to fertilizer and shape to type
p=plot_ordination(beta_div_non_agglo, non_agglo_CAP, shape ="old_or_new", color = "fertilizer")
p = p + geom_point(size = 7, alpha = 5)
plot_CAP=p + scale_colour_manual(name = "Fertilizer", labels = fert_names, values = colorblind_Palette)+
  scale_shape_discrete(name = "Age", labels = old_or_new_names)+
  theme_bw()+
  theme(axis.title.y = element_text(color="Black", size=23),
        axis.title.x = element_text(color="Black", size=23),
        axis.text.y = element_text(size = 23),
        axis.text.x = element_text(size = 23),
        legend.position="none",
        axis.line = element_line(linewidth=1),
          axis.ticks = element_line(linewidth=1)
        )
plot_CAP
library(cowplot)
# Extract the legend
legend_Dw_Cul_alt <- get_legend(plot_CAP)

# adonis
BC <- phyloseq::distance(beta_div_non_agglo, "bray")
#Microhabitat effect
design <- read.delim("design_2.tsv", sep = "\t", header=TRUE, row.names=1)
design <- design[sample_names(beta_div_non_agglo), ]
design
adonis2(BC ~ old_or_new*fertilizer*type , data= design, permutations = 10000)

# Fig 6A
library(UpSetR)
sam_d=import_qiime_sample_data("../input/rare_final_clean_ps_samData.tsv")
meta_dm <- sample_data(sam_d)
new_rare_final_clean_ps=phyloseq(rare_final_clean_ps@otu_table,rare_final_clean_ps@tax_table,meta_dm)
upsetda <- get_upset(obj=new_rare_final_clean_ps, factorNames="fertilizer")
upset_plot=upset(upsetda, sets=unique(as.vector(sample_data(new_rare_final_clean_ps)$fertilizer)), 
                 sets.bar.color = "#56B4E9",
                 order.by = "freq", 
                 empty.intersections = "on",
                 nintersects=NA, point.size = 4, line.size = 1, mb.ratio = c(0.7, 0.3), text.scale=c(2,2,2,2,2,2), 
                 number.angles = 25,mainbar.y.label = "Intersection Size", sets.x.label = "Set size")
pdf("upsetPlot2.pdf",width=15,height=5, onefile=FALSE)
upset_plot
dev.off()

# Fig 6B
library(microViz)
supp.labs <- c("Mineral fertilizer", "No fertilizer", "Organic 2", "Organic 4", "Organic 6")
names(supp.labs) <- c("mineral", "no_fertilizer", "organic_2", "organic_4", "organic_6")
mv_plot=rare_final_clean_ps%>%
  tax_fix(unknowns = c("uncultured"))%>% 
  comp_barplot(tax_level = "Phylum", n_taxa = 10, merge_other = TRUE, label=NULL,#bar_width = 0.9,bar_outline_colour=NA,
               sample_order = c("B_min_C4","B_min_C5","B_min_C6","B_noF_C1","B_noF_C2","B_noF_C3","B_o2_C10","B_o2_C11","B_o2_C12","B_o4_C16","B_o4_C17","B_o4_C18","B_o6_C22","B_o6_C23","B_o6_C24","L_min_C52","L_min_C53","L_min_C54","L_noF_C49","L_noF_C50","L_noF_C51","L_o2_C58","L_o2_C59","L_o2_C60","L_o4_C64","L_o4_C65","L_o4_C66","L_o6_C70","L_o6_C71","L_o6_C72","S_min_C28","S_min_C29","S_min_C30","S_noF_C27","S_o2_C34","S_o2_C35","S_o2_C36","S_o4_C40","S_o4_C41","S_o4_C42","S_o6_C46","S_o6_C47","S_o6_C48","F_min_C76","F_min_C77","F_min_C78","F_noF_C73","F_noF_C74","F_noF_C75","F_o2_C82","F_o2_C83","F_o2_C84","F_o4_C88","F_o4_C89","F_o4_C90","F_o6_C94","F_o6_C95","F_o6_C96","FW_min_C124","FW_min_C125","FW_min_C126","FW_noF_C121","FW_noF_C122","FW_noF_C123","FW_o2_C130","FW_o2_C131","FW_o2_C132","FW_o4_C136","FW_o4_C137","FW_o4_C138","FW_o6_C142","FW_o6_C143","FW_o6_C144","I_min_C101","I_min_C102","I_noF_C97","I_noF_C98","I_o2_C106","I_o2_C107","I_o2_C108","I_o4_C112","I_o4_C113","I_o6_C118","I_o6_C119","I_o6_C120","R_min_C198","R_noF_C193","R_noF_C194","R_noF_C195","R_o2_C202","R_o2_C203B","R_o2_C204","R_o4_C208","R_o4_C210","R_o6_C214","R_o6_C215","R_o6_C216"),
  )+
  facet_grid(
    rows = vars(fertilizer),labeller = labeller(fertilizer = supp.labs),
    scales = "free", space = "free") +# these options are critically important!
  theme(axis.ticks.y = element_blank(), strip.text = element_text(size=12), legend.text=element_text(size=12), 
        axis.text=element_text(family="ArialMT", size=12), axis.title=element_text(family="ArialMT", size=12), 
        legend.title =element_text(family="ArialMT", size=12),legend.position = "left")
mv_plot
# relative abundance of phyla per treatment
phylum_table <- tax_glom(rare_final_clean_ps, taxrank = "Phylum")
ps.rel = transform_sample_counts(rare_final_clean_ps, function(x) x/sum(x)*100)
glom <- tax_glom(ps.rel, taxrank = "Phylum")
ps.melt <- psmelt(glom)
ps.melt$Phylum <- as.character(ps.melt$Phylum)

ps.melt <- ps.melt %>%
  group_by(fertilizer, Phylum) %>%
  mutate(median=median(Abundance))
keep <- unique(ps.melt$Phylum[ps.melt$median > 1])
ps.melt$Phylum[!(ps.melt$Phylum %in% keep)] <- "< 1%"
ps.melt_sum <- ps.melt %>%
  group_by(fertilizer,Phylum) %>%
  summarise(Abundance=sum(Abundance))

# merging bars
rare_final_clean_ps %>%
  tax_fix(unknowns = c("uncultured"))%>% 
  ps_select(fertilizer, fertilizer_cultivar) %>% # avoids lots of phyloseq::merge_samples warnings
  phyloseq::merge_samples(group = "fertilizer_cultivar", fun=mean) %>%
  comp_barplot(tax_level = "Phylum", n_taxa = 10, bar_width = 0.8, 
               sample_order=c("b_m","l_m","s_m","r_m","fw_m","f_m","i_m","b_noF","l_noF","s_noF","r_noF","fw_nF","f_noF","i_noF","b_o2","l_o2","s_o2","r_o2","fw_o2","f_o2","i_o2","b_o4","l_o4","s_o4","r_o4","fw_o4","f_o4","i_o4", "b_o6","l_o6","s_o6","r_o6","fw_o6","f_o6","i_o6")
  ) +
  coord_flip() + theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

# Differential expression
# Fig 7
head (rare_ps_genus@sam_data)
ds=phyloseq_to_deseq2(rare_ps_genus, ~fertilizer)
ds=DESeq(ds)
# no fertilizer vs. mineral
noFvsM=results(ds,contrast=c("fertilizer","no_fertilizer","mineral"), alpha=0.05)
noFvsM = noFvsM[order(noFvsM$padj, na.last=NA), ]
noFvsM_sig = noFvsM[(noFvsM$padj < 0.05), ]
noFvsM_sig
noFvsM_sig_full = cbind(as(noFvsM_sig, "data.frame"), as(phyloseq::tax_table(rare_ps_genus)[rownames(noFvsM_sig), ], "matrix"))
nrow(noFvsM_sig)

# no fertilizer vs. organic_2
noFvsO2=results(ds,contrast=c("fertilizer","no_fertilizer","organic_2"), alpha=0.05)
noFvsO2 = noFvsO2[order(noFvsO2$padj, na.last=NA), ]
noFvsO2_sig = noFvsO2[(noFvsO2$padj < 0.05), ]
noFvsO2_sig
noFvsO2_sig_full = cbind(as(noFvsO2_sig, "data.frame"), as(phyloseq::tax_table(rare_ps_genus)[rownames(noFvsO2_sig), ], "matrix"))
nrow(noFvsO2_sig)

# no fertilizer vs. organic_4
noFvsO4=results(ds,contrast=c("fertilizer","no_fertilizer","organic_4"), alpha=0.05)
noFvsO4 = noFvsO4[order(noFvsO4$padj, na.last=NA), ]
noFvsO4_sig = noFvsO4[(noFvsO4$padj < 0.05), ]
noFvsO4_sig
noFvsO4_sig_full = cbind(as(noFvsO4_sig, "data.frame"), as(phyloseq::tax_table(rare_ps_genus)[rownames(noFvsO4_sig), ], "matrix"))
nrow(noFvsO4_sig)

# no fertilizer vs. organic_6
noFvsO6=results(ds,contrast=c("fertilizer","no_fertilizer","organic_6"), alpha=0.05)
noFvsO6 = noFvsO6[order(noFvsO6$padj, na.last=NA), ]
noFvsO6_sig = noFvsO6[(noFvsO6$padj < 0.05), ]
noFvsO6_sig
noFvsO6_sig_full = cbind(as(noFvsO6_sig, "data.frame"), as(phyloseq::tax_table(rare_ps_genus)[rownames(noFvsO6_sig), ], "matrix"))
nrow(noFvsO6_sig)

# mineral vs. organic_2
MvsO2=results(ds,contrast=c("fertilizer","mineral","organic_2"), alpha=0.05)
MvsO2 = MvsO2[order(MvsO2$padj, na.last=NA), ]
MvsO2_sig = MvsO2[(MvsO2$padj < 0.05), ]
MvsO2_sig
MvsO2_sig_full = cbind(as(MvsO2_sig, "data.frame"), as(phyloseq::tax_table(rare_ps_genus)[rownames(MvsO2_sig), ], "matrix"))
nrow(MvsO2_sig)

# mineral vs. organic_4
MvsO4=results(ds,contrast=c("fertilizer","mineral","organic_4"), alpha=0.05)
MvsO4 = MvsO4[order(MvsO4$padj, na.last=NA), ]
MvsO4_sig = MvsO4[(MvsO4$padj < 0.05), ]
MvsO4_sig
MvsO4_sig_full = cbind(as(MvsO4_sig, "data.frame"), as(phyloseq::tax_table(rare_ps_genus)[rownames(MvsO4_sig), ], "matrix"))
nrow(MvsO4_sig)

# mineral vs. organic_6
MvsO6=results(ds,contrast=c("fertilizer","mineral","organic_6"), alpha=0.05)
MvsO6 = MvsO6[order(MvsO6$padj, na.last=NA), ]
MvsO6_sig = MvsO6[(MvsO6$padj < 0.05), ]
MvsO6_sig
MvsO6_sig_full = cbind(as(MvsO6_sig, "data.frame"), as(phyloseq::tax_table(rare_ps_genus)[rownames(MvsO6_sig), ], "matrix"))
nrow(MvsO6_sig)

# organic_2 vs. organic_4
O2vsO4=results(ds,contrast=c("fertilizer","organic_2","organic_4"), alpha=0.05)
O2vsO4 = O2vsO4[order(O2vsO4$padj, na.last=NA), ]
O2vsO4_sig = O2vsO4[(O2vsO4$padj < 0.05), ]
O2vsO4_sig
O2vsO4_sig_full = cbind(as(O2vsO4_sig, "data.frame"), as(phyloseq::tax_table(rare_ps_genus)[rownames(O2vsO4_sig), ], "matrix"))
nrow(O2vsO4_sig)

# organic_2 vs. organic_6
O2vsO6=results(ds,contrast=c("fertilizer","organic_2","organic_6"), alpha=0.05)
O2vsO6 = O2vsO6[order(O2vsO6$padj, na.last=NA), ]
O2vsO6_sig = O2vsO6[(O2vsO6$padj < 0.05), ]
O2vsO6_sig
O2vsO6_sig_full = cbind(as(O2vsO6_sig, "data.frame"), as(phyloseq::tax_table(rare_ps_genus)[rownames(O2vsO6_sig), ], "matrix"))
nrow(O2vsO6_sig)

# organic_4 vs. organic_6
O4vsO6=results(ds,contrast=c("fertilizer","organic_4","organic_6"), alpha=0.05)
O4vsO6 = O4vsO6[order(O4vsO6$padj, na.last=NA), ]
O4vsO6_sig = O4vsO6[(O4vsO6$padj < 0.05), ]
O4vsO6_sig
O4vsO6_sig_full = cbind(as(O4vsO6_sig, "data.frame"), as(phyloseq::tax_table(rare_ps_genus)[rownames(O4vsO6_sig), ], "matrix"))
nrow(O4vsO6_sig)

nrow(noFvsM_sig)
nrow(noFvsO2_sig)
nrow(noFvsO4_sig)
nrow(noFvsO6_sig)
nrow(MvsO2_sig)
nrow(MvsO4_sig)
nrow(MvsO6_sig)
nrow(O2vsO4_sig)
nrow(O2vsO6_sig)
nrow(O4vsO6_sig)

length(noFvsM_sig$log2FoldChange[noFvsM_sig$log2FoldChange<0])
length(noFvsM_sig$log2FoldChange[noFvsM_sig$log2FoldChange>0])
length(noFvsO2_sig$log2FoldChange[noFvsO2_sig$log2FoldChange<0])
length(noFvsO2_sig$log2FoldChange[noFvsO2_sig$log2FoldChange>0])
length(noFvsO4_sig$log2FoldChange[noFvsO4_sig$log2FoldChange<0])
length(noFvsO4_sig$log2FoldChange[noFvsO4_sig$log2FoldChange>0])
length(noFvsO6_sig$log2FoldChange[noFvsO6_sig$log2FoldChange<0])
length(noFvsO6_sig$log2FoldChange[noFvsO6_sig$log2FoldChange>0])
length(MvsO2_sig$log2FoldChange[MvsO2_sig$log2FoldChange<0])
length(MvsO2_sig$log2FoldChange[MvsO2_sig$log2FoldChange>0])
length(MvsO4_sig$log2FoldChange[MvsO4_sig$log2FoldChange<0])
length(MvsO4_sig$log2FoldChange[MvsO4_sig$log2FoldChange>0])
length(MvsO6_sig$log2FoldChange[MvsO6_sig$log2FoldChange<0])
length(MvsO6_sig$log2FoldChange[MvsO6_sig$log2FoldChange>0])
length(O2vsO4_sig$log2FoldChange[O2vsO4_sig$log2FoldChange<0])
length(O2vsO4_sig$log2FoldChange[O2vsO4_sig$log2FoldChange>0])
length(O2vsO6_sig$log2FoldChange[O2vsO6_sig$log2FoldChange<0])
length(O2vsO6_sig$log2FoldChange[O2vsO6_sig$log2FoldChange>0])
length(O4vsO6_sig$log2FoldChange[O4vsO6_sig$log2FoldChange<0])
length(O4vsO6_sig$log2FoldChange[O4vsO6_sig$log2FoldChange>0])


# write.csv(noFvsM_sig_full,file="noFvsM.csv")
# write.csv(noFvsO2_sig_full,file="noFvsO2.csv")
# write.csv(noFvsO4_sig_full,file="noFvsO4.csv")
# write.csv(noFvsO6_sig_full,file="noFvsO6.csv")
# write.csv(MvsO2_sig_full,file="MvsO2.csv")
# write.csv(MvsO4_sig_full,file="MvsO4.csv")
# write.csv(MvsO6_sig_full,file="MvsO6.csv")
# write.csv(O2vsO4_sig_full,file="O2vsO4.csv")
# write.csv(O2vsO6_sig_full,file="O2vsO6.csv")
# write.csv(O4vsO6_sig_full,file="O4vsO6.csv")
# 
nof_m=read.csv(file="noFvsM.csv")
nof_o2=read.csv(file="noFvsO2.csv")
nof_o4=read.csv(file="noFvsO4.csv")
nof_o6=read.csv(file="noFvsO6.csv")
m_o2=read.csv(file="MvsO2.csv")
m_o4=read.csv(file="MvsO4.csv")
m_o6=read.csv(file="MvsO6.csv")
o2_o4=read.csv(file="O2vsO4.csv")
o2_o6=read.csv(file="O2vsO6.csv")
o4_o6=read.csv(file="O4vsO6.csv")


total=merge(nof_m, nof_o2, by="X", all=TRUE)
total=merge(total, nof_o4, by="X", all=TRUE)
total=merge(total, nof_o6, by="X", all=TRUE)
total=merge(total, m_o2, by="X", all=TRUE)
total=merge(total, m_o4, by="X", all=TRUE)
total=merge(total, m_o6, by="X", all=TRUE)
total=merge(total, o2_o4, by="X", all=TRUE)
total=merge(total, o2_o6, by="X", all=TRUE)
total=merge(total, o4_o6, by="X", all=TRUE)

# write.csv(total, file="total.csv")
total_edit=read.csv(file="09.total_edited.csv")    
total_edit_df=as.data.frame(total_edit)
colnames(total_edit_df)
colnames(total_edit_df)=c("ID", "nof_m","nof_o2","nof_o4","nof_o6","m_o2", "m_o4","m_o6", "o2_o4", "o2_o6","o4_o6")

tax_tab=read.csv(file="tax_tab.csv")
colnames(tax_tab)
colnames(tax_tab)=c("ID","Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus" ,  "Species")
head (tax_tab)

total_binary=read.csv(file="09.total_edited_v3.csv")
total_binary[is.na(total_binary)]=""
total_edit_tax=merge(total_binary, tax_tab, by="ID")

nwk <- read.tree("tree.tre")
tips_to_keep=c("'03e58d7cca7e86b7a0f1882b794daaf7'","'0b9cb58df7bf94f4a39c63eec8a40a16'","'1b78072b77acd0ad0a696f710faed818'","'26806e97b60d87f8e699950a0924c620'","'272045d665c4203b3e4530e461953a0c'","'36fc13f12732167cf880e50cb8b9f3c0'","'428bf41da247f0c5a097bef0acf5ccd9'","'4bbd43412cb3b83e4bd7bb8a49e6ec5d'","'4f30775eaacdf9a258b0ce78114dd309'","'50a54ebbff32e375978e5c72416749b9'","'530437e9bcf3bc07e44fcd275a63b1e1'","'5424adb0fd084805e10fa977ca6dee02'","'5d28dc37909b222f704ce5f617902d8a'","'6b43f05bd4d20fb8ac54ce8f8acf0764'","'6fbc11ee6181120518f4c1ec2c14da6c'","'7863e499e498fe835c5c5a74c19e0568'","'7a29003cfb55d8d33b8b0718234ffe05'","'83b3017444d6d1073dde897542f8bf61'","'9c161abbc647a811281dddcf1a0e9d53'","'a9b4afff69b9010904a19568b80ca4eb'","'b20738cc5d13918176859602c10fb841'","'c19eb04f6823c90ded0cd253c4d41382'","'c6f61b95201848ba396fcf8f9f393bd8'","'d97c1bb6b21bed11d2d424df1e8a25ab'","'dde70f26699a7091ec3666372358b497'","'e90dbd62bb0a32729bc98d3ff18bfd2d'","'efbef087504c53a3dc3abd8b5d82188d'","'f0edbb81dae082b85925da3c07085614'")
nwk$tip.label
pruned_tree = keep.tip(nwk,tips_to_keep)

