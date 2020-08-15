# -----------------------------------------------------------------------------#
# Dog Oral Bacteria
# Initial exploratory analyses of processed 16S reads
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     vegan v 2.5.6
#                     ggpubr v 0.4.0
#                     patchwork v 1.0.1
#                     corncob v 0.1.0
#                     broom v 0.7.0
#                     lme4 v 1.1.23
#                     car v 3.0.8
# -----------------------------------------------------------------------------#


# PACKAGES, SCRIPTS, AND SETUP ####
library(phyloseq); packageVersion("phyloseq")
library(tidyverse); packageVersion("tidyverse")
library(vegan); packageVersion("vegan")
library(purrr); packageVersion("purrr")
library(ggpubr); packageVersion("ggpubr")
library(patchwork); packageVersion("patchwork")
library(corncob); packageVersion("corncob")
library(broom); packageVersion("broom")
library(lme4); packageVersion("lme4")
library(car); packageVersion("car")


# library(microbiome)

source("./R/palettes.R")
source("./R/plot_bar2.R")
source("./R/bbdml_helper.R")

#Set ggplot theme
theme_set(theme_bw())

# IMPORT DATA ####
ps <- readRDS("./output/16S_clean_ps_object_w_tree.RDS")


# Investigate SILVA assignments at each taxon level
ps_sp <- ps
phy <- !is.na(tax_table(ps_sp)[,2])
cla <- !is.na(tax_table(ps_sp)[,3])
ord <- !is.na(tax_table(ps_sp)[,4])
fam <- !is.na(tax_table(ps_sp)[,5])
gen <- !is.na(tax_table(ps_sp)[,6])
spp <- !is.na(tax_table(ps_sp)[,7])
assignments_sponge <- data.frame(Phylum=phy, Class=cla,Order=ord,Family=fam,Genus=gen,Species=spp)

assignments_sponge %>% pivot_longer(1:6) %>% mutate(name=factor(name,levels = c("Phylum","Class","Order","Family","Genus","Species"))) %>%
  ggplot(aes(x=name,fill=value)) + geom_bar() + scale_fill_manual(values=c("Gray","Black")) +
  labs(x="Taxonomic level",y="Count",fill="Unambiguous\nassignment")
ggsave("./output/figs/SILVA_Taxonomic_Assignment_Efficiency_at_Each_Taxonomic_Rank.png",dpi=300)
rm(phy,cla,ord,fam,gen,spp,assignments_sponge,ps_sp)

# remove "NA" Phylum taxa
ps <- subset_taxa(ps,!is.na(tax_table(ps)[,2]))

# Look at available metadata
glimpse(sample_data(ps))

#distribution of taxa and sample sums
summary(taxa_sums(ps))
summary(sample_sums(ps))

# quick alpha div plots
plot_richness(ps,x="Pre_or_PostDentalChew",measures = "Shannon") + 
  labs(y="Shannon diversity")
ggsave("./output/figs/16S_Shannon_diversity_dotplot_by_PrePost.png",dpi=300)


# Compare Treat taxa with oral taxa



# Calculate alpha diversity measures and add to metadata
ps@sam_data$Shannon <- estimate_richness(ps, measures="Shannon")$Shannon
ps@sam_data$Richness <- specnumber(otu_table(ps))

# Merge samples for plotting ####
# merge based on sponge species and acidification
ps_merged <- merge_samples(ps,"Pre_or_PostDentalChew")

# Diversity BarPlots ####
# stacked boxplots x-axis=Acidification, y-axis relative-abundance

# Phylum
ps_merged %>%
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Phylum") + scale_fill_manual(values = pal.discrete) + 
  labs(y="Relative abundance",x="Sample type") +
  theme(axis.text.x = element_text(angle = 60,hjust=1),
        axis.title = element_text(face="bold",size=16),
        axis.text = element_text(face="bold",size=12),
        strip.background = element_blank(),
        strip.text = element_text(size=12,face="bold.italic"))

ggsave("./output/figs/16S_Phylum_Diversity_BarChart_by_SampleType.png",dpi=300,height=8,width=12)

# Beta-diversity distances and ordinations ####
unifrac.dist <- UniFrac(ps,weighted = TRUE,normalized = TRUE,parallel = TRUE)

glimpse(sample_data(ps))
ordu = ordinate(ps, "PCoA","unifrac", weighted=TRUE)
plot_ordination(ps, ordu, color="Pre_or_PostDentalChew") +
  geom_point(size=3,alpha=.5) + scale_color_manual(values=pal.discrete) +
  labs(caption = "MDS/PCoA on weighted-UniFrac distance")




# Model alpha diversity ####
# as a function of acidification, species, bleaching status

# Consider, first, the data without Seawater samples...not sure how to treat those...
meta <- microbiome::meta(ps)
glimpse(meta)
meta2 <- meta[meta$Pre_or_PostDentalChew != "Control",]
mod1 <- glm(data = meta2, 
            Richness ~ Pre_or_PostDentalChew)
summary(mod1)

mod2 <- glm(data = meta2, 
            Shannon ~ Pre_or_PostDentalChew)
summary(mod2)

# ANOVA of diversity
mod3 <- aov(data = meta2,
            Shannon ~ Pre_or_PostDentalChew)
summary(mod3)

mod3b <- aov(data = meta2,
             Richness ~ Pre_or_PostDentalChew)
summary(mod3b)

# Dental chews have no clear effect on richness or Shannon diversity
sink("./output/16S_Diversity_linear_models.txt")
print("Richness")
summary(mod1)
summary(mod3b)
print("Shannon")
summary(mod2)
summary(mod3)
sink(NULL)

# Paired differences ####
#remove treat controls
ps2 <- ps %>% subset_samples(Pre_or_PostDentalChew != "Control")
# paired t-test of diversity with individual dogs as grouping variable

# Set up new dataframe
Treatment <- ps2@sam_data$Pre_or_PostDentalChew
Subject <- map_chr(str_split(ps2@sam_data$SampleID,"-"),1)
Age_months <- ps2@sam_data$DogAgeMonths
Gender <- ps2@sam_data$DogGender
Size <- ps2@sam_data$DogSize
Shannon <- ps2@sam_data$Shannon
Richness <- ps2@sam_data$Richness
RegularChewUse <- ps2@sam_data$DogDentalChewUse
DentalIssues <- str_trunc(ps2@sam_data$DogOralIssues,width = 1,ellipsis = "")

meta3 <- data.frame(Treatment=Treatment,
                    Subject=Subject,
                    Age_months=Age_months,
                    Sex=Gender,
                    Size=Size,
                    Shannon=Shannon,
                    Richness=Richness,
                    RegularChewUse=RegularChewUse,
                    DentalIssues=DentalIssues)

mod4 <- lmer(data=meta3,
             Shannon ~ Treatment*RegularChewUse + (1|Subject))
mod5 <- lmer(data=meta3,
             Richness ~ Treatment*RegularChewUse + (1|Subject))

sink("./output/lmer_models.txt")
summary(mod4)
Anova(mod4,test.statistic = "F")
summary(mod5)
Anova(mod5,test.statistic = "F")
sink(NULL)



# Differential abundance of bacterial taxa ####
# Including treat controls
ps_genus <- tax_glom(ps,"Genus")
ps_genus <- ps_genus %>% subset_samples(Pre_or_PostDentalChew != "Control")
set.seed(123)
da_analysis <- differentialTest(formula = ~ Pre_or_PostDentalChew, #abundance
                                phi.formula = ~ Pre_or_PostDentalChew, #dispersion
                                formula_null = ~ 1, #mean
                                phi.formula_null = ~ 1,
                                test = "Wald", boot = FALSE,
                                data = ps_genus,
                                fdr_cutoff = 0.05)


length(da_analysis$significant_models)

bbdml_obj <- multi_bbdml(da_analysis,
                         ps_object = ps_genus,
                         mu_predictor = "Pre_or_PostDentalChew",
                         phi_predictor = "Pre_or_PostDentalChew",
                         taxlevels = 2:7)
length(bbdml_obj)

sink("./output/BBDML_results.txt")
bbdml_obj
sink(NULL)

plot_multi_bbdml(bbdml_list = bbdml_obj,
                 color = "Pre_or_PostDentalChew",
                 pointsize = 3)

plot_1 <- bbdml_plot_1 + theme(legend.position = "none") + labs(y="")
plot_2 <- bbdml_plot_2 + theme(legend.position = "none")   
plot_3 <- bbdml_plot_3 + theme(legend.position = "bottom") + labs(y="")   
  
plot_1 / plot_2 / plot_3

ggsave("./output/figs/DA_Genus-Level_Taxa.png",dpi=300,height = 6,width = 12)


# lmer models on these 3 taxa

taxa_sums(ps_genus)
# add columns to meta3
Staphylococcus <- subset_taxa(ps_genus,Genus == "Staphylococcus") %>% otu_table()
meta3$Staphylococcus <- Staphylococcus@.Data
meta3$Staphylococcus <- meta3$Staphylococcus/sum(meta3$Staphylococcus) %>% as.numeric()
Acinetobacter <- subset_taxa(ps_genus,Genus == "Acinetobacter") %>% otu_table()
meta3$Acinetobacter <- Acinetobacter@.Data
meta3$Acinetobacter <- meta3$Acinetobacter/sum(meta3$Acinetobacter) %>% as.numeric()
Undibacterium <- subset_taxa(ps_genus,Genus == "Undibacterium") %>% otu_table()
meta3$Undibacterium <- Undibacterium@.Data
meta3$Undibacterium <- meta3$Undibacterium/sum(meta3$Undibacterium) %>% as.numeric()

mod.staphy <- lmer(data=meta3,
             Staphylococcus ~ Treatment*RegularChewUse + (1|Subject))
Anova(mod.staphy,test.statistic = "F")

mod.acineto <- lmer(data=meta3,
                   Acinetobacter ~ Treatment*RegularChewUse + (1|Subject))
Anova(mod.acineto,test.statistic = "F")
mod.undi <- lmer(data=meta3,
                   Undibacterium ~ Treatment*RegularChewUse + (1|Subject))
Anova(mod.undi,test.statistic = "F")


# Paired t.test ####
# keep only subjects with pre and post measurements
meta4 <- meta3[meta3$Subject %in% meta3[duplicated(meta3$Subject),]$Subject,]
str(meta4)

Pre <- meta4 %>% filter(Treatment == "Pre") %>% select(Staphylococcus)
Pre <- Pre$Staphylococcus

Post <- meta4 %>% filter(Treatment == "Post") %>% select(Staphylococcus)
Post <- Post$Staphylococcus


t.test(x=Pre,y=Post,paired=TRUE)
