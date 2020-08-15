# -----------------------------------------------------------------------------#
# Dog Oral Bacteria
# Processing trimmed Reads
# Author: Geoffrey Zahn
# Software versions:  R v 3.6.3
#                     tidyverse v 1.3.0
#                     dada2 v 1.14.1
#                     decontam v 1.6.0
#                     phyloseq v 1.30.0
#                     purrr v 0.3.4
#                     Biostrings 2.54.0
# -----------------------------------------------------------------------------#


# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dada2); packageVersion("dada2")
library(decontam); packageVersion("decontam")
library(phyloseq); packageVersion("phyloseq")
library(purrr); packageVersion("purrr")
library(Biostrings); packageVersion("Biostrings")
library(readxl)

source("./R/palettes.R")
source("./R/plot_bar2.R")


#################################################################################
#                               Main workflow                                   #
# Filter and trim, denoise, sample inferrence, chimera and contaminant removal, # 
# taxonomic assignment, phylogenetic tree, combine sequence table and metadata  #
#################################################################################

# PARSE FILE PATHS ####

# File parsing - For this, we will use only the forward illumina reads - make sure to move fwd reads into their own directory for simplest processing
path <- "./seqs" # CHANGE to the directory containing your demultiplexed fastq files
filtpath <- file.path(path, "filtered") # Filtered files go into the filtered/ subdirectory
if(!file_test("-d", filtpath)) dir.create(filtpath) # make directory for filtered fqs if not already present
fns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R1_001.fastq.gz"))
rns <- sort(list.files(file.path(path), full.names = TRUE, pattern = "R2_001.fastq.gz"))


sample.names <- unlist(map(strsplit(basename(fns), "_"), 1))

# visualize a couple of fwd read quality profiles to help select reasonable filtration parameters
plotQualityProfile(fns[1:2])
plotQualityProfile(rns[1:2])

# FILTER AND TRIM ####
filts_f <- file.path(path, "filtered", paste0(sample.names, "_FWD_filt.fastq.gz"))
filts_r <- file.path(path, "filtered", paste0(sample.names, "_REV_filt.fastq.gz"))

out <- filterAndTrim(fns, filts_f, rns, filts_r,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, truncLen = c(300,200),
                     compress=TRUE, multithread=4) # On Windows set multithread=FALSE
head(out,10)
nrow(out)
saveRDS(out, "./output/trackreads.RDS")

# Some input samples had no reads pass the filter!


# sanity check  comparison of before and after filtration
plotQualityProfile(c(fns[1:2],filts_f[1:2]))
plotQualityProfile(c(rns[1:2],filts_r[1:2]))

# LEARN ERROR RATES ####
# Since some samples may have had zero reads pass QC, reassign filts
filts_f <- sort(list.files(filtpath, full.names = TRUE,pattern = "FWD"))
filts_r <- sort(list.files(filtpath, full.names = TRUE,pattern = "REV"))
length(fns);length(filts_f)
length(rns);length(filts_r) 

# learn errors
set.seed(123)
errF <- learnErrors(filts_f, multithread=TRUE, MAX_CONSIST = 20,verbose = 2,randomize = TRUE)
saveRDS(errF,"./output/errF.RDS")
set.seed(123)
errR <- learnErrors(filts_r, multithread=TRUE, MAX_CONSIST = 20,verbose = 2,randomize = TRUE)
saveRDS(errR,"./output/errR.RDS")

# sanity check for error model
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# DEREPLICATION ####
derepF <- derepFastq(filts_f, verbose=TRUE)
derepR <- derepFastq(filts_r, verbose=TRUE)


# Name the derep-class objects by the sample names
# If some samples were removed (no reads passed QC), reassign sample.names
if(length(derepF) != length(sample.names)){
  sample.names <- unlist(map(strsplit(basename(filts_f), "_filt"), 1))
}


if(identical(unlist(map(strsplit(basename(filts_f), "FWD_filt"), 1)),unlist(map(strsplit(basename(filts_r), "REV_filt"), 1)))){
  names(derepF) <- sample.names
  names(derepR) <- sample.names
} else {
  stop("Make sure fwd and rev files are in same order!")
}  



# SAMPLE INFERRENCE ####
dadaFs <- dada(derepF, err=errF, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")
dadaRs <- dada(derepR, err=errR, multithread=TRUE, selfConsist = TRUE, verbose=TRUE, pool = "pseudo")



mergers <- mergePairs(dadaFs, filts_f, dadaRs, filts_r, verbose=TRUE)

# MAKE SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)

# REMOVE CHIMERAS ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim,"./output/seqtab.nochim.RDS")
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# reassign "out" to remove any missing reads
out = out[as.data.frame(out)$reads.out > 0,]

# TRACK READS THROUGH PIPELINE ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track = as.data.frame(track)
track$total.loss.proportion = (track[,1]-track$nonchim)/track[,1]
head(track)

write.csv(track, file = "./output/16S_read_counts_at_each_step.csv", row.names = TRUE)

row.names(seqtab)
# IMPORT METADATA ####
meta = read_xlsx("./data/DM_Metadata.xlsx")
row.names(meta) <- as.character(meta$SampleID)

# reorder metadata to match seqtab
df <- data.frame(seqtab_rows=row.names(seqtab.nochim),
                 SampleID=row.names(seqtab.nochim))
df2 <- left_join(df,meta,by="SampleID")
row.names(df2) <- df2$SampleID
identical(row.names(df2),row.names(seqtab.nochim))


# Remove all seqs with fewer than 100 nucleotides ####
keeper_esvs <- nchar(names(as.data.frame(seqtab.nochim))) > 99
seqtab.nochim <- seqtab.nochim[,keeper_esvs]

# remove singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]
# remove samples with fewer than 1000 reads
seqtab.nochim <- seqtab.nochim[rowSums(seqtab.nochim) >= 1000,]
# remove newly singleton taxa
seqtab.nochim <- seqtab.nochim[,colSums(seqtab.nochim) > 1]

# save cleaned up seqtab
saveRDS(seqtab.nochim,"./output/seqtab.nochim.clean.RDS")

# Find and remove contaminants
meta$Control <- meta$SampleSource == "NegativeControl"
meta <- meta[meta$SampleID %in% row.names(seqtab.nochim),] 
contams = isContaminant(seqtab.nochim, neg = meta$Control, normalize = TRUE)
table(contams$contaminant)
write.csv(contams, file = "./output/likely_contaminants.csv", row.names = TRUE)
seqtab.nochim = seqtab.nochim[,(which(contams$contaminant != TRUE))]
seqtab.nochim = seqtab.nochim[meta$Control == FALSE,]
meta = meta[meta$Control == FALSE,]



# ASSIGN TAXONOMY ####
taxa <- assignTaxonomy(seqtab.nochim, "./taxonomy/rdp_train_set_16.fa.gz", multithread=4)

# Save intermediate taxonomy file
saveRDS(taxa, file = "./output/RDP_Taxonomy_from_dada2.RDS")

# add_species
taxa <- addSpecies(taxa, "./taxonomy/rdp_species_assignment_16.fa.gz")

# Save completed taxonomy file
saveRDS(taxa, file = "./output/RDP_Taxonomy_from_dada2_sp.RDS")


# inspect taxonomy
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# re-load point
# seqtab.nochim <- readRDS("./output/dada2_seqtable.RDS")
# taxa <- readRDS("./output/RDP_Taxonomy_from_dada2.RDS")
# meta <- read_delim("./data/Metadata.csv",delim = ",")
# row.names(meta) <- as.character(meta$SampleID)



# Hand off to Phyloseq ####
otu <- otu_table(seqtab.nochim,taxa_are_rows = FALSE)
tax <- tax_table(taxa)
met <- sample_data(meta)
row.names(met) <- meta$SampleID


ps <- phyloseq(otu,met,tax)
saveRDS(ps,"./output/ps_not-cleaned.RDS")


# Find non-bacteria ####
ps_nonbact <- subset_taxa(ps, Kingdom != "Bacteria")

ps %>% transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Kingdom")
ggsave("./output/figs/16S_Kingdom-Level_Taxonomic_Proportions.png",dpi=300)

ps_nonbact %>% 
  transform_sample_counts(function(x){x/sum(x)}) %>%
  plot_bar2(fill="Kingdom")


# REMOVE NON-BACTERIA, CHLOROPLASTS, MITOCHONDRIA, and empty samples/taxa ####
ps <- subset_taxa(ps, Kingdom == "Bacteria")
ps <- subset_taxa(ps,Class != "Chloroplast")
ps <- subset_taxa(ps, taxa_sums(ps) > 0)
ps <- subset_samples(ps, sample_sums(ps) > 0)

# Save DNA sequences apart from rownames (from subsetted ps object)
seqs <- taxa_names(ps)
seqs <- DNAStringSet(seqs)
saveRDS(seqs,"./output/16S_ASV_reference_sequences.RDS")

# Save RDS object for Phyloseq
saveRDS(ps, file = "./output/clean_phyloseq_object.RDS")


