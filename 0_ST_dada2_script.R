#!/usr/bin/env Rscript
#title: "Dada2_16S_EO4"
#original author: "Jacob Agerbo Rasmussen"
#created: "8/18/2020"
#edited by: "Sigri Kløve"

#set dependencies
library("dada2")
#packageVersion("dada2") # 1.11.5 when this was put together
#list.files() # make sure what we think is here is actually here

#r set create data and check QC
## first we're setting a few variables we're going to use ##
# one with all sample names, by scanning our "samples" file we made earlier
path <- "EO5"
# one holding the file names of all the forward reads
fnFs <- sort(list.files(path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE))
# and one with the reverse
fnRs <- sort(list.files(path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_L001_R1_001.fastq.gz"), `[`, 1)

#Quality trimming/filtering
pdf("QC_profile.pdf", height = 15, width = 15)
plotQualityProfile(fnFs[1:8]) #problems to show all together
plotQualityProfile(fnRs[1:8]) #more dramatic drop in reverse reads
dev.off()

# Description for filtering, using DADA2
#INPUT: forward_reads, reverse_reads
#OUTPUT: filtered_forward_reads, filtered_reverse_reads:
# maxEE=Maximum amount of estimated errors that you expect to see and read
#rm.phix: removes any reads that match the PhiX bacteriophage genome, which is typically added to Illumina sequencing runs for quality monitoring
#truncLen: parameter setting the minimum size to trim the forward and reverse reads to in order to keep the quality scores roughly above 30 overall
#minLen: is setting the minimum length reads we want to keep after trimming
#no minLen value in dada2 tutorial
#maxN: additional filtering default parameter that is removing any sequences containing any Ns
#truncq= truncates your reads at the first base pair with a quality score equal to or less than 2
#compress=TRUE: FASTQ files gzipped
#multithread=TRUE: FASTQ files process in parallel

#Filtering
filtFs <- file.path(path, "filtered_dec24", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered_dec24", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

filtered_out = filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                    truncLen=c(280,250), #decide on numbers based on Quality Profile
                    trimLeft = c(21,25), #Remove primers if they are not already removed 
                    maxN=0, 
                    maxEE=c(2,5), 
                    truncQ=2, 
                    rm.phix=TRUE, 
                    compress=TRUE, 
                    multithread=TRUE)

### Store filtered reads as R dataset
saveRDS(filtered_out, "filtered_out.RData")

pdf("QC_profile_filtered.pdf", height = 15, width = 15)
plotQualityProfile(filtFs[1:8]) #problems to show all together
plotQualityProfile(filtRs[1:8]) #more dramatic drop in reverse reads
dev.off()

#r filtering based on error rates
#Generating an error model of our data by learning the specific error-signature of our dataset
err_forward_reads <- learnErrors(filtFs, multithread=TRUE)
err_reverse_reads <- learnErrors(filtRs, multithread=TRUE)
#The developers have incorporated a plotting function to visualize how well the estimated error rates match up with the observed:
#The red line is what is expected based on the quality score, the black line represents the estimate, and the black dots represent the observed.
#In generally speaking, you want the observed (black dots) to track well with the estimated (black line)
#In geneal, error rate decreases as Q scores increases
#QC filtered reads

pdf("plotErrors.pdf", height = 15, width = 15)
plotErrors(err_forward_reads, nominalQ=TRUE)
plotErrors(err_reverse_reads, nominalQ=TRUE)
dev.off()

#Dereplication:keep/process one, and just attach the identical sequences to it
#When DADA2 dereplicates sequences, it also generates a new quality-score profile of each unique sequence
#based on the average quality scores of each base of all of the sequences that were replicates of it.
#the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow


#QC filtered reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

#Inferring of ASVs
#It does this by incorporating the consensus quality profiles and abundances of each unique sequence, and then figuring out if each sequence more likely to be of biological origin or more likely to be spurious.

#Inferring ASVs
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE)

#maxMismatch is by default 0 but you can make it say 1 or 2 if you’re finding that a lot of your forward and reverse reads are not merging
#minOverlap option has been change because they were sequenced at 250 bp

#Merge reads
#Merging forward and reverse reads
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, minOverlap=12, verbose=TRUE)
### Store filtered amplicons as R dataset
saveRDS(merged_amplicons, "merged_amplicons.RData")

###Generate ASV table
#Generating a count table
seqtab <- makeSequenceTable(merged_amplicons)
table(nchar(getSequences(seqtab)))
#to see the breakdown of the size of these amplicons. On the top row is the size of the merged reads, and on the botton is the frequency

#Chimera identification
seqtab.nochim <- removeBimeraDenovo(seqtab, multithread=TRUE, verbose=TRUE)
saveRDS(seqtab.nochim, "seqtab.RData")

### Create a Summary
# this is one quick way to look at sequences that have been lost, to know whether they held a lot in terms of abundance
sum(seqtab.nochim)/sum(seqtab)

#Overview of counts throughout:quick way to pull out how many reads were dropped at various points of the pipeline
# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))
write.table(summary_tab, "summary_reads_table.txt", sep="\t", quote=F)

#Assigning taxonomy(before this step, download the database, at: https://zenodo.org/record/1172783#.XzvzQJMzbmE)
#                  There are different DADA2-formatted databases available in DADA2 website
                   
#Assign Taxanomy
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz", tryRC=T, multithread=TRUE)
saveRDS(taxa, "taxa.RData")

                   
#Extracting the standard goods from DADA2
#Write out tables for further processing}
#giving to seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
                   
#making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")
                   
#count table:
asv_tab <- t(seqtab.nochim)
write.table(asv_tab, "ASVs_counts.csv", sep=",", quote=F, col.names=NA)
                   
#tax table:
asv_tax <- taxa
write.table(asv_tax, "ASVs_taxonomy.csv", sep=",", quote=F, col.names=NA)

