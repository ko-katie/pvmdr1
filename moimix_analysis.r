install.packages("BiocManager")
BiocManager::install("remotes") 
BiocManager::install("bahlolab/moimix", build_vignettes = TRUE)
BiocManager::install("SeqArray")

library(SeqArray)

#convert vcf.gz file into .gds file
seqVCF2GDS("all_samples_variants_filtered2.vcf.gz", "all_samples_variants_filtered2.gds")

#read .gds file into R and check that it matches original VCF file
samples <- seqOpen("all_samples_variants_filtered2.gds")
seqSummary(samples)
samples.id <- seqGetData(samples, "sample.id")

library(moimix)

#perform moimix analysis
coords <- getCoordinates(samples)
head(coords)

#Filter out variants on apicoplast and mitochondria because they skew the results
seqSetFilter(samples, variant.id = coords$variant.id[coords$chromosome != "PvP01_API_v1" & coords$chromosome !="PvP01_MIT_v1"])

#estimate Fws
fws_all <- getFws(samples)
hist(fws_all, breaks=20, xlim=c(0,1.0))

sum(fws_all < 0.95)

fws_results <- data.frame(fws_all)
rownames(fws_results) <- samples.id
write.csv(fws_results, file="fws_results.csv")
