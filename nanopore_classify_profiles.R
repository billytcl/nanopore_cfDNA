#nanopore_classify_profiles.R
#this file is meant to be used within Rstudio

#this script takes in two methylation profiles (tumor and immune), merges them, and then takes another set of methylation profiles (eg. plasma) and computes the number of matching reads to a tumor profile

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(readxl)

#load reference methylation profiles
#these should be sorted by position
immune_profile <- fread("/mnt/ix1/Seq_Runs/20220512_PRM_1224/megalodon/modified_bases.5mC.sorted.bed") #P6199
tumor_profile <- fread("/mnt/ix1/Projects/M083_20211228_nanopore_methylation_GI/longitudinal_pts/00_samples/P6199/P6199_primary_modified_bases.5mC.sorted.bed")

#merge and filter on V5>0 as specified in BedMethyl formatting -- these are zero coverage sites
merge_profile <- inner_join(immune_profile %>% filter(V5 > 0), tumor_profile %>% filter(V5 > 0), by=c("V1","V2","V3")) %>% select(V1,V2,V3,V10.x,V11.x,V10.y,V11.y)

merge_profile_filter <- merge_profile %>% 
  filter(V10.x > 4 & V10.y > 4)

#get list of dumped bam files
pat <- "P6199"
files_tumor <- Sys.glob(paste("/mnt/ix1/Projects/M083_20211228_nanopore_methylation_GI/longitudinal_pts/00_samples/",pat,"/*dump.bed",sep=""))
files_tumor <- files_tumor[grep(pattern = "_normal_|_primary_", x=files_tumor, invert=T)]

#classify dumped files
out <- list()
for (i in 1:length(files_tumor)) {
  print(i)
  file <- files_tumor[i]
  sample <- basename(file)
  sample <- gsub(pattern = "_mod_mappings.dump.bed", replacement = "", x = sample)
  
  file_profile <- fread(file)
  
  print("intersecting")
  profiles_intersect <- file_profile %>% inner_join(merge_profile_filter, by=c("V1","V2","V3"))
  
  #compute a matching probability for each site
  profiles_intersect <- profiles_intersect %>% mutate(
    test_t = case_when(
      (V4 == 0) ~ (100-V11.y),
      (V4 == 1) ~ (V11.y)
    ),
    test_p = case_when(
      (V4 == 0) ~ (100-V11.x),
      (V4 == 1) ~ (V11.x)
    )
  )
  
  print("summarizing")
  #group by read and compute an average score across the read
  profiles_intersect_summary <- profiles_intersect %>% group_by(V9) %>% summarize(s=sum(test_t)/n(), sp=sum(test_p)/n()) %>% as.data.table()
  
  #normalize the score
  profiles_intersect_summary$ratio <- profiles_intersect_summary$s/(profiles_intersect_summary$s + profiles_intersect_summary$sp)
  
  #thresholding based on immune < 0.1 and cancer > 0.9
  nt <- profiles_intersect_summary %>% filter(ratio > 0.9) %>% as.data.frame() %>% nrow()
  np <- profiles_intersect_summary %>% filter(ratio < 0.1) %>% as.data.frame() %>% nrow()
  
  #output result
  out[[i]] <- c(sample, nt, np, nrow(profiles_intersect_summary), length(unique(file_profile$V9)))
}

#transform into dataframe and ensure type
out_df <- as.data.frame(do.call("rbind",out))
out_df$V2 <- as.numeric(out_df$V2)
out_df$V3 <- as.numeric(out_df$V3)
out_df$V4 <- as.numeric(out_df$V4)
out_df$V5 <- as.numeric(out_df$V5)
out_df$i <- 1:nrow(out_df)

#output a flat file because this takes a long time
fwrite(out_df, paste(pat,"_longitudinal_dump_reads_table_2.txt", sep=""), col.names=F, row.names=F, sep="\t", quote=F)

#use a separate excel file containing draw dates -- not necessarily needed for plotting
out_df <- out_df %>% separate(V1, into=c("pt","sample"), sep="_")
draw_dates <- read_excel("draw_dates.xlsx", sheet = "P6199")

out_df <- out_df %>% left_join(draw_dates, by=c("sample"="FromSample")) %>% select(pt,sample,V2,V3,V4,V5,i,SampleDate)
out_df <- out_df %>% mutate(cumldate = difftime(SampleDate, first(SampleDate), units="days"))
out_df$f <- out_df$V2/(out_df$V4)

#ensure sample column naming
colnames(out_df) <- c("pt","sample","V2","V3","V4","V5","i","SampleDate","cumldate","f")


#plotting
out_df_melt <- melt(out_df, measure.vars=c("f","V4"))
out_df_melt$variable <- factor(out_df_melt$variable, levels=c("V4","f"))

labels <- c(
  "V4"="cfDNA burden (reads)",
  "f" = "Frac. reads confidently classified as tumor-specific"
)

ggplot(out_df_melt, aes(x=cumldate,y=value, fill=variable)) + 
  facet_wrap(~variable, scales="free", labeller=as_labeller(labels), nrow=2) + 
  geom_point(pch=21, color="black", size=5) +
  geom_line(linetype="dashed", size=1, alpha=0.2) +
  theme_classic(40) +
  xlab("Cumulative draw date (days)") +
  theme(legend.position="none",
        strip.background=element_blank(),
        axis.title.y = element_blank())

ggsave(paste(pat,"_longitudinal_dumpreads.png", sep=""), dpi=300, width=18, height = 9, units="in")

