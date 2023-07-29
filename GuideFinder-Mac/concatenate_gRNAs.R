# read in files and add PAM column
Sepi6_NGG_guides <- read_csv("S_epi_6/Sepi6_NGG_Prom_AllUseableGuides.csv") %>%
  mutate(PAM = "NGG")
Sepi6_NGG_lowered_guides <- read_csv("S_epi_6/Sepi6_NGG_Prom_GuidesUsingLowerThresholds.csv") %>%
  mutate(PAM = "NGG")
# Sepi6_NAN_guides <- read_csv("S_epi_6/Sepi6_NAN_Prom_1_AllUseableGuides.csv") %>%
#   mutate(PAM = "NAN")
# Sepi6_NAN_lowered_guides <- read_csv("S_epi_6/Sepi6_NAN_Prom_1_GuidesUsingLowerThresholds.csv") %>%
#   mutate(PAM = "NAN")


# take out PAM sequence from guides (in all useable guides files)
Sepi6_NGG_guides <- Sepi6_NGG_guides %>%
  mutate(Guide = substr(Guide,1,(nchar(Guide) - 3)))

# Sepi6_NAN_guides <- Sepi6_NAN_guides %>%
#   mutate(Guide = substr(Guide,1,(nchar(Guide) - 3)))

# put all guides in one data frame
# give each guide its own ID
# total 8209 guides; 2298 genes
Sepi6_all_guides <- rbind(Sepi6_NGG_guides, Sepi6_NGG_lowered_guides) %>% 
  arrange(ID, Location) %>%
  mutate(guide_num = c(1:length(ID))) %>%
  select(-...1)

# checking for duplicates
# duplicate guides most likely result from promoter sequence overlapping with gene upstream
Sepi6_NGG_guides %>% group_by(Guide) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head()
Sepi6_NGG_lowered_guides %>% group_by(Guide) %>% summarize(count=n()) %>% arrange(desc(count)) %>% head()

# 52 duplicate guides
duplicate_guides <- Sepi6_all_guides %>% group_by(Guide) %>% summarize(count=n()) %>% filter(count>1)
duplicate_guides_all <- Sepi6_all_guides %>% filter(Guide %in% duplicate_guides$Guide) %>% arrange(ID, PAM, Location)

# read in CDS and promoter data - just need a data table that says whether genes are on plus or minus strand
Sepi6_allCDS_Prom <- read_csv("S_epi_6/Sepi6_allCDS_andProm.csv")

#Sepi6_completeGuidesList <- rbind(read_csv("S_epi_6/Sepi6_NGG_Prom_CompleteGuidesList.csv"), read_csv("S_epi_6/Sepi6_NAN_CompleteGuidesList.csv"))

# sort through duplicate guides
# creates a column to store T or F based on whether to keep the duplicate guide
duplicate_guides_all <- duplicate_guides_all %>%
  arrange(Guide) %>%
  mutate(Keep = NA)

# filter out guides that shouldn't have passed off-target BLAST search
# found by manually looking through list for which pairs of guides are not from neighboring genes
take_out <- which(duplicate_guides_all$ID %in% c(1849, 1397))
duplicate_guides_all$Keep[take_out] <- FALSE

# filter out guides where promoter overlaps with gene body of neighboring gene
# make sure duplicate guides is sorted by guide so that pairs are together
# if genes are on minus strand, keep second guide
for (i in c(1:length(duplicate_guides_all$ID))[c(T,F)]) {
  if (!is.na(duplicate_guides_all$Keep[i])) {next}
  
  ID_first <- duplicate_guides_all$ID[i]
  strand <- Sepi6_allCDS_Prom[which(Sepi6_allCDS_Prom$ID==ID_first),"Strand"]
  
  if (strand=="+") {
    duplicate_guides_all$Keep[i] <- T
    duplicate_guides_all$Keep[i+1] <- F
  } else if (strand=="-") {
    duplicate_guides_all$Keep[i] <- F
    duplicate_guides_all$Keep[i+1] <- T
  }
}

# filter out 27/52 duplicate guides
duplicates_to_drop <-  duplicate_guides_all[!duplicate_guides_all$Keep,]

# filter out guides to drop
# renumber guides
# total 8182 guides
Sepi6_all_guides_filtered <-  Sepi6_all_guides %>%
  filter(!(guide_num %in% duplicates_to_drop$guide_num)) %>%
  arrange(ID, Location) %>%
  mutate(guide_num = c(1:length(ID)))

# between 1 and 24 guides per gene
# 2291 genes represented
num_guides_per_gene <- Sepi6_all_guides_filtered %>% group_by(GeneProduct) %>% summarize(count=n())

# average of 3.571366 genes per guide
mean_guides_per_gene <- mean(num_guides_per_gene$count)

# Twist oligo prep - add primer/etc sequence around guide
Sepi6_Twist_guides <- Sepi6_all_guides_filtered %>%
  mutate(twist_guide = paste0("GCTGTTTTGAATGGTCCCGGTCTCGAAAC", Guide, "GTTTTTGAGACCGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGG")) %>%
  mutate(ID_loc = paste0(guide_num, "_", PAM, "_", ID, "_", Location)) %>%
  select(SequenceName = ID_loc, Sequence = twist_guide)

# Twist guide file
write.csv(Sepi6_Twist_guides, file="Sepi6_Twist_oligos.csv", row.names=F, quote=F)

# Agilent oligo prep
Sepi6_Agilent_guides <- Sepi6_all_guides_filtered %>%
  mutate(ID_loc = paste0(guide_num, "_", PAM, "_", ID, "_", Location)) %>%
  mutate(Replication = 1) %>%
  mutate(Sequence = paste0("GCTGTTTTGAATGGTCCCGGTCTCGAAAC", Guide, "GTTTTTGAGACCGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGG")) %>%
  select(SequenceName = ID_loc, Sequence, Replication)

# one column Agilent file
write.table(Sepi6_Agilent_guides[,2], file="Sepi6_Agilent_oligos_onecol.txt", row.names=F, col.names=T, quote=F)

# three column Agilent file
write.csv(Sepi6_Agilent_guides, file="Sepi6_Agilent_oligos_threecol.csv", row.names=F, quote=F)
