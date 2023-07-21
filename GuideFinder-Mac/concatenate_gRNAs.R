# read in files and add PAM column
Sepi6_NGG_guides <- read_csv("S_epi_6/Sepi6_NGG_Prom_AllUseableGuides.csv") %>%
  mutate(PAM = "NGG")
Sepi6_NGG_lowered_guides <- read_csv("S_epi_6/Sepi6_NGG_Prom_GuidesUsingLowerThresholds.csv") %>%
  mutate(PAM = "NGG")
Sepi6_NAN_guides <- read_csv("S_epi_6/Sepi6_NAN_Prom_AllUseableGuides.csv") %>%
  mutate(PAM = "NAN")
Sepi6_NAN_lowered_guides <- read_csv("S_epi_6/Sepi6_NAN_Prom_GuidesUsingLowerThresholds.csv") %>%
  mutate(PAM = "NAN")

# take out PAM sequence from guides (in all useable guides files)
Sepi6_NGG_guides <- Sepi6_NGG_guides %>%
  mutate(Guide = substr(Guide,1,(nchar(Guide) - 3)))

Sepi6_NAN_guides <- Sepi6_NAN_guides %>%
  mutate(Guide = substr(Guide,1,(nchar(Guide) - 3)))

Sepi6_all_guides <- rbind(Sepi6_NGG_guides, Sepi6_NAN_guides,
                          Sepi6_NGG_lowered_guides, Sepi6_NAN_lowered_guides) %>%
  select(-...1)

# between 1 and 129 guides per gene
num_guides_per_gene <- Sepi6_all_guides %>% group_by(GeneProduct, PAM) %>% summarize(count=n())

# average of 20.05021 genes per guides
mean_guides_per_gene <- mean(num_guides_per_gene$count)

# Twist oligo prep - add primer/etc sequence around guide
Sepi6_Twist_guides <- Sepi6_all_guides %>%
  mutate(twist_guide = paste0("GCTGTTTTGAATGGTCCCGGTCTCGAAAC", Guide, "GTTTTTGAGACCGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGG")) %>%
  mutate(ID_loc = paste0(ID, "_", PAM, "_", Location)) %>%
  select(SequenceName = ID_loc, Sequence = twist_guide)

# Twist guide file
write.csv(Sepi6_Twist_guides, file="Sepi6_Twist_oligos.csv", row.names=F, quote=F)

# Agilent oligo prep
Sepi6_Agilent_guides <- Sepi6_all_guides %>%
  mutate(ID_loc = paste0(ID, "_", PAM, "_", Location)) %>%
  mutate(Replication = 1) %>%
  mutate(Sequence = paste0("GCTGTTTTGAATGGTCCCGGTCTCGAAAC", Guide, "GTTTTTGAGACCGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGG")) %>%
  select(SequenceName = ID_loc, Sequence, Replication)

# one column Agilent file
write.table(Sepi6_Agilent_guides[,2], file="Sepi6_Agilent_oligos_onecol.txt", row.names=F, col.names=T, quote=F)

# three column Agilent file
write.csv(Sepi6_Agilent_guides, file="Sepi6_Agilent_oligos_threecol.csv", row.names=F, quote=F)
