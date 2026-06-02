# Prepare protein interanction networks biogrid and intact 
# latter is tiny for mouse alone
# Very conservative approach, only mouse, only direct and physical interactions
# Save all metadata and entrez ids only of the interactions for using with epimod

library("tidyverse")
library("here")
library("biomaRt")

results_dir <- here("results")
dirOut <- file.path(results_dir,"1-output")
if(!dir.exists(dirOut)){dir.create(dirOut, recursive = TRUE)}

# PPI from biogrid #----
# https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-5.0.255/BIOGRID-ORGANISM-5.0.255.mitab.zip
# Extracted mouse file
biogrid <- read.delim("data/BIOGRID-ORGANISM-Mus_musculus-5.0.255.mitab.txt", stringsAsFactors = FALSE)

# Filter interactions and species
#table(biogrid$Interaction.Types)
biogrid <- biogrid %>% 
  dplyr::filter(Interaction.Types %in% c("psi-mi:MI:0407(direct interaction)", "psi-mi:MI:0915(physical association)")) %>%
  dplyr::filter(Taxid.Interactor.A %in% c("taxid:10090") & Taxid.Interactor.B %in% c("taxid:10090")) # could consider to keep also humans (taxid:9606) in, but then have to do uniprot/entrez/genename mapping separately

# Get Entrez ids biogrid
extract_entrez <- function(x) {
  sub(".*locuslink:([0-9]+).*", "\\1", x)
}

biogrid$A <- extract_entrez(biogrid$X.ID.Interactor.A)
biogrid$B <- extract_entrez(biogrid$ID.Interactor.B)

saveRDS(biogrid, file = file.path(dirOut, "biogrid_mouse_ppi.Rds"))

# Keep only valid rows
biogrid_clean <- biogrid[
  grepl("^[0-9]+$", biogrid$A) &
    grepl("^[0-9]+$", biogrid$B),
  c("A", "B")
] %>%
  distinct()

saveRDS(biogrid_clean, file = file.path(dirOut, "biogrid_mouse_ppi_entrez.Rds"))

# PPI from IntAct #----
# https://www.ebi.ac.uk/intact/interactomes
# miTab 2.7
intact <- read.delim("data/intact_2.7_mouse.txt", stringsAsFactors = FALSE)

#table(intact$Interaction.type.s.)

intact_ppi <- intact %>% 
  dplyr::filter(Taxid.interactor.A == "taxid:10090(mouse)|taxid:10090(Mus musculus)" & Taxid.interactor.B == "taxid:10090(mouse)|taxid:10090(Mus musculus)") %>% # mouse only
  dplyr::filter(Type.s..interactor.A == "psi-mi:MI:0326(protein)" & Type.s..interactor.B == "psi-mi:MI:0326(protein)") %>%
  dplyr::filter(Interaction.type.s. %in% c("psi-mi:MI:0407(direct interaction)","psi-mi:MI:0915(physical association)")) %>%
  droplevels() # tiny sset

# Entrez IDs
intact_ppi$A_uniprot <- gsub("uniprotkb:","",intact_ppi$X.ID.s..interactor.A)
intact_ppi$B_uniprot <- gsub("uniprotkb:","",intact_ppi$ID.s..interactor.B)

# intact_ppi$A <- mapIds(
#   org.Mm.eg.db,
#   keys = intact_ppi$A_uniprot,
#   column = "ENTREZID",
#   keytype = "UNIPROT",
#   multiVals = "first"   # or "list" if multiple mappings exist
# )
# # too many NA
# sum(is.na(intact_ppi$A))
# [1] 23239

saveRDS(intact_ppi, file = file.path(dirOut, "intact_mouse_ppi.Rds"))

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "mmusculus_gene_ensembl"
)

#attrs <- listAttributes(ensembl)

# Map A UniProt IDs to Entrez
mapping <- getBM(
  attributes = c("uniprotswissprot", "entrezgene_id"),
  filters = "uniprotswissprot",
  values = intact_ppi$A_uniprot,
  mart = ensembl
)

colnames(mapping) <- c("A_uniprot","A")
intact_ppi <- left_join(intact_ppi, mapping, by = "A_uniprot", relationship = "many-to-many") # not a one-to-one mapping

# Map B
mapping <- getBM(
  attributes = c("uniprotswissprot", "entrezgene_id"),
  filters = "uniprotswissprot",
  values = intact_ppi$B_uniprot,
  mart = ensembl
)

colnames(mapping) <- c("B_uniprot","B")
intact_ppi <- left_join(intact_ppi, mapping, by = "B_uniprot", relationship = "many-to-many") # not a one-to-one mapping
# sum(is.na(intact_ppi$A))
# [1] 2694
# Even more NAs

# Keep only valid interactions
intact_clean <- intact_ppi %>%
  dplyr::select(A,B) %>%
  na.omit() %>%
  distinct()

saveRDS(intact_clean, file = file.path(dirOut, "intact_mouse_ppi_entrez.Rds"))

