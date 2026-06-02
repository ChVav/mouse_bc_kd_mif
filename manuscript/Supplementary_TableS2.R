# Write Excel table of final core gene sets used for PPR analyses
library(openxlsx)
library(tibble)
library(dplyr)
library(stringr)
library(here)

# Output file
out_file <- here("results/TableS1_cancerhubgenes.xlsx")

# Final gene sets used
core_sets <- list(
  "PR regulatory network" = c(
    "Pgr","Ncoa1","Ncoa2","Ncoa3","Ncor1","Ncor2","Pelp1","Kdm4b",
    "Ep300","Crebbp","Med1","Kdm6b","Hsp90aa1","Hsp90ab1","Fkbp4","Fkbp5"
  ),
  
  "RANKL–NF-κB" = c(
    "Tnfsf11","Tnfrsf11a","Tnfrsf11b","Traf6","Nfkb1","Rela","Chuk",
    "Ikbkb","Ikbkg","Nfkb2","Relb","Map3k14","Map3k7","Tab2","Traf2",
    "Traf3","Nfkbia","Birc2","Birc3"
  ),
  
  "Wnt/β-catenin signaling" = c(
    "Wnt4","Rspo1","Lgr4","Lgr5","Fzd2","Fzd7","Lrp5","Lrp6","Dvl2",
    "Ctnnb1","Tcf7","Tcf7l2","Lef1","Axin2","Porcn","Gsk3b","Apc",
    "Csnk1a1","Rnf43","Znrf3"
  ),
  
  "PRL–JAK2–STAT5 signaling" = c(
    "Prlr","Jak2","Stat5a","Stat5b","Cish","Socs2","Socs3","Elf5"
  ),
  
  "CCND1/cell-cycle" = c(
    "Ccnd1","Cdk4","Cdk6","Rb1","E2f1","E2f2","E2f3","Cdkn1b","Cdkn1a"
  ),
  
  "PI3K–AKT–mTOR signaling" = c(
    "Pik3ca","Pik3cb","Pik3r1","Pik3r2","Akt1","Pdpk1","Pten","Tsc1",
    "Tsc2","Mtor","Rheb","Rptor","Rictor","Rps6kb1","Eif4ebp1"
  ),
  
  "NF-κB survival" = c(
    "Bcl2","Bcl2l1","Mcl1","Birc5","Bcl2a1","Xiap","Traf1"
  )
)

# Create table
tab <- bind_rows(lapply(names(core_sets), function(gs) {
  tibble(
    `Gene set` = gs,
    Genes = paste(core_sets[[gs]], collapse = ", ")
  )
}))

# Write Excel
wb <- createWorkbook()
addWorksheet(wb, "Core gene sets")
writeData(wb, sheet = 1, x = tab)

# Basic formatting
header_style <- createStyle(
  textDecoration = "bold",
  halign = "center",
  border = "bottom"
)

addStyle(wb, 1, header_style, rows = 1, cols = 1:2, gridExpand = TRUE)
setColWidths(wb, 1, cols = 1, widths = 30)
setColWidths(wb, 1, cols = 2, widths = 120)
freezePane(wb, 1, firstRow = TRUE)

saveWorkbook(wb, out_file, overwrite = TRUE)
