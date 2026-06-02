
library(tidyverse)
library(igraph)
library(tidygraph)
library(ggraph)

# Make names readable #-----
normalize_name <- function(x) {
  
  # dictionary of exceptions
  exceptions <- c(
    # ions / elements
    "na" = "Na",
    "cl" = "Cl",
    
    # receptors / signaling
    "gpcr" = "GPCR",
    "fgfr1" = "FGFR1",
    "fgfr1b" = "FGFR1B",
    "fgfr1c" = "FGFR1C",
    "fgfr2" = "FGFR2",
    "fgfr2b" = "FGFR2B",
    "fgfr2c" = "FGFR2C",
    "fgfr3" = "FGFR3",
    "fgfr3b" = "FGFR3B",
    "fgfrl1" = "FGFRL1",
    "foxo" = "FOXO",
    "notch3" = "NOTCH3",
    "tgfbr3" = "TGFBR3",
    "ntrk2" = "NTRK2",
    "erbb2" = "ERBB2",
    "erbb4" = "ERBB4",
    "frs2" = "FRS2",
    "frs3" = "FRS3",
    "fyn" = "FYN",
    "irf7" = "IRF7",
    "ptk"= "PTK",
    "traf6"= "TRAF6",
    "tlr7"= "TLR7",
    "fzd" = "FZD",
    "tcr" = "TCR",
    "alk" = "ALK",
    "rnd1" = "Rnd1",
    "rnd2" = "Rnd2",
    "rnd3" = "Rnd3",
    "stat5" = "STAT5",
    "ras" = "Ras",
    "akt1" = "AKT1",
    "e17k" = "E17K",
    "bcr" = "BCR",
    "mapk" = "MAPK",
    
    # proteins / genes
    "rac1" = "RAC1",
    "irs" = "IRS",
    "grb7" = "GRB7",
    "grb2" = "GRB2",
    "abl" = "ABL",
    "mitf" = "MITF",
    "neurog3" = "NEUROG3",
    "cdh11" = "CDH11",
    "shc" = "SHC",
    "frs" = "FRS",
    "pi3k" = "PI3K",
    "akt" = "AKT",
    "cd22" = "CD22",
    "cd28" = "CD28",
    "plk1" = "PLK1",
    "pten" = "PTEN",
    "mtor" = "mTOR",
    "bmal1" = "BMAL1",
    "arntl" = "ARNTL",
    "fbxl7" = "FBXL7", 
    "aurka" = "AURKA",
    "tak1" = "TAK1",
    "irak1" = "IRAK1",
    "hur elavl1" = "HuR ELAVL1",
    "runx2" = "RUNX2",
    
    # hormones / signaling molecules
    "igf" = "IGF",
    "igfbps" = "IGFBPs",
    "tnfs" = "TNFs",
    "igf1r" = "IGF1R",
    "ret" = "RET",
    
    # transcription factors
    "pou5f1" = "POU5F1",
    "oct4" = "OCT4",
    "sox2" = "SOX2",
    "nanog" = "NANOG",
    "ap 2" = "AP-2",
    "tfap2" = "TFAP2",
    "nfe2l2" = "NFE2L2",
    "mecp2" = "MECP2",
    
    # receptor classes
    "gaba" = "GABA",
    "trka" = "TRKA",
    
    # complexes
    "b wich" = "B-WICH",
    "ikk" = "IKK",
    
    # some more
    "gaps" = "GAPs",
    "cd163" = "CD163",
    "srebf" = "SREBF",
    "srebp" = "SREBP",
    "tnfr2" = "TNFR2",
    "nf kb" = "NF-κB",
    "pparalpha" = "PPARα",
    "mrna" = "mRNA",
    "axin" = "AXIN",
    "notch4" = "NOTCH4",
    "dvl" = "DVL",
    "uch" = "UCH",
    "sema4d" = "SEMA4D",
    "tysnd1" = "TYSND1",
    "snrna" = "snRNA",
    "fgfr1" = "FGFR1",
    "fgfr4" = "FGFR4",
    "smad2" = "SMAD2",
    "smad3" = "SMAD3",
    "smad4" = "SMAD4",
    "robo" = "ROBO",
    "slit" = "SLIT",
    "spry" = "SPRY",
    "fgf" = "FGF",
    "hsf1" = "HSF1",
    "odc" = "ODC",
    "met" = "MET",
    "ptpn11" = "PTPN11",
    "b5" = "B5",
    "sars cov 1" = "SARS-CoV-1",
    "slc2a4" = "SLC2A4",
    "glut4" = "GLUT4",
    "mapk6" = "MAPK6",
    "mapk4" = "MAPK4",
    "flt3" = "FLT3",
    "abc" = "ABC",
    "netrin 1" = "Netrin-1",
    "runx3" = "RUNX3",
    "yap1" = "YAP1",
    "ecd" = "ECD",
    "egfrviii" = "EGFRvIII",
    "adme" = "ADME",
    "hdr" = "HDR",
    "hrr" = "HRR",
    "crmps" = "CRMPs",
    "sema3a" = "Sema3A",

    # general abbreviations
    "dna" = "DNA",
    "rna" = "RNA",
    "atp" = "ATP",
    "trna" = "tRNA",
    "mrna" = "mRNA",
    "rrna" = "rRNA",
    "tca" = "TCA",
    "ecm" = "ECM",
    "golgi" = "Golgi",
    "gtpase" = "GTPase",
    "orc1" = "ORC1",
    "gpi" = "GPI",
    "hs gag" = "HS-GAG",
    "hh np" = "Hh-Np",
    "micrornas" = "microRNAs",
    "pi" = "PI",
    "RNA polymerase ii" = "RNA polymerase II",
    "snrna" = "snRNA",
    "24 hydroxycholesterol" = "24-hydroxycholesterol",
    "hiv" = "HIV",
    "anti inflammatory" = "anti-inflammatory",
    "transcription factors" = "TFs",
    "o linked" = "o-linked",
    "cell cell" = "cell-cell",
    "phase ii" = "phase II",
    "hdacs" = "HDACs",
    "regulation of gene expression" = "gene expression regulation",
    "rhov" = "RHO", # RHOV does not make sense??
    "type i" = "type I",
    "cd209 dc sign" = "CD209 DC-SIGN",
    "ip3" = "IP3",
    "ip4" = "IP4",
    "chst3" = "CHST3",
    "sedcjd" = "SEDCJD",
    "p14 arf" = "p14-ARF",
    "respiratory electron transport" = "RET",
    "gli1" = "GLI1",
    "c type lectin receptors clrs" = "C-type lectin receptors",
    "ub specific" = "Ub-specific",
    "fc epsilon receptor fceri" = "FCERI",
    "hdl" = "HDL",
    
    # Full labels
    "negative regulation of tcf dependent signaling by wnt ligand antagonists" = "negative regulation of TCF signaling by Wnt antagonists",
    "interferon alpha beta signaling" = "interferon alpha/beta signaling",
    "keap1 NFE2L2" = "KEAP1-NFE2L2",
    "raf independent mapk1 3 activation" = "RAF-independent MAPK1/3 activation",
    "tp53 regulates transcription of several additional cell death genes whose specific roles in p53 dependent apoptosis remain uncertain" = "TP53-regulated apoptosis genes (uncertain role)",
    "regulation of tp53 activity through association with co factors" = "cofactor regulation of TP53",
    "defective c1galt1c1 causes tnps" = "defective C1GALT1C1 causes TNPS",
    "regulation of MITF m dependent genes involved in extracellular matrix focal adhesion and epithelial to mesenchymal transition"= "MITF-M regulation of ECM and EMT",
    "signaling by type 1 insulin like growth factor 1 receptor IGF1R" = "IGF1R signaling",
    "acetylcholine inhibits contraction of outer hair cells" = "ACh inhibits outer hair cell contraction",
    "regulation of MITF m dependent genes involved in apoptosis" = "MITF-M regulation of apoptosis",
    "regulation of mRNA stability by proteins that bind au rich elements" = "ARE-binding proteins and mRNA stability",
    "synthesis of bile acids and bile salts via 24-hydroxycholesterol" = "bile acid synthesis via 24-hydroxycholesterol",
    "defective chsy1 causes tpbs" = "defective CHSY1 causes TPBS",
    "gene expression regulation in endocrine committed NEUROG3 progenitor cells" = "gene expression regulation in NEUROG3 endocrine progenitors",
    "regulation of insulin like growth factor IGF transport and uptake by insulin like growth factor binding proteins IGFBPs" = "regulation of IGF transport and uptake by IGFBPs"
  )
  
  x <- x %>%
    tolower() %>%
    gsub("reactome_", "",.) %>%
    gsub("_", " ", .) %>%
    gsub("\\s+", " ", .) %>%
    trimws()
  
  # replace exception words
  for (w in names(exceptions)) {
    x <- gsub(paste0("\\b", w, "\\b"), exceptions[w], x)
  }
  
  # capitalize first letter only
  x <- sub("^([a-z])", "\\U\\1", x, perl = TRUE)
  
  x
  
}

# Helper: define prioritized themes and their regex patterns #-----
default_theme_patterns <- tribble(
  ~theme,                                   ~pattern,
  
  # High-priority axis for MPA context
  "Steroid hormone metabolism / RANKL axis", "(progesterone|progestin|androgen|estrogen|steroid(ogenesis)?|pregnenolone biosynthesis|CYP1[179]|CYP19|CYP21|HSD(3|17)B|RANKL|TNFSF11|\\bRANK\\b|TNFRSF11A|OPG|TNFRSF11B)",
  "FGFR / RTK signaling",                   "(\\bFGFR\\b|FGF|IGF1R|ERBB|EGFR|MET\\b|FLT3|RET|TRKA|NTRK1|IGFBPS|NTRK2)",
  "PI3K / AKT / PTEN",                      "(\\bPI3K\\b|AKT\\b|PTEN\\b|PIP3|SLC2A4|GLUT4)",
  "MAPK signaling",                         "(\\bMAPK\\b|RAF\\b|TAK1|ERK|MAPK1/3|MAPK6|MAPK4|Ras by GAPs)",
  "TGF-beta / SMAD signaling",              "(\\bTGF\\b|SMAD\\b|TGFBR\\b|TGFBR3|SMAD2|SMAD3|SMAD4)",
  "Notch / Hedgehog / Wnt",                 "(NOTCH|HEDGEHOG|SHH|PTCH|SMO|WNT|DVL|AXIN|\\bbeta\\s*catenin)",
  "Developmental / lineage programs",       "(gastrulation|somitogenesis|mesoderm|cardiogenesis|kidney development|neuronal system|developmental cell lineages|neural plate|germ cells|lineage|formation of .* mesoderm|sensory perception|beta cell|pancreatic precursor|testis differentiation|adipocyte|adipogenesis|white adipocyte|PPAR\\s?γ|PPARG)",
  "Transcription factors / chromatin",      "(POU5F1|OCT4|SOX2|NANOG|MITF|RUNX2|TFAP2|MECP2|HSF1|FOXO|RUNX3|YAP1|HDACs|NEUROG3)",
  "Apoptosis / TP53 / senescence",          "(TP53|p53|apoptosis|senescence|proapoptotic|hypoxia|\\bHIF\\b)",
  "Immune / inflammation",                  "(interferon|interleukin|TNFR2|TCR\\b|BCR\\b|antigen processing|inflammasome|HIV|SARS|KEAP1\\-NFE2L2|NRF2|CD163|innate immune system|antimicrobial peptides|IRAK1|\\bIKK\\b|\\bTLR\\b|\\bSIRP\\b|chemical stress|TNFS)",
  "GPCR signaling",                         "(\\bGPCR\\b|chemokine|ACh\\b|muscarinic|adrenergic|incretin)",
  "Proteostasis / UPS / autophagy",         "(proteasome|ubiquitin|ubiquitination|neddylation|sumoylation|\\bUCH\\b|autophagy|chaperone|heat stress|protein repair)",
  "RNA processing / translation",           "(mRNA editing|deadenylation|ARE\\-|HuR|ELAVL1|RNA polymerase|snRNA|mRNA stability|translation|metabolism of RNA)",
  "DNA replication / repair / cell cycle",  "(DNA repair|DNA replication|cell cycle)",
  "Metabolism (lipid / AA / vitamin)",      "(fatty acid|beta oxidation|bile acid|phenylalanine|histidine|lysine catabolism|ornithine decarboxylase|pantothenate|vitamin|glucocorticoid|sulfide oxidation|glycerophospholipid|cofactor|SREBF|SREBP|PPAR|synthesis of PI|phase ii conjugation|glycosphingolipid|sialic acid|polyamines|chylomicron|creatine)",
  "Mitochondria / mitophagy",               "(mitochondria|mitophagy|respiration|electron transport|peroxisom)",
  "ECM / adhesion / GAG",                   "(extracellular matrix|ECM|elastic fibre|glycosaminoglycan|HS\\-GAG|heparan sulfate|chondroitin|dermatan|o[- ]?linked glycosylation|cell\\-cell junction|junction organization|adhesion|mucopolysaccharidoses|diseases associated with glycosaminoglycan|TNPS|CHSY1)",
  "Cytoskeleton / guidance / Rho",          "(rho gtpase|rho\\b|eph\\b|ephrin|sema\\b|semaphorin|slit\\b|robo\\b|reelin|axon guidance|cell repulsion|muscle contraction|myosin|actin cytoskeleton|netrin\\-1|RAC1|Sema3A)",
  "Transport / ion channels",               "(ion channel|sodium calcium exchanger|sodium channel|\\bna channel\\b|gap junction|monocarboxylate|ABC family|transmembrane transporter|\\bSLC\\b|transport of fatty acids|stimuli sensing|transport of bile salts)",
  "Circadian / clock",                      "(circadian|BMAL1|ARNTL|clock)",
  "Hormone / cytokine receptors",           "(growth hormone receptor|prolactin receptor|erythropoietin|insulin receptor|IGF1R|physiological factors)",
  "Disease / host interaction",             "(defects in vitamin)"
)

# Manual overrides 
manual_overrides <- c(
  "neuronal system" = "Transport / ion channels",
  "keap1-nfe2l2 pathway" = "Immune / inflammation",
  "uch proteinases" = "Proteostasis / UPS / autophagy",
  "phase 1 inactivation of fast na channels" = "Transport / ion channels",
  "phase ii conjugation of compounds" = "Metabolism (lipid / AA / vitamin)",
  "signal regulatory protein family interactions" = "Immune / inflammation"
)

# Core function: add primary theme and all matching themes
add_pathway_theme <- function(df, name_col = "gs_name_norm",
                              theme_patterns = default_theme_patterns,
                              overrides = manual_overrides) {
  
  stopifnot(name_col %in% names(df))
  # Prepare lowercase name for matching
  nm <- rlang::sym(name_col)
  
  df2 <- df %>%
    mutate(
      .name_lc = str_to_lower(!!nm)
    )
  
  # Apply manual overrides first
  df2 <- df2 %>%
    mutate(
      theme_manual = recode(.name_lc, !!!overrides, .default = NA_character_)
    )
  
  # Find all matching themes per row
  matches <- theme_patterns %>%
    mutate(pattern_i = row_number()) %>%
    pmap(function(theme, pattern, pattern_i) {
      # Return a logical vector of matches for this pattern
      str_detect(df2$.name_lc, regex(pattern, ignore_case = TRUE))
    }) %>%
    set_names(theme_patterns$theme) %>%
    as_tibble()
  
  # Collect all themes that matched (could be none or several)
  df3 <- bind_cols(df2, matches) %>%
    pivot_longer(cols = all_of(theme_patterns$theme),
                 names_to = "theme_candidate", values_to = "matched") %>%
    group_by(across(all_of(names(df)))) %>%
    summarize(
      theme_all = paste0(theme_candidate[matched], collapse = "; "),
      .groups = "drop"
    ) %>%
    mutate(theme_all = ifelse(theme_all == "", NA_character_, theme_all))
  
  # Prioritized primary theme:
  # - use manual override if present
  # - otherwise first theme in the ordered list that matched
  df_final <- df2 %>%
    left_join(df3, by = names(df)) %>%
    mutate(
      theme_primary = coalesce(
        theme_manual,
        map_chr(.name_lc, function(x) {
          idx <- which(map_lgl(seq_len(nrow(theme_patterns)), function(i) {
            str_detect(x, regex(theme_patterns$pattern[i], ignore_case = TRUE))
          }))
          if (length(idx) == 0) NA_character_ else theme_patterns$theme[min(idx)]
        })
      )
    ) %>%
    dplyr::select(-.name_lc, -theme_manual)
  
  df_final
}

# Build pseudo-DAG for pathway pruning #-----
build_pseudo_dag <- function(gene_sets){
  
  # Create a list: pathway -> genes
  pathway_genes <- gene_sets %>%
    dplyr::select(gs_name, gene_primary) %>%
    group_by(gs_name) %>%
    summarise(genes = list(unique(gene_primary))) %>%
    deframe()
  
  # Define edges by subset relationships
  edges <- data.frame(from = character(), to = character(), stringsAsFactors = FALSE)
  
  paths <- names(pathway_genes)
  n <- length(paths)
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) next
      # genes i and j
      genes_i <- pathway_genes[[paths[i]]]
      genes_j <- pathway_genes[[paths[j]]]
      
      # if i is subset of j, create edge i -> j
      if (all(genes_i %in% genes_j)) {
        edges <- rbind(edges, data.frame(from = paths[i], to = paths[j], stringsAsFactors = FALSE))
      }
    }
  }
  
  # build graph
  g_pseudo <- igraph::graph_from_data_frame(edges, directed = TRUE)
  
  # check if DAG
  if (!igraph::is_dag(g_pseudo)) {
    stop("DAG check failed: g_pseudo contains cycles. Cannot proceed.")
  } 
  
  message("DAG successfully built.")
  
  return(g_pseudo)
}

# Pathway pruning #----
prune_leaves <- function(ids, graph) {
  keep <- ids
  for (i in ids) {
    if (!(i %in% V(graph)$name)) next
    # get descendant vertices by name
    des <- igraph::subcomponent(graph, i, mode="out")
    des_names <- V(graph)$name[des]  # convert numeric vertex IDs to names
    # remove i if any descendant is also significant
    if (any(des_names %in% ids & des_names != i)) keep <- setdiff(keep, i)
  }
  keep
}

prune_specific <- function(ids, graph) {
  keep <- ids
  for (i in ids) {
    if (!(i %in% V(graph)$name)) next
    # get ancestor vertices by name
    anc <- igraph::subcomponent(graph, i, mode = "in")
    anc_names <- V(graph)$name[anc]
    # remove i if any ancestor is also significant
    if (any(anc_names %in% ids & anc_names != i)) keep <- setdiff(keep, i)
  }
  keep
}

prune_paths_from_dag <- function(g_dag, sig_paths, approach = "general"){
  
  missing_paths <- setdiff(sig_paths, V(g_dag)$name)
  message(paste0(length(missing_paths), " not in provided DAG.\n", missing_paths))
  
  sig_paths_mapped <- ifelse(
    sig_paths %in% V(g_pseudo)$name,
    sig_paths,
    NA
  )
  sig_paths_mapped <- na.omit(sig_paths_mapped)
  
  if(approach == "general"){
    pruned <- prune_leaves(sig_paths_mapped, g_dag)
  } else if(approach == "specific"){
    pruned <- prune_specific(sig_paths_mapped, g_dag)
  }
  
  pruned_all <- c(pruned, missing_paths)
  
  return(pruned_all)
}

# Vizualize DAG/subDAG #-------
plot_dag <- function(g_dag, paths_sub = NULL, layout = "kk"){
  
  # Subset graph if subset provided
  if(!is.null(paths_sub)){
    subg <- induced_subgraph(g_dag, vids = intersect(V(g_dag)$name, paths_sub))   
  } else{
    subg <- g_dag
  }
  
  # Make plot
  graph_tbl <- as_tbl_graph(subg)
  lay <- create_layout(graph_tbl, layout = layout)
  p <- ggraph(lay) +
    geom_edge_link(
      alpha = 0.4,
      width = 0.4,
      arrow = arrow(length = unit(3, "mm")),
      end_cap = circle(2, "mm"),
      color = "darkgrey"
    ) +
    geom_node_point(size = 2, color = "darkgrey") +
    geom_node_text(
      aes(label = normalize_name(V(subg)$name)),
      size = 3,
      repel = TRUE,
      max.overlaps = Inf
    ) +
    theme_void()
  
  return(p)
}