This repository is provided for reproducibility of the accompanying manuscript: <br>
Vavourakis CD*, Weber DD*, Aminzadeh-Gohari S*, Tevini J, Barret JE, Herzog C, Catalano L, Stefan VE, Redl E, Alkasalias T, Felder TK, Lang R, Kofler B, Widschwendter M. Ketogenic diet and mifepristone suppress mammary tumorigenesis and reverse progestin-driven epigenetic field defects.  <br>
<br><br>

**Please note**: <br>
* EpiMod functions from West et al (2013, *Scientific reports*, https://doi.org/10.1038/srep01630) were revised to work with igraph package version 2.2.2. Additionally, functionallity was added to remove highly dense, non-informative network regions, and to output ggplot objects. Licensing status of the original implementation is unknown. Users should cite the original work. <br> 
* TCA analysis was performed as described by Rahmani et al (2019, *Nature Communications*, https://doi.org/10.1038/s41467-019-11052-9). <br>
<br><br>

**Instructions**:

Data needed:
* Raw methylation data (Mouse methylation array) is available through the NCBI GEO repository GSE2362360. <br>
* Beta matrix creation, CpG annotation and cell deconvolution was done as described by Barrett et al (2025, *Communications Medicine*, https://doi.org/10.1038/s43856-025-00779-w). <br>
* Precalculated DNA methylation scores, cell type estimates and mouse phenotypic data needed to further reproduce the analysis results are available here 10.5281/zenodo.20513216 (will be released upon manuscript acceptance). <br>
<br><br>

To start please, download the raw data and create a file called beta_final.Rdata as described by Barrett et al and store it in /data.<br>
For the final beta matrix: Rows = CpGs, Columns = samples. <br>
<br><br>

In /data also store: <br>
* mouse phenotypic data (from Zenodo) <br>
* PPI from Biogrid (extract from https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-5.0.255/BIOGRID-ORGANISM-5.0.255.mitab.zip) <br>
* PPI from IntAct (intact_2.7_mouse.txt; miTab 2.7 from https://www.ebi.ac.uk/intact/interactomes) <br>

In /src download: <br>
* vignette_analysis_tca.R from https://github.com/cozygene/TCA