suppressMessages(library(spatialGE))

# Remove after importing flat file for the mSigDB data
library(msigdbr)

# args from command line:
args <- commandArgs(TRUE)

# note, -k --kclusters can be integer or character ('dtc') for valid options
option_list <- list(
    make_option(
        c('-f', '--input_file'),
        help='Path to the count matrix input.'
    ),
    make_option(
        c('-c', '--coordinates_file'),
        help='Path to the barcode spatial coordinates input.'
    ),
    make_option(
        c('-s', '--sample_name'),
        help='Sample name'
    ),
    make_option(
        c('-o','--output_file_prefix'),
        help='The prefix for the output file'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))


# load the data
spat <- STlist(rnacounts="./s117d", samples=c("s117d"))
# normalize
spat <- transform_data(spat, method='sct')

# Load the gene sets for human from mSigDB
# Returns a tibble (4,331,807 x 15)
#   gs_cat gs_subcat      gs_name        gene_symbol entrez_gene ensembl_gene   
#   <chr>  <chr>          <chr>          <chr>             <int> <chr>          
# 1 C3     MIR:MIR_Legacy AAACCAC_MIR140 ABCC4             10257 ENSG00000125257
# 2 C3     MIR:MIR_Legacy AAACCAC_MIR140 ABRAXAS2          23172 ENSG00000165660
# 3 C3     MIR:MIR_Legacy AAACCAC_MIR140 ACTN4                81 ENSG00000130402
geneSets <- msigdbr(species="Homo sapiens")
# Filter the gene sets for Hallmark genes
# Just an example to shrink the test space
# Optional
geneSets <- geneSets[geneSets['gs_cat'] == "H",]

# Convert gene set dataframe to list
# The result is a named list. The names of the list are the names of each gene set
# The contents of each list element are the gene names within each gene set
geneSets <- split(x=geneSets[['gene_symbol']], f=geneSets[['gs_name']])
# Example: geneSets$HALLMARK_WNT_BETA_CATENIN_SIGNALING
# [1] "ADAM17" "AXIN1"  "AXIN2"  "CCND2"  "CSNK1E" "CTNNB1" "CUL1"   "DKK1"  
# [9] "DKK4"   "DLL1"   "DLL1"   "DVL2"   "FRAT1"  "FZD1"   "FZD8"   "GNAI1" 
#[17] "HDAC11" "HDAC2"  "HDAC5"  "HEY1"   "HEY2"   "JAG1"   "JAG2"   "KAT2A" 
#[25] "LEF1"   "MAML1"  "MAML1"  "MYC"    "NCOR2"  "NCSTN"  "NKD1"   "NOTCH1"
#[33] "NOTCH4" "NOTCH4" "NOTCH4" "NOTCH4" "NOTCH4" "NOTCH4" "NOTCH4" "NUMB"  
#[41] "PPARD"  "PSEN2"  "PTCH1"  "RBPJ"   "SKP2"   "TCF7"   "TP53"   "WNT1"  
#[49] "WNT5B"  "WNT6"  


# Perform enrichment analysis (exports as tibble)
# Options to potentially make transparent
# Below are default values
# reps = number of random samples to be extracted
# num_sds = number of standard deviations to set min gene set expression threshold
# min_genes = min number of genes in a gene set in the data set to include that gene set
# min_units = min number of spots with high expression of apathway for that gene set to
#             be considered in the analysis
stenrich_df <- STenrich(
    spat, 
    gene_sets = geneSets, 
    reps = 1000, 
    num_sds = 1, 
    min_genes = 5,
    min_units = 20, 
    pval_adj_method = "BH",
    seed = 12345
)

# Need to pull specific sample out of tibble
# Dropped first column that only outlined the sample name
# If we want to put the gene in the gene set in the tibble, we
# will have to join with the geneSets object.
write.table(
    stenrich_df$s117d[, 2:dim(stenrich_df$s117d)[2]],
    opt$output_file_prefix,
    sep="\t",
    row.names=F,
    quote=F
)