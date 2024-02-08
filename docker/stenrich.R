suppressMessages(library(spatialGE))
suppressMessages(library(reshape2))

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
        c('-g', '--gmt_file'),
        help='GMT formatted gene set input'
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

# load the GMT
gmt <- read.table(
    opt$gmt_file,
    header=F,
    sep="\t",
    fill=T
)
# Munge the GMT into the format required by spatialGE
# Clunky but works
geneSets <- t(gmt[, 3:dim(gmt)[2]])
colnames(geneSets) <- gmt[, 1]
geneSets <- melt(geneSets)
geneSets <- geneSets[geneSets$value != "", ] # Drop the empty gene fields
# If you are not sending individual pathways files for the
# various collections, a filter will be needed.

# Convert dataframe to named list
geneSets <- split(x=geneSets[['value']], f=geneSets[['Var2']])

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