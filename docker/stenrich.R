suppressMessages(library(spatialGE))
suppressMessages(library(reshape2))
suppressMessages(library(rjson))
suppressMessages(library(optparse))
suppressMessages(library(dplyr))
source('/usr/local/bin/prep_stlist.R')

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
        c('-n', '--normalization'),
        help='Normalization method of `log` or `sct`'
    ),
    make_option(
        c('-g', '--gmt_file'),
        help='GMT formatted gene set input'
    ),
    make_option(
        c('-m', '--map_file'),
        help='TSV file for mapping gene identifiers'
    ),
    make_option(
        c('-i', '--gene_ids'),
        help='The gene identifier system used.'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Sanity checks on the inputs:
# Check that the file was provided:
if (is.null(opt$input_file)){
    message('Need to provide a count matrix with the -f/--input_file arg.')
    quit(status=1)
}

if (is.null(opt$coordinates_file)){
    message('Need to provide a count matrix with the -c/--coordinates_file arg.')
    quit(status=1)
}

if (is.null(opt$map_file)){
    message('Need to provide a gene ID mapping file with the -m/--map_file arg.')
    quit(status=1)
}

if (is.null(opt$gene_ids)){
    message('Need to provide the gene identifier with the -i/--gene_ids arg.')
    quit(status=1)
} else {
    gene_ids <- toupper(opt$gene_ids)
}

if (is.null(opt$gmt_file)){
    message('Need to provide a GMT file with the -g/--gmt_file arg.')
    quit(status=1)
}

# transform the name of the normalization scheme:
if (is.null(opt$normalization)){
    message('Need to provide a normalization scheme with the -n/--normalization arg.')
    quit(status=1)
} else if(tolower(opt$normalization) == 'sctransform'){
    norm_scheme <- 'sct'
} else if(tolower(opt$normalization) == 'log'){
    norm_scheme <- 'log'
} else {
    message('We only accept `log` or `SCTransform` for the normalization scheme.')
    quit(status=1)
}

# change the working directory to co-locate with the counts file:
working_dir <- dirname(opt$input_file)
setwd(working_dir)

# dataframe which maps from one system to another:
gene_mapping_df <- read.table(opt$map_file, sep='\t', header=T)

# prepare an STList instance. Note that we are potentially re-mapping
# the original gene identifiers to symbols such that they work with
# the MSigDB files:
spat_list <- prep_stlist(opt$input_file, 
                         opt$coordinates_file,
                         opt$sample_name,
                         gene_mapping_df,
                         gene_ids,
                         'SYMBOL')
spat <- spat_list$spat

# normalize
spat <- transform_data(spat, method=opt$normalization)

# load the GMT. Note that fill=T handles the fact that we have a 
# "jagged" table since each gene set can have a different number of genes 
# (and hence each line has a different number of columns/entries)
gmt <- read.table(
    opt$gmt_file,
    header=F,
    sep="\t",
    fill=T
)
# Munge the GMT into the format required by spatialGE
# Clunky but works
# The first two cols have the name of the gene set and a URL. The rest of the
# row contains the gene names. The next two lines create sets of genes in a 
# column-wise fashion. For smaller gene sets, the "bottom" of the column
# is populated by emptry strings
geneSets <- t(gmt[, 3:dim(gmt)[2]])
colnames(geneSets) <- gmt[, 1]
geneSets <- melt(geneSets)
geneSets <- geneSets[geneSets$value != "", ] # Drop the empty gene fields
# geneSets looks like:
# > head(geneSets)
#   Var1                  Var2 value
# 1   V3 HALLMARK_ADIPOGENESIS ABCA1
# 2   V4 HALLMARK_ADIPOGENESIS ABCB8
# 3   V5 HALLMARK_ADIPOGENESIS ACAA2
# 4   V6 HALLMARK_ADIPOGENESIS ACADL
# 5   V7 HALLMARK_ADIPOGENESIS ACADM
# 6   V8 HALLMARK_ADIPOGENESIS ACADS


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
sample_df <- stenrich_df[[opt$sample_name]]
print(head(sample_df))
sample_df <- sample_df[, 2:dim(sample_df)[2]]

# map the geneSet object from symbols (Default) back to 
# the identifier system used in the original counts.
remapped_geneSets <- lapply(geneSets, function(pathway_genes){
    m <- merge(gene_mapping_df[,c('SYMBOL', gene_ids)],
                   pathway_genes, by.x='SYMBOL', by.y=1) %>% distinct()
    unique(m[, gene_ids])
})
q = apply(
    sample_df,
    1,
    function(r){
        list(
            pathway=r[['gene_set']], 
            pval=as.numeric(r[['p_value']]), 
            padj=as.numeric(r[['adj_p_value']]), 
            size=as.numeric(r[['size_gene_set']]),
            genes=remapped_geneSets[[r['gene_set']]]
        )
    }
)

results_json_str <- toJSON(q)
output_filename = 'stenrich_results.json'
results_json_file <- paste(working_dir, output_filename, sep='/')
write(results_json_str, results_json_file)

# for WebMEV compatability, need to create an outputs.json file.
json_str = paste0(
       '{"stenrich_results":"', results_json_file, '"}'
)
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)