FROM ghcr.io/web-mev/base-spatialge-docker:sha-b9c9d80863af47a4ac85f420b96fd1ff57e2cde7

USER root
RUN mkdir -p /opt/resources && \
  chmod a+r /opt/resources && \
  chown $MAMBA_USER:$MAMBA_USER /opt/resources

USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# in addition to the packages in the base image, we need
# the annotation packages to enable creation of gene identifier mapping files
RUN micromamba install -y \
    -c conda-forge \
    -c bioconda \
    bioconductor-org.hs.eg.db=3.18.0 \
    bioconductor-org.mm.eg.db=3.18.0

# Download msigdb database files for human and mouse:
# Human:
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/h.all.v2023.2.Hs.symbols.gmt > /opt/resources/h.all.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c2.cp.reactome.v2023.2.Hs.symbols.gmt > /opt/resources/c2.cp.reactome.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt > /opt/resources/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c5.go.bp.v2023.2.Hs.symbols.gmt > /opt/resources/c5.go.bp.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c5.go.mf.v2023.2.Hs.symbols.gmt > /opt/resources/c5.go.mf.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c5.go.cc.v2023.2.Hs.symbols.gmt > /opt/resources/c5.go.cc.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c6.all.v2023.2.Hs.symbols.gmt > /opt/resources/c6.all.v2023.2.Hs.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/human/gene_sets/c7.immunesigdb.v2023.2.Hs.symbols.gmt > /opt/resources/c7.immunesigdb.v2023.2.Hs.symbols.gmt

# Mouse
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/mh.all.v2023.2.Mm.symbols.gmt > /opt/resources/mh.all.v2023.2.Mm.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/m2.cp.reactome.v2023.2.Mm.symbols.gmt > /opt/resources/m2.cp.reactome.v2023.2.Mm.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/m2.cp.wikipathways.v2023.2.Mm.symbols.gmt > /opt/resources/m2.cp.wikipathways.v2023.2.Mm.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/m5.go.bp.v2023.2.Mm.symbols.gmt > /opt/resources/m5.go.bp.v2023.2.Mm.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/m5.go.mf.v2023.2.Mm.symbols.gmt > /opt/resources/m5.go.mf.v2023.2.Mm.symbols.gmt
RUN curl ftp://ftp.broadinstitute.org/pub/gsea/msigdb/mouse/gene_sets/m5.go.cc.v2023.2.Mm.symbols.gmt > /opt/resources/m5.go.cc.v2023.2.Mm.symbols.gmt

ADD create_gene_mappings.R /usr/local/bin

# Create the gene mapping files:
RUN Rscript /usr/local/bin/create_gene_mappings.R human /opt/resources/human_genes.tsv
RUN Rscript /usr/local/bin/create_gene_mappings.R mouse /opt/resources/mouse_genes.tsv

ADD stenrich.R /usr/local/bin
ADD run.sh /usr/local/bin