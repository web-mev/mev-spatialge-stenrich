#!/bin/bash

RAW_COUNTS=$1
COORDS=$2
SAMPLE_ID=$3
NORMALIZATION_METHOD=$4
GENE_SETS_DB=$5
ORGANISM=$6
GENE_IDS=$7
XPOS_COL=$8
YPOS_COL=$9
OUTPUT_FILENAME=$10

declare -A HUMAN_SETS_MAP
HUMAN_SETS_MAP["MSigDB Hallmark"]="/opt/resources/h.all.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["Reactome"]="/opt/resources/c2.cp.reactome.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["WikiPathways"]="/opt/resources/c2.cp.wikipathways.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["Ontology - Biological process"]="/opt/resources/c5.go.bp.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["Ontology - Molecular function"]="/opt/resources/c5.go.mf.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["Ontology - Cellular component"]="/opt/resources/c5.go.cc.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["(Human only) Oncogenic signatures"]="/opt/resources/c6.all.v2023.2.Hs.symbols.gmt"
HUMAN_SETS_MAP["(Human only) ImmuneSigDB"]="/opt/resources/c7.immunesigdb.v2023.2.Hs.symbols.gmt"

declare -A MOUSE_SETS_MAP
MOUSE_SETS_MAP["MSigDB Hallmark"]="/opt/resources/mh.all.v2023.2.Mm.symbols.gmt"
MOUSE_SETS_MAP["Reactome"]="/opt/resources/m2.cp.reactome.v2023.2.Mm.symbols.gmt"
MOUSE_SETS_MAP["WikiPathways"]="/opt/resources/m2.cp.wikipathways.v2023.2.Mm.symbols.gmt"
MOUSE_SETS_MAP["Ontology - Biological process"]="/opt/resources/m5.go.bp.v2023.2.Mm.symbols.gmt"
MOUSE_SETS_MAP["Ontology - Molecular function"]="/opt/resources/m5.go.mf.v2023.2.Mm.symbols.gmt"
MOUSE_SETS_MAP["Ontology - Cellular component"]="/opt/resources/m5.go.cc.v2023.2.Mm.symbols.gmt"

if [[ $ORGANISM == "Human" ]]
then
    GENE_SET_FILE=${HUMAN_SETS_MAP[${GENE_SETS_DB}]}
    GENE_MAP_FILE=/opt/resources/human_genes.tsv
elif [[ $ORGANISM == "Mouse" ]]
then
    GENE_SET_FILE=${MOUSE_SETS_MAP[${GENE_SETS_DB}]}
    GENE_MAP_FILE=/opt/resources/mouse_genes.tsv
else
    echo "Not a valid organism choice: $ORGANISM" >&2
    exit 1;
fi

if [ ! -z "$GENE_SET_FILE" ]
then
    Rscript /usr/local/bin/stenrich.R \
        -f $RAW_COUNTS \
        -c $COORDS \
        -s $SAMPLE_ID \
        -n $NORMALIZATION_METHOD \
        -g $GENE_SET_FILE \
        -m $GENE_MAP_FILE \
        -i $GENE_IDS \
        -o $OUTPUT_FILENAME \
        -x $XPOS_COL \
        -y $YPOS_COL
else
    echo "The choice of \"$GENE_SETS_DB\" was not a valid choice for your organism ($ORGANISM). Note that some gene sets are organism specific."
    exit 1;
fi
