# Make combined data table with Rob and my analyses
# 

library(dplyr)
library(GenomicRanges)
library(magrittr)
library(ParasiteXML)
source('david-reader.R')

elt2.data = read_ELT2_binding_data(as_genomic_ranges = FALSE) # from david-reader.R
rob = read_rob_all_merged() %>% dplyr::select(-starts_with("pvalue."),-starts_with("lfcSE."))

merge = right_join(rob, elt2.data, by = "WBGeneID")

# fix a name 
merge %<>% dplyr::rename(L3.log_chip_signal_mean = L3_log.chip_signal_mean, chrom=seqnames)

# select columns
alldata = merge %>% dplyr::select(
                chrom,
                start,
                end,
                strand,
                wikigene_name,
                WBGeneID,
                ends_with("log_chip_signal_mean"), 
                starts_with("rlogc."),
                starts_with("log2FoldChange."),
                ends_with("_int_exp"),
                ends_with("_bound"),
                din.status.description)

# there are some gene IDs missing from Rob's wiki_gene_name,
# use external_gene_id from biomart
library(biomaRt)
library(ParasiteXML)

bm = getParamart()
QUERY='<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "parasite_mart" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
	<Dataset name = "wbps_gene" interface = "default" >
		<Filter name = "biotype" value = "protein_coding"/>
		<Filter name = "species_id_1010" value = "caelegprjna13758"/>
		<Attribute name = "wbps_gene_id" />
		<Attribute name = "external_gene_id" />
	</Dataset>
</Query>'

library(dplyr)
R_query = format_BM_from_XML(QUERY)
mart = getParamart()
genes = runWithMart(R_query, mart)
genes %<>% dplyr::rename(WBGeneID=wbps_gene_id)
rownames(genes) <- genes$WBGeneID
alldata$name = genes[alldata$WBGeneID, 'external_gene_id']
alldata %<>% relocate(name, .after = WBGeneID)
alldata$wikigene_name = NULL

# sort by chromosome, start point
alldata %<>% arrange(chrom,start) 
write_tsv(alldata, "../03_output/promoter_data/all.bedplus")
write_rds(alldata, "../03_output/promoter_data/all.df.rds", compress="xz")
