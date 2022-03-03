### GO data and processing functions

# WORMGO
#```{r GO Functions, echo=FALSE, cache=TRUE}
# get the annotations from PARASITE

# if (! "paramart" %in% ls()) {
#   system.time({paramart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)})}
# if (! "WORMGO" %in% ls()) {
#   # go to https://parasite.wormbase.org/biomart/martview/ to figure out the values to use for this
#   system.time({WORMGO = biomaRt::getBM(
#     mart = paramart,
#     filter = c("species_id_1010","biotype"),
#     value = list(species_id_1010="caelegprjna13758",biotype="protein_coding"),
#     attributes = c(
#       "wbps_gene_id",
#       "external_gene_id",
#       "go_accession",
#       "go_name_1006",
#       "go_linkage_type"
#     )
#   )}) 
# }

C_elegans_query = function(paramart) {

  WORMGO = biomaRt::getBM(
        mart = paramart,
        filter = c("species_id_1010","biotype"),
        value = list(species_id_1010="caelegprjna13758",biotype="protein_coding"),
        attributes = c(
          "wbps_gene_id",
          "external_gene_id",
          "go_accession",
          "go_name_1006",
          "go_linkage_type"
        )
      )
  WORMGO
}

# create an object where you can access all the GO terms that are assigned to a specific gene

geneID2GO = function(WORMGO) {
  by(WORMGO$go_accession, 
     WORMGO$wbps_gene_id, 
     function(x) as.character(x))
}

# all.genes <- unique(as.character(WORMGO$wbps_gene_id))
# unique.ap.wbid = unique(annotatedPeaks$ap$feature)

# set_genelist=list(fore=unique.ap.wbid, back=all.genes)
# dataset = mkGO(set_genelist$fore, set_genelist$back)
# BP.tab = GOSummary(dataset$BP) %>% filter(fisher < .05) %>% arrange(fisher) %>% mutate(DB="BP")
# CC.tab = GOSummary(dataset$CC) %>% filter(fisher < .05) %>% arrange(fisher) %>% mutate(DB="CC")
# MF.tab = GOSummary(dataset$MF) %>% filter(fisher < .05) %>% arrange(fisher) %>% mutate(DB="MF")

GOSummary<- function(GOdata, topNodes = 200) {
  library(topGO)
  resultClassic <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultElim <- topGO::runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  topGO::GenTable(
    object=GOdata, 
    classicFisher = resultClassic,
    elim=resultElim,    
    orderBy="classicFisher",
    topNodes = topNodes
  ) -> tab
  # not sure where the conversion to char is happening. convert back
  tab$fisher = as.numeric(tab$classicFisher)
  tab$elim = as.numeric(tab$elim)
  return(tab)
}

runGO = function(foreground_genes, background_genes, WORMGO, topNodes = 200)
{
  go = mkGO(foreground_genes, background_genes, WORMGO)
  go$BP.result = GOSummary(go$BP, topNodes)
  go$MF.result = GOSummary(go$MF, topNodes)
  go$CC.result = GOSummary(go$CC, topNodes)
  go
}

mkGO = function(foreground_genes, background_genes, WORMGO) {
  library(topGO)
  # create a TRUE/FALSE list the size of background genes, TRUE if in foreground
  geneList = factor(as.integer(background_genes %in% foreground_genes))
  names(geneList) = background_genes
  
  geneID2GO = geneID2GO(WORMGO)
  
  BP.go = new("topGOdata", ontology='BP'
              , allGenes = geneList
              , annot = topGO::annFUN.gene2GO
              , gene2GO = geneID2GO)
  MF.go = new("topGOdata", ontology='MF'
              , allGenes = geneList
              , annot = topGO::annFUN.gene2GO
              , gene2GO = geneID2GO)
  CC.go = new("topGOdata", ontology='CC'
              , allGenes = geneList
              , annot = topGO::annFUN.gene2GO
              , gene2GO = geneID2GO)
  list(BP=BP.go,CC=CC.go,MF=MF.go)
}