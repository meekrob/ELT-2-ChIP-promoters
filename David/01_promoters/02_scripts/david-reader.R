library(magrittr)
library(readr)

###### Intestine gene categories
read_rob_intestine_gene_categories = function(){
rob.dir = normalizePath('../../../Rob/03_emb_L1_L3_intestine_RNAseq/03_output')
fileroot = file.path(rob.dir,'intestine_gene_categories')

intestine.gene.categories.fnames = list(LE='embryo_intestine_gene_categories.csv',
                                        L1='L1_intestine_gene_categories.csv',
                                        L3='L3_intestine_gene_categories.csv')

intestine.gene.categories = lapply(intestine.gene.categories.fnames, 
      function(f) { read_csv(file.path(fileroot, f)) })

intestine.gene.categories$LE %<>%  dplyr::rename(embryo_altHyp=altHyp,embryo_int_exp=intestine_expression)
intestine.gene.categories$L1 %<>%  dplyr::rename(L1_altHyp=altHyp,L1_int_exp=intestine_expression)
intestine.gene.categories$L3 %<>%  dplyr::rename(L3_altHyp=altHyp,L3_int_exp=intestine_expression)

# return single, merged dataframe
intestine.gene.categories$LE %>% full_join(intestine.gene.categories$L1, by = 'WBGeneID') %>%
                                        full_join(intestine.gene.categories$L3, by = 'WBGeneID')
}
read_rob_dineen_sets = function() {
  # wd: David/01_promoters/02_scripts
  robdir = normalizePath("../../../Rob")
  # Dineen results
  elt2_regulated_gene_sets <- read.table(file.path(robdir,
                                                   "05_elt2_RNAseq/03_output/elt2_regulated_gene_sets.csv"), 
                                         sep=",", header=T)
  # rob rerun of dineen
  res_elt2D_v_wt <- read.table(file.path(robdir,"05_elt2_RNAseq/03_output/res_elt2D_v_wt.csv"),sep=',',header=T)
  
  ELT2.din = full_join(res_elt2D_v_wt, elt2_regulated_gene_sets, by='WBGeneID')
  return(ELT2.din)
}

read_rob_ashr_shrunk_rlogc = function() {
  # wd: David/01_promoters/02_scripts
  
  # Also including Rlog normalized counts:  # via from David/01_promoters/02_scripts
  robdir = normalizePath("../../../Rob")
  rob.counts.path = file.path(robdir, '03_emb_L1_L3_intestine_RNAseq/03_output/rlog_counts/GFPplus_samples_rlog_counts.tsv')
  
  rob.counts = read.table(rob.counts.path, header=T) %>% dplyr::rename(
    embryo_rep1.rlogc=embryo_GFPplus_rep1,   
    embryo_rep2.rlogc=embryo_GFPplus_rep2,   
    embryo_rep3.rlogc=embryo_GFPplus_rep3,   
    L1_rep1.rlogc=L1_GFPplus_rep1,
    L1_rep3.rlogc=L1_GFPplus_rep3,     
    L3_rep1.rlogc=L3_GFPplus_rep1,      
    L3_rep2.rlogc=L3_GFPplus_rep2, 
    L3_rep3.rlogc=L3_GFPplus_rep3 
  ) %>% rowwise() %>% mutate(
    rlogc.embryo = mean(c(embryo_rep1.rlogc,embryo_rep2.rlogc,embryo_rep3.rlogc)),
    rlogc.L1     = mean(c(L1_rep1.rlogc,L1_rep3.rlogc)),
    rlogc.L3     = mean(c(L3_rep1.rlogc,L3_rep2.rlogc,L3_rep3.rlogc))
  )
  
  
  rob.dir = normalizePath('../../../Rob/03_emb_L1_L3_intestine_RNAseq/03_output')
  rob.shrunk.files = list(embryo='pairwise_shrunk_DE_results/res_embryoGFPplus_vs_embryoGFPminus_ashr_shrunk.csv',
                          L1='pairwise_shrunk_DE_results/res_L1GFPplus_vs_L1GFPminus_ashr_shrunk.csv',
                          L3='pairwise_shrunk_DE_results/res_L3GFPplus_vs_L3GFPminus_ashr_shrunk.csv')
  
  
  shrunk = lapply(rob.shrunk.files, function(f)
  {
    read.csv( file.path(rob.dir,f) )
  })
  
  merged = shrunk$embryo %>% full_join(shrunk$L1, by = "WBGeneID", suffix=c(".embryo", ".L1"))
  merged %<>% full_join( shrunk$L3 %>% dplyr::rename(
    baseMean.L3 = baseMean,
    log2FoldChange.L3 = log2FoldChange,
    lfcSE.L3 = lfcSE,
    pvalue.L3 = pvalue,
    padj.L3 = padj
  ))
  
  merged %<>% full_join(rob.counts %>% dplyr::select(WBGeneID,starts_with("rlogc.")), by = "WBGeneID")
  
  
  return(merged)
}

read_rob_all_merged = function() {
  a = read_rob_ashr_shrunk_rlogc()
  b = read_rob_dineen_sets() %>% 
    dplyr::select(WBGeneID, status, description, log2FoldChange) %>% 
    dplyr::rename(din.status=status, din.status.description=description, din.log2FoldChange = log2FoldChange)
  
  #return
  full_join(a,b,by='WBGeneID') %>% full_join(read_rob_intestine_gene_categories(), by="WBGeneID")
}

read_ELT2_binding_data = function(as_genomic_ranges=FALSE) {
  LE_promoter_tsv = "../03_output/LE.promoters.hilo.tsv"
  
  if(!file.exists(LE_promoter_tsv)) {
    stop("%s can't be read. Working directory must be .../ELT-2-ChIP-revision/David/01_promoters/02_scripts", LE_promoter_tsv)
  }
  
  LE_tsv = read.table(LE_promoter_tsv, header=T)
  L1_promoter_tsv = "../03_output/L1.promoters.hilo.tsv"
  L1_tsv = read.table(L1_promoter_tsv, header=T)
  L3_promoter_tsv = "../03_output/L3.promoters.hilo.tsv"
  L3_tsv = read.table(L3_promoter_tsv, header=T)
  
  LE_tsv$stage = "LE"
  L1_tsv$stage = "L1"
  L3_tsv$stage = "L3"
  
  cbound = LE_tsv %>% dplyr::arrange(wbps_gene_id) %>% 
    dplyr::select(seqnames, start, end, width, strand, wbps_gene_id, log_chip_signal_mean,
                  log_chip_signal_max,
                  IDR_logTEN_max,
                  IDR_logTEN_mean,
                  IDR_logTEN_sum,
                  class
    ) %>% 
    dplyr::rename(LE_wbps_gene_id = wbps_gene_id,
                  LE.log_chip_signal_mean=log_chip_signal_mean,
                  LE.log_chip_signal_max=log_chip_signal_max,
                  LE.IDR_logTEN_max=IDR_logTEN_max,
                  LE.IDR_logTEN_mean=IDR_logTEN_mean,
                  LE.IDR_logTEN_sum=IDR_logTEN_sum,
                  LE.class = class
    ) %>%
    cbind(L1_tsv %>% dplyr::arrange(wbps_gene_id) %>% 
            dplyr::select(wbps_gene_id, log_chip_signal_mean,
                          log_chip_signal_max,
                          IDR_logTEN_max,
                          IDR_logTEN_mean,
                          IDR_logTEN_sum,
                          class
            ) %>%
            dplyr::rename(L1_wbps_gene_id = wbps_gene_id,
                          L1.log_chip_signal_mean=log_chip_signal_mean,
                          L1.log_chip_signal_max=log_chip_signal_max,
                          L1.IDR_logTEN_max=IDR_logTEN_max,
                          L1.IDR_logTEN_mean=IDR_logTEN_mean,
                          L1.IDR_logTEN_sum=IDR_logTEN_sum,
                          L1.class=class
            )
    ) %>%
    cbind(L3_tsv %>% dplyr::arrange(wbps_gene_id) %>% 
            dplyr::select(wbps_gene_id, log_chip_signal_mean,
                          log_chip_signal_max,
                          IDR_logTEN_max,
                          IDR_logTEN_mean,
                          IDR_logTEN_sum,
                          class
            ) %>%
            dplyr::rename(L3_wbps_gene_id=wbps_gene_id,
                          L3_log.chip_signal_mean=log_chip_signal_mean,
                          L3.log_chip_signal_max=log_chip_signal_max,
                          L3.IDR_logTEN_max=IDR_logTEN_max,
                          L3.IDR_logTEN_mean=IDR_logTEN_mean,
                          L3.IDR_logTEN_sum=IDR_logTEN_sum,
                          L3.class = class
            )
    )
  
  stopifnot(all(cbound$LE_wbps_gene_id == cbound$L1_wbps_gene_id) &&
              all(cbound$LE_wbps_gene_id == cbound$L3_wbps_gene_id))
  
  cbound = cbound %>% 
    dplyr::select(-LE_wbps_gene_id,-L1_wbps_gene_id) %>%
    dplyr::rename(WBGeneID=L3_wbps_gene_id) 
  
  cbound = cbound %>% mutate(LE_bound=is.finite(LE.IDR_logTEN_max),
                             L1_bound=is.finite(L1.IDR_logTEN_max), 
                             L3_bound=is.finite(L3.IDR_logTEN_max))
  
  if(as_genomic_ranges) {
    return(GenomicRanges::makeGRangesFromDataFrame(cbound,keep.extra.columns = T))
  }
  return(cbound)
}
