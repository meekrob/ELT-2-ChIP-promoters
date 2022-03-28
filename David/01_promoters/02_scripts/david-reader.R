read_rob_dineen_sets = function() {
  # wd: David/01_promoters/02_scripts
  robdir = normalizePath("../../../Rob")
  # Dineen results
  elt2_regulated_gene_sets <- read.table(file.path(robdir,
                                                   "05_elt2_RNAseq/03_output/elt2_regulated_gene_sets.csv"), 
                                         sep=",", header=T)
  # rob rerun of dineen
  res_elt2D_v_wt <- read.table(file.path(robdir,"05_elt2_RNAseq/03_output/res_elt2D_v_wt.csv"),sep=',',header=T)
  
  ELT2.din = inner_join(res_elt2D_v_wt, elt2_regulated_gene_sets, by='WBGeneID')
  return(ELT2.din)
}

read_rob_ashr_shrunk = function() {
  # wd: David/01_promoters/02_scripts
  ##### Rlog normalized counts:  # via from David/01_promoters/02_scripts
  robdir = normalizePath("../../../Rob")
  rob.counts.path = file.path(robdir, '03_emb_L1_L3_intestine_RNAseq/03_output/rlog_counts/GFPplus_samples_rlog_counts.tsv')
  
  rob.counts = read.table(rob.counts.path, header=T) 
  
  
  rob.dir = normalizePath('../../../Rob/03_emb_L1_L3_intestine_RNAseq/03_output')
  rob.shrunk.files = list(LE='pairwise_shrunk_DE_results/res_embryoGFPplus_vs_embryoGFPminus_ashr_shrunk.csv',
                          L1='pairwise_shrunk_DE_results/res_L1GFPplus_vs_L1GFPminus_ashr_shrunk.csv',
                          L3='pairwise_shrunk_DE_results/res_L3GFPplus_vs_L3GFPminus_ashr_shrunk.csv')
  
  
  shrunk = lapply(rob.shrunk.files, function(f)
  {
    read.csv( file.path(rob.dir,f) )
  })
  
  shrunk$LE = shrunk$LE %>% inner_join(rob.counts, by = "WBGeneID") %>% mutate(LE.rlog.rep1 = rob.counts$embryo_GFPplus_rep1,
                                   LE.rlog.rep2 = rob.counts$embryo_GFPplus_rep2,
                                   LE.rlog.rep3 = rob.counts$embryo_GFPplus_rep3)
  
  shrunk$L1 = shrunk$LE %>% mutate(L1.rlog.rep1 = rob.counts$L1_GFPplus_rep1,
                                   L1.rlog.rep3 = rob.counts$L1_GFPplus_rep3)
  
  
  return(shrunk)
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