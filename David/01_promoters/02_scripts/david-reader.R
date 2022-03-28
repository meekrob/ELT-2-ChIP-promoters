read_ELT2_binding_data = function(as_genomic_ranges=FALSE) {
  LE_promoter_tsv = "../03_output/LE.promoters.hilo.tsv"
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