## functions to read KDA results both locally and from synapse

# read_nw <- function(synID){
#   nw = read.graph( synGet( synID )$path,format = 'graphml' )
#   return(nw)
#   }

# from synapse
read_kd <- function(synID_res_dir){
  
  synDir <- synGetChildren(synID_kda_dir)$asList() %>% tibble(f = .) %>% unnest_wider(f) 
  synDir <- synDir %>% filter( name == 'kda' ) %>% pull(id) %>% 
    synGetChildren(.) %>% .$asList() %>% tibble(f = .) %>% unnest_wider(f) %>% 
    bind_rows(synDir,.)

  queryNodes <- synDir %>% filter( grepl('queryList', name) ) %>% 
    pull(id) %>% synGet() %>% .$path %>% read_tsv(., col_types = cols(), col_names = 'gene')
  kd.ef0 <- synDir %>% filter( grepl('edgeFactor_0.results.txt', name) ) %>% 
    pull(id) %>% synGet() %>% .$path %>% read_tsv(.,col_types = cols())
  kd.ef1 <- synDir %>% filter( grepl('edgeFactor_1.results.txt', name) ) %>% 
    pull(id) %>% synGet() %>% .$path %>% read_tsv(.,col_types = cols())
  kd.top <- synDir %>% filter( grepl('edgeFactor_1.tophits.txt', name) ) %>% 
    pull(id) %>% synGet() %>% .$path %>% read_tsv(.,col_types = cols()) %>% 
    mutate(kd = paste0(NODE,'-',MODULE))

  kd = bind_rows(
    kd.ef0 %>% mutate(ef = 'ef0_fdr'),
    kd.ef1 %>% mutate(ef = 'ef1_fdr')
  ) %>% 
    mutate(FDR = -log10(FDR)) %>%
    pivot_wider(id_cols = c(MODULE, NODE, MEMBER), names_from = ef, values_from = FDR, values_fill = 0) %>%
    mutate(delta_fdr = ef1_fdr - ef0_fdr) %>% arrange(desc(delta_fdr)) %>% 
    mutate(
      node_status = if_else(NODE %in% queryNodes$gene, 'query', 'added'), 
      top = 0
      )
  
  if( any( kd.top$FDR < 0.05 ) ){
    kd$top[ which(paste0(kd$NODE,'-',kd$MODULE) %in% kd.top$kd[kd.top$FDR < 0.05]) ] <- 1
  }
  
  return(kd)
}

# from local directory
read_kd_local <- function(path){
  
  resDir <- tibble( name = list.files(path, full.names = T) )
  resDir <- tibble( name = list.files(paste0(path,'/kda'), full.names = T) ) %>% 
    bind_rows(resDir, .)

  queryNodes <- resDir %>% filter( grepl('queryList', name) ) %>% 
    pull(name) %>% read_tsv(., col_types = cols(), col_names = 'gene')
  kd.ef0 <- resDir %>% filter( grepl('edgeFactor_0.results.txt', name) ) %>% 
    pull(name) %>% read_tsv(.,col_types = cols())
  kd.ef1 <- resDir %>% filter( grepl('edgeFactor_1.results.txt', name) ) %>% 
    pull(name) %>% read_tsv(.,col_types = cols())
  kd.top <- resDir %>% filter( grepl('edgeFactor_1.tophits.txt', name) ) %>% 
    pull(name) %>% read_tsv(.,col_types = cols()) %>% 
    mutate(kd = paste0(NODE,'-',MODULE))

  kd = bind_rows(
    kd.ef0 %>% mutate(ef = 'ef0_fdr'),
    kd.ef1 %>% mutate(ef = 'ef1_fdr')
  ) %>% 
    mutate(FDR = -log10(FDR)) %>%
    pivot_wider(id_cols = c(MODULE, NODE, MEMBER), names_from = ef, values_from = FDR, values_fill = 0) %>%
    mutate(delta_fdr = ef1_fdr - ef0_fdr) %>% arrange(desc(delta_fdr)) %>% 
    mutate(
      node_status = if_else(NODE %in% queryNodes$gene, 'query', 'added'), 
      top = 0
    )
  
  if( any( kd.top$FDR < 0.05 ) ){
    kd$top[ which(paste0(kd$NODE,'-',kd$MODULE) %in% kd.top$kd[kd.top$FDR < 0.05]) ] <- 1
  }
  
  return(kd)
}
