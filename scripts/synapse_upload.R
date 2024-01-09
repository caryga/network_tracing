# setup -------------------------------------------------------------------

library(synapser)
library(tidyverse)

synLogin()

# function to read directory on synapse into a tibble
read_syn_dir <- function(parent_id){
  synapser::synGetChildren(parent_id)$asList() %>%
    tibble(f = .) %>% unnest_wider(f)
}

# biodomains --------------------------------------------------------------

# enumerate results directory on synapse
parent_id = 'syn51739831' # biodomains

# enumerate results directory on synapse
synLogin()
syn.dir <- read_syn_dir(parent_id)

# top results directory
local.dir <- tibble(top = list.dirs(here::here('results','biodomain'), recursive = F) ) %>%
  mutate(dom = basename(top)) %>%
  rowwise() %>% mutate( nw_f = list.files(top) %>% list(),
                        kda_f = list.files(paste0(top,'/','kda')) %>% list(),
                        kda_f_2 = list.files(paste0(top,'/','kda_subdom_modules')) %>% list()
                        )

# prepare synapse directories and upload data
for( i in 1:nrow(local.dir) ){

  if( local.dir$dom[i] %in% syn.dir$name ){
    syn.i <- which(syn.dir$name == local.dir$dom[i])
    id <- syn.dir$id[syn.i]
  } else{
    foo <- synStore( Folder( local.dir$dom[i] , parent = parent_id) )
    id <- foo$properties$id
  }

  # # graph files
  # foo <- map(
  #   1:length(local.dir$nw_f[[i]]),
  #   ~ synStore( File(
  #     paste0(local.dir$top[i],'/',local.dir$nw_f[[i]][.x]),
  #     parent = id
  #   )))

  # # kda dir on synapse?
  # if( 'kda' %in% read_syn_dir(id)$name ){
  #   kda_id <-  read_syn_dir(id) %>% filter(name == 'kda') %>% pull(id)
  # } else {
  #   foo2 <- synStore( Folder( 'kda' , parent = id) )
  #   kda_id <- foo2$properties$id
  # }
  # 
  # # kda files
  # foo3 <- map(
  #   1:length(local.dir$kda_f[[i]]),
  #   ~ synStore( File(
  #     paste0(local.dir$top[i],'/kda/',local.dir$kda_f[[i]][.x]),
  #     parent = kda_id
  #   )))
  
  # kda dir on synapse?
  if( 'kda_subdom_modules' %in% read_syn_dir(id)$name ){
    kda_id <-  read_syn_dir(id) %>% filter(name == 'kda_subdom_modules') %>% pull(id)
  } else {
    foo2 <- synStore( Folder( 'kda_subdom_modules' , parent = id) )
    kda_id <- foo2$properties$id
  }

  # kda files
  foo3 <- map(
    1:length(local.dir$kda_f_2[[i]]),
    ~ synStore( File(
      paste0(local.dir$top[i],'/kda_subdom_modules/',local.dir$kda_f_2[[i]][.x]),
      parent = kda_id
    )))

}

# subdomains --------------------------------------------------------------


# enumerate results directory on synapse
parent_id = 'syn52314422' # subdomains

cat(paste0(
  'uploading files from: ', here::here('results','subdomain'),
  '\n','to: ', parent_id,
  '\n'
))

# top results directory
local.dir <- tibble(top = list.dirs(here::here('results','subdomain'), recursive = F) ) %>% 
  mutate(dom = basename(top)) %>% 
  rowwise() %>% mutate( 
    subdom = list.files(top) %>% list()
  ) %>% 
  unnest_longer(subdom) %>% 
  rowwise() %>% mutate(
    nw_f = paste0(top,'/',subdom) %>% list.files() %>% str_subset('kda',negate = T) %>% list(), 
    kda_f = paste0(top,'/',subdom,'/kda') %>% list.files() %>% list() 
  )

cat(paste0(
  '# directories = ',nrow(local.dir),'\n'))

# prepare synapse directories and upload data
for( i in 1:nrow(local.dir) ){
  
  cat(paste0(
    'starting subdomain #', i,': ' ,
      local.dir$dom[i], ' / ', local.dir$subdom[i],' ... '
    ))
  
  # enumerate results directory on synapse
  syn.dir <- read_syn_dir(parent_id) %>% 
    filter(grepl('Folder', type)) %>% 
    select(bd_name = name, bd_id = id) %>% 
    rowwise() %>% 
    mutate(subdom = read_syn_dir(bd_id) %>% list()) 
  
  cat('which bd directory / sd directory ... ','\n')
  
  if( local.dir$dom[i] %in% syn.dir$bd_name ){
    bd_idx <- which(syn.dir$bd_name == local.dir$dom[i])
    bd_id <- syn.dir$bd_id[bd_idx]
  } else {
    foo <- synStore( Folder( local.dir$dom[i] , parent = parent_id) )
    bd_id <- foo$properties$id
  }
  
  if( local.dir$subdom[i] %in% syn.dir$subdom[[bd_idx]]$name ){
    sd_idx <- which(syn.dir$subdom[[bd_idx]]$name == local.dir$subdom[i])
    sd_id <- syn.dir$subdom[[bd_idx]]$id[sd_idx]
  } else {
    foo <- synStore( Folder( local.dir$subdom[i] , parent = bd_id) )
    sd_id <- foo$properties$id
  }
  
  cat('upload files ... ','\n')
  
  # graph files
  foo <- map(
    1:length(local.dir$nw_f[[i]]),
    ~ synStore( File( 
      paste0(local.dir$top[i],'/',local.dir$subdom[i],'/',local.dir$nw_f[[i]][.x]),
      parent = sd_id
    )))
  
  # kda dir on synapse?
  if( 'kda' %in% read_syn_dir(sd_id)$name ){
    
    kda_id <-  read_syn_dir(sd_id) %>% filter(name == 'kda') %>% pull(id)
    
  } else {
    
    foo2 <- synStore( Folder( 'kda' , parent = sd_id) )
    kda_id <- foo2$properties$id
    
  }
  
  # kda files
  foo3 <- map(
    1:length(local.dir$kda_f[[i]]),
    ~ synStore( File(
      paste0(local.dir$top[i],'/',local.dir$subdom[i],'/kda/',local.dir$kda_f[[i]][.x]),
      parent = kda_id
    )))

  }

# # cell type specific biodomain NWs ----------------------------------------
# 
# # enumerate results directory on synapse
# parent_id = 'syn52383741' # biodomain by cell type
# 
# # enumerate results directory on synapse
# syn.dir <- synapser::synGetChildren(parent_id)$asList() %>%
#   tibble(f = .) %>% unnest_wider(f)
# 
# # get the local data dirs
# base <- paste0(here::here(),'/results/biodomain_by_celltype')
# nm1 <- c('path','d1','d2', 'nw')
# dirs <- tibble(path = list.dirs(base, recursive = T)) %>%
#   mutate(d = str_split_fixed(path, 'biodomain_by_celltype',2) %>% .[,2]) %>%
#   filter(!grepl('kda',d), d!='') %>%
#   mutate(d = str_split(d, '\\/')) %>%
#   unnest_wider(d, names_repair = ~nm1) %>% suppressMessages %>%
#   select(-path, -d1) %>%
#   filter(!is.na(nw)) %>%
#   rowwise() %>%
#   mutate( data = list.files( paste0(base,'/',d2,'/',nw),recursive = F) %>% list(),
#           kda = list.files( paste0(base,'/',d2,'/',nw,'/','kda'), recursive = F) %>% list())
# 
# # prepare synapse directories
# top <- dirs %>% select(d2) %>% distinct() %>%
#   mutate(foo = synStore(Folder(d2, parent = parent_id)) %>% list(),
#          id = pluck(foo, 'properties') %>% pluck(., 'id'))
# 
# dirs <- left_join(dirs, top %>% select(-foo)) %>% relocate(id, .after = 'd2')
# 
# mid <- dirs %>%
#   mutate(foo = synStore(Folder(nw, parent = id)) %>% list(),
#          id_nw = pluck(foo, 'properties') %>% pluck(., 'id'),
#          foo2 = synStore(Folder('kda', id_nw)) %>% list(),
#          id_kda = pluck(foo2, 'properties') %>% pluck(.,'id'))
# 
# dirs <- left_join(dirs, mid %>% select(nw, id_nw, id_kda))
# 
# # store the data
# nw_dir <- dirs %>% unnest_longer(data) %>% filter(data != 'kda') %>%
#   rowwise() %>%
#   mutate(foo = synStore( File( paste0(base,'/',d2,'/',nw,'/',data), parent = id_nw)) %>% list())
# 
# kda_dir <- dirs %>% unnest_longer(kda) %>%
#   rowwise() %>%
#   mutate(foo = synStore(File(paste0(base,'/',d2,'/',nw,'/kda/',kda), parent = id_kda)) %>% list())
# 
# ##
# # remove duplicates across celltype directories
# read_syn_dir <- function(parent_id){
#   synapser::synGetChildren(parent_id)$asList() %>%
#     tibble(f = .) %>% unnest_wider(f)
# }
# 
# parent_id = 'syn52383741' # biodomain by cell type
# syn.dir <- read_syn_dir(parent_id)
# 
# dirs = map_dfr(1:nrow(syn.dir), 
#                ~ read_syn_dir(syn.dir$id[.x]) %>% 
#                  mutate(parent = syn.dir$name[.x]) %>% 
#                  relocate(parent)) %>%
#   rowwise() %>% 
#   mutate(nw.dir = read_syn_dir(id) %>% list(),
#          kd.dir = read_syn_dir(id) %>% filter(name == 'kda') %>% pull(id) %>% 
#            read_syn_dir(.) %>% list(), 
#          celltype = str_split_fixed(parent, '_', 2)[,2] %>% str_to_sentence()) 
# 
# dirs <- dirs %>% 
#   select(parent, celltype, n_name = name, n_id = id, nw.dir, kd.dir) %>% 
#   pivot_longer(cols = c('nw.dir', 'kd.dir'), names_repair = 'minimal') %>% 
#   rename(dir.name = name) %>% 
#   unnest(value) 
# 
# ids_to_remove <- dirs %>% 
#   rowwise() %>% 
#   filter((grepl(celltype, name)|name == 'kda'|grepl('queryList',name))==F) %>% 
#   pull(id)
# 
# foo <- map(ids_to_remove, ~synDelete(.x))
# 
# ##
# # upload re-filtered MM graphs
# dom = 'Mitochondrial_Metabolism'
# mm.syndir = dirs %>% filter(name == dom) %>% 
#   mutate(celltype = str_split_fixed(parent, '_',2)[,2]) %>% 
#   relocate(celltype)
# 
# mm.graphs <- list.files('results/biodomain_by_celltype', 
#                         pattern = '^Mitochondrial_Metabolism.*graphml', 
#                         recursive = T, 
#                         full.names = T)
# 
# bd_genes <- biodom %>% 
#   filter(Biodomain == str_replace_all(dom, '_',' '),
#          GO_ID != 'GO:0043227' # membrane-bound organelle w/ 12k genes annotated
#   ) %>%
#   pull(symbol) %>% 
#   unlist() %>% 
#   unique() %>% 
#   .[!is.na(.)]
# 
# foo <- map(
#   mm.graphs,
#   ~{
#     trace.nw <- igraph::read_graph(.x, format = 'graphml') 
#     trace_nw_genes <- names(igraph::V(trace.nw))
#     trace_nw_filter <- setdiff(trace_nw_genes, bd_genes)
#     trace.nw.bdFilt <- igraph::delete_vertices(
#       trace.nw,
#       v=igraph::V(trace.nw)[ names(igraph::V(trace.nw)) %in% trace_nw_filter ]
#     ) # PREVIOUS: igraph::induced_subgraph( trace.nw, v )
#     igraph::write_graph(
#       trace.nw.bdFilt,
#       paste0(dirname(.x), '/bdFiltered_', basename(.x)),
#       format = "graphml"
#     )
#     celltype = str_extract_all(.x, 'astro|micro|exc|inh') %>% unlist()
#     id = mm.syndir$id[which(mm.syndir$celltype == celltype)]
#     synStore(File(paste0(dirname(.x), '/bdFiltered_', basename(.x)), id))
#   }
# )

# pseudotime --------------------------------------------------------------

# # enumerate results directory on synapse
# parent_id = 'syn52392218' # pseudotime
# 
# # # enumerate results directory on synapse
# # syn.dir <- synapser::synGetChildren(parent_id)$asList() %>%
# #   tibble(f = .) %>% unnest_wider(f)
# 
# # get the local data dirs
# base <- paste0(here::here(),'/results/pseudotime')
# nm1 <- c('path','d1','nw')
# dirs <- tibble(path = list.dirs(base, recursive = T)) %>%
#   mutate(d = str_split_fixed(path, 'pseudotime',2) %>% .[,2]) %>%  
#   filter(!grepl('kda',d), d!='') %>% 
#   mutate(d = str_split(d, '\\/')) %>% 
#   unnest_wider(d, names_repair = ~nm1) %>% suppressMessages %>% 
#   select(-path, -d1) %>% 
#   filter(!is.na(nw)) %>% 
#   rowwise() %>% 
#   mutate( nw_data = list.files( paste0(base,'/',nw),recursive = F) %>% list(),
#           kda_data = list.files( paste0(base,'/',nw,'/','kda'), recursive = F) %>% list())
# 
# # prepare synapse directories
# top <- dirs %>% select(nw) %>% distinct() %>% 
#   mutate(foo = synStore(Folder(nw, parent = parent_id)) %>% list(),
#          id_nw = pluck(foo, 'properties') %>% pluck(., 'id'),
#          foo2 = synStore(Folder('kda', id_nw)) %>% list(), 
#          id_kda = pluck(foo2, 'properties') %>% pluck(.,'id'))
# 
# dirs <- left_join(dirs, top %>% select(-foo, -foo2)) %>% relocate(id_nw,id_kda, .after = 'nw')
# 
# # store the data
# nw_dir <- dirs %>% unnest_longer(nw_data) %>% filter(nw_data != 'kda') %>% 
#   rowwise() %>% 
#   mutate(foo = synStore( File( paste0(base,'/',nw,'/',nw_data), parent = id_nw)) %>% list())
# 
# kda_dir <- dirs %>% unnest_longer(kda_data) %>% 
#   rowwise() %>% 
#   mutate(foo = synStore(File(paste0(base,'/',nw,'/kda/',kda_data), parent = id_kda)) %>% list())

# EOF ---------------------------------------------------------------------
