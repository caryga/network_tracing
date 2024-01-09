

# setup -------------------------------------------------------------------

library(synapser)
library(tidyverse)

synLogin()


# kda.dir <- map_dfr(
#   syn.dir %>% filter(grepl('Folder', type)) %>% pull(id),
#   ~ {
#     id = synGetChildren(.x)$asList() %>% tibble(f = .) %>% unnest_wider(f) %>%
#       filter(grepl('Folder', type)) %>% pull(id)
#     synGetChildren(id)$asList() %>% tibble(f=.) %>% unnest_wider(f)
#     }
# )

# local dirs to push ------------------------------------------------------

# # enumerate results directory on synapse
# # parent_id <- 'syn51117833'
# parent_id = 'syn52314422' # subdomains
# 
# # enumerate results directory on synapse
# syn.dir <- synapser::synGetChildren(parent_id)$asList() %>%
#   tibble(f = .) %>% unnest_wider(f)


# subdomains --------------------------------------------------------------

# # enumerate results directory on synapse
# parent_id = 'syn52314422' # subdomains
# 
# # enumerate results directory on synapse
# syn.dir <- synapser::synGetChildren(parent_id)$asList() %>%
#   tibble(f = .) %>% unnest_wider(f)
# 
# # top results directory
# top <- list.dirs(paste0(here::here(),'/results/subdomain'), recursive = F)
# 
# # get the data dirs
# d <- map( top, ~ list.dirs( .x, recursive = F ) %>% str_remove_all( paste0(.x,'/') ) ) %>%
#   setNames(., top %>% str_remove_all('^.*results/subdomain/'))
# 
# # prepare synapse directories and upload data
# for( i in 1:length(d) ){
# 
#   # top-level directory
#   foo <- synStore( Folder( d[i] %>% names() , parent = parent_id) )
# 
#   for( j in 1:length(d[[i]]) ){
# 
#     # graph directory
#     foo2 <- synStore( Folder( d[[i]][j], parent = foo$properties$id) )
# 
#     # list graph files
#     f <- list.files( paste0( top[i],'/',d[[i]][j], '/' ) ) %>% str_subset('kda',negate=T)
# 
#     # upload graph files
#     for( k in 1:length(f) ){
#       foo3 <- synStore( File(
#         paste0( top[i],'/',d[[i]][j], '/', f[k] ),
#         parent=foo2$properties$id
#       ))
#     }
# 
#     ### KDA RESULTS ###
# 
#     # is there a kda directory?
#     if( !dir.exists(paste0( top[i],'/',d[[i]][j], '/','kda')) ){ next }
# 
#     # create kda sub-directory
#     foo3 <- synStore( Folder('kda', parent = foo2$properties$id) )
# 
#     # list kda files
#     f <- list.files( paste0( top[i],'/',d[[i]][j], '/kda') )
# 
#     # store kda files
#     for(k in 1:length(f)){
#       foo4 <- synStore( File(
#         paste0(paste0( top[i],'/',d[[i]][j], '/kda/', f[k])),
#         parent=foo3$properties$id
#       ))
#     }
# 
#   }
# }

# re-upload re-filtered MM networks
synLogin()
parent_id = 'syn52314422' # subdomains
syn.dir <- synapser::synGetChildren(parent_id)$asList() %>%
  tibble(f = .) %>% unnest_wider(f)

mm.syndir <- synGetChildren(
  syn.dir$id[syn.dir$name == 'Mitochondrial_Metabolism'])$asList() %>%
  tibble(f = .) %>% unnest_wider(f)

mm.dir <- tibble( path = list.dirs('results/subdomain/Mitochondrial_Metabolism', recursive = F) ) %>%
  rowwise() %>%
  mutate(base = basename(path),
        files = list.files(path,full.names = T) %>% list()) %>%
  unnest_longer(files) %>% filter(grepl('bdFilt', files))

for(i in 1:nrow(mm.syndir)){
  idx = which(mm.dir$base == mm.syndir$name[i])
  foo <- synStore( File( mm.dir$files[idx], parent = mm.syndir$id[i] ) )
}

for(i in 1:nrow(mm.dir)){
  idx = which(mm.syndir$name == mm.dir$base[i])
  kda.syndir <- synGetChildren(mm.syndir$id[idx])$asList() %>%
    tibble(f = .) %>% unnest_wider(f) %>% 
    filter(grepl('kda', name)) %>% pull(id)
  files <- list.files(paste0(mm.dir$path[i],'/kda'), full.names = T)
  foo <- map(files, ~ synStore( File( .x, parent = kda.syndir ) ))
}

# cell type specific biodomain NWs ----------------------------------------

# enumerate results directory on synapse
parent_id = 'syn52383741' # biodomain by cell type

# enumerate results directory on synapse
syn.dir <- synapser::synGetChildren(parent_id)$asList() %>%
  tibble(f = .) %>% unnest_wider(f)

# get the local data dirs
base <- paste0(here::here(),'/results/biodomain_by_celltype')
nm1 <- c('path','d1','d2', 'nw')
dirs <- tibble(path = list.dirs(base, recursive = T)) %>%
  mutate(d = str_split_fixed(path, 'biodomain_by_celltype',2) %>% .[,2]) %>%
  filter(!grepl('kda',d), d!='') %>%
  mutate(d = str_split(d, '\\/')) %>%
  unnest_wider(d, names_repair = ~nm1) %>% suppressMessages %>%
  select(-path, -d1) %>%
  filter(!is.na(nw)) %>%
  rowwise() %>%
  mutate( data = list.files( paste0(base,'/',d2,'/',nw),recursive = F) %>% list(),
          kda = list.files( paste0(base,'/',d2,'/',nw,'/','kda'), recursive = F) %>% list())

# prepare synapse directories
top <- dirs %>% select(d2) %>% distinct() %>%
  mutate(foo = synStore(Folder(d2, parent = parent_id)) %>% list(),
         id = pluck(foo, 'properties') %>% pluck(., 'id'))

dirs <- left_join(dirs, top %>% select(-foo)) %>% relocate(id, .after = 'd2')

mid <- dirs %>%
  mutate(foo = synStore(Folder(nw, parent = id)) %>% list(),
         id_nw = pluck(foo, 'properties') %>% pluck(., 'id'),
         foo2 = synStore(Folder('kda', id_nw)) %>% list(),
         id_kda = pluck(foo2, 'properties') %>% pluck(.,'id'))

dirs <- left_join(dirs, mid %>% select(nw, id_nw, id_kda))

# store the data
nw_dir <- dirs %>% unnest_longer(data) %>% filter(data != 'kda') %>%
  rowwise() %>%
  mutate(foo = synStore( File( paste0(base,'/',d2,'/',nw,'/',data), parent = id_nw)) %>% list())

kda_dir <- dirs %>% unnest_longer(kda) %>%
  rowwise() %>%
  mutate(foo = synStore(File(paste0(base,'/',d2,'/',nw,'/kda/',kda), parent = id_kda)) %>% list())

##
# remove duplicates across celltype directories
read_syn_dir <- function(parent_id){
  synapser::synGetChildren(parent_id)$asList() %>%
    tibble(f = .) %>% unnest_wider(f)
}

parent_id = 'syn52383741' # biodomain by cell type
syn.dir <- read_syn_dir(parent_id)

dirs = map_dfr(1:nrow(syn.dir), 
               ~ read_syn_dir(syn.dir$id[.x]) %>% 
                 mutate(parent = syn.dir$name[.x]) %>% 
                 relocate(parent)) %>%
  rowwise() %>% 
  mutate(nw.dir = read_syn_dir(id) %>% list(),
         kd.dir = read_syn_dir(id) %>% filter(name == 'kda') %>% pull(id) %>% 
           read_syn_dir(.) %>% list(), 
         celltype = str_split_fixed(parent, '_', 2)[,2] %>% str_to_sentence()) 

dirs <- dirs %>% 
  select(parent, celltype, n_name = name, n_id = id, nw.dir, kd.dir) %>% 
  pivot_longer(cols = c('nw.dir', 'kd.dir'), names_repair = 'minimal') %>% 
  rename(dir.name = name) %>% 
  unnest(value) 

ids_to_remove <- dirs %>% 
  rowwise() %>% 
  filter((grepl(celltype, name)|name == 'kda'|grepl('queryList',name))==F) %>% 
  pull(id)

foo <- map(ids_to_remove, ~synDelete(.x))

##
# upload re-filtered MM graphs
dom = 'Mitochondrial_Metabolism'
mm.syndir = dirs %>% filter(name == dom) %>% 
  mutate(celltype = str_split_fixed(parent, '_',2)[,2]) %>% 
  relocate(celltype)

mm.graphs <- list.files('results/biodomain_by_celltype', 
                        pattern = '^Mitochondrial_Metabolism.*graphml', 
                        recursive = T, 
                        full.names = T)

bd_genes <- biodom %>% 
  filter(Biodomain == str_replace_all(dom, '_',' '),
         GO_ID != 'GO:0043227' # membrane-bound organelle w/ 12k genes annotated
  ) %>%
  pull(symbol) %>% 
  unlist() %>% 
  unique() %>% 
  .[!is.na(.)]

foo <- map(
  mm.graphs,
  ~{
    trace.nw <- igraph::read_graph(.x, format = 'graphml') 
    trace_nw_genes <- names(igraph::V(trace.nw))
    trace_nw_filter <- setdiff(trace_nw_genes, bd_genes)
    trace.nw.bdFilt <- igraph::delete_vertices(
      trace.nw,
      v=igraph::V(trace.nw)[ names(igraph::V(trace.nw)) %in% trace_nw_filter ]
    ) # PREVIOUS: igraph::induced_subgraph( trace.nw, v )
    igraph::write_graph(
      trace.nw.bdFilt,
      paste0(dirname(.x), '/bdFiltered_', basename(.x)),
      format = "graphml"
    )
    celltype = str_extract_all(.x, 'astro|micro|exc|inh') %>% unlist()
    id = mm.syndir$id[which(mm.syndir$celltype == celltype)]
    synStore(File(paste0(dirname(.x), '/bdFiltered_', basename(.x)), id))
  }
)

# pseudotime --------------------------------------------------------------

# enumerate results directory on synapse
parent_id = 'syn52392218' # pseudotime

# # enumerate results directory on synapse
# syn.dir <- synapser::synGetChildren(parent_id)$asList() %>%
#   tibble(f = .) %>% unnest_wider(f)

# get the local data dirs
base <- paste0(here::here(),'/results/pseudotime')
nm1 <- c('path','d1','nw')
dirs <- tibble(path = list.dirs(base, recursive = T)) %>%
  mutate(d = str_split_fixed(path, 'pseudotime',2) %>% .[,2]) %>%  
  filter(!grepl('kda',d), d!='') %>% 
  mutate(d = str_split(d, '\\/')) %>% 
  unnest_wider(d, names_repair = ~nm1) %>% suppressMessages %>% 
  select(-path, -d1) %>% 
  filter(!is.na(nw)) %>% 
  rowwise() %>% 
  mutate( nw_data = list.files( paste0(base,'/',nw),recursive = F) %>% list(),
          kda_data = list.files( paste0(base,'/',nw,'/','kda'), recursive = F) %>% list())

# prepare synapse directories
top <- dirs %>% select(nw) %>% distinct() %>% 
  mutate(foo = synStore(Folder(nw, parent = parent_id)) %>% list(),
         id_nw = pluck(foo, 'properties') %>% pluck(., 'id'),
         foo2 = synStore(Folder('kda', id_nw)) %>% list(), 
         id_kda = pluck(foo2, 'properties') %>% pluck(.,'id'))

dirs <- left_join(dirs, top %>% select(-foo, -foo2)) %>% relocate(id_nw,id_kda, .after = 'nw')

# store the data
nw_dir <- dirs %>% unnest_longer(nw_data) %>% filter(nw_data != 'kda') %>% 
  rowwise() %>% 
  mutate(foo = synStore( File( paste0(base,'/',nw,'/',nw_data), parent = id_nw)) %>% list())

kda_dir <- dirs %>% unnest_longer(kda_data) %>% 
  rowwise() %>% 
  mutate(foo = synStore(File(paste0(base,'/',nw,'/kda/',kda_data), parent = id_kda)) %>% list())

# EOF ---------------------------------------------------------------------
