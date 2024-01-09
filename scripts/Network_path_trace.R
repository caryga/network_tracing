# This Rscript performs two broad functions:
# 1. self-trace leading edge genes from enriched biodomain term


# setup -------------------------------------------------------------------

cat('
##################################
##        PATHWAY TRACING       ##
##################################
    ')

# Package names
packages <- c('synapser','igraph','tidyverse')

# Load packages
suppressPackageStartupMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)

# source path tracing functions
source(paste0(here::here(), '/scripts/igraph_NW_exp_functions.R'))

# # source Mergeomics
# source('https://raw.githubusercontent.com/jessicading/mergeomics/master/Mergeomics_Version_1.99.0.R')

theme_set(theme_bw())

cat('\npackages loaded:\n', packages, '\n')

# base network
synLogin()
net <- igraph::read_graph(synGet('syn51110930')$path, format = 'graphml') 

# biological domain annotations
biodom <- full_join(
  # biodomains
  readRDS(synGet('syn25428992')$path),
  # domain labels
  read_csv(synGet('syn26856828')$path,
           col_types = cols()),
  by = c('Biodomain'='domain')
) %>%
  mutate(Biodomain = case_when(Biodomain == 'none' ~ NA_character_, T ~ Biodomain))

domains <- biodom %>% pull(Biodomain) %>% unique() %>% sort() %>% .[!is.na(.)]


# read arguments ----------------------------------------------------------

# arguments
# 1. query filepath
# 2. trace directed 
# 3. node filter
# 4. edge filter
# 5. cell type

# parse args and establish settings
args <- commandArgs(trailingOnly = TRUE)

cat('\narguments:\n',args, '\n')

# parse args
query.file <- args[ which( grepl('\\.txt|\\.tsv', args) ) ] 
directed <- any( grepl('dir', args) )
filt_nodes <- any( grepl('node', args) )
filt_edges <- any( grepl('edge', args) )
noTrace <- any( grepl('noTrace', args) )
filt_cellType <- any( grepl('Exc|Inh|Astro|Micro', args) )
if(filt_cellType){ cellType <- args[ which( grepl('Exc|Inh|Astro|Micro', args) ) ] }

# report args
full_path = normalizePath(query.file) %>% dirname()
working_path = normalizePath(query.file) %>% dirname() %>% basename()
cat('\n\nMove to dir: ',  full_path , '\n')
setwd(full_path)
cat('Query file to trace: ', query.file, '\n')
cat('File found in current path?: ', basename(query.file) %in% list.files(), '\n')
cat('Name of working directory: ', working_path, '\n\n')

cat('PATH TRACE OPTIONS\nTrace directed edges?: ', directed, '\n')
cat('Filter nodes for brain expression?: ', filt_nodes, '\n')
cat('Filter edges for PMID evidence?: ', filt_edges, '\n')
cat('Filter nodes for expression in certain cells?: ', filt_cellType, '\n')
if(filt_cellType){ cat('Trace NW in: ', cellType, '\n') }
# cat('Skip pathway tracing and only run wKDA?: ', noTrace, '\n\n')

# filter base network -----------------------------------------------------

cat('\n\n','Filtering Pathway Commons network based on specifications...','\n')

# nodes & edges
if( filt_edges & filt_nodes ){
  tmp <- igraph::delete_vertices(
    net,
    igraph::V(net)[ igraph::V(net)$brain_exp == 0 ]
  )
  nw <- igraph::subgraph.edges(
    tmp,
    igraph::E(tmp)[ igraph::E(tmp)$n_edge_evidence > 1 ],
    delete.vertices = T 
  )
  filt = '_filt_node_edge'
} else if( filt_nodes ){
  nw <- igraph::delete_vertices(
    net,
    igraph::V(net)[ igraph::V(net)$brain_exp == 0 ]
  )
  filt = '_filt_node'
} else if( filt_edges ){
  nw <- igraph::subgraph.edges(
    net,
    igraph::E(net)[ igraph::E(net)$n_edge_evidence > 1 ],
    delete.vertices = T 
  )
  filt = '_filt_edge'
} else { 
  nw <- net
  filt = '_filt_none'
  }

# cell type
if( filt_cellType ){
  nw <- igraph::delete_vertices(
    nw,
    V(nw)[ which( vertex_attr(nw, cellType) %in% c(0, NaN) ) ]
  )
  filt = paste0(filt,'_',cellType)
}

# directionality
if( directed ){
  nw <- igraph::subgraph.edges( 
    nw, 
    igraph::E(nw)[ igraph::E(nw)$directed == 1 ],
    delete.vertices = T 
  ) %>% 
    as.directed(mode = 'arbitrary')
  directionality = '_directed'  
} else { 
  directionality = '_undirected'
}

cat( 'Base network filtered.', '\n')

# gene list to trace ------------------------------------------------------

# Read file
file_name = query.file %>% basename()

cat('\n','Tracing paths for genes in the file: ', query.file %>% basename(), 
    '\n At the path ', getwd(), ' \n\n')

query = read_tsv( query.file %>% basename(), col_names = 'gene')

# Intersect with network nodes
input.gene.list <- query$gene
if( length(setdiff( input.gene.list, names(V(nw)))) > 0 ){
  cat('\nQuery genes missing from filtered NW object: ',
      setdiff( input.gene.list, names(V(nw)) ) %>% sort(), sep = '\n')
  cat('\n')
  
  input.gene.list = intersect( 
    query$gene, names(V(nw))
    )
}

cat('Number of queried genes to trace: ', 
    length(input.gene.list), ' / ', length(query$gene),
    ' (', length(input.gene.list)/length(query$gene)*100, '%)',
    ' \n')

# path tracing ------------------------------------------------------------

cat('\n','Beginning trace... \n')
Sys.time()

# ID target biodomain from path name
dom = full_path %>% 
  str_split('\\/') %>% 
  unlist() %>% 
  str_replace_all(., '_',' ') %>% 
  str_subset(., paste0(domains, collapse = '|'))

bd_genes <- biodom %>% 
  filter(Biodomain == dom,
         GO_ID != 'GO:0043227' # membrane-bound organelle w/ 12k genes annotated
  ) %>%
  pull(symbol) %>% 
  unlist() %>% 
  unique() %>% 
  .[!is.na(.)]

# remove redundant edges
snw <- igraph::simplify(nw, 
                        remove.multiple = TRUE,
                        remove.loops = FALSE,
                        edge.attr.comb = list( 
                          interaction = 'concat',
                          edge = 'random',
                          occurrance = 'concat',
                          n_edge = 'max',
                          n_edge_types = 'max',
                          n_edge_evidence = 'max',
                          n_source = 'max',
                          sources = 'concat',
                          n_evidence = 'sum',
                          evidence_pmid = 'concat',
                          n_pathways = 'max',
                          pathway_names = 'concat',
                          directed = 'max')
                        )

# calculate edge weights for Dijkstra
edges <- e_info(snw)
verts <- v_info(snw)

node_weights = verts %>% 
  select(name, trs = TargetRiskScore) %>% 
  mutate(trs = if_else(is.nan(trs), 0, trs),
         bd_gene = if_else(name %in% bd_genes, 5, 0)) 

edge_weights = edges %>% 
  select(edge) %>% distinct() %>% 
  mutate(head = str_split_fixed(edge, ':',2)[,1], 
         tail = str_split_fixed(edge, ':',2)[,2]) %>% 
  left_join(., node_weights %>% select(head = name, h.trs = trs, h.bd = bd_gene)) %>% 
  left_join(., node_weights %>% select(tail = name, t.trs = trs, t.bd = bd_gene)) %>% 
  mutate(
    ht.bd = case_when(
      h.bd + t.bd == 10 ~ 5, 
      h.bd + t.bd ==  5 ~ 3,
      h.bd + t.bd ==  0 ~ 0),
    # weight = (10-h.trs-h.bd)+(10-t.trs-t.bd)
    weight = (5-h.trs)+(5-t.trs)+(5-ht.bd)
    )

# p = edge_weights %>%
#   mutate(set = case_when(h.bd == 5 & t.bd == 5 ~ 'both bd',
#                          h.bd == 5 & t.bd == 0 | h.bd == 0 & t.bd == 5 ~ 'one bd',
#                          h.bd == 0 & t.bd == 0 ~ 'neither bd')) %>%
#   ggplot(aes(weight, fill = set))+
#   geom_density(alpha = .3)
  # geom_histogram(position = 'dodge')
# ggsave(here::here())


weights = e_info(nw) %>% select(edge) %>% left_join(., edge_weights, by = 'edge')
  igraph::edge_attr(nw, 'dijkstra_dist', index = igraph::E(nw)) <- weights$weight

##
# trace paths in parallel

future::plan(strategy = 'multisession', workers = 10)
trace <- furrr::future_map(
  input.gene.list,
  ~ short_paths(
    tnet = snw,
    target = .x,
    targets = input.gene.list,
    sentinals = input.gene.list,
    
    #conduct a breadth first search
    edge_weights = NULL, 
    
    # #specify weights to perform Dijkstra's algorithm
    # edge_weights = edge_weights$weight, 
    cores = 1)
)
future::plan(strategy = 'sequential')

trace1 <- map_dfr(
  trace, 
  ~ tibble(
    target = list(.x$Inter), 
    sentinal = list(.x$Sentinal), 
    nodes = .x$Nodes %>% list())
  ) 
trace1$source = input.gene.list

trace2 = trace1 %>% unnest(nodes) %>% 
  group_by(name, status) %>% 
  summarise(min_path_length = min(min_path)) %>% 
  mutate(in_bd = map_lgl(name, ~ .x %in% bd_genes))

##
# Filter NW obj for traced nodes

# trace_filt <- unlist(trace) %>% unique()
trace_filt <- trace2$name

trace.nw <- igraph::induced_subgraph(
  nw,
  v=igraph::V(nw)[ names(igraph::V(nw)) %in% trace_filt ]
)

##
# Annotate nodes that are queried or added by trace

status = v_info(trace.nw) %>% select(name) %>% left_join(., trace2, by = 'name')
  igraph::vertex_attr(trace.nw, 'node_inBD', index = igraph::V(trace.nw)) <- status$in_bd
  igraph::vertex_attr(trace.nw, 'node_status', index = igraph::V(trace.nw)) <- status$status
  igraph::vertex_attr(trace.nw, 'min_path_length', index = igraph::V(trace.nw)) <- status$min_path_length

##
# save NW trace and filtered NW

  igraph::write_graph(
  trace.nw,
  paste0( full_path, '/',
          working_path, directionality, filt,
          '.graphml'),
  format = "graphml"
)

# ##
# # Remove redundant edges
# nw.simple <- igraph::simplify(
#   trace.nw,
#   remove.multiple = TRUE,
#   remove.loops = FALSE,
#   edge.attr.comb = list( 
#     interaction = 'concat',
#     edge = 'random',
#     occurrance = 'concat',
#     n_edge = 'max',
#     n_edge_types = 'max',
#     n_edge_evidence = 'max',
#     n_source = 'max',
#     sources = 'concat',
#     n_evidence = 'sum',
#     evidence_pmid = 'concat',
#     n_pathways = 'max',
#     pathway_names = 'concat',
#     directed = 'max'
#   )
# )

cat('\n','Trace complete and network saved. \n')
Sys.time()

# # generate plots ----------------------------------------------------------
# 
# cat('\n\nGenerating plots...\n')
# 
# # calculate network stats
# nw.stats <- tibble(  
#   network = c('full','pre-trace','traced','simplified'),
#   n_nodes = c( 
#     net %>% V %>% length,
#     nw %>% V %>% length,
#     trace.nw %>% V %>% length,
#     nw.simple %>% V %>% length
#     ),
#   n_edges = c(
#     net %>% E %>% length,
#     nw %>% E %>% length,
#     trace.nw %>% E %>% length,
#     nw.simple %>% E %>% length
#     ),
#   avg_path_length = c(
#     net %>% average.path.length, 
#     nw %>% average.path.length, 
#     trace.nw %>% average.path.length, 
#     nw.simple %>% average.path.length
#     ),
#   assortativity_coef = c(
#     net %>% assortativity(., types1 = V(.)),
#     nw %>% assortativity(., types1 = V(.)),
#     trace.nw %>% assortativity(., types1 = V(.)),
#     nw.simple %>% assortativity(., types1 = V(.))
#     ),
#   connected_components = c(
#     net %>% no.clusters, 
#     nw %>% no.clusters,
#     trace.nw %>% no.clusters, 
#     nw.simple %>% no.clusters
#     )
# )
# 
# write_csv(nw.stats,
#           paste0( full_path, '/',
#                   working_path, directionality, filt,
#                   '_netStats.csv' )
#           )
# 
# # plot network stats
# nw.stats %>% 
#   pivot_longer(cols = -network, names_to = 'properties', values_to = 'val') %>% 
#   mutate(properties = factor(properties, 
#                              levels = c('n_nodes',
#                                         'n_edges', 
#                                         'avg_path_length',
#                                         'assortativity_coef',
#                                         'connected_components')) 
#          , network = factor(network, levels = c('full',
#                                                 'pre-trace',
#                                                 'traced',
#                                                 'biodom_filtered',
#                                                 'simplified'))) %>% 
#   ggplot(aes(network, val)) +
#   geom_bar(stat = 'identity', position = 'dodge')+
#   theme(legend.position = 'top',
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
#   labs(y = '',x = '', subtitle = 'base network properties')+
#   facet_wrap(~properties, scales = 'free_y', ncol = 2)
# 
# ggsave(
#   paste0( full_path, '/',
#           working_path, directionality, filt,
#           '_netStats.pdf' )
# )

# # synapse upload ----------------------------------------------------------
# 
# parent_id <- 'syn51117833'
# 
# d <- list.dirs('results', recursive = F) %>% 
#   str_subset('input_gene_lists|wKDA', negate=T)
# 
# for(j in 4:5){
#   foo <- synStore( Folder(d[j] %>% str_remove_all('results/'), parent = parent_id) )
#   f <- list.files(d[j]) %>% str_subset('kda',negate=T)
#   for(i in 1:length(f)){
#     foo2 <- synStore( File(
#       paste0(d[j],'/',f[i]),
#       parent=foo$properties$id
#     ))
#   }
#   foo3 <- synStore( Folder('kda', parent = foo$properties$id) )
#   f <- list.files( paste0(d[j],'/','kda')) 
#   for(i in 1:length(f)){
#     foo4 <- synStore( File(
#       paste0(d[j],'/kda/',f[i]),
#       parent=foo3$properties$id
#     ))
#   }
# }

# EOF ##
