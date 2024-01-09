
# PACKAGES ----------------------------------------------------------------
library(synapser)
library(igraph)
library(tidygraph)
library(ggraph)
library(tidyverse)

source(here::here('..','treatAD_biodomains','scripts',
                  'biodom_tally.R'))
source(here::here('..','treatAD_biodomains','scripts',
                  'annotate_enr_term_biodomains.R'))

theme_set(theme_bw())


# FUNCTIONS ---------------------------------------------------------------
read_kd <- function(path) {
  
  # detect if KDA dir path is on synapse
  if (substr(path, 1, 3) == 'syn') {
    
    synDir <- synGetChildren(path)$asList() %>% tibble(f = .) %>% unnest_wider(f)
    
    kd.ef0 <- synDir %>% 
      filter(grepl('edgeFactor_0.results.txt', name)) %>% 
      pull(id) %>% synGet() %>% .$path %>% 
      read_tsv(col_types = cols())
    
    ef0.top <- synDir %>% 
      filter(grepl('edgeFactor_0.tophits.txt', name)) %>% 
      pull(id) %>% synGet() %>% .$path %>% 
      read_tsv(col_types = cols()) %>%
      mutate(kd = paste0(NODE, '-', MODULE))
   
    kd.ef1 <- synDir %>% 
      filter(grepl('edgeFactor_1.results.txt', name)) %>% 
      pull(id) %>% synGet() %>% .$path %>% 
      read_tsv(col_types = cols())
    
    ef1.top <- synDir %>% 
      filter(grepl('edgeFactor_1.tophits.txt', name)) %>% 
      pull(id) %>% synGet() %>% .$path %>% 
      read_tsv(col_types = cols()) %>%
      mutate(kd = paste0(NODE, '-', MODULE))

  # or read from a local path
  } else {
    
    kd.pth <- list.files(path, full.names = T)
    
    if(length(kd.pth) < 13){
      warning(paste0('\nNot full kda result set:\n', path))
      return( kd = tibble(NULL))
    }
    
    kd.ef0 <- kd.pth %>% 
      str_subset('edgeFactor_0.results.txt') %>% 
      read_tsv(col_types = cols())
    
    ef0.top <- kd.pth %>% 
      str_subset('edgeFactor_0.tophits.txt') %>% 
      read_tsv(col_types = cols()) %>%
      mutate(kd = paste0(NODE, '-', MODULE))
    
    kd.ef1 <- kd.pth %>% 
      str_subset('edgeFactor_1.results.txt') %>% 
      read_tsv(col_types = cols())
    
    ef1.top <- kd.pth %>% 
      str_subset('edgeFactor_1.tophits.txt') %>% 
      read_tsv(col_types = cols()) %>%
      mutate(kd = paste0(NODE, '-', MODULE))
    
  }
  
  # Join up the results
  kd = bind_rows(
    kd.ef0 %>% mutate(ef = 'ef0_fdr'),
    kd.ef1 %>% mutate(ef = 'ef1_fdr')
    ) %>%
    mutate(FDR = -log10(FDR)) %>%
    pivot_wider(
      id_cols = c(MODULE, NODE, MEMBER),
      names_from = ef,
      values_from = FDR,
      values_fill = 0
    ) %>%
    mutate(delta_fdr = ef1_fdr - ef0_fdr) %>% 
    arrange(desc(delta_fdr)) %>%
    mutate(top_ef1 = 0, top_ef0 = 0)
  
  # Indicate which driver-module pairs are top for each analysis
  if(any(ef1.top$FDR < 0.05)) {
    kd$top_ef1[which(paste0(kd$NODE, '-', kd$MODULE) %in%
                       ef1.top$kd[ef1.top$FDR < 0.05])] <- 1
  }
  if(any(ef0.top$FDR < 0.05)) {
    kd$top_ef0[which(paste0(kd$NODE, '-', kd$MODULE) %in%
                       ef0.top$kd[ef0.top$FDR < 0.05])] <- 1
  }
  
  kd <- kd %>% rowwise() %>%
    mutate(top_kd = case_when(
      (top_ef1 == 1 & top_ef0 == 1) ~ 'both',
      (top_ef1 == 1 & top_ef0 == 0) ~ 'ef1',
      (top_ef1 == 0 & top_ef0 == 1) ~ 'ef0',
      (top_ef1 == 0 & top_ef0 == 0) ~ 'neither'
    ))
  
  return(kd)
}

integrate_kd <- function(path) {
  
  nodes <- list.files(dirname(path), full.names = T) %>%
    str_subset('graphml') %>%
    str_subset('bdFiltered', negate = T) %>%
    read_graph(format = 'graphml') %>%
    v_info() %>%
    mutate(n_nodes = nrow(.))
  
  nw_terms <- biodom %>%
    group_by(GOterm_Name, symbol, n_symbol) %>%
    summarise(n_nw = length(intersect(symbol, nodes$name)), .groups = 'keep') %>%
    mutate(pct_term = 100 * (n_nw / n_symbol)) %>%
    distinct() %>% 
    ungroup()
  
  kd.res <- read_kd(path)
  
  if (nrow(kd.res) == 0) {
    kd.res = tibble(NODE = NA_character_, MODULE = NA_character_)
  }
  
  kd.res <- kd.res %>%
    left_join(., nodes %>% select(NODE = name
                                  , n_nodes
                                  , node_status
                                  # , node_inBD
                                  # , min_path_length
                                  )
              , by = 'NODE') %>%
    left_join(.,
              biodom %>% 
                select(MODULE = GOterm_Name
                       , Biodomain:TopSD
                       , abbr:color)
              ,by = 'MODULE') %>%
    left_join(.,
              nw_terms %>% 
                select(MODULE = GOterm_Name
                       , n_symbol:pct_term) %>% distinct()
              ,by = 'MODULE') %>%
    mutate(nw = basename(dirname(path))) %>%
    mutate(pct_nodes = (n_nw / n_nodes) * 100) %>%
    select(
      nw
      , n_nodes
      , contains('_fdr')
      , NODE
      , node_status
      # , node_inBD
      # , min_path_length
      , contains('top_')
      , contains('MEMBER')
      , MODULE
      , Biodomain
      , Subdomain
      , subdomain_idx
      , TopSD
      , pct_nodes
      , pct_term
      , n_nw
      , n_symbol
      , abbr
      , label
      , color
    )
  }

simp <- list(
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
  directed = 'max',
  dijkstra_dist = 'min'
)



# CRISPR BRAIN DATA -------------------------------------------------------

synLogin()

# CRISPRbrain data
syn.dir <- synGetChildren('syn51183357')$asList() %>% 
  tibble(f = .) %>% unnest_wider(f)

# Simple Screens
simple <- map_dfr(
  syn.dir %>% with(which( !grepl('RNA|CROP',name) )),
  ~ { 
    id = syn.dir$id[.x]
    n = syn.dir$name[.x]
    synGet(id)$path %>% 
      read_csv( col_types = cols() ) %>% 
      mutate(cell = n %>% str_split_fixed(.,'-(?!M)', 2) %>% .[,1],
             mode = n %>% str_extract('CRISPR[ain]'),
             expt = n %>% str_remove_all('-CRISPR[ain].csv') %>% 
               str_split_fixed('-(?!M)',2) %>% .[,2] ) %>% 
      relocate(cell, mode, expt)
  }
  ) %>% 
  rename(pval = `P Value`, pheno = `Phenotype`, gene_score = `Gene Score`)

# Add extended phenotype descriptor
simple <- tibble(expt = simple$expt %>% unique() %>% sort(),
                 phenotype = c('Expansion (CD34 Staining)',
                               'Reactive Oxygen Species (CellRox Intensity)',
                               'Proliferation (CFSE Staining)',
                               'Survival (14 day)', 'Survival (21 day)', 'Survival (28 day)',
                               'Labile Iron (FeRhoNox Intensity)',
                               'Immune Activation (CD38 levels)',
                               'Lysosome Exocytosis (cell surface LAMP1), +IL1a +TNF +C1q',
                               'Lysosome Exocytosis (cell surface LAMP1), Vehicle',
                               'Peroxidized Lipids (Liperfluo intensity)',
                               'Lysosome (LysoTracker intensity)',
                               'Lysosome Mass or pH (LysoTracker staining), +IL1a +TNF +C1q',
                               'Lysosome Mass or pH (LysoTracker staining), Vehicle',
                               'Survival, no antioxidants',
                               'Phagocytosis (pHRodo-rat synaptosomes)',
                               'Phagocytosis (pHRodo-labeled synaptosomes), +IL1a +TNF +C1q',
                               'Phagocytosis (pHRodo-labeled synaptosomes), Vehicle',
                               'Survival, PSAP KO & no antioxidants',
                               'Survival, PSAP KO',
                               'Survival',
                               'Survival-Proliferation',
                               'Inflammatory activation (cell-surface VCAM1 levels), +IL1a +TNF +C1q')
                 ) %>% left_join(simple, . , by = 'expt')

# Group related phenotypes
simple <- simple %>% 
  mutate(
    pheno2 = case_when(
      grepl('rolif|xpansion',phenotype) ~ 'proliferation',
      grepl('urvival', phenotype) ~ 'survival',
      grepl('hagocytosis',phenotype) ~ 'phagocytosis',
      grepl('LysoTracker',phenotype) ~ 'lysotracker',
      grepl('ysosome',phenotype) ~ 'lysosome_exocytosis',
      grepl('ctivation',phenotype) ~ 'inflammation',
      grepl('FeRhoNox',phenotype) ~ 'FeRhoNox',
      grepl('CellRox',phenotype) ~ 'CellRox',
      grepl('Liperfluo',phenotype) ~ 'LiperFluo'
    )
  )

# Simplify and call hits vs non-hits per experiment group
simple1 <- simple %>% 
  filter( cell %in% c("Glutamatergic Neuron","iAstrocyte","iTF-Microglia" )) %>% 
  mutate(cb_hit = if_else(`Hit Class` == 'Non-Hit', 0, 1)) %>% 
  select(Gene
         , pheno2
         # , expt
         , cb_hit) %>% 
  group_by(Gene
           , pheno2
           # , expt
           ) %>% 
  mutate(cb_hit = if_else( any(cb_hit == 1), 1, 0)) %>% 
  # renamge(pheno2 = expt) %>% 
  distinct()

# KD RESULTS DIRS ---------------------------------------------------------

results_dir <- here::here('results')

kd_paths <- tibble(
  dir = c(
    'biodomain',
    'biodomain_dijkstra',
    'subdomain',
    'subdomain_dijkstra' ),
  kd = list(
    c('kda_genetics', 'kda_smallMod_test', 'kda_subdom_modules'),
    #'kda',
    c('kda_', 'kda_full_add', 'kda_full_mult', 'kda_smallMod_test'),
    c('kda'),
    c('kda_full_add', 'kda_full_mult', 'kda_smallMod_test') )
  ) %>%  
  mutate(nw_dirs = map(
    dir,
    ~ if (grepl('sub', .x)) {
      list.dirs(here::here(results_dir, .x), recursive = F) %>% 
        list.dirs(recursive = F)
      } else {
        list.dirs(here::here(results_dir, .x), recursive = F)
        })) %>%
  relocate(kd, .after = nw_dirs) %>%
  unnest(kd)

# PROCESS RESULTS ---------------------------------------------------------

# p = list()
cb_lr = list()
kd_res = list()

for (i in 1:nrow(kd_paths)) {
  
  ###
  # READ KD
  ###
  
  kd = map_dfr(
    kd_paths$nw_dirs[[i]],
    ~ paste0(.x, '/', kd_paths$kd[i]) %>% integrate_kd()
    )
  
  kd_res[[i]] <- kd
  
  ###
  # CBRAIN OVERLAP
  ###
  
  cbrain_lr = map_dfr(
    kd$nw %>% unique(),
    ~ {
      # load NW and get node degree
      nw <- kd_paths$nw_dirs[[i]] %>%
        str_subset(.x) %>%
        list.files(recursive = F, full.names = T) %>%
        str_subset('graphml') %>%
        str_subset('bdFilt', negate = T) %>%
        read_graph(format = 'graphml')
      
      snw <-igraph::simplify(nw,
                             remove.multiple = T,
                             edge.attr.comb = simp)
      
      nodes <- v_info(nw) %>%
        mutate(n_nodes = nrow(.),
               node_degree = degree(snw))
      
      # join CB phenos, NW nodes, and KDA results
      tmp <- inner_join( simple1, 
                  nodes %>% select(Gene = name, node_degree),
                  by = 'Gene') %>%
        left_join(
          .,
          kd %>%
            filter(nw == .x
                   , ef1_fdr >= -log10(5e-2)            ##    TUNING
                   # , top_kd %in% c('ef1','both')        ##    KDA
                   # , delta_fdr > 0                      ##    RESULT
                   # , pct_term >= 50                     ##    FILTER
                   ) %>% 
            select(Gene = NODE, node_status) %>% 
            distinct(),
          by = 'Gene'
          )
      
      # ID drivers, perform logistic regression
      tmp %>% 
        mutate(
          key.driver = if_else(is.na(node_status), 0, 1),
         # node_degree = scale(node_degree, center = F, scale = T) %>% c(),
         # node_degree = max(node_degree) - node_degree
         ) %>%
        group_by(pheno2) %>%
        mutate(
          lr = glm(cb_hit ~ key.driver + node_degree,
                   family = binomial(link = 'logit')) %>% broom::tidy() %>%
            filter(term == 'key.driver') %>%
            mutate(nobs = glm(cb_hit ~ key.driver + node_degree,
                              family = binomial(link = 'logit')) %>%
                     broom::glance() %>% pull(nobs)
                   ) %>% list() ) %>%
        unnest(lr)  %>%
        select(-Gene, -cb_hit, -key.driver, -node_degree, -node_status) %>%
        distinct() %>%
        ungroup() %>%
        mutate(nw = .x)
      
    })
  
  cb_lr[[i]] <- cbrain_lr
                               
   ###
   # PLOT
   ###
  
   # p[[i]] = ggplot(cbrain_lr, aes(statistic, nw)) +
   #   geom_vline(xintercept = 0, lwd = .2) +
   #   theme(
   #     strip.text.y = element_text(angle = 0),
   #     legend.position = 'bottom'
   #   ) +
   #   scale_size_continuous('# observations', limits = c(0, max(cbrain_lr$nobs))) +
   #   viridis::scale_fill_viridis(
   #     direction = -1,
   #     na.value = 'white',
   #     limits = c(-log10(0.05), max(-log10(
   #       cbrain_lr$p.value
   #     )))
   #   ) +
   #   geom_point(
   #     shape = 21,
   #     aes(fill = -log10(p.value), size = nobs),
   #     alpha = .5
   #   ) +
   #   labs(x = '',
   #        title = paste0(kd_paths$dir[i], ' -- ', kd_paths$kd[i]))+
   #   facet_wrap(~ pheno2)
   
}

kd_paths$kd_results = kd_res
kd_paths$cb_logit = cb_lr

kd_paths$set = paste0(kd_paths$dir, ' -- ', kd_paths$kd)

tmp <- kd_paths %>% select(set, cb_logit) %>% unnest(cb_logit)
tmp1 <- tmp %>% 
  filter(pheno2 == 'survival'
         , grepl('biodomain',set)) 


p = map(
  tmp$pheno2 %>% unique, 
  ~ {
    tmp1 <- tmp %>% filter(pheno2 == .x, grepl('biodomain',set)) 
    
    ggplot(tmp1, aes(statistic, set)) +
      geom_vline(xintercept = 0, lwd = .2) +
      theme(
        strip.text.y = element_text(angle = 0),
        legend.position = 'bottom'
      ) +
      scale_size_continuous(
        '# observations', 
        limits = c(0, max(tmp1$nobs))
      ) +
      viridis::scale_fill_viridis(
        direction = -1,
        na.value = 'white',
        limits = c(-log10(0.05), max(-log10(tmp1$p.value)))
      ) +
      geom_point(
        shape = 21,
        aes(fill = -log10(p.value), size = nobs),
        alpha = .5
      ) +
      labs(x = ''
           , title = .x
      )+
      facet_wrap(~ nw)
  }
)

