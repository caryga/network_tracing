# initialize biodomain NW trace queries and directories


# setup -------------------------------------------------------------------

# Package names
packages <- c('synapser','tidyverse')

# Load packages
suppressPackageStartupMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)

synLogin()

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

subdomains <- biodom %>% select(Biodomain, subdomain_idx, Subdomain) %>% 
  distinct() %>% filter(subdomain_idx != 0)

# # enriched biodomain terms
# trs.enr <- read_csv(synGet('syn45824995')$path, col_types = cols()) %>% 
#   mutate(leadingEdge_genes = str_split(leadingEdge_genes, '\\|'))
# gen.enr <- read_csv(synGet('syn45824969')$path, col_types = cols()) %>% 
#   mutate(leadingEdge_genes = str_split(leadingEdge_genes, '\\|'))
# omic.enr <- read_csv(synGet('syn45824835')$path, col_types = cols()) %>% 
#   mutate(leadingEdge_genes = str_split(leadingEdge_genes, '\\|'))
# 
# enr.bd <-  trs.enr

enr.scores <- readRDS('/projects/carter-lab/caryg/treatAD_biodomains/results/risk_score_biodom_enrichments.rds')
enr.bd <- enr.scores$anno[[3]] %>% 
  mutate(leadingEdge_genes = str_split(core_enrichment, '\\/'))

# biodomain ---------------------------------------------------------------

# initialize directories and leading edge gene list files
setwd('/projects/carter-lab/caryg/network_tracing/results/biodomain')
for( bd in domains ){
  
  # pull gene list from leading edge genes
  gl <- enr.bd %>% 
    filter(
      Biodomain == bd,
      p.adjust < 0.01,
      NES > 1.7
    ) %>% 
    pull(leadingEdge_genes) %>% 
    unlist() %>% 
    unique()

  cat(paste0(bd, ': ', length(gl), ' leading edge genes\n'))
  
  if(length(gl) == 0){next}
  
  bd_filename <- bd %>% str_replace_all(.,' ','_')
  if(!dir.exists(bd_filename)){ dir.create(bd_filename) }
  
  # write query list to file
  write_tsv(tibble(x = gl), paste0(bd_filename, '/queryList_', bd_filename, '.tsv'), col_names = F)
  
}


# sub-domain --------------------------------------------------------------

# read curated list of sub-domain terms 
subdom <- readxl::read_xlsx(
  paste0(here::here(), '/data/Subendophenotype Trackers.xlsx'),
  sheet = 1, 
  col_names = T
)

subdom$Biodomain[which(subdom$Biodomain == 'Epigenetics')] <- 'Epigenetic'

# join subdomain and enrichment results table
# rename sub-domain with top term for that sub-domain
subdom <- full_join( 
    enr.bd, 
    subdom %>% select(Biodomain, ID = GO_ID, Subdomain, TopSD),
    by = c('Biodomain','ID')) %>% 
  group_by(Biodomain, Subdomain) %>% 
  mutate(sub_dom = pathway[TopSD == 1]) %>% 
  ungroup()

# initialize directories and leading edge gene list files
setwd('/projects/carter-lab/caryg/network_tracing/results/subdomain')
domains <- subdom %>% filter(!is.na(sub_dom)) %>% pull(Biodomain) %>% unique() %>% sort()
for( bd in domains ){
  
  # full domain gene list
  gl1 <- subdom %>% 
    filter(
      Biodomain == bd,
      padj < 0.01,
      NES > 1.7
    ) %>% 
    pull(leadingEdge_genes) %>% 
    unlist() %>% 
    unique()
  
  cat('\n\nfull domain: \n', 
      paste0(bd, ': ', length(gl1), ' leading edge genes\n\n'))
  
  # list sub-domains with signif terms
  sds <- subdom %>% 
    filter(Biodomain == bd, padj < 0.01, NES > 1.7, !is.na(sub_dom)) %>% 
    pull(sub_dom) %>% unique() %>% sort()
  for( sd in sds ){
    
    # pull gene list from leading edge genes
    gl <- subdom %>% 
      filter(
        sub_dom == sd,
        padj < 0.01,
        NES > 1.7
      ) %>% 
      pull(leadingEdge_genes) %>% 
      unlist() %>% 
      unique()
    
    cat(paste0(sd, ': ', length(gl), ' leading edge genes\n'))
    
    if(length(gl) == 0){next}

    bd_filename <- bd %>% str_replace_all(.,' ','_')
    if(!dir.exists(bd_filename)){ dir.create(bd_filename) }
    sd_filename <- sd %>% str_replace_all(.,' ','_')
    if(!dir.exists( paste0(bd_filename,'/', sd_filename) )){ dir.create( paste0(bd_filename,'/', sd_filename) ) }

    # write query list to file
    write_tsv(tibble(x = gl), paste0(bd_filename,'/', sd_filename, '/queryList_', sd_filename, '.tsv'), col_names = F)
    
  }
  
}

# # M42 query list ----------------------------------------------------------
# 
# system('wget -q --output-document=tmt_modules_supplement.xlsx https://www.biorxiv.org/content/biorxiv/early/2021/08/13/2021.04.05.438450/DC2/embed/media-2.xlsx ')
# tmt.modules <- readxl::read_xlsx('tmt_modules_supplement.xlsx', sheet = 4, skip = 2)
# tmt.modules <- tmt.modules %>%
#   mutate(prot = str_split_fixed(UniqueID, '\\|',2)[,1],
#          mod = str_split_fixed(kMEtableSortVector, '\\|',2)[,1]) %>% 
#   relocate(prot, mod, .after = UniqueID)
# 
# gl = tmt.modules %>% filter(grepl('M42', mod)) %>% pull(prot) %>% unique()
# write_tsv(tibble(x = gl), 'results/modules/m42/queryList_M42.tsv', col_names = F)


# pseudotime query lists --------------------------------------------------

synapser::synLogin()

# ROSMAP txomic PT
rna.f.de <- read_csv( synapser::synGet('syn39989123')$path, col_types = cols() ) %>% rename(gene = gene_names)
rna.f.enr <- read_tsv( synapser::synGet('syn47728345')$path, col_types = cols() )
rna.m.de <- read_csv( synapser::synGet('syn39990047')$path, col_types = cols() ) %>% rename(gene = gene_names)
rna.m.enr <- read_tsv( synapser::synGet('syn47728065')$path, col_types = cols() )

# ROSMAP proteomic PT
prot.f.de <- read_csv( synapser::synGet('syn40616521')$path, col_types = cols() ) %>% rename(gene = gene_short_name)
prot.f.enr <- read_csv( synapser::synGet('syn50881293')$path, col_types = cols() )
prot.m.de <- read_csv( synapser::synGet('syn40621972')$path, col_types = cols() ) %>% rename(gene = gene_short_name)
prot.m.enr <- read_csv( synapser::synGet('syn50881295')$path, col_types = cols() )

# pull gene list from leading edge genes
input.gene.list <- 
  bind_rows( prot.f.de %>% mutate(s = 'f'), 
             prot.m.de %>% mutate(s = 'm')
    ) %>% 
  pivot_wider(id_cols = c(gene_names, gene, state), 
              values_from = c(pvalue, effect), 
              names_from = s) %>% 
  mutate(
    pvalue_f = if_else( !is.finite(pvalue_f), 0, pvalue_f ),
    pvalue_f = if_else( pvalue_f == 0 , min(pvalue_f[pvalue_f>0], na.rm =T), pvalue_f),
    pvalue_m = if_else( !is.finite(pvalue_m), 0, pvalue_m ),
    pvalue_m = if_else( pvalue_m == 0 , min(pvalue_m[pvalue_m>0], na.rm =T), pvalue_m)
  ) %>% 
  rowwise() %>% 
  mutate(
    sig = case_when( 
      (pvalue_f <= 0.01 & pvalue_m > 0.01) ~ 'f_only',
      (pvalue_m <= 0.01 & pvalue_f > 0.01) ~ 'm_only',
      (pvalue_f > 0.01 & pvalue_m > 0.01) ~ 'neither',
      (pvalue_f <= 0.01 & pvalue_m <= 0.01) ~ 'both'),
    mn = mean(-log10(pvalue_f), -log10(pvalue_m), na.rm = T)
  ) %>%
  mutate( state = case_when(
    state %in% c(2,3,4) ~ 'early', 
    state %in% c(5,6) ~ 'mid',
    state %in% c(7) ~ 'late',
    T ~ as.character(state)) ) 
  

for( s in unique(input.gene.list$state) ){
  input.gene.list %>% 
    filter(state == s
           , sig == 'both'
           , sign(effect_f) == sign(effect_m) ) %>% 
    select(gene) %>% 
    distinct() %>% 
    write_tsv(., paste0('results/pseudotime/',s,'/queryList_',s,'.tsv'), col_names = F)
}
