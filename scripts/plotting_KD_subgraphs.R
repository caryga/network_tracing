dom <- 'Synapse'
kd.res <- kd %>% filter(nw == dom)
nw = read_nw(which(doms == dom))
snw = igraph::simplify(nw, remove.multiple = T)

# ID KDs: sig & top ef1, delta > 0 
kd.nodes <- kd.res %>% 
  filter(ef1_fdr > -log10(0.05)
         , top == 1
         , delta_fdr > 0
  ) %>% 
  select(NODE, node_status, MODULE, Biodomain, MEMBER, top, contains('fdr')) %>% 
  arrange(desc(delta_fdr), desc(ef1_fdr)) %>%
  group_by(NODE) %>% 
  summarise(n_biodom = length(unique(Biodomain))
         , n_term = length(unique(MODULE))
         , ef1_max = max(ef1_fdr)
         , delta_max = max(delta_fdr)) 


x = kd.nodes %>% arrange(desc(n_term), desc(delta_max)) %>% slice(1:20)
tmp <- kd.res %>% 
  filter(NODE %in% x$NODE
         , ef1_fdr > -log10(0.05)
         , top == 1
         , delta_fdr > 0) %>% 
  group_by(NODE, Biodomain, color, node_status) %>% 
  summarise(n_GO = length(unique(MODULE)),
            ef1_max = max(ef1_fdr),
            delta_max = max(delta_fdr)) %>% 
  mutate(NODE = factor(NODE, level = x$NODE)) %>% 
  arrange(NODE)

ggplot(tmp, aes(delta_max, NODE))+
  geom_point(shape = 21, aes(fill = color, size = n_GO), alpha = .5, color = 'black')+
  scale_fill_identity()+scale_y_discrete(limits = rev)


x = kd %>% 
  filter(ef1_fdr > -log10(0.05)
         , top == 1
         , delta_fdr > 0) %>% 
  select(NODE, node_status, MODULE, Biodomain, MEMBER, top, contains('fdr'), nw) %>% 
  arrange(desc(delta_fdr), desc(ef1_fdr)) %>%
  group_by(NODE) %>% 
  summarise(n_biodom = length(unique(Biodomain))
            , n_term = length(unique(MODULE))
            , n_nw = length(unique(nw))
            , ef1_max = max(ef1_fdr)
            , delta_max = max(delta_fdr)) %>% 
  arrange(desc(n_term), desc(delta_max)) %>% 
  slice(1:20)

tmp <- kd %>% 
  filter(NODE %in% x$NODE
         , top == 1
         , delta_fdr > 0
         , ef1_fdr > -log10(0.05)  
         # , MODULE %in% enr.bd$pathway[enr.bd$padj < 0.05]
  ) %>% 
  group_by(NODE, Biodomain, color, node_status) %>% 
  summarise(
    n_nw = length(unique(nw)),
    n_terms = length(unique(MODULE)),
    # n_biodom = length(unique(Biodomain)),
    delta_max = max(delta_fdr),
    ef1_max = max(ef1_fdr)
  ) %>% 
  ungroup() %>% 
  mutate(NODE = factor(NODE, level = x$NODE)) %>% 
  arrange(NODE)

ggplot(tmp, aes(delta_max, NODE))+
  geom_jitter(shape = 21, aes(fill = color, size = n_nw), alpha = .5, color = 'black')+
  scale_fill_identity()+scale_y_discrete(limits = rev)


# KD nodes ----------------------------------------------------------------

# subnetwork of just KD Nodes
vert_num <- match(kd.nodes$NODE, V(snw)$name)

g2 <- induced_subgraph(snw, vert_num)
unlinked <- names(degree(g2)[which(degree(g2)==0)])
g2 <- delete.vertices(g2, unlinked)

# KD evidence
x = tibble( node = V(g2) %>% names ) %>%
  left_join(., kd.nodes, by = c('node'='NODE'))
for(i in 2:ncol(x)){
  igraph::vertex_attr(g2, names(x)[i], index = igraph::V(g2)) <- x %>% pull(i)
}

# plot
ggraph(g2, layout = 'lgl') +
  geom_edge_link(alpha = .3, color = 'grey75') +
  geom_node_point( color = "black", alpha = .5, 
                   aes(size = n_GO, fill = delta_max, shape = node_status )  )+
  geom_node_text(aes(label = name), size = 2) +
  scale_shape_manual( values = c(21,23) )+
  scale_size_continuous( range = c(2,10) )+
  viridis::scale_fill_viridis()+
  theme_graph(background = 'white')

# KD & neighbors ----------------------------------------------------------

# subnetwork of KD and neighbors
vert_neigh <- neighbors(snw, vert_num)

g3 <- induced_subgraph(snw, c(vert_num, vert_neigh))
unlinked <- names(degree(g3)[which(degree(g3)==0)])
g3 <- delete.vertices(g3, unlinked)

# KD evidence
x = tibble( node = V(g3) %>% names ) %>% 
  left_join(., kd.nodes, by = c('node'='NODE')) %>% 
  mutate(kd = if_else(is.na(n_GO), 0, 1))
for(i in 2:ncol(x)){
  igraph::vertex_attr(g3, names(x)[i], index = igraph::V(g3)) <- x %>% pull(i)
}

ggraph(g3, layout = 'graphopt') +
  geom_edge_link(alpha = .3, color = 'grey90') +
  geom_node_point( shape = 21, color = "black", alpha = .3, fill = 'grey90', 
                   aes(filter= .N()$kd == 0)  )+
  geom_node_point( shape = 23, color = "black", alpha = .6, 
                   aes(size = delta_max, fill = ef1_max, filter = .N()$kd == 1) )+
  scale_shape_manual(name='Driver Status', breaks=c('key driver', 'non-driver'),
                     values=c('key driver'=23, 'non-driver'=21))+
  geom_node_text(aes(label = name, filter = .N()$kd == 1), size = 2) +
  viridis::scale_fill_viridis()+
  theme_graph(background = 'white')

