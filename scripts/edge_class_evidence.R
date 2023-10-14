edge <- e_info(undir_node_edge)

edge <- edge %>% 
  mutate(edge_class = case_when(
    interaction %in% c('in-complex-with','interacts-with') ~ 'binding',
    interaction %in% c('controls-state-change-of','controls-phosphorylation-of') ~ 'modification',
    interaction %in% c('controls-expression-of') ~ 'expression'
  ))

ggplot(data = subset(edge, edge_class %in% c('binding','modification')), 
       aes(n_evidence, fill = edge_class))+ 
  geom_histogram(position = 'dodge')+scale_x_log10()

edge %>% 
  filter(edge_class %in% c('binding','modification')) %>% 
  pivot_wider(id_cols = edge, names_from = edge_class, values_from = n_edge_evidence) %>%
  ggplot(aes(modification, binding))+geom_point()
