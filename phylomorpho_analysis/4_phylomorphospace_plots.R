##### Phylomorphospace plots ####
#Gen Morinaga: this is a script I wrote to make a phylomorphospace in ggplot2. This can be done normally using phytools::phylomorphospace() but I prefer using ggplot2.
gg_phylomorpho<-function(tree, dat){
  anc_dat<-data.frame(apply(dat, 2, fastAnc, tree=tree))
  colnames(anc_dat)<-colnames(dat)
  tips_nodes<-dat[tree$tip.label,]
  tips_nodes<-rbind(tips_nodes, anc_dat)
  tips_nodes$node_id<-1:(tree$Nnode + length(tree$tip.label))
  
  edge.dat<-data.frame(tree$edge) 
  colnames(edge.dat)<-c('node1', 'node2')
  edge_coords<-left_join(left_join(edge.dat, tips_nodes, by = c('node1' = 'node_id')), tips_nodes, by = c('node2'= 'node_id'))
  edge_coords<-merge(
    merge(edge.dat, tips_nodes, by.x = 'node1', by.y = 'node_id'),
    tips_nodes, by.x = 'node2', by.y = 'node_id'
  )
  code<-'Paste the following to produce a nice enough plot. First geom_point (for original data) needs to be modified before plotting. ggplot() + geom_segment(data = x$edge_coords, aes(x = trait1.x, xend = trait1.y, y = trait2.x, yend = trait2.y)) + geom_point(data = x$original, aes(x = trait1, y = trait2)) + geom_point(data = x$anc_dat, aes(x = trait1, y = trait2))'
  pm_data<-list(dat, anc_dat, edge_coords, code)
  names(pm_data)<-c('original', 'anc_dat', 'edge_coords', 'code')
  pm_data
}

#estimating the phylomorphospace coordinates for scales
x<-gg_phylomorpho(aedini.tree, data.frame(adult.mca$ind$coord) %>% select(Dim.1, Dim.2))

#some minor housekeeping steps to make the plot a little more informative by identifying which nodes belong to which clade
min(sapply(1:nrow(data.frame(t(combn(clade_a, 2)))), function (x) findMRCA(aedini.tree, tips = as.character(data.frame(t(combn(clade_a, 2)))[x,]))))
max(sapply(1:nrow(data.frame(t(combn(clade_a, 2)))), function (x) findMRCA(aedini.tree, tips = as.character(data.frame(t(combn(clade_a, 2)))[x,]))))
min(sapply(1:nrow(data.frame(t(combn(clade_b, 2)))), function (x) findMRCA(aedini.tree, tips = as.character(data.frame(t(combn(clade_b, 2)))[x,]))))
max(sapply(1:nrow(data.frame(t(combn(clade_b, 2)))), function (x) findMRCA(aedini.tree, tips = as.character(data.frame(t(combn(clade_b, 2)))[x,]))))
min(sapply(1:nrow(data.frame(t(combn(psor, 2)))), function (x) findMRCA(aedini.tree, tips = as.character(data.frame(t(combn(psor, 2)))[x,]))))
max(sapply(1:nrow(data.frame(t(combn(psor, 2)))), function (x) findMRCA(aedini.tree, tips = as.character(data.frame(t(combn(psor, 2)))[x,]))))

#clade A - 150:234
#clade B - 235:288
#psorophora - 289:293

##Identify which tip labels belong to which clade
pdf("aediniTree.tipLabels.pdf", width=12, height=18)
plot.phylo(aedini.tree, show.tip.label = TRUE)
tiplabels(cex=0.6)
nodelabels(cex=0.5)
dev.off()

#clade A - 1:86
#clade B - 87:141
#psorophora - 142:147

x$edge_coords<-x$edge_coords %>%
  mutate(clade = case_when(node2 %in% 1:86 ~ 'Clade A',
                           node2 %in% 150:234 ~ 'Clade A',
                           node2 %in% 87:141 ~ 'Clade B',
                           node2 %in% 289:293 ~ 'Clade B',
                           node2 %in% 142:147 ~ 'Psorophora',
                           node2 %in% 283:293 ~'Psorophora',
                           .default = NA),
         clade = factor(clade, levels = c('Clade A', 'Clade B', 'Psorophora')))

##first find variance of dimensions
adult.mca$eig #dim1=8.436, dim2=5.67

#creating the plot
pdf("phylomorphospace_dim1_dim2.pdf", width=10, height=12)
(scale.setae.plot<-ggplot() + 
    geom_segment(data = x$edge_coords, aes(x = Dim.1.x, xend = Dim.1.y, y = Dim.2.x, yend = Dim.2.y, color = clade), lwd = 1.25, alpha = 0.5) + 
    geom_point(data = data.frame(adult.mca$ind$coord) %>% 
                 rownames_to_column('species') %>% 
                 left_join(all.dat.hab5.dim %>% 
                             select(species, hab, clade)), 
               aes(x =  Dim.1, y = Dim.2, fill = hab, shape = clade), size = 3) +
    geom_point(data = x$anc_dat, aes(x = Dim.1, y = Dim.2), size = 1, shape = 21, fill = 'white') +
    labs(x = 'Dim 1 (8.444% Variance)', y = 'Dim 2 (5.67% of Variance)') +
    ggtitle('Scale + Setae MCA') +
    scale_fill_manual(values = c("#00306FFF", "aquamarine4", "lightgreen", "#FFEA46FF", "azure3"),'Larval Habitat')+
    scale_shape_manual(values=c(21, 22, 23))+
    scale_color_manual(values = c("#60136EFF", "#CB4149FF", "black"), na.value = '#FCFFA4FF', labels = c(' Clade A', 'Clade B', 'Psorophora', 'Root'),
                       'Clade') +
    theme_bw() +
    theme(aspect.ratio = 1))
dev.off()

