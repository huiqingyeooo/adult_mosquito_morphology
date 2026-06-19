##### MCA visualization #####
##### scatterplot #####
# scatterplot for adults
ggplot(adult.mca.dim)+
  geom_point(aes(x=Dim.1, y=Dim.2, shape=clade, fill= hab), size = 3, alpha=0.7)+
  scale_fill_manual(values = c("#00306FFF","#FFEA46FF"), 'Larval Habitat') +
  scale_shape_manual(values=c(21, 22, 23))+
  xlab(paste0("Dimension 1 (", round(adult.mca$eig[1,2],2),"%)")) +
  ylab(paste0("Dimension 2 (", round(adult.mca$eig[2,2],2),"%)")) +
  guides(fill = guide_legend(override.aes = list(colour = c("#00306FFF","#FFEA46FF")))) +
  ggtitle("Adult MCA")+
  theme_bw() +
  theme(aspect.ratio = 1)

### Visualizing adult scatterplots ###
# scatterplot for larvae
ggplot(larva.mca.dim)+
  geom_point(aes(x=Dim.1, y=Dim.2, shape=clade, fill= hab), size = 3, alpha=0.7)+
  scale_fill_manual(values = c("#00306FFF","#FFEA46FF"), 'Larval Habitat') +
  scale_shape_manual(values=c(21, 22, 23))+
  xlab(paste0("Dimension 1 (", round(larva.mca$eig[1,2],2),"%)")) +
  ylab(paste0("Dimension 2 (", round(larva.mca$eig[2,2],2),"%)")) +
  guides(fill = guide_legend(override.aes = list(colour = c("#00306FFF","#FFEA46FF")))) +
  ggtitle("Larval MCA")+
  theme_bw() +
  theme(aspect.ratio = 1)


##### plots to visualize points with 01,02 and 23 categories #####
#dim 1 vs dim 2 
(mca1.raw<- ggplot(all.dat.hab.raw.dim)+
   geom_point(aes(x=Dim.1, y=Dim.2, shape=clade, fill= hab), size = 3, alpha=0.7)+
   scale_fill_manual(values = c("#00306FFF", "aquamarine4", "lightgreen", "#FFEA46FF", "azure3","red","purple","cornflowerblue"), 'Larval Habitat') +
   scale_shape_manual(values=c(21, 22, 23))+
   scale_color_manual("black")+
   labs(x = 'Dim 1 (8.46% Variance)', y = 'Dim 2 (5.65% of Variance)') +
   guides(fill = guide_legend(override.aes = list(colour = c("#00306FFF", "aquamarine4", "lightgreen", "#FFEA46FF", "azure3","red","purple","cornflowerblue")))) +
   theme_bw() +
   theme(aspect.ratio = 1))
#dim 1 vs dim 3
(mca2.raw <- ggplot(all.dat.hab.raw.dim)+
    geom_point(aes(x=Dim.1, y=Dim.3, shape=clade, fill= hab), size = 3, alpha=0.7)+
    scale_fill_manual(values = c("#00306FFF", "aquamarine4", "lightgreen", "#FFEA46FF", "azure3","red","purple","cornflowerblue"), 'Larval Habitat') +
    scale_shape_manual(values=c(21, 22, 23))+
    scale_color_manual("black")+
    labs(x = 'Dim 1 (8.46% Variance)', y = 'Dim 3 (4.17% of Variance)') +
    guides(fill = guide_legend(override.aes = list(colour = c("#00306FFF", "aquamarine4", "lightgreen", "#FFEA46FF", "azure3","red","purple","cornflowerblue")))) +
    theme_bw() +
    theme(aspect.ratio = 1))

#dim 2 vs dim 3
(mca3.raw <- ggplot(all.dat.hab.raw.dim)+
    geom_point(aes(x=Dim.2, y=Dim.3, shape=clade, fill= hab), size = 3, alpha=0.7)+
    scale_fill_manual(values = c("#00306FFF", "aquamarine4", "lightgreen", "#FFEA46FF", "azure3","red","purple","cornflowerblue"), 'Larval Habitat') +
    scale_shape_manual(values=c(21, 22, 23))+
    scale_color_manual("black")+
    labs(x = 'Dim 2 (5.65% Variance)', y = 'Dim 3 (4.17% of Variance)') +
    guides(fill = guide_legend(override.aes = list(colour = c("#00306FFF", "aquamarine4", "lightgreen", "#FFEA46FF", "azure3","red","purple","cornflowerblue")))) +
    theme_bw() +
    theme(aspect.ratio = 1))

pdf("mca_habRaw_dim1-3.pdf", width=15, height=5)
ggarrange(mca1.raw,mca2.raw,mca3.raw,
          ncol=3, nrow=1, 
          labels=c("A","B","C"),
          common.legend = TRUE, legend="right")
dev.off()

# plot with 5 habitats and species labels for comparison 
(mca1.spp <- ggplot(all.dat.hab5.dim)+
    geom_point(aes(x=Dim.1, y=Dim.2, shape=clade, fill=clade), size = 3, alpha=0.7)+
    scale_fill_manual(values = c("#60136EFF", "#CB4149FF", "black"), na.value = '#FCFFA4FF',labels= c('Clade A', 'Clade B', 'Psorophora', 'Root')) +
    scale_shape_manual(values=c(21, 22, 23))+
    labs(x = 'Dim 1 (8.46% Variance)', y = 'Dim 3 (4.17% of Variance)') +
    geom_text_repel(aes(x = Dim.1, y = Dim.3, label = species), size = 3) +
    theme_bw() +
    theme(aspect.ratio = 1))

(mca2.spp <- ggplot(all.dat.hab5.dim)+
    geom_point(aes(x=Dim.1, y=Dim.3, shape=clade, fill=clade), size = 3, alpha=0.7)+
    scale_fill_manual(values = c("#60136EFF", "#CB4149FF", "black"), na.value = '#FCFFA4FF',labels= c('Clade A', 'Clade B', 'Psorophora', 'Root')) +
    scale_shape_manual(values=c(21, 22, 23))+
    labs(x = 'Dim 1 (8.46% Variance)', y = 'Dim 3 (4.17% of Variance)') +
    geom_text_repel(aes(x = Dim.1, y = Dim.3, label = species), size = 3) +
    theme_bw() +
    theme(aspect.ratio = 1))

(mca3.spp <- ggplot(all.dat.hab5.dim)+
    geom_point(aes(x=Dim.2, y=Dim.3, shape=clade, fill=clade), size = 3, alpha=0.7)+
    scale_fill_manual(values = c("#60136EFF", "#CB4149FF", "black"), na.value = '#FCFFA4FF',labels= c('Clade A', 'Clade B', 'Psorophora', 'Root')) +
    scale_shape_manual(values=c(21, 22, 23))+
    labs(x = 'Dim 2 (5.65% Variance)', y = 'Dim 3 (4.17% of Variance)') +
    geom_text_repel(aes(x = Dim.2, y = Dim.3, label = species), size = 3) +
    theme_bw() +
    theme(aspect.ratio = 1))

# for figure with posthoc plot
# mca plot
mca1.raw

#boxplot (posthoc)
(mod1.d1.boxplot<-ggplot(phylo.dat, aes(x = hab, y = Dim.1, fill = hab)) +
    geom_boxplot(width = 0.15, outliers = FALSE) +
    geom_signif(y_position = c(1.1), xmin = c('Open water'), xmax = c('Container'))+
                #annotation = c('P < 0.001'), tip_length = 0.01) +
    geom_jitter(data = phylo.dat, aes(x = hab, y = Dim.1, fill = hab), shape = 21, width = 0.1, alpha=0.5) +
    scale_fill_manual(values = c("#00306FFF", "#FFEA46FF"), 'Larval Habitat') +
    labs(x = 'Larval Habitat', y = 'MCA Dimension 1') +
    theme_bw() +
    coord_flip()+
    theme(aspect.ratio = 0.1, legend.position='none',
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()))

top.good <- mod1.d1.boxplot+mca1.raw+plot_layout(nrow = 2, ncol = 1, height = c(1,7))
ggsave2(filename = 'topgood.png', plot = top.good, width = 9, height = 9)
