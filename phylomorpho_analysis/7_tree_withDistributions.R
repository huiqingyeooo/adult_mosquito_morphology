######## Plot tree with habitat + distribution #####
library(treeio)
library(ggtree)
library(ggplot2)

aedini.tree2<-geiger::treedata(aedini.tree, data = phylo.dat, warnings = F, sort = T)$phy

distribution<-as.data.frame(phylo.dat %>%
  select(neo,nea,pal,ind,aus,afr)%>%
  mutate(neo = str_replace_all(neo, "1", "2")) %>%
  mutate(pal = str_replace_all(pal, "1", "3")) %>%
  mutate(ind = str_replace_all(ind, "1", "4")) %>%
  mutate(aus = str_replace_all(aus, "1", "5")) %>%
  mutate(afr = str_replace_all(afr, "1", "6")) %>%
  mutate(across(c("pal","afr","ind","neo","nea","aus"), as.factor)))
str(distribution)

# plot tree with distribution heatmap
p<-ggtree(aedini.tree)
p
gheatmap(p,distribution, width=.4, offset=22, colnames=F)
#  scale_fill_manual(values=c("#ffffff","#EE6C4D","#FF9F1C","#293241","#3D5A80","#98C1D9","#84D2F6"))+
#  geom_tiplab(size=0.9)+
#  scale_x_ggtree()+
#  xlab("Time")+
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

# Group tree based on clades
clades<-list(clade_abp = c(clade_a, clade_b, psor),clade_a = clade_a, clade_b = clade_b, psorophora = psor)
aedini.tree.clades<-groupOTU(aedini.tree, clades)

# Plot tree with habitats + distribution heatmap
(tree.plot<-ggtree(aedini.tree.clades, aes(color = group), lwd = 0.7) %<+% 
    all.dat.hab5 +
    geom_tippoint(aes(fill = hab), shape = 21, size = 3, color = 'black')+
    #scale_color_manual(values = c("#60136EFF",'#FCFFA4FF', "#CB4149FF", "black"), 'Clade',
    scale_color_manual(values = c("#60136EFF",'black', "#CB4149FF", "black"), 'Clade',
                       labels = c('clade_abp' = 'Root', 'clade_a' = 'Clade A', 'clade_b' = 'Clade B',
                                  'psorophora' = 'Psorophora'))+
    theme(aspect.ratio = 1.7))

#pdf("tree_combined.pdf", width=12, height=18)
gheatmap(tree.plot,distribution, width=.4, offset=45, colnames=F)+
  scale_fill_manual(values = c("#00306FFF", "aquamarine4", "#2EC4B6", "#FFEA46FF", "azure3", #colours for larval habitat
                               "#F5F5F5","#EE6C4D","#FF9F1C","#293241","#3D5A80","#98C1D9","#84D2F6" #colours for distribution
                               ), 'Larval Habitat') +
  geom_tiplab(size=2.5, offset=1.2, color="black")+
  scale_x_ggtree()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))
#dev.off()
