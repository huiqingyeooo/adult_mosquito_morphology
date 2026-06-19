##### larval dimension 1 vs adult dimension 1 ####
library(ggtree)
library(ggtreeExtra)
library(ggplot2)

# MCA
mca.dim1<-data.frame(adult.mca$ind$coord) %>% select(Dim.1) %>% 
  rownames_to_column('species')%>%
  rename(dim.1.adult = Dim.1) %>%
  left_join(
    data.frame(larva.mca$ind$coord) %>% 
      select(Dim.1) %>%
      rename(dim.1.larva = Dim.1) %>%
      rownames_to_column('species'), by="species") %>%
  left_join(all.dat.hab2 %>% select(species,hab), by="species") %>%
  mutate(hab = as.factor(hab)) %>%
  rename(hab2 = hab) %>%
  left_join(all.dat.hab5 %>% select(species,hab,clade), by="species") %>%
  mutate(hab = as.factor(hab)) %>%
  rename(hab5 = hab)

nrow(adult.dat); nrow(larva.dat);nrow(mca.dim1)

##### MCA biplots ###### 
# determine impt variables and visualise
# dimension 1
adult.biplot <- get_mca_var(adult.mca) 
adult.contributions <- data.frame(adult.biplot$contrib)
threshold <- quantile(adult.contributions[, "Dim.1"], 0.9) #threshold for top 10% biplot variables
top_vars <- rownames(adult.contributions[adult.contributions[, "Dim.1"] >= threshold, ])

(bp1<-fviz_mca_var(adult.mca, select.var = list(name = top_vars), 
                   geom = c("point", "text"), label = "var",
                   repel = TRUE) +
    ggtitle("Biplot of top 10% adult Dim 1 variables")+
    theme_bw() +
    theme(aspect.ratio = 1))

larva.biplot <- get_mca_var(larva.mca)
larva.contributions <- data.frame(larva.biplot$contrib)
threshold <- quantile(larva.contributions[, "Dim.1"], 0.9) #threshold for top 10% biplot variables
top_vars <- rownames(larva.contributions[larva.contributions[, "Dim.1"] >= threshold, ])

(bp2<-fviz_mca_var(larva.mca, select.var = list(name = top_vars), 
                   geom = c("point", "text"), label = "var",
                   repel = TRUE) +
    ggtitle("Biplot of top 10% larval Dim 1 variables")+
    theme_bw() +
    theme(aspect.ratio = 1))


# adult dim 1 vs larval dim 1. two habitats
ggplot(mca.dim1, aes(x=dim.1.adult, y=dim.1.larva,color=hab2,fill=hab2))+
  geom_point(aes(shape=clade),size=3, alpha=0.9, color='black')+
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("#00306FFF","#FFEA46FF"), 'Larval Habitat') +
  scale_fill_manual(values = c("#00306FFF","#FFEA46FF"), 'Larval Habitat') +
  scale_shape_manual(values=c(21, 22, 23))+
  theme_bw()+
  theme(aspect.ratio = 1)

# adult dim 1 vs larval dim 1. five habitats. with species labels.
ggplot(mca.dim1, aes(x=dim.1.adult, y=dim.1.larva,color=hab5, fill=hab5, label=species))+
  geom_point(size=3, alpha=0.9, pch=21, color='black')+
  geom_text_repel(size=1.5,color="black")+
  #geom_smooth(method = "lm")+
  scale_fill_manual(values = c("#00306FFF", "aquamarine4", "lightgreen", "#FFEA46FF", "azure3"), 'Larval Habitat') +
  theme_bw()
  theme(aspect.ratio = 1)


##### phylolm #####
# dataset for phylolm
phylo.dat <- mca.dim1 %>%
  column_to_rownames('species') 
#%>%
#  filter(!grepl("Brackish", hab5)) %>%
#  filter(!grepl("Freshwater rock", hab5)) %>%
# filter(!grepl("Crabholes", hab5))

# ensure tips in tree matches species in dataset
aedini.tree<-geiger::treedata(mcmc.tree, data =phylo.dat, warnings = T, sort = T)$phy
aedini.tree<-drop.tip(aedini.tree, tip = out_group)

# phylolm models
mod0 <- phylolm(dim.1.larva ~ 1, data = phylo.dat, phy = aedini.tree, model='lambda')
mod1 <- phylolm(dim.1.larva ~ dim.1.adult, data = phylo.dat, phy = aedini.tree, model='lambda')
mod2 <- phylolm(dim.1.larva ~ dim.1.adult+hab2, data = phylo.dat, phy = aedini.tree, model='lambda')
mod3 <- phylolm(dim.1.larva ~ dim.1.adult*hab2, data = phylo.dat, phy = aedini.tree, model='lambda')
model.sel(mod0, mod1, mod2, mod3)
summary(mod2)

model0 <- phylolm(dim.1.adult ~ 1, data = phylo.dat, phy = aedini.tree, model='lambda')
model1 <- phylolm(dim.1.adult ~ dim.1.larva, data = phylo.dat, phy = aedini.tree, model='lambda')
model2 <- phylolm(dim.1.adult ~ dim.1.larva+hab2, data = phylo.dat, phy = aedini.tree, model='lambda')
model3 <- phylolm(dim.1.adult ~ dim.1.larva*hab2, data = phylo.dat, phy = aedini.tree, model='lambda')
model4 <- phylolm(dim.1.adult ~ hab2, data = phylo.dat, phy = aedini.tree, model='lambda')
model.sel(model0, model1, model2, model3, model4)

summary(model2)
summary(model3)

# getting r2 
library(phylolm)
library(rr2)
m_full    <- phylolm(dim.1.adult ~ dim.1.larva + hab2, data = phylo.dat, phy = aedini.tree, model = "lambda")
m_larva   <- phylolm(dim.1.adult ~ dim.1.larva, data = phylo.dat, phy = aedini.tree, model = "lambda")
m_habitat <- phylolm(dim.1.adult ~ hab2, data = phylo.dat, phy = aedini.tree, model = "lambda")

# Total R² (vs intercept-only)
R2(mod = m_full, phy=aedini.tree)

# Partial R² — unique contribution of each predictor
R2(mod = m_full, mod.r = m_habitat, phy=aedini.tree)   # unique to larva
R2(mod = m_full, mod.r = m_larva, phy=aedini.tree)     # unique to habitat
  
### make plot
# Group tree based on clades
clades<-list(clade_abp = c(clade_a, clade_b, psor),clade_a = clade_a, clade_b = clade_b, psorophora = psor)
aedini.tree.clades<-groupOTU(aedini.tree, clades)

# visualize the tree 
(p <- ggtree(aedini.tree.clades, aes(color=group),lwd=0.7)%<+% 
  phylo.dat + 
  geom_tiplab(size=3, color="black", offset=5) + xlim_tree(500)+
  scale_color_manual(values = c("#60136EFF",'black', "#CB4149FF", "black"), 'Clade',
                   labels = c('clade_abp' = 'Root', 'clade_a' = 'Clade A', 'clade_b' = 'Clade B',
                              'psorophora' = 'Psorophora'))+
  theme(aspect.ratio = 1.7))

pdf("larval_vs_adult_dim1.pdf", width=12, height=18)
p + geom_fruit(data=mca.dim1, geom=geom_bar,
               mapping = aes(y=species,x=mca.dim1$dim.1.larva, fill=mca.dim1$hab5),
               stat="identity", orientation="y",width=0.8, offset=2.8, pwidth=0.7, colour="black",linewidth=0.1) +
  geom_fruit(data=mca.dim1, geom=geom_bar,
             mapping = aes(y = species,x = mca.dim1$dim.1.adult,fill = mca.dim1$hab5),
             stat = "identity", orientation = "y",width=0.8, offset=0.9, pwidth=0.7,colour="black",linewidth=0.1) +
  scale_fill_manual(values = c("#00306FFF", "aquamarine4", "lightgreen", "#FFEA46FF", "azure3"), 'Larval habitat')+
  ggplot2::annotate("text", x = 430, y = 150, label = "Adult dimension 1", fontface = "bold", size=5)+
  ggplot2::annotate("text", x = 300, y = 150, label = "Larva dimension 1", fontface = "bold", size=5)+
  ggplot2::annotate("text", x = 180, y = 151, label = " ", fontface = "bold")+
  guides(fill = guide_legend(override.aes = list(size = 0.8))) +
  theme_tree(legend.text = element_text(size = 8))
dev.off()

# visualize tree with heatmap and without tip labels
(p.small <- ggtree(aedini.tree.clades, aes(color=group),lwd=0.5)%<+% 
    phylo.dat + 
    #geom_tiplab(size=3, color="black", offset=5)+ 
    #xlim_tree(500)+
    scale_color_manual(values = c("#60136EFF",'black', "#CB4149FF", "black"), 'Clade',
                       labels = c('clade_abp' = 'Root', 'clade_a' = 'Clade A', 'clade_b' = 'Clade B',
                                  'psorophora' = 'Psorophora'))+
    theme(aspect.ratio = 1))

pdf("larval_vs_adult_dim1_woLabels.pdf", width=8, height=12)
p.small + geom_fruit(data=mca.dim1, geom=geom_bar,
               mapping = aes(y=species,x=mca.dim1$dim.1.larva, fill=mca.dim1$hab5),
               stat="identity", orientation="y",width=1, offset=0.6, pwidth=0.6, colour="black",linewidth=0.1,
               axis.params = list(axis = "x",text.size = 2)) +
  geom_fruit(data=mca.dim1, geom=geom_bar,
             mapping = aes(y = species,x = mca.dim1$dim.1.adult,fill = mca.dim1$hab5),
             stat = "identity", orientation = "y",width=1, offset=0.6, pwidth=0.6,colour="black",linewidth=0.1,
             axis.params = list(axis = "x",text.size = 2)) +
  scale_fill_manual(values = c("#00306FFF", "aquamarine4", "lightgreen", "#FFEA46FF", "azure3"), 'Larval habitat')+
  #ggplot2::annotate("text", x = 225, y = 150, label = "Adult dimension 1", fontface = "bold", size=3)+
  #ggplot2::annotate("text", x = 125, y = 150, label = "Larva dimension 1", fontface = "bold", size=3)+
  guides(fill = guide_legend(override.aes = list(size = 0.8))) +
  theme_tree(legend.text = element_text(size = 8))+
  theme(aspect.ratio = 1.2)
dev.off()
