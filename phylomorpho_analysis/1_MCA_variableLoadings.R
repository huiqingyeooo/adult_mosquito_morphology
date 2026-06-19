##### MCA #####
adult.mca <-adult.dat %>% MCA(X = ., ncp = length(scale.setae.col), graph = F)
larva.mca <-larva.dat %>% MCA(X = ., ncp = length(larva.col), graph = F)

# Generate matrix with dimension coordinates, species, habitats and clades
adult.mca.dim <- data.frame(adult.mca$ind$coord) %>%
  rownames_to_column('species') %>%
  left_join(all.dat.hab2 %>% select(species, hab, clade))

larva.mca.dim <- data.frame(larva.mca$ind$coord) %>%
  rownames_to_column('species') %>%
  left_join(all.dat.hab2 %>% select(species, hab, clade))

##### MCA scree plot ####
# Find number of dimensions that explains >= 1.5% of variance
length(adult.mca$eig[,2][adult.mca$eig[,2] >= 1.5]) #19
# Scree plot
data.frame(eigenvalue = adult.mca$eig[, 1]) %>%
  mutate(per_var = eigenvalue/sum(eigenvalue) * 100, cum_var = cumsum(per_var),dim = 1:nrow(.)) %>%
  ggplot(aes(x = dim, y = per_var, label = round(cum_var, 2))) +
  geom_vline(xintercept = length(adult.mca$eig[,2][adult.mca$eig[,2] >= 1.5]), lty = 2) +
  geom_col(fill = 'grey65', color = 'black') + geom_line(color = 'red') +geom_hline(yintercept = 1.5, lty = 2) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 25),expand = c(0, 0)) + scale_y_continuous(breaks = scales::pretty_breaks(n = 8),expand = c(0, 0.25)) +
  ylab('% Variance') + geom_text(nudge_y = 0.5, size = 2, angle = 45) +
  theme_void() +
  ggtitle('Scree plot for Scales + Setae') +
  theme(aspect.ratio = 0.75,axis.text.x = element_blank(),axis.title.x = element_blank(),axis.title.y = element_text(size = 9, angle = 90),axis.text.y = element_text(size= 8))

##### MCA biplots ######
# adult
# determine impt variables and visualise
# dimension 1
biplot.var <- get_mca_var(adult.mca) 
contributions <- data.frame(biplot.var$contrib)
threshold <- quantile(contributions[, "Dim.1"], 0.9) #threshold for top 10% biplot variables
top_vars <- rownames(contributions[contributions[, "Dim.1"] >= threshold, ])

(bp1<-fviz_mca_var(adult.mca, select.var = list(name = top_vars), 
                   geom = c("point", "text"), label = "var",
                   repel = TRUE) +
    ggtitle("Biplot of top 10% adult Dim 1 variables")+
    theme_bw() +
    theme(aspect.ratio = 1))

# dimension 2
threshold2 <- quantile(contributions[, "Dim.2"], 0.9) #threshold for top 10% biplot variables
top_vars2 <- rownames(contributions[contributions[, "Dim.2"] >= threshold2, ])
(bp2<-fviz_mca_var(adult.mca, select.var = list(name = top_vars2), 
                   geom = c("point", "text",title=""), label = "var",
                   repel = TRUE) +
    ggtitle("Biplot of top 10% Dim 2 variables")+
    theme_bw() +
    theme(aspect.ratio = 1))

# dimension 3
threshold3 <- quantile(contributions[, "Dim.3"], 0.9) #threshold for top 10% biplot variables
top_vars3 <- rownames(contributions[contributions[, "Dim.3"] >= threshold3, ])
(bp3<-fviz_mca_var(adult.mca, select.var = list(name = top_vars3), 
                   geom = c("point", "text"), label = "var",
                   repel = TRUE, axes = c(2, 3),arrows=TRUE) +
    ggtitle("Biplot of top 10% Dim 3 variables")+
    theme_bw() +
    theme(aspect.ratio = 1))

pdf("mca_adult_biplot_imptChr_dim1-3.pdf", width=15, height=6)
ggarrange(bp1, bp2, bp3,
          ncol=3, nrow=1, 
          labels=c("A","B","C"),
          common.legend = TRUE, legend="right")
dev.off()


# Generate and export biplot variable list
contributions[contributions[, "Dim.1"] >= threshold, ]$Dim.1
contributions[contributions[, "Dim.2"] >= threshold2, ]$Dim.2
contributions[contributions[, "Dim.3"] >= threshold3, ]$Dim.3

# match top 10% characters with dimension values
coords<-data.frame(biplot.var$coord)%>%rownames_to_column('characters')

dim.values1<- as.data.frame(top_vars)%>%
  rename(characters=top_vars)%>% #rename column top_vars to characters
  left_join(coords%>%select(characters,Dim.1))%>%
  mutate(values=Dim.1, dim="dim1")%>% #rename column Dim.1 to values
  select(characters,values,dim)

dim.values2<- as.data.frame(top_vars2)%>%
  rename(characters=top_vars2)%>%
  left_join(coords%>%select(characters,Dim.2))%>%
  mutate(values=Dim.2, dim="dim2")%>% #rename column Dim.1 to values
  select(characters,values,dim)

dim.values3<- as.data.frame(top_vars3)%>%
  rename(characters=top_vars3)%>%
  left_join(coords%>%select(characters,Dim.3))%>%
  mutate(values=Dim.3, dim="dim3")%>% #rename column Dim.1 to values
  select(characters,values,dim)

combined_list <- as.data.frame(rbind(dim.values1,dim.values2,dim.values3))
write.csv(combined_list, "mca_adult_biplot_imptChr_dim1-3.csv")

# larvae
# determine impt variables and visualise
# dimension 1
biplot.var <- get_mca_var(larva.mca) 
contributions <- data.frame(biplot.var$contrib)
threshold <- quantile(contributions[, "Dim.1"], 0.9) #threshold for top 10% biplot variables
top_vars <- rownames(contributions[contributions[, "Dim.1"] >= threshold, ])

(bp1<-fviz_mca_var(larva.mca, select.var = list(name = top_vars), 
                   geom = c("point", "text"), label = "var",
                   repel = TRUE) +
    ggtitle("Biplot of top 10% adult Dim 1 variables")+
    theme_bw() +
    theme(aspect.ratio = 1))

# dimension 2
threshold2 <- quantile(contributions[, "Dim.2"], 0.9) #threshold for top 10% biplot variables
top_vars2 <- rownames(contributions[contributions[, "Dim.2"] >= threshold2, ])
(bp2<-fviz_mca_var(larva.mca, select.var = list(name = top_vars2), 
                   geom = c("point", "text",title=""), label = "var",
                   repel = TRUE) +
    ggtitle("Biplot of top 10% Dim 2 variables")+
    theme_bw() +
    theme(aspect.ratio = 1))

# dimension 3
threshold3 <- quantile(contributions[, "Dim.3"], 0.9) #threshold for top 10% biplot variables
top_vars3 <- rownames(contributions[contributions[, "Dim.3"] >= threshold3, ])
(bp3<-fviz_mca_var(larva.mca, select.var = list(name = top_vars3), 
                   geom = c("point", "text"), label = "var",
                   repel = TRUE, axes = c(2, 3),arrows=TRUE) +
    ggtitle("Biplot of top 10% Dim 3 variables")+
    theme_bw() +
    theme(aspect.ratio = 1))

pdf("mca_larva_biplot_imptChr_dim1-3.pdf", width=15, height=6)
ggarrange(bp1, bp2, bp3,
          ncol=3, nrow=1, 
          labels=c("A","B","C"),
          common.legend = TRUE, legend="right")
dev.off()


# Generate and export biplot variable list
contributions[contributions[, "Dim.1"] >= threshold, ]$Dim.1
contributions[contributions[, "Dim.2"] >= threshold2, ]$Dim.2
contributions[contributions[, "Dim.3"] >= threshold3, ]$Dim.3

# match top 10% characters with dimension values
coords<-data.frame(biplot.var$coord)%>%rownames_to_column('characters')

dim.values1<- as.data.frame(top_vars)%>%
  rename(characters=top_vars)%>% #rename column top_vars to characters
  left_join(coords%>%select(characters,Dim.1))%>%
  mutate(values=Dim.1, dim="dim1")%>% #rename column Dim.1 to values
  select(characters,values,dim)

dim.values2<- as.data.frame(top_vars2)%>%
  rename(characters=top_vars2)%>%
  left_join(coords%>%select(characters,Dim.2))%>%
  mutate(values=Dim.2, dim="dim2")%>% #rename column Dim.1 to values
  select(characters,values,dim)

dim.values3<- as.data.frame(top_vars3)%>%
  rename(characters=top_vars3)%>%
  left_join(coords%>%select(characters,Dim.3))%>%
  mutate(values=Dim.3, dim="dim3")%>% #rename column Dim.1 to values
  select(characters,values,dim)

combined_list <- as.data.frame(rbind(dim.values1,dim.values2,dim.values3))
write.csv(combined_list, "mca_larva_biplot_imptChr_dim1-3.csv")


##### Finding important variables based on high R2 #####
adult.dimdesc<-bind_rows(lapply(1:19, function(x) data.frame(dimdesc(adult.mca, axes = x)[[1]]$quali) %>% rownames_to_column('var')), .id = 'dim') %>%
    mutate(var = factor(var),dim = factor(dim, levels = 1:19))

# heatmap
ggplot() +
    geom_tile(data = adult.dimdesc, aes(x = dim, y= var, fill = R2)) +
    geom_text(data = adult.dimdesc %>% filter(R2 > 0.3), aes(x = dim, y = var), label = '*') +
    ggtitle('Heat map of Variable importance for Scale+Setae Variables') +
    scale_fill_viridis(option = 'viridis') +
    labs(x = 'Dimension', y = 'Variable') +
    theme_minimal() +
    theme(aspect.ratio = 1.75)

# Finds and returns list of variables and characters with “high” cutoff of R2 ≥0.3
# axes refers to dimensions here
scale.setae.var1<-data.frame(dimdesc(adult.mca, axes = 1)[[1]]$quali) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  filter(R2 >= 0.3) %>% #all characters with R2>=0.3 have significant p-values
  rownames_to_column('var') %>%
  pull(var)

scale.setae.cat1<-data.frame(dimdesc(adult.mca, axes = 1)[[1]]$category) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  rownames_to_column('var_cat') %>%
  separate(var_cat, into = c('var', 'cat'), sep = '=') %>%
  filter(var %in% scale.setae.var1) %>%
  pull(cat)

scale.setae.var2<-data.frame(dimdesc(adult.mca, axes = 2)[[1]]$quali) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  filter(R2 >= 0.3) %>%
  rownames_to_column('var') %>%
  pull(var)

scale.setae.cat2<-data.frame(dimdesc(adult.mca, axes = 2)[[1]]$category) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  rownames_to_column('var_cat') %>%
  separate(var_cat, into = c('var', 'cat'), sep = '=') %>%
  filter(var %in% scale.setae.var2) %>%
  pull(cat)

scale.setae.var3<-data.frame(dimdesc(adult.mca, axes = 3)[[1]]$quali) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  filter(R2 >= 0.3) %>%
  rownames_to_column('var') %>%
  pull(var)

scale.setae.cat3<-data.frame(dimdesc(adult.mca, axes = 3)[[1]]$category) %>%
  mutate(p.adj = p.adjust(p.value, method = 'fdr')) %>%
  rownames_to_column('var_cat') %>%
  separate(var_cat, into = c('var', 'cat'), sep = '=') %>%
  filter(var %in% scale.setae.var2) %>%
  pull(cat)

# Generate a matrix of highly correlated variables with the column imp.dim
all.scale.setae.var<-data.frame(adult.mca$var$coord)
important.scale.setae.var<-all.scale.setae.var %>%
  rownames_to_column('cat')%>%
  filter(cat %in% scale.setae.cat1 | cat %in% scale.setae.cat2 | cat %in% scale.setae.cat3)%>% # filters by list of variables impt in dimensions 1 to 3
  select(cat:Dim.3)%>% # selects the first 3 dimensions of the matrix
  # create a column 'imp.dim' and label rows with the dimension numbers
  mutate(imp.dim = case_when(cat %in% scale.setae.cat1 & cat %in% scale.setae.cat2 & !(cat %in% scale.setae.cat3) ~ 'Dim.1&2',
                             cat %in% scale.setae.cat1 & !(cat %in% scale.setae.cat2) & !(cat %in% scale.setae.cat3) ~ 'Dim.1',
                             !(cat %in% scale.setae.cat1) & cat %in% scale.setae.cat2 & !(cat %in% scale.setae.cat3) ~ 'Dim.2',
                             !(cat %in% scale.setae.cat1) & !(cat %in% scale.setae.cat2) & cat %in% scale.setae.cat3 ~ 'Dim.3',
                             cat %in% scale.setae.cat1 & !(cat %in% scale.setae.cat2) & cat %in% scale.setae.cat3 ~ 'Dim.1&3',
                             !(cat %in% scale.setae.cat1) & cat %in% scale.setae.cat2 & cat %in% scale.setae.cat3 ~ 'Dim.2&3',
                             cat %in% scale.setae.cat1 & cat %in% scale.setae.cat2 & cat %in% scale.setae.cat3 ~ 'Dim.1&2&3'),
         imp.dim = factor(imp.dim, levels = c('Dim.1', 'Dim.2', 'Dim.3','Dim.1&2','Dim.1&3','Dim.2&3','Dim.1&2&3'))
  )

# Plot characters for dimensions 1,2,and 3
p1<-ggplot() +
  geom_point(data = all.scale.setae.var, aes(x = Dim.1, y = Dim.2), shape = 1, size = 3) +
  geom_point(data = important.scale.setae.var, aes(x = Dim.1, y = Dim.2, fill = imp.dim), shape = 21, size = 3) +
  scale_fill_manual(values = c("#482173FF", "#2D708EFF", "#85D54AFF",'grey',"#FDE725FF")) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  ggtitle("Scales + Setae: Dim.1 vs Dim.2")

p2<-ggplot() +
  geom_point(data = important.scale.setae.var, aes(x = Dim.1, y = Dim.2, fill = imp.dim), shape = 21, size = 3) +
  geom_text_repel(data = important.scale.setae.var, aes(x = Dim.1, y = Dim.2, label = cat, color = imp.dim), size = 2.5) +
  scale_fill_manual(values = c("#482173FF", "#2D708EFF", "#85D54AFF",'grey',"#FDE725FF"))  +
  scale_color_manual(values = c("#482173FF", "#2D708EFF", "#85D54AFF",'grey',"#FDE725FF"))  +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = 'none')+
  ggtitle("Scales + Setae: Dim.1 vs Dim.2 \n(Important characters)")

p3<-ggplot() +
  geom_point(data = all.scale.setae.var, aes(x = Dim.1, y = Dim.3), shape = 1, size = 3) +
  geom_point(data = important.scale.setae.var, aes(x = Dim.1, y = Dim.3, fill = imp.dim), shape = 21, size = 3) +
  scale_fill_manual(values = c("#482173FF", "#2D708EFF", "#85D54AFF",'grey',"#FDE725FF")) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  ggtitle("Scales + Setae: Dim.1 vs Dim.3")

p4<-ggplot() +
  geom_point(data = important.scale.setae.var, aes(x = Dim.1, y = Dim.3, fill = imp.dim), shape = 21, size = 3) +
  geom_text_repel(data = important.scale.setae.var, aes(x = Dim.1, y = Dim.3, label = cat, color = imp.dim), size = 2.5) +
  scale_fill_manual(values = c("#482173FF", "#2D708EFF", "#85D54AFF",'grey',"#FDE725FF"))  +
  scale_color_manual(values = c("#482173FF", "#2D708EFF", "#85D54AFF",'grey',"#FDE725FF"))  +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = 'none')+
  ggtitle("Scales + Setae: Dim.1 vs Dim.3 \n(Important characters)")

p5<-ggplot() +
  geom_point(data = all.scale.setae.var, aes(x = Dim.2, y = Dim.3), shape = 1, size = 3) +
  geom_point(data = important.scale.setae.var, aes(x = Dim.2, y = Dim.3, fill = imp.dim), shape = 21, size = 3) +
  scale_fill_manual(values = c("#482173FF", "#2D708EFF", "#85D54AFF",'grey',"#FDE725FF")) +
  theme_bw() +
  theme(aspect.ratio = 1)+
  ggtitle("Scales + Setae: Dim.2 vs Dim.3")

p6<-ggplot() +
  geom_point(data = important.scale.setae.var, aes(x = Dim.2, y = Dim.3, fill = imp.dim), shape = 21, size = 3) +
  geom_text_repel(data = important.scale.setae.var, aes(x = Dim.2, y = Dim.3, label = cat, color = imp.dim), size = 2.5) +
  scale_fill_manual(values = c("#482173FF", "#2D708EFF", "#85D54AFF",'grey',"#FDE725FF"))  +
  scale_color_manual(values = c("#482173FF", "#2D708EFF", "#85D54AFF",'grey',"#FDE725FF"))  +
  theme_bw() +
  theme(aspect.ratio = 1,
        legend.position = 'none')+
  ggtitle("Scales + Setae: Dim.2 vs Dim.3 \n(Important characters)")

p1 #dim1 vs dim2
p2 #dim1 vs dim2 (impt characters)
p3 #dim1 vs dim3
p4 #dim1 vs dim3 (impt characters)
p5 #dim2 vs dim3
p6 #dim2 vs dim3 (impt characters)

pdf("mca_adult_R2_imptChr_dim1-3.pdf", width=15, height=10)
ggarrange(p1,p3,p5,p2,p4,p6, ncol=3, nrow=2, 
          labels=c("A","B","C","D","E","F"),
          common.legend = TRUE, legend="right")
dev.off()

print.data.frame(important.scale.setae.var)
write.csv(important.scale.setae.var, "mca_adult_R2_imptChr_dim1-3.csv")

  
