##### Phylogenetic regression (phylolm) #####
# models 
mod0.d1<-phylolm(Dim.1 ~ 1, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod1.d1<-phylolm(Dim.1 ~ hab, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod2.d1<-phylolm(Dim.1 ~ pal+afr+ind+neo+nea+aus, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod3.d1<-phylolm(Dim.1 ~ oldWorld+newWorld, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod4.d1<-phylolm(Dim.1 ~ hab+pal+afr+ind+neo+nea+aus, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod5.d1<-phylolm(Dim.1 ~ hab+oldWorld+newWorld, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod6.d1<-phylolm(Dim.1 ~ hab*oldWorld+hab*newWorld, data = phylo.dat, phy = aedini.tree, model = 'lambda')

model.sel(mod0.d1,mod1.d1,mod2.d1,mod3.d1,mod4.d1,mod5.d1,mod6.d1) # mod1.d1 is the best model
mod1.d1 %>% summary
#Estimate    StdErr t.value   p.value    
#(Intercept)  -0.162647  0.119022 -1.3665    0.1739    
#habContainer  0.374532  0.065488  5.7191 5.904e-08 ***

# mod1.d1 posthoc
mod1.d1.rg<-qdrg(formula = mod1.d1$formula, data = phylo.dat, coef = mod1.d1$coefficients, vcov = mod1.d1$vcov, df = nrow(phylo.dat) - (length(mod1.d1$coefficients) + 1) - 1)
mod1.d1.emm<-data.frame(mod1.d1.rg, specs = 'hab')
mod1.d1.contrasts<-data.frame(pairs(mod1.d1.rg, adjust = 'fdr'))

(mod1.d1.boxplot<-ggplot(phylo.dat, aes(x = hab, y = Dim.1, fill = hab)) +
    geom_boxplot(width = 0.15, outliers = FALSE) +
    geom_signif(y_position = c(1.1), xmin = c('Open water'), xmax = c('Container'), annotation = c('P < 0.001'), tip_length = 0.01) +
    geom_jitter(data = phylo.dat, aes(x = hab, y = Dim.1, fill = hab), shape = 21, width = 0.1, alpha=0.5) +
    scale_fill_manual(values = c("#00306FFF", "#FFEA46FF"), 'Larval Habitat') +
    labs(x = 'Larval Habitat', y = 'MCA Dimension 1') +
    theme_bw() +
    coord_flip()+
    theme(aspect.ratio = 0.15, legend.position='none',
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank()))

pdf("mod1.d1_posthoc.pdf", width=6, height=5)
mod1.d1.boxplot
dev.off()

# dimensions 2
mod0.d2<-phylolm(Dim.2 ~ 1, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod1.d2<-phylolm(Dim.2 ~ hab, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod2.d2<-phylolm(Dim.2 ~ pal+afr+ind+neo+nea+aus, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod3.d2<-phylolm(Dim.2 ~ oldWorld+newWorld, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod4.d2<-phylolm(Dim.2 ~ hab+pal+afr+ind+neo+nea+aus, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod5.d2<-phylolm(Dim.2 ~ hab+oldWorld+newWorld, data = phylo.dat, phy = aedini.tree, model = 'lambda')

model.sel(mod0.d2,mod1.d2,mod2.d2,mod3.d2,mod4.d2,mod5.d2) # mod1.d2, mod0.d2, mod3.d2
mod1.d2 %>% summary
mod5.d2 %>% summary
mod0.d2 %>% summary

# dimensions 3
mod0.d3<-phylolm(Dim.3 ~ 1, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod1.d3<-phylolm(Dim.3 ~ hab, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod2.d3<-phylolm(Dim.3 ~ pal+afr+ind+neo+nea+aus, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod3.d3<-phylolm(Dim.3 ~ oldWorld+newWorld, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod4.d3<-phylolm(Dim.3 ~ hab+pal+afr+ind+neo+nea+aus, data = phylo.dat, phy = aedini.tree, model = 'lambda')
mod5.d3<-phylolm(Dim.3 ~ hab+oldWorld+newWorld, data = phylo.dat, phy = aedini.tree, model = 'lambda')

model.sel(mod0.d3,mod1.d3,mod2.d3,mod3.d3,mod4.d3,mod5.d3) #mod0.d3
mod0.d3 %>% summary
mod3.d3 %>% summary

# dataset for phylolm (exclude individuals without blood preference data)
phylo.dat.2 <- phylo.dat %>% drop_na() %>%
  rownames_to_column('species') #147 to 63 species
aedini.bf.tree<-geiger::treedata(mcmc.tree, data = phylo.dat.2 %>% column_to_rownames('species'), warnings = F, sort = T)$phy
phylo.dat2<- phylo.dat.2 %>% column_to_rownames("species")
str(phylo.dat2)

#models
m0.d1<-phylolm(Dim.1 ~ 1, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m1.d1<-phylolm(Dim.1 ~ hab, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m2.d1<-phylolm(Dim.1 ~ amp_prop+ave_prop+mam_prop+rep_prop, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m3.d1<-phylolm(Dim.1 ~ oldWorld+newWorld, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m4.d1<-phylolm(Dim.1 ~ pal+afr+ind+neo+nea+aus, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')

m5.d1<-phylolm(Dim.1 ~ hab+amp_prop+ave_prop+mam_prop+rep_prop, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m6.d1<-phylolm(Dim.1~ hab+oldWorld+newWorld, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m7.d1<-phylolm(Dim.1 ~ hab+pal+afr+ind+neo+nea+aus, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')

m8.d1<-phylolm(Dim.1 ~ hab+amp_prop+ave_prop+mam_prop+rep_prop+oldWorld+newWorld, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m9.d1<-phylolm(Dim.1 ~ hab+amp_prop+ave_prop+mam_prop+rep_prop+pal+afr+ind+neo+nea+aus, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m10.d1<-phylolm(Dim.1~ hab*oldWorld+hab*newWorld, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')

model.sel(m0.d1,m1.d1,m2.d1,m3.d1,m4.d1,m5.d1,m6.d1,m7.d1,m8.d1,m9.d1,m10.d1) #m1.d1, m6.d1
m1.d1 %>% summary
m6.d1 %>% summary

m0.d2<-phylolm(Dim.2 ~ 1, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m1.d2<-phylolm(Dim.2 ~ hab, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m2.d2<-phylolm(Dim.2 ~ amp_prop+ave_prop+mam_prop+rep_prop, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m3.d2<-phylolm(Dim.2 ~ oldWorld+newWorld, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m4.d2<-phylolm(Dim.2 ~ pal+afr+ind+neo+nea+aus, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')

m5.d2<-phylolm(Dim.2 ~ hab+amp_prop+ave_prop+mam_prop+rep_prop, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m6.d2<-phylolm(Dim.2~ hab+oldWorld+newWorld, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m7.d2<-phylolm(Dim.2 ~ hab+pal+afr+ind+neo+nea+aus, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')

m8.d2<-phylolm(Dim.2 ~ hab+amp_prop+ave_prop+mam_prop+rep_prop+oldWorld+newWorld, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m9.d2<-phylolm(Dim.2 ~ hab+amp_prop+ave_prop+mam_prop+rep_prop+pal+afr+ind+neo+nea+aus, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')

model.sel(m0.d2,m1.d2,m2.d2,m3.d2,m4.d2,m5.d2,m6.d2,m7.d2,m8.d2,m9.d2) #m8.d2, m5.d2
m8.d2 %>% summary
m5.d2 %>% summary

m0.d3<-phylolm(Dim.3 ~ 1, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m1.d3<-phylolm(Dim.3 ~ hab, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m2.d3<-phylolm(Dim.3 ~ amp_prop+ave_prop+mam_prop+rep_prop, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m3.d3<-phylolm(Dim.3 ~ oldWorld+newWorld, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m4.d3<-phylolm(Dim.3 ~ pal+afr+ind+neo+nea+aus, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')

m5.d3<-phylolm(Dim.3 ~ hab+amp_prop+ave_prop+mam_prop+rep_prop, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m6.d3<-phylolm(Dim.3~ hab+oldWorld+newWorld, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m7.d3<-phylolm(Dim.3 ~ hab+pal+afr+ind+neo+nea+aus, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')

m8.d3<-phylolm(Dim.3 ~ hab+amp_prop+ave_prop+mam_prop+rep_prop+oldWorld+newWorld, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')
m9.d3<-phylolm(Dim.3 ~ hab+amp_prop+ave_prop+mam_prop+rep_prop+pal+afr+ind+neo+nea+aus, data = phylo.dat2, phy = aedini.bf.tree, model = 'lambda')

model.sel(m0.d3,m1.d3,m2.d3,m3.d3,m4.d3,m5.d3,m6.d3,m7.d3,m8.d3,m9.d3) #m2.d3, m0.d3, m3.d3
m2.d3 %>% summary
m0.d3 %>% summary
m3.d3 %>% summary


##### Calculate phylogenetic signal #####
phylosig(aedini.tree, x = all.dat.hab2.dim %>% pull(Dim.1, name = species), method = 'K', test = T)
phylosig(aedini.tree, x = all.dat.hab5.dim %>% pull(Dim.1, name = species), method = 'K', test = T)

phylosig(aedini.tree, x = all.dat.hab2.dim %>% pull(Dim.1, name = species), method = 'lambda', test = T)
phylosig(aedini.tree, x = all.dat.hab2.dim %>% pull(Dim.1, name = species), method = 'lambda', test = T)
