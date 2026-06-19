##### Load libraries and setwd #####
library(tidyverse)
library(FactoMineR)
library(factoextra)
library(phytools)
library(ggtree)
library(clipr)
library(geiger)
library(viridis)
library(cowplot)
library(patchwork)
library(ggrepel)
library(ggsignif)
library(emmeans)
library(phylolm)
library(ggpubr)
library(MuMIn)
rm(list=ls())
setwd("/path/to/working/dir")

##### Load datasets #####
larv.dat<-read.csv('reinert_aedini_larvae_char.csv', header = T)
adult.dat<-read.csv('reinert_aedini_adult_char.csv', header = T)
mcmc.tree<-read.tree('S2017S2023genbank_subgen.chronos.RHK.treefile')
#hostpref.biogeog.dat<-read.csv('Aedini_hostPref_biogeo_updated.csv')

# joins larval with adult characters (270 rows)
all.dat<-larv.dat %>%
  select(-c(X336)) %>% #remove larval habitats 336 as it is present in adult.dat
  left_join(adult.dat, by="species") %>%
  mutate(across(X5:X336, ~ factor(gsub('(\\[)([0-9]+)(\\])', '\\2', .x))), #removes [ ] in dataset
         X162 = case_when(X162 == '{' ~ '-', .default = X162),
         X336 = as.factor(case_when(species == 'Aedes_Cornetius_cozi' ~ '3', #update larval habitats
                                    species == 'Armigeres_Leicesteria_longipalpis' ~ '3',
                                    .default = X336)))

# removes rows not represented in tree (153 rows)
all.dat2<-data.frame(geiger::treedata(mcmc.tree, data = all.dat %>% column_to_rownames('species'), warnings = T, sort = T)$data) %>%
  mutate(across(X5:X336, ~ as.factor(.x)))%>%
  rownames_to_column('species')

# exclude species without larval habitat data
all.dat2<-all.dat2 %>% filter(!(species %in% c('Aedes_OchlerotatusChrysoconops_fulvus', 
                                               'Aedes_OchlerotatusChrysoconops_pallens',
                                               'Aedes_Aedimorphus_eritreae')))
# exclude outgroups
out_group<-c('Culex_Culex_quinquefasciatus', 'Orthopodomyia_signifera', 'Culiseta_inornata', 'Mansonia_titillans')
all.dat2<-all.dat2 %>% filter(!(species %in% out_group))

##### Adult characters #####
# Determine if there is missing data (?) in the adult dataset
all.dat2 %>%
  select(species, X141:X335) %>% #selects for adult characters
  pivot_longer(cols = -species, names_to = 'var', values_to = 'val') %>% #long form: spp, character, state
  filter(val == '?') #filters for rows with ?

# categories for adult characters
body.col<-paste0('X', c(144:255)) #columns that code different characters 
rat.col<-paste0('X', c(159:161, 267, 272, 280, 283, 315, 322))
fem.col<-paste0('X', c(256:285))
mal.col<-paste0('X', c(285:335))
scale.col<-paste0('X', c(141:145, 147, 148, 151:153, 156, 164, 165, 167, 168, 173:183, 186:195, 197, 198, 201:203, 205, 206, 208, 210:213, 217:222, 226, 230:235, 237:245,251:253))
setae.col<-paste0('X', c(149:152, 155, 162, 169:172, 184:186, 196, 199, 200, 204, 207, 209, 214, 215, 217:219, 223:225, 254))

# find characters in both scale and setae list
ex1<-scale.col[scale.col %in% setae.col]

# look for invariant characters in scale dataset
# breakdown of steps
#a<-apply(all.dat2 %>% select(species, all_of(scale.col)), 2, function(a) length(unique(a)) == 1) #find character where unique=1, indicating invariance
#a1<-colnames(all.dat2 %>% select(species,all_of(scale.col))) #get colnames in the dataframe
#ex2<-a1[a]
ex2<-colnames(all.dat2 %>% select(species,all_of(scale.col)))[apply(all.dat2 %>% select(species, all_of(scale.col)), 2, function(a) length(unique(a)) == 1)]

scale.col<-scale.col[!(scale.col %in% ex1)]
scale.col<-scale.col[!(scale.col %in% ex2)]
scale.col ## SCALE DATASET
length(scale.col)

# look for invariant characters in setae dataset
ex3<-colnames(all.dat2 %>% select(species, all_of(setae.col)))[apply(all.dat2 %>% select(species, all_of(setae.col)), 2, function(a) length(unique(a)) == 1)]
setae.col<-setae.col[!(setae.col %in% ex1)]
setae.col<-setae.col[!(setae.col %in% ex3)]
setae.col ##SETAE DATASET
length(setae.col)

# Character column list for scales and setae traits. Some characters are repeated in both scale and setae, thus unique() was applied to obtain a non-redundant character list
scale.setae.col<-unique(c(scale.col,setae.col))

# Exclude invariant columns
ex4<-colnames(all.dat2 %>% select(species, all_of(scale.setae.col)))[apply(all.dat2 %>% select(species, all_of(scale.setae.col)), 2, function(a) length(unique(a)) == 1)]
scale.setae.col<-scale.setae.col[!(scale.setae.col %in% ex4)]
scale.setae.col ##SCALE+SETAE DATASET
length(scale.setae.col) #89 characters for adults

##### Larval characters #####
# Determine if there is missing data (?) in the larva dataset
all.dat2 %>%
  select(species, X5:X101) %>% #selects for adult characters
  pivot_longer(cols = -species, names_to = 'var', values_to = 'val') %>% #long form: spp, character, state
  filter(val == '?')%>% #filters for rows with ?
  filter(!grepl('flavifrons', species)) %>% filter(!grepl('nivalis', species)) %>% #remove flavifrons and nivalis
  group_by(var) %>% summarise(count = n())

# characters with ? and ordered characters
exclude.col<-paste0('X', c(23,25,26,39,40,52,54,56,94))
ordered.col<-paste0('X', c(13,14,17,20,62,67,72,84))

# select for larval characters, remove flavifrons and nivalis
larva.dat <- all.dat2 %>%
  select(species, X5:X101) %>% #selects for adult characters
  select(-any_of(exclude.col)) %>% #exclude characters with ?
  select(-any_of(ordered.col)) %>% #exclude ordered characters
  filter(!grepl('flavifrons', species)) %>% filter(!grepl('nivalis', species)) #remove flavifrons and nivalis
ncol(larva.dat)

larva.col<-colnames(larva.dat%>%select(X5:X101))
length(larva.col) # 80 characters for larvae

# check for invariant characters
larva.ex<-colnames(larva.dat %>%
                     select(species,all_of(larva.col)))[apply(larva.dat %>%
                     select(species,all_of(larva.col)), 2, function(a) length(unique(a)) == 1)]
# no invariant characters

##### Adult + larva dataset #####
adult.dat<-all.dat2 %>% 
  select(species, all_of(scale.setae.col)) %>%
  filter(!grepl('flavifrons', species)) %>% filter(!grepl('nivalis', species)) %>% #remove flavifrons and nivalis
  column_to_rownames('species')

larva.dat <- all.dat2 %>%
  select(species, X5:X101) %>% #selects for adult characters
  select(-any_of(exclude.col)) %>% #exclude characters with ?
  select(-any_of(ordered.col)) %>% #exclude ordered characters
  filter(!grepl('flavifrons', species)) %>% filter(!grepl('nivalis', species)) %>% #remove flavifrons and nivalis
  column_to_rownames('species')

###### Aedini tree #####
aedini.tree<-geiger::treedata(mcmc.tree, data = adult.dat, warnings = T, sort = T)$phy
aedini.tree<-drop.tip(aedini.tree, tip = out_group)

# lists taxa in clades A, B and Psorophora
clade_b<-extract.clade(aedini.tree, node = findMRCA(aedini.tree, tips = c('Aedes_Aedes_cinereus','Aedes_Fredwardsius_vittatus')))$tip.label
clade_a<-extract.clade(aedini.tree, node = findMRCA(aedini.tree, tips = c('Aedes_Abraedes_papago', 'Opifex_Nothoskusea_chathamicus')))$tip.label
psor<-extract.clade(aedini.tree, node = findMRCA(aedini.tree, tips = aedini.tree$tip.label[grepl('Psorophora', aedini.tree$tip.label)]))$tip.label

##### Habitat categories #####
# Add a column for clade and another column for habitat type
# Two categories
all.dat.hab2<-all.dat2 %>%
  mutate(clade = case_when(species %in% clade_a ~ 'Clade A',
                           species %in% clade_b ~ 'Clade B',
                           species %in% psor ~ 'Psorophora'),
         hab = factor(case_when(X336 == '01' ~ '0',
                                X336 == '02' ~ '0',
                                X336 == '23' ~ '3',
                                .default = X336)),
         hab = factor(case_when(hab == '0' ~ 'Open water', #freshwater ground pools
                                hab == '1' ~ 'Open water', #brackish-water ground and rock pools
                                hab == '2' ~ 'Container', #freshwater rock pools
                                hab == '3' ~ 'Container', #freshwater containers
                                hab == '4' ~ 'Container'), #crabholes
                      levels = c('Open water', 'Container')))

# Five categories
all.dat.hab5<-all.dat2 %>%
  mutate(clade = case_when(species %in% clade_a ~ 'Clade A',
                           species %in% clade_b ~ 'Clade B',
                           species %in% psor ~ 'Psorophora'),
         hab = factor(case_when(X336 == '01' ~ '0',
                                X336 == '02' ~ '0',
                                X336 == '23' ~ '3',
                                .default = X336)),
         hab = factor(case_when(hab == '0' ~ 'Freshwater ground pools', #freshwater ground pools
                                hab == '1' ~ 'Brackish ground and rock pools', #brackish-water ground and rock pools
                                hab == '2' ~ 'Freshwater rock pools', #freshwater rock pools
                                hab == '3' ~ 'Freshwater containers', #freshwater containers
                                hab == '4' ~ 'Crabholes'), #crabholes
                      levels = c('Freshwater ground pools', 'Brackish ground and rock pools','Freshwater rock pools','Freshwater containers','Crabholes')))
all.dat.hab5 %>% group_by(hab) %>% count()

# Raw categories (includes species with two categories eg. 01)
all.dat.hab.raw<-all.dat2 %>%
  mutate(clade = case_when(species %in% clade_a ~ 'Clade A',
                           species %in% clade_b ~ 'Clade B',
                           species %in% psor ~ 'Psorophora'),
         hab = factor(case_when(X336 == '01' ~ '01',
                                X336 == '02' ~ '02',
                                X336 == '23' ~ '23',
                                .default = X336)),
         hab = factor(case_when(hab == '0' ~ 'Freshwater ground pools', #freshwater ground pools
                                hab == '1' ~ 'Brackish ground and rock pools', #brackish-water ground and rock pools
                                hab == '2' ~ 'Freshwater rock pools', #freshwater rock pools
                                hab == '3' ~ 'Freshwater containers', #freshwater containers
                                hab == '4' ~ 'Crabholes', #crabholes
                                hab == '01' ~ '01',
                                hab == '02' ~ '02',
                                hab == '23' ~ '23'), 
                      levels = c('Freshwater ground pools', 'Brackish ground and rock pools','Freshwater rock pools','Freshwater containers','Crabholes','01','02','23')))
all.dat.hab.raw %>% group_by(hab) %>% count()

