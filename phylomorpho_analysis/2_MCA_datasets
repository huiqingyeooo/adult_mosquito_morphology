# Generate datasets with dimension coordinates, species, habitats and clades

all.dat.hab.raw.dim <- data.frame(adult.mca$ind$coord) %>%
  rownames_to_column('species') %>%
  left_join(all.dat.hab.raw %>% select(species, hab, clade))

all.dat.hab5.dim <- data.frame(adult.mca$ind$coord) %>%
  rownames_to_column('species') %>%
  left_join(all.dat.hab5 %>% select(species, hab, clade))

all.dat.hab2.dim <- data.frame(adult.mca$ind$coord) %>%
  rownames_to_column('species') %>%
  left_join(all.dat.hab5 %>% select(species, hab, clade))

# Read in distribution data
dat<-read.csv('aedini_hostPref_biogeo_updated.csv', header = T)

bg.dat <- dat %>% rename(species = RHKspp_in_phylotree) %>%
  left_join(all.dat.hab2, by="species") %>%
  select(-contains("X"))

# Dataset for phylolm (all individuals)
phylo.dat <-data.frame(Dim.1 = adult.mca$ind$coord[, 1],
                       Dim.2 = adult.mca$ind$coord[, 2],
                       Dim.3 = adult.mca$ind$coord[, 3])%>%
  rownames_to_column('species') %>%
  left_join(bg.dat, by="species") %>%
  select(-contains("X"), -"Distribution", -"unep2024",-"motw_2021")%>%
  mutate(across(c("clade","hab","pal","afr","ind","neo","nea","aus","oldWorld","newWorld"), as.factor))%>%
  column_to_rownames('species')
str(phylo.dat)
