salloc -N1 -n1 -c4 --mem=5gb -t 4:00:00
export PATH=/work/soghigian_lab/apps/iqtree-latest/bin:$PATH

# collapse nodes below 90 support value
iqtree2 -t 12862_2017_1092_MOESM14_ESM.nwk -minsupnew 90
mv 12862_2017_1092_MOESM14_ESM.nwk.collapsed 12862_2017_1092_MOESM14_ESM.nwk.collapsed90.treefile

# rename tip labels
cp 12862_2017_1092_MOESM14_ESM.nwk.collapsed90.treefile S2017.collapsed90.treefile
while read -r pattern replacement; do
sed -i "s/$pattern/$replacement/" S2017.collapsed90.treefile
done < tipLabels_S2017.txt
