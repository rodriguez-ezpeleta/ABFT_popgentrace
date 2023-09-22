sed -i '1d' populations.treemix
gzip populations.treemix
for f in `seq 0 10`; do
treemix -i populations.treemix.gz -m $f -o stem_m${f}
done


##### LOOK AT TREEMIX RESULTS:
#####   A)EXPLORE LIKELIHOODS FROM *.lik FILES, WHICH ARE LIKELIHOOD ASSOCIATED TO EACH TREE + MIGRATION EVENTS ######
#####            for f in `seq 1 10`; do cat stem_m${f}.llik; done
#####   B)PLOT MIGRATION EVENTS IN R
#####

