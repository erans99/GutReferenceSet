import pandas
from GutReferenceSet.Build_Species_Set import config

NUM_IN_PART = 100
compute_centrality = False
out_reps = True


if __name__ == '__main__':
    clusts = pandas.read_csv(config.clusters_file, index_col=0)

    metadata = pandas.read_csv(config.metadata_filtered_file, index_col=0)
    clusts.set_index('fasta', inplace=True)
    metadata = metadata.join(clusts[['cluster_s', 'cluster_g', 'cluster_f']])
    metadata.to_csv(config.metadata_filtered_with_clusts)

    res = []
    for species, species_df in metadata.groupby('cluster_s'):
        res.append([int(species), len(species_df), species_df[config.quality_score_col].idxmax()])
    res = pandas.DataFrame(res, columns=['cluster_s', 'num_ASMs', 'representative'])
    res = res.set_index('representative')
    res.to_csv(config.reps_file)

    tax = clusts.loc[res.index].set_index('cluster_s')
    tax['family_name'] = 'fBin__' + tax['cluster_f'].astype(str)
    tax['genus_name'] = tax['family_name'] + '|gBin__' + tax['cluster_g'].astype(str)
    tax['species_name'] = tax['genus_name'] + '|sBin__' + tax.index.astype(str)
    tax.to_csv(config.species_taxonomy)
