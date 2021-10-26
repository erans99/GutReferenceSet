import os
import numpy

build_dir = "XXXX" # work directory, in it should already exist full_metadata.csv

Min_Completeness = 70
Max_Contamination = 5

Coef_Completeness = 1
Coef_Contamination = -5
Coef_log_N50 = 15

quality_score_col = f'Score_{Coef_Completeness}_{Coef_Contamination}_{Coef_log_N50}'

species_level_mash_distance = 0.05

quality_dir = os.path.join(build_dir, "quality")
dists_dir = os.path.join(build_dir, "dists")
cluster_dir = os.path.join(dists_dir, "between_person")

metadata_file = os.path.join(build_dir, "full_metadata.csv")
metadata_with_quality_file = os.path.join(build_dir, "full_metadata_w_quality.csv")
metadata_with_filter_file = os.path.join(build_dir, "full_metadata_w_filter.csv")
metadata_filtered_file = os.path.join(build_dir, "full_metadata_filtered.csv")
metadata_filtered_with_clusts = os.path.join(build_dir, "full_metadata_filtered_w_clusts.csv")

quality_stats = os.path.join(build_dir, "quality_stats.png")

cluster_sketch_dir = os.path.join(cluster_dir, "sketches")
cluster_dist_dir = os.path.join(cluster_dir, "dists")

MEM_DTYPE = numpy.float16
SIZE_IN_MEM = 2  # number of bytes the config.MEM_DTYPE takes in memory

dists_mat16_index_file = os.path.join(cluster_dir, "dists_mat16_index.csv")
dists_mat16_file = os.path.join(cluster_dir, "dists_mat16.dat")
dists_mat16_triu_file = os.path.join(cluster_dir, "dists_mat16_triu.dat")

clusters_file = os.path.join(build_dir, "clusters_avg_dist.csv")
clusters_links = os.path.join(build_dir, "clusters_avg_dist_link.csv")
reps_file = os.path.join(build_dir, "chosen_representatives.csv")
species_taxonomy = os.path.join(build_dir, "taxonomy.csv")

naming_dir = os.path.join(build_dir, "GTDB")
reps_GTDB_names = os.path.join(build_dir, "GTDB_naming.csv")

singles_dir = os.path.join(build_dir, "singles")
phy_path = os.path.join(build_dir, "PhyPhlan")
phy_names = os.path.join(phy_path, "df_out_PhyPhlan.csv")
