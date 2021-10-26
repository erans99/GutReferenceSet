# Need 300G free memory space to run on 240k samples.
import numpy as np
import fastcluster
from scipy.cluster.hierarchy import fcluster
import pandas
import os
import time
from GutReferenceSet.Build_Species_Set import config


if __name__ == "__main__":
    method = 'average'
    th = [['s', 0.05], ['g', 0.15], ['f', 0.30]]

    for out_file in [config.clusters_file, config.clusters_links]:
        assert not os.path.exists(out_file), "%s exists. Exiting" % out_file
    for inp_file in [config.dists_mat16_triu_file, config.dists_mat16_index_file]:
        assert os.path.exists(inp_file), "%s doesn't exist. Exiting" % inp_file

    res = pandas.read_csv(config.dists_mat16_index_file, index_col=0)[['fasta']]
    num_files = len(res)
    print("Start clustering ", time.ctime())
    size_of_upper = int(num_files*(num_files-1)/2)
    print("Size %d*%d (condensed %d)" % (num_files, num_files, size_of_upper))
    upper_mem = np.memmap(config.dists_mat16_triu_file, dtype=config.MEM_DTYPE, mode='r', shape=(size_of_upper,))
    upper = np.array(upper_mem)
    del upper_mem
    print("Read into memory", time.ctime())

    link = fastcluster.linkage(upper, method=method, preserve_input=False)
    del upper
    pandas.DataFrame(data=link).to_csv(config.clusters_links)
    print("Ended creating linkage", time.ctime())

    for level, threshold in th:
            res['cluster_' + level] = fcluster(link, threshold, 'distance')
            print("End cutting level %s" % level, time.ctime())
    res.to_csv(config.clusters_file)
