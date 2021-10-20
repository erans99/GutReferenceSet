# note: this code assumes all the single fasta files already exist in the singles_dir directory

import os
import pandas
import numpy
from GutReferenceSet.Utils.config import phy_path, phy_names, python_exe_dir, mash_exe, Segata_desc, \
    phylophlan_path
from GutReferenceSet.Build_Species_Set.config import singles_dir

NUM_TH = 50
CREATE_TREE = True
PHYPHLAM_NAME = True


if __name__ == '__main__':
    os.makedirs(phy_path, exist_ok=True)
    os.chdir(phy_path)

    if PHYPHLAM_NAME and (not os.path.exists(phy_names)):
        if not os.path.exists(os.path.join(phy_path, "output.tsv")):
            os.environ["PATH"] += os.pathsep + os.path.dirname(mash_exe)

            com = "ln -s  %s/phylophlan_ ." % phylophlan_path
            os.system(com)

            com = "ln -s  %s/phylophlan_databases ." % phylophlan_path
            os.system(com)

            assert os.path.exists(singles_dir), \
                "Single fasta directory doesn't exist"
            com = "%s/phylophlan_metagenomic -i %s -o output --nproc %d -n 1 -d SGB.Jan19 --verbose" % \
                  (python_exe_dir, singles_dir, NUM_TH)
            os.system(com)

        if os.path.exists(Segata_desc):
            Segata4930 = pandas.read_csv(Segata_desc, index_col=0).index

        df = pandas.read_csv(os.path.join(phy_path, "output.tsv"), index_col=0, delim_whitespace=True)
        col = df.columns[0]
        df = df.iloc[:, 0].str.split(':', expand=True)
        df.columns = col.split(':')
        df.avg_dist = df.avg_dist.astype(float)
        df['num_SGB'] = df['[u|k]_[S|G|F]GBid'].str.split("_").str[-1].astype(int)
        if os.path.exists(Segata_desc):
            df['is_human'] = numpy.isin(df['num_SGB'], Segata4930)
        df['is_SGB_dist'] = (df['avg_dist'] < 0.05)

        df['num_rep'] = df.index.str.split("_").str[-1].astype(int)
        df.sort_values('num_rep', inplace=True)
        df.to_csv(phy_names)

    if CREATE_TREE:
        assert not os.path.exists(os.path.join(phy_path, "output_tree")), "Tree directory exists!"
        print("Creating tree. This takes a day and a half, dont be discouraged!")
        com = ("%s/phylophlan -i %s -o output_tree -d phylophlan -f supermatrix_aa.cfg --nproc %d " + \
               "--diversity high --fast --genome_extension .fa --verbose") % \
              (python_exe_dir, singles_dir, NUM_TH)
        os.system(com)
