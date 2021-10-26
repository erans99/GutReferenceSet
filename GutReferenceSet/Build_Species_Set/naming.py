# this runs on 50 threads of the computer

import os
import pandas
from GutReferenceSet.Build_Species_Set import config
from GutReferenceSet.Utils.config import gtdb_exe, gtdb_info

NUM_TH = 50
SUBDIR = "out_GTDB"

if __name__ == '__main__':
    if not os.path.exists(os.path.join(config.naming_dir, SUBDIR, "gtdbtk.bac120.summary.tsv")):
        os.makedirs(config.naming_dir, exist_ok=True)
        reps = pandas.read_csv(config.reps_file, index_col=0)
        batch_file = os.path.join(config.naming_dir, "%s_fastas.csv" % SUBDIR)
        reps[['cluster_s']].to_csv(batch_file, sep='\t', header=False)
        os.environ["GTDBTK_DATA_PATH"] = gtdb_info
        com = "%s classify_wf --batchfile %s --out_dir %s --cpus %d" % (gtdb_exe, batch_file,
                                                                        os.path.join(config.naming_dir, SUBDIR), NUM_TH)
        os.system(com)

    if os.path.exists(os.path.join(config.naming_dir, SUBDIR, "gtdbtk.bac120.summary.tsv")):
        bac_inp = os.path.join(config.naming_dir, SUBDIR, "gtdbtk.bac120.summary.tsv")
        df_bac = pandas.read_csv(bac_inp, index_col=0, sep='\t')
        print("Bacteria:", df_bac.shape)
    else:
        df_bac = pandas.DataFrame()

    if os.path.exists(os.path.join(config.naming_dir, SUBDIR, "gtdbtk.ar122.summary.tsv")):
        arc_inp = os.path.join(config.naming_dir, SUBDIR, "gtdbtk.ar122.summary.tsv")
        df_arc = pandas.read_csv(arc_inp, index_col=0, sep='\t')
        print("Archaea:", df_arc.shape)
    else:
        df_arc = pandas.DataFrame()

    all_reps = pandas.read_csv(config.reps_file, index_col=0)
    if len(df_bac) + len(df_arc) != len(all_reps):
        print("Len named (%d) and len of representatives (%d) don't match" % (len(df_bac) + len(df_arc), len(all_reps)))
    df = pandas.concat([df_bac, df_arc])
    df.sort_index(inplace=True)

    for i, tax in enumerate('dpcofgs'):
        df[tax] = df.classification.str.split(';').str[i]
    GTDB_warn = df[df.warnings ==
                   'Genome not assigned to closest species as it falls outside its pre-defined ANI radius'].index
    df['GTDB_warn'] = False
    df.loc[GTDB_warn, 'GTDB_warn'] = True
    df['GTDB_no_s'] = (df.s == 's__')
    df.index = "Rep_" + df.index.astype(str)

    df.to_csv(config.reps_GTDB_names)
