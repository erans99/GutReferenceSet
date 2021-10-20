import os
import pandas as pd
from GutReferenceSet.Utils.config import checkm, prodigal, hmmer, pplacer

results_file = os.path.join('storage', 'bin_stats_ext.tsv')


def run(input_dir, output_dir, threads=2):
    if os.path.exists(output_dir):
        raise Exception(f'quality output dir should not exist {output_dir}')  # this will crash checkM
    os.mkdir(output_dir)

    os.environ['PATH'] = os.environ['PATH'] + os.pathsep + os.pathsep.join([prodigal, hmmer, pplacer])
    cmd = f'{checkm} lineage_wf -t {threads} -x fa {input_dir} {output_dir}'
    os.system(cmd)


def read_results(output_dir):
    if not os.path.exists(os.path.join(output_dir, results_file)):
        raise Exception(f'quality was not generated for {os.path.join(output_dir)}')

    results = pd.read_csv(os.path.join(output_dir, results_file), sep='\t', index_col=0, header=None).iloc[:, 0]
    results = results.apply(lambda x: eval(x))
    results = pd.DataFrame(index=results.index, data=list(results))
    results['quality_dir'] = output_dir
    return results
