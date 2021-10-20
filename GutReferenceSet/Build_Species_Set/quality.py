import os
import pandas as pd
from GutReferenceSet.Utils import checkm, compression
from GutReferenceSet.Utils.data_frame import simple_statistics
from GutReferenceSet.Build_Species_Set.config import metadata_file, metadata_with_quality_file, quality_dir

genome_params = ['Genome size', 'Completeness', 'Contamination']  # ['# genomes']
marker_params = ['marker lineage', '# markers', '# marker sets']
base_pairs_params = ['GC', 'GC std', '# ambiguous bases']  # ['GCN0', 'GCN1', 'GCN2', 'GCN3', 'GCN4', 'GCN5+']
scafold_params = ['N50 (scaffolds)', '# scaffolds', 'Longest scaffold', 'Mean scaffold length']
contig_params = ['N50 (contigs)', '# contigs', 'Longest contig', 'Mean contig length']
genes_params = ['Coding density', '# predicted genes']  # ['Translation table']
unknown = ['0', '1', '2', '3', '4', '5+']

params2keep = genome_params + base_pairs_params + contig_params + genes_params + ['quality_dir']


def get_params(b):

    params = []
    for b_start, b_end in zip(b[:-1], b[1:]):
        in_dir = os.path.join(quality_dir, f'files_{b_start}_{b_end}')
        out_dir = os.path.join(in_dir, 'checkM')
        params.append([b_start, b_end, in_dir, out_dir])

    return params


if __name__ == '__main__':
    run_checkm = True
    get_results = True
    zip_files = False

    trds = 2
    max_r = 150
    mem_def = '35G'
    batch_size = 100

    metadata = pd.read_csv(metadata_file, index_col=0)
    run = metadata.copy().sort_values(by='Method', ascending=False)
    run = run.loc[run['DoNotTake'].isna(), 'AssemblyName']
    batches = range(0, run.shape[0] + batch_size, batch_size)

    if run_checkm:
        for batch_start, batch_end, input_dir,  output_dir in get_params(batches):

            # create dirs with links
            if os.path.exists(input_dir):
                raise Exception(f'dir already exist {input_dir}')
            os.makedirs(input_dir)
            for index, value in run.iloc[batch_start:batch_end].items():
                os.system(f'ln -s {index} {os.path.join(input_dir, os.path.splitext(value)[0])}.fa')

            # run checkM
            checkm.run(input_dir, output_dir, trds)

    if get_results:

        all_results = []
        for _, _, _, output_dir in get_params(batches):

            # decompress results if needed
            results_file = os.path.join(output_dir, checkm.results_file)
            if (not os.path.exists(results_file)) and os.path.exists(results_file + '.gz'):
                compression.unzip_files(zipped_paths=[results_file])

            # read results
            all_results.append(checkm.read_results(output_dir=output_dir))

        # summarize results
        all_results = pd.concat(all_results)
        metadata = metadata[list(set(metadata.columns) - set(params2keep))]
        metadata['temp'] = metadata['AssemblyName'].str.replace('.fa', '').str.replace('.fna', '')
        metadata = metadata.join(all_results[params2keep], how='left', on='temp')
        metadata = metadata.drop('temp', axis=1)
        metadata.to_csv(metadata_with_quality_file)

    if zip_files:

        # compress results
        files2zip = [output_dir for _, _, _, output_dir in get_params(batches)]
        compression.zip_files(unzipped_paths=files2zip)

    print(f'ran on {run.shape[0]} assemblies')

    simple_statistics(metadata[params2keep])
