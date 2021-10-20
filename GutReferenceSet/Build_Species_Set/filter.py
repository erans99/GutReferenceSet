import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from GutReferenceSet.Utils.data_frame import simple_statistics
from GutReferenceSet.Utils.mash import create_dists, make_sketches, read_dists
from GutReferenceSet.Build_Species_Set.config import \
    dists_dir, quality_stats, species_level_mash_distance, \
    metadata_with_quality_file, metadata_with_filter_file, metadata_filtered_file, \
    Min_Completeness, Max_Contamination, Coef_Completeness, Coef_Contamination, Coef_log_N50, quality_score_col


thresholds = {
    'Completeness': ('>=', Min_Completeness),
    'Contamination': ('<=', Max_Contamination)
}

person_defining_cols = ['Source', 'RegistrationCode']


def get_check_quality_condition(df):
    return df['DoNotTake'].isna()


def filter_by_quality(df):

    for k, (sign, value) in thresholds.items():
        check_quality = get_check_quality_condition(df)

        if df.loc[check_quality, k].isna().sum() == 0:
            if sign == '>=':
                good_quality = df[k] >= value
            elif sign == '<=':
                good_quality = df[k] <= value
            else:
                raise Exception(f'Unknown sign "{sign}"')

            df.loc[check_quality & (~good_quality), 'DoNotTake'] = 'Quality'
            print(f'checked {sum(check_quality)} and lost {sum(check_quality & (~good_quality))} on {k}')

        else:
            raise Exception(f'Some genomes do not have {k}')

    return df


def compute_quality_stats(df):
    df['l_N50'] = np.log10(df['N50 (contigs)'])
    df['l_length'] = np.log10(df['Genome size'])

    cols = ['Completeness', 'Contamination', 'l_N50', 'l_length', '# contigs']

    st = {}
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    for i, c in enumerate(cols):
        st[c] = [df[c].mean(), df[c].std(), df[c].quantile(0.05), df[c].quantile(0.95)]
        ax = axes[int(i / 3)][i % 3]
        ax.hist(df[c], bins=100)
        h = ax.get_ylim()[1]
        ax.vlines(st[c][0], 0, h, color='r', ls='-')
        ax.vlines(st[c][0] + st[c][1], 0, h, color='r', ls='--')
        ax.vlines(st[c][0] - st[c][1], 0, h, color='r', ls='--')
        ax.set_xlabel("%s (mean %g +- 1 std %g)" % (c, st[c][0], st[c][1]))
    plt.savefig(quality_stats)
    plt.close()
    st = pd.DataFrame.from_dict(st, orient='index', columns=['mean', 'std', 'q05', 'q95'])
    st['diff_qaunt'] = st['q95'] - st['q05']
    st['Completeness_std_ratio'] = st.loc['Completeness', 'std'] / st['std']
    st['Completeness_diff_quant_ratio'] = st.loc['Completeness', 'diff_qaunt'] / st['diff_qaunt']
    st.to_csv(quality_stats.replace('.png', '.csv'))


def calculate_quality_score(df):
    df[quality_score_col] = \
        Coef_Completeness*df['Completeness'] + \
        Coef_Contamination*df['Contamination'] + \
        Coef_log_N50*np.log10(df['N50 (contigs)'])
    # the contamination minus comes from the config

    return df


def mash_from_best(person_df):
    person_df = person_df.sort_values(quality_score_col, ascending=False)

    source, registration = person_df.loc[person_df.index[0], person_defining_cols]
    person = f'{source}_{registration}'.replace('|', '_')
    sketch_name = os.path.join(dists_dir, 'within_person', f'{person}.msh')
    mash_name = os.path.join(dists_dir, 'within_person', f'{person}_{person}.csv')

    make_sketches(out_file=sketch_name, fs=person_df.index.tolist(), save_file_list=True)
    create_dists(f1=sketch_name, f2=sketch_name, out_name=mash_name, num_threads=2, pr=False)
    dists = read_dists(out_name=mash_name)

    included_genomes = [person_df.index[0]]
    for considered_genome in person_df.index[1:]:
        if dists.loc[considered_genome, included_genomes].min() > species_level_mash_distance:
            included_genomes.append(considered_genome)

    person_df.loc[person_df.index.difference(included_genomes), 'DoNotTake'] = 'Duplicated species within person'

    return person_df


def filter_by_person(df):
    multiple_genomes = df[person_defining_cols].duplicated(keep=False)

    df.loc[multiple_genomes] = df.loc[multiple_genomes].groupby(person_defining_cols).apply(mash_from_best)\
        .droplevel(person_defining_cols)

    return df


if __name__ == '__main__':
    batch_size = 1000

    metadata = pd.read_csv(metadata_with_quality_file, index_col=0)
    metadata = filter_by_quality(metadata)
    check_quality_condition = get_check_quality_condition(metadata)
    compute_quality_stats(metadata.loc[check_quality_condition].copy())
    metadata = calculate_quality_score(metadata)

    os.makedirs(os.path.join(dists_dir, 'within_person'), exist_ok=True)

    res = []
    run = metadata.loc[check_quality_condition].sort_values(person_defining_cols).copy()
    end = 0
    while end < run.shape[0]:
        start = end
        end = start + batch_size
        if end >= run.shape[0]:
            end = run.shape[0]
        else:
            end_person = run.iloc[end-1][person_defining_cols].values
            while (end < run.shape[0]) and (run.iloc[end][person_defining_cols].values == end_person).all():
                end = end + 1
        res.append(filter_by_person(run.iloc[start:end]))

    res = pd.concat(res)
    if res.shape[0] != sum(check_quality_condition):
        raise Exception('Not all jobs returned properly')
    metadata = res.loc[metadata.index]
    metadata.to_csv(metadata_with_filter_file)

    print(f'ran on {metadata.shape[0]} assemblies')
    simple_statistics(metadata[[quality_score_col, 'DoNotTake']])

    check_quality_condition = get_check_quality_condition(metadata)
    metadata = metadata.loc[check_quality_condition]
    metadata.to_csv(metadata_filtered_file)
