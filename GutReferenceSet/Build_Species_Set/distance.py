import pandas
import numpy
import os
import time
from GutReferenceSet.Build_Species_Set import config
from GutReferenceSet.Utils.mash import make_sketches, merge_sketches, create_dists, read_dists

NUM_SKETCHES_PER_JOB = 1000
NUM_MERGE_SKETCH_FILES = 5  # to merge together before creating dists
NUM_SKETCHES_PER_MERGED_FILE = NUM_SKETCHES_PER_JOB * NUM_MERGE_SKETCH_FILES
CONDENSE_BATCH_SIZE = 3 * 10 ** 8


def get_dists_file_name(f1, f2, out_dir):
    return os.path.join(out_dir, ("%s__%s.csv" % (os.path.splitext(os.path.basename(f1))[0],
                                                  os.path.splitext(os.path.basename(f2))[0])))


def unite_horizontal_blocks(sketch_files, i, n):
    if (i + 1) * NUM_SKETCHES_PER_MERGED_FILE > n:
        n_block_rows = n % NUM_SKETCHES_PER_MERGED_FILE
    else:
        n_block_rows = NUM_SKETCHES_PER_MERGED_FILE
    mem_offset = i * NUM_SKETCHES_PER_MERGED_FILE * n
    start_col = i * NUM_SKETCHES_PER_MERGED_FILE
    f1 = sketch_files[i]
    print("Started rows of %s" % f1, time.ctime())
    print("offset %d-%d" % (mem_offset * config.SIZE_IN_MEM, (mem_offset + n_block_rows * n) * config.SIZE_IN_MEM))
    mem_block = numpy.memmap(config.dists_mat16_file, mode="r+", dtype=config.MEM_DTYPE,
                             offset=mem_offset * config.SIZE_IN_MEM, shape=(n_block_rows, n))
    for f2 in sketch_files[i:]:
        block = read_dists(get_dists_file_name(f1, f2, config.cluster_dist_dir)).T.values
        if f1 == f2:
            block = numpy.triu(block)
        mem_block[:, start_col:start_col + block.shape[1]] = block
        start_col += block.shape[1]
    del mem_block
    print("Finished rows of %s" % f1, time.ctime())
    return


def create_upper_condense_matrix(n, mem_offset, mem_length, row_numbers):
    print("Started condensing %d to %d" % (mem_offset, mem_offset + mem_length), time.ctime())
    print("offset %d-%d" % (mem_offset * config.SIZE_IN_MEM, (mem_offset + mem_length) * config.SIZE_IN_MEM))

    dists_mat16 = numpy.memmap(config.dists_mat16_file, mode="r", dtype=config.MEM_DTYPE, shape=(n, n))
    mem_block = numpy.memmap(config.dists_mat16_triu_file, mode="r+", dtype=config.MEM_DTYPE,
                             offset=(mem_offset * config.SIZE_IN_MEM), shape=mem_length)
    writing_pos = 0
    for reading_row in row_numbers:
        mem_block[writing_pos:writing_pos + n - (reading_row + 1)] = dists_mat16[reading_row, (reading_row + 1):]
        writing_pos += n - (reading_row + 1)
    del dists_mat16
    del mem_block
    print("Finished condensing", time.ctime())
    return


def unite_all_dists(sketch_files, n):
    if os.path.exists(config.dists_mat16_file):
        assert not os.path.exists(config.dists_mat16_file), "%s exist" % config.dists_mat16_file
    if os.path.exists(config.dists_mat16_triu_file):
        assert not os.path.exists(config.dists_mat16_triu_file), "%s exist" % config.dists_mat16_triu_file

    # the 2 following lines creates the needed memory map, but then we do not need the pointer to it
    dists_mat16 = numpy.memmap(config.dists_mat16_file, mode="w+", dtype=config.MEM_DTYPE, shape=(n, n))
    del dists_mat16

    horizontal_args = []
    for i in range(len(sketch_files)):
        horizontal_args.append((sketch_files, i, n))

    # need to jump by 2, because consecutive files cause memory mis-writes (due to memory cache)
    # thus finish the evens, and then do the odds
    for start in [0, 1]:
        for args in horizontal_args[start::2]:
            unite_horizontal_blocks(*args)
        print("Sent %s half of unite jobs" % (start+1))

    print("Creating upper condensed matrix", time.ctime())

    # the 2 following lines creates the needed memory map, but then we do not need the pointer to it
    dists_mat16_triu = numpy.memmap(config.dists_mat16_triu_file, mode="w+", dtype=config.MEM_DTYPE,
                                    shape=int(n * (n - 1) / 2))
    del dists_mat16_triu

    mem_offset = 0
    mem_length = 0
    row_numbers = []
    condense_args = []
    for i in range(n):
        mem_length += n - i - 1
        row_numbers.append(i)
        if (mem_length > CONDENSE_BATCH_SIZE) or (i == (n-1)):
            condense_args.append((n, mem_offset, mem_length, row_numbers))
            mem_offset += mem_length
            mem_length = 0
            row_numbers = []

    # need to jump by 2, because consecutive files cause memory mis-writes (due to memory cache)
    # thus finish the evens, and then do the odds
    for start in [0, 1]:
        for args in condense_args[start::2]:
            create_upper_condense_matrix(*args)
        print("Sent %s half of upper condense jobs" % (start + 1))


if __name__ == '__main__':
    fasta_files = pandas.read_csv(config.metadata_filtered_file, index_col=0).index

    raw_sketch_dir = os.path.join(config.cluster_sketch_dir, "raw")
    os.makedirs(raw_sketch_dir, exist_ok=True)
    os.makedirs(config.cluster_dist_dir, exist_ok=True)

    raw_sketch_files = []
    for i in range(0, len(fasta_files), NUM_SKETCHES_PER_JOB):
        raw_sketch_files.append(os.path.join(raw_sketch_dir, "raw_part_%d_%d.msh" %
                                             (i, min(i + NUM_SKETCHES_PER_JOB, len(fasta_files)))))
        make_sketches(raw_sketch_files[-1], fasta_files[i:i + NUM_SKETCHES_PER_JOB], True)

    print("Merging sketches", time.ctime())
    merged_sketch_files = []
    for i in range(0, len(raw_sketch_files), NUM_MERGE_SKETCH_FILES):
        merged_sketch_files.append(os.path.join(config.cluster_sketch_dir, "merged_part_%d_%d.msh" %
                                                (i * NUM_SKETCHES_PER_JOB,
                                                 min((i + NUM_MERGE_SKETCH_FILES) * NUM_SKETCHES_PER_JOB,
                                                     len(fasta_files)))))
        merge_sketches(merged_sketch_files[-1], raw_sketch_files[i:i + NUM_MERGE_SKETCH_FILES], True, True)

    print("Creating dist matrices", time.ctime())
    file_list = []
    for i, fn1 in enumerate(merged_sketch_files):
        file_list.append(pandas.read_csv(fn1.replace(".msh", ".csv"), header=None, index_col=None))
        file_list[-1]['sketch_file'] = fn1
        for fn2 in merged_sketch_files[i:]:
            dists_file_name = get_dists_file_name(fn1, fn2, config.cluster_dist_dir)
            create_dists(fn1, fn2, dists_file_name, 2, True,)

    file_list = pandas.concat(file_list, ignore_index=True)
    file_list.columns = ['fasta', 'sketch_file']
    assert (file_list['fasta'].values == fasta_files).all(), "Input and output file list don't agree"
    file_list.to_csv(config.dists_mat16_index_file)

    print("Uniting", time.ctime())
    unite_all_dists(merged_sketch_files, len(fasta_files))

    print("Done", time.ctime())

