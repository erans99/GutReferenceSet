import os
import time
import pandas

from GutReferenceSet.Utils.config import mash_exe


def make_sketches(out_file, fs, save_file_list=True):
    '''
    :param out_file: name of the mash output file ('.msh' will be added, unless already exists)
    :param fs: files list to mash
    :param save_file_list: bool, whether to output also the file list (as csv, in same location as the output msh file)
    :return: None
    '''
    if os.path.splitext(out_file)[1] == ".msh":
        out_file = os.path.splitext(out_file)[0]
    out_csv = out_file + ".csv"
    assert not os.path.exists(out_file + ".msh"), "Sketch %s file exists. Exiting" % (out_file + ".msh")
    assert not os.path.exists(out_csv), "File list %s file exists. Exiting" % out_csv

    pandas.Series(fs).to_csv(out_csv, index=False, header=False)
    com = "%s sketch -s 1e4 -l %s -o %s" % (mash_exe, out_csv, out_file)
    res = os.system(com)
    if res != 0:
        raise Exception("Mash on %s did not succeed" % out_csv)
    if not save_file_list:
        com = "rm %s" % out_csv
        os.system(com)


def merge_sketches(out_file, fs, save_file_list=True, remove_inputs=True):
    '''
    :param out_file: name of the mash output file ('.msh' will be added, unless already exists)
    :param fs: files list to mash
    :param save_file_list: bool, whether to output also the file list (as csv, in same location as the output msh file)
    :param remove_inputs: bool, whether or remove the input .msh and .csv files
    :return: None
    '''
    if os.path.splitext(out_file)[1] == ".msh":
        out_file = os.path.splitext(out_file)[0]
    out_csv = out_file + ".csv"
    assert not os.path.exists(out_file + ".msh"), "Sketch %s file exists. Exiting" % (out_file + ".msh")
    assert not os.path.exists(out_csv), "File list %s file exists. Exiting" % out_csv

    com = "%s paste %s" % (mash_exe, out_file)
    file_list = []
    for f in fs:
        com += " %s" % f
        if save_file_list:
            file_list.append(pandas.read_csv(f.replace(".msh", ".csv"), header=None))
    res = os.system(com)
    if res != 0:
        raise Exception("Mash on %s did not succeed" % out_csv)

    if save_file_list:
        assert not os.path.exists(out_csv), "File list %s file exists. Exiting" % out_csv
        pandas.concat(file_list, ignore_index=True).to_csv(out_csv, index=False, header=False)

    if remove_inputs:
        com = "rm -f %s"
        for f in fs:
            os.system(com % f)
            os.system(com % f.replace(".msh", ".csv"))


def create_dists(f1, f2, out_name, num_threads=2, pr=False):
    assert not os.path.exists(out_name), "Distances file %s file exists. Exiting" % out_name
    if pr:
        print("Create dists %s %s" % (os.path.basename(f1), os.path.basename(f2)), time.ctime())
    com = "%s dist -t -p %d %s %s > %s" % (mash_exe, num_threads, f1, f2, out_name)
    os.system(com)
    if pr:
        print("Done", time.ctime())


def read_dists(out_name):
    return pandas.read_csv(out_name, delimiter='\t', index_col=0, header=0)
