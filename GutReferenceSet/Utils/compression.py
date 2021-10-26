import os


def unzip_files(zipped_paths):
    print("Unzipping %d files" % len(zipped_paths))
    for f in zipped_paths:
        com = 'gunzip -r %s' % f
        os.system(com)
    print("Done")


def zip_files(unzipped_paths):
    print("Zipping %d files" % len(unzipped_paths))
    for f in unzipped_paths:
        com = 'gzip -r %s' % f
        os.system(com)
    print("Done")
