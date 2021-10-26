import argparse
import itertools
import os
import multiprocessing
import random

min_run = 21221
max_run = 72455

base_dir = '/global/cscratch1/sd/skohn/dyb30/raw_adtime/EH3'

def rename(filenames):
    for filename in filenames:
        if len(filename) == 5:  # subdirectory
            continue
        if filename[0] == 'e':  # 'events_ad?_?????_????.root'
            run_num_str = filename[11:16]
        elif filename[0] == 'm':  # 'muons_?????_????.root'
            run_num_str = filename[6:11]
        subdir = run_num_str[0:3] + '00'
        new_location = os.path.join(base_dir, subdir, filename)
        if os.path.isfile(new_location):
            continue
        os.makedirs(os.path.join(base_dir, subdir), exist_ok=True)
        os.link(  # hard link
            os.path.join(base_dir, filename),
            new_location
        )
    return


def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks

    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"

    Taken from python itertools recipes.
    """
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def remove(file_nlink_pairs):
    for filename_nlinks in file_nlink_pairs:
        if filename_nlinks is None:
            continue
        filename, n_hardlinks = filename_nlinks
        if not os.path.isfile(filename):
            continue
        if n_hardlinks == 2:
            if random.random() > 0.9999:
                print(f'Removing {os.path.abspath(filename)}')
            os.remove(os.path.abspath(filename))
        elif n_hardlinks > 2:
            print(f'anomalous n_hardlinks = {n_hardlinks} for {os.path.abspath(filename)}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--rename', action='store_true')
    parser.add_argument('--remove-old', action='store_true')
    args = parser.parse_args()
    if args.rename:
        with multiprocessing.Pool() as pool:
            file_list = os.listdir(base_dir)
            split_by_1000s = [file_list[i : i + 1000] for i in range(0, len(file_list), 1000)]
            pool.map(rename, split_by_1000s)
    if args.remove_old:
        dir_entry_list = os.scandir(base_dir)
        file_nlink_list = [
            (entry.path, entry.stat().st_nlink) for entry in dir_entry_list
            if entry.name[0] in 'em'
        ]
        split_by_1000s = grouper(file_nlink_list, 1000)
        with multiprocessing.Pool() as pool:
            pool.map(remove, split_by_1000s)
