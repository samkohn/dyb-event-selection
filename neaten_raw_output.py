import os
import multiprocessing

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

with multiprocessing.Pool() as pool:
    file_list = os.listdir(base_dir)
    split_by_1000s = [file_list[i : i + 1000] for i in range(0, len(file_list), 1000)]
    pool.map(rename, split_by_1000s)
