#!usr/bin/env python
# -*- coding: utf-8 -*-


import os


def get_files(dir_name):
    # create a list of file and sub directories
    # names in the given directory
    files = os.listdir(dir_name)
    all_files = []
    # Iterate over all the entries
    for entry in files:
        # Create full path
        full_path = os.path.join(dir_name, entry)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(full_path):
            all_files += get_files(full_path)
        else:
            all_files.append(full_path)
    return all_files


def get_names(path):
    return [(file.split('/')[3][:-3]) for file in get_files(path)]


def get_dir_names(filenames):
    dirn = set()
    subd = set()
    subsub = set()
    for filename in filenames:
        dir_splt = filename.split('/')
        dirn.add(dir_splt[0])
        subd.add(dir_splt[1])
        subsub.add(dir_splt[2])
    return next(iter(dirn)), next(iter(subd)), next(iter(subsub))


def get_dir_name(filename):
    subd, subsub = set(), set()
    names = filename.split('/')
    for i in range(len(names)):
        subd.add(names[1])
        subsub.add(names[2])
    return next(iter(subd)), next(iter(subsub))


def make_me_a_folder(folder_name):
    os.getcwd()
    try:
        os.makedirs(folder_name)
    except OSError:
        print("Directory already exists")
    return os.path.exists(folder_name)


def get_full_name(dir_out, sud_dir, ssub):
    return os.path.join(dir_out, sud_dir, ssub)


def get_fasta_files(dir_name):
    infiles = []
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            input_files = os.path.join(path, name)
            if input_files.endswith('fa.gz') or input_files.endswith('.fa') \
                    or input_files.endswith('.fasta') or input_files.endswith('.fa.gz') \
                    or input_files.endswith('.fna') or input_files.endswith('.fna.gz'):
                infiles.append(input_files)
    return infiles

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z