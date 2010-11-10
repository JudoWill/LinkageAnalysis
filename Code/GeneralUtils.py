import csv
import os



def make_mapping_dict(in_file):
    
    mapping_dict = {}
    with open(in_file) as handle:
        for row in csv.DictReader(handle, delimiter = '\t'):
            if row['name'] == 'None':
                mapping_dict[row['key']] = None
            else:
                mapping_dict[row['key']] = row['name']
    return mapping_dict


def mapping_func(mapping_dict, name):
    return mapping_dict.get(name.lower(), None)

def gi_from_path(path):

    fname = path.split(os.sep)[-1]
    gi = fname.split('.')[0]
    return gi
