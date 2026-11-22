#!/usr/bin/python3

import sys

def combine_files(args):
    table = {}
    for tsv_file in args:
        with open(tsv_file) as f:
            for line in f:
                row = line.rstrip().split('\t')
                if row[0] in table.keys():
                    table[row[0]] = table[row[0]]+row[1:]
                else:
                    table[row[0]] = row[1:]
    return table


if __name__ == "__main__":
    args = sys.argv[1:]
    table = combine_files(args)
    for key, val in table.items():
        print('\t'.join([key]+val))
