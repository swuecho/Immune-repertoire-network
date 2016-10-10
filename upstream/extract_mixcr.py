#!/usr/bin/python

import os, sys, re
import string, csv
import glob
import time
import argparse


def trim_name_suffix(name):
    return name.split('*')[0]


def find_barcode(mixcr_file):
    return mixcr_file.rstrip(".txt")


def mixcr_to_csv(mixcr_file, trb_or_igh):
    mixcr_file_base = os.path.basename(mixcr_file)
    mixcr_file_dir = os.path.dirname(mixcr_file)
    barcode = find_barcode(mixcr_file_base)
    output_filename = os.path.join(
        os.path.dirname(mixcr_file_dir),
        '_'.join(['mixcr', trb_or_igh, 'extracted']), mixcr_file_base.replace(
            ".txt", ".csv"))
    output_fh = open(output_filename, 'wb')
    writer = csv.writer(output_fh)
    header_info = ['vgene', 'jgene', 'dgene', 'cdr3', 'count', 'barcode']
    writer.writerow(header_info)

    input_fh = open(mixcr_file, 'r')
    mixcr_reader = csv.reader(input_fh, delimiter='\t')
    # This skips the first row of the CSV file.
    # csvreader.next() also works in Python 2.
    #
    # in rare case, file size is zero
    if os.path.getsize(mixcr_file):
        next(mixcr_reader)
        print("!!! valid cdr3 in sample %s" % barcode)
    else:
        print("!!! NO valid cdr3 in sample %s" % barcode)

    for row in mixcr_reader:
        vgene, jgene, dgene, count, cdr3 = row
        if vgene and jgene and len(cdr3) > 7 and not (
                '*' in cdr3 or '~' in cdr3 or '_' in cdr3):
            result = [trim_name_suffix(vgene), trim_name_suffix(jgene),
                      trim_name_suffix(dgene), cdr3, count, barcode]
            writer.writerow(result)
        else:
            print(row)

    output_fh.close()
    input_fh.close()


def extract_and_move_file_in_dir(work_dir):
    start = time.time()
    names = glob.glob(work_dir + "/mixcr_trb/*.txt")
    for mixcr_result_x in names:
        mixcr_to_csv(mixcr_result_x, 'trb')

    names = glob.glob(work_dir + "/mixcr_igh/*.txt")
    for mixcr_result_x in names:
        mixcr_to_csv(mixcr_result_x, 'igh')

    print("total time used %f" % (time.time() - start))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("work_dir")
    args = parser.parse_args()
    print(args.work_dir)

    work_dir = args.work_dir

    extract_and_move_file_in_dir(work_dir)
