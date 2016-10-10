#!/usr/bin/python
import os, sys
import time
import glob
import argparse
import hcvb_lab


def run_mixcr_bcr(input_file_name):
    mixcr = '/home/hwu/app/mixcr-1.6/mixcr'
    output_file_name = input_file_name.replace(
        "raw_fastq", "mixcr_igh").replace('.fastq', '.txt')
    align_file_name = output_file_name.replace('.txt', '.vdjca')
    assemle_file_name = output_file_name.replace('.txt', '.clns')
    align_command = "%s align --chains IGH  %s %s;" % (mixcr, input_file_name,
                                                       align_file_name)
    assemle_command = "%s assemble  %s %s;" % (mixcr, align_file_name,
                                               assemle_file_name)
    export_command = " %s exportClones -vHit -jHit -dHit -count -aaFeature CDR3 %s %s;" % (
        mixcr, assemle_file_name, output_file_name)
    command = align_command + assemle_command + export_command
    print(input_file_name)
    os.system(command)
    #print(command)


def run_mixcr_tcr(input_file_name):
    mixcr = '/home/hwu/app/mixcr-1.6/mixcr'
    output_file_name = input_file_name.replace(
        "raw_fastq", "mixcr_trb").replace('.fastq', '.txt')
    align_file_name = output_file_name.replace('.txt', '.vdjca')
    assemle_file_name = output_file_name.replace('.txt', '.clns')
    align_command = "%s align --chains TRB  %s %s;" % (mixcr, input_file_name,
                                                       align_file_name)
    assemle_command = "%s assemble  %s %s;" % (mixcr, align_file_name,
                                               assemle_file_name)
    export_command = " %s exportClones -vHit -jHit -dHit -count -aaFeature CDR3 %s %s;" % (
        mixcr, assemle_file_name, output_file_name)
    command = align_command + assemle_command + export_command
    print(input_file_name)
    os.system(command)


def run_mixcr_for_dir(work_dir):
    # Make sure the output dir is created
    print("mixcr start")
    start = time.time()
    fastq = os.path.join(work_dir, "raw_fastq/*.fastq")
    names = glob.glob(fastq)
    for name in names:
        barcode_name = os.path.basename(name)
        if hcvb_lab.is_bcr(barcode_name):
            run_mixcr_bcr(name)
        else:
            run_mixcr_tcr(name)

    print("mixcr done")
    print("total time used %f" % (time.time() - start))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("work_dir")
    args = parser.parse_args()
    print(args.work_dir)
    work_dir = args.work_dir
    run_mixcr_for_dir(work_dir)
