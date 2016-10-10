#
# Author: Hao Wu <echowuhao@gmail.com>
#
# create barcode file in the dir
# python barcode.py /home/hwu/analysis/109-Stimulation_61515a_Proliferating_T_and_B_cells_dupli_1_8-18-2015_149_119
# will create a barcode.json file in the dir /home/hwu/analysis/109-Stimulation_61515a_Proliferating_T_and_B_cells_dupli_1_8-18-2015_149_119
#
import os, glob
import os, sys
import json
import logging
from os.path import join, dirname, basename
import time
import csv

from extract_mixcr import extract_and_move_file_in_dir
from run_mixcr import run_mixcr_for_dir


def bam_to_fastq(bam_and_fastq_name):
    command = "/usr/local/bin/bamToFastq -i %s  -fq %s" % bam_and_fastq_name
    os.system(command)


#TODO: produce barcode name
def rawname_to_fastq_name(bam_name, work_dir):
    serial = bam_name.split('_')[1]
    barcode_dict = {}
    with open(os.path.join(work_dir, 'barcode.json')) as f:
        barcode_dict = json.load(f)
    barcode_name = barcode_dict.get(serial, False)
    if barcode_name:
        return barcode_dict.get(serial) + '.fastq'
    else:
        print("!!!! ========================== !!!")
        print("!!!! %s do not have barcode !!!" % (serial))
        print("!!!! ========================== !!!")


def transform_bam_to_fastq(work_dir):
    raw_bam_dir = os.path.basename(work_dir)

    # we striped the Auto_user part in post_sync.py
    if os.path.basename(work_dir)[0].isdigit():
        raw_bam_dir = 'Auto_user_SN2-' + raw_bam_dir

    raw_bam_path = os.path.join(RAW_DATA_ROOT, raw_bam_dir)

    print("!! Generating fastq files for %s" % raw_bam_dir)
    start = time.time()
    # Make sure the output dir is created
    bam = os.path.join(raw_bam_path, 'download_links', "IonXpress*.bam")
    # remove basecaller file
    names = [name for name in glob.glob(bam) if "basecaller" not in name]
    for name in names:
        fastq_name = rawname_to_fastq_name(basename(name), work_dir)
        if fastq_name:
            fastq_path = join(ANALYSIS_ROOT, work_dir, 'raw_fastq', fastq_name)
            # check if fastq file already exists
            if not os.path.exists(fastq_path):
                bam_and_fastq_path = (name, fastq_path)
                print(fastq_name)
                bam_to_fastq(bam_and_fastq_path)
    print("total time used %f" % (time.time() - start))


def analysis_dir_to_raw_data_dir(analysis_dir):
    """we striped the Auto_user part in post_sync.py
    so add it back
    """
    raw_data_dir = os.path.basename(analysis_dir)
    if os.path.basename(raw_data_dir)[0].isdigit():
        raw_data_dir = 'Auto_user_SN2-' + raw_data_dir
    return raw_data_dir


def create_barcode_json(work_dir):
    expmeta_json_file = 'ion_params_00.json'
    expmeta_path = os.path.join(RAW_DATA_ROOT,
                                analysis_dir_to_raw_data_dir(work_dir),
                                expmeta_json_file)
    #read expmeta json
    barcode_dict = {}
    with open(expmeta_path) as f:
        barcode_dict = json.load(f)
    barcode_info = barcode_dict['barcodeInfo']

    barcode_json = {}

    for barcode in barcode_info:
        sample_name = barcode_info[barcode]['sample']
        if sample_name != 'none':
            barcode_json[barcode[-3:]] = sample_name

    with open(os.path.join(work_dir, 'barcode.json'), 'w') as f:
        f.write(json.dumps(barcode_json, sort_keys=True, indent=4))


if __name__ == '__main__':

    RAW_DATA_ROOT = '/home/hwu/raw_data'
    ANALYSIS_ROOT = '/home/hwu/analysis'
    LOG_ROOT = '/home/hwu/analysis_logs'
    LOG_FILENAME = os.path.join(LOG_ROOT, 'upstream.log')
    logging.basicConfig(filename=LOG_FILENAME,
                        level=logging.DEBUG,
                        stream=sys.stdout,
                        format='%(asctime)s %(message)s')

    logging.getLogger().addHandler(logging.StreamHandler())
    logging.info('%s start', __file__)

    # TODO: better doc
    # this file scan the orignal data dir 
    # prepare dir for anlaysis and crate a link to the original dir

    # raw data if store in
    #/home/hwu/raw_data/Auto_user_SN2-48-Template-TCR_sequencing_84_055/plugin_out/FastaCreator_out.79

    # the output is in

    #/home/hwu/analysis/SN2-48-Template-TCR_sequencing_84_055/

    # basically, remove the directory and remove the redudancy in file and dir names

    counter = 0
    new_dirs = []
    for root, dirs, files in os.walk(RAW_DATA_ROOT):
        # when top is FastaCreat_out*
        # we are in the dir of fasta file
        if basename(root).startswith('FastaCreator_out'):
            counter = counter + 1
            # /home/hwu/raw_data/Auto_user_SN2-48-Template-TCR_sequencing_84_055/plugin_out/FastaCreator_out.79
            # Auto_user_SN2-48-Template-TCR_sequencing_84_055
            parent_dir = basename(dirname(dirname(root))).lstrip(
                'Auto_user_SN2-')
            abs_dir = join(ANALYSIS_ROOT, parent_dir)
            # make new dir
            if not os.path.exists(abs_dir):
                new_dirs.append(abs_dir)
                os.makedirs(abs_dir)
                logging.info("%s created", abs_dir)
            # create a link to original dir
            raw_data_dir = join(abs_dir, 'raw_fasta')

            if not os.path.exists(raw_data_dir):
                os.symlink(root, raw_data_dir)

            for analysis_dir in ['raw_fastq',
                                 'mixcr_igh',
                                 'mixcr_igh_extracted',
                                 'mixcr_igh_export_align',
                                 'shm',
                                 'mixcr_trb',
                                 'mixcr_trb_extracted', ]:
                abs_analysis_dir = join(abs_dir, analysis_dir)
                if not os.path.exists(abs_analysis_dir):
                    os.makedirs(abs_analysis_dir)

    logging.info("Total number of fasta dir is %s", counter)
    logging.info("number of new fasta dir is %s", len(new_dirs))
    for new_created_dir in new_dirs:
        logging.info("creating barcode.json file in %s", new_created_dir)
        create_barcode_json(new_created_dir)
        transform_bam_to_fastq(new_created_dir)
        run_mixcr_for_dir(new_created_dir)
        extract_and_move_file_in_dir(new_created_dir)
