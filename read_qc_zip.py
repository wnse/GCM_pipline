import os
import re
import json
import logging
import argparse
from zipfile import ZipFile
from zipfile import is_zipfile
import pandas as pd


def unzip_file(f, des_dir='./'):
    if is_zipfile(f):
        # logging.info(f)
        with ZipFile(f) as myzip:
            # for i in myzip.namelist():
            #     print(i)
            try:
                myzip.extractall(path=des_dir)
            except Exception as e:
                logging.error(e)

def read_summary(f):
    if os.path.isfile(f):
        try:
            df = pd.read_csv(f, sep='\t', header=None)
            df = df.set_index(1)
            df.index = df.index.str.lower().str.replace(' ','_').str.replace('.','_', regex=False)
            return df[0].to_dict()
        except Exception as e:
            logging.error(f'read_summary {e}')

def get_qc_file_path(qc_dir):
    qc_path = {}
    for d_name in os.listdir(qc_dir):
        if os.path.isdir(os.path.join(qc_dir, d_name)) and re.match(r'\S+_fastqc', d_name):
            sample = re.match(r'(\S+)_fastqc', d_name).group(1)
            tmp_qc_path = {}
            tmp_qc_path['name'] = sample
            tmp_qc_path['path'] = {}
            for item in os.scandir(os.path.join(qc_dir, d_name)):
                if item.is_file():
                    key = item.name.lower().replace(' ','_').replace('.','_')
                    tmp_qc_path['path'][key] = os.path.join(os.path.split(qc_dir)[1], d_name, item.name)
                    if item.name == 'summary.txt':
                        tmp_qc_path['summary'] = read_summary(os.path.join(qc_dir, d_name, item.name))
                if item.is_dir():
                    for tmp_file in os.listdir(os.path.join(qc_dir, d_name, item.name)):
                        key = tmp_file.lower().replace(' ','_').replace('.','_')
                        tmp_qc_path['path'][key] = os.path.join(os.path.split(qc_dir)[1], d_name, item.name, tmp_file)
            qc_path[sample] = (tmp_qc_path)
    return qc_path

def get_qc_data(qc_zip_dir, unzip_dir):
    if os.path.isdir(qc_zip_dir):
        for f in os.listdir(qc_zip_dir):
            unzip_file(os.path.join(qc_zip_dir, f), des_dir=unzip_dir)
        if os.path.isdir(unzip_dir):
            qc_data = get_qc_file_path(unzip_dir)
            return qc_data
    return None


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s',datefmt='%Y-%m-%d %H:%M:%S')
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-d', '--qc_zip_dir', required=True, help='QC dir')
    parser.add_argument('-o', '--qc_unzip_dir', default='./', help='unzip qc dir')
    args = parser.parse_args()
    #
    # f = './query5_1_fastqc.zip'
    # path = '../test_qc'
    # unzip_file(f, des_dir=path)

    # f = '../test_qc/query5_1_fastqc/summary.txt'
    # print(read_summary(f))
    # tmp_dict = get_png_path(path)
    # print(json.dumps(tmp_dict, indent=2))
    qc_data = get_qc_data(args.qc_zip_dir, args.qc_unzip_dir)
    print(json.dumps(qc_data, indent=2))
