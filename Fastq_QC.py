import json
import multiprocessing
import os
import shutil
import argparse
import gzip
import logging
import re
from Bio import SeqIO
import pandas as pd

from mkdir import mkdir
import run_docker
from unzip_file import file_type
from read_qc_zip import get_qc_data
from post_status import post_url
from post_status import copy_file
from post_status import write_status

def cal_q20(f):
    ft = file_type(f)
    if ft == 'FQ':
        open_file = open
    elif ft == 'gz':
        open_file = gzip.open
    else:
        return 0, 0, 0, 0
    total_reads = 0
    total_base = 0
    q20 = 0
    q30 = 0
    gc = 0
    with open_file(f, 'rt') as handle:
        for rec in SeqIO.parse(handle, 'fastq'):
            a = rec.letter_annotations['phred_quality']
            q20 += sum(i > 20 for i in a)
            q30 += sum(i > 30 for i in a)
            gc += sum(i in ['G','C', 'g', 'c'] for i in rec.seq)
            total_base += len(a)
            total_reads += 1
    return total_reads, total_base, q20/total_base, q30/total_base, gc/total_base

def cal_avg_len(f):
    ft = file_type(f)
    if ft == 'FQ':
        open_file = open
    elif ft == 'gz':
        open_file = gzip.open
    else:
        return 0
    total = 0
    reads = 0
    with open_file(f, 'rt') as handle:
        for rec in SeqIO.parse(handle, 'fastq'):
            total += len(rec.seq)
            reads += 1
    return total/reads

def run_cal_q20(fq_path_list, name_list=[]):
    # q20 = []
    # q30 = []
    # total_reads = []
    # total_base = []
    # gc = []
    total_info = {}
    logging.info('q20 start')
    for i, file in enumerate(fq_path_list):
        if name_list:
            name = name_list[i]
        else:
            name = get_file_name(file)
        if name:
            total_reads_tmp, total_base_tmp, q20_tmp, q30_tmp, gc_tmp = cal_q20(file)
        # q20.append(q20_tmp)
        # q30.append(q30_tmp)
        # total_reads.append(total_reads_tmp)
        # total_base.append(total_base_tmp)
        # gc.append(gc_tmp)
            avg_tmp = total_base_tmp / total_reads_tmp
            total_info[name] = dict(zip(['total_reads', 'total_bases', 'Q20', 'Q30', 'GC', 'average_length'],
                                        [str(total_reads_tmp), str(total_base_tmp), str(q20_tmp), str(q30_tmp), str(gc_tmp), str(avg_tmp)]))
    logging.info('q20 end')
    # return total_reads, total_base, q20, q30, gc
    return total_info

def run_fastqc(fq_path_list, outdir, threads):
    logging.info('fastqc')
    qc_out, qc_dir = run_docker.fastqc_docker(fq_path_list, outdir, threads)
    logging.info(f'{qc_out}')
    return qc_out, qc_dir

def run_trim_pe(fq1, fq2, outdir, threads):
    logging.info('trimmomatic')
    trime_cmd_out, fq_trime_1, fq_trime_2 = run_docker.trim_docker_pe(fq1, fq2, outdir, threads)
    logging.info(f'{trime_cmd_out}')
    return trime_cmd_out, fq_trime_1, fq_trime_2

def run_musket(fq1, fq2, outdir, threads):
    logging.info('musket')
    musket_cmd_out, fq_cor_1, fq_cor_2 = run_docker.musket_docker(fq1, fq2, outdir, threads)
    logging.info(f'{musket_cmd_out}')
    return musket_cmd_out, fq_cor_1, fq_cor_2

def run_jellyfish(fq, outdir, threads):
    logging.info('jefflyfish count')
    jellyfish_count_out, jf_file = run_docker.jellyfish_docker_count(fq, outdir, threads)
    logging.info(f'{jellyfish_count_out}')
    if jf_file:
        logging.info('jefflyfish histo')
        jellyfish_histo_out, histo_file = run_docker.jellyfish_docker_histo(jf_file, outdir, threads)
        logging.info(f'{jellyfish_histo_out}')
        return jellyfish_histo_out, histo_file
    return jellyfish_count_out, 0

def run_genomescope(histo_file, outdir, read_length):
    logging.info('genomescope')
    genomescope_cmd_out, summary_file, outpng = run_docker.genomescope_docker(histo_file, outdir, read_length=int(read_length))
    logging.info(f'{genomescope_cmd_out}')
    if summary_file:
        try:
            with open(summary_file, 'rt') as h:
                res = {}
                for l in h.readlines():
                    if len(res) > 0 or re.search(r'property', l):
                        tmp = (re.split(r' [ ]+', l.strip()))
                        res[tmp[0]] = tmp[1:]
            return genomescope_cmd_out, res, outpng
        except Exception as e:
            logging.error(f'run_genomescope {e}')
    return genomescope_cmd_out, 0, 0

def get_file_name(file):
    f = os.path.split(file)[1]
    if re.match('(\S+)\.fq', f):
        return re.match('(\S+)\.fq', f).group(1)
    elif re.match('(\S+)\.fastq', f):
        return re.match('(\S+)\.fastq', f).group(1)
    else:
        return None

def Fastq_QC(fq1, fq2, outdir, threads):
    mkdir(outdir)
    out_file_list = {}
    out_file_list['json'] = {}
    out_file_list['file'] = []
    fastqc_sta_clean = {}
    fastqc_sta_raw = {}
    fastqc_summary = {}
    fq_cor_1, fq_cor_2 = ('', '')
    qc_out, qc_dir = run_fastqc([fq1, fq2], outdir, threads)
    qc_out_dir = os.path.join(outdir, 'fastqc_unzip_dir')
    if qc_out:
        logging.error(f'fastqc error {qc_out}')
    try:
        fastqc_summary = get_qc_data(qc_dir, qc_out_dir)
        out_file_list['file'].append(qc_out_dir)
    except Exception as e:
        logging.error(e)
    try:
        fastqc_sta_raw = run_cal_q20([fq1, fq2])
    except Exception as e:
        logging.error(e)

    trime_out, fq_trime_1, fq_trime_2 = run_trim_pe(fq1, fq2, outdir, threads)
    if trime_out:
        logging.info(f'trime error {trime_out}')
    if fq_trime_1 and fq_trime_2:
        musket_out, fq_cor_1, fq_cor_2 = run_musket(fq_trime_1, fq_trime_2, outdir, threads)
        if musket_out:
            logging.info(f'musket error {musket_out}')
        if fq_cor_1 and fq_cor_2:
            jellyfish_out, histo_file = run_jellyfish(fq_cor_1, outdir, threads)
            read_length = cal_avg_len(fq_cor_1)
            fastqc_sta_clean = run_cal_q20([fq_cor_1, fq_cor_2], name_list=[get_file_name(fq1), get_file_name(fq2)])
            if jellyfish_out:
                logging.info(f'jellyfish error {jellyfish_out}')
        else:
            jellyfish_out, histo_file = run_jellyfish(fq_trime_1, outdir, threads)
            read_length = cal_avg_len(fq_trime_1)
            fastqc_sta_clean = run_cal_q20([fq_trime_1, fq_trime_2], name_list=[get_file_name(fq1), get_file_name(fq2)])
    else:
        jellyfish_out, histo_file = run_jellyfish(fq1, outdir, threads)
        read_length = cal_avg_len(fq1)

    if histo_file:
        logging.info(f'read length:  {read_length}')
        genomescope_out, genomescope_res , genomescope_png = run_genomescope(histo_file, outdir, read_length)
        if genomescope_res:
            df_genomescope = pd.DataFrame.from_dict(genomescope_res).set_index('property').T
            genomescope_out_file = os.path.join(outdir, 'qc_genome_size_predicted.csv')
            df_genomescope.to_csv(genomescope_out_file)
            out_file_list['file'].append(genomescope_out_file)
            out_file_list['file'].append(genomescope_png)
            df_genomescope['value'] = df_genomescope['min'].astype(str) + '-' + df_genomescope['max'].astype(str)
            df_genomescope.index = df_genomescope.index.str.lower().str.replace(' ','_').str.replace('.','_',regex=False)
            out_file_list['json'].update({'genomescope_kmer_png':os.path.split(genomescope_png)[1]})
            out_file_list['json'].update({'genomescope_predict':df_genomescope['value'].astype(str).to_dict()})

    fastq_sta = []
    for name, summary in fastqc_summary.items():
        tmp = summary
        if name in fastqc_sta_raw.keys():
            tmp['raw_sta'] = fastqc_sta_raw[name]
        if name in fastqc_sta_clean.keys():
            tmp['clean_sta'] = fastqc_sta_clean[name]
        fastq_sta.append(tmp)

    out_file_list['json'].update({'fastqc_info': fastq_sta})
    return [fq_cor_1, fq_cor_2], out_file_list


if __name__ == '__main__':
    cpu_num = multiprocessing.cpu_count()
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', nargs='+', required=True, help='R1 fastq file, R2 fastq file')
    parser.add_argument('-o', '--outdir', default='./', help='output dir')
    parser.add_argument('-n', '--name', default='test', help='sample name')
    parser.add_argument('-t', '--threads', default=cpu_num, help='threads for fastqc')
    parser.add_argument('-tID', '--taskID', default='', help='task ID for report status')
    args = parser.parse_args()
    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)

    outdir = args.outdir
    mkdir(outdir)
    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    status_report = os.path.join(outdir, 'status_report.txt')
    tmp_dir = os.path.join(outdir, 'tmp_dir')
    mkdir(tmp_dir)
    taskID = args.taskID
    fq_cor_1, fq_cor_2 = ('', '')
    try:
        ## Fastqc
        s = f'Fastqc\tR\t'
        write_status(status_report, s)
        [fq_cor_1, fq_cor_2], out_file_list_tmp = Fastq_QC(args.input[0], args.input[1], tmp_dir, threads=args.threads)
        copy_file(out_file_list_tmp['file'], outdir)
        with open(os.path.join(outdir, 'Fastqc.json'), 'w') as H:
            json.dump(out_file_list_tmp['json'], H, indent=2)
        s = f'Fastqc\tD\t'
    except Exception as e:
        logging.error(f'Fastqc {e}')
        s = f'Fastqc\tE\t'
    try:
        write_status(status_report, s)
        post_url(taskID, 'Fastqc')
    except Exception as e:
        logging.error(f'Fastqc status {e}')

    # try:
        # clean_fq_list, out_file_list = Fastq_QC(args.input[0], args.input[1], tmp_dir, args.threads)
        # if out_file_list['file']:
        #     copy_file(out_file_list['file'], outdir)
        # if out_file_list['json']:
        #     print(json.dumps(out_file_list['json'], indent=2))

        # for i in out_file_list['dir']:
        #     if os.path.isfile(i) or os.path.isdir(i):
        #         try:
        #             des_file = shutil.move(i, args.outdir)
        #             logging.info(des_file)
        #         except Exception as e:
        #             logging.error(e)
        #     else:
        #         logging.error(f' not exitst {i}')
    # except Exception as e:
    #     logging.error(e)
