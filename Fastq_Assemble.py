import os
import multiprocessing
import shutil
import argparse
import logging
import sys
import re
import json
import pandas as pd
from Bio import SeqIO

import Fastq_QC
import run_docker
from mkdir import mkdir
from plot_gc_cover import fasta2gcCover
from plot_gc_cover import plot_len_dis
from post_status import post_url
from post_status import post_pid
from post_status import copy_file
from post_status import write_status

def run_spades(fq1, fq2, outdir, threads):
    logging.info('spades')
    spades_out, scaffolds_fasta = run_docker.spades_docker(fq1, fq2, outdir, threads=threads)
    logging.info(spades_out)
    logging.info(scaffolds_fasta)
    return scaffolds_fasta

def filter_fa_by_len(fa, outfile, cutoff=500, name='scaffolds'):
    filter_fa = []
    max_len = 0
    seq_count = 0
    name_no = 1
    try:
        for seq in SeqIO.parse(fa, 'fasta'):
            tmp_len = len(seq.seq)
            if tmp_len > max_len:
                max_len = tmp_len
            if tmp_len >cutoff:
                seq.id = name + '_' + str(name_no)
                seq.description = seq.id + ':' + str(len(seq.seq))
                filter_fa.append(seq)
                name_no += 1
        if filter_fa:
            seq_count = SeqIO.write(filter_fa, outfile, 'fasta')
    except Exception as e:
        logging.error(f'filter_fa_by_len {e}')
    return seq_count, max_len

def run_bwa_align(fq1, fq2, scaffolds_fasta, outdir, threads=1):
    logging.info('bwa align')
    cmd_out = run_docker.bwa_align_docker(fq1, fq2, scaffolds_fasta, outdir, threads)
    if cmd_out[0]:
        logging.error(cmd_out)
        return cmd_out, 0, 0
    bwa_align_out, out_bam_file, out_genomecov_file, out_flagstat_file = cmd_out
    logging.info(bwa_align_out)
    gc_cover_png = os.path.join(outdir, 'align.genomecov.depth.png')
    try:
        logging.info('plot gc&depth')
        genome_info_file = fasta2gcCover(scaffolds_fasta, out_genomecov_file, outfile=gc_cover_png)
    except Exception as e:
        logging.error(f'run_bwa_align {e}')
    return bwa_align_out, out_bam_file, gc_cover_png, genome_info_file

def run_picard_CIAM(sortbam, outdir):
    logging.info('picard CollectInsertSizeMetrics')
    cmd_out = run_docker.picard_CIAM_docker(sortbam, outdir)
    if cmd_out[0]:
        logging.error(cmd_out)
        return cmd_out, 0, 0
    picard_out, outTxt, outPdf = cmd_out
    outpng = outPdf + '.png'
    try:
        check = 0
        size_lst = []
        with open(outTxt, 'rt') as h:
            for i in h.readlines():
                if check:
                    if re.match(r'\d+', i):
                        tmp = i.strip().split()
                        size_lst.extend([int(tmp[0])] * int(tmp[1]))
                if re.match('insert_size', i):
                    check = 1
        plot_len_dis(size_lst, outpng, xlabel='Insert Length')
    except Exception as e:
        logging.error(e)
    logging.info(picard_out)
    return picard_out, outTxt, outpng

def run_checkm(fa, outdir, threads):
    try:
        scaffolds_fasta_dir = os.path.join(outdir, 'scaffolds_fasta')
        mkdir(scaffolds_fasta_dir)
        shutil.copy(fa, scaffolds_fasta_dir)
        logging.info('checkm lineage_wf')
        checkm_lineage_out, checkm_outfile = run_docker.checkm_lineage_docker(scaffolds_fasta_dir, outdir, threads)
        logging.info(checkm_lineage_out)

        content = []
        checkm_outfile_csv = os.path.join(outdir, 'checkm_out.csv')
        if checkm_lineage_out == 0:
            with open(checkm_outfile, 'rt') as h:
                for i in h.readlines():
                    if not re.match('^----', i):
                        content.append(re.split(r' [ ]+', i.strip()))
        else:
            return checkm_lineage_out, 0
        df = pd.DataFrame(content).T.set_index(0).T
        df.to_csv(checkm_outfile_csv, index=False)
        return checkm_lineage_out, checkm_outfile_csv
    except Exception as e:
        logging.error(f'run_checkm {e}')

def Fastq_Assemble(fq1, fq2, outdir, threads, sampleTag='test', fastqc=None):
    mkdir(outdir)
    out_file_list = {}
    out_file_list['json'] = {}
    out_file_list['file'] = []
    # if fastqc:
    fq_cor_1, fq_cor_2 = fq1, fq2
    # else:
    #     # [fq_cor_1, fq_cor_2], qc_out_file = Fastq_QC.Fastq_QC(fq1, fq2, outdir, threads=threads)
    #     # out_file_list['json'].update(qc_out_file)
    #     [fq_cor_1, fq_cor_2], out_file_list_tmp = Fastq_QC.Fastq_QC(fq1, fq2, outdir, threads=threads)
    #     fastqc_json = os.path.join(outdir, 'Fastqc.json')
    #     with open(fastqc_json, 'w') as H:
    #         json.dump(out_file_list_tmp['json'], H, indent=2)
    #     out_file_list['file'].extend(out_file_list_tmp['file'])
    #     out_file_list['file'].append(fastqc_json)
    # fq_cor_1, fq_cor_2 = (fq1, fq2)
    # fq_cor_1 = os.path.join(outdir, 'musket_dir', 'trime_corrected.1.fastq') #test
    # fq_cor_2 = os.path.join(outdir, 'musket_dir', 'trime_corrected.2.fastq') #test
    if fq_cor_1 and fq_cor_2:
        scaffolds_fasta = run_spades(fq_cor_1, fq_cor_2, outdir, threads=threads)
        # scaffolds_fasta = os.path.join(outdir, 'spades_21_33_55', 'scaffolds.fasta') #test
        if scaffolds_fasta:
            scaffolds_fasta_file = os.path.join(outdir, sampleTag+'.scaffolds.fasta')
            scaffolds_count, scaffolds_max_len = filter_fa_by_len(scaffolds_fasta, scaffolds_fasta_file, name=sampleTag)
            logging.info(f'scaffolds {scaffolds_count} maxLen {scaffolds_max_len}')
            if scaffolds_count and scaffolds_max_len:
                try:
                    checkm_out, checkm_out_file = run_checkm(scaffolds_fasta_file, outdir, threads)
                    if checkm_out_file:
                        df_tmp = pd.read_csv(checkm_out_file).astype(str)
                        df_tmp.columns = df_tmp.columns.astype(str).str.lower().str.replace(' ','_').str.replace('.','_', regex=False)
                        df_tmp.columns = df_tmp.columns.str.replace('+','', regex=False).str.replace('#','', regex=False)
                        df_tmp.columns = 'c_' + df_tmp.columns
                        tmp_dict = {'genome_checkm_evaluate':df_tmp.loc[0].to_dict()}
                        out_file_list['json'].update(tmp_dict)
                        out_file_list['file'].append(checkm_out_file)
                except Exception as e:
                    logging.error(f'Fastq_Assemble checkM {e}')
                try:
                    bwa_align_out, out_bam_file, gc_cover_png, genome_info_file = run_bwa_align(fq1, fq2, scaffolds_fasta_file, outdir, threads)
                    if bwa_align_out == 0:
                        out_file_list['file'].append(gc_cover_png)
                        out_file_list['file'].append(genome_info_file)
                        df_tmp = pd.read_csv(genome_info_file, header=None, index_col=0)[1].astype(str).to_dict()
                        out_file_list['json'].update({'genome_info_summary':df_tmp})
                        out_file_list['json'].update({'gc_cover_png':os.path.split(gc_cover_png)[1]})
                        try:
                            picard_CIAM_out, inser_size, insert_size_pdf = run_picard_CIAM(out_bam_file, outdir)
                            out_file_list['file'].append(insert_size_pdf)
                            out_file_list['json'].update({'insert_len_size_pdf': os.path.split(insert_size_pdf)[1]})
                        except Exception as e:
                            logging.error(f'Fastq_Assemble picard {e}')
                except Exception as e:
                    logging.error(f'Fastq_Assemble bwa {e}')
                out_file_list['file'].append(scaffolds_fasta_file)
                out_file_list['json'].update({'scaffolds_fasta_file':os.path.split(scaffolds_fasta_file)[1]})
            else:
                logging.info('0 scaffolds')
        else:
            logging.info('assemble error')
    else:
        logging.info(f'Fastq_Assemble {fq1} or {fq2} error')
    return out_file_list

def exit_now(input):
    if len(input) != 2:
        sys.exit('input fastq must be R1 & R2')

if __name__ == '__main__':
    cpu_num = multiprocessing.cpu_count()
    if cpu_num > 4:
        cpu_num = cpu_num - 2
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', nargs='+', required=True, help='R1 fastq file, R2 fastq file')
    parser.add_argument('-o', '--outdir', default='./', help='output dir')
    parser.add_argument('-n', '--name', default='test', help='sample name')
    parser.add_argument('-t', '--threads', default=cpu_num, help='threads for fastqc')
    parser.add_argument('-tID', '--taskID', default='', help='task ID for report status')
    parser.add_argument('-Nqc', '--Nqc', action='store_true')
    parser.add_argument('-debug', '--debug', action='store_true')
    args = parser.parse_args()

    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)
    taskID = args.taskID
    post_pid(taskID)
    try:
        exit_now(args.input)
        outdir = args.outdir
        mkdir(outdir)
        logfile = os.path.join(outdir, 'log')
        logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        tmp_dir = os.path.join(outdir, 'tmp_dir')
        mkdir(tmp_dir)
        status_report = os.path.join(outdir, 'status_report.txt')

        fq1 = args.input[0]
        fq2 = args.input[1]
        if not args.Nqc:
            try:
                ## Fastqc
                s = f'Fastqc\tR\t'
                write_status(status_report, s)
                [fq_cor_1, fq_cor_2], out_file_list_tmp = Fastq_QC.Fastq_QC(args.input[0], args.input[1], tmp_dir,threads=args.threads)
                copy_file(out_file_list_tmp['file'], outdir)
                with open(os.path.join(outdir, 'Fastqc.json'), 'w') as H:
                    json.dump(out_file_list_tmp['json'], H, indent=2)
                s = f'Fastqc\tD\t'
                fq1 = fq_cor_1
                fq2 = fq_cor_2
            except Exception as e:
                logging.error(f'Fastqc {e}')
                s = f'Fastqc\tE\t'
            try:
                write_status(status_report, s)
            except Exception as e:
                logging.error(f'Fastqc status {e}')

        ## Assembly
        scaffolds_fasta_file = ''
        try:
            s = f'Assembly\tR\t'
            write_status(status_report, s)
            out_file_list_tmp = Fastq_Assemble(fq1, fq2, tmp_dir, args.threads, sampleTag=args.name)
            scaffolds_fasta_file = out_file_list_tmp['file'][-1]
            copy_file(out_file_list_tmp['file'], outdir)
            with open(os.path.join(outdir, 'Assembly.json'), 'w') as H:
                json.dump(out_file_list_tmp['json'], H, indent=2)
            s = f'Assembly\tD\t'
        except Exception as e:
            logging.error(f'Assembly {e}')
            s = f'Assembly\tE\t'
        try:
            write_status(status_report, s)
            try:
                post_url(taskID, '2', 'http://localhost/task/getTaskRunningStatus/')
            except Exception as e:
                logging.error(f'post_url getTaskRunningStatus {e}')
            # post_url(taskID, 'Assembly')
        except Exception as e:
            logging.error(f'Assembly status {e}')

        if not args.debug:
            try:
                shutil.rmtree(tmp_dir)
            except Exception as e:
                logging.error(e)

        # for i in out_file_list:
        #     if os.path.isfile(i) or os.path.isdir(i):
        #         try:
        #             des_file = shutil.move(i, args.outdir)
        #             logging.info(des_file)
        #         except Exception as e:
        #             logging.error(e)
        #     else:
        #         logging.error(f' not exitst {i}')
    except Exception as e:
        logging.error(e)

