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

# def run_spades(pe_fq_list, se_fq, outdir, threads):
def run_spades(outdir, pe_fq, se_fq, pacbio_clr, nanopore, threads):
    logging.info('spades')
    spades_out, scaffolds_fasta = run_docker.spades_docker(outdir, pe_fq, se_fq, pacbio_clr, nanopore, threads=threads)
    logging.info(spades_out)
    logging.info(scaffolds_fasta)
    return scaffolds_fasta

def run_canu(outdir, pacbio, nanopore):
    logging.info('canu')
    spades_out, scaffolds_fasta = run_docker.canu_docker(outdir, pacbio, nanopore)
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

def run_bwa_align(fq_list, scaffolds_fasta, outdir, threads=1, nanopore=None, pacbio=None):
    logging.info('bwa align')
    cmd_out = run_docker.bwa_align_docker(fq_list, scaffolds_fasta, outdir, threads, nanopore=nanopore, pacbio=pacbio)
    if cmd_out[0]:
        logging.error(cmd_out)
        return cmd_out, 0, 0, 0
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

def Fastq_Assemble(fq_list, outdir, threads, sampleTag='test', pacbio=None, nanopore=None, fastqc=None):
    mkdir(outdir)
    out_file_list = {}
    out_file_list['json'] = {}
    out_file_list['file'] = []
    scaffolds_fasta = ''
    # if fastqc:
    # fq_cor_1, fq_cor_2 = None, None
    # if len(fq_list) == 2:
    #     fq_cor_1, fq_cor_2 = fq_list
    # else:
    #     fq_cor_1 = fq_list[0]
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

    pe_fq, se_fq = [], []
    nanopore_tmp = None 
    pacbio_tmp = None 
    scaffolds_fasta = None
    if fq_list and fq_list[0]:
        for fqs in fq_list:
            if len(fqs) == 2:
                pe_fq.append(fqs)
            elif len(fqs) == 1:
                se_fq.extend(fqs)
            else:
                logging.error(f'fastq file number not correct: {fqs}')

        # if pacbio_ccs:
        #     se_fq.extend(pacbio_ccs)
        try:
            scaffolds_fasta = run_spades(outdir, pe_fq, se_fq, pacbio, nanopore, threads=threads)
            # logging.info(f'Assebly for {pe_fq} | {se_fq} | {pacbio} | {nanopore}')
        except Exception as e:
            logging.error(f'Fastq_Assemble Assembly {e}')
    elif (not fq_list) and nanopore or pacbio:
        if nanopore:
            nanopore_tmp = nanopore[0]
        if pacbio:
            pacbio_tmp = pacbio[0]
        try:
            scaffolds_fasta = run_canu(outdir, pacbio_tmp, nanopore_tmp)
        except Exception as e:
            logging.error(f'canu Assemble Assembly {e}')        

    # if fq_cor_1 and fq_cor_2 and pacbio_ccs:
    #     scaffolds_fasta = run_spades(fq_list, pacbio_fq, outdir, threads=threads)
    # elif fq_cor_1 and fq_cor_2:
    #     scaffolds_fasta = run_spades(fq_list, None, outdir, threads=threads)
    # elif fq_cor_1:
    #     scaffolds_fasta = run_spades(None, fq_cor_1, outdir, threads=threads)
    #     # scaffolds_fasta = os.path.join(outdir, 'spades_21_33_55', 'scaffolds.fasta') #test
    # else:
    #     logging.info(f'Assebly for {fq_list} & {pacbio_fq}')
    if os.path.isfile(scaffolds_fasta):
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
                if fq_list and fq_list[0]:
                    bwa_align_out, out_bam_file, gc_cover_png, genome_info_file = run_bwa_align(fq_list[0], scaffolds_fasta_file, outdir, threads)
                elif nanopore or pacbio:
                    bwa_align_out, out_bam_file, gc_cover_png, genome_info_file = run_bwa_align(None, scaffolds_fasta_file, outdir, threads, nanopore_tmp, pacbio_tmp)
                if bwa_align_out == 0:
                    out_file_list['file'].append(gc_cover_png)
                    out_file_list['file'].append(genome_info_file)
                    df_tmp = pd.read_csv(genome_info_file, header=None, index_col=0)[1].astype(str).to_dict()
                    out_file_list['json'].update({'genome_info_summary':df_tmp})
                    out_file_list['json'].update({'gc_cover_png':os.path.split(gc_cover_png)[1]})
                    if fq_list and (len(fq_list[0]) == 2):
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
    parser.add_argument('-i', '--input', nargs='+', action='append', help='R1 fastq file, R2 fastq file')
    parser.add_argument('-pacbio', '--pacbio', action='append', help='pacbio fastq file')
    # parser.add_argument('-pacbioccs', '--pacbioccs', action='append', help='pacbio ccs fastq file')
    parser.add_argument('-nanopore', '--nanopore', action='append', help='nanopore fastq file')
    parser.add_argument('-o', '--outdir', default='./', help='output dir')
    parser.add_argument('-n', '--name', default='test', help='sample name')
    parser.add_argument('-t', '--threads', default=cpu_num, help='threads for fastqc')
    parser.add_argument('-tID', '--taskID', default='', help='task ID for report status')
    parser.add_argument('-Nqc', '--Nqc', action='store_true')
    parser.add_argument('-debug', '--debug', action='store_true')
    args = parser.parse_args()

    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)
    outdir = args.outdir
    mkdir(outdir)
    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    tmp_dir = os.path.join(outdir, 'tmp_dir')
    mkdir(tmp_dir)
    status_report = os.path.join(outdir, 'status_report.txt')
    taskID = args.taskID
    post_pid(taskID)
    try:
        # fq1 = args.input[0]
        # fq2 = args.input[1]
        fq_list_total = args.input
        fq_cor_list_total = []

        if not args.Nqc:
            s = f'Fastqc\tR\t'
            write_status(status_report, s)
            # for fq_list in fq_list_total:

            try:
                ## Fastqc
                if fq_list_total and fq_list_total[0]:
                    fq_list = fq_list_total[0]
                    fq_cor_list, out_file_list_tmp = Fastq_QC.Fastq_QC(fq_list, tmp_dir, threads=args.threads)
                    fq_cor_list_total.append(fq_cor_list)
                    fq_cor_list_total.extend(fq_list_total[1:])
                else:
                    fq_cor_list = []
                    out_file_list_tmp = {}
                    out_file_list_tmp['json'] = {}
                    out_file_list_tmp['file'] = {}
                with open(os.path.join(outdir, 'Fastqc.json'), 'w') as H:
                    json.dump(out_file_list_tmp['json'], H, indent=2)
                s = f'Fastqc\tD\t'
                # fq1 = fq_cor_1
                # fq2 = fq_cor_2
                # fq_list = fq_cor_list
            except Exception as e:
                logging.error(f'Fastqc {e}')
                s = f'Fastqc\tE\t'
            try:
                copy_file(out_file_list_tmp['file'], outdir)
            except Exception as e:
                logging.error(f'copyfile {e}')
            try:
                post_url(taskID, 'Fastqc')
                write_status(status_report, s)
            except Exception as e:
                logging.error(f'Fastqc status {e}')

        ## Assembly
        scaffolds_fasta_file = ''
        try:
            s = f'Assembly\tR\t'
            write_status(status_report, s)
            # out_file_list_tmp = Fastq_Assemble(fq_list, tmp_dir, args.threads, sampleTag=args.name)
            logging.info(f'GII data: {fq_cor_list_total}')
            logging.info(f'GIII pacbio: {args.pacbio}')
            logging.info(f'GIII nanopore: {args.nanopore}')

            out_file_list_tmp = Fastq_Assemble(fq_cor_list_total, tmp_dir, args.threads, sampleTag=args.name, pacbio=args.pacbio, nanopore=args.nanopore)
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
                post_url(taskID, 'Assembly')
                post_url(taskID, '2', 'http://localhost/task/getTaskRunningStatus/')
            except Exception as e:
                logging.error(f'post_url getTaskRunningStatus {e}')
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

