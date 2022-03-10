import os
import sys
import shutil
import json
import logging
import argparse
import multiprocessing
import datetime
import pandas as pd

from mkdir import mkdir
from Fastq_QC import Fastq_QC
from Fastq_Assemble import Fastq_Assemble
from Fasta_Taxonomy import Fasta_Taxonomy
from Fasta_FactorAnno import FactorsAnno
from Fasta_Typing import Fasta_Typing
from Fasta_Annotation import FunctionAnno
from post_status import post_url
from post_status import post_pid
from post_status import copy_file
from post_status import write_status

def exit_now_fqPE(input):
    if len(input) != 2:
        sys.exit('input fastq must be R1 & R2')

def exit_now_fa(input):
    if len(input) != 1:
        sys.exit('input fasta should be 1')



if __name__ == '__main__':
    cpu_num = multiprocessing.cpu_count()
    if cpu_num > 4:
        cpu_num = cpu_num - 2
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', nargs='+', required=True, help='R1 fastq file, R2 fastq file')
    parser.add_argument('-o', '--outdir', default='./', help='output dir')
    parser.add_argument('-n', '--name', default='test', help='sample name')
    parser.add_argument('-t', '--threads', default=cpu_num, help='threads for analysis')
    parser.add_argument('-tID', '--taskID', default='', help='task ID for report status')
    parser.add_argument('-type', '--type', default='fastqPE', choices=['fqPE', 'fa'], help='data type for analysis')
    parser.add_argument('-d16S', '--db_16S', default='/Bio/tax-20200810_DB/database/best.16s',
                        help='16S database, pre blastn index')
    parser.add_argument('-i16S', '--info_16S', default='/Bio/tax-20200810_DB/database/16sdb.info_temp',
                        help='16S database taxonomy info')
    parser.add_argument('-dg', '--db_genome', default='/Bio/tax-20200810_DB/database/genome.msh',
                        help='genome database, mash index')
    parser.add_argument('-ig', '--info_genome', default='/Bio/tax-20200810_DB/database/all.genome.best.info.final',
                        help='genome database taxonomy info')
    parser.add_argument('-dgf', '--db_genome_fa', default='/Bio/tax-20200810_DB/fasta',
                        help='genome database fasta file dir')
    parser.add_argument('-con_db_path', '--con_db_path', default=f'{os.path.join(bin_dir, "conf_DB_path.json")}',
                        help='conf_DB_path.json genome database path info')
    parser.add_argument('-cgmlst_db_path', '--cgmlst_db_path', default=f'{os.path.join(bin_dir, "cgMLST_db.csv")}',
                        help='cgMLST database path for schema_seed, training file and core gene')
    parser.add_argument('-db_list', '--db_list', default=['kegg','cog', 'NR', 'MetaCyc', 'cazy', 'PHI', 'SwissProt'],
                        help='database for annotation', nargs='+')
    parser.add_argument('-debug', '--debug', action='store_true')
    args = parser.parse_args()

    # if args.type == 'fastqPE':
    #     exit_now_fqPE(args.input)
    # else:
    #     exit_now_fa(args.input)

    if os.path.isfile(args.con_db_path):
        with open(args.con_db_path, 'rt') as h:
            con_db_path = json.load(h)

    db_anno_list = args.db_list
    outdir = args.outdir
    taskID = args.taskID
    status_report = os.path.join(outdir, 'status_report.txt')
    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)
    tmp_dir = None
    scaffolds_fasta_file = None

    try:
        mkdir(outdir)
        logfile = os.path.join(outdir, 'log')
        logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        logging.info(f'{vars(args)}')
        logging.info('Analysis Start')
        tmp_dir = os.path.join(outdir, 'tmp_dir')
        mkdir(tmp_dir)
    except Exception as e:
        logging.error(e)
        sys.exit(e)

    post_pid(taskID)
    fq_cor_list = []
    if args.type == 'fastqPE':
        fq_cor_1, fq_cor_2 = ('', '')
        try:
            ## Fastqc
            s = f'Fastqc\tR\t'
            write_status(status_report, s)
            fq_cor_list, out_file_list_tmp = Fastq_QC(args.input, tmp_dir, threads=args.threads)
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

        ## Assembly
        scaffolds_fasta_file = ''
        if fq_cor_list:
            logging.info(f'{fq_cor_list}')
        else:
            fq_cor_list = args.input
        try:
            s = f'Assembly\tR\t'
            write_status(status_report, s)
            out_file_list_tmp = Fastq_Assemble(fq_cor_list, tmp_dir, args.threads, sampleTag=args.name, fastqc=1)
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
            post_url(taskID, 'Assembly')
        except Exception as e:
            logging.error(f'Assembly status {e}')

    if args.type == 'fa':
        scaffolds_fasta_file = args.input[0]
    species = None
    faa_file = None
    if scaffolds_fasta_file:
        try:
            ## Taxonomy
            s = f'Taxonomy\tR\t'
            write_status(status_report, s)
            out_file_list_tmp = Fasta_Taxonomy(scaffolds_fasta_file, tmp_dir, args.db_16S, args.info_16S, args.db_genome, args.info_genome, args.db_genome_fa, args.threads, args.name)
            species_file = out_file_list_tmp['file'][-1]
            species = ' '.join(pd.read_csv(species_file, index_col=0, header=None).loc[['family','species'],1].to_list())
            copy_file(out_file_list_tmp['file'], outdir)
            with open(os.path.join(outdir, 'Taxonomy.json'), 'w') as H:
                json.dump(out_file_list_tmp['json'], H, indent=2)
            s = f'Taxonomy\tD\t'
        except Exception as e:
            logging.error(f'Taxonomy {e}')
            s = f'Taxonomy\tE\t'

        try:
            write_status(status_report, s)
            post_url(taskID, 'Taxonomy')
        except Exception as e:
            logging.error(f'Taxonomy status {e}')

        try:
            ## FactorAnno
            s = f'FactorAnno\tR\t'
            write_status(status_report, s)
            out_file_list_tmp = FactorsAnno(scaffolds_fasta_file, tmp_dir, con_db_path['FactorsDB'], args.threads)
            faa_file = out_file_list_tmp['file'][0]
            out_file_list_tmp['file'].remove(faa_file)
            copy_file(out_file_list_tmp['file'], outdir)
            with open(os.path.join(outdir, 'FactorAnno.json'), 'w') as H:
                json.dump(out_file_list_tmp['json'], H, indent=2)
            s = f'FactorAnno\tD\t'
        except Exception as e:
            logging.error(f'FactorAnno {e}')
            s = f'FactorAnno\tE\t'

        try:
            write_status(status_report, s)
            post_url(taskID, 'FactorAnno')
        except Exception as e:
            logging.error(f'FactorAnno status {e}')

        try:
            s = f'Typing\tR\t'
            write_status(status_report, s)
        except Excception as e:
            logging.error(f'Typing status {e}')

        if species:
            try:
                ## Typing
                out_file_list_tmp = Fasta_Typing(scaffolds_fasta_file, tmp_dir, args.threads, args.cgmlst_db_path, species)
                copy_file(out_file_list_tmp['file'], outdir)
                with open(os.path.join(outdir, 'Typing.json'), 'w') as H:
                    json.dump(out_file_list_tmp['json'], H, indent=2)
                s = f'Typing\tD\t'
            except Exception as e:
                logging.error(f'Typing {e}')
                s = f'Typing\tE\t'
        else:
            s = f'Typing\tE\t'
        try:
            write_status(status_report, s)
            post_url(taskID, 'Typing')
        except Exception as e:
            logging.error(f'Typing status {e}')

        try:
            ## GeneAnno
            s = f'GeneAnno\tR\t'
            write_status(status_report, s)
            out_file_list_tmp = FunctionAnno(scaffolds_fasta_file, tmp_dir, con_db_path['GeneAnnoDB'], con_db_path['GeneAnnoInfo'], args.threads, db_anno_list, faafile=faa_file)
            copy_file(out_file_list_tmp['file'], outdir)
            with open(os.path.join(outdir, 'GeneAnno.json'), 'w') as H:
                json.dump(out_file_list_tmp['json'], H, indent=2)
            s = f'GeneAnno\tD\t'
        except Exception as e:
            logging.error(f'GeneAnno {e}')
            s = f'GeneAnno\tE\t'
        try:
            write_status(status_report, s)
            post_url(taskID, 'GeneAnno')
        except Exception as e:
            logging.error(f'GeneAnno status {e}')

    try:
        post_url(taskID, '2', 'http://localhost/task/getTaskRunningStatus/')
    except Exception as e:
        logging.error(f'post_url getTaskRunningStatus {e}')

    if not args.debug:
        try:
            shutil.rmtree(tmp_dir)
        except Exception as e:
            logging.error(e)
