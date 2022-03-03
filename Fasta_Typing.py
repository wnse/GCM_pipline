import os
import re
import multiprocessing
import logging
import argparse
import shutil
import pandas as pd
import numpy as np
import json

from mkdir import mkdir
import run_docker
from Fasta_Taxonomy import Fasta_Taxonomy
from post_status import post_url
from post_status import post_pid
from post_status import copy_file
from post_status import write_status


def run_chewBBACA(fa, outdir, schema_dir, training_file, cgmlstfile, threads=1):
    out_file_list = []
    mkdir(outdir)
    scaffolds_fasta_dir = os.path.join(outdir, 'scaffolds_fasta')
    mkdir(scaffolds_fasta_dir)
    try:
        shutil.copy(fa, scaffolds_fasta_dir)
    except Exception as e:
        logging.error(e)
    logging.info('chewBBACA')
    cmd_out = run_docker.chewBBACA_docker(scaffolds_fasta_dir, outdir, schema_dir, training_file=training_file, threads=threads)
    if cmd_out[0]:
        logging.error(cmd_out)
        return cmd_out
    chewbbaca_out, alleles_file = cmd_out
    logging.info(f'{chewbbaca_out}')
    out_file_list.append(alleles_file)
    df_cg = pd.read_csv(cgmlstfile, index_col=0, sep='\t', low_memory=False)
    df_a = pd.read_csv(alleles_file, index_col=0, sep='\t')
    df_alle = df_a.loc[df_a.index[0]]
    mask = df_alle.str.isnumeric().isna()
    find_alle_in_cg = len(np.intersect1d(df_alle[mask].index, df_cg.columns))
    summaryfile = os.path.join(outdir, 'cgMLST_summary.csv')
    out_file_list.append(summaryfile)
    with open(summaryfile, 'w') as h:
        print(f'cg_num,cg_miss_num,cg_miss_ratio', file=h)
        print(f'{df_cg.shape[1]},{df_cg.shape[1] - find_alle_in_cg},{round(1 - find_alle_in_cg / df_cg.shape[1] , 2)}', file=h)

    # df_cg_tmp = df_cg
    # if df_cg.shape[0] > 2000:
    #     df_cg_tmp = df_cg.head(2000)
    # df_total = df_cg_tmp.append(df_alle)[df_cg.columns]
    #
    # cg_total_file = os.path.join(outdir, 'cgMLST_total.csv')
    # df_total.to_csv(cg_total_file, sep='\t')
    # tree_file = run_r_tree(cg_total_file, outdir)
    # out_file_list.append(tree_file)

    return out_file_list

def run_SeqSero2(fa, outdir, threads=1):
    mkdir(outdir)
    logging.info('SeqSero2')
    cmd_out = run_docker.SeqSero2_docker(fa, outdir, threads)
    if cmd_out[0]:
        logging.error(cmd_out)
        return cmd_out, 0
    SeqSero2_out, outfile = cmd_out
    logging.info(f'{SeqSero2_out}')
    return SeqSero2_out, outfile

def run_mlst(fa, outdir, threads=1):
    mkdir(outdir)
    logging.info('mlst')
    cmd_out = run_docker.mlst_docker(fa, outdir, threads)
    if cmd_out[0]:
        logging.error(cmd_out)
        return cmd_out, 0
    mlst_out, outfile = cmd_out
    logging.info(f'{mlst_out}')
    return mlst_out, outfile

def Fasta_Typing(fa, outdir, threads, schema_path, species):
    mkdir(outdir)
    out_file_list = {}
    out_file_list['json'] = {}
    out_file_list['file'] = []
    db_dict = pd.read_csv(schema_path, index_col=0).fillna('').to_dict(orient='index')
    logging.info(species)
    try:
        species_used = None
        if species in db_dict:
            species_used = species
        else:
            if ' '.join(species.split()[-2:]) in db_dict:
                species_used = ' '.join(species.split()[-2:])
            if species.split()[1] in db_dict:
                species_used = species.split()[1]
            elif species.split()[0] in db_dict:
                species_used = species.split()[0]

        if species_used in db_dict:
            logging.info(db_dict[species_used])
            # (cgmlst_csv, cgmlst_summary, cgmlst_tree) = run_chewBBACA(fa, outdir, db_dict[species]['schema_seed'], db_dict[species]['training_file'], db_dict[species]['cgFile'], threads)
            (cgmlst_csv, cgmlst_summary) = run_chewBBACA(fa, outdir, db_dict[species_used]['schema_seed'], db_dict[species_used]['training_file'], db_dict[species_used]['cgFile'], threads)

            # out_file_list['file'].extend([cgmlst_csv, cgmlst_summary, cgmlst_tree])
            # out_file_list['json'].update({'cgmlst_tree': os.path.split(cgmlst_tree)[1]})
            out_file_list['file'].extend([cgmlst_csv, cgmlst_summary])
            out_file_list['json'].update({'cgmlst_csv': os.path.split(cgmlst_csv)[1]})
            df_tmp = pd.read_csv(cgmlst_summary).astype(str).loc[0].to_dict()
            out_file_list['json'].update({'cgmlst_summary':df_tmp})
        else:
            logging.info(f'{species} not in cgMLST database')
            # cgmlst_summary = os.path.join(outdir, 'cgMLST_summary.csv')
            # with open(cgmlst_summary, 'a') as h:
            #     print(f'{species} not in cgMLST Database', file=h)
            # out_file_list.append(cgmlst_summary)
    except Exception as e:
        logging.error(f'Fasta_Typing cgMLST {e}')
    try:
        serotype, serotype_result = run_SeqSero2(fa, outdir, threads)
        out_file_list['file'].append(serotype_result)
        tmp_dict = {}
        with open(serotype_result, 'rt') as h:
            for i in h.readlines():
                try:
                    tmp = i.strip().split(r':')
                    if len(tmp)>1 and (not (tmp[0] == '')):
                        tmp[1] = re.sub('\t','',tmp[1])
                        key = re.sub(r'\(\S+\)', '', tmp[0]).lower().replace(' ', '_')
                        tmp_dict[key] = ':'.join(tmp[1:])
                except Exception as e:
                    logging.error(f'Fasta_Typing serotype {e}')
        out_file_list['json'].update({'serotype_summary':tmp_dict})
    except Exception as e:
        logging.error(f'Fasta_Typing serotype {e}')
    try:
        mlst_out, mlst_out_file = run_mlst(fa, outdir, threads)
        out_file_list['file'].append(mlst_out_file)
        df_tmp = {}
        with open(mlst_out_file, 'rt') as h:
            tmp = h.readlines()[0].split()
            df_tmp = {'tax':tmp[1], 'ST':tmp[2], 'st_info':tmp[3:]}
        # df_tmp = pd.read_csv(mlst_out_file, header=None, sep='\t').astype(str)
        # df_tmp = df_tmp.rename(columns={0: 'name', 1: 'tax', 2: 'ST', 3: 'dnaE', 4: 'gyrB', 5: 'recA', 6: 'dtdS', 7: 'pntA', 8: 'pyrC', 9: 'tnaA'})
        # df_tmp = df_tmp.loc[0].to_dict()
        out_file_list['json'].update({'mlst_summary':df_tmp})
    except Exception as e:
        logging.error(f'Fasta_Typing mlst {e}')
    return out_file_list


if __name__ == '__main__':
    cpu_num = multiprocessing.cpu_count()
    if cpu_num > 4:
        cpu_num = cpu_num - 2
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='assembly scaffolds fasta file')
    parser.add_argument('-o', '--outdir', default='./', help='output dir')
    parser.add_argument('-s', '--species', default=None, help='species in cgMLST list')
    parser.add_argument('-dp', '--cgmlst_db_path', default=f'{os.path.join(bin_dir, "cgMLST_db.csv")}',
                        help='cgMLST database path for schema_seed, training file and core gene')
    parser.add_argument('-t', '--threads', default=cpu_num, help='threads')
    parser.add_argument('-tID', '--taskID', default='', help='task ID for report status')
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
    parser.add_argument('-debug', '--debug', action='store_true')
    parser.add_argument('-n', '--name', default='test', help='sample name')
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
    post_pid(taskID)
    scaffolds_fasta_file = args.input
    species = args.species
    if not species:
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
        except Exception as e:
            logging.error(f'Taxonomy status {e}')

    if species:
        try:
            ## Typing
            s = f'Typing\tR\t'
            write_status(status_report, s)
            out_file_list_tmp = Fasta_Typing(scaffolds_fasta_file, tmp_dir, args.threads, args.cgmlst_db_path, species)
            copy_file(out_file_list_tmp['file'], outdir)
            with open(os.path.join(outdir, 'Typing.json'), 'w') as H:
                json.dump(out_file_list_tmp['json'], H, indent=2)
            s = f'Typing\tD\t'
        except Exception as e:
            logging.error(f'Typing {e}')
            s = f'Typing\tE\t'
        try:
            write_status(status_report, s)
            try:
                post_url(taskID, '2', 'http://localhost/task/getTaskRunningStatus/')
            except Exception as e:
                logging.error(f'post_url getTaskRunningStatus {e}')
            # post_url(taskID, 'Typing')
        except Exception as e:
            logging.error(f'Typing status {e}')

    if not args.debug:
        try:
            shutil.rmtree(tmp_dir)
        except Exception as e:
            logging.error(e)

    # out_file_list = Fasta_Typing(args.input, tmp_dir, args.threads, args.cgmlst_db_path, args.species)
    #
    # for i in out_file_list:
    #     if os.path.isfile(i) or os.path.isdir(i):
    #         try:
    #             des_file = shutil.copy(i, args.outdir)
    #             logging.info(des_file)
    #         except Exception as e:
    #             logging.error(e)
    #     else:
    #         logging.error(f' not exitst {i}')

