import os
import re
import multiprocessing
import logging
import argparse
import shutil
import pandas as pd
import numpy as np

from mkdir import mkdir
import run_docker

def run_r_tree(cg_mlst_file, outdir):
    logging.info('r cgmlst tree')
    outfile = os.path.join(outdir, 'cgMLST.tree')
    r_script = os.path.join(bin_dir, 'R_cgMLST_tree.r')
    logfile = os.path.join(outdir, 'r.log')
    r_out = run_docker.r_docker('Rscript', r_script, cg_mlst_file, outfile, f'>{logfile} 2>&1')
    logging.info(f'{r_out}')
    if r_out:
        return 0
    return outfile

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

    df_cg_tmp = df_cg
    if df_cg.shape[0] > 2000:
        df_cg_tmp = df_cg.head(2000)
    df_total = df_cg_tmp.append(df_alle)[df_cg.columns]

    cg_total_file = os.path.join(outdir, 'cgMLST_total.csv')
    df_total.to_csv(cg_total_file, sep='\t')
    tree_file = run_r_tree(cg_total_file, outdir)
    out_file_list.append(tree_file)

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
    db_dict = pd.read_csv(schema_path, index_col=0).to_dict(orient='index')
    try:
        if species in db_dict:
            logging.info(db_dict[species])
            (cgmlst_csv, cgmlst_summary, cgmlst_tree) = run_chewBBACA(fa, outdir, db_dict[species]['schema_seed'], db_dict[species]['training_file'], db_dict[species]['cgFile'], threads)
            out_file_list['file'].extend([cgmlst_csv, cgmlst_summary, cgmlst_tree])
            out_file_list['json'].update({'cgmlst_csv': os.path.split(cgmlst_csv)[1]})
            out_file_list['json'].update({'cgmlst_tree': os.path.split(cgmlst_tree)[1]})
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
                    if not (tmp[0] == ''):
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

cpu_num = multiprocessing.cpu_count()
bin_dir = os.path.split(os.path.realpath(__file__))[0]
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='assembly scaffolds fasta file')
    parser.add_argument('-o', '--outdir', default='./', help='output dir')
    parser.add_argument('-s', '--species', required=True, help='species in cgMLST list')
    parser.add_argument('-dp', '--db_path', default=f'{os.path.join(bin_dir, "cgMLST_db.csv")}',
                        help='cgMLST database path for schema_seed, training file and core gene')
    parser.add_argument('-t', '--threads', default=cpu_num, help='threads')
    args = parser.parse_args()

    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)
    outdir = args.outdir
    mkdir(outdir)
    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    tmp_dir = os.path.join(outdir, 'tmp_dir')
    mkdir(tmp_dir)

    out_file_list = Fasta_Typing(args.input, tmp_dir, args.threads, args.db_path, args.species)

    for i in out_file_list:
        if os.path.isfile(i) or os.path.isdir(i):
            try:
                des_file = shutil.copy(i, args.outdir)
                logging.info(des_file)
            except Exception as e:
                logging.error(e)
        else:
            logging.error(f' not exitst {i}')