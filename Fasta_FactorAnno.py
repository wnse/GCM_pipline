import os
import multiprocessing
import argparse
import logging
import json
import re
import pandas as pd
from Bio import SeqIO

import run_docker
from mkdir import mkdir
from Fastq_Assemble import run_checkm
from plot_gc_cover import plot_len_dis
from post_status import post_url
from post_status import copy_file
from post_status import write_status


def get_gff_info(gfffile, outfile, outpng):
    total_info = []
    with open(gfffile, 'rt') as h:
        for i in h.readlines():
            if not re.match('^#', i):
                total_info.append(i.strip().split('\t'))
    df = pd.DataFrame(total_info, columns=['id', 'version', 'type', 'start', 'end', 'score', 'dir', 'none', 'info'])
    df['start'] = df['start'].astype(float)
    df['end'] = df['end'].astype(float)
    df['len'] = df['end'] - df['start'] + 1
    try:
        plot_len_dis(df['len'], outpng, xlabel='Gene Length')
    except Exception as e:
        logging.error(f'plot len dist {e}')
    total_gene = df.shape[0]
    avg_inter = ((df.groupby('id')['end'].max() - df.groupby('id')['start'].min() + 1 ) - df.groupby('id')['len'].sum()).sum() / total_gene
    partial_gene = df['info'].str.split(';').str[1].str.split('=').str[1].astype(float).astype(bool).sum() / total_gene
    avg_len = df['len'].mean()
    total_len = df['len'].sum()
    df_tmp = {'total_gene_num':total_gene, 'total_gene_len':total_len, 'avg_gene_len':avg_len, 'avg_inter_len':avg_inter, 'partial_gene_ratio':partial_gene}
    pd.DataFrame.from_dict(df_tmp, orient='index').to_csv(outfile, header=None)

def run_prodigal(fa, outdir):
    logging.info('prodigal')
    gene_info_file = os.path.join(outdir, 'prodigal_gene_info.csv')
    out_len_dis_png = os.path.join(outdir, 'prodigal_gene_len.png')
    cmd_out = run_docker.prodigal_docker(fa, outdir)
    if cmd_out[0]:
        logging.error(cmd_out)
        return 0, 0
    prodigal_out, gene_dir, faafile, gfffile = cmd_out
    try:
        get_gff_info(gfffile, gene_info_file, out_len_dis_png)
    except Exception as e:
        logging.error(e)
    return gene_dir, faafile, gene_info_file, out_len_dis_png

def run_diamond(faa, db, db_path, outdir, threads, id=40, query_cover=40, subject_cover=40):
    logging.info(f'{db}')
    logging.info(f'{db_path}')
    outfile = os.path.join(outdir, db + '.diamond.out')
    diamond_out = run_docker.diamond_docker(faa, outfile, db_path, threads=threads, id=id, query_cover=query_cover, subject_cover=subject_cover)
    logging.info(f'{diamond_out}')
    if diamond_out:
        return diamond_out
    return outfile


def run_rgi(fa, outdir, threads):
    logging.info('rgi')
    cmd_out = run_docker.rgi_docker(fa, outdir, threads)
    if cmd_out[0]:
        logging.error(cmd_out)
        return 0
    RGI_out, rgi_out_file = cmd_out
    logging.info(f'{RGI_out}')
    if os.path.isfile(rgi_out_file):
        df = pd.read_csv(rgi_out_file, sep='\t')
        out_list = ['Contig', 'Best_Hit_ARO', 'Drug Class', 'AMR Gene Family', 'Best_Identities',
                    'Percentage Length of Reference Sequence']
        outfile = os.path.join(outdir, 'FactorAnno_rgi.csv')
        df[out_list].to_csv(outfile, index=False)
        return outfile
    return rgi_out_file

def process_vfdb(diamond_out, vfdb_fa, outfile):
    try:
        logging.info('process_vfdb')
        df_diamond = pd.read_csv(diamond_out, sep='\t', header=None)
        pattern = re.compile('(VFG\\d+).*\\((.*)\\) (.*) \\[(.+) \\((.*)\\)\\] \\[(.+)\\]')
        vfdb_info = {rec.id:(pattern.match(rec.description).group(2, 4, 3)) for rec in SeqIO.parse(vfdb_fa, "fasta")}
        vfdb_seqLen = {rec.id: len(rec.seq) for rec in SeqIO.parse(vfdb_fa, "fasta")}
        df_diamond['coverage'] = 100 * (df_diamond[9] - df_diamond[8]) / df_diamond[1].map(vfdb_seqLen)
        df_diamond['gene'] = df_diamond[1].map(vfdb_info).str[0]
        df_diamond['VFs'] = df_diamond[1].map(vfdb_info).str[1]
        df_diamond['gene_product'] = df_diamond[1].map(vfdb_info).str[2]
        out_list = [0, 1, 'gene', 'VFs', 'gene_product', 2, 'coverage']
        df_diamond[out_list].rename(columns={0:'Contig', 1:'Best_Hit', 2:'Identities'}).to_csv(outfile, index=False)
    except Exception as e:
        logging.error(f'process_vfdb {e}')

def FactorsAnno(fa, outdir, db_path, threads, checkm=False, faafile=None):
    mkdir(outdir)
    out_file_list = {}
    out_file_list['json'] = {}
    out_file_list['file'] = []
    if checkm:
        _, checkm_outfile_csv = run_checkm(fa, outdir, threads)
        if checkm_outfile_csv:
            out_file_list['file'].append(checkm_outfile_csv)
            df_tmp = pd.read_csv(checkm_out_file).astype(str)
            df_tmp.index = df_tmp.index.astype(str).str.lower().str.replace(' ', '_').str.replace('.', '_', regex=False)
            tmp_dict = {'check_evaluate': df_tmp.loc[0].to_dict()}
            out_file_list['json'].update(tmp_dict)
    if not faafile:
        prodigal_dir, faafile, gene_info_file, out_len_dis_png= run_prodigal(fa, outdir)
        if prodigal_dir:
            out_file_list['file'].append(faafile)
            out_file_list['file'].append(prodigal_dir)
            out_file_list['file'].append(gene_info_file)
            out_file_list['file'].append(out_len_dis_png)
            df_tmp = pd.read_csv(gene_info_file, header=None, index_col=0)[1].astype(str).to_dict()
            out_file_list['json'].update({'gene_info_summary':df_tmp})
            out_file_list['json'].update({'gene_predicted_dir':os.path.split(prodigal_dir)[1]})
            out_file_list['json'].update({'gene_len_dis_png': os.path.split(out_len_dis_png)[1]})
    if faafile:
        if "VFDB" in db_path:
            try:
                vfdb_path, vfdb_fa = db_path["VFDB"].split()
                vfdb_diamond_out = run_diamond(faafile, "VFDB", vfdb_path, outdir, threads, id=40, query_cover=40, subject_cover=40)
                out_file_tmp = os.path.join(outdir, 'FactorAnno_VFs.csv')
                process_vfdb(vfdb_diamond_out, vfdb_fa, out_file_tmp)
                out_file_list['file'].append(out_file_tmp)
                out_file_list['json'].update({'vfdb_csv':os.path.split(out_file_tmp)[1]})
            except Exception as e:
                logging.error(f'FactorsAnno VFDB {e}')
        try:
            RGI_out = run_rgi(fa, outdir, threads)
            if RGI_out:
                out_file_list['file'].append(RGI_out)
                out_file_list['json'].update({'rgi_csv': os.path.split(RGI_out)[1]})
        except Exception as e:
            logging.error(f'FactorsAnno RGI {e}')
    return out_file_list



cpu_num = multiprocessing.cpu_count()
bin_dir = os.path.split(os.path.realpath(__file__))[0]
if __name__ == '__main__':
    cpu_num = multiprocessing.cpu_count()
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='genome fasta file')
    parser.add_argument('-o', '--outdir', default='./', help='output dir')
    parser.add_argument('-n', '--name', default='test', help='sample name')
    parser.add_argument('-t', '--threads', default=cpu_num, help='threads')
    parser.add_argument('-tID', '--taskID', default='', help='task ID for report status')
    parser.add_argument('-c', '--checkm', action='store_true', help='evaluate genome by checkm')
    parser.add_argument('-dp', '--db_path', default=f'{os.path.join(bin_dir, "conf_DB_path.json")}', help='conf_DB_path.json genome database path info')
    args = parser.parse_args()

    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)
    outdir = args.outdir
    mkdir(outdir)
    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    taskID = args.taskID

    if os.path.isfile(args.db_path):
        with open(args.db_path, 'rt') as h:
            db_path = json.load(h)
        logging.info(db_path['FactorsDB'])
        tmp_dir = os.path.join(outdir, 'tmp_dir')
        mkdir(tmp_dir)

        try:
            s = f'FactorAnno\tR\t'
            write_status(status_report, s)
            out_file_list_tmp = FactorsAnno(args.input, tmp_dir, db_path['FactorsDB'], args.threads, args.checkm)
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

        # for i in out_file_list:
        #     if os.path.isfile(i) or os.path.isdir(i):
        #         try:
        #             des_file = shutil.copy(i, args.outdir)
        #             logging.info(des_file)
        #         except Exception as e:
        #             logging.error(e)
        #     else:
        #         logging.error(f' not exitst {i}')
    else:
        logging.error(f'not file {args.db_path}')