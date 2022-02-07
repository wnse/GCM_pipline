import os
import shutil
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

def run_rnammer(fa, outdir):
    logging.info('rnammer')
    rnammer_out = run_docker.rnammer_docker(fa, outdir, m='lsu,tsu,ssu')
    logging.info(rnammer_out)
    return rnammer_out

def run_blastn(query, outfile, threads):
    logging.info('blastn')
    db = '/BioBin/ref/bacteria.16s'
    img_name = 'registry.servicemgr.bjwsws:5000/blast:16s'
    tmp_out = outfile + '.blast.out'
    blastn_out = run_docker.blastn_docker(query, db, tmp_out, threads=threads, max_target_seqs=1, img_name=img_name)
    logging.info(blastn_out)
    return blastn_out

def run_trf(fa, outdir):
    logging.info('trf')
    trf_out = run_docker.trf_docker(fa, outdir)
    logging.info(trf_out)
    return trf_out

def run_pilercr(fa, outdir):
    logging.info('pilercr')
    outfile = os.path.join(outdir, 'pilercr.out')
    pilercr_out = run_docker.pilercr_docker(fa, outfile)
    logging.info(pilercr_out)
    return pilercr_out

def run_trnascan(fa, outdir):
    logging.info('tRNAscan se')
    trnascan_out = run_docker.trnascan_docker(fa, outdir)
    logging.info(trnascan_out)
    return trnascan_out

def run_Rfam(fa, outdir, db_path):
    logging.info('Rfam')
    rfam_out = run_docker.Rfam_docker(fa, outdir, db_path)
    logging.info(rfam_out)
    return rfam_out

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
        plot_len_dis(df['len'], outpng)
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

def run_pfam(faa, outdir, db_path, threads):
    logging.info('pfam')
    outfile = os.path.join(outdir, 'pfam.out')
    pfam_out = run_docker.pfam_docker(faa, outfile, db_path, threads=threads)
    logging.info(f'{pfam_out}')
    return pfam_out

def run_diamond(faa, db, db_path, outdir, threads, id=40, query_cover=40, subject_cover=40):
    logging.info(f'{db} {db_path}')
    outfile = os.path.join(outdir, db + '.diamond.out')
    diamond_out = run_docker.diamond_docker(faa, outfile, db_path, threads=threads, id=id, query_cover=query_cover, subject_cover=subject_cover)
    logging.info(f'{diamond_out}')
    if diamond_out:
        return diamond_out
    return outfile

def run_antismash(fa, outdir, threads):
    logging.info('antismash')
    cmd_out = run_docker.antismash_docker(fa, outdir, threads)
    if cmd_out[0]:
        logging.error(cmd_out)
        return cmd_out
    antismash_out, antismash_dir = cmd_out
    logging.info(antismash_out)
    return antismash_dir


def process_kegg(diamond_out, db_info, outfile, col_list=[0,1]):
    try:
        df_diamond = pd.read_csv(diamond_out, sep='\t', header=None, index_col=0)
        if df_diamond.empty:
            return 0
        subject_list = set(df_diamond[1].to_list())
        subject_info = []
        with open(db_info, 'rt') as h:
            for i in h:
                if not re.match('^#', i):
                    info = i.strip().split('\t')
                    if info[0] in subject_list:
                        subject_info.append([info[col] for col in col_list if len(info) > col])
        df_subject = pd.DataFrame(subject_info).set_index(0)
        df_merge = pd.merge(df_diamond, df_subject, left_on=1, right_index=True, how='left')
        df_merge.to_csv(outfile, header=False, sep='\t')
        return len(df_diamond.index.unique())
    except Exception as e:
        logging.error(e)

def process_phi(diamond_out, db_info, outfile):
    try:
        df_diamond = pd.read_csv(diamond_out, sep='\t', header=None, index_col=0)
        df_subject = pd.read_csv(db_info, header=None, sep='#', index_col=1)
        df_merge = pd.merge(df_diamond, df_subject, left_on=1, right_index=True, how='left')
        df_merge.to_csv(outfile, header=False, sep='\t')
        return df_merge.shape[0]
    except Exception as e:
        logging.error(e)

def process_cazy(diamond_out, db_info, outfile):
    try:
        df_diamond = pd.read_csv(diamond_out, sep='\t', header=None, index_col=0)
        # df_subject = pd.read_csv(db_info, header=None, sep='#', index_col=1)
        # df_merge = pd.merge(df_diamond, df_subject, left_on=1, right_index=True, how='left')
        df_diamond['ref'] = df_diamond[1].str.split('|').str[-1]
        df_diamond.to_csv(outfile, header=False, sep='\t')
        return len(df_diamond.index.unique())
    except Exception as e:
        logging.error(e)

def process_cog(diamond_out, db_info, outfile):
    try:
        gi_file, class_file = db_info.split()
        df_diamond = pd.read_csv(diamond_out, sep='\t', header=None)
        df_diamond['gi'] = df_diamond[1].str.split('|').str[1]
        df_gi_cog = pd.read_csv(gi_file, header=None, dtype='str')
        df_class = pd.read_csv(class_file, sep='\t', header=None)
        df_merge = pd.merge(df_diamond, df_gi_cog[[0,6]].rename(columns={0:'gi', 6:'cog'}), left_on='gi', right_on='gi', how='left')
        df_merge = df_merge.drop_duplicates(subset=0, keep='first')
        df_merge = pd.merge(df_merge, df_class.rename(columns={0:'cog'}), left_on='cog', right_on='cog', how='left')
        df_merge.drop('gi', axis=1).to_csv(outfile, header=None, sep='\t', index=False)
        return len(df_diamond[0].unique())
    except Exception as e:
        logging.error(e)
        return 0

def StructureAnno(fa, outdir, db_path, threads, checkm=False):
    out_file_list = []
    if checkm:
        _, checkm_outfile_csv = run_checkm(fa, outdir, threads)
        if checkm_outfile_csv:
            out_file_list.append(checkm_outfile_csv)
    rnammer_out = run_rnammer(fa, outdir)
    if rnammer_out:
        logging.info(f'{rnammer_out}')
    # trf_out = run_trf(fa, outdir)
    pilercr_out = run_pilercr(fa, outdir)
    # trnascan_out = run_trnascan(fa, outdir)
    # rfam_out = run_Rfam(fa, outdir, db_path)

def FunctionAnno(fa, outdir, db_path, db_info, threads, db_list, faafile=None):
    mkdir(outdir)
    out_file_list = {}
    out_file_list['json'] = {}
    out_file_list['file'] = []
    if not faafile:
        prodigal_dir, faafile, gene_info_file, out_len_dis_png = run_prodigal(fa, outdir)
        if prodigal_dir:
            out_file_list['file'].append(faafile)
            out_file_list['file'].append(prodigal_dir)
            out_file_list['file'].append(gene_info_file)
            out_file_list['file'].append(out_len_dis_png)
            df_tmp = pd.read_csv(gene_info_file, header=None, index_col=0)[1].astype(str).to_dict()
            out_file_list['json'].update({'gene_info_summary':df_tmp})
            out_file_list['json'].update({'gene_len_dis_png': os.path.split(out_len_dis_png)[1]})
            out_file_list['json'].update({'gene_predicted_dir':os.path.split(prodigal_dir)[1]})

    # antisamsh_dir = run_antismash(fa, outdir, threads)
    # out_file_list.append(antisamsh_dir)
    if faafile:
        gene_count = len(list(SeqIO.parse(faafile, 'fasta')))
        anno_count = {}
        for db in db_list:
            if db in db_path:
                if db == 'pfam':
                    # next
                    pfam_out = run_pfam(faafile, outdir, db_path['pfam'], threads)
                else:
                    diamond_out = run_diamond(faafile, db, db_path[db], outdir, threads)
                    outfile = ''
                    if os.path.isfile(diamond_out):
                        if db in ['kegg', 'MetaCyc'] and db in db_info:
                            outfile = os.path.join(outdir, str(db)+'_anno_table.csv')
                            anno_count[db] = process_kegg(diamond_out, db_info[db], outfile)
                            out_file_list['file'].append(outfile)
                        if db in ['NR'] and db in db_info:
                            outfile = os.path.join(outdir, str(db)+'_anno_table.csv')
                            anno_count[db] = process_kegg(diamond_out, db_info[db], outfile, col_list=[0, 1, 2, 3])
                            out_file_list['file'].append(outfile)
                        if db in ['SwissProt'] and db in db_info:
                            outfile = os.path.join(outdir, str(db)+'_anno_table.csv')
                            anno_count[db] = process_kegg(diamond_out, db_info[db], outfile, col_list=[0, 1, 2, 5])
                            out_file_list['file'].append(outfile)
                        if db in ['PHI'] and db in db_info:
                            outfile = os.path.join(outdir, str(db)+'_anno_table.csv')
                            anno_count[db] = process_phi(diamond_out, db_info[db], outfile)
                            out_file_list['file'].append(outfile)
                        if db in ['cazy'] and db in db_info:
                            outfile = os.path.join(outdir, str(db)+'_anno_table.csv')
                            anno_count[db] = process_cazy(diamond_out, db_info[db], outfile)
                            out_file_list['file'].append(outfile)
                        if db == 'cog' and db in db_info:
                            outfile = os.path.join(outdir, 'cog_anno_table.csv')
                            anno_count[db] = process_cog(diamond_out, db_info[db], outfile)
                            out_file_list['file'].append(outfile)
                        if outfile:
                            out_file_list['json'].update({f'{db}_gene_anno_csv': os.path.split(outfile)[1]})
            else:
                logging.info(f'No database path for {db}')

        anno_summary = os.path.join(outdir, 'Gene_Anno_summary.txt')
        df = pd.DataFrame.from_dict(anno_count, orient='index').rename_axis(index='db')
        df.columns = ['gene_count']
        df['ratio'] = df['gene_count'] / gene_count
        df.to_csv(anno_summary)
        df_tmp = [v for i, v in df.astype(str).reset_index().to_dict(orient='index').items()]
        out_file_list['file'].append(anno_summary)
        out_file_list['json'].update({'gene_anno_summary':df_tmp})
    return out_file_list


if __name__ == '__main__':
    cpu_num = multiprocessing.cpu_count()
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='genome fasta file')
    parser.add_argument('-o', '--outdir', default='./', help='output dir')
    parser.add_argument('-n', '--name', default='test', help='sample name')
    parser.add_argument('-t', '--threads', default=cpu_num, help='threads')
    parser.add_argument('-c', '--checkm', action='store_true', help='evaluate genome by checkm')
    parser.add_argument('-db', '--db', default=['kegg','cog', 'NR', 'MetaCyc', 'cazy', 'PHI', 'SwissProt'], help='database for annotation', nargs='+')
    parser.add_argument('-dp', '--db_path', default=f'{os.path.join(bin_dir, "conf_DB_path.json")}', help='conf_DB_path.json genome database path info')
    args = parser.parse_args()

    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)
    outdir = args.outdir
    mkdir(outdir)
    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    if os.path.isfile(args.db_path):
        with open(args.db_path, 'rt') as h:
            db_path = json.load(h)
        db_anno_list = args.db
        logging.info(db_anno_list)
        logging.info(db_path)
        tmp_dir = os.path.join(outdir, 'tmp_dir')
        mkdir(tmp_dir)
        status_report = os.path.join(outdir, 'status_report.txt')

        try:
            s = f'GeneAnno\tR\t'
            write_status(status_report, s)
            out_file_list_tmp = FunctionAnno(args.input, tmp_dir, db_path['GeneAnnoDB'], db_path['GeneAnnoInfo'],
                                         args.threads, db_anno_list)
            copy_file(out_file_list_tmp['file'], outdir)

            with open(os.path.join(outdir, 'GeneAnno.json'), 'w') as H:
                json.dump(out_file_list_tmp['json'], H, indent=2)
            s = f'GeneAnno\tD\t'
        except Exception as e:
            logging.error(f'GeneAnno {e}')
            s = f'GeneAnno\tE\t'

        try:
            write_status(status_report, s)
            # post_url(taskID, 'GeneAnno')
        except Exception as e:
            logging.error(f'GeneAnno status {e}')

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
        logging.info(f'{args.db_path} not exists')

        # stru_anno = StructureAnno(args.fasta, args.outdir, args.db_path, args.threads, checkm=args.checkm)
        # print(stru_anno)



