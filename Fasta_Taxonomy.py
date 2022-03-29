import os
import multiprocessing
import shutil
import argparse
import logging
import re
import json
import pandas as pd
from Bio import SeqIO

from mkdir import mkdir
import run_docker
from post_status import post_url
from post_status import post_pid
from post_status import copy_file
from post_status import write_status

def filter_fa_by_len(fa, outfa, len_cutoff=1200):
    tmp_len_list = [len(seq.seq) for seq in SeqIO.parse(fa, 'fasta')]
    if tmp_len_list:
        max_len = max(tmp_len_list)
    else:
        return 0
    if max_len > len_cutoff:
        max_len = len_cutoff
    res_fa = [seq for seq in SeqIO.parse(fa, 'fasta') if len(seq.seq) >= max_len]
    res_count = SeqIO.write(res_fa, outfa, 'fasta')
    return res_count

def run_rnammer(fa, outdir):
    logging.info('rnammer')
    cmd_out = run_docker.rnammer_docker(fa, outdir)
    if cmd_out[0]:
        logging.error(cmd_out)
        return cmd_out, 0
    rnammer_out, fa_file = cmd_out
    logging.info(rnammer_out)
    if os.path.isfile(fa_file):
        filter_fa_file = os.path.join(outdir, 'RNAmmer.filter.fasta')
        filter_count = filter_fa_by_len(fa_file, filter_fa_file)
        logging.info(f'filtered 16S sequences {filter_count}')
        if filter_count:
            return rnammer_out, filter_fa_file
        return rnammer_out, 0
    return rnammer_out, 0

def run_blastn(query, db, outfile, threads):
    logging.info('blastn')
    blastn_out = run_docker.blastn_docker(query, db, outfile, threads=threads)
    logging.info(blastn_out)
    return blastn_out

def merge_blast_with_tax(blast_16S_out, db_16S_info, blast_16S_out_tax, cutoff=10):
    logging.info('merge_blast_with_tax')
    df_blast = pd.read_csv(blast_16S_out, sep='\t', header=None)
    df_db = pd.read_csv(db_16S_info, sep='\t', header=None)
    df_blast.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    df_merge = pd.merge(df_blast.head(cutoff), df_db.set_index(0)[[2,13,14]], left_on='sseqid', right_on=0, how='left')
    df_merge = df_merge.rename(columns={2:'slen',13:'species', 14:'lineage'})
    df_merge['q_pos'] = df_merge['qstart'].astype(str) + '-' + df_merge['qend'].astype(str)
    df_merge['s_pos'] = df_merge['sstart'].astype(str) + '-' + df_merge['send'].astype(str)
    df_merge['completness'] = df_merge['length']/df_merge['slen']
    df_merge.loc[df_merge['completness']>100, 'completness'] = 100
    out_list = ['sseqid', 'completness', 'length' , 'q_pos', 's_pos', 'pident', 'bitscore', 'species', 'lineage']
    df_merge[out_list].to_csv(blast_16S_out_tax, index=False)

def Fasta_16S_Taxonomy(fa, outdir, db_16S, db_16S_info, threads):
    mkdir(outdir)
    rnammer_out, fa_16S = run_rnammer(fa, outdir)
    if fa_16S:
        blast_16S_out = os.path.join(outdir, '16S.blast')
        blastn_out = run_blastn(fa_16S, db_16S, blast_16S_out, threads)
        blast_16S_out_tax = os.path.join(outdir, 'tax.16S.csv')
        if blastn_out == 0:
            merge_blast_with_tax(blast_16S_out, db_16S_info, blast_16S_out_tax, cutoff=10)
        return blast_16S_out_tax
    return 0

def run_mash_dist(query, db, outfile, evalue=0.05, threads=1):
    logging.info('mash dist')
    mash_dist_out = run_docker.mash_dist_docker(query, db, outfile, evalue=evalue, threads=threads)
    logging.info(mash_dist_out)
    return mash_dist_out

def merge_msh_with_tax(mash_out, db_info, mash_out_tax, cutoff=10):
    logging.info('merge_msh_with_tax')
    df_mash = pd.read_csv(mash_out, sep='\t', header=None)
    df_db = pd.read_csv(db_info, sep='\t')
    # df_db['Species'] = df_db['Tax'].str.split('|').str[-1].str.split('__').str[-1]
    # df_db.loc[df_db['Species'].isin(['', '-']), 'Species'] = df_db.loc[df_db['Species'].isin(['', '-']), 'Species ID']
    df_mash.columns = ['reference', 'query', 'mash-distance', 'p-value', 'matching-hashes']
    df_mash['reference'] = df_mash['reference'].str.strip('.fna')
    df_merge = pd.merge(df_mash.head(cutoff), df_db, left_on='reference', right_on='SampleID', how='left')
    out_list = ['reference', 'mash-distance', 'Species ID', 'Strain id', 'Tax']
    df_merge[out_list].to_csv(mash_out_tax, index=False)

def cal_ANI(fa1, fa2, outdir, threads=1):
    res_ani_total = {}
    logging.info(f'{os.path.split(fa1)[1]} vs {os.path.split(fa2)[1]}')
    logging.info('cal ANI: fastANI')
    outfile = os.path.join(outdir, 'fastANI.out')
    tmp_out = run_docker.fastANI_docker(fa1, fa2, outfile, threads=threads)
    logging.info(f'{tmp_out}')
    if tmp_out == 0:
        try:
            with open(outfile, 'rt') as h:
                ani_value = h.readline().strip().split('\t')[2]
                res_ani_total['fastANI'] = ani_value
        except Exception as e:
            logging.error(e)

    logging.info('cal ANI: OrthoANI')
    outfile = os.path.join(outdir, 'OrthoANI.out')
    tmp_out = run_docker.OrthoANI_docker(fa1, fa2, outfile, threads=threads)
    logging.info(f'{tmp_out}')
    if tmp_out == 0:
        try:
            with open(outfile, 'rt') as h:
                res_ani_total['OrthoANI'] = [i.split(' ')[2] for i in h.readlines() if re.match('OrthoANI :', i)][-1]
        except Exception as e:
            logging.error(e)

    logging.info('cal ANI: OAU')
    outfile = os.path.join(outdir, 'OAU.out')
    tmp_out = run_docker.OAU_docker(fa1, fa2, outfile, threads=threads)
    logging.info(f'{tmp_out}')
    if tmp_out == 0:
        try:
            with open(outfile, 'rt') as h:
                res_ani_total['OAU'] = h.readlines()[-1].split('\t')[1]
        except Exception as e:
            logging.error(e)

    # outfile = os.path.join(outdir, 'ANI.out')
    # tmp_out = run_docker.ANI_docker(fa1, fa2, outfile, threads=threads)
    # logging.info(f'{outfile} {tmp_out}')
    return res_ani_total

def Fasta_genome_Taxonomy(fa, outdir, db_genome, db_info, db_genome_fa_dir, threads):
    mkdir(outdir)
    mash_dist_out_file = os.path.join(outdir, 'genome.mash')
    mash_dist_out = run_mash_dist(fa, db_genome, mash_dist_out_file, threads=threads)
    if mash_dist_out == 0:
        mash_dist_out_tax = os.path.join(outdir, 'tax.genome.csv')
        try:
            merge_msh_with_tax(mash_dist_out_file, db_info, mash_dist_out_tax, cutoff=3)
        except Exception as e:
            logging.error(f'merge_msh_with_tax, {e}')
        df_genome_tax = pd.read_csv(mash_dist_out_tax, index_col=0)
        if os.path.isfile(mash_dist_out_tax):
            ref_list = pd.read_csv(mash_dist_out_tax)['reference']
            ref_ANI = {}
            cal_ani_dir = os.path.join(outdir, 'cal_ani_dir')
            mkdir(cal_ani_dir)
            for ref in ref_list:
                ref_path = os.path.join(db_genome_fa_dir, ref+'.fna')
                if os.path.isfile(ref_path):
                    fa_link = os.path.join(cal_ani_dir, os.path.split(fa)[1])
                    ref_link = os.path.join(cal_ani_dir, os.path.split(ref_path)[1])
                    shutil.copy(fa, fa_link)
                    shutil.copy(ref_path, ref_link)
                    ref_ANI[ref] = cal_ANI(fa_link, ref_link, outdir, threads=threads)
            genome_tax_ani = os.path.join(outdir, 'genome.ani')
            df_ani = pd.DataFrame.from_dict(ref_ANI, orient='index')
            df_ani.to_csv(genome_tax_ani)
            df_genome_tax_ani = pd.merge(df_genome_tax, df_ani, left_index=True, right_index=True, how='left')
            df_genome_tax_ani.to_csv(mash_dist_out_tax)
        return mash_dist_out_tax
    else:
        return 0

def merge_tax_genome_16s(tax_file_genome, tax_file_16s, outfile, name):
    logging.info('merge_tax_genome_16s')
    family_genome, genus_genome, species_genome, strain_genome, ani_genome, family_16s, genus_16s, species_16S = [''] * 8
    try:
        if tax_file_16s:
            df_16s = pd.read_csv(tax_file_16s)
            species_16S = df_16s.loc[0,'species']
            genus_16s = df_16s.loc[0,'lineage'].split('|')[-2].split('_')[-1]
            family_16s = df_16s.loc[0,'lineage'].split('|')[-3].split('_')[-1]
    except Exception as e:
        logging.error(f'merge tax 16S {e}')
    try:
        if tax_file_genome:
            df_genome = pd.read_csv(tax_file_genome)
            species_genome = df_genome.loc[0, 'Tax'].split('|')[-1].split('__')[-1]
            if species_genome in ('', '-'):
                species_genome = df_genome.loc[0,'Species ID']
            genus_genome = df_genome.loc[0, 'Tax'].split('|')[-2].split('_')[-1]
            family_genome = df_genome.loc[0, 'Tax'].split('|')[-3].split('_')[-1]
            strain_genome = df_genome.loc[0, 'Strain id']
            ani_genome = ''
            if 'fastANI' in df_genome.columns:
                ani_genome = df_genome.loc[0, 'fastANI']
            if not ani_genome:
                ani_genome = df_genome.loc[0, 'OrthoANI']
    except Exception as e:
        logging.error(f'merge tax genome {e}')

    logging.info(f'family: {family_genome}')
    logging.info(f'genus: {genus_genome}')
    logging.info(f'species: {species_genome}, related_strain: {strain_genome}')
    if genus_genome != genus_16s:
        logging.info(f'Family inconsist of genome {family_genome} and 16S {family_16s}')
    if genus_genome != genus_16s:
        logging.info(f'Genus inconsist of genome {genus_genome} and 16S {genus_16s}')
    if species_genome != species_16S:
        logging.info(f'Species inconsist of genome {species_genome} and 16S {species_16S}')
    tax_summary_dict = {'sample':name, 'family':family_genome, 'genus':genus_genome, 'species':species_genome,
                        'related_strain':strain_genome, 'ANI':ani_genome}
    pd.DataFrame.from_dict(tax_summary_dict, orient='index').T.set_index('sample').T.to_csv(outfile)

def Fasta_Taxonomy(fa, outdir, db_16s, info_16s, db_genome, info_genome, db_genome_fa, threads, name):
    out_file_list = {}
    out_file_list['json'] = {}
    out_file_list['file'] = []
    tax_file_16s = Fasta_16S_Taxonomy(fa, outdir, db_16s, info_16s, threads)
    if os.path.isfile(tax_file_16s):
        try:
            out_file_list['file'].append(tax_file_16s)
            df_tmp = pd.read_csv(tax_file_16s).astype(str)
            df_tmp.columns = df_tmp.columns.str.lower().str.replace(' ','_').str.replace('.','_',regex=False).str.replace('-','_')
            df_tmp = [v for i, v in df_tmp.to_dict(orient='index').items()]
            out_file_list['json'].update({'tax_16S':df_tmp})
        except Exception as e:
            logging.error(f'Fasta_Taxonomy {e}')
    else:
        logging.info(f'Fasta_Taxonomy not exists {tax_file_16s}')
    tax_file_genome = Fasta_genome_Taxonomy(fa, outdir, db_genome, info_genome, db_genome_fa, threads)
    if os.path.isfile(tax_file_genome):
        try:
            out_file_list['file'].append(tax_file_genome)
            df_tmp = pd.read_csv(tax_file_genome).astype(str)
            df_tmp.columns = df_tmp.columns.str.lower().str.replace(' ', '_').str.replace('.', '_', regex=False).str.replace('-', '_')
            df_tmp = [v for i, v in df_tmp.to_dict(orient='index').items()]
            out_file_list['json'].update({'tax_genome': df_tmp})
        except Exception as e:
            logging.error(f'Fasta_Taxonomy {e}')
    else:
        logging.info(f'Fasta_Taxonomy not exists {tax_file_genome}')
    merge_file = os.path.join(outdir, 'tax_summary.csv')
    merge_tax_genome_16s(tax_file_genome, tax_file_16s, merge_file, name)
    df_tmp = pd.read_csv(merge_file, index_col=0).fillna('')[name].to_dict()
    out_file_list['json'].update({'tax_summary': df_tmp})
    out_file_list['file'].append(merge_file)
    return out_file_list

if __name__ == '__main__':
    cpu_num = multiprocessing.cpu_count()
    if cpu_num > 4:
        cpu_num = cpu_num - 2
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='genome fasta file')
    parser.add_argument('-o', '--outdir', default='./', help='output dir')
    parser.add_argument('-n', '--name', default='test', help='sample name')
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
    args = parser.parse_args()


    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)
    outdir = args.outdir
    mkdir(outdir)
    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    status_report = os.path.join(outdir, 'status_report.txt')
    taskID = args.taskID
    post_pid(taskID)
    if os.path.isfile(args.db_16S + '.nsq') and os.path.isfile(args.info_16S) and os.path.isfile(args.db_genome) and os.path.isfile(args.info_genome):
        tmp_dir = os.path.join(outdir, 'tmp_dir')
        mkdir(tmp_dir)
        try:
            ## Taxonomy
            s = f'Taxonomy\tR\t'
            write_status(status_report, s)
            out_file_list_tmp = Fasta_Taxonomy(args.input, tmp_dir, args.db_16S, args.info_16S, args.db_genome, args.info_genome, args.db_genome_fa, args.threads, args.name)
            species_file = out_file_list_tmp['file'][-1]
            species = pd.read_csv(species_file, index_col=0, header=None).loc['species',1]
            copy_file(out_file_list_tmp['file'], outdir)
            with open(os.path.join(outdir, 'Taxonomy.json'), 'w') as H:
                json.dump(out_file_list_tmp['json'], H, indent=2)
            s = f'Taxonomy\tD\t'
        except Exception as e:
            logging.error(f'Taxonomy {e}')
            s = f'Taxonomy\tE\t'

        try:
            write_status(status_report, s)
            try:
                post_url(taskID, 'Taxonomy')
                post_url(taskID, '2', 'http://localhost/task/getTaskRunningStatus/')
            except Exception as e:
                logging.error(f'post_url getTaskRunningStatus {e}')
            # post_url(taskID, 'Taxonomy')
        except Exception as e:
            logging.error(f'Taxonomy status {e}')

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
    else:
        logging.error(f'something missed {args.db_16S} {args.info_16S} {args.db_genome} {args.info_genome}')




