import os
import json
import logging
import multiprocessing
import argparse
import pandas as pd

import run_docker
from mkdir import mkdir
from post_status import post_url


def get_value(d, key_list, name):
    if type(d) is dict:
        if key_list[0] in d:
            k_tmp = key_list[1:]
            #print(k_tmp)
            if k_tmp:
                return get_value(d[key_list[0]], k_tmp, name)
            else:
                return {name : d[key_list[0]]}
        else:
            logging.info(f'{key_list[0]} not in {d.keys()}')
    elif type(d) is list:
        tmp_out = []
        for i, d_tmp in enumerate(d):
            tmp_out.append(get_value(d_tmp, key_list, name+f'_{i+1}'))
        return tmp_out

def get_all_value(json_dict, config_file, step='Fastqc'):
    with open(config_file, 'rt') as h:
        con = json.load(h)
    df = []
    for i in con[step].keys():
        tmp_list = con[step][i]['key']
        df.append(get_value(json_dict, tmp_list, i))
    df_dict = {}
    for i in df:
        if type(i) is dict:
            df_dict.update(i)
        elif type(i) is list:
            for j in i:
                df_dict.update(j)
    return df_dict

# def merge_json(json_list, name_list, step, config_file):
def merge_json(dir_list, name_list, step, config_file):
    total_json = {}
    total_json['batch'] = {}
    total_json['samples'] = []
    total_df = pd.DataFrame()
    json_list, name_list_tmp, file_fail_list, name_fail_list = get_file_list(dir_list, name_list, step+'.json')
    for i, j in enumerate(json_list):
        with open(j, 'rt') as h:
            try:
                tmp = json.load(h)
                if config_file:
                    tmp_dict = get_all_value(tmp, config_file, step)
                    tmp_df = pd.DataFrame.from_dict(tmp_dict, orient='index')
                    tmp_df.columns = [name_list_tmp[i]]
                    total_df = pd.concat([total_df, tmp_df], axis=1)
            except Exception as e:
                logging.error(f'{j} {e}')
            tmp['sampleID'] = name_list_tmp[i]
            tmp['samplePath'] = os.path.split(j)[0]
            total_json['samples'].append(tmp)
    for i, j in enumerate(file_fail_list):
        tmp = {}
        tmp['sampleID'] = name_fail_list[i]
        tmp['samplePath'] = os.path.split(j)[0]
        total_json['samples'].append(tmp)
    return total_json, total_df

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

def make_cg_tree(cg_file_list, name_list, outdir):
    logging.info('merge_cg_tree')
    df_total = pd.DataFrame()
    for i, cg_file in enumerate(cg_file_list):
        df_tmp = pd.read_csv(cg_file, index_col=0, sep='\t')
        try:
            df_tmp.index = [name_list[i]]
        except Exception as e:
            logging.error(f'make_cg_tree {e}')
        df_total = df_total.append(df_tmp)

    cg_total_file = os.path.join(outdir, 'cgMLST_total.csv')
    df_total = df_total.rename_axis(index='FILE')
    df_total.to_csv(cg_total_file, sep='\t')
    tree_file = run_r_tree(cg_total_file, outdir)
    return tree_file, cg_total_file


def get_rgi_csv(f, cutoff=90):
    df = pd.read_csv(f)
    if 'Best_Identities' in df.columns:
        df = df[df['Best_Identities']>cutoff]
    elif 'Identities' in df.columns:
        df = df[df['Identities']>cutoff]
    # df['sampleID'] = name
    return df

def merge_factor_csv(file_list, name_list, col='Drug Class'):
    df_total = pd.DataFrame()
    for i, f in enumerate(file_list):
        if os.path.isfile(f):
            df_tmp = get_rgi_csv(f)
            df_tmp = df_tmp[col].str.split('; ').apply(pd.Series).stack().value_counts()
            df_tmp.name = name_list[i]
            df_total = pd.concat([df_total, df_tmp], axis=1)
    # df_sta = df_total.groupby([col, 'sampleID'])['Contig'].count().unstack().fillna(0).rename_axis(index=None,columns=None)
    df_sta = df_total.fillna(0)
    for i in name_list:
        if i not in df_sta.columns:
            df_sta[i] = 0
    return df_sta.T

def get_file_list(dir_list, name_list, file_name):
    file_list = []
    name_list_tmp = []
    file_fail_list = []
    name_fail_list = []
    for i, d in enumerate(dir_list):
        f = os.path.join(d, file_name)
        if os.path.isfile(f):
            if not os.path.getsize(f):
                logging.info(f'{f} is empty')
                file_fail_list.append(f)
                name_fail_list.append(name_list[i])
            else:
                file_list.append(f)
                name_list_tmp.append(name_list[i])
        else:
            logging.info(f'{name_list[i]} {f} not exists')
            file_fail_list.append(f)
            name_fail_list.append(name_list[i])
    return file_list, name_list_tmp, file_fail_list, name_fail_list

def merge(dir_list, name_list, outdir='./', type='Fastqc', cofig_file=None):
    if type in ['Fastqc', 'Assembly', 'Taxonomy', 'FactorAnno', 'Typing', 'GeneAnno']:
        # json_list, name_list_tmp = get_file_list(dir_list, name_list, type+'.json')
        total_json_dict, total_df = merge_json(dir_list, name_list, type, cofig_file)
        if not total_df.empty:
            total_df.to_csv(os.path.join(outdir, type+'.csv'))

        if type in ['FactorAnno']:
            filename = 'FactorAnno_VFs.csv'
            VF_list, name_list_tmp, _, _ = get_file_list(dir_list, name_list, filename)
            out_file = os.path.join(outdir, filename)
            df_out = merge_factor_csv(VF_list, name_list_tmp, col='VFs')
            df_out.to_csv(out_file)
            total_json_dict['batch'].update({'vfdb_csv':filename})

            filename = 'FactorAnno_rgi.csv'
            rgi_list, name_list_tmp, _, _ = get_file_list(dir_list, name_list, filename)
            out_file = os.path.join(outdir, filename)
            df_out = merge_factor_csv(rgi_list, name_list_tmp)
            df_out.to_csv(out_file)
            total_json_dict['batch'].update({'rgi_csv': filename})

        if type in ['GeneAnno']:
            type_factor = 'FactorAnno'
            # json_list, name_list_factor_tmp = get_file_list(dir_list, name_list, type_factor + '.json')
            total_json_dict_tax, total_df_factor = merge_json(dir_list, name_list, type_factor, cofig_file)
            if not total_df_factor.empty:
                total_df_factor.to_csv(os.path.join(outdir, type_factor + '.csv'))
            for i, d in enumerate(total_json_dict_tax['samples']):
                total_json_dict['samples'][i].update(d)

        if type in ['Typing']:
            cgfile_list, name_list_tmp, _, _ = get_file_list(dir_list, name_list, "results_alleles.tsv")

            type_tax = 'Taxonomy'
            # json_list, name_list_tax_tmp = get_file_list(dir_list, name_list, type_tax + '.json')
            total_json_dict_tax, total_df_tax = merge_json(dir_list, name_list, type_tax, cofig_file)
            if not total_df_tax.empty:
                total_df_tax.to_csv(os.path.join(outdir, type_tax + '.csv'))
            for i, d in enumerate(total_json_dict_tax['samples']):
                total_json_dict['samples'][i].update(d)

            out_tree_file, out_cgMLST_file = None, None
            if cgfile_list:
                try:
                    out_tree_file, out_cgMLST_file =  make_cg_tree(cgfile_list, name_list_tmp, outdir)
                except Exception as e:
                    logging.error(f'make_cg_tree {e}')
            else:
                logging.info(f'no file for cgMLST')
            if out_cgMLST_file:
                total_json_dict['batch'].update({'cgMLST_file': os.path.split(out_cgMLST_file)[1]})
            if out_tree_file:
                total_json_dict['batch'].update({'cgMLST_tree': os.path.split(out_tree_file)[1]})

        outfile = os.path.join(outdir, type+'.json')
        with open(outfile, 'w') as H:
            json.dump(total_json_dict, H, indent=2)

        return outfile
    else:
        logging.info(f'{type} not supported')

if __name__ == '__main__':
    cpu_num = multiprocessing.cpu_count()
    bin_dir = os.path.split(os.path.realpath(__file__))[0]
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=True, nargs="+", help='dir list')
    parser.add_argument('-n', '--name', required=True, nargs="+", help='sample name list')
    parser.add_argument('-s', '--step', required=True, nargs="+", help='types for merge')
    parser.add_argument('-t', '--threads', default=cpu_num, help='threads')
    parser.add_argument('-o', '--outdir', default='./', help='output dir')
    parser.add_argument('-tID', '--taskID', default='', help='task ID for report status')
    parser.add_argument('-con', '--config', default=os.path.join(bin_dir, 'merge_json_config.json'), help='task ID for report status')
    parser.add_argument('-debug', '--debug', action='store_true')
    args = parser.parse_args()

    if len(args.input) != len(args.name):
        sys.exit('input dirs muster equal names')

    os.environ['NUMEXPR_MAX_THREADS'] = str(args.threads)
    outdir = args.outdir
    mkdir(outdir)
    logfile = os.path.join(outdir, 'log')
    logging.basicConfig(level=logging.INFO, filename=logfile, format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    dir_list, name_list = args.input, args.name
    for step in args.step:
        try:
            merge(dir_list, name_list, outdir=args.outdir, type=step, cofig_file=args.config)
        except Exception as e:
            logging.error(f'merge {e}')
    taskID = args.taskID
    try:
        post_url(taskID, '2', 'http://localhost/task/getMergeStatus/')
    except Exception as e:
        logging.error(f'post_url getTaskRunningStatus {e}')

