import os
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import logging

def cal_gc(seq, winsize=500, step=20):
    gc_ratio = []
    bin = int(len(seq)/step)
    tmp_list = [1 if i in ['g','G','c','C'] else 0 for i in seq]
    total_gc_ratio = sum(tmp_list)/len(tmp_list)
    for win in range(bin):
        win_pos = step * win
        gc_tmp = sum(tmp_list[win_pos:win_pos+winsize])
        gc_ratio.append(round(gc_tmp/winsize*100, 2))
    return total_gc_ratio, gc_ratio

def get_depth(genomecov_file, seq_id, seq_len, winsize=500, step=20):
    # try:
    total_cover = {}
    with open(genomecov_file, 'rt') as h:
        for i in h.readlines():
            tmp = i.strip().split('\t')
            if tmp[0] == seq_id:
                total_cover[int(tmp[1])-1] = float(tmp[2])

    cover_list = [total_cover[i] if i in total_cover else 0 for i in range(seq_len)]
    cover_avg = []
    for win in range(int(seq_len/step)):
        win_pos = step * win
        cover_tmp = sum([i for i in cover_list[win_pos:win_pos+winsize]])
        cover_avg.append(cover_tmp/winsize)
    return cover_avg
    # except Exception as e:
        # logging.info(e)

def plot_gcCover(gc_ratio, cover_avg, outpng):
    df = pd.DataFrame([gc_ratio, cover_avg], index=['GC content(%)', 'Sequencing depth (x)']).T
    sns.jointplot(x='GC content(%)', y='Sequencing depth (x)', data=df,
                  xlim=(20, 60), ylim=(0,200), alpha=0.2, s=10
                  )
    plt.savefig(outpng, dpi=300)
    plt.close()

def fasta2gcCover(fa, genomecov_file, outfile='test_gc_cover.png', winsize=500, step=20):
    genome_info = os.path.join(os.path.split(outfile)[0], 'assemble_summary.csv')

    len_list = []
    total_gc_ratio = []
    gc_ratio = []
    cover_avg = []
    n_num = 0
    try:
        check = 0
        for seq in SeqIO.parse(fa, 'fasta'):
            if check > 10:
                break
            len_list.append(len(seq.seq))
            total_gc_ratio_tmp , gc_ratio_tmp = cal_gc(seq.seq)
            total_gc_ratio.append(total_gc_ratio_tmp)
            gc_ratio.extend(gc_ratio_tmp)
            cover_avg.extend(get_depth(genomecov_file, seq.id, len(seq.seq)))
            n_num += len([i for i in seq.seq if i in['n','N']])
            check += 1
        plot_gcCover(gc_ratio, cover_avg, outfile)
    except Exception as e:
        logging.error(f'fasta2gcCover {e}')
    if total_gc_ratio:
        total_gc_ratio_avg = sum(total_gc_ratio)/len(total_gc_ratio)
        if len_list:
            df = pd.DataFrame(len_list).describe()
            df.loc['genome_size'] = sum(len_list)
            df.loc['GC'] = total_gc_ratio_avg
            df.loc['cover_depth']  = sum(cover_avg)/len(cover_avg)
            df.loc['n_number'] = n_num
            df = df.rename({'count':'scaffolds_num', 'mean':'len_avg', 'std':'len_std', 'min':'len_min',
                            'max':'len_max', '25%':'N25', '50%':'N50', '75%':'N75'})
            df.to_csv(genome_info, header=None)
    return genome_info

def plot_len_dis(lst, outfile, xlabel=''):
    try:
        g = sns.displot(lst)
        # g.set(xlim=(0,xmax))
        plt.xlabel(xlabel)
        plt.savefig(outfile, dpi=300)
        plt.close()
    except Exception as e:
        logging.error(f'plot_len_dis {e}')



# logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
# fa = 'test.scaffolds.fasta'
# genomecov_file = 'align.genomecov.depth.txt'
# fasta2gcCover(fa, genomecov_file)



