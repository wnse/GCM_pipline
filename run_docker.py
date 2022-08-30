import os
import re
import unzip_file
import logging
from mkdir import mkdir

def run_docker(img_name, all_data_path, cmd):
    s = set(os.path.split(i)[0] for i in all_data_path if i)
    volume = []
    for i in s:
        volume.append(f"-v {i}:{i}")
    cmd_run = f'sudo docker run --rm {" ".join(volume)} {img_name} {cmd}'
    logging.info(f'{cmd_run}')
    return os.system(cmd_run)
    # return f"{cmd_run}\n\n"

def fastqc_docker(fq_path_list, outdir, threads):
    img_name = 'registry.servicemgr.bjwsws:5000/fastqc:latest'
    all_data_path = fq_path_list + [outdir]
    qc_dir = os.path.join(outdir, 'fastqc_dir')
    mkdir(qc_dir)
    logfile = os.path.join(outdir, 'fastqc.log')
    cmd = (f'/usr/local/bin/fastqc -t {threads} -o {qc_dir} {" ".join(fq_path_list)} '
           f'>{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0
    return cmd_out, qc_dir


def trim_docker_pe(fq1, fq2, outdir, threads):
    img_name = 'registry.servicemgr.bjwsws:5000/trimmomatic_sickle:latest'
    all_data_path = [fq1, fq2, outdir]
    trim_dir = os.path.join(outdir, 'trimmomatic_dir')
    mkdir(trim_dir)
    fq1_trime = os.path.join(trim_dir, "trim_1.fastq")
    fq2_trime = os.path.join(trim_dir, "trim_2.fastq")
    fq1_unpaired = os.path.join(trim_dir, "forward_unpaired.fastq")
    fq2_unpaired = os.path.join(trim_dir, "reverse_unpaired.fastq")
    logfile = os.path.join(outdir, "trimmomatic.log")
    cmd = (f'java -jar /software/applications/Trimmomatic/0.36/trimmomatic-0.36.jar  PE -threads {threads} '
           f'{fq1} {fq2} {fq1_trime} {fq1_unpaired} {fq2_trime} {fq2_unpaired} '
           f'ILLUMINACLIP:/software/applications/Trimmomatic/0.36/adapters/TruSeq3-PE.fa:2:30:10 '
           f'SLIDINGWINDOW:5:20 MINLEN:20 '
           f'>{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0, 0
    return cmd_out, fq1_trime, fq2_trime

def trim_docker_se(fq, outdir, threads):
    img_name = 'registry.servicemgr.bjwsws:5000/trimmomatic_sickle:latest'
    all_data_path = [fq, outdir]
    trim_dir = os.path.join(outdir, 'trimmomatic_dir')
    mkdir(trim_dir)
    fq_trime = os.path.join(trim_dir, "trim.fastq")
    logfile = os.path.join(outdir, "trimmomatic.log")
    cmd = (f'java -jar /software/applications/Trimmomatic/0.36/trimmomatic-0.36.jar  SE -threads {threads} '
           f'{fq} {fq_trime} '
           f'ILLUMINACLIP:/software/applications/Trimmomatic/0.36/adapters/TruSeq3-SE.fa:2:30:10 '
           f'SLIDINGWINDOW:5:20 MINLEN:20 '
           f'>{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0
    return cmd_out, fq_trime

def musket_docker(fq_list, outdir, threads):
    img_name = 'registry.servicemgr.bjwsws:5000/srctools:latest'
    all_data_path = fq_list + [outdir]
    musket_dir = os.path.join(outdir, 'musket_dir')
    mkdir(musket_dir)
    logfile = os.path.join(outdir, "musket.log")
    cmd = (f'/BioBin/musket-1.1/bin/musket -inorder -p {threads} -omulti {os.path.join(musket_dir, "trime_corrected")} '
           f'{" ".join(fq_list)} >{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0, 0
    fq_cor_list = []
    for i, f in enumerate(fq_list):
        fq_cor = os.path.join(musket_dir, f"trime_corrected.{i+1}.fastq")
        # fq_cor_2 = os.path.join(musket_dir, "trime_corrected.2.fastq")
        os.system(f'sudo mv {os.path.join(musket_dir, f"trime_corrected.{i}")} {fq_cor}')
        # os.system(f'sudo mv {os.path.join(musket_dir, "trime_corrected.1")} {fq_cor_2}')
        fq_cor_list.append(fq_cor)
    return cmd_out, fq_cor_list


def unzip_fq_file(file, outdir):
    ft = unzip_file.file_type(file)
    if ft == 'FQ':
        return file
    if ft == 'gz':
        if os.path.splitext(file)[1] == '.gz':
            file_name = os.path.join(outdir, 'tmp.fq')
            unzip_file.un_gz(file, file_name)
            return file_name
        else:
            return f"unkown file type {ft}"
    return "unkown file type"

def jellyfish_docker_count(fq, outdir, threads):
    # file = unzip_fq_file(fq, outdir)
    img_name = 'registry.servicemgr.bjwsws:5000/jellyfish:latest'
    all_data_path = [fq, outdir]
    jellyfish_dir = os.path.join(outdir, 'jellyfish_dir')
    mkdir(jellyfish_dir)
    outfile = os.path.join(jellyfish_dir, 'jellyfish.jf')
    logfile = os.path.join(outdir, "jellyfish_count.log")
    cmd = (f'jellyfish count -C -m 17 -s 1000000000 -t {threads} {fq} -o {outfile} >{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0
    return cmd_out, outfile

def jellyfish_docker_histo(jf_file, outdir, threads):
    img_name = 'registry.servicemgr.bjwsws:5000/jellyfish:latest'
    all_data_path = [jf_file, outdir]
    jellyfish_dir = os.path.join(outdir, 'jellyfish_dir')
    mkdir(jellyfish_dir)
    outfile = os.path.join(jellyfish_dir, 'jellyfish.histo')
    logfile = os.path.join(os.path.split(outfile)[0], 'jellyfish_histo.log')
    cmd = (f'jellyfish histo -t {threads} {jf_file} > {outfile} 2>{logfile}')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd, 0
    return cmd, outfile

def genomescope_docker(histo_file, outdir, read_length):
    img_name = 'registry.servicemgr.bjwsws:5000/genomescope:latest'
    all_data_path = [histo_file, outdir]
    genomescope_dir = os.path.join(outdir, 'genomescope_dir')
    mkdir(genomescope_dir)
    outfile = os.path.join(genomescope_dir, "summary.txt")
    outpng = os.path.join(genomescope_dir, 'plot.png')
    logfile = os.path.join(outdir, 'genomescope.log')
    cmd = (f'genomescope.R {histo_file} 17 {read_length} {genomescope_dir} >{logfile} 2>&1')
    cmd_out =  run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd, 0, 0
    return cmd, outfile, outpng



def canu_docker(outdir, pacbio, nanopore):
    all_data_path = [outdir]
    canu_outdir = os.path.join(outdir,'canu_dir')
    if pacbio:
        all_data_path.append(pacbio)
        x_para = f'-pacbio {pacbio}'
    elif nanopore:
        all_data_path.append(nanopore)
        x_para = f'-nanopore {nanopore}'
    logfile = os.path.join(outdir, "canu.log")
    img_name = 'staphb/canu:latest'
    cmd = (f'canu -p canu -d {canu_outdir} genomeSize=4m minInputCoverage=0 stopOnLowCoverage=0 {x_para} >{logfile} 2>&1')
    cmd_out =  run_docker(img_name, all_data_path, cmd)
    outfile = os.path.join(canu_outdir, 'canu.contigs.fasta')
    if cmd_out:
        return cmd, 0
    return cmd, outfile


# def spades_docker(PE_fq_list=None, SE_fq = None, outdir='./', threads=1, k_mer=[21,33,55]):
def spades_docker(outdir, PE_fq_list, SE_fq_list=[], pacbio_list=None, nanopore_list=None, threads=1, k_mer=[21,33,55]):
    k_mer_str = 'spades_' + '_'.join([str(i) for i in k_mer])
    k_mer_dir = os.path.join(outdir, k_mer_str)
    mkdir(k_mer_dir)
    log_file = os.path.join(outdir, k_mer_str+".log")
    img_name = 'registry.servicemgr.bjwsws:5000/assembletools:latest'
    all_data_path = [outdir]
    pe_input = ''
    se_input = ''
    if PE_fq_list:
        if len(PE_fq_list) == 1:
            pe_fq = PE_fq_list[0]
            if os.path.isfile(pe_fq[0]) and os.path.isfile(pe_fq[1]):
                pe_input += f' -1 {pe_fq[0]} -2 {pe_fq[1]}'
                all_data_path = all_data_path + pe_fq     
        else: 
            for pe_fq  in PE_fq_list:
                if len(pe_fq) == 2:
                    if os.path.isfile(pe_fq[0]) and os.path.isfile(pe_fq[1]):
                        pe_input += f' --pe1-1 {pe_fq[0]} --pe1-2 {pe_fq[1]}'
                        all_data_path = all_data_path + pe_fq
                    else:
                        logging.error(f'fastq file path not correct: {pe_fq}')

                else:
                    logging.error(f'fastq file number not correct: {pe_fq}')
        
    if SE_fq_list:
        for se_fq in SE_fq_list:
            if os.path.isfile(se_fq):
                se_input += f' -s {se_fq}'
                all_data_path = all_data_path + [se_fq]
            else:
                logging.error(f'fastq file path not correct: {se_fq}')
    pacbio_fq = ''
    nanopore_fq = ''
    if pacbio_list:
        for fq in pacbio_list:
            if os.path.isfile(fq):
                pacbio_fq += f' --pacbio {fq}'
                all_data_path = all_data_path + [fq]
            else:
                logging.error(f'fastq file path not correct: {fq}')
    if nanopore_list:
        for fq in nanopore_list:
            if os.path.isfile(fq):
                nanopore_fq += f' --nanopore {fq}'
                all_data_path = all_data_path + [fq]
            else:
                logging.error(f'fastq file path not correct: {fq}')

    if pe_input or se_input:
        cmd = (f'python /BioBin/SPAdes-3.13.0-Linux/bin/spades.py --careful --sc --disable-gzip-output '
               f'{pe_input} {se_input} {pacbio_fq} {nanopore_fq} '
               f'-k {",".join([str(i) for i in k_mer])} -o {k_mer_dir} -t {threads} '
               f'>{log_file} 2>&1')
    #
    # cmd = (f'python /BioBin/SPAdes-3.13.0-Linux/bin/spades.py --careful --sc --disable-gzip-output '
    #        f'-1 {fq1} -2 {fq2} -k {",".join([str(i) for i in k_mer])} -o {k_mer_dir} -t {threads} '
    #        f'>{log_file} 2>&1')
        cmd_out = run_docker(img_name, all_data_path, cmd)
        print(cmd_out)
        if cmd_out:
            return cmd_out, 0
        scaffolds_fasta = os.path.join(k_mer_dir, "scaffolds.fasta")
        if os.path.isfile(scaffolds_fasta):
            return cmd_out, scaffolds_fasta
        return cmd_out, 0
    return 0, 0

def bwa_align_docker(fq_list, scaffolds_fasta, outdir, threads=1, nanopore=None, pacbio=None):
    bwa_align_dir = os.path.join(outdir, 'bwa_align_dir')
    mkdir(bwa_align_dir)
    out_bam_file = os.path.join(bwa_align_dir, 'align.sorted.bam')
    out_genomecov_file = os.path.join(bwa_align_dir, 'align.genomecov.depth.txt')
    out_flagstat_file = os.path.join(bwa_align_dir, 'align.flagstat.txt')
    log_file = os.path.join(outdir, 'bwa_align.log')
    sam_file = out_bam_file + ".sam"

    if fq_list:
        all_data_path = fq_list + [scaffolds_fasta, outdir]
        img_name = 'registry.servicemgr.bjwsws:5000/bwa:latest'
        cmd = (f'/BioBin/bwa/bwa index {scaffolds_fasta} >>{log_file} 2>>{log_file}')
        cmd_out = run_docker(img_name, all_data_path, cmd)
        if cmd_out:
            return "bwa index error"

        cmd = (f'/BioBin/bwa/bwa mem {scaffolds_fasta} {" ".join(fq_list)} -t {threads} > {sam_file} 2>>{log_file}')
        cmd_out = run_docker(img_name, all_data_path, cmd)
        if cmd_out:
            return "bwa mem error"
    elif nanopore or pacbio:
        data = None
        if pacbio:
            x_para = 'map-pb'
            data = pacbio
        if nanopore:
            x_para = 'map-ont'
            data = nanopore

        all_data_path = [data, scaffolds_fasta, outdir]
        img_name = 'staphb/minimap2:latest'
        cmd = (f'minimap2 -ax {x_para} {scaffolds_fasta} {data} -t {threads} -o {sam_file} 2>>{log_file}')
        cmd_out = run_docker(img_name, all_data_path, cmd)
        if cmd_out:
            return "minimap2 error"
    else:
        return 'no data input'

    img_name = 'registry.servicemgr.bjwsws:5000/bwa:latest'
    cmd = (f'/BioBin/samtools/samtools sort {sam_file} -o {out_bam_file} 2>>{log_file}')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return "bwa mem error"

    cmd = (f'/BioBin/bedtools2/bin/bedtools genomecov -ibam {out_bam_file} -d  >{out_genomecov_file} 2>>{log_file}')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return "bedtools genomecov error"

    cmd = (f'/BioBin/samtools/samtools flagstat {out_bam_file} >{out_flagstat_file} 2>>{log_file}')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return "samtools flagsta error"
    return cmd_out, out_bam_file, out_genomecov_file, out_flagstat_file

def picard_CIAM_docker(sortBam, outdir):
    img_name = 'registry.servicemgr.bjwsws:5000/picard:latest'
    all_data_path = [sortBam, outdir]
    picard_dir = os.path.join(outdir, 'picard_dir')
    mkdir(picard_dir)
    outTxt = os.path.join(picard_dir, 'insert_size_metrics.txt')
    outPdf = os.path.join(picard_dir, 'insert_size_histogram.pdf')
    log_file = os.path.join(outdir, 'picard.log')
    cmd = (f'java -jar /usr/picard/picard.jar CollectInsertSizeMetrics '
           f'I={sortBam} O={outTxt} H={outPdf} >{log_file} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return "picard CIAM error"
    return cmd_out, outTxt, outPdf



def checkm_lineage_docker(inputdir, outputdir, threads):
    img_name = 'registry.servicemgr.bjwsws:5000/checkm:latest'
    all_data_path = [inputdir, outputdir]
    checm_dir = os.path.join(outputdir, 'checkm_dir')
    mkdir(checm_dir)
    out_file = os.path.join(outputdir, 'checkm_out.txt')
    log_file = os.path.join(outputdir, 'checkm.log')
    cmd = (f'checkm lineage_wf -t {threads} -x fasta -f {out_file} {inputdir} {checm_dir} >{log_file} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0
    return cmd_out, out_file


def rnammer_docker(fa, outputdir, s='bac', m='ssu'):
    img_name = 'registry.servicemgr.bjwsws:5000/rnammer:v1.2'
    all_data_path = [fa, outputdir]
    rnammer_dir = os.path.join(outputdir, 'rnammer_dir')
    mkdir(rnammer_dir)
    xml_file = os.path.join(rnammer_dir, 'RNAmmer.xml')
    fasta_file = os.path.join(rnammer_dir, 'RNAmmer.fasta')
    report_file = os.path.join(rnammer_dir, 'RNAmmer.report')
    gff_file = os.path.join(rnammer_dir, 'RNAmmer.gff')
    log_file = os.path.join(outputdir, 'RNAmmer.log')
    cmd = (f'rnammer -S {s} -multi -m {m} -xml {xml_file} -f {fasta_file} '
           f'-h {report_file} -gff {gff_file} {fa} >{log_file} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return "rnammer error"
    return cmd_out, fasta_file

def blastn_docker(query, db, out, evalue=1e-5, max_hsps=1, max_target_seqs=500, outfmt=6, threads=1, img_name='registry.servicemgr.bjwsws:5000/blast:2.9.0'):
    all_data_path = [query, db, out]
    log_file = os.path.join(os.path.split(out)[0], 'blastn.log')
    cmd = (f'/BioBin/blast/bin/blastn -query {query} -db {db} -out {out} -outfmt {outfmt} -evalue {evalue} '
           f'-max_hsps {max_hsps} -max_target_seqs {max_target_seqs} -num_threads {threads} >{log_file} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return "blastn error"
    return cmd_out

def mash_dist_docker(query, db, out, evalue=0.05, threads=1):
    img_name = 'mash:latest'
    all_data_path = [query, db, out]
    log_file = os.path.join(os.path.split(out)[0], 'mash_dist.log')
    cmd = (f'mash dist -v {evalue} -p {threads} {db} {query} | sort -gk3 >{out}  2>{log_file}')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return "mash dist error"
    return cmd_out

def fastANI_docker(fa1, fa2, out, threads=1):
    img_name = 'registry.servicemgr.bjwsws:5000/fastani:20200713'
    all_data_path = [fa1, fa2, out]
    log_file = os.path.join(os.path.split(out)[0], 'fastANI.log')
    cmd = (f'/usr/local/bin/fastANI -q {fa1} -r {fa2} -t {threads} -o {out} >{log_file} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return "fastANI error"
    return cmd_out

def OrthoANI_docker(fa1, fa2, out, threads=1):
    img_name = 'registry.servicemgr.bjwsws:5000/orthoani:latest'
    all_data_path = [fa1, fa2, out]
    log_file = os.path.join(os.path.split(out)[0], 'OrthoANI.log')
    cmd = (f'java -jar /BioBin/OrthoANI/OAT_cmd.jar -fasta1 {fa1} -fasta2 {fa2} -num_threads {threads} '
           f'-blastplus_dir /BioBin/ncbi-blast-2.10.0+/bin -method ani >{out} 2>{log_file}')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return "OrthoANI error"
    return cmd_out

def OAU_docker(fa1, fa2, out, threads=1):
    img_name = 'registry.servicemgr.bjwsws:5000/aou:20200713'
    all_data_path = [fa1, fa2, out]
    log_file = os.path.join(os.path.split(out)[0], 'OAU.log')
    cmd = (f'java -jar /BioBin/OAU.jar -f1 {fa1} -f2 {fa2} -o {out} -n {threads} '
           f'-u /BioBin/usearch11.0.667_i86linux32 >{log_file} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return "OAU error"
    return cmd_out

def ANI_docker(fa1, fa2, out):
    img_name = 'registry.servicemgr.bjwsws:5000/animaster:20200817'
    all_data_path = [fa1, fa2, out]
    log_file = os.path.join(os.path.split(out)[0], 'ANI.log')
    cmd = (f'perl /BioBin/ANI/ANI.pl --fd /BioBin/ANI/blast-2.2.23/bin/formatdb '
           f'-bl /BioBin/ANI/blast-2.2.23/bin/blastall --qr {fa1} --sb {fa2} --od {os.path.split(out)[0]} >{out} 2>{log_file}')
    # cmd_out = run_docker(img_name, all_data_path, cmd)
    # if cmd_out != 0:
    #     return "ANI error"
    return run_docker(img_name, all_data_path, cmd)

def trf_docker(fa, outdir):
    img_name = 'registry.servicemgr.bjwsws:5000/repeatmasker:latest'
    all_data_path = [fa, outdir]
    trf_dir = os.path.join(outdir, 'trf')
    mkdir(trf_dir)
    log_file = os.path.join(outdir, 'trf.log')
    cmd = (f'sh -c "cd {trf_dir}; /usr/local/bin/trf {fa} 2 7 7 80 10 50 500 -f -d -m >{log_file} 2>&1"')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    # if cmd_out != 0:
    #     return "trf error"
    return cmd_out

def pilercr_docker(fa, outfile):
    img_name = 'registry.servicemgr.bjwsws:5000/metatools:lite'
    all_data_path = [fa, outfile]
    log_file = os.path.join(os.path.split(outfile)[0], 'pilercr.log')
    cmd = (f'/BioBin/pilercr1.06/pilercr -in {fa} -out {outfile} >{log_file} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    return cmd_out

def trnascan_docker(fa, outdir):
    img_name = 'registry.servicemgr.bjwsws:5000/trnascan:v2.0'
    all_data_path = [fa, outdir]
    log_file = os.path.join(outdir, 'trnascan.log')
    trnascan_dir = os.path.join(outdir, 'trnascan')
    mkdir(trnascan_dir)
    confile = '/usr/local/bin/tRNAscan-SE.conf'
    outfile = os.path.join(trnascan_dir, 'tRNAscan.out')
    stafile = os.path.join(trnascan_dir, 'tRNAscan.sta')
    cmd = (f'tRNAscan-SE -qQ -Y -o {outfile} -m {stafile} -c {confile} -B {fa} -out >{log_file} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    return cmd_out

def Rfam_docker(fa, outdir, db_path):
    img_name = 'registry.servicemgr.bjwsws:5000/infernal:latest'
    all_data_path = [fa, outdir, db_path]
    logfile = os.path.join(outdir, 'esl-seqstat.log')
    outfile = os.path.join(outdir, 'esl-seqstat.out')
    cmd = f'esl-seqstat {fa} >{outfile} 2>{logfile}'
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out == 0 and os.path.isfile(outfile):
        with open(outfile, 'rt') as h:
            seqstat = [i.strip().split(':')[1] for i in h.readlines() if re.match('Total\s+#\s+residues:', i)][-1]
        if seqstat:
            num = round(float(seqstat)*2/1000000, 4)
            outfile = os.path.join(outdir, 'cmscan.out')
            logfile = os.path.join(outdir, 'cmscan.log')
            cmd = (f'cmscan -Z {num} --cut_ga --rfam --nohmmonly --fmt 2 '
                   f'--clanin {os.path.join(db_path, "Rfam", "Rfam.clanin")} --tblout {outfile} '
                   f'{os.path.join(db_path, "Rfam", "Rfam.cm")} {fa} >{logfile} 2>&1')
            cmd_out = run_docker(img_name, all_data_path, cmd)
            return cmd_out
    return cmd_out

def prodigal_docker(fa, outdir):
    img_name = 'registry.servicemgr.bjwsws:5000/metatools:lite'
    all_data_path = [fa, outdir]
    prodigal_dir = os.path.join(outdir, 'prodigal_dir')
    mkdir(prodigal_dir)
    logfile = os.path.join(outdir, 'prodigal.log')
    outAAfile = os.path.join(prodigal_dir, 'scaffolds_gene.faa')
    outNAfile = os.path.join(prodigal_dir, 'scaffolds_gene.fna')
    outGFFfile = os.path.join(prodigal_dir, 'scaffolds_gene.gff')

    cmd = (f'/BioBin/Prodigal/bin/prodigal -f gff -g 11 -p single -m '
           f'-a {outAAfile} -d {outNAfile} -o {outGFFfile} -i {fa} >{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0, 0
    return cmd_out, prodigal_dir, outAAfile, outGFFfile

def pfam_docker(faa, outfile, db_path, threads=1):
    img_name = 'pfam:latest'
    all_data_path = [faa, outfile, db_path]
    logfile = os.path.join(os.path.split(outfile)[0], 'pfam.log')
    cmd = (f'perl /usr/local/bin/pfam_scan.pl -fasta {faa} -dir {db_path} '
           f'-outfile {outfile} -cpu {threads} >{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    return cmd_out

def diamond_docker(faa, outfile, db_path, threads=1, id=40, query_cover=40, subject_cover=40):
    img_name = 'registry.servicemgr.bjwsws:5000/diamond:new'
    all_data_path = [faa, outfile, db_path]
    logfile = outfile + '.log'
    cmd = (f'diamond blastp -e 1e-5 -k 1 --max-hsps 1 '
           f'--id {id} --query-cover {query_cover} --subject-cover {subject_cover} '
           f'-p {threads} -d {db_path} -q {faa} -o {outfile} >{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    return cmd_out

def antismash_docker(fa, outdir, threads=1):
    img_name = 'registry.servicemgr.bjwsws:5000/antismash:v4.0.2'
    antisamsh_dir = os.path.join(outdir, 'antismash')
    all_data_path = [fa, outdir]
    logfile = os.path.join(outdir, 'antismash.log')
    cmd = (f'antismash --taxon bacteria --input-type nucl -c {threads} '
           f'--clusterblast {fa} --outputfolder {antisamsh_dir} >{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0, 0
    return cmd_out, antisamsh_dir

def rgi_docker(fa, outdir, threads=1):
    img_name = 'rgi:latest'
    all_data_path = [fa, outdir]
    rgi_dir = os.path.join(outdir, 'rgi_dir')
    mkdir(rgi_dir)
    logfile = os.path.join(outdir, 'rgi.log')
    outfile = os.path.join(rgi_dir, 'rgi.out')
    cmd = (f'sh -c "export PATH=$PATH:/rgi/bin; '
           f'rgi main -t contig --include_loose -a DIAMOND -d chromosome '
           f'-i {fa} -o {outfile} -n {threads} >{logfile} 2>&1"')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0
    return cmd_out, outfile+'.txt'

def chewBBACA_docker(fadir, outdir, schema_dir, training_file='', threads=1):
    img_name = 'chewbbaca:v1'
    all_data_path = [fadir, outdir, schema_dir, training_file]
    logfile = os.path.join(outdir, 'chewbbaca.log')
    chewbbaca_dir = os.path.join(outdir, 'chewbbaca')
    mkdir(chewbbaca_dir)
    ptf = ''
    if training_file:
        ptf = f'--ptf {training_file}'
    cmd = (f'chewBBACA.py AlleleCall -i {fadir} -o {chewbbaca_dir} -g {schema_dir} {ptf} --cpu {threads} >{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    outfile = os.path.join(chewbbaca_dir, os.listdir(chewbbaca_dir)[0], 'results_alleles.tsv')
    if cmd_out:
        return cmd_out, 0
    return cmd_out, outfile

def SeqSero2_docker(fa, outdir, threads=1):
    img_name = 'seqsero2:1.2.1'
    all_data_path = [fa, outdir]
    logfile = os.path.join(outdir, 'SeqSero2.log')
    SeqSero2_dir = os.path.join(outdir, 'SeqSero2')
    outfile = os.path.join(SeqSero2_dir, 'SeqSero_result.txt')
    mkdir(SeqSero2_dir)
    cmd = (f'/NGStools/SeqSero2-1.2.1/bin/SeqSero2_package.py -m k -t 4 '
           f'-i {fa} -d {SeqSero2_dir} -p {threads} >{logfile} 2>&1')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0
    return cmd_out, outfile

def mlst_docker(fa, outdir, threads=1):
    img_name = 'mlst:2.19'
    all_data_path = [fa, outdir]
    logfile = os.path.join(outdir, 'mlst.log')
    mlst_out = os.path.join(outdir, 'Typing_mlst.out')
    cmd = (f'/Bio/mlst-2.19.0/bin/mlst {fa} --threads {threads} >{mlst_out} 2>{logfile}')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    if cmd_out:
        return cmd_out, 0
    return cmd_out, mlst_out

def r_docker(*args):
    img_name = 'rbase'
    all_data_path = args[:-1]
    # logfile = os.path.join(outdir, 'r.log')
    cmd = (f'{" ".join(args)}')
    cmd_out = run_docker(img_name, all_data_path, cmd)
    return cmd_out


