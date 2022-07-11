

##  GCM增加三代测序数据的组装命令

> 增加了pacbio和nanopore的参数，pacbio的输入可以是fastq格式或者fasta格式

```shell
## 原来的命令
python /Bio/script/GCM_pipline.py -i {fastq1} {fastq2} -o {outdir} -n {name} -tID {taskID}

python /Bio/script/Fastq_Assemble.py -i {fastq1} {fastq2} -o {outdir} -n {name} -tID {taskID}

## 增加的命令
## 增加的命令
python /Bio/script/GCM_pipline.py -i {fastq1} {fastq2} -i {fastq} -o {outdir} -n {name} -tID {taskID} -pacbio {pacbio_fq_or_fa} -nanopore {nanopore_fq}


python /Bio/script/Fastq_Assemble.py -i {fastq1} {fastq2} -i {fastq} -o {outdir} -n {name} -tID {taskID} -pacbio {pacbio_fq_or_fa} -nanopore {nanopore_fq}

python /Bio/script/Fastq_Assemble.py -i {fastq1} {fastq2} -i {fastq} -o {outdir} -n {name} -tID {taskID} -pacbio {pacbio_fq_or_fa} -nanopore {nanopore_fq}

## 示例
python /Bio/script/Fastq_Assemble.py -i /Bio/test_data/test3.fq -i /Bio/test_data/test1.fq /Bio/test_data/test2.fq -pacbio /Bio/test_data/testQ5.scaffolds.fasta -nanopore /Bio/test_data/test4.fq -o /home/imcas/yangkai_test/test/test_assemble -n name -tID taskID 

python /Bio/script/Fastq_Assemble.py -i /Bio/test_data/test3.fq -i /Bio/test_data/test1.fq /Bio/test_data/test2.fq -pacbio /Bio/test_data/testQ5.scaffolds.fasta -nanopore /Bio/test_data/test4.fq -o /home/imcas/yangkai_test/test/test_assemble -n name -tID taskID 
```







