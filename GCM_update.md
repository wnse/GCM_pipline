## 细菌流程使用方法

### 全流程（WholeGenome）

- /Bio/script/GCM_pipline.py

- 两种使用方法
  
  1. 命令行中参数不加-type fa
     
     - 产出结果：Fastqc.json/Assembly.json/Taxonomy.json/Typing.json/FactorAnno.json/GeneAnno.json
     - 数据格式：.fastq/.fq/.fastq.gz/.fq.gz
     - 具体参数见以下说明：`单独二代数据`，`单独三代数据`，`二代+三代数据`
  
  2. 命令行中参数加-type fa
     
     - 产出结果：Taxonomy.json/Typing.json/FactorAnno.json/GeneAnno.json
     - 数据格式：.fasta/.fa
     - 参数：

```shell
python /Bio/script/GCM_pipline.py 
-i {fasta} 
-o {outdir} 
-n {name} 
-tID {taskID}
-type fa
```

### 数据质控（Fastqc）

- /Bio/script/Fastq_QC.py
- 产出结果：Fastqc.json
- 数据格式：.fastq/.fq/.fastq.gz/.fq.gz
- 参数

```shell
## 双端测序数据
python /Bio/script/Fastq_QC.py 
-i {fastq1} {fastq2} 
-o {outdir} 
-n {name} 
-tID {taskID}

## 单端测序数据
python /Bio/script/Fastq_QC.py 
-i {fastq}
-o {outdir} 
-n {name} 
-tID {taskID}
```

### 基因组拼接（Assembly）

- /Bio/script/Fastq_Assemble.py
- 产出结果：Fastqc.json/Assembly.json
- 具体参数见以下说明：`单独二代数据`，`单独三代数据`，`二代+三代数据`

### 物种鉴定（Taxonomy）

- /Bio/script/Fasta_Taxonomy.py
- 产出结果：Taxonomy.json
- 数据格式：.fasta/.fa
- 参数

```shell
python /Bio/script/Fasta_Taxonomy.py 
-i {fasta}  
-o {outdir}
-n {name}
-tID {taskID}
```

### 分型分析（Typing）

- /Bio/script/Fasta_Typing.py
- 产出结果：Taxonomy.json/Typing.json
- 数据格式：.fasta/.fa
- 参数

```shell
python /Bio/script/Fasta_Typing.py 
-i {fasta}  
-o {outdir}
-n {name}
-tID {taskID}
```

### 耐药毒力基因分析（FactorAnno）

- /Bio/script/Fasta_FactorAnno.py
- 产出结果：FactorAnno.json
- 数据格式：.fasta/.fa
- 参数

```shell
python /Bio/script/Fasta_FactorAnno.py 
-i {fasta}  
-o {outdir}
-n {name}
-tID {taskID}
```

### 基因注释（GeneAnno）

- /Bio/script/Fasta_Annotation.py
- 产出结果：FactorAnno.json/GeneAnno.json
- 数据格式：.fasta/.fa
- 参数

```shell
python /Bio/script/Fasta_Annotation.py 
-i {fasta}  
-o {outdir}
-n {name}
-tID {taskID}
```

### 单独二代数据

- 适用全流程和基因组拼接
- 数据格式：.fastq/.fq/.fastq.gz/.fq.gz
- 参数

```shell
## 双端测序数据
python /Bio/script/GCM_pipline.py 
-i {fastq1} {fastq2} 
-o {outdir} 
-n {name} 
-tID {taskID}


## 单端测序数据
python /Bio/script/GCM_pipline.py 
-i {fastq}
-o {outdir} 
-n {name} 
-tID {taskID}


## 一个样本多个二代测序数据
## 双端测序数据和单端测序数据皆可
python /Bio/script/GCM_pipline.py 
-i {fastq1_1} {fastq1_2} 
-i {fastq2_1} {fastq2_2}
-i {fastq3}
-o {outdir} 
-n {name} 
-tID {taskID}
```

### 单独三代数据

- 适用全流程和基因组拼接
- 数据格式：.fastq/.fq/.fastq.gz/.fq.gz/.fasta/.fa/.fasta.gz/.fa.gz
- 参数

```shell
## pacbio测序数据
python /Bio/script/GCM_pipline.py 
-pacbio {fastq}/{fasta}
-o {outdir} 
-n {name} 
-tID {taskID}

## nanopore测序数据
python /Bio/script/GCM_pipline.py 
-nanopore {fastq}/{fasta}
-o {outdir} 
-n {name} 
-tID {taskID}

## pacbio + nanopore测序数据
python /Bio/script/GCM_pipline.py 
-pacbio {fastq}/{fasta}
-nanopore {fastq}/{fasta}
-o {outdir} 
-n {name} 
-tID {taskID}
```

### 二代+三代数据

- 适用全流程和基因组拼接
- 二代数据格式：.fastq/.fq/.fastq.gz/.fq.gz
- 三代数据格式：.fastq/.fq/.fastq.gz/.fq.gz/.fasta/.fa/.fasta.gz/.fa.gz
- 参数

```shell
## 单独二代参数 + 单独三代参数
## 即 
## -i 为二代数据
## -pacbio 和 -nanopore 分别为三代pacbio与nanopore数据
## 可根据输入数据选择
python /Bio/script/GCM_pipline.py 
-i {fastq1} {fastq2} 
-pacbio {fastq}/{fasta}
-nanopore {fastq}/{fasta}
-o {outdir} 
-n {name} 
-tID {taskID}
```
