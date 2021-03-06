****
## 目录
* [Tools](#Tools)
  * [数据下载](#数据下载)
    * [sra-tools](#sra-tools)
    * [gdc-tranfer-tool](#gdc-tranfer-tool)
    * [aspera](#aspera)
  * [辅助工具](#辅助工具)
    * [GNU-parallel](#GNU-parallel)
  * [质量控制](#质量控制)
    * [fastqc](#fastqc)
    * [cutadapter](#cutadapter)
    * [trim-galore](#trim-garole)
    * [trimmomatic](#trimmomatic)
    * [fastx-toolkit](#fastx-toolkit)
  * [基因组比对](#基因组比对)
    * [bowtie](#bowtie)
    * [bowtie2](#bowtie2)
    * [bwa](#bwa)
    * [STAR](#STAR)
    * [Kallisto](Kallisto)
  * [峰值探测](#峰值探测)
    * [MACS2](#MACS2)
  * [motif分析](#motif分析)
    * [Homer](#Homer)
  * [可视化](#可视化)
    * [deeptools*](#deeptools*)
    * [igv](#igv)


# Tools

### 数据下载

##### sra-tools

&emsp;&emsp;常用命令：`prefetch`、`fastq-dump`

&emsp;&emsp;SRA是NIH的高通量测序数据的主要档案，是国际核苷酸序列数据库协作（INSDC）的一部分，包括NCBI SRA，欧洲生物信息学研究所（EBI）和日本DNA数据库（DDBJ），提交给三个组织中任何一个的数据在它们之间共享。[SRA-tools-download-documentation](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)

* `prefetch`

  `prefetch`就是将SRA file从数据库下载到本地的基本命令，也是SRAtools中最常用的一个命令，有以下三种下载形式：

  ```bash
  prefetch [options] <path/SRA file | path/kart file> [<path/file> ...]
  prefetch [options] <SRA accession>
  prefetch [options] --list <kart_file>
  ```

  其中有几个参数需要注意：

  * **`-X`** 该参数表示下载数据的上限（单位是KB），默认是20G, 超过20G的不下载
  * **`-O`**该参数表示下载文件的路径

* `fastq-dump`

  默认下载的数据是SRR格式，我们需要将其解压缩出来，所以需要用到`fastq-dump`，有以下两种用法：

  ```bash
  fastq-dump [options] <path/file> [<path/file> ...]
  fastq-dump [options] <accession>  # 这个相当于下载、解压一步到位
  ```

  以下几个参数需要注意：

  * **`-O`**该参数表示下载文件的路径
  * **`--split-files`** 如果是双端的话，必须注明该参数

##### gdc-tranfer-tool

&emsp;&emsp;常用命令：`./gdc-client download 22a29915-6712-4f7a-8dba-985ae9a1f005`

&emsp;&emsp;Genomic Data Commons（GDC）是美国国家癌症研究所（NCI）的研究项目。GDC的使命是为**癌症研究界(cancer research community)**提供统一的数据存储库，以便在癌症基因组研究中共享数据，以支持精准医学。[gdc-documentation](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/)

##### aspera

```shell
# 安装与使用
wget https://download.asperasoft.com/download/sw/connect/3.8.1/ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz
tar -xvf ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.tar.gz
./ibm-aspera-connect-3.8.1.161274-linux-g2.12-64.sh  # installation dir:~/.aspera
export PATH=$PATH:~/.aspera/connect/bin/  # temp
ascp -i ~/asperaweb_id_dsa.openssh anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByExp/sra/SRX/SRX189/SRX189773/SRR576933/SRR576933.sra ./
```

&emsp;&emsp;The Aspera transfer platform is the the most advanced software solution for file transfer, synchronization and streaming of digital assets, allowing users and enterprises secure high speed movement of all of their data over any distance, to any environment, with none of the waiting.[中文参考](https://www.jianshu.com/p/a6ac81456c01)

### 辅助工具

##### GNU-parallel

&emsp;&emsp;GNU Parallel是一个shell工具，为了在一台或多台计算机上并行的执行计算任务。通常的输入是文件列表、主机列表、用户列表、URL列表或者表格列表；一个计算任务也可以是一个从管道读取的一条命令。GNU Parallel会把输入分块，然后通过管道并行的执行。GNU Parallel保证它的输出与顺序执行计算任务时是一样的，这样就可以方便的把GNU Parallel的输出做为其它程序的输入。[官方参考](https://www.gnu.org/software/parallel/), [中文参考1](https://blog.csdn.net/huozhanfeng/article/details/38497707),[中文参考2](https://www.jianshu.com/p/c5a2369fa613)

```shell
(wget -O - pi.dk/3 || curl pi.dk/3/) | bash  # installation dir :~/bin
man parallel_tutorial  # man手册
```

### 质量控制

##### fastqc

&emsp;&emsp;FastQC是一款基于Java的软件，一般都是在linux环境下使用命令行运行，它可以快速多线程地对测序数据进行质量评估（Quality Control）。[中文参考](https://zhuanlan.zhihu.com/p/20731723)，[Fastqc官方](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/)

functionality:

- **Import of data from BAM, SAM or FastQ files (any variant)**
- Providing a quick overview to tell you in which areas there may be problems
- Summary graphs and tables to quickly assess your data
- Export of results to an HTML based permanent report
- Offline operation to allow automated generation of reports without running the interactive application

```shell
fastqc [-o output_dir] [--(no)extract] [-f fastq|bam|sam] [-c contaminant file] seqfile1 .. seqfileN

# 主要是包括前面的各种选项和最后面的可以加入N个文件
# -o --outdir FastQC生成的报告文件的储存路径，生成的报告的文件名是根据输入来定的
# --extract 生成的报告默认会打包成1个压缩文件，使用这个参数是让程序不打包
# -t --threads 选择程序运行的线程数，每个线程会占用250MB内存，越多越快咯


# -c --contaminants 污染物选项，输入的是一个文件，格式是Name [Tab] Sequence，里面是可能的污染序列，如果有这个选项，FastQC会在计算时候评估污染的情况，并在统计的时候进行分析，一般用不到
# -a --adapters 也是输入一个文件，文件的格式Name [Tab] Sequence，储存的是测序的adpater序列信息，如果不输入，目前版本的FastQC就按照通用引物来评估序列是否有adapter的残留
# -q --quiet 安静运行模式，一般不选这个选项的时候，程序会实时报告运行的状况。
```

 &emsp;举例说明：

```shell
fastqc -o qc_output -t 5 ExampleData_hg19/Pat1_1.fq
```



##### cutadapter

&emsp;&emsp;cutadapt软件是最常用的去adapter的工具。它是基于Python编写的一个Python包.[中文参考](https://zhuanlan.zhihu.com/p/20776942),[官方文档](https://cutadapt.readthedocs.io/en/stable/guide.html)，主要包括三部分内容：

- read modifications
- cut adapters
- Filtering of processed reads

```shell
# cut adapters
# 去掉3‘端AAAAAAA和5’端的adapter TTTTTTT ,若输入输出压缩文件，则直接output.fq.gz 
# -b表示both,3'和5'端, -j 表示核心数
$ cutadapt -a AAAAAAA -g TTTTTTT -j 4 -o output.fastq input.fastq
#For paired-end reads:
$ cutadapt -a ADAPT1 -A ADAPT2 [options] -o out1.fastq -p out2.fastq in1.fastq in2.fastq

# read modifications
$ cutadapt -u 5 -u -5 -o trimmed.fastq input_reads.fastq

# cutadapt软件可以使用-q参数进行reads质量的过滤。基本原理就是，一般reads头和尾会因为测序仪状态或者是反应时间的问题造成测序质量差，比较粗略的一个过滤办法就是-q进行过滤。需要特别说明的是，这里的-q对应的数字和phred值是不一样的，它是软件根据一定的算法计算出来的。
# 3‘端进行一个简单的过滤,--quality-base=33是指序列使用的是phred33计分系统
$ cutadapt -q 10 --quality-base=33 -o output.fastq input.fastq 
# 3‘端 5’端都进行过滤,3'的阈值是10，5‘的阈值是15
$ cutadapt -q 10,15 --quality-base=33 -o output.fastq input.fastq

# 下面是根据长度过滤
[--minimum-length N or -m N] # 当序列长度小于N的时候，reads扔掉
[--too-short-output FILE] # 上面参数获得的这些序列不是直接扔掉，而是输出到一个文件中
[--maximum-length N or -M N] # 当序列长度大于N的时候，reads扔掉
[--too-long-output FILE] # 上面参数获得的这些序列不是直接扔掉，而是输出到一个文件中
```



##### trim-galore

&emsp;&emsp;<u>Trim Galore是对FastQC和Cutadapt的包装</u>。**适用于所有高通量测序，可以自动检测adapter**，包括RRBS(Reduced Representation Bisulfite-Seq ), Illumina、Nextera 和smallRNA测序平台的双端和单端数据。[推荐一篇有趣的文章](https://www.jianshu.com/p/7a3de6b8e503)，[trim_galore](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)的主要功能包括:

- 去除低质量碱基

- 去除3' 末端的adapter

  Illumina: AGATCGGAAGAGC

  Small RNA: TGGAATTCTCGG

  Nextera: CTGTCTCTTATA

  ```shell
  trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 \
              --paired fq1 fq2  \
              --gzip -o output_dir
   # 质量过滤           
  --quality：设定Phred quality score阈值，默认为20。
  --phred33：：选择-phred33或者-phred64，表示测序平台使用的Phred quality score。
  # 去除adpater
  --adapter：输入adapter序列。也可以不输入，Trim Galore!会自动寻找可能性最高的平台对应的adapter。自动搜选的平台三个，也直接显式输入这三种平台，即--illumina、--nextera和--small_rna。
  --stringency：设定可以忍受的前后adapter重叠的碱基数，默认为1。可以适度放宽，因为后一个adapter几乎不可能被测序仪读到。
  # 长度过滤
  --length：设定输出reads长度阈值，小于设定值会被抛弃。
  # 双端过滤
  --paired：对于双端测序结果，一对reads中，如果有一个被剔除，那么另一个会被同样抛弃，而不管是否达到标准。
  --retain_unpaired：对于双端测序结果，一对reads中，如果一个read达到标准，但是对应的另一个要被抛弃，达到标准的read会被单独保存为一个文件。
  
  --gzip和 --dont_gzip：清洗后的数据zip打包或者不打包。
  ```


##### trimmomatic

&emsp;&emsp;Trimmomatic是**针对Illumina高通量测序平台**设计的**接头去除和低质量reads清洗**软件。软件中包括有Illumina平台常见接头序列，可以很方便处理单端和双端RNA-seq数据。Trimmomatic也支持自己设计要去除的接头序列文件，**目前的HiSeq系列和MiSeq系列用的都是TruSeq3，TruSeq2是以前GA2系列的测序仪所用的**。[Trimmomatic](https://www.jianshu.com/p/a8935adebaae)，[官方文档](http://www.usadellab.org/cms/?page=trimmomatic)

```shell
# PE
java -jar path/to/trimmomatic-0.36.jar PE -phred33 -trimlog logfile input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:/path/to/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36

# SE
java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


# 参数讨论

# ILLUMINACLIP:/path/to/adapters/TruSeq3-PE.fa:2:30:10
TruSeq3-PE.fa是接头序列，2是比对时接头序列时所允许的最大错配数；30指的是要求PE的两条read同时和PE的adapter序列比对，匹配度加起来超30%，那么就认为这对PE的read含有adapter，并在对应的位置需要进行切除【注】。10和前面的30不同，它指的是，我就什么也不管，反正只要这条read的某部分和adpater序列有超过10%的匹配率，那么就代表含有adapter了，需要进行去除；

# SLIDINGWINDOW
滑动窗口长度的参数, SLIDINGWINDOW:5:20代表窗口长度为5，窗口中的平均质量值至少为20，否则会开始切除；
# LEADING
规定read开头的碱基是否要被切除的质量阈值；
#TRAILING，
规定read末尾的碱基是否要被切除的质量阈值；
# MINLEN
规定read被切除后至少需要保留的长度，如果低于该长度，会被丢掉。

# -threads 线程数
# -summary <statsSummaryFile> 信息
# -trimlog <logfile> 日志
# -quiet 安静模式
# -validatePairs PE模式下验证双端
```



##### fastx-toolkit

&emsp;&emsp;fastx Toolkit是包含处理fastq/fasta文件的一系列的工具，它是基于java开发，PGM 测序数据一般用fastx-toolkit，当然其他数据也可以。[中文参考-简](https://zhuanlan.zhihu.com/p/20776942)，[官方文档](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)，包括以下工具：

> FASTQ-to-FASTA---格式转换
> FASTQ/A Quality Statistics---基本统计信息
> FASTQ Quality chart
> FASTQ/A Nucleotide Distribution chart
> FASTQ/A Clipper---reads过滤和adapter裁剪
> FASTQ/A Renamer
> FASTQ/A Trimmer---快速序列切割
> FASTQ/A Collapser
> FASTQ/A Artifacts Filter
> FASTQ Quality Filter
> FASTQ/A Reverse Complement
> FASTA Formatter
> FASTA nucleotides changer
> FASTA Clipping Histogram
> FASTX Barcode Splitter

&emsp;&emsp;&emsp;记录几个主要工具：

```shell
$ fastx_clipper [-a ADAPTER] [-k -v -z] [-i INFILE] [-o OUTFILE]
[-h] = 获得帮助信息.
[-a ADAPTER] = Adapter序列信息. 默认的是CCTTAAGG
[-z] = 调用GZip软件，输出的文件自动经过压缩.

[-k] = 报告adater的序列信息
[-v] = 报告序列总数
[-D]= Debug output.

[-l N] = 如果1条reads小于N就抛弃，默认5.
[-d N] = 保留adapter并保留后面的Nbp，如果设置-d 0等于没有用这个参数.
[-c] = 只保留包含adapter的序列
[-C] = 只保留不包含adapter的序列
[-n] = 如果reads中有N，保留reads.（默认是有N的序列删除）
```

​	

```shell
$ fastx_trimmer [-h] [-f N] [-l N] [-z] [-v] [-i INFILE] [-o OUTFILE]
[-h] = 获得帮助信息.
[-f N] = 序列中从第几个碱基开始保留. 默认是1.
[-l N] = 序列最后保留到多少个碱基，默认是整条序列全部保留.
[-z] = 调用GZip软件，输出的文件自动经过压缩.
```



### 基因组比对

&emsp;&emsp;BWA和Bowtie都是基于BWT转换算法构建的快速比对算法，Bowtie2则是对Bowtie算法的一个改进，如允许deletion的出现。二代测序数据一般长250bp,有相对较高的精度0.1%（Q30），目前比较常用的就是Bowtie2。

##### bowtie

&emsp;&emsp;bowtie1 2009年出现的工具，对于测序长度在50bp以下的序列效果不错，而bowtie2主要针对的是长度在50bp以上的测序的；另外很重要一点，**bowtie不支持 gap open**，**建议一般使用bowtie2。**

##### bwa

&emsp;&emsp;BWA也是一个基因组比对的软件包。它由三种算法组成：BWA-backtrack，BWA-SW和BWA-MEM。第一种算法设计用于Illumina测序序列，最多长100bp；而其余两种用于较长序列，范围为70bp至1Mbp。BWA-MEM和BWA-SW具有相似的功能，例如支持long reads和split alignment，但最新的BWA-MEM通常被推荐用于高质量查询，因为它更快，更准确。对于70-100bp Illumina读数，BWA-MEM还具有比BWA-backtrack更好的性能。[官方文档](http://bio-bwa.sourceforge.net/bwa.shtml)

* build index

  ```bash
  $ bwa index -a bwtsw ref.fa
  
  # 生成以下index文件
  $ ll
  -rw-r--r-- 1 liunianping qukun  3088286508 Sep  2 20:06 hg38.fa.bwt
  -rw-r--r-- 1 liunianping qukun   772071602 Sep  2 20:06 hg38.fa.pac
  -rw-r--r-- 1 liunianping qukun         954 Sep  2 20:06 hg38.fa.ann
  -rw-r--r-- 1 liunianping qukun       16843 Sep  2 20:06 hg38.fa.amb
  -rw-r--r-- 1 liunianping qukun  1544143256 Sep  2 20:20 hg38.fa.sa
  ```

  这里需要注意参数`-a`, 即指定build index的算法，默认是`is`，比对人类基因组、鼠基因组一般都是使用`bwtsw`算法

* mapping

  `mem`算法可用于70bp-1Mbp序列的比对，简单的来讲，`mem` 使用的 MEMs(maximal exact matches） 进行seedling alignments, 再使用 SW(affine-gap Smith-Waterman 算法）进行seeding extending.
  mem 进行局部比对，因此，对于一条序列的不同区域可能会产生多种最优匹配结果， 这对于long reads 来说尤为重要。 有些软件如 Picard’s markDuplicates 跟 mem 的这种剪接性比对不兼容, 在这种情况下，可以使用 –M 选项来将 shorter split hits 标记为次优。。目前bwa团队还更新了比`bwa mem`更先进的方法，参考[github repository](https://github.com/lh3/bwa)基本使用形式为：`bwa mem [options] ref.fa reads.fq [mates.fq]`

  

  > [特别说明]([http://starsyi.github.io/2016/05/24/BWA-%E5%91%BD%E4%BB%A4%E8%AF%A6%E8%A7%A3/](http://starsyi.github.io/2016/05/24/BWA-命令详解/))：
  >
  > * 如果 mates.fq 缺省，且参数 –p 未设定，那么 reads.fq 被认为是 single-end;
  > * 如果 mates.fq 存在，且参数 –p 未设定，那么 `mem` 命令会认为 read.fq 和 mates.fq 中的 i-th reads 组成一个read对 (a read pair)，这个模式是常用的 paired-end mode。
  > * 如果参数 –p 被设定，那么， mem 命令会认为 read.fq 中的 第 2i-th 和 第 (2i + 1)-th 的 reads 组成一个 read 对 （a read pair），这种方式也被成为交错式的（interleaved paired-end)。 在这种情况下，即使有 mates.fq，也会被忽略。

  

  ```bash
  $ bwa mem ref.fa reads.fq > mem-se.sam
  $ bwa mem ref.fa read1.fq read2.fq > mem-pe.sam
  $ bwa mem -t 4 -M -R "\@RG\tID:{library}\tLB:{library}\tPL:Illumina\tPU:{sample}\tSM:{sample}\"  ref.fa read1.fastq read2.fastq > mem-pe.sam 2> ./mem-pe.log
  ```

  

  - `-t INT`  线程数，默认是1。
  - `-M`      将 shorter split hits 标记为次优，以兼容 Picard’s markDuplicates 软件。
  - `-p`      若无此参数：输入文件只有1个，则进行单端比对；若输入文件有2个，则作为paired reads进行比对。若加入此参数：则仅以第1个文件作为输入(输入的文件若有2个，则忽略之)，该文件必须是read1.fq和read2.fa进行reads交叉的数据。
  - `-R STR`  完整的read group的头部，可以用 '\t' 作为分隔符， 在输出的SAM文件中被解释为制表符TAB. read group 的ID，会被添加到输出文件的每一个read的头部。
  - `-T INT`  当比对的分值比 INT 小时，不输出该比对结果，这个参数只影响输出的结果，不影响比对的过程。
  - `-a`      将所有的比对结果都输出，包括 single-end 和 unpaired paired-end的 reads，但是这些比对的结果会被标记为次优。

  

##### bowtie2

&emsp;&emsp;2012年出现的比对工具，用于将测序基因比对到参考基因组，Bowtie 2 支持gapped, local和 paired-end比对，[bowtie2官方文档](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

```shell
# 代码示例
$ bowtie2-build files_str_sepWith_comma ./hg19_index/bt2_hg19_index
$ bowtie2 -p 6 --phred33 --very-sensitive -x ./hg19_index/bt2_hg19_index -U ./sra/SRR363795.fastq -S test_95.sam


# 必须参数
-x <bt2-idx> 由bowtie2-build所生成的索引文件的前缀
-S <hit> 所生成的SAM格式的文件前缀。默认是输入到标准输出。

-1 <m1> 双末端测寻对应的文件1。可以为多个文件，并用逗号分开；多个文件必须和 -2 
<m2> 中制定的文件一一对应。比如:"-1 flyA_1.fq,flyB_1.fq -2 flyA_2.fq,flyB
_2.fq". 测序文件中的reads的长度可以不一样。
-2 <m2> 双末端测寻对应的文件2.
-U <r> 非双末端测序对应的文件。可以为多个文件，并用逗号分开。测序文件中的reads的
长度可以不一样。

# 可选参数
# 输入形式选项
-q # 输入的文件为FASTQ格式文件，此项为默认值。
-qseq # 输入的文件为QSEQ格式文件。
-f # 输入的文件为FASTA格式文件,此时 --ignore-quals被默认加上。
-r # 输入的文件中，每一行代表一条序列，没有序列名和测序质量等,此时 --ignore-quals被默认加上。
-c <seq_str_sepWith_comma> # 紧跟其后为比对的reads序列，序列间用逗号隔开,此时 --ignore-quals被默认加上。
-s/--skip <int> # input的reads中，跳过前<int>个reads或者pairs。
-u/--qupto <int> # 只比对前<int>个reads或者pairs
-5/--trim5 <int> # 剪掉5'端<int>长度的碱基，再用于比对。(default: 0).
-3/--trim3 <int> # 剪掉3'端<int>长度的碱基，再用于比对。(default: 0).
--phred33 # 输入的碱基质量等于ASCII码值加上33.
--phred64 # 输入的碱基质量等于ASCII码值加上64.
--solexa-quals # 将Solexa的碱基质量转换为Phred。在老的GA Pipeline版本中得以
运用,Default: off.
--int-quals # 输入文件中的碱基质量为用“ ”分隔的数值，而不是ASCII码,Default: off.
–local # local alignment
```

&emsp;&emsp;举一个实例，来源[官方文档](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

1. 比对单端测序数据

   - 建立索引

     ```shell
     $ BT2_HOME='path/to/bowtie2-2.3.4.3/'
     $ BT2_HOME/bowtie2-build $BT2_HOME/example/reference/lambda_virus.fa lambda_virus
     ```

     ![](./imgs/bt2_buid.png)

   - 比对

     ```shell
     # 默认全局比对，可以加上--local，进行局部比对
     $ bowtie2 -x lambda_virus -U $BT2_HOME/example/reads/reads_1.fq -S eg1.sam
     ```

     ![](./imgs/bt2-align.png)

2. 比对双端测序数据

   - 比对

     ```shell
     $BT2_HOME/bowtie2 -x $BT2_HOME/example/index/lambda_virus -1 $BT2_HOME/example/reads/reads_1.fq -2 $BT2_HOME/example/reads/reads_2.fq -S eg2.sam
     ```

   - SAM-to-BAM

     ```shell
     $ samtools view -bS eg2.sam > eg2.bam
     ```

   - BAM-to-sortedBAM

     ```shell
     $ samtools sort eg2.bam -o eg2.sorted.bam
     ```

     &emsp;&emsp;Sorted BAM is a useful format because the alignments are (a) compressed, which is convenient for long-term storage, and (b) sorted, which is conveneint for variant discovery. To generate variant calls in VCF format.

   - downstream analysis

     ```shell
     $ samtools mpileup -uf $BT2_HOME/example/reference/lambda_virus.fa eg2.sorted.bam | bcftools view -Ov - > eg2.raw.bcf
     
     $ bcftools view eg2.raw.bcf
     ```

##### STAR

&emsp;&emsp;ENCODE计划御用比对软件，也是比较权威，并且比对速度极快，一般用于RNA-Seq测序数据的比对。详细请参考[官方文档](https://github.com/alexdobin/STAR/)，或则[中文参考](http://www.bioinfo-scrounger.com/archives/288)以下为核心内容：

- 安装

  ```shell
  # 最新版本请参考官方文档
  wget https://github.com/alexdobin/STAR/archive/2.7.0c.tar.gz
  tar -xzf 2.7.0c.tar.gz
  cd STAR-2.7.0c
  ```

- 建立索引

  ```shell
  # 参数如下
  --runThreadN NumberOfThreads
  --runMode genomeGenerate
  --genomeDir /path/to/genomeDir
  --genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 ... 
  --sjdbGTFfile /path/to/annotations.gtf
  --sjdbOverhang ReadLength-1 
  ```

- 比对到基因组

  ```shell
  # 参数如下
  --runThreadN NumberOfThreads
  --genomeDir /path/to/genomeDir
  --readFilesIn /path/to/read1 [/path/to/read2]
  --readFilesCommand UncompressionCommand
  --outFileNamePrefix /path/to/output/dir/prefix
  --outSAMtype BAM SortedByCoordinate  # 输出默认排序的bam文件
  --genomeLoad LoadAndKeep  # 共享内存
  
  # 需要注意的是：
  
  # for single-end
   --readFilesIn sample1.fq,sample2.fq,sample3.fq
  # for pair-end
   --readFilesIn sample1read1.fq,sample2read1.fq sample1read2.fq,sample2read2.fq 
   
  # outSamType
  --outSAMtype BAM Unsorted
  --outSAMtype BAM Unsorted SortedByCoordinate
  ```

- 输出文件

  - Log Files

    包括 Log.out， Log.progress.out，Log.final.out

  - SAM

  - Unsorted and sorted-by-coordinate BAM

  - Splice junctions

- two-pass比对

  > &emsp;&emsp;For the most sensitive novel junction discovery,I would recommend running STAR in the 2-pass mode. **It does not increase the number of detected novel junctions, but allows to detect more splices reads mapping to novel junctions.** The basic idea is **to run 1st pass of STAR mapping with the usual parameters, then collect the junctions detected in the  rst pass, and use them as "annotated" junctions for the 2nd pass mapping.** 



##### Kallisto

> 论文：Near-optimal probabilistic RNA-seq quantification（http://dx.doi.org/10.1038/nbt.3519），它的优势是快速、耗费内存小，可以在普通的台式机上几分钟之内完成人或其他物种的转录组二代测序比对任务。其最大特点是不需要完整的参考基因组，只需要对某物种的全转录本序列建立索引，然后再进行pseudo alignment.他是一款 RNA-seq 数据量化软件，即从RNA-seq 的测序数据中计算出每一个基因的表达量。传统的 RNA-seq 数据分析思路分两步，第一步是把利用 RNA-seq 方法得到的测序数据先比对到参考基因组序列上（tophat2, bowtie2, HISAT, STAR等软件）；第二步是从比对结果中计算表达量，可以理解为数每一个基因的 reads 数量（Cufflinks，HTseq-count 等软件）。



### 峰值探测

##### MACS2

&emsp;&emsp;MACS2是peak calling最常用的工具，这是MACS2的主要功能，因为MACS2的目的就是找peak，其他功能都是可有可无，唯独`callpeak`不可取代。[官方文档](<https://github.com/taoliu/MACS/>)

基本用法：

```shell
# Example for regular peak calling:
macs2 callpeak -t ChIP.bam -c Control.bam -f BAM -g hs -n test -B -q 0.01
# Example for broad peak calling:
macs2 callpeak -t ChIP.bam -c Control.bam --broad -g hs --broad-cutoff 0.1
```

参数如下：

![](./imgs/macs2.png)

- `-t`  实验组
- `-c` 对照组
- `-f` 文件类型，支持BAM SAM BED等，常见文件类型可以自动识别，除了"BAMPE"和"BEDPE"
- `-g` 基因组的大小，人类基因组使用'hs'【The default hs -- 2.7e9 is recommended for UCSC human hg18 assembly】, 鼠基因组使用'mm'
- `--outdir` 输出文件的文件路径
- `-n`  即name，取一个有区分性有意义的名字，用于作为结果输出文件的前缀
- `-B` 指定输出bedgraph格式的文件
- `--verbose`  verbose模式，如果是0表示不输出过程信息
- `-s/--tsize` 测序的读长，如果不指定，将会自动推测 
- `--qvalue / --pvalue` q value即minimum FDR，默认cutoff值是0.05，如果指定了pvalue, 那么qvalue将会被忽略。
- `--broad / --broad-cutoff` peak有narrow peak和broad peak之分，cutoff和qvalue相似 
- `--nomodel` 这个参数说明不需要MACS去构建模型，双端测序不需要这个值。
- `--extsize` 当设置了`--nomodel`时，MACS会用`--extsize`这个参数从5'->3'方向扩展reads修复fragments。比如说你的转录因子结合范围200bp，就设置这个参数是200。
- `--shift`  当设置了`--nomodel`，MACS用这个参数从5' 端移动剪切，然后用`--extsize`延伸，如果`--shift`是负值表示从3'端方向移动。建议ChIP-seq数据集这个值保持默认值为0，对于检测富集剪切位点如DNAsel数据集设置为EXTSIZE的一半。

输出结果解读：

&emsp;假设我们对ATAC-seq数据进行 callpeak：`macs2 callpeak -t norm.pe.q10.sort.rmdup.shift.bed -f BED -B -g hs --outdir output -q 0.01 -n myatac --nomodel --shift 0`，结果将会得到以下文件：

- `myatac_peaks.xls`

  ![](./imgs/macs0001.png)

- `myatac_peaks.narrowPeak`

  前四行不言而喻，分别是 染色体、起始位置、结束位置、peak名称

  - 5th: integer score for display calculated as `int(-10*log10qvalue)`. Please note that currently this value might be out of the [0-1000] range defined in [UCSC Encode narrowPeak format](https://genome.ucsc.edu/FAQ/FAQformat.html#format12)
  - 7th: fold-change
  - 8th: -log10pvalue
  - 9th: -log10qvalue
  - 10th: relative summit position to peak start

- `myatac_summits.bed`

  前四行同上，第五行是-log10pvalue

  > If you want to find the motifs at the binding sites, this file is recommended

- bedgraph文件—> `myatac_control_lambda.bdg`和 `myatac_treat_pileup.bdg`

  这两个文件可以直接导入UCSC Browser
  
 > 更多使用疑问，可以参考[google 讨论组](https://groups.google.com/g/macs-announcement),如[pileup](https://groups.google.com/g/macs-announcement/c/c__KCotHsok)

### 差异分析

##### Deseq2

&emsp;差异基因表达分析常用的R包有DESeq2，edgeR，limma等几种，目前最常用的就是DESeq2

- 使用方式1-一步执行

  1. 构造dds的对象

     -  准备数据

     ```R
     > library("pasilla")
     > pasCts <- system.file("extdata",
     +                       "pasilla_gene_counts.tsv",
     +                       package="pasilla", mustWork=TRUE)
     > pasAnno <- system.file("extdata",
     +                        "pasilla_sample_annotation.csv",
     +                        package="pasilla", mustWork=TRUE)
     > cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
     > coldata <- read.csv(pasAnno, row.names=1)
     > data.frame(coldata[,c('condition')],row.names = row.names(coldata)) -> coldata
     > names(coldata)<-c('group_info')
     > coldata
                  group_info
     treated1fb      treated
     treated2fb      treated
     treated3fb      treated
     untreated1fb  untreated
     untreated2fb  untreated
     untreated3fb  untreated
     untreated4fb  untreated
     > head(cts,1)
                 untreated1 untreated2 untreated3 untreated4 treated1 treated2 treated3
     FBgn0000003          0          0          0          0        0        0        1
     > 
     > # 处理cts和coldata一致性问题
     > rownames(coldata) <- sub("fb", "", rownames(coldata))
     > cts <- cts[, rownames(coldata)]
     > coldata
                group_info
     treated1      treated
     treated2      treated
     treated3      treated
     untreated1  untreated
     untreated2  untreated
     untreated3  untreated
     untreated4  untreated
     > head(cts,1)
                 treated1 treated2 treated3 untreated1 untreated2 untreated3 untreated4
     FBgn0000003        0        0        1          0          0          0          0
     > all(rownames(coldata) == colnames(cts))
     [1] TRUE
     ```

     - 构造dds对象

     ```R
     > suppressMessages(library("DESeq2"))
     > dds <- DESeqDataSetFromMatrix(countData = cts,
     +                               colData = coldata,
     +                               design = ~ group_info)
     > dds
     ```

     

  2. 用`DESeq()`函数进行normalization处理

      

     ```R
     > dds <- DESeq(dds)
     using pre-existing size factors
     estimating dispersions
     found already estimated dispersions, replacing these
     gene-wise dispersion estimates
     mean-dispersion relationship
     final dispersion estimates
     fitting model and testing
     ```

     

  3. 用`results()`函数来提取差异比较结果

     ```R
     > resultsNames(dds)
     [1] "Intercept"                       "group_info_untreated_vs_treated"
     > res<-results(dds, name='group_info_untreated_vs_treated')
     > head(res,3)
     > summary(res)
     
     # 按pvalue排序
     > resOrdered <- res[order(res$pvalue),]
     > head(resOrdered)
     
     # 可视化结果
     > plotMA(res, ylim=c(-2,2))
     ```

- 使用方式2-分步执行

  首先准备数据部分以及构建dds对象同上，我就不重复这一步工作了

  1. estimation of size factors(`dds <- estimateSizeFactors(dds)`)

      

     ```R
     > dds <- estimateSizeFactors(dds)
     > count_norm=counts(dds, normalized=TRUE)  # 获取normalized之后的数据
     > count_raw=counts(dds, normalized=F)  # 获取原始数据
     > write.table(count_norm,file='a.txt',row.names=T,col.names=T,quote=F,sep='\t')
     
     > count_mat_norm=as.matrix(count_norm)
     > log2_count_mat_norm=log2(count_mat_norm+1) 
     > dimnames(log2_count_mat_norm)=list(row.names(counts),colnames(counts))
     > write.table(log2_count_mat_norm,file='b.txt',row.names=T,
                   col.names=T,quote=F,sep='\t')
     ```

     

  2. estimation of spersions(`dds <- estimateDispersions(dds)`)

       

     ```R
     > dds=estimateDispersions(dds,fitType="local")
     gene-wise dispersion estimates
     mean-dispersion relationship
     final dispersion estimates
     ```

     

  3. Negative Binomial GLM fitting and Wald statistics(`dds <- nbinomWaldTest(dds)`)

  ```R
  > dds <- nbinomWaldTest(dds)
  > resultsNames(dds)
  [1] "Intercept"                       "group_info_untreated_vs_treated"
  ```

&emsp;剩下的步骤也和上面方式1相同~

### motif分析

##### Homer

&emsp;&emsp;HOMER最初是为了使发现ChIP-Seq peaks富集motif的过程自动化，更一般地，HOMER分析富集motif的基因组位置，不仅限于ChIP-Seq peaks。使用homer，所有用户真正需要的是包含基因组坐标的文件（HOMER peak file or  BED file），然后HOMER会处理后续过程。[homer motif analysis manual](<http://homer.ucsd.edu/homer/ngs/peakMotifs.html>)



- 安装

  [官方安装参考](<http://homer.ucsd.edu/homer/introduction/install.html>)

  ```shell
  cd /path/to/install_homer
  wget http://homer.ucsd.edu/homer/configureHomer.pl
  perl /Users/chucknorris/homer/configureHomer.pl -install
  
  # 安装你需要的参考基因组，这一步取决于网速
  perl /path-to-homer/configureHomer.pl -install hg19
  perl /path-to-homer/configureHomer.pl -install hg38
  ```

- 基本使用

  ```shell
  # By default this will perform de novo motif discovery as well as check the enrichment of known motifs.
  findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]
  findMotifsGenome.pl ERpeaks.txt hg19 Output/motif/ -size 200 -mask
  ```

- 参数解析

  -  `-mask`  : optional and tells the program to use the repeat-masked sequence
  -  `-size`  : 必须参数, 默认是200
  -  `-bg <peak/BED file>` ： 设置背景数据
  -  `-noweight / -nlen 0` : 关闭自动对数据的normalization

  > `-size`参数是最重要的参数之一，也是许多人容易混淆的地方。如果希望根据你的Peak本身宽度来搜寻motif，请使用`-size given`，但是对于转录因子的peak，motif一般位于peak中心+/- 50-75的位置，所以用一个确切的 `-size`也许会更好。


- 输入文件`<peak/BED file>`

  如果输入文件是peak文件，至少需要5列，分别如下：

  - Column1: Unique Peak ID
  - Column2: chromosome
  - Column3: starting position
  - Column4: ending position
  - Column5: Strand (+/- or 0/1, where 0="+", 1="-")

  如果输入是BED文件，至少需要六列，分别如下：

  - Column1: chromosome
  - Column2: starting position
  - Column3: ending position
  - Column4: Unique Peak ID
  - Column5: not used
  - Column6: Strand (+/- or 0/1, where 0="+", 1="-")

  当然，如果你的输入文件多余上述的规定，也不要紧，程序会自动忽略。

- 输出文件
  - knownResults.html : formatted output of known motif finding.
  - homerResults.html : formatted output of *de novo* motif finding.
  - homerMotifs.motifs : these are the output files from the de novo motif finding, separated by motif length, and represent separate runs of the algorithm.
  - homerMotifs.all.motifs : Simply the concatenated file composed of all the homerMotifs.motifs<#> files
  - motifFindingParameters.txt : command used to execute findMotifsGenome.pl
  - knownResults.txt : text file containing statistics about known motif enrichment (open in EXCEL).
  - seq.autonorm.tsv : autonormalization statistics for lower-order oligo normalization.

- 程序执行步骤

  Homer通过以下步骤来[搜寻高质量的motif](<http://homer.ucsd.edu/homer/ngs/peakMotifs.html>)：

  1. Verify peak/BED file
  2. Extract sequences from the genome corresponding to the regions in the input file, filtering sequences that are >70% "N"
  3. Calculate GC/CpG content of peak sequences
  4. Preparse the genomic sequences of the selected size to serve as background sequences.
  5. Randomly select background regions for motif discovery
  6. Autonormalization of sequence bias
  7. Check enrichment of known motifs
  8. *de novo* motif finding

Homer这个程序集中还包括很多其他的脚本，这里再添加几个脚本的用法：

- `scanMotifGenomeWide.pl`

  > To make it easier to predict motif sites across the genome, HOMER contains a program called **scanMotifGenomeWide.pl** to assist with scanning large FASTA files.

  * 基本用法

    ```shell
    scanMotifGenomeWide.pl <motif file> <genome> [options] # 默认输出到stdout
    scanMotifGenomeWide.pl ctcf.motif hg19 -bed > pu1.sites.mm9.bed
    ```

  * 其他参数

    * **-bed** : 指定输出格式为bed格式
    * **-keepAll** :默认情况下输出结果不包括自身，此选项会保留
    * **-mask** :不在RepeatMasked sequence (lower case sequence in FASTA files)中搜寻motif
    * **-5p** : 基于motif 序列的5'端确定其位置
    * **-homer1** or **-homer2** : 指定使用哪个版本 (-homer2 is default)

  * 输出文件

    * Tab delimited text file (default):
      1. Site ID (motif name + number)
      2. chr
      3. start
      4. end
      5. strand
      6. log-odds score
      7. sequence

    * BED (tab) format (use **-bed**):
      1. chr
      2. start
      3. end
      4. motif name
      5. log-odds score (will be floored to an integer)
      6. strand

### 可视化

##### deeptools*

&emsp;&emsp;deeptools就是用来对单个或者多个比对好的bam文件进行信息统计并进行可视化分析的，所以ChIP-seq和RNA-seq及其它类型的二代测序结果都是可以借以分析。[deeptools-简书](https://www.jianshu.com/p/7cc5df9f7900)

##### igv

&emsp;[基因组数据可视化工具](<https://software.broadinstitute.org/software/igv/home>)(Integrative Genomics Viewer Tools)，其可以对`sorted BAM，BED，GTF/GFF，bigwig/wig，bedgraph`等数据进行可视化，是最常用可视化工具之一。无论是全基因组测序，RNA-Seq，ChIP-Seq等等都可以进行可视化。基本操作如下：(导入时同时需要sorted bam和bam index文件，具体操作请见[File_Format](FileFormat.md))

1. 最上面的Human hg19 是选择参考基因组信息；
2. 中间的空白部分是载入数据以后的可视化部分；
3. 中间的空白栏是一个搜索栏，可以搜索一个基因或具体的坐标比如 chr1:10000-12000 ；
4. 最下方的注释信息为RefSeq Genes；

