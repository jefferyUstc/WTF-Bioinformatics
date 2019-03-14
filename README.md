## 目录

- [Tools](Tools)
- [How_To_Read_Paper](ReadPaper)
- [File_Format](FileFormat)
- [Common-Sense](#Common-Sense)
  * [bio-database](#bio-database)
  * [sequence-representation](#sequence-representation)
  * [suggestion](#suggestion)

# Common-Sense

### bio-database

- UCSC
- ENSEMBL
- NCBI

### sequence-representation

- Amino-Acid

  > ```text
  > A  alanine               P  proline       
  > B  aspartate/asparagine  Q  glutamine      
  > C  cystine               R  arginine      
  > D  aspartate             S  serine      
  > E  glutamate             T  threonine      
  > F  phenylalanine         U  selenocysteine      
  > G  glycine               V  valine        
  > H  histidine             W  tryptophan        
  > I  isoleucine            Y  tyrosine
  > K  lysine                Z  glutamate/glutamine
  > L  leucine               X  any
  > M  methionine            *  translation stop
  > N  asparagine            -  gap of indeterminate length
  > ```

  &emsp;

- Nucleotide

  > ```text
  > A  adenosine          C  cytidine             G  guanine
  > T  thymidine          N  A/G/C/T (any)        U  uridine 
  > K  G/T (keto)         S  G/C (strong)         Y  T/C (pyrimidine) 
  > M  A/C (amino)        W  A/T (weak)           R  G/A (purine)        
  > B  G/T/C              D  G/A/T                H  A/C/T      
  > V  G/C/A              -  gap of indeterminate length
  > ```

&emsp;

### suggestion

1. 建议分析测序数据之前先搞清楚以下两个方面：
   - 原始数据是通过哪种测序平台产生的；它们的错误率分布是怎么样的；是否有一定的偏向性和局限性；是否会显著受GC含量的影响等；
   - 评估它们有可能影响哪些方面的分析；
2. 从以下方面认识原始测序数据：
   - read各个位置的碱基质量值分布
   - 碱基的总体质量值分布
   - read各个位置上碱基分布比例，目的是为了分析碱基的分离程度
   - GC含量分布
   - read各位置的N含量
   - read是否还包含测序的接头序列
   - read重复率，这个是实验的扩增过程所引入的