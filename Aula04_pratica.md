# Bioinformática Instrumental - Módulo II
# Identificação de Variantes Genéticas

## Aula04: Pré-processamento de dados genômicos  

Essa é uma etapa obrigatória que precede toda a descoberta de variantes. 
Envolve o pré-processamento dos dados brutos da sequência (geralmente obtidos no formato FASTQ) para produzir arquivos BAM prontos para análise. O principal passo é o alinhamento das sequências de leituras com um genoma de referência, bem como algumas operações de "limpeza" dos dados para corrigir vieses técnicos e torná-los mais adequados para chamada de variantes.

Utilizaremos como referência para a aula as boas práticas de pré-processamento preconizadas pelo [GATK](https://software.broadinstitute.org/gatk/) (_Genome Analysis Toolkit_):


>*The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.* McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303
 
>*From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline.* Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M, 2013 CURRENT PROTOCOLS IN BIOINFORMATICS 43:11.10.1-11.10.33  


**Fluxograma de trabalho desta aula prática:**  

![preprocessing](https://github.com/cmasotti/BioinfoInstrumental-Aula04/blob/master/bestPractices_preprocessing_smaller.jpg)


A seguir, executaremos o passo-a-passo do pré-processamento dos dados:

### PASSO 1: ACESSO AO SERVIDOR REMOTO BIOINFO
Para relembrar como fazer a conexão ao servidor remoto via Putty, reveja os passos [aqui](https://github.com/cmasotti/BioinfoInstrumental-Aula04/blob/master/Acesso_servidor_remoto.pdf).

### PASSO 2: ORGANIZAR ESTRUTURA DE DIRETÓRIOS NO SERVIDOR REMOTO

No prompt da linha de comando no servidor remoto, criar diretórios para pré-processamento.

```bash  
aluno30@5b6864eb3f67:~$
aluno30@5b6864eb3f67:~$ mkdir preprocessing
aluno30@5b6864eb3f67:~$ cd preprocessing
aluno30@5b6864eb3f67:~/preprocessing$ mkdir hg38
aluno30@5b6864eb3f67:~/preprocessing$ mkdir mapping
aluno30@5b6864eb3f67:~/preprocessing$ mkdir alignment_metrics
aluno30@5b6864eb3f67:~/preprocessing$ mkdir markduplicates
aluno30@5b6864eb3f67:~/preprocessing$ mkdir references
aluno30@5b6864eb3f67:~/preprocessing$ mkdir bqsr
```  

Confira os diretórios criados:
```bash  
aluno30@5b6864eb3f67:~/preprocessing$ ls
alignment_metrics  bqsr  hg38  mapping  markduplicates  references
```  

### PASSO 3: PREPARAR O GENOMA DE REFERÊNCIA (Hg38)
 - Genoma de referência é um arquivo muito grande - trabalharemos com o genoma já pronto para a análise. 
 - Peculiaridades do hg38 em relação ao hg19 [HumanGenomeReferenceBuilds](https://gatkforums.broadinstitute.org/gatk/discussion/11010/human-genome-reference-builds-grch38-hg38-b37-hg19)
 - hg38.fa pode ser baixado a partir do repositório [GATK-Bundle](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/). 
 - Uma vez disponível, é preciso indexar o hg38.fa de acordo com as exigências do mapper/caller que utilizaremos. O BWA requer um conjunto diferente de arquivos de índice para alinhamento e o comando abaixo cria cinco dos seis arquivos necessários (**.fai, .pac, .bwt, .ann, .amb, .sa**).
  
```bash   
aluno30@5b6864eb3f67:~/preprocessing$ bwa index

Usage:   bwa index [options] <in.fasta>

Options: -a STR    BWT construction algorithm: bwtsw, is or rb2 [auto]
         -p STR    prefix of the index [same as fasta name]
         -b INT    block size for the bwtsw algorithm (effective with -a bwtsw) [10000000]
         -6        index files named as <in.fasta>.64.* instead of <in.fasta>.* 

Warning: `-a bwtsw' does not work for short genomes, while `-a is' and
         `-a div' do not work not for long genomes.
```  
 - Confira os índices para análise:
```bash   
aluno30@5b6864eb3f67:~/preprocessing$ ls -ltr /mnt/dados/aula4/hg38/
```  
 >**hg38.25chrs.fa** (download GATK-budle)  
 >hg38.fa.fai, hg38.fa.amb, hg38.fa.ann, hg38.fa.bwt, hg38.fa.pac, hg38.fa.sa   

 - É preciso criar o arquivo hg38.dict,um dicionário das sequências de referência FASTA (Picard [CreateSequenceDicttionary](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.3.0/picard_sam_CreateSequenceDictionary.php))  

**Faça um link simbólico para o genoma referência e seus índices na sua pasta "preprocessing" e confira:**
```bash   
aluno30@5b6864eb3f67:~/preprocessing/hg38$ ln -s /mnt/dados/aula4/hg38/* .  
aluno30@5b6864eb3f67:~/preprocessing/hg38$ ls
```  

### PASSO 4: SALVAR E CONFERIR ARQUIVOS DE DADOS BRUTOS DE SEQUENCIAMENTO
Criar link simbólico para os arquivos FASTQs de dados públicos de câncer no diretório **/mapping**.
```bash   
aluno30@5b6864eb3f67:~/preprocessing$ cd mapping
aluno30@5b6864eb3f67:~/preprocessing/mapping$ ln -s /mnt/dados/aula4/raw/*fastq .
aluno30@5b6864eb3f67:~/preprocessing/mapping$ ls
```   

Explore os arquivos fastqs com head/tail:
```bash   
aluno30@5b6864eb3f67:~/preprocessing/mapping$ head TCGA-BH-A1F0-01A_BRCA_R1.fastq 
@HWI-ST467_110093165:5:2:11191:90620/1
CAGGCACAGCCAAGAGGGCTGAAGAAATGGTAGAACGGAGCAGCTGGTGATGTGTGGGCCCACCGGCCCCAGGCTCCTGTCTCCCCCCAGGGGGGTGGTG
+
9;;9?CCBEDDEDEBECCDEDBEECDEDDDBAECEE>DCEDDFEEDD?DABDAD?BAD?ADED<6B;7B,ADA:AB?0D@BDD?CAAD?E##########
```   

### PASSO 5: MAPEAMENTO

Para o correto maepamento, é preciso atribuir corretamente os nomes dos reads, ou "read groups".
>__Por que atribuir corretamente os read groups (RG)?__  
>Para diferenciar não apenas amostras, mas também características técnicas de artefatos. Com essas informações em mãos, podemos mitigar os efeitos desses artefatos durante as etapas de marcação de reads duplicados (PASSO X) e BQSR (PASSO XX).[RG_required_by_GATK](https://software.broadinstitute.org/gatk/documentation/article?id=6472)  

Buscamos as informações das amostras aqui analisadas no repositório de dados do TCGA:  

>BREAST CANCER SAMPLE (TCGA-BH-A1F0)  
> - **WXS primary tumor** [TCGA-BH-A1F0-01A-11D-A135-09](https://portal.gdc.cancer.gov/files/68ada300-f0a2-447a-aa47-865770a80125)  
>  - **WXS adjacent normal tissue** [TCGA-BH-A1F0-01A-11D-A135-09](https://portal.gdc.cancer.gov/files/68ada300-f0a2-447a-aa47-865770a80125)  

![info @RG](https://github.com/cmasotti/BioinfoInstrumental-Aula04/blob/master/RG.png)  

**Executar a linha de comando a seguir (piped command line) para as duas amostras TCGA:**   
>TCGA-BH-A1F0-01A (WXS primary tumor)
```bash   
aluno30@5b6864eb3f67:~/preprocessing/mapping$ bwa mem -M -t4 -R '@RG\tID:74ed7812-25ef-40ff-aca8-dea5ccb39851\tSM:TCGA-BH-A1F0-01A\tPL:ILLUMINA\t' ../hg38/hg38.fa TCGA-BH-A1F0-01A_BRCA_R1.fastq TCGA-BH-A1F0-01A_BRCA_R2.fastq | samtools view -@4 -Sb - -O BAM -o TCGA-BH-A1F0-01A_BRCA.bam   
```  
>TCGA-BH-A1F0-11B (WXS normal tissue)
```bash   
aluno30@5b6864eb3f67:~/preprocessing/mapping$ bwa mem -M -t4 -R '@RG\tID:3ac135b5-f024-4534-a513-7adb9f04cc00\tSM:TCGA-BH-A1F0-11B\tPL:ILLUMINA\t' ../hg38/hg38.fa TCGA-BH-A1F0-11B_BRCA_R1.fastq TCGA-BH-A1F0-11B_BRCA_R2.fastq | samtools view -@4 -Sb - -O BAM -o TCGA-BH-A1F0-11B_BRCA.bam   
```   

Observe o que a segunda parte da linha de comando faz. 
Por que convertemos .sam para .bam?

Para prosseguir, também precisamos ordenar:
```bash   
aluno30@5b6864eb3f67:~/preprocessing/mapping$ samtools sort -@4 TCGA-BH-A1F0-01A_BRCA.bam -O BAM -o TCGA-BH-A1F0-01A_BRCA_sorted.bam  
aluno30@5b6864eb3f67:~/preprocessing/mapping$ samtools sort -@4 TCGA-BH-A1F0-11B_BRCA.bam -O BAM -o TCGA-BH-A1F0-11B_BRCA_sorted.bam  
```  
... e indexar nossos dados mapeados:
```bash   
aluno30@5b6864eb3f67:~/preprocessing/mapping$ samtools index -@2 TCGA-BH-A1F0-01A_BRCA_sorted.bam  
aluno30@5b6864eb3f67:~/preprocessing/mapping$ samtools index -@2 TCGA-BH-A1F0-11B_BRCA_sorted.bam  
```   

### PASSO 6: ANÁLISE DA QUALIDADE DO MAPEAMENTO  
Para avaliar a qualidade do mapeamento dos reads, executamos um script do pacote de ferramentas Picard

