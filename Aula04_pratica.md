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

Preparando o **genoma de referência (hg38)**:
 - Genoma de referência é um arquivo muito grande - trabalharemos com o genoma já pronto para a análise. 
 - Peculiaridades do hg38 em relação ao hg19 [HumanGenomeReferenceBuilds](https://gatkforums.broadinstitute.org/gatk/discussion/11010/human-genome-reference-builds-grch38-hg38-b37-hg19)
 - hg38.fa pode ser baixado a partir do repositório [GATK-Bundle]([ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/). 
  - Uma vez disponível, é preciso indexar o hg38.fa de acordo com as exigências do mapper/caller que utilizaremos (BWA-MEM/GATK, vide abaixo).
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



![preprocessing]()  



