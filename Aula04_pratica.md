# Bioinformática Instrumental - Módulo II
# Identificação de Variantes Genéticas

## Aula04: Pré-processamento de dados genômicos  

Essa é uma etapa obrigatória que precede toda a descoberta de variantes. 
Envolve o pré-processamento dos dados brutos da sequência (geralmente obtidos no formato FASTQ) para produzir arquivos BAM prontos para análise. O principal passo é o alinhamento das sequências de leituras com um genoma de referência, bem como algumas operações de "limpeza" dos dados para corrigir vieses técnicos e torná-los mais adequados para chamada de variantes.

Utilizaremos como referência para a aula as boas práticas de pré-processamento preconizadas pelo [GATK](https://software.broadinstitute.org/gatk/) (_Genome Analysis Toolkit_):

![referências bibliográficas]()

*The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data.* McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA, 2010 GENOME RESEARCH 20:1297-303
 
*From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline.* Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M, 2013 CURRENT PROTOCOLS IN BIOINFORMATICS 43:11.10.1-11.10.33  


**Fluxograma de trabalho desta aula prática:**  

![preprocessing](https://github.com/cmasotti/BioinfoInstrumental-Aula04/blob/master/bestPractices_preprocessing_smaller.jpg)


A seguir, executaremos o passo-a-passo do pré-processamento dos dados:

### PASSO 1: ACESSO AO SERVIDOR REMOTO BIOINFO
Para relembrar como fazer a conexão ao servidor remoto via Putty, reveja os passos [aqui](https://github.com/cmasotti/BioinfoInstrumental-Aula04/blob/master/Acesso_servidor_remoto.pdf).


```bash  
$test  
``` 


![preprocessing]()  



