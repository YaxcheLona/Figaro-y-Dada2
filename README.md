
**Figaro_Y_Dada2**



##



Figaro es una herramienta bioinformatica que nos permite optimizar los parametros de corte de mirobioma rRNA. Se encarga de obtener las tasas de error de los archivos FASTQ para utilizarlos en pipelines de alta resolución de Dada2 y Deblur.


En este caso, los archivos que se analizarán, deben de cumplir con ciertos requisitos, por ejemplo:

-Debe de haber un directorio espécifico con las secuencias FASTQ

-Los archivos forward y reverse tienen que estar juntos

-las lecturas deben de provenir de la misma secuenciación y del mismo metodo de preparación de librerías

**Instalación**


+Se usará la instalación por línea de comandos y se clonará FIGARO en tu repositorio

```
git clone https://github.com/Zymo-Research/figaro.git

```
+Desciende en el directorio de FIGARO

```
cd figaro

```
+Instala los paquetes de Python, porque este programa requiere te una interpretación python3

```
pip3 install -r requirements.txt

```

**Correr el programa**

Figaro analizará analizará el directorio con todos los archivos FASTQ con el siguiente comando:

```
python3 figaro.py -i /path/to/fastq/directory -o /path/to/output/files -a 450 -f 20 -r 20

```
**Dada2**

Dada2 tambien puede usarse como paqueteria de R, donde se puede usar Rstudios (server) para su uso
a traves del servidor en la siguiente liga:http://132.248.15.30:4607/

**1. Descarga**

- [ ] Descargar Dada2 directo de github (en caso de existir una versión reciente cambiarla)

```bash
install.packages("devtools")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16")```

**2. Carga y lectura de datos**

- [ ] Cargamos Dada2 y observamos la version 

```bash
library(dada2)
packageVersion("dada2")
```

- [ ] Cargar las librerias 

```bash
path <- "Direccion/en/tu/directorio"
```

- [ ] Observamos el directorio cargado 

```bash
list.files(path)
```

- [ ] Transformamos los archivos en objetos los cuales nos permiten separar en Forward y Reverse 

```bash
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
```

- [ ] Simplificamos los nombres de los archivos

```bash
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
sample.names <- sapply(strsplit(basename(fnRs), "_R"), `[`, 1)
```
- [ ] Observamos los archivos y nombres de muestras

```bash
cat("Los archivos de las secuencias son:\n")
print(fnFs)
cat("\nLos nombres de las muestras son:\n")
print(sample.names)
```
- [ ] Graficamos la calidad de los archivos en R1 y R2

```bash
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```
**3. Filtrado**

- [ ] Asignar nombres a los archivos para los fastq filtrados 

```bash
library(dada2)
packageVersion("dada2")
```

- [ ] Cargar las librerias 

```bash
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

- [ ] Filrado
En esta parte se utilizaran los resultados de Figaro.
Donde: truncLen= Indica la longitud final a la que se recortarán las secuencias
       MaxN= Número máximo permitido de bases ambiguas (base N)
       maxEE= Máximo error esperado permitido en una secuencia antes de truncarla
       truncQ= Calidad mínima para las bases en la secuencia antes de truncarlas.  
```bash
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(300, 250),
                     maxN = 0,
                     maxEE = c(2, 2),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE)
```

- [ ] Observamos los primeros resultados despues del filtrado 

```bash
head(out)
```
**4. Errores**
El algoritmo de Dada2 trabaja con los errores esperados y aprende de ellos 

- [ ] Aprender de la tasa de error 

```bash
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
```

- [ ] Graficamos los errores
Donde la linea roja nos muestra los errores esperados, los errores estimados se observan en linea negra y los errores
observados se muestran con puntos

```bash
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)
```
**5. Variantes** 

- [ ] Aplicamos el algoritmo central de inferencia de muestras a los datos de secuencias filtrados y recortados.

```bash
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
```

- [ ] Observamos los objetos despues de aplicarles el algoritmo

```bash
head(dadaFs[[1]])
```
- [ ] Fusionamos las lecturas

```bash
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
```

- [ ] Inspeccionando nuestras muestras fusionadas

```bash
head(mergers[[2]])
```
**6. Tabla de Variantes** 

- [ ] Construir la tabla de secuencia (Una version nueva de la OTU table)

```bash
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

- [ ] Inspeccionando los largos de la secuencia de distribución

```bash
table(nchar(getSequences(seqtab)))
```
**7. Remover quimeras** 

Donde el resultado "sum(seqtab.nochim)/sum(seqtab)", nos mostrara el porcentaje total de quimeras 
en las muestras 

- [ ] Remover quimeras

```bash
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method ="consensus", 
                                    multithread=TRUE, 
                                    verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
```
**8. Resumen del analisis** 

- [ ] Resumen

```bash
getN <- function(x) sum(getUniques(x))

track<- cbind (out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers,getN), 
               rowSums(seqtab.nochim))


colnames(track)<- c("input", 
                    "filtered",
                    "denoisedF", 
                    "denoisedR", 
                    "merged",
                    "nochim")

rownames(track) <- sample.names

head(track)
```

**9. Asignacion taxonomica** 

- [ ] Cargamos librerias de Silva

```bash
taxa <- assignTaxonomy(seqtab.nochim,
"~/Directorio/donde/este/Silva",
                       multithread = TRUE)
taxa <- addSpecies(taxa, 
  "~/Directorio/donde/este/complemento/de/Silva")
```
- [ ] Eliminando nombre de las filas de secuencias para visualización

```bash
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```
**10. Inferencias para analisis con phyloseq**

- [ ] Cargamos la paqueteria necesaria

```bash
cran_packages <- c("knitr", "qtl", "bookdown", "ggplot2", "grid", "gridExtra", 
                   "tidyverse", "xtable", "devtools", "kableExtra", "remotes")
bioc_packages <- c("phyloseq", "DECIPHER", "phangorn", "ggpubr", 
                   "DESeq2", "genefilter", "philr")
git_packages <- c("btools", "fantaxtic")

sapply(c(cran_packages, bioc_packages, git_packages), 
       require, character.only = TRUE)
```
- [ ] Creamos objeto en phyloseq para inferencias

```bash
file_path <- "Direccion/donde/se/encuentren/meta/datos"
csv_content <- readLines(file_path)
csv_content_clean <- gsub("\uFFFD", "", csv_content)
samdf <- read.csv(text = csv_content_clean, header = TRUE, row.names = 4, 
                  sep = ",")
samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```
