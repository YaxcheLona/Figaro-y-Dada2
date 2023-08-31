
Figaro_Y_Dada2



##



**Qiime2** (Quantitative Insights Into Microbial Ecology) es un *pipeline* desarrollado para el análisis de metataxonomía ([Bolyen et al., 2019](https://www.nature.com/articles/s41587-019-0209-9)). Contiene herramientas para limpiar secuencias, agrupar, asignar taxonomía, reconstruir filogenias, inferir métricas de diversidad, abundancia diferencial, etc. Es de código abierto, posee una [interfaz gráfica](https://view.qiime2.org/) amigable, [mucha documentación](https://docs.qiime2.org/2022.11/plugins/available/diversity/), [tutoriales](https://docs.qiime2.org/2023.7/tutorials/) y [foros de ayuda](https://forum.qiime2.org/).



Recordemos el flujo de análisis que se pueden hacer dentro de este *pipeline*.

<img src="https://github.com/DianaOaxaca/02.Amplicones_16S_Qiime2/blob/main/qiime2_wf.jpg" style="zoom:67%;" />



#### Los datos

Para este taller trabajaremos con datos de amplicones de la región V3-V4 del 16S rRNA de muestras de tres tiempos de fermentación del **pulque**, estos se obtuvieron con una plataforma ILLUMINA MiSeq (2 x 300 pb) y están en formato FASTQ. Los datos fueron depositados en NCBI y ENA bajo el BioProject **PRJEB13870** del artículo **[Deep microbial community profiling along the fermentation process of pulque, a biocultural resource of Mexico](https://www.sciencedirect.com/science/article/pii/S0944501320304614#sec0010)** (Rocha-Arriaga et al., 2020).

Recordemos como es una secuencia en formato **fastq**:

**FASTQ**: Normalmente se compone de cuatro lineas por secuencia:

```
@HWI-D00525:129:C6ACWANXX:5:1106:11467:38087
TGTCAAGACAGGCACGCACGTTCTGAATCAGCCGACGGTCGGTACGTAAGGCCATGTATATCGTGGCGTCCTTGTAAGTGATTTCCTTGCGTCCG
+
CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGFGGGGGFGGGGGGEDGGGG
```

- L1 Comienza con '@' seguido del identificador de la secuencia y una descripción opcional
- L2 La secuencia de nucleótidos
- L3 Comienza con un '+' opcionalmente incluye el identificador de la secuencia
- L4 Indica los valores de calidad de la secuencia, debe contener el mismo número de símbolos que el número de nucleótidos. (Actualmente se usa phred+33, puedes ver explicaciones más detalladas [aquí](https://en.wikipedia.org/wiki/Phred_quality_score). 





> **NOTA**: Cada programa tiene una ayuda y un manual de usuario, es **importante** revisarlo y conocer cada parámetro que se ejecute. En terminal se puede consultar el manual con el comando `man` y también se puede consultar la ayuda con `-h` o `--help`, por ejemplo `qiime tools import --help`.
>
> La presente práctica sólo es una representación del flujo de trabajo para el análisis de amplicones, sin embargo, no sustituye a los manuales de cada programa y el flujo puede variar dependiendo del tipo de datos y pregunta de investigación. Una explicación mucho más detallada de cada paso se encuentra en el *[overview](https://docs.qiime2.org/2023.7/tutorials/overview/)*  y en el *[moving-pictures](https://docs.qiime2.org/2023.7/tutorials/moving-pictures/)* de QIIME2.



### Manos a la obra

QIIME2 se puede correr en un ambiente **Conda**, las instrucciones de instalación para cada sistema operativo las puedes encontrar [aquí.](https://docs.qiime2.org/2023.7/install/native/#install-qiime-2-within-a-conda-environment)



- [ ] Entra al servidor

```bash
ssh -X -Y USUARIO@132.248.15.30 -p NÚMERO_PUERTO -o ServerAliveInterval=60
```

- [ ] Dirígete al espacio en `/botete/` en el que siempre estarás trabajando

```bash
cd /botete/USUARIO
```

- [ ] Inicia conda

```bash
/opt/anaconda3/bin/conda init bash
```

- [ ] Sal del servidor
- [ ] Vuelve a entrar y entra a tu espacio en `/botete/usuario`
- [ ] Activa el ambiente de Qiime2

```bash
conda activate /botete/diana/.conda/envs/qiime2-2022.11
```

- [ ] Vamos a crear el directorio y espacio de trabajo

```bash
mkdir -p 02.Amplicones_16S_Qiime2/{data,src,results}
cd 02.Amplicones_16S_Qiime2/
```

- [ ] Crea una liga simbólica a los datos, copia el archivo de metadata y los scripts.

```bash
# Liga simbólica
ln -s /botete/diana/Hackeando_las_comunidades_microbianas_v1/02.Amplicones_16S_Qiime2/data/*.gz data/
# copiar el metadata
cp /botete/diana/Hackeando_las_comunidades_microbianas_v1/02.Amplicones_16S_Qiime2/data/metadata.tsv data/
# copiar los scripts
cp /botete/diana/Hackeando_las_comunidades_microbianas_v1/02.Amplicones_16S_Qiime2/src/* src/
# Copiar la base de datos
cp /botete/diana/Hackeando_las_comunidades_microbianas_v1/02.Amplicones_16S_Qiime2/data/*.qza data/

```

Ya que copiaste los archivos necesarios, lista el contenido de tus directorio `data` , veamos que contiene el archivo `metadata.tsv`

Recordemos que Qiime2 requiere de un archivo `manifest` que contenga la ubicación e información de los archivos `fastq`. Así que ejecuta el `script 02` para crear tu archivo manifest:

```bash
bash src/02.create_manifest.sh
```

- [ ] Observa que contiene el archivo `manifest.csv`

Ahora si tenemos todos los datos de entrada para comenzar con el *pipeline*.

1. **Importar los datos**

   Observa el contenido de `03.import_data.sh` y ve a la ayuda de qiime para conocer cada uno de los *plugins* que estamos poniendo. 

   Para conocer más el tipo de datos que se pueden importar, puedes visitar [esta página.](https://docs.qiime2.org/2023.7/tutorials/importing/)

```bash
bash src/03.import_data.sh
```

- [ ] Veamos como se ven las secuencias que estamos analizando. Puedes descargar el archivo `results/01.demux.qzv` que acabas de generar o puedes descargarlo desde este repositorio y verlo en [QIIME2view](https://view.qiime2.org/).

2. **Remover adaptadores**

   Como no hemos hecho ningún preprocesamiento a los datos y sabemos que nuestras secuencias tienen adaptadores, los removeremos con `cutadapt` dentro de QIIME2.

   - [ ] Primero observa la ayuda de este paso de limpieza.

     ```bash
     qiime cutadapt trim-paired --help
     ```

   - [ ] Ahora si removamos los adaptadores, puedes modificar, agregar parámetros al script para mejorar.

   ```bash
   bash src/04.clean_adapter.sh
   ```

   Observemos el archivo `02.demux_clean_adapter.qzv` , ¿qué ocurrió con la calidad?

3. **Eliminación de ruido y agrupación con Dada2**

   Antes de comenzar el *denoising* , veamos la ayuda.

   ```bash
   qiime dada2 denoise-paired --i-demultiplexed-seqs --help
   ```

   Como habrás notado, es necesario tomar decisiones basadas en la calidad de nuestras lecturas para definir los valores de truncado.

   - [ ] Edita el script `src/05.denoising.sh` y modifica los valores de truncado del *denoising* y/o los que consideres que pueden mejorar los resultados.

   ```bash
   bash src/05.denoising.sh
   ```

   - [ ] **Opcional**: Si realizaste más de una versión extra a las 3 versiones propuestas, puedes obtener sus estadisticos con el siguiente script.

   ```bash
   #bash src/06.get_denoising_info.sh
   ```

   - [ ] O puedes obtener los estadísticos de la versión 4 que creaste con la siguiente línea

   ```bash
   qiime tools export --input-path results/03.denoising-stats_v4.qza --output-path results/03.denoising-stats_v4
   ```

   - [ ] Observa  los estadísticos de los resultados que obtuviste y compáralos con los de las versiones v1-v3, para ello, copia los resultados de las versiones 1 a 3 y de la nueva versión 4 que ya hiciste, a un nuevo directorio de estadísticos de las versiones de denoising.

   ```bash
   mkdir -p results/stats_versions
   cp /botete/diana/Hackeando_las_comunidades_microbianas_v1/02.Amplicones_16S_Qiime2/results/03.denoising-stats_v*/*stats_v* results/stats_versions/
   cp results/03.denoising-stats_v4/stats.tsv results/03.denoising-stats_v4/stats_v4.tsv
   ```

   - [ ] Vamos a compararlos

   ```bash
   cat results/stats_versions/*.tsv
   ```

   Antes de conocer el resultado que obtuviste, al comparar las tres versiones, seleccionamos la v3 por ser con la que se recuperó un mayor número de ASVs. Así que en adelante trabajaremos con esta versión.

   

4. **Asignación taxonómica**

   Utilizaremos `sklearn` para realizar la asignación taxonómica, por lo tanto utilizaremos una base de datos preentrenada. La puedes encontrar [aquí](https://docs.qiime2.org/2022.11/data-resources/).

   Recuerda que puedes modificar el script con la versión que obtuviste o cambiando algún parámetro.

   ```bash
   bash src/07.tax_assign.sh
   ```

   

### Git

```bash
git config --global user.name "tuNombre"
git config --global user.email "correo@gmail.com
#entra al directorio de trabajo
git init
nano .gitignore
git add .gitignore
git add *
git commit -m "create repository"
git branch -M main
git remote add origin https://github.com/DianaOaxaca/02.Amplicones_16S_Qiime2.git

git push -u origin main
#nombre de usuario
#public-key
```



#### Referencias

Bolyen, E., Rideout, J.R., Dillon, M.R. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852–857 (2019). https://doi.org/10.1038/s41587-019-0209-9

Rocha-Arriaga, C., Espinal-Centeno, A., Martinez-Sanchez, S., Caballero-Pérez, J., Alcaraz, L. D., & Cruz-Ramírez, A. (2020). Deep microbial community profiling along the fermentation process of pulque, a biocultural resource of Mexico. *Microbiological Research*, *241*, 126593
