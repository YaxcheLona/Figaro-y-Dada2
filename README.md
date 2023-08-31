
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

