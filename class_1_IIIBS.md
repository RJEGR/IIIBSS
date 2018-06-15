> Group notes at: **http://132.248.34.204:9001/p/EscuelaVerano2018**
>
> details in https://swcarpentry.github.io/shell-novice-es/04-pipefilter/



### 1st chapter:  

Caracteristicas:

 - Multiprocesor: varias tareas al mismo tiempo

 - Multiuser: Varios usuarios conectados simultaneamente y ejecutando procesos.

   

   Porque Unix ?:

   Comandos incorporan facil solucion para manejar datosflujos de trabajos automatizadoscombinación de programas existentes ie. cosas complejas a “pocas” teclas de distanciaProtocolos de transferencia segura de archivos:pyzillaputtymova xtermprocesos y archivosProcesos:

   

   Programa en ejecucion de atributos:nombretamanoUIDEstadoTiempo de CPUetc ..Archivos:Son elementos de almacenamiento de tipo:Binario (información codificada).exe , .com.jpg, .png.mov, .mpg, .mp3.doc, .xls, .pptTexto plano:tipo de documento que contiene texto sin otros atributos tipográficos (color, subrayado, negrita, etc.).Editor de texto plano: vi, emacs, pico, nano ;)

Lets figure out the frequency of regulated genes by proteins (ie. number of proteins P than regulate genes G)

```
grep -v "^#" network_tf_gene.txt | cut -f1,2 | sort -u | cut -f2 | sort | uniq -c | sort -n | cut -c 1-8 | uniq -c >FreqTFGene.txt
```

| Number of genes regulated | `#` Protein |
| ------------------------: | :---------: |
|                       764 |      1      |
|                       409 |      2      |
|                       298 |      3      |
|                       137 |      4      |
|                       103 |      5      |
|                        36 |      6      |
|                        22 |      7      |
|                        12 |      8      |
|                         4 |      9      |
|                        14 |     10      |
|                         1 |     11      |
|                         6 |     13      |
|                         1 |     14      |



# Loops

```
for i in *.dat; do echo $i; done
```

## ciclos anidados

```
for species in cubane etgane methane
do
	for temperature in 25 30 37 40
	do
		mkdir $species-temperature
	done
done
```



