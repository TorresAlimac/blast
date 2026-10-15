# Sesión 08: BLAST local y candidatos LGR en Deroceras laeve

## Objetivo

En esta práctica comparé proteínas de referencia de la vía WNT con el proteoma de la babosa Deroceras laeve. El objetivo fue encontrar proteínas similares, seleccionar candidatos relacionados con LGR y ubicarlos en el genoma. Después alineé sus secuencias y realicé una búsqueda de vuelta para comparar a qué referencia se parecía más cada candidato.

## Programas y datos

Trabajé en el servidor ken, en /mnt/data/bioinfo1/ctorres/sesion08/.

- BLASTP: 2.16.0+.
- MAFFT: v7.526 (2024/Apr/26).
- También utilicé SeqKit y los comandos cut, sort, uniq y awk.

Los datos se copiaron de /mnt/data/bioinfo1/jmiranda/sesion08/datos/

| Archivo                  | Contenido |
|--------------------------|-----------|
| proteoma_dlaeve_v1.2.faa | Proteoma de D. laeve, versión 1.2, con 43,905 secuencias |
| consultas_wnt.faa        | 72 proteínas de referencia correspondientes a 18 genes y cuatro especies |
| mapa_mrna.tsv            | Correspondencia entre identificadores, cromosomas, posiciones y anotaciones |

## Preparación y búsqueda

Construí una base de proteínas con makeblastdb y activé -parse_seqids para poder recuperar después las secuencias por su identificador. El número de secuencias coincidió entre makeblastdb, el conteo de encabezados del FASTA y blastdbcmd: 43,905. La base contiene 15,023,582 residuos en total.

```bash
module load blast+ anaconda3 seqkit
source activate prokka-1.14.6
blastp -version
mafft --version

mkdir -p datos
cp /mnt/data/bioinfo1/jmiranda/sesion08/datos/* datos/

makeblastdb -in datos/proteoma_dlaeve_v1.2.faa -dbtype prot -parse_seqids \
  -out proteoma_dlaeve -title "Proteoma D. laeve v1.2"

grep -c '^>' datos/proteoma_dlaeve_v1.2.faa
blastdbcmd -db proteoma_dlaeve -info

blastp -query datos/consultas_wnt.faa -db proteoma_dlaeve \
  -evalue 1e-5 -outfmt "6 std qcovs" -num_threads 1 -out hits.tsv
```

La salida contiene las 12 columnas estándar de BLAST y una columna adicional de cobertura de la consulta. Obtuve 7,935 filas de alineamientos, en las que aparecen 632 proteínas distintas de la babosa. Una misma proteína puede aparecer en varias filas, por consultas diferentes o por varios segmentos de alineamiento local.

## Ejercicios 1 y 2: genes y número de alineamientos

Los 18 genes de referencia aparecen en hits.tsv.

| Gen de referencia | Filas de alineamientos |
|-------------------|------------------------|
| CELSR1            | 2513                   |
| LGR6              | 1481                   |
| LGR5              | 1475                   |
| RXFP1             | 1401                   |
| FSHR              | 218                    |                 
| PRICKLE1          | 152                    |             
| PRICKLE2          | 147                    |
| DVL2              | 81                     |
| FZD7              | 78                     |
| WNT10A            | 76                     |
| WNT11             | 76                     |
| WNT6              | 76                     |
| FZD5              | 70                     |
| TCF7L2            | 36                     |
| CTNNB1            | 27                     |
| PIWI              | 16                     |
| VANGL1            | 8                      |
| RSPO2             | 4                      |

CELSR1 produjo más filas y RSPO2 menos. Esto no significa que CELSR1 sea más importante. El número depende de las proteínas similares presentes, de los dominios compartidos y de los segmentos que BLAST logra alinear. Por ejemplo, las cuatro filas de RSPO2 corresponden a una sola proteína de la babosa comparada con cuatro referencias.

## Ejercicio 3: mejor alineamiento por bit score

El mejor alineamiento de toda la tabla fue el de CELSR1 de Bradybaena similaris, accession BFZ08600.1, contra dlaeveM00000596688. Tuvo 78.197 % de identidad, un bit score de 4728 y 99 % de cobertura de la consulta. El E-value se reportó como 0.0; esto representa un valor extremadamente pequeño en la salida, no que se haya demostrado una probabilidad exactamente igual a cero.

## Ejercicio 4: filtro por identidad y cobertura

Al pedir identidad de al menos 50 % y cobertura de consulta de al menos 80 %, quedaron 90 de las 7,935 filas. La fracción que pasó fue aproximadamente 0.0113, equivalente a 1.13 %.

## Ejercicios 5 y 6: mejor hit por consulta

Ordené los alineamientos por consulta y por bit score descendente, y conservé la primera fila de cada consulta en mejores.tsv. Obtuve 72 filas. La comparación de identificadores no mostró consultas sin hit, por lo que las 72 referencias tuvieron al menos una coincidencia con los parámetros usados.

Las cuatro referencias de RSPO2 tuvieron como mejor hit a dlaeveM00000555743.

| Referencia RSPO2     | Identidad (%) | E-value  | Cobertura de consulta (%) |
|----------------------|---------------|----------|---------------------------|
| Humano               | 41.000        | 1.90e-44 | 80                        |
| Bradybaena similaris | 75.152        | 9.59e-95 | 56                        |
| Magallana gigas      | 48.598        | 4.23e-70 | 72                        |
| Pomacea canaliculata | 56.621        | 2.14e-84 | 73                        |

En LGR5, las referencias no apuntaron todas a la misma proteína.

| Referencia LGR5      | Mejor hit en la babosa | Identidad (%) | E-value   | Cobertura de consulta (%) |
|----------------------|------------------------|---------------|-----------|---:|
| Humano               | dlaeveM00000441839     | 33.908        | 1.11e-143 | 89 |
| Bradybaena similaris | dlaeveM00000727719     | 83.243        | 0.0       | 67 |
| Magallana gigas      | dlaeveM00000727719     | 50.166        | 0.0       | 81 |
| Pomacea canaliculata | dlaeveM00000727719     | 58.758        | 0.0       | 79 |

La referencia de Bradybaena similaris tuvo la identidad más alta para ese candidato. Su cercanía evolutiva con la babosa ayuda a detectar similitudes, aunque también hay que revisar cobertura y puntaje. La cobertura más alta no siempre corresponde a la identidad más alta.

## Ejercicio 7: selección de candidatos

En la búsqueda inicial, 239 proteínas distintas de la babosa tuvieron algún alineamiento con las referencias de LGR5. Al pedir E-value menor que 1e-30 quedaron 12. Al agregar cobertura de consulta de al menos 60 % quedaron seis, cuyos identificadores guardé en lgr_candidatos_ids.txt.

```bash
awk -F'\t' '$1 ~ /^LGR5\|/ && $11<1e-30 && $13>=60 {print $2}' hits.tsv \
  | sort -u > lgr_candidatos_ids.txt
```

El filtro de cobertura eliminó coincidencias que no cubrían una proporción suficiente de la referencia. Una región conservada puede dar un E-value bajo sin que el resto de la proteína tenga la misma arquitectura que LGR5.

## Ejercicio 8: ubicación de los candidatos

Relacioné los identificadores de los seis candidatos con mapa_mrna.tsv para recuperar cromosoma, posición y símbolo anotado. Las posiciones se expresan en megabases y están redondeadas a un decimal.

| Candidato          | Cromosoma | Posición (Mb) | Símbolo anotado |
|--------------------|-----------|---------------|-----------------|
| dlaeveM00000132191 | chr2      | 15.7          | GPRGPH          |
| dlaeveM00000134217 | chr2      | 17.7          | GPRGPH          |
| dlaeveM00000308285 | chr4      | 45.9          | TSHR            |
| dlaeveM00000351276 | chr5      | 18.7          | GPRGPH          |
| dlaeveM00000441839 | chr6      | 30.9          | LGR4            |
| dlaeveM00000727719 | chr10     | 23.7          | LGR4            |

### El caso de dlaeveM00000909458

Este gen no apareció en lgr_candidatos_ids.txt. Contra LGR5 humano se reportaron dos segmentos de alineamiento local: uno con 26.593 % de identidad y E-value de 8.39e-14, y otro con 30.734 % de identidad y E-value de 2.60e-12. La cobertura de consulta reportada fue de 43 %. No pasó ni el umbral de E-value menor que 1e-30 ni la cobertura mínima de 60 %.

Esta proteína tiene repeticiones ricas en leucina (LRR), pero carece del dominio de siete hélices transmembranales (7TM). Por eso puede coincidir con una parte de LGR5 sin ser un receptor LGR completo. La ausencia de 7TM es información proporcionada en el material; no la determiné con el filtro de BLAST.

## Alineamiento de la familia

Recuperé las seis proteínas candidatas con blastdbcmd y las junté con las cuatro referencias de LGR5. Después ejecuté MAFFT.

```bash
seqkit grep -nrp '^LGR5\|' datos/consultas_wnt.faa > referencias_lgr5.faa
blastdbcmd -db proteoma_dlaeve -entry_batch lgr_candidatos_ids.txt > candidatos_lgr.faa
cat candidatos_lgr.faa referencias_lgr5.faa > familia_lgr.faa
mafft --auto familia_lgr.faa > familia_lgr.aln.faa
seqkit stats familia_lgr.faa familia_lgr.aln.faa
```

El conjunto contiene 10 secuencias. Antes del alineamiento sus longitudes iban de 892 a 1,682 aminoácidos. Después, todas quedaron con 2,099 posiciones, incluyendo los gaps. El alineamiento tiene, por lo tanto, 10 secuencias y 2,099 columnas. Los guiones permiten colocar regiones comparables en las mismas columnas; no son aminoácidos nuevos agregados a las proteínas.

## Bonus: búsqueda de vuelta

Construí una base con las 72 referencias y busqué contra ella las seis proteínas candidatas. Conservé el mejor alineamiento por bit score para cada una.

```bash
makeblastdb -in datos/consultas_wnt.faa -dbtype prot -out consultas
blastp -query candidatos_lgr.faa -db consultas -evalue 1e-5 \
  -outfmt "6 std qcovs" -num_threads 1 -out vuelta.tsv
sort -t$'\t' -k1,1 -k12,12nr vuelta.tsv \
  | awk -F'\t' '!v[$1]++' | cut -f1,2,3,13
```

En los seis casos, el mejor hit fue una referencia de Bradybaena similaris.

| Candidato          | Mejor referencia | Identidad (%) | Cobertura del candidato (%) |
|--------------------|------------------|---------------|-----------------------------|
| dlaeveM00000132191 | RXFP1            | 53.315        | 94                          |
| dlaeveM00000134217 | RXFP1            | 73.164        | 99                          |
| dlaeveM00000308285 | FSHR             | 67.825        | 100                         |
| dlaeveM00000351276 | RXFP1            | 70.838        | 99                          |
| dlaeveM00000441839 | LGR5             | 56.807        | 82                          |
| dlaeveM00000727719 | LGR5             | 83.243        | 99                          |

En esta búsqueda las consultas son las proteínas de la babosa. Por eso la cobertura se refiere al candidato y no a la referencia LGR5 usada en la primera búsqueda. El mejor hit está limitado al conjunto de 72 referencias utilizado.

## Interpretación

La búsqueda inicial dejó seis candidatos parecidos a LGR5, pero no todos se clasificaron de la misma manera en la búsqueda de vuelta. Dos, dlaeveM00000441839 y dlaeveM00000727719, tuvieron como mejor referencia LGR5 y están anotados como LGR4, en chr6 a 30.9 Mb y chr10 a 23.7 Mb. Son los dos candidatos principales al grupo LGR4/5/6. Los otros tres candidatos prefirieron RXFP1 y uno prefirió FSHR. Así, seis coincidencias filtradas no equivalen a seis copias de LGR5. Para afirmar que los dos candidatos principales son receptores completos hace falta confirmar que tienen tanto LRR como 7TM; el alineamiento y BLAST por sí solos no verifican ambos dominios. El gen dlaeveM00000909458 se excluye porque, según la práctica, carece de 7TM y porque sus alineamientos no cumplen los filtros utilizados, aunque presenta similitud significativa en una región.








































