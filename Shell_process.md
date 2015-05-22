### Data analyses for Pan-Core Genome project
1\. Data clean-up

Change fasta file with a clean title line for each of the sequence. Then the modified fasta files will be fed to **"Prokka"** to perform prokayotic genome annotations.Shell scripts below:

```bash
    #!/bin/bash
    # change the first line of fasta file
    var=">Contig1"
    file=$1
    #sed was used to change the first line of each fasta file
    sed "1s/.*/$var/" $file > ${file%%fasta}1.fa
```

If the multi-fasta file has a title more than 20 characters long, prokka will recommend you change a more readable name to use as an input.

2\. Prokaryotic genome annotation with **Prokka**

Prokka will take fasta file as an input to predict CDS, rRNA, CRISPR, and ncRNA from prokaryotic genome (also works for Archae and virus). The script for genome annotation is below:

```bash
    #!/bin/bash
    FILE=$1
    prokka --outdir ~/Files/OrthoMCL_turf/Fasta/${FILE%%.fa} \
    --force --prefix ${FILE%%.fa} --locustag ${FILE%%.fa} \
    --kingdom Bacteria --addgenes $FILE
    # copy all the file to
    cp ~/Files/OrthoMCL_turf/Fasta/*/*.gff ../GFF/
```

3\. Run pan-genome process with **Roary**
**Roary** is a recently developed tool to study pan and core genomes of prokaryotic species. Function-wise is very similar to what have been developed form the OrthoMCL method. But their clustering algorithm prior to all-to-all blastp reduce the computation time for the orthologs analysis. Script to run **Roary** is below:

```bash
roary -p 11 -e *.gff  # this step may take about 40-50 min

# query for union and intersection between groups of bacteria

# intersections for all turf pathogens
query_pan_genome -a intersection -g clustered_proteins COLB1.gff INDB2.gff INV.gff KL3.gff MDB1.gff NCT3.gff QH1.gff QHB1.gff Sa2.gff SF12.gff SH7.gff MOR.gff
cat pan_genome_results | cut -f1 -d : | sort > gene_set1_intersect
# intersections for all maize pathogens
query_pan_genome -a intersection -g clustered_proteins AA38.gff AA78-5.gff Aa99-2.gff
cat pan_genome_results | cut -f1 -d : | sort > gene_set2_intersect

# union for all turf pathogens
query_pan_genome -a union -g clustered_proteins COLB1.gff INDB2.gff INV.gff KL3.gff MDB1.gff NCT3.gff QH1.gff QHB1.gff Sa2.gff SF12.gff SH7.gff MOR.gff
cat pan_genome_results | cut -f1 -d : | sort > gene_set1_union
# union for all maize pathogens
query_pan_genome -a union -g clustered_proteins AA38.gff AA78-5.gff Aa99-2.gff
cat pan_genome_results | cut -f1 -d : | sort > gene_set2_union

```

After extract the gene list from the intersection and union from the data, genes unique to maize and genes unique to turf can be extracted from the gene lists. Extraction code below:

```bash
comm -23 gene_set1_intersect gene_set2_union

comm -23 gene_set1_intersect gene_set2_union
```
