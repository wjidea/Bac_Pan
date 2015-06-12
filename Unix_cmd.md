## Data analyses for Pan-Core Genome project
**1\. Data clean-up**

Change fasta file with a clean title line for each of the sequence. Then the
modified fasta files will be fed to **"Prokka"** to perform prokayotic genome
annotations.Shell scripts below:

```bash
    #!/bin/bash
    change the first line of fasta file
    var=">Contig1"
    file=$1
    #sed was used to change the first line of each fasta file
    sed "1s/.*/$var/" $file > ${file%%fasta}1.fa
```

If the multi-fasta file has a title more than 20 characters long, prokka will
recommend you change a more readable name to use as an input.

**2\. Prokaryotic genome annotation with _Prokka_**

Prokka will take fasta file as an input to predict CDS, rRNA, CRISPR, and ncRNA
from prokaryotic genome (also works for Archae and virus). The script for genome
annotation is below:

```bash
    #!/bin/bash
    # save this script as a shell script file
    FILE=$1
    prokka --outdir ~/Files/OrthoMCL_turf/Fasta/${FILE%%.fa} \
    --force --prefix ${FILE%%.fa} --locustag ${FILE%%.fa} \
    --kingdom Bacteria --addgenes --cpus 11 $FILE

```

```bash
    # this step may take about 1 hour
    for FILE in *.fasta; do ./prokka_run.sh $FILE; done;
    # copy all the file to
    cp ~/Files/OrthoMCL_turf/Fasta/*/*.gff ../GFF/

```

**3\. Run pan-genome process with _Roary_**

**Roary** is a recently developed tool to study pan and core genomes of prokaryotic
species. Function-wise is very similar to what have been developed form the OrthoMCL
method. But their clustering algorithm prior to all-to-all blastp reduce the
computation time for the orthologs analysis. Script to run **Roary** is below:

```bash
    roary -p 11 -e *.gff  # this step may take about 40-50 min

    # query for union and intersection between groups of bacteria

    # intersections for all turf pathogens
    query_pan_genome -a intersection -g clustered_proteins COLB1.gff INDB2.gff \
      INV.gff KL3.gff MDB1.gff NCT3.gff QH1.gff QHB1.gff Sa2.gff SF12.gff SH7.gff MOR.gff
    cat pan_genome_results | cut -f1 -d : | sort > gene_set1_intersect
    # intersections for all maize pathogens
    query_pan_genome -a intersection -g clustered_proteins AA38.gff AA78-5.gff Aa99-2.gff
    cat pan_genome_results | cut -f1 -d : | sort > gene_set2_intersect

    # union for all turf pathogens
    query_pan_genome -a union -g clustered_proteins COLB1.gff INDB2.gff INV.gff \
      KL3.gff MDB1.gff NCT3.gff QH1.gff QHB1.gff Sa2.gff SF12.gff SH7.gff MOR.gff
    cat pan_genome_results | cut -f1 -d : | sort > gene_set1_union
    # union for all maize pathogens
    query_pan_genome -a union -g clustered_proteins AA38.gff AA78-5.gff Aa99-2.gff
    cat pan_genome_results | cut -f1 -d : | sort > gene_set2_union

```

After extract the gene list from the intersection and union from the data, genes
 unique to maize and genes unique to turf can be extracted from the gene lists.
 Extraction code below:

```bash
    # filter the genes present in all turf pathogens but the maize
    comm -23 gene_set1_intersect gene_set2_union
    # filter the genes present in all maize pathogens but the turf
    comm -23 gene_set2_intersect gene_set1_union
```

Based on the whole genome SNP tree of Acidovorax spp., we found three majoe clusters:
 - Group1: KL3, INV, QHB1, MDB1, NCT3
 - Group2: SH7, QH1, MOR, COLB1, INDB2, SF12
 - Group3: Sa2  
We are intersted looking into the genes shared among those groups

```bash
    # intersections for group1 turf pathogens
    query_pan_genome -a intersection -g clustered_proteins KL3.gff INV.gff \
                     QHB1.gff MDB1.gff NCT3.gff
    cat pan_genome_results | cut -f1 -d : | sort > gene_group1_intersect
    wc -l < gene_group1_intersect
    # intersections for group2 turf pathogens
    query_pan_genome -a intersection -g clustered_proteins SH7.gff QH1.gff MOR.gff \
                     COLB1.gff INDB2.gff SF12.gff
    cat pan_genome_results | cut -f1 -d : | sort > gene_group2_intersect
    wc -l < gene_group2_intersect
    # genes in group3
    cat Sa2.gff | grep -P "\tCDS\t"

```
Pan genome results with genes shared within groups:
| Groups | Intersection |
|--------|--------------|
| 1      | 1444         |
| 2      |              |
| 3      | 3935         |


###TODO
- [x] list of core genes shared within each host of plants
- [x] check the groups within Aa turf
- [x] find gene and gene contents among those bacteria
- [ ] Phylogenomics analysis of the CDSs -- in progress



### Dependencies
1. [Roary](https://github.com/sanger-pathogens/Roary) - Thanks team Roary for
developing such nice tool
2. [Prokka](http://www.vicbioinformatics.com/software.prokka.shtml) - I wish
eukaryotic genome annotation can be as easy as in prokka
