### Data analyses for Pan-Core Genome project
1\. Data clean-up

Change fasta file with a clean title line for each of the sequence. Then the modified fasta files will be fed to **"Prokka"** to perform prokayotic genome annotations.Shell scripts below:

```bash
    #!/bin/bash
    # change the first line of fasta file
    var=">Contig1"
    file=$1
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
```
