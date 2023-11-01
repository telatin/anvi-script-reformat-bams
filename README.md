

### Test files

The **test** directory contains a simple reference (FASTA) and a paired-end dataset (FASTQ).

```text
├── bin
│   └── anvi-reformat-bam
└── test
    ├── input
    │   ├── reads_R1.fq.gz
    │   ├── reads_R2.fq.gz
    │   └── ref.fa
    └── output
        └── formatted.bam
```

Running the `Makefile` will:

* index the reference (`bwa index``), 
* align the reads (`bwa mem``), 
* reformat the fasta file (`script-reformat-fasta``),
* and reformat the BAM file (using this repository's scripts).

```text
└── test
    ├── REFORMAT.fa
    ├── REPORT_FILE.txt
    ├── input
    │   ├── raw.bam
    │   ├── raw.bam.csi
    │   ├── reads_R1.fq.gz
    │   ├── reads_R2.fq.gz
    │   ├── ref.fa
    │   ├── ref.fa.amb
    │   ├── ref.fa.ann
    │   ├── ref.fa.bwt
    │   ├── ref.fa.pac
    │   └── ref.fa.sa
    └── output
        └── formatted.bam
```