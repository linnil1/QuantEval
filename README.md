# QuantEval
## About
I rewrite the [original Code](https://github.com/dn070017/QuantEval) to make it cleaner.

Two main features are implemented in this repo:

1. Run the pipeline of quantification on de-novo transcriptome assembly.
2. Get the quantification results by building a ambiguity cluster based on connected components.

Both methods are explained in the following study.
> Ping-Han Hsieh, Yen-Jen Oyang and Chien-Yu Chen. Effect of de novo transcriptome assembly on transcript quantification. Scientific Reports volume 9, Article number: 8304 (2019).

## Requirement
I test successfully with below requirments

* Python:
   * Python3.8.5
   * Packages: pandas(1.1.0), numpy(1.19.1), biopython(1.77), scipy(1.52), PyYAML (5.3.1), ipython(7.17.0)
   * Install it with `pip3 install -r requirments.txt`
* Pipeline tools:
    * All tools are ran in docker(Podman 1.6.4)
    * wget
* Generate figures and table:(Unreviewed by linnil1)
   * R (3.3.0)
   * R pacakges: gridExtra, grid, stats, tidyverse, plyr, ggplot2, reshape2

## Usage
### Run pipelines
**It takes a long time and large space to run.**
```
ipython pipeline.ipy -- --method=all
```

All the configuration (include tool version, parameters and data source) are written in [`metadata/meta.yml`](https://github.com/linnil1/QuantEval/blob/master/metadata/meta.yaml).

Feel free to modify it.

The result of one of species(yeast) and one of dataset(lower depth simulation)
should look like this:

```
data/yeast/
├── chromosome.fasta
├── fastqc
│   ├── simlow_r1_fastqc.html
│   ├── simlow_r1_fastqc.zip
│   ├── simlow_r2_fastqc.html
│   └── simlow_r2_fastqc.zip
├── mRNA.blast.self.tsv
├── mRNA.fasta
├── mRNA.gtf
├── simlow.mRNA.answer.tsv
├── simlow.mRNA.fasta -> mRNA.fasta
├── simlow.mRNA.blast.self.tsv-> mRNA.blast.self.tsv
├── simlow.mRNA.kallisto.tsv
├── simlow.mRNA.rsem.tsv
├── simlow.mRNA.salmon.tsv
├── simlow.mRNA.transrate.csv
├── simlow_r1.fastq
├── simlow_r1.trim.fastq
├── simlow_r2.fastq
├── simlow_r2.trim.fastq
├── simlow.rnaspades.blast.contig_to_mrna.tsv
├── simlow.rnaspades.blast.mrna_to_contig.tsv
├── simlow.rnaspades.blast.self.tsv
├── simlow.rnaspades.fasta
├── simlow.rnaspades.kallisto.tsv
├── simlow.rnaspades.rsem.tsv
├── simlow.rnaspades.salmon.tsv
├── simlow.rnaspades.transrate.csv
├── simlow.transabyss.blast.contig_to_mrna.tsv
├── simlow.transabyss.blast.mrna_to_contig.tsv
├── simlow.transabyss.blast.self.tsv
├── simlow.transabyss.fasta
├── simlow.transabyss.kallisto.tsv
├── simlow.transabyss.rsem.tsv
├── simlow.transabyss.salmon.tsv
├── simlow.transabyss.transrate.csv
├── simlow.trinity.blast.contig_to_mrna.tsv
├── simlow.trinity.blast.mrna_to_contig.tsv
├── simlow.trinity.blast.self.tsv
├── simlow.trinity.fasta
├── simlow.trinity.kallisto.tsv
├── simlow.trinity.rsem.tsv
├── simlow.trinity.salmon.tsv
├── simlow.trinity.transrate.csv
├── simulation
│   ├── flux_simulator_clean_sorted.gtf
│   ├── flux_simulator.gtf
│   ├── flux_simulator_yeast_low.lib
│   └── flux_simulator_yeast_low.pro
└── trimmed
    ├── simlow_r1.unpaired.fastq
    └── simlow_r2.unpaired.fastq
```

### Run QuantEval Main Program
In brief, TPM/counts will be calculated based on each components, compared to normal method that show tpm/counts for each transcript individually.

Three files are needed before running QuantEval
* Abundance files(Of course)
* Self-Blast result(In tsv, aka blast format 6)
* Transrate result

#### Abundance files
Run quantification and get the abundance file(e.g. kallisto)
```bash
# contig.fasta is the file you assembled by de-novo tools
# read_r1.fastq and read_r2.fastq is the original read files
# Abundance file: kallisto_abundance.tsv
kallisto index -i ./kallisto/kallisto.index -k 31 contig.fasta
kallisto quant -i ./kallisto/kallisto.index -o ./kallisto read_r1.fastq read_r2.fastq
mv ./kallisto/abundance.tsv kallisto_abundance.tsv
```

#### Self-Blast
Run pairwise BLASTn
```bash
blastn -db contig.fasta -query contig.fasta -outfmt 6 -evalue 1e-5 -perc_identity 95 -out ./contig.self.tsv
```

#### Transrate
Run TransRate to get contig information
```bash
transrate --assembly contig.fasta --output ./transrate/ --left read_r1.fastq --right read_r2.fastq
mv ./transrate/contigs.csv transrate.csv
```

#### Main
After all input filess are ready, we can run the main program.
```bash
# Preprocess the abundance tsv to simplified table
# add --merge can merge multiple tsv together
# Only tsv of quantifier: kallisto, rsem, salmon are allowed
python3 QuantEval.py pre_abundance
    --abundance kallisto_abundance.tsv
    --quantifier kallisto
    --merge
    --output    abundance.tsv

# main
python3 QuantEval.py contig
    --abundance abundance.tsv
    --transrate transrate.csv
    --blast     contig.self.tsv
    --output    quanteval.contig.tsv
```

#### QuantEval Output format
`quanteval.contig.tsv`

| column | description |
|--------|-------------|
| match_name | alignment (ref.contig.strand) |
| contig_name | target contig name |
| ref_name | target ref name |
| accuracy | accuracy of the alignment |
| recovery | recovery of the alignment |
| ***contig/ref***\_length | length of ***contig/ref*** |
| ***contig/ref***\_tr\_***transrate_score*** | ***transrate score*** of ***contig/ref*** |
| ***contig/ref***\_xprs\_***tpm/count***_***quantifier*** | quantification result of ***contig/ref*** |
| ***contig/ref***\_component | label of connected component of ***contig/ref*** |
| ***contig/ref***\_component_size | number of sequences in the connected component of ***contig/ref*** |
| ***contig/ref***\_component_contribute_xprs_***tpm/count***\_***quantifier*** | proportion of ***TPM/read count (RPEA)*** in the connected component of ***contig/ref*** |
| ***contig/ref***\_component_relative_xprs_***tpm/count***\_***quantifier*** | ***TPM/count*** of ***contig/ref*** / highest ***TPM/count*** in the same connected component |
| ***contig/ref***\_component_max_xprs_***tpm/count***\_***quantifier*** | highest ***TPM/count*** of ***contig/ref*** in the same connected component |
| ***contig/ref***\_component_avg_xprs_***tpm/count***\_***quantifier*** | average ***TPM/count*** of ***contig/ref*** in the same connected component |
| ***contig/ref***\_component_tot_xprs_***tpm/count***\_***quantifier*** | total ***TPM/count*** of ***contig/ref*** in the same connected component |
| ref_gene_contribute_xprs_***tpm/count***\_***quantifier*** | proportion of ***TPM/read count*** in the gene of ref |
| ref_gene_relative_xprs_***tpm/count***\_***quantifier*** | ***TPM/count*** of ref / highest ***TPM/count*** in the same gene |
| ref_gene_max_xprs_***tpm/count***\_***quantifier*** | highest ***TPM/count*** of ref in the same gene |
| ref_gene_avg_xprs_***tpm/count***\_***quantifier*** | average ***TPM/count*** of ref in the same gene |
| ref_gene_tot_xprs_***tpm/count***\_***quantifier*** | total ***TPM/count*** of ref in the same gene |
| length_difference | the difference of length between contig and reference |
| xprs_***tpm/count***\_error_***quantifier*** | quantificaion error for the estimated abundance of contig |
> Note that this is the superset of the output fields (the match mode). The content of output will be different depends on reference/contig/match mode (e.g. one can only find the columns start with ***contig*** from contig mode), but one can find all the description on the table above.
