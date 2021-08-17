# About

This repository provides supplementary code for the paper:

Wingfield, B., Lapsley, C., McDowell, A. et al. Variations in the oral
microbiome are associated with depression in young adults. Sci Rep 11, 15009
(2021). https://doi.org/10.1038/s41598-021-94498-6

The 16S sequence data for this study have been deposited in the European
Nucleotide Archive (ENA) at EMBL-EBI under accession number PRJEB44565
(https://www.ebi.ac.uk/ena/browser/view/PRJEB44565)

# Running the pipeline

This pipeline reproduces the plots and key findings from the paper described
above. To run this pipeline you'll need to have the following installed:

* Nextflow
* Docker
* R and [qiime2R](https://github.com/jbisanz/qiime2R)

Running the pipeline should be simple:

`nextflow run main.nf`

The `results/` directory will contain all the summary figures. 

