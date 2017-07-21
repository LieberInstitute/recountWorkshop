# recount workshop: Learn to leverage 70,000 human RNA-seq samples for your projects

## Level of difficulty

beginner

## Abstract

The recount2 project re-processed RNA sequencing (RNA-seq) data on over 70,000 human RNA-seq samples spanning a variety of tissues, cell types and disease conditions. Researchers can easily access these data via the `recount` Bioconductor package, and can quickly import gene, exon, exon-exon junction and base-pair coverage data for uniformly processed data from the SRA, GTEx and TCGA projects in R for analysis. This workshop will cover different use cases of the `recount` package, including downloading and normalizing data, processing and cleaning relevant phenotype data, performing differential expression (DE) analyses, and creating reports for exploring the results using other Bioconductor packages. The workshop will also cover how to use the base-pair coverage data for annotation-agnostic DE analyses and for visualizing coverage data for features of interest. After taking this workshop, attendees will be ready to enhance their analyses by leveraging RNA-seq data from 70,000 human samples.

[recount2 website](https://jhubiostatistics.shinyapps.io/recount/), [recount package](http://bioconductor.org/packages/recount), [paper](http://www.nature.com/nbt/journal/v35/n4/full/nbt.3838.html).

## Topic

RNA-seq

## Expected outcomes

Learn how to search projects, download data, explore the metadata, add more phenotype information, and prepare the data for a DE analysis. Then perform the DE analysis with `DESeq2` and explore the results using `regionReport`.
Participant prerequisites: basic familiarity with packages such as `GenomicRanges` and `DESeq2`. Functions from those packages used in the workshop will be briefly described.

## Installation

```{r}
## You will need R 3.4.x and Bioconductor 3.6 (currently Bioc-devel)
library('devtools')
install_github('LieberInstitute/recountWorkshop')
```

## Feedback and bug reports

For feedback and bug reports, please create a [new issue](https://github.com/LieberInstitute/recountWorkshop/issues). Remember to check the [posting guidelines](http://www.bioconductor.org/help/support/posting-guide/) from the Bioconductor support website before submitting your new post. Thank you.
