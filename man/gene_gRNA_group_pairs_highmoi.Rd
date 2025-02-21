% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_highmoi_old.R
\docType{data}
\name{gene_gRNA_group_pairs_highmoi}
\alias{gene_gRNA_group_pairs_highmoi}
\title{Gene-gRNA group pairs}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 120 rows and 3 columns.
}
\usage{
data(gene_gRNA_group_pairs_highmoi)
}
\description{
The pairs of genes and gRNA groups that we seek to test for association. Columns \code{gene_id} (required), \code{gRNA_group} (required), and \code{pair_type} (optional).
}
\section{Overview}{

This object is part of the example data. The data are taken from the paper "A genome-wide framework for mapping gene regulation via cellular genetic screens" by Gasperini et al., 2019. The authors used a CRISPRi-based assay to target 5,000+ putative enhancers in a population of K562 cells. The authors additionally targeted 200+ gene transcription start sites (TSSs) and designed a library of 50 non-targeting gRNAs to serve as negative controls. Genes, gRNAs, and cells are all down-sampled to reduce the size of the data.


This object is part of the example data. The data are taken from the paper "A genome-wide framework for mapping gene regulation via cellular genetic screens" by Gasperini et al., 2019. The authors used a CRISPRi-based assay to target 5,000+ putative enhancers in a population of K562 cells. The authors additionally targeted 200+ gene transcription start sites (TSSs) and designed a library of 50 non-targeting gRNAs to serve as negative controls. Genes, gRNAs, and cells are all down-sampled to reduce the size of the data.
}

\keyword{datasets}
