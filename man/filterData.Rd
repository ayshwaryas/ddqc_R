% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ddqc.R
\name{filterData}
\alias{filterData}
\title{Filter the Seurat object}
\usage{
filterData(data, df.qc)
}
\arguments{
\item{data}{Seurat object}

\item{df.qc}{result of ddqc.metrics}
}
\value{
Filtered Seurat object
}
\description{
This function filters Seurat object based on df.qc
}
