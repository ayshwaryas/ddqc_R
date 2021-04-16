.clusterData <- function(data, norm.factor=10000, n.pieces=50, res=1, random.seed=29) {
  set.seed(random.seed)
  data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = norm.factor)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x = data)
  data <- ScaleData(data, features = all.genes)
  data <- RunPCA(data, npcs=n.pieces, features = VariableFeatures(data))
  data <- FindNeighbors(data, dims = 1:n.pieces)
  data <- FindClusters(data, resolution = res)
  return(data)
}


.metricFilter <- function(data, df.qc, param=2, metric.name, do.upper.co=FALSE, do.lower.co=FALSE,
                          lower.bound=10E10, upper.bound=-10E10) {
  passed.qc <- vector(mode="logical", length=length(colnames(data)))
  names(passed.qc) <- colnames(data)

  df.qc[[metric.name]] <- data[[metric.name]][[metric.name]]
  df.qc[[paste0(metric.name, ".upper.co")]] <- NaN
  df.qc[[paste0(metric.name, ".lower.co")]] <- NaN
  df.qc[[paste0(metric.name, ".passed.qc")]] <- FALSE

  for (cl in levels(data$seurat_clusters)) {
    idx <- data$seurat_clusters == cl
    values = data[ ,idx][[metric.name]][[metric.name]]

    median.v <- median(values)
    mad.v <- mad(values)
    lower.co <- min(median.v - param * mad.v, lower.bound)
    upper.co <- max(median.v + param * mad.v, upper.bound)

    qc.pass.cl <- vector(mode="logical", length=length(values))
    qc.pass.cl <- TRUE

    if (do.lower.co) {
      qc.pass.cl <- qc.pass.cl & (values >= lower.co)
      df.qc[idx, ][[paste0(metric.name, ".lower.co")]] <- lower.co
    }
    if (do.upper.co) {
      qc.pass.cl <- qc.pass.cl & (values <= upper.co)
      df.qc[idx, ][[paste0(metric.name, ".upper.co")]] <- upper.co
    }

    df.qc[idx, ][[paste0(metric.name, ".passed.qc")]] <- qc.pass.cl
    passed.qc[idx] <- qc.pass.cl
  }

  return(df.qc)
}


.ddqcBoxplot <- function(df.qc, metric.name, h.line.x=0, do.log=False) {
  plt.data <- data.frame(metric=df.qc[[metric.name]], clusters=df.qc$cluster_labels)
  colnames(plt.data) <- c("metric", "clusters")

  if (do.log) {
    plt.data$metric <- log2(plt.data$metric)
    axis.labels <- labs(y=paste0("log2(", metric.name, ")"))
  } else {
    axis.labels <- labs(y=metric.name)
  }
  plt.data$clusters = with(plt.data, reorder(clusters, -metric, mean))

  horizontal_line <- NULL
  if (h.line.x > 0){
    horizontal_line <- geom_hline(yintercept=h.line.x, color="red", size=0.5)
  }

  boxplot <- ggplot(plt.data, aes(x=clusters, y=metric)) + geom_boxplot() + axis.labels + horizontal_line
  print(boxplot)
}


#' Perform initial quality control
#'
#' This function takes a Seurat object and does initial QC on it
#'
#' @param data Seurat object
#' @param basic.n.genes lower nFeature_RNA filtering. 100 by default
#' @param basic.percent.mt upper percent.mt filtering. 80 by default
#' @param mt.prefix gene regular expression used to calculate percent.mt in a cell
#' "MT-" by default
#' @param rb.prefix gene regular expression used to calculate percent.rb in a cell
#' @param basic.percent.mt upper percent.mt filtering. 80 by default
#' @return Filltered Seurat object
#' @export
initialQC <- function(data, basic.n.genes=100, basic.percent.mt=80, mt.prefix="MT-", rb.prefix="^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA") {
  data[["percent.mt"]] <- PercentageFeatureSet(data, features=grep(mt.prefix, rownames(data$RNA), ignore.case=TRUE))
  data[["percent.rb"]] <- PercentageFeatureSet(data, features=grep(rb.prefix, rownames(data$RNA), ignore.case=TRUE))

  data <- subset(data, subset = nFeature_RNA >= basic.n.genes & percent.mt <= basic.percent.mt)
  return(data)
}


#' Calculate which cells are passing ddqc
#'
#' This function takes a Seurat object after InitialQC, and then performs ddqc on it
#' Returns a data.frame that tells which cells have passed ddqc and additional information
#'
#' @param data Seurat object
#' @param res clustering resolution. 1 by default
#' @param threshold MAD multiplier for ddqc. 2 by default
#' @param do.counts whether to consider nCount_RNA for ddqc. TRUE by default
#' @param do.genes whether to consider nFeature_RNA for ddqc. TRUE by default
#' @param do.mito whether to consider percent.mt for ddqc. TRUE by default
#' @param do.ribo whether to consider percent.rb for ddqc. TRUE by default
#' @param n.genes.lower.bound bound for lower nFeature_RNA cluster-level threshold. 200 by default
#' @param percent.mito.upper.bound bound for upper percent.mt cluster-level threshold. 10 by default
#' @param random.state random seed for clustering results reproducibility. 29 by default
#'
#' @return data.frame with ddqc statistics
#' @export
ddqc.metrics <- function(data, res=1, threshold=2, do.counts=TRUE, do.genes=TRUE, do.mito=TRUE, do.ribo=FALSE,
                         n.genes.lower.bound=200, percent.mito.upper.bound=10, random.state=29) {
  data <- .clusterData(data, res=res, random.seed = random.state)

  df.qc <- data.frame("cluster_labels"=data$seurat_clusters, row.names=colnames(data))
  passed.qc <- vector(mode="logical", length=length(data$seurat_clusters))
  passed.qc <- TRUE

  if (do.counts) {
    df.qc <- .metricFilter(data, df.qc, threshold, "nCount_RNA", do.lower.co=TRUE)
    passed.qc <- passed.qc & df.qc$nCount_RNA.passed.qc
  }
  if (do.genes) {
    df.qc <- .metricFilter(data, df.qc, threshold, "nFeature_RNA", do.lower.co=TRUE,
                           lower.bound=n.genes.lower.bound)
    passed.qc <- passed.qc & df.qc$nFeature_RNA.passed.qc
  }
  if (do.mito) {
    df.qc <- .metricFilter(data, df.qc, threshold, "percent.mt", do.upper.co=TRUE,
                           upper.bound=percent.mito.upper.bound)
    passed.qc <- passed.qc & df.qc$percent.mt.passed.qc
  }
  if (do.ribo) {
    df.qc <- .metricFilter(data, df.qc, threshold, "percent.rb", do.upper.co=TRUE)
    passed.qc <- passed.qc & df.qc$percent.rb.passed.qc
  }
  df.qc[["passed.qc"]] <- passed.qc

  .ddqcBoxplot(df.qc, "nFeature_RNA", log2(200), TRUE)
  .ddqcBoxplot(df.qc, "percent.mt", 10, FALSE)

  return(df.qc)
}


#' Filter the Seurat object
#'
#' This function filters Seurat object based on df.qc
#'
#' @param data Seurat object
#' @param df.qc result of ddqc.metrics

#' @return Filtered Seurat object
#' @export
filterData <- function(data, df.qc) {
  data[["passed.qc"]] <- df.qc$passed.qc
  subset(data, subset = passed.qc)
  data[["passed.qc"]] <- NULL
  return(data)
}
