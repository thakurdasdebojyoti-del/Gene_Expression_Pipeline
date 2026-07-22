# plumber.R — GEO differential-expression pipeline exposed as a REST API.
# Run:  Rscript -e "plumber::pr_run(plumber::plumb('plumber.R'), host='0.0.0.0', port=8000)"
# This is the BACKEND. It must run on a host that runs R (container on AWS/Render/
# Railway/Fly) — NOT on Vercel. The frontend (index.html) calls these endpoints.

suppressPackageStartupMessages({
  library(plumber); library(GEOquery); library(limma); library(Biobase)
})
options(timeout = 600)

REGISTRY <- list(
  GSE2034="Breast cancer (classic DE example)", GSE45827="Breast cancer subtypes",
  GSE10072="Lung cancer vs normal",             GSE53757="Kidney cancer (RNA-seq)",
  GSE14520="Liver cancer",                       GSE39582="Colorectal cancer",
  GSE48350="Alzheimer's disease",                GSE55235="Rheumatoid arthritis",
  GSE15641="Influenza infection",                GSE11121="Type 2 Diabetes"
)

# ---- helpers ---------------------------------------------------------------
row_var   <- function(x) rowSums((x - rowMeans(x))^2) / (ncol(x) - 1)
needs_log2 <- function(m){ q <- as.numeric(quantile(m,c(0,.25,.5,.75,.99,1),na.rm=TRUE))
  (q[5] > 100) || (q[6]-q[1] > 50 && q[2] > 0) }

auto_group <- function(pheno){
  ctrl_p <- "normal|control|healthy|non[- ]?tumou?r|adjacent|benign|non[- ]?relapse|no relapse|uninfected|mock"
  case_p <- "tumou?r|cancer|carcinoma|disease|patient|case|relapse|malignant|affected|infected|diabet"
  cols <- pheno[, vapply(pheno, function(x) is.character(x)||is.factor(x), logical(1)), drop=FALSE]
  for (col in names(cols)){
    v <- tolower(as.character(cols[[col]])); ctrl <- grepl(ctrl_p,v); case <- grepl(case_p,v) & !ctrl
    if (sum(ctrl)>=3 && sum(case)>=3 && sum(ctrl)+sum(case)==length(v)){
      g <- factor(ifelse(ctrl,"Control","Case"), levels=c("Control","Case")); attr(g,"col") <- col; return(g)
    }
  }
  NULL
}
finite0 <- function(x){ x[!is.finite(x)] <- 0; x }

# ---- core pipeline (returns a list plumber serializes to JSON) --------------
run_pipeline <- function(gse, logfc=1, padj=0.05){
  gse <- toupper(trimws(gse))
  gset <- getGEO(gse, GSEMatrix=TRUE, AnnotGPL=FALSE)
  eset <- gset[[1]]
  mat  <- exprs(eset); meta <- pData(eset)
  title <- as.character(experimentData(eset)@title)
  plat  <- annotation(eset)

  if (nrow(mat)==0 || all(is.na(mat)))
    return(list(gse=gse, title=title, platform=plat, status=
      "no series matrix (RNA-seq — counts are in supplementary files; use getGEOSuppFiles + DESeq2)"))

  mat <- mat[rowMeans(is.na(mat)) <= 0.10, , drop=FALSE]
  if (needs_log2(mat)){ mat[mat<=0] <- NaN; mat <- log2(mat) }
  mat <- mat[complete.cases(mat), , drop=FALSE]
  mat_norm <- normalizeQuantiles(mat)
  v <- row_var(mat_norm); mat_norm <- mat_norm[v > quantile(v,.25) & v > 0, , drop=FALSE]

  grp0 <- auto_group(meta)
  grouped <- !is.null(grp0)
  grp <- if (grouped) factor(grp0[match(colnames(mat_norm), rownames(meta))],
                             levels=c("Control","Case")) else factor(rep("sample", ncol(mat_norm)))

  # PCA (all samples) -------------------------------------------------------
  pc  <- prcomp(t(mat_norm), scale.=TRUE)
  ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  pca <- data.frame(sample=colnames(mat_norm),
                    PC1=finite0(pc$x[,1]), PC2=finite0(pc$x[,2]),
                    group=as.character(grp), stringsAsFactors=FALSE)

  # heatmap: top-40 variable genes, row z, columns ordered by group ---------
  ord   <- order(grp)
  topv  <- head(order(row_var(mat_norm), decreasing=TRUE), 40)
  M     <- mat_norm[topv, ord, drop=FALSE]
  Z     <- finite0(t(scale(t(M))))
  heat  <- list(genes=rownames(M), groups=as.character(grp)[ord],
                z=lapply(seq_len(nrow(Z)), function(i) as.numeric(Z[i,])))

  # differential expression (only if groups derived) ------------------------
  de_out <- list(); n_up <- 0; n_down <- 0; group_counts <- list()
  if (grouped){
    group_counts <- as.list(table(grp))
    design <- model.matrix(~ grp); colnames(design) <- c("Intercept","Case_vs_Control")
    fit <- eBayes(lmFit(mat_norm, design))
    de  <- topTable(fit, coef="Case_vs_Control", number=Inf, sort.by="P")
    sig <- subset(de, adj.P.Val < padj & abs(logFC) > logfc)
    n_up <- sum(sig$logFC > 0); n_down <- sum(sig$logFC < 0)
    top <- head(de[order(de$P.Value), ], 2000)          # cap payload; full run via batch script
    de_out <- data.frame(gene=rownames(top), logFC=round(top$logFC,4),
                         `adj.P.Val`=signif(top$adj.P.Val,4), check.names=FALSE)
  }

  list(gse=gse, title=title, platform=plat,
       n_samples=ncol(mat_norm), n_genes=nrow(mat_norm),
       grouped=grouped, group_col=if(grouped) attr(grp0,"col") else "",
       group_counts=group_counts, n_up=n_up, n_down=n_down, status="ok",
       de=de_out, pca=pca, pc_var=as.numeric(ve[1:2]), heatmap=heat)
}

# ---- CORS (so the Vercel frontend can call this cross-origin) --------------
#* @filter cors
function(req, res){
  res$setHeader("Access-Control-Allow-Origin", "*")
  if (identical(req$REQUEST_METHOD, "OPTIONS")){
    res$setHeader("Access-Control-Allow-Methods", "GET,POST,OPTIONS")
    h <- req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS; if (is.null(h)) h <- "Content-Type"
    res$setHeader("Access-Control-Allow-Headers", h)
    res$status <- 200; return(list())
  }
  forward()
}

#* Health check
#* @get /health
function(){ list(status="ok", service="geo-de-api", datasets=length(REGISTRY)) }

#* Preset series list (feeds the frontend dropdown)
#* @get /datasets
function(){ lapply(names(REGISTRY), function(k) list(gse=k, desc=REGISTRY[[k]])) }

#* Inspect a series' phenotype columns (to set groups manually)
#* @get /inspect
function(gse){
  eset <- getGEO(toupper(trimws(gse)), GSEMatrix=TRUE, AnnotGPL=FALSE)[[1]]
  ph <- pData(eset)
  keep <- vapply(ph, function(x){ u <- length(unique(x)); u>1 && u<=12 }, logical(1))
  lapply(ph[, keep, drop=FALSE], function(x) as.list(table(x)))
}

#* Run the full pipeline for one series
#* @post /analyze
#* @serializer unboxedJSON
function(gse="", logfc=1, padj=0.05){
  tryCatch(
    run_pipeline(gse, as.numeric(logfc), as.numeric(padj)),
    error = function(e) list(gse=toupper(trimws(gse)), status=paste("ERROR:", conditionMessage(e)))
  )
}
