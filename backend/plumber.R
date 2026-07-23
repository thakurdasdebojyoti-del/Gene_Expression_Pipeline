# plumber.R — Expression Console: ONE R server doing three jobs.
#
#   1. serves the frontend   (./www/index.html)   at  /
#   2. real authentication   (/auth/*)            — hashed passwords + session tokens
#   3. the analysis API      (/analyze, ...)      — limma differential expression
#
# Because the page and the API share an origin there is NO CORS problem, and the
# frontend needs no API URL configured.
#
# Local run:  Rscript start.R    ->  http://localhost:8000
#
# USER STORAGE: users live in USERS_FILE (JSON). Inside a container that path is
# EPHEMERAL — it resets on redeploy. Mount a volume at /data, or swap
# load_users()/store_users() for a database, before relying on it.

suppressPackageStartupMessages({
  library(plumber); library(GEOquery); library(limma); library(Biobase)
  library(jsonlite); library(sodium)
})
options(timeout = 600)

# ============================================================================ #
#  AUTH                                                                         #
# ============================================================================ #
USERS_FILE <- Sys.getenv("USERS_FILE", "/data/users.json")
dir.create(dirname(USERS_FILE), recursive = TRUE, showWarnings = FALSE)

load_users <- function(){
  if (!file.exists(USERS_FILE)) return(list())
  tryCatch(jsonlite::fromJSON(USERS_FILE, simplifyVector = FALSE), error = function(e) list())
}
store_users <- function(u){
  tryCatch(write(jsonlite::toJSON(u, auto_unbox = TRUE, pretty = TRUE), USERS_FILE),
           error = function(e) warning("could not persist users: ", conditionMessage(e)))
}

# passwords hashed with libsodium (scrypt) — plaintext is never stored
hash_pw   <- function(pw) sodium::password_store(pw)
verify_pw <- function(pw, hash) isTRUE(tryCatch(sodium::password_verify(hash, pw),
                                                error = function(e) FALSE))

# seed a demo account on first boot so the app is usable immediately
if (length(load_users()) == 0) {
  u <- list()
  u[["researcher@lab.edu"]] <- list(email="researcher@lab.edu", name="Researcher",
    rid="8824", plan="Premium Plan", hash=hash_pw("demo1234"), created=as.character(Sys.time()))
  store_users(u)
  message(">>> seeded demo account:  researcher@lab.edu  /  demo1234")
}

# in-memory session tokens (cleared on restart)
SESSIONS  <- new.env(parent = emptyenv())
TOKEN_TTL <- 60 * 60 * 12   # 12 hours

new_token <- function(email){
  tok <- paste0(sample(c(letters, LETTERS, 0:9), 44, replace = TRUE), collapse = "")
  assign(tok, list(email = email, expires = as.numeric(Sys.time()) + TOKEN_TTL), envir = SESSIONS)
  tok
}
token_user <- function(tok){
  if (is.null(tok) || !nzchar(tok) || !exists(tok, envir = SESSIONS, inherits = FALSE)) return(NULL)
  s <- get(tok, envir = SESSIONS)
  if (s$expires < as.numeric(Sys.time())) { rm(list = tok, envir = SESSIONS); return(NULL) }
  s$email
}
req_token <- function(req){
  h <- req$HTTP_AUTHORIZATION
  if (is.null(h)) return(NULL)
  sub("^[Bb]earer ", "", h)
}
require_auth <- function(req, res){
  em <- token_user(req_token(req))
  if (is.null(em)) { res$status <- 401; return(NULL) }
  em
}

# ============================================================================ #
#  ANALYSIS                                                                     #
# ============================================================================ #
REGISTRY <- list(
  GSE2034="Breast cancer (classic DE example)", GSE45827="Breast cancer subtypes",
  GSE10072="Lung cancer vs normal",             GSE53757="Kidney cancer (RNA-seq)",
  GSE14520="Liver cancer",                      GSE39582="Colorectal cancer",
  GSE48350="Alzheimer's disease",               GSE55235="Rheumatoid arthritis",
  GSE15641="Influenza infection",               GSE11121="Type 2 Diabetes"
)

row_var    <- function(x) rowSums((x - rowMeans(x))^2) / (ncol(x) - 1)
needs_log2 <- function(m){ q <- as.numeric(quantile(m, c(0,.25,.5,.75,.99,1), na.rm=TRUE))
  (q[5] > 100) || (q[6]-q[1] > 50 && q[2] > 0) }
finite0    <- function(x){ x[!is.finite(x)] <- 0; x }

# Derive a clean 2-level Control/Case factor from the series metadata.
# Fires only when EVERY sample lands unambiguously in one of two groups of >= 3.
auto_group <- function(pheno){
  ctrl_p <- "normal|control|healthy|non[- ]?tumou?r|adjacent|benign|non[- ]?relapse|no relapse|uninfected|mock"
  case_p <- "tumou?r|cancer|carcinoma|disease|patient|case|relapse|malignant|affected|infected|diabet"
  cols <- pheno[, vapply(pheno, function(x) is.character(x)||is.factor(x), logical(1)), drop=FALSE]
  for (col in names(cols)){
    v <- tolower(as.character(cols[[col]])); ctrl <- grepl(ctrl_p,v); case <- grepl(case_p,v) & !ctrl
    if (sum(ctrl)>=3 && sum(case)>=3 && sum(ctrl)+sum(case)==length(v)){
      g <- factor(ifelse(ctrl,"Control","Case"), levels=c("Control","Case"))
      attr(g,"col") <- col; return(g)
    }
  }
  NULL
}

run_pipeline <- function(gse, logfc = 1, padj = 0.05){
  gse  <- toupper(trimws(gse))
  eset <- getGEO(gse, GSEMatrix = TRUE, AnnotGPL = FALSE)[[1]]
  mat  <- exprs(eset); meta <- pData(eset)
  title <- as.character(experimentData(eset)@title); plat <- annotation(eset)

  if (nrow(mat) == 0 || all(is.na(mat)))
    return(list(gse=gse, title=title, platform=plat, status=
      "no series matrix (RNA-seq — counts are in supplementary files; use getGEOSuppFiles + DESeq2)"))

  mat <- mat[rowMeans(is.na(mat)) <= 0.10, , drop=FALSE]
  if (needs_log2(mat)){ mat[mat<=0] <- NaN; mat <- log2(mat) }
  mat <- mat[complete.cases(mat), , drop=FALSE]
  mat_norm <- normalizeQuantiles(mat)
  v <- row_var(mat_norm); mat_norm <- mat_norm[v > quantile(v,.25) & v > 0, , drop=FALSE]

  grp0    <- auto_group(meta)
  grouped <- !is.null(grp0)
  grp <- if (grouped) factor(grp0[match(colnames(mat_norm), rownames(meta))],
                             levels=c("Control","Case")) else factor(rep("sample", ncol(mat_norm)))

  pc  <- prcomp(t(mat_norm), scale.=TRUE)
  ve  <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  pca <- data.frame(sample=colnames(mat_norm), PC1=finite0(pc$x[,1]), PC2=finite0(pc$x[,2]),
                    group=as.character(grp), stringsAsFactors=FALSE)

  ord  <- order(grp)
  topv <- head(order(row_var(mat_norm), decreasing=TRUE), 40)
  M    <- mat_norm[topv, ord, drop=FALSE]
  Z    <- finite0(t(scale(t(M))))
  heat <- list(genes=rownames(M), samples=colnames(M), groups=as.character(grp)[ord],
               z=lapply(seq_len(nrow(Z)), function(i) as.numeric(Z[i,])))

  de_out <- list(); n_up <- 0; n_down <- 0; group_counts <- list()
  if (grouped){
    group_counts <- as.list(table(grp))
    design <- model.matrix(~ grp); colnames(design) <- c("Intercept","Case_vs_Control")
    fit <- eBayes(lmFit(mat_norm, design))
    de  <- topTable(fit, coef="Case_vs_Control", number=Inf, sort.by="P")
    sig <- subset(de, adj.P.Val < padj & abs(logFC) > logfc)
    n_up <- sum(sig$logFC > 0); n_down <- sum(sig$logFC < 0)
    top <- head(de[order(de$P.Value), ], 3000)
    de_out <- data.frame(gene=rownames(top), logFC=round(top$logFC,4),
                         p=signif(top$P.Value,4), adjP=signif(top$adj.P.Val,4),
                         baseMean=round(2^top$AveExpr,2), stat=round(top$t,3),
                         check.names=FALSE)
  }

  list(gse=gse, title=title, platform=plat,
       n_samples=ncol(mat_norm), n_genes=nrow(mat_norm),
       grouped=grouped, group_col=if(grouped) attr(grp0,"col") else "",
       groups=group_counts, n_up=n_up, n_down=n_down, status="ok", real=TRUE,
       de=de_out, pca=pca, pc_var=as.numeric(ve[1:2]), heat=heat)
}

# ============================================================================ #
#  ROUTES                                                                       #
# ============================================================================ #

#* Serve the frontend from ./www at the site root
#* @assets ./www /
list()

#* CORS — needed only if the page is hosted on a different origin
#* @filter cors
function(req, res){
  res$setHeader("Access-Control-Allow-Origin", "*")
  if (identical(req$REQUEST_METHOD, "OPTIONS")){
    res$setHeader("Access-Control-Allow-Methods", "GET,POST,OPTIONS")
    h <- req$HTTP_ACCESS_CONTROL_REQUEST_HEADERS
    res$setHeader("Access-Control-Allow-Headers",
                  if (is.null(h)) "Content-Type,Authorization" else h)
    res$status <- 200; return(list())
  }
  forward()
}

#* Health check (public)
#* @get /health
#* @serializer unboxedJSON
function(){ list(status="ok", service="expression-console", users=length(load_users())) }

#* Create an account
#* @post /auth/register
#* @serializer unboxedJSON
function(req, res, email="", password="", name=""){
  email <- tolower(trimws(email))
  if (!grepl("^[^@[:space:]]+@[^@[:space:]]+\\.[^@[:space:]]+$", email)){
    res$status <- 400; return(list(error="Enter a valid email address")) }
  if (nchar(password) < 8){
    res$status <- 400; return(list(error="Password must be at least 8 characters")) }
  u <- load_users()
  if (!is.null(u[[email]])){
    res$status <- 409; return(list(error="An account with that email already exists")) }
  rid <- sprintf("%04d", sample(1000:9999, 1))
  u[[email]] <- list(email=email, name=if(nzchar(name)) name else strsplit(email,"@")[[1]][1],
                     rid=rid, plan="Free Plan", hash=hash_pw(password),
                     created=as.character(Sys.time()))
  store_users(u)
  list(token=new_token(email),
       user=list(email=email, name=u[[email]]$name, rid=rid, plan="Free Plan"))
}

#* Sign in
#* @post /auth/login
#* @serializer unboxedJSON
function(req, res, email="", password=""){
  email <- tolower(trimws(email))
  rec <- load_users()[[email]]
  # same message whether or not the account exists — don't leak registered emails
  if (is.null(rec) || !verify_pw(password, rec$hash)){
    res$status <- 401; return(list(error="Incorrect email or password"))
  }
  list(token=new_token(email),
       user=list(email=rec$email, name=rec$name, rid=rec$rid, plan=rec$plan))
}

#* Validate a token / fetch the signed-in user
#* @get /auth/me
#* @serializer unboxedJSON
function(req, res){
  em <- require_auth(req, res); if (is.null(em)) return(list(error="Not signed in"))
  rec <- load_users()[[em]]
  list(user=list(email=rec$email, name=rec$name, rid=rec$rid, plan=rec$plan))
}

#* Sign out
#* @post /auth/logout
#* @serializer unboxedJSON
function(req){
  tok <- req_token(req)
  if (!is.null(tok) && exists(tok, envir=SESSIONS, inherits=FALSE)) rm(list=tok, envir=SESSIONS)
  list(ok=TRUE)
}

#* Curated dataset list (public)
#* @get /datasets
#* @serializer unboxedJSON
function(){ lapply(names(REGISTRY), function(k) list(gse=k, desc=REGISTRY[[k]])) }

#* Inspect a series' phenotype columns — for assigning groups manually
#* @get /inspect
#* @serializer unboxedJSON
function(req, res, gse=""){
  em <- require_auth(req, res); if (is.null(em)) return(list(error="Not signed in"))
  ph <- pData(getGEO(toupper(trimws(gse)), GSEMatrix=TRUE, AnnotGPL=FALSE)[[1]])
  keep <- vapply(ph, function(x){ n <- length(unique(x)); n>1 && n<=12 }, logical(1))
  lapply(ph[, keep, drop=FALSE], function(x) as.list(table(x)))
}

#* Run the differential expression pipeline (sign-in required)
#* @post /analyze
#* @serializer unboxedJSON
function(req, res, gse="", logfc=1, padj=0.05){
  em <- require_auth(req, res); if (is.null(em)) return(list(error="Not signed in"))
  tryCatch(run_pipeline(gse, as.numeric(logfc), as.numeric(padj)),
           error=function(e) list(gse=toupper(trimws(gse)),
                                  status=paste("ERROR:", conditionMessage(e))))
}
