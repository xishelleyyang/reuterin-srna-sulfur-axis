#!/usr/bin/env Rscript
# Numerical engine: native limma/edgeR; controlled by the Python run_all entry point.
# No synthetic biological data. Seeds affect rotation/permutation tests only.
Sys.setlocale("LC_ALL", "C") # GenePix legacy headers contain non-UTF-8 bytes.
suppressPackageStartupMessages({library(limma); library(edgeR)})
a <- commandArgs(TRUE)
stopifnot(length(a) >= 2)
DATA <- normalizePath(a[1]); OUT <- a[2]
pilot <- length(a) > 2 && a[3] == "--pilot"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
OUT <- normalizePath(OUT)
NROT <- if(pilot) 199L else 9999L
SEEDS <- if(pilot) 42L else c(42L, 20260912L)
write_out <- function(x, name, row.names=FALSE) write.csv(x, file.path(OUT,name), row.names=row.names, na="")
ft <- read.delim(gzfile(file.path(DATA,"mg1655_feature_table.txt.gz")),check.names=FALSE)
ft <- ft[ft[[1]] == "gene" & !is.na(ft$symbol) & ft$symbol != "",]
b2s <- setNames(ft$symbol, ft$locus_tag)
SETS <- list(
 "ISC operon"=c("iscR","iscS","iscU","iscA","hscB","hscA","fdx"),
 "SUF operon"=c("sufA","sufB","sufC","sufD","sufS","sufE"),
 "RyhB targets"=c("sodB","fumA","acnA","sdhC","sdhD","sdhA","sdhB","bfr","ftnA","iscS","iscU","iscA"),
 "Fur iron uptake"=c("fhuA","fhuE","fhuF","fecA","fepA","fepB","entC","entB","entE","fes","cirA","fiu","tonB","exbB","exbD"),
 "OxyR regulon"=c("ahpC","ahpF","katG","dps","grxA","trxC","gor","fur","sufA","sufB","sufC","sufD","sufS","sufE"),
 "Cysteine regulon"=c("cysK","cysM","cysE","cysA","cysW","cysU","cysJ","cysI","cysH","cysD","cysN","cysC","cysP","sbp"),
 "Glutathione/glyoxalase"=c("gshA","gshB","gor","ggt","gloA","gloB","gloC","yajG","grxA","grxB"),
 "Electrophile detoxification"=c("nemA","nemR","yqhD","yqhC","dkgA","dkgB"))
HPX_DELETED <- c("katE","katG","ahpC")
stopifnot(all(unique(unlist(SETS)) %in% ft$symbol))

make_de <- function(fit, bnum, symbols, coef=1) {
 fit <- eBayes(fit)
 tt <- topTable(fit,coef=coef,number=Inf,sort.by="none")
 data.frame(feature_id=rownames(tt), bnum=bnum, symbol=symbols,
   logFC=tt$logFC,AveExpr=tt$AveExpr,t=tt$t,P.Value=tt$P.Value,
   adj.P.Val=tt$adj.P.Val,B=tt$B,
   moderated_SE=fit$stdev.unscaled[,coef]*sqrt(fit$s2.post),
   residual_df=fit$df.residual, stringsAsFactors=FALSE)
}
# -------- Microarray: identical upstream processing and native avereps MAList semantics.
soft <- readLines(gzfile(file.path(DATA,"GPL7445_family.soft.gz")))
hdr <- grep("^ID\t",soft)[1]; en <- grep("^!platform_table_end",tolower(soft))[1]
plat <- read.delim(textConnection(soft[hdr:(en-1)]),comment.char="#",check.names=FALSE)
bcol <- names(plat)[vapply(plat,function(x) sum(grepl("^b[0-9]{4}$",x))>1000,logical(1))][1]
id2b <- setNames(plat[[bcol]],plat$ID)
files <- file.path(DATA,paste0("GSM493",569:573,".gpr.gz"))
rg <- read.maimages(files,source="genepix",columns=list(R="F635 Median",G="F532 Median",Rb="B635 Median",Gb="B532 Median"),annotation=c("ID","Name"),wt.fun=function(x)as.numeric(x$Flags>=0))
ma <- normalizeWithinArrays(backgroundCorrect(rg,method="normexp",offset=50),method="loess")
ma$M[,2] <- -ma$M[,2]
bn <- unname(id2b[ma$genes$ID]); ok <- !is.na(bn) & bn != ""
M <- avereps(ma$M[ok,],ID=bn[ok]); A <- avereps(ma$A[ok,],ID=bn[ok])
W <- avereps(ma$weights[ok,],ID=bn[ok]); stopifnot(identical(rownames(M),rownames(W)))
colnames(M) <- colnames(A) <- colnames(W) <- paste0("rep",1:5)
sym <- unname(b2s[rownames(M)])
# Explicit sensitivity: exclude bad-flag probes from the within-gene mean.
masked_probe <- ma$M[ok,]; masked_probe[ma$weights[ok,] == 0] <- NA_real_
Mm <- avereps(masked_probe,ID=bn[ok]); colnames(Mm) <- colnames(M)
primary_fit <- lmFit(M,design=matrix(1,5,1),weights=W)
DE <- list(Reuterin=make_de(primary_fit,rownames(M),sym))
DE$Reuterin$AveExpr <- rowMeans(A)
DE$Reuterin$n_arrays_positive_weight <- rowSums(W>0)
DE$Reuterin$test_status <- ifelse(is.finite(DE$Reuterin$P.Value),"tested","not_estimable_zero_weight")
write_out(DE$Reuterin,"DE_Reuterin.csv")
write_out(data.frame(bnum=rownames(M),symbol=sym,M,check.names=FALSE),"M_Reuterin.csv")
write_out(data.frame(bnum=rownames(W),W,check.names=FALSE),"W_Reuterin.csv")
alt <- list(
 unweighted=make_de(lmFit(M,design=matrix(1,5,1)),rownames(M),sym),
 masked_probe_mean=make_de(lmFit(Mm,design=matrix(1,5,1),weights=W),rownames(Mm),sym))
for(nm in names(alt)) write_out(alt[[nm]],paste0("DE_Reuterin_",nm,".csv"))
# Native roast requires strictly positive precision weights. Do not replace zero
# quality weights by arbitrary epsilon values. Primary rotation tests use genes
# measurable with positive weight in all five arrays; DE retains partial cases.
valid_partial <- !is.na(sym) & rowSums(W>0)>=2 & rowSums(!is.finite(M))==0
valid <- !is.na(sym) & rowSums(W>0)==5 & rowSums(!is.finite(M))==0
stopifnot(!anyDuplicated(sym[valid]))
write_out(data.frame(bnum=rownames(M),symbol=sym,eligible_rotation_primary=valid,
 eligible_unweighted_sensitivity=valid_partial,n_positive_arrays=rowSums(W>0)),"array_rotation_eligibility.csv")
models <- list(Reuterin=list(X=M[valid,],W=W[valid,],symbols=sym[valid],design=matrix(1,5,1),coef=1,trend=FALSE))
Mm2 <- Mm[valid,]; Mm2[!is.finite(Mm2)] <- 0 # zero-weight placeholder only, not an observed value
sensitivity_models <- list(
 array_unweighted=list(Reuterin=list(X=M[valid_partial,],W=NULL,symbols=sym[valid_partial],design=matrix(1,5,1),coef=1,trend=FALSE)),
 array_masked_probe=list(Reuterin=list(X=Mm2,W=W[valid,],symbols=sym[valid],design=matrix(1,5,1),coef=1,trend=FALSE)),
 rnaseq_trend=list())

# -------- RNA-seq: frozen contrasts, TMM/voom primary, trend and edgeR QL sensitivity.
d <- read.delim(file.path(DATA,"GSE126176_counts_genes_formatted_metadata.txt.gz"),check.names=FALSE)
annot <- unique(d[,c("Gene","Transcript","GeneName")]); t2g <- setNames(annot$Gene,annot$Transcript)
meta_rows <- list(); count_summaries <- list(); blocked_models <- list()
run_rnaseq <- function(X,group,bnum,label,blocks=NULL) {
 stopifnot(all(is.finite(X)),all(X>=0),all(abs(X-round(X))<1e-7))
 group <- factor(group,levels=unique(c(if(label=="Hpx") "WT" else "W",group)))
 design <- model.matrix(~group)
 y <- DGEList(counts=X,group=group)
 keep <- filterByExpr(y,design=design)
 y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method="TMM")
 v <- voom(y,design,plot=FALSE)
 fit <- lmFit(v,design)
 bb <- bnum[keep]; ss <- unname(b2s[bb])
 de <- make_de(fit,bb,ss,coef=2)
 de$excluded_from_geneset <- if(label=="Hpx") ss %in% HPX_DELETED else FALSE
 de$test_status <- "tested"
 write_out(de,paste0("DE_",label,".csv"))
 logcpm <- cpm(y,log=TRUE,prior.count=0.5)
 trfit <- eBayes(lmFit(logcpm,design),trend=TRUE)
 tr <- topTable(trfit,coef=2,number=Inf,sort.by="none")
 write_out(data.frame(feature_id=rownames(tr),bnum=bb,symbol=ss,tr),paste0("DE_",label,"_trend.csv"))
 yy <- estimateDisp(y,design,robust=TRUE)
 ql <- glmQLFTest(glmQLFit(yy,design,robust=TRUE),coef=2)
 qt <- topTags(ql,n=Inf,sort.by="none")$table
 write_out(data.frame(feature_id=rownames(qt),bnum=bb,symbol=ss,qt),paste0("DE_",label,"_edgeR_QL.csv"))
 use <- !is.na(ss) & ss != ""
 stopifnot(!anyDuplicated(ss[use]))
 prim <- list(X=v$E[use,],W=v$weights[use,],symbols=ss[use],design=design,coef=2,trend=FALSE)
 trm <- list(X=logcpm[use,],W=NULL,symbols=ss[use],design=design,coef=2,trend=TRUE)
 if(!is.null(blocks)) {
   block <- factor(blocks)
   desb <- model.matrix(~block+group)
   vb <- voom(y,desb,plot=FALSE)
   db <- make_de(lmFit(vb,desb),bb,ss,coef=ncol(desb))
   write_out(db,paste0("DE_",label,"_blocked.csv"))
   blocked_models[[label]] <<- list(X=vb$E[use,],W=vb$weights[use,],symbols=ss[use],design=desb,coef=ncol(desb),trend=FALSE)
 }
 meta_rows[[label]] <<- data.frame(contrast=label,sample=colnames(X),group=as.character(group),block=if(is.null(blocks)) "not_used" else blocks,raw_library_size=colSums(X),filtered_library_size=y$samples$lib.size,TMM_factor=y$samples$norm.factors)
 count_summaries[[label]] <<- data.frame(contrast=label,n_input_features=nrow(X),n_filtered_features=nrow(y),n_mapped_features=sum(use),n_unmapped_features=sum(!use),n_biological_samples=ncol(X))
 list(de=de,primary=prim,trend=trm)
}
cc <- as.matrix(d[,grepl("-aft$",names(d))]); rownames(cc) <- d$Gene
for(pair in list(c("H","HOCl"),c("F","Ferrate"))) {
 sel <- substr(colnames(cc),2,2) %in% c("W",pair[1])
 rr <- run_rnaseq(cc[,sel],substr(colnames(cc)[sel],2,2),d$Gene,pair[2],substr(colnames(cc)[sel],1,1))
 DE[[pair[2]]] <- rr$de; models[[pair[2]]] <- rr$primary; sensitivity_models$rnaseq_trend[[pair[2]]] <- rr$trend
}
hfiles <- sort(list.files(DATA,pattern="^GSM50.*8h.*txt.gz$",full.names=TRUE))
hfiles <- hfiles[grepl("minus",hfiles)]
stopifnot(length(hfiles)==4)
hd <- lapply(hfiles,function(p){z<-read.delim(gzfile(p),header=FALSE,col.names=c("tx","count"));z$tx<-sub("transcript:","",z$tx);aggregate(count~tx,z,sum)})
stopifnot(all(vapply(hd,function(x)identical(x$tx,hd[[1]]$tx),logical(1))))
HX <- do.call(cbind,lapply(hd,`[[`,"count")); rownames(HX)<-hd[[1]]$tx; colnames(HX)<-basename(hfiles)
rr <- run_rnaseq(HX,ifelse(grepl("HPX",colnames(HX)),"HPX","WT"),unname(t2g[rownames(HX)]),"Hpx")
DE$Hpx<-rr$de;models$Hpx<-rr$primary;sensitivity_models$rnaseq_trend$Hpx<-rr$trend
write_out(do.call(rbind,meta_rows),"sample_normalization.csv")
write_out(do.call(rbind,count_summaries),"rnaseq_feature_accounting.csv")

# -------- Matched model, curated membership, full 32-test BH family, two seeds.
score_sets <- function(ms,seed,setnames=names(SETS)) {
 rows<-list()
 for(ct in names(ms)) {
   z<-ms[[ct]]; rownames(z$X)<-z$symbols
   if(!is.null(z$W)) rownames(z$W)<-z$symbols
   exc<-if(ct=="Hpx") HPX_DELETED else character(0)
   # "Measured" = has an individually testable DE call (Table S4 basis), independent of
   # whether it also clears the stricter roast() complete-positive-weight requirement.
   measured_ct<-unique(DE[[ct]]$symbol[DE[[ct]]$test_status=="tested" & !is.na(DE[[ct]]$symbol)])
   for(sn in setnames) {
     intended<-SETS[[sn]]; base<-setdiff(intended,exc)
     used<-intersect(base,z$symbols)
     not_measured<-setdiff(base,measured_ct)
     roast_excluded<-setdiff(setdiff(base,not_measured),used)
     stopifnot(length(used)>=3)
     # Stable per-contrast/set seed: changing unrelated membership cannot shift RNG stream.
     testseed<-seed+match(ct,c("Reuterin","HOCl","Ferrate","Hpx"))*100L+match(sn,names(SETS))
     set.seed(testseed)
     r<-roast(z$X,index=used,design=z$design,contrast=z$coef,weights=z$W,trend=z$trend,nrot=NROT,set.statistic="mean")
     df<-DE[[ct]]; fc<-setNames(df$logFC,df$symbol)
     # Use model's own fit for mean logFC so sensitivity rows are not mixed with primary.
     fitz<-lmFit(z$X,design=z$design,weights=z$W)
     fcz<-setNames(fitz$coefficients[,z$coef],z$symbols)
     rows[[length(rows)+1]]<-data.frame(contrast=ct,gene_set=sn,n_genes=length(used),n_samples=ncol(z$X),
       raw_mean_logFC=mean(fcz[used]),roast_p_up=r$p.value["Up","P.Value"],roast_p_mixed=r$p.value["Mixed","P.Value"],
       intended_genes=paste(intended,collapse=";"),used_genes=paste(used,collapse=";"),
       genotype_excluded_genes=paste(intersect(intended,exc),collapse=";"),
       not_measured_genes=paste(not_measured,collapse=";"),
       roast_excluded_partial_weight_genes=paste(roast_excluded,collapse=";"),seed=seed,nrot=NROT)
   }
 }
 ans<-do.call(rbind,rows)
 ans$roast_q_up<-p.adjust(ans$roast_p_up,"BH")
 ans$roast_q_mixed<-p.adjust(ans$roast_p_mixed,"BH")
 ans$BH_family_size<-nrow(ans)
 ans
}
primary<-score_sets(models,SEEDS[1]);write_out(primary,"genesets_primary.csv")
if(length(SEEDS)>1) write_out(score_sets(models,SEEDS[2]),"genesets_seed_check.csv")
# Each sensitivity gets its own full 32-test family; never mix selected p-values.
for(nm in names(sensitivity_models)) {
 mm<-models;mm[names(sensitivity_models[[nm]])]<-sensitivity_models[[nm]]
 write_out(score_sets(mm,42L),paste0("genesets_sensitivity_",nm,".csv"))
}
set3<-c("ISC operon","SUF operon","Electrophile detoxification")
bp<-score_sets(blocked_models,42L,set3)
up<-score_sets(models[c("HOCl","Ferrate")],42L,set3)
bp$model<-"replicate_block_voom";up$model<-"unblocked_voom"
write_out(rbind(up,bp),"batch_sensitivity.csv")
counts<-do.call(rbind,lapply(names(DE),function(nm){z<-DE[[nm]];data.frame(contrast=nm,n_rows=nrow(z),n_tested=sum(is.finite(z$P.Value)),FDR05=sum(z$adj.P.Val<0.05,na.rm=TRUE),FDR10=sum(z$adj.P.Val<0.1,na.rm=TRUE))}))
write_out(counts,"de_summary.csv")
write_out(data.frame(symbol=names(b2s),canonical_symbol=unname(b2s)),"reference_symbol_map.csv")
writeLines(capture.output(sessionInfo()),file.path(OUT,"R_sessionInfo.txt"))
saveRDS(list(sets=SETS,models=models,DE=DE,metadata=meta_rows),file.path(OUT,"analysis_state.rds"))
print(counts)
print(primary[,c("contrast","gene_set","n_genes","raw_mean_logFC","roast_q_up")])
cat("PILOT:",pilot,"NROT:",NROT,"Output:",OUT,"\n")
