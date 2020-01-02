#################################
# @author: Mariano Alvarez
# @author: Alessandro Vasciaveo
# @copyright: 2019
# -------------------------------

source("sources/examples/als-tables-generation/viper-utility-functions.R")

viper_analysis <- list()
viper_analysis$data_dir <- file.path("data/examples/als-tables-generation")
viper_analysis$output_dir <- file.path("output/")

# Marina analysis
library(viper)
library(limma)
library(aanot)
library(atools)
load( file.path( viper_analysis$data_dir , "als_sod1_raw_expset.rda" ) )
annot <- t(sapply(strsplit(readLines( file.path( viper_analysis$data_dir , "samples.txt" )) , "\t" ), function(x) x))
groups <- annot[match(colnames(dset), annot[, 1]), 2:4]
groups <- cbind(groups, unlist(apply(groups, 1, paste, collapse="-"), use.names=F))
colnames(groups) <- c("treat", "time", "batch", "index")
groups <- as.data.frame(groups)

# Processing only the 72h time point
filtro <- groups$time %in% c(0, 72)
dset <- dset[, filtro]
groups <- groups[filtro, ]
set.seed(1)
d1 <- DEtransform(dset)
# d1 <- limma::voom(dset)$E

#Absolute expression
absexp <- rnaseqExpressionFilter(d1, cores=4)
absexp <- rowMeans(absexp)

# Computing the 72h signatures
treat <- factor(paste(groups$treat, groups$time, sep="."))
batch <- factor(groups$batch)
design <- model.matrix(~0+treat+batch)
colnames(design) <- gsub("treat", "", colnames(design))
fit <- lmFit(d1, design)
cont <- makeContrasts("72h-0h"=(Tg.72-Tg.0)-(WT.72-WT.0+NTg.72-NTg.0)/2, levels=design)
fit1 <- eBayes(contrasts.fit(fit, cont))
ss <- qnorm(fit1$p.value/2, lower.tail=F)[, 1]*sign(fit1$t)[, 1]

# membrane
load( file.path( viper_analysis$data_dir , "mus_brain_GSE10415_mem-regulon.rda" ) , verbose = TRUE ) 
data(desc, package="aanot")
load( file.path( viper_analysis$data_dir , "mouse-GO.rda" ) , verbose = TRUE)
pmem <- c(gomus[["GO:0005886"]], "94185")
regul <- pruneRegulon(regul[names(regul) %in% pmem], 100 )

res <- msviper(ss, regul, minsize=50 ) #, adaptive.size=F, ges.filter=T)
# res <- msviper(ss, regul, nullmodel = .ges.null.model , minsize=25, adaptive.size=F, ges.filter=T)
tmp <- cbind(GeneID=names(res$es$nes), Symbol=aanot::entrez2gene(names(res$es$nes)), NES=round(res$es$nes, 2), "p-value"=signif(res$es$p.value, 3), FDR=signif(p.adjust(res$es$p.value, "fdr"), 3), Bonferroni=signif(p.adjust(res$es$p.value, "bonferroni"), 3))[order(res$es$p.value), ]
# tmp <- cbind(GeneID=names(res$es$nes), Symbol=aanot::entrez2gene(names(res$es$nes)), NES=round(res$es$nes, 2), "p-value"=signif(res$es$p.value, 3), FDR=signif(p.adjust(res$es$p.value, "fdr"), 3), Bonferroni=signif(p.adjust(res$es$p.value, "bonferroni"), 3))[order(res$es$p.value), ]
tmp <- rbind(colnames(tmp), tmp)

active_receptors_filename <- file.path( viper_analysis$output_dir,"active-receptors.txt")
cat(unlist(apply(tmp, 1, paste, collapse="\t"), use.names=FALSE), sep="\n", file = active_receptors_filename )

# # Crossing receptors with ligands
# # ligands <- sapply(strsplit(readLines( file.path( viper_analysis$data_dir , "putative_ligands_mouse.txt") )[-1], "\t"), function(x) x[1])
# # ligands <- sapply(strsplit(readLines("sources/analysis-for-paper/new-data/putative_ligands_mouse_new_data_may2019.txt")[-1], ","), function(x) x[1])
# ligands <- sapply(strsplit(readLines(active_receptors_filename)[-1], ","), function(x) x[3]) # Get the ENTREZ ID
# viper_analysis$ligands <- ligands
# 
# load( file.path( viper_analysis$data_dir , "string_mus_from_homo.rda" ) , verbose = TRUE )
# cr <- string[match(names(res$es$nes)[order(res$es$p.value)], names(string))]
# cr <- cr[!is.na(names(cr))]
# cr <- lapply(cr, function(x, genes) {
#   pos <- which(x$geneid %in% genes)
#   if (length(pos)==0) return(NULL)
#   tmp <- x$score[pos]
#   names(tmp) <- x$geneid[pos]
#   return(tmp[order(tmp, decreasing=T)])
# }, genes=ligands)
# cr <- cr[sapply(cr, length)>0]
# names(cr) <- as.vector(aanot::entrez2gene(names(cr)))
# cr <- lapply(cr, function(x) {names(x) <- as.vector(aanot::entrez2gene(names(x))); x})
# cr1 <- lapply(cr, function(x) x[x>800])
# cr1 <- cr1[sapply(cr1, length)>0]
# 
# # Dump the whole database
# res1 <- msviperAnnot(res, aanot::entrez2gene(unique(c(rownames(res$signature), names(res$regulon)))))
# pos <- sapply(cr, length)
# tmp <- cbind(Receptor=rep(names(cr), pos), Activity=round(rep(res1$es$nes[match(names(cr), names(res1$es$nes))], pos), 3), "p-value"=rep(signif(p.adjust(res1$es$p.value, "bonferroni"), 3)[match(names(cr), names(res1$es$nes))], pos), Ligand=unlist(lapply(cr, names), use.names=F), "String score"=unlist(cr, use.names=F))
# cat(unlist(apply(rbind(colnames(tmp), tmp), 1, paste, collapse="\t"), use.names=F), sep="\n", file=file.path( viper_analysis$output_dir ,"_string_all_interactions.txt") )
# 
# tmp1 <- tmp[as.numeric(tmp[, 5])>850, ] # Applying filter from string
# # message(">>> Not Applying STRING filter")
# # tmp1 <- tmp # NOT  Applying filter from string
# 
# cat(unlist(apply(rbind(colnames(tmp1), tmp1), 1, paste, collapse="\t"), use.names=F), sep="\n", file=file.path( viper_analysis$data_dir ,"_string_reliable_interactions.txt"))
# 
# tmp2 <- tapply(1:nrow(tmp1), tmp1[, 1], function(i, tmp1) c(tmp1[i[1], 1:3], Ligands=paste(tmp1[i, 4], collapse=", ")), tmp1=tmp1)
# tmp2 <- t(sapply(tmp2, function(x) x))
# tmp2 <- tmp2[order(as.numeric(tmp2[, 3])), ]
# names(absexp) <- as.vector(aanot::entrez2gene(names(absexp)))
# pval <- as.numeric(tmp2[, 3])
# tmp2 <- cbind(tmp2[, 1:3], FDR=signif(p.adjust(pval, "fdr"), 2), Bonferroni=signif(p.adjust(pval, "bonferroni"), 2), Expression=round(absexp[match(tmp2[, 1], names(absexp))], 2), Ligands=tmp2[, 4])
# cat(unlist(apply(rbind(colnames(tmp2), tmp2), 1, paste, collapse="\t"), use.names=F), sep="\n", file= file.path( viper_analysis$data_dir ,"_receptors-manuscript.txt") )
# cat(unlist(apply(rbind(colnames(tmp2), tmp2), 1, paste, collapse="\t"), use.names=F), sep="\n", file=file.path( viper_analysis$data_dir ,"_string-filtered-receptors.txt") )
# 
# # pdf("data/for-paper/from-mariano/DR6-clean-latest/_marina-manuscript.pdf", w=9, h=7, useD=F, pointsize=18)
# # par(mai=c(.1, .1, 1, 1))
# # res2 <- res1
# # res2$es$p.value <- p.adjust(res2$es$p.value, "bonferroni")
# # #res2$regulon <- pruneRegulon(res2$regulon, 200)
# # plot(res2, tmp2[as.numeric(tmp2[, 3])<.01 & as.numeric(tmp2[, 4])>.8, 1], hybrid=T)
# # dev.off()
# # 
# # pdf("data/for-paper/from-mariano/DR6-clean-latest/_marina-manuscript01.pdf", w=9, h=3, useD=F, pointsize=18)
# # par(mai=c(.1, .1, 1, 1))
# # res2 <- res1
# # plot(res2, c("Gria2", "Tnfrsf21", "Itgb1", "Ptprc", "Cd74"), hybrid=T)
# # dev.off()
# 
