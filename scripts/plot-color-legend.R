#!/usr/bin/env Rscript

args = commandArgs(TRUE)
colors = read.table(args[[1]], stringsAsFactors=FALSE, header=FALSE)
colnames(colors) = c("label","hex")

pdf(args[[2]])
plot(c(0,nrow(colors)), c(0,nrow(colors)), col=0, xlab='', ylab='', main=NULL, axes=FALSE)
legend("topleft", as.character(colors$label), fill=as.character(colors$hex), ncol=3, bty='n', cex=2)
invisible(dev.off())

