#!/usr/bin/env Rscript

args = commandArgs(TRUE)

magnitude.scale.int = function(x) {
  s = 1
  if (x >= 1e9 - 1) {
    s = 1e9
  } else if (x >= 1e6 - 1) {
    s = 1e6
  } else if (x > 1e3 - 1) {
    s = 1e3
  } else if (x > 1e2 - 1) {
    s = 1e2
  } else if (x > 1e1 - 1) {
    s = 1e1
  }
  return(s)
}

magnitude.scale.str = function(x) {
  s = ""
  if (x >= 1e9 - 1) {
    s = "G"
  } else if (x >= 1e6 - 1) {
    s = "M"
  } else if (x > 1e3 - 1) {
    s = "k"
  }
  return(s)
}

if (length(args) < 3) {
    cat("plot-genomic-intervals.R <chr.sizes> <out.pdf> <file> [...]\n", file=stderr())
    quit(status=1)
}

chrom.sizes = read.table(args[[1]], stringsAsFactors=FALSE, header=FALSE)
chrom.names = chrom.sizes[,1]
chrom.sizes = chrom.sizes[,2]

names(chrom.sizes) = chrom.names

max.x = max(chrom.sizes)
max.y = length(chrom.names) * 0.5 - 0.5
mag.x = magnitude.scale.int(max.x)
lab.x = magnitude.scale.str(max.x)

chrom.offsets = rev(seq(0, max.y, 0.5))

names(chrom.offsets) = chrom.names


pdf(args[[2]]);
par(mar=c(5.1,5.1,1.1,1.1))
plot(c(0,max.x),c(0,max.y), xlab='Locus index', ylab='', main="", axes=FALSE, col=0, cex.lab=1.5)
axis(1)  #, at=seq(0, max.x, 0.5*mag.x), labels=seq(0, max.x, 0.5*mag.x)/mag.x)
axis(2, at=chrom.offsets, labels=chrom.names, las=2)

segments(0, chrom.offsets, chrom.sizes, chrom.offsets, lwd=2, col="black", lend=0)

sign = c(1,-1)

jitr = list()
for (i in 3:length(args)) {
    features.args = strsplit(args[[i]], ":")[[1]];

    feature.table = read.table(args[[i]], stringsAsFactors=FALSE, header=FALSE)
    feature.table = feature.table[feature.table[,1] %in% chrom.names,]

    features = as.character(feature.table[,ncol(feature.table)])
    for (feature in unique(features)) {
        if (! (feature %in% jitr)) {
            jitr[[feature]] = jitter(0)
        }
    }
    
    segments(
        feature.table[,3],
	chrom.offsets[feature.table[,1]] - 0.125 + unlist(jitr[features]),
	feature.table[,3],
	chrom.offsets[feature.table[,1]] + 0.125 + unlist(jitr[features]),
	col=as.character(features),
	lwd=0.5
    )
}
invisible(dev.off())

