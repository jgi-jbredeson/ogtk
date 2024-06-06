#!/usr/bin/env Rscript

args = commandArgs(TRUE)

d = read.table(args[[1]], stringsAsFactors=FALSE, header=TRUE)
M.obs = as.matrix(d[,2:ncol(d)])
rownames(M.obs) = d[,1]


# r.sizes = read.table(args[[2]], stringsAsFactors=FALSE, header=FALSE)
# r.names = r.sizes[,1]
# r.sizes = r.sizes[,2]
# names(r.sizes) = r.names

# c.sizes = read.table(args[[3]], stringsAsFactors=FALSE, header=FALSE)
# c.names = c.sizes[,1]
# c.sizes = c.sizes[,2]
# names(c.sizes) = c.names

# # subset and reorder
# r.sizes = r.sizes[rownames(M.obs)]
# c.sizes = c.sizes[colnames(M.obs)]

r.sizes = apply(M.obs, 1, sum)
names(r.sizes) = rownames(M.obs)

c.sizes = apply(M.obs, 2, sum)
names(c.sizes) = colnames(M.obs)

r.exp = r.sizes / sum(r.sizes, na.rm=TRUE)
c.exp = c.sizes / sum(c.sizes, na.rm=TRUE)


M.exp = sum(M.obs) * (r.exp %*% t(c.exp))
M.fit = M.obs >= 3 & M.exp > 0 & M.obs >= M.exp

pdf(sprintf("%s.obs.pdf", args[[2]]))
image((seq(nrow(M.obs))-1)/(nrow(M.obs)-1), (seq(ncol(M.obs))-1)/(ncol(M.obs)-1), z=M.obs, axes=FALSE, xlab="", ylab="")
axis(1, at=(seq(nrow(M.obs))-1)/(nrow(M.obs)-1), labels=rownames(M.obs), cex.lab=0.5, las=2)
axis(2, at=(seq(ncol(M.obs))-1)/(ncol(M.obs)-1), labels=colnames(M.obs), cex.lab=0.5, las=2)
dev.off()

pdf(sprintf("%s.pass.pdf", args[[2]]))
image((seq(nrow(M.fit))-1)/(nrow(M.fit)-1), (seq(ncol(M.fit))-1)/(ncol(M.fit)-1), z=M.fit, axes=FALSE, xlab="", ylab="")
axis(1, at=(seq(nrow(M.fit))-1)/(nrow(M.fit)-1), labels=rownames(M.fit), cex.lab=0.5, las=2)
axis(2, at=(seq(ncol(M.fit))-1)/(ncol(M.fit)-1), labels=colnames(M.fit), cex.lab=0.5, las=2)
dev.off()

rownames(M.exp) = row.names = rownames(M.obs)
colnames(M.exp) = col.names = colnames(M.obs)

# write.table(M.exp, file=sprintf("%s.exp.matrix", args[[2]]), row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

L.i = c()
L.j = c()
L.obs = c()
L.exp = c()
L.fit = c()
for (i in 1:nrow(M.obs)) {
  for (j in 1:ncol(M.obs)) {
    if (is.na(M.obs[i,j])) {
      next
    }
    L.i = append(L.i, row.names[i])
    L.j = append(L.j, col.names[j])
    L.obs = append(L.obs, M.obs[i,j])
    L.exp = append(L.exp, M.exp[i,j])
    L.fit = append(L.fit, M.fit[i,j])
  }
}

write.table(
  data.frame(
    INDV1=L.i,
    INDV2=L.j,
    N_OBS=L.obs,
    N_EXP=L.exp,
    PASS=L.fit,
    stringsAsFactors=FALSE
  ),
  file=sprintf("%s.tsv", args[[2]]),
  row.names=FALSE,
  col.names=TRUE,
  sep="\t",
  quote=FALSE
)
