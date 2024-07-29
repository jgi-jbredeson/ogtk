#!/usr/bin/env Rscript

# require(viridis)

args = commandArgs(TRUE)

M.obs = as.matrix(read.table(args[[1]], stringsAsFactors=FALSE, row.names=1, header=TRUE))

r.sizes = apply(M.obs, 1, sum)
c.sizes = apply(M.obs, 2, sum)

r.exp = r.sizes / sum(r.sizes, na.rm=TRUE)
c.exp = c.sizes / sum(c.sizes, na.rm=TRUE)


M.exp = sum(M.obs) * (r.exp %*% t(c.exp))
M.ind = matrix(0, nrow(M.exp), ncol(M.exp))
M.ind[(M.obs >= 3) & (M.exp > 0) & (M.obs >= M.exp)] = 1

kwb = colorRampPalette(c("black","white","blue"))

M.diff = M.obs - M.exp
z.lim = 0.25 * max(abs(M.diff))
M.diff[M.diff < -z.lim] = -z.lim
M.diff[M.diff > +z.lim] = +z.lim

pdf(sprintf("%s.pdf", args[[2]]))
image(
  (seq(nrow(M.obs))-1)/(nrow(M.obs)-1),
  (seq(ncol(M.obs))-1)/(ncol(M.obs)-1),
  z=M.diff,
  zlim=z.lim*c(-1.0,1.0),
  axes=FALSE, xlab="", ylab="", col=kwb(length(unique(c(M.diff)))))
axis(1, at=(seq(nrow(M.obs))-1)/(nrow(M.obs)-1), labels=rownames(M.obs), cex.lab=0.5, las=2)
axis(2, at=(seq(ncol(M.obs))-1)/(ncol(M.obs)-1), labels=colnames(M.obs), cex.lab=0.5, las=2)


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
    if (M.ind[i,j] > 0) {
      points((i-1)/(nrow(M.obs)-1), (j-1)/(ncol(M.obs)-1), pch=8, cex=0.5, col="cyan2")  #blue4")
    }
    L.i = append(L.i, row.names[i])
    L.j = append(L.j, col.names[j])
    L.obs = append(L.obs, M.obs[i,j])
    L.exp = append(L.exp, M.exp[i,j])
    L.fit = append(L.fit, M.ind[i,j])
  }
}
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=NULL, border="black", lwd=1.5)
invisible(dev.off())

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
