aa <- get.seq("P21554")
blast <- blast.pdb(aa)
hits <- plot.blast(blast, cutoff = 50) #DEfault cutoff = 368, hits = 12 
#head(hits$pdb.id)
files <- get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE)
pdbs <- pdbaln(files, fit=TRUE, web.args=list(email="anmayroy@gmail.com"))
ids <- basename.pdb(pdbs$id)
plot(pdbs, labels=ids)

# Calculate sequence conservation
cons <- conserv(pdbs, method="entropy22")

# SSE annotations
sse <- pdbs2sse(pdbs, ind=1, rm.gaps=FALSE)

# Plot conservation per residue
plotb3(cons, sse=sse, ylab="Sequence entropy")

anno <- pdb.annotate(ids)
print(unique(anno$source))

# find invariant core
core <- core.find(pdbs)

# superimpose all structures to core
pdbs$xyz = pdbfit(pdbs, core)

# Perform PCA
pc.xray <- pca(pdbs)

# Calculate RMSD
rd <- rmsd(pdbs)

# Structure-based clustering
hc.rd <- hclust(dist(rd))
grps.rd <- cutree(hc.rd, k=3) 
plot(pc.xray, 1:2, col="grey50", bg=grps.rd, pch=21, cex=1)

# Visualize first principal component
mktrj(pc.xray, pc=1, file="pc_1.pdb")
library(ggplot2)
library(ggrepel) 

df <- data.frame(x=pc.xray$z[,1], y=pc.xray$z[,2])
col <- as.factor(grps.rd)

p <- ggplot(df, aes(x, y)) +
  geom_point(aes(col=col), size=2) +
  xlab("PC1") +
  ylab("PC2") +
  scale_color_discrete(name="Clusters") + 
  geom_text_repel(aes(label=ids))

p

