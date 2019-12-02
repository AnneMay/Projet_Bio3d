## Préparation
## Set up de l'environnement de travail 
#setwd("C:/Users/anmay/OneDrive - Universite de Montreal/Documents/DESS_Bio-INFO/BIF7105/Bio3d/script/")
## Importation des library
# install.packages("bio3d", dependencies=TRUE)
library(bio3d)
library(tidyverse)

##Référence: https://bioboot.github.io/bimm143_W18/class-material/lecture11-BIMM143_W18.pdf

# List of bio3d functions with brief description
help(package=bio3d)

pdb <- read.pdb("6MSM")
pdb

attributes(pdb)
head(pdb$atom)
# Note that individual $atom records can also be accessed like this
pdb$atom$elety[1:2] 

# Which allows us to do the following
plot.bio3d(pdb$atom$b[pdb$calpha], sse=pdb, typ="l", ylab="B-factor")
##Residue temperature factors for PDB ID 4q21 with secondary structure element (SSE) annotation 
##in marginal regions plotted with function plot.bio3d()

pdb$xyz 
dim(pdb$xyz) 
pdb$xyz[ 1, atom2xyz(1:2) ] 
ca.inds <- atom.select(pdb, "calpha")
ca.inds 
head( pdb$atom[ca.inds$atom, ] )
a.inds <- atom.select(pdb, chain="A")
ca.inds <- atom.select(pdb, "calpha", chain="A")
cab.inds <- atom.select(pdb, elety=c("CA","CB"), chain="A", resno=10:20)
b.inds <- atom.select(pdb, "back")
backpdb <- trim.pdb(pdb, b.inds)
write.pdb(backpdb, file="46MSM_back.pdb")

backpdb <- trim.pdb(pdb, "backbone")
backpdb <- atom.select(pdb, "backbone", value=TRUE)
write.pdb(backpdb, resno=backpdb$atom$resno+10)
write.pdb(backpdb, chain="B")

# Test de multiples structures 
blast <- blast.pdb(pdb)
hits <- plot.blast(blast) #188 hits, cutoff = 40
files <- get.pdb(hits, split=TRUE)
pdbs <- pdbaln(files, fit=TRUE, web.args=list(email="anmayroy@gmail.com"))
#Visualisation des donneés 
pdbs
pdbs$ali[1:5, 1000:1008] 
gaps <- gap.inspect(pdbs$ali)
head(gaps$f.inds)
pdbs$ali[, gaps$f.inds]

#PCA
seqID <- seqidentity(pdbs) 
rd <- rmsd(pdbs)
hc <- hclust(as.dist(rd))
grps <- cutree(hc, k=7)
hclustplot(hc, k=7)
pc <- pca(pdbs)
plot(pc, col=grps)

#NMA 
modes <- nma(pdbs)
plot(modes, pdbs, col=grps, spread=TRUE)


# Test Pratique 

aa <- get.seq("6MSM_A")
blast <- blast.pdb(aa)
hits <- plot.blast(blast, cutoff = 300) #DEfault cutoff = 368, hits = 12 
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
grps.rd <- cutree(hc.rd, k=4) 
plot(pc.xray, 1:2, col="grey50", bg=grps.rd, pch=21, cex=1)

# Visualize first principal component
mktrj(pc.xray, pc=1, file="pc_1.pdb")


library(ggplot2)
library(ggrepel) 

df <- data.frame(x=pc.xray$z[,1], y=pc.xray$z[,2])
col <- as.factor(grps)

p <- ggplot(df, aes(x, y)) +
  geom_point(aes(col=col), size=2) +
  xlab("PC1") +
  ylab("PC2") +
  scale_color_discrete(name="Clusters") + 
  geom_text_repel(aes(label=ids))

p ##On voit ici que la strucure 3GD7_A (et 6GJQ_A) cause la majorité de la variation. On devrait reprendre l'analyse en exculant cette valeur

pdbs2 <- trim(pdbs, row.inds = c(1:17, 19:20))
ids2 <- basename.pdb(pdbs2$id)

#PCA
seqID <- seqidentity(pdbs) 
rd <- rmsd(pdbs)
hc <- hclust(as.dist(rd))
grps <- cutree(hc, k=7)
hclustplot(hc, k=7)
pc <- pca(pdbs)
plot(pc, col=grps) ##On voit ici que la strucure 3GD7_A et 6GJQ_A cause la majorité de la variation. On devrait reprendre l'analyse en exculant cette valeur

#NMA 
modes <- nma(pdbs)
plot(modes, pdbs, col=grps, spread=TRUE)

