### Référence : Tutoriel sur le NMA 
#### http://thegrantlab.org/bio3d/tutorials/normal-mode-analysis

##Set-up
library(bio3d)

## Protéine d'intérêt: 6MSM 
pdb <- read.pdb("6MsM") #Modèle: 9703 atomes, 2 chaines (A et B); Protéine: 9466 atomes

# Example 1: Basic Normal Mode Analysis
## Example 1A: Normal mode calculation

# Calcul des modes de la protéine 
modes <- nma(pdb)
# Somaire des modes (nombre total de modes = 3N où N est le nombre d'atomes de la structure)
## 6 premiers modes = mouvement de structure rigide avec une valeur de 0 
print(modes)
# Aperçu rapide des résultats du NMA 
plot(modes, sse=pdb)

## Example 1B: Specifying a force field

# Calculate modes with various force fields
modes.a <- nma(pdb, ff="calpha")
modes.b <- nma(pdb, ff="anm")
modes.c <- nma(pdb, ff="pfanm")
modes.d <- nma(pdb, ff="reach")
modes.e <- nma(pdb, ff="sdenm")

# Root mean square inner product (RMSIP)
rab <- rmsip(modes.a, modes.b)
rac <- rmsip(modes.a, modes.c)
plot(rab, xlab="ANM", ylab="C-alpha FF") 
plot(rac, xlab="PFANM", ylab="C-alpha FF") 
## Analysis of mode similarity between modes obtained from the *ANM* and *calpha* force fields 
## by calculating mode overlap and root mean square inner product (RMSIP) with function **rmsip()**. 
## An RMSIP value of *1* depicts identical directionality of the two mode subspaces.
### Permet de faire des comparaison entre les différents champs de forces utilisé pour le NMA, 
### selon la question de recherche on pourrait préférer un type par rapport au autres. 

## Example 1C: Normal mode analysis of the GroEL subunit (Autre sous-unité dans notre cas?)

# Download PDB, calcualte normal modes of the open subunit
pdb.full   <- read.pdb("6MSM")
pdb.open   <- trim.pdb(pdb.full, atom.select(pdb.full, chain="A"))
modes      <- nma(pdb.open)

# Make a PDB trajectory
mktrj(modes, mode=7, file = "6MSM_mktrj.pdb")

# Vector field representation (see Figure 3.)
pymol(modes, mode=7, file = "6MSM_vectField.pdb")

# Calculate the cross-correlation matrix
cm <- dccm(modes)

# Plot a correlation map with plot.dccm(cm)
plot(cm, sse=pdb.open, contour=F, col.regions=bwr.colors(20), at=seq(-1,1,0.1) )

# View the correlations in the structure (see Figure 5.)
pymol(cm, pdb.open, type="launch")

# Deformation energies
defe <- deformation.nma(modes)
defsums <- rowSums(defe$ei[,1:3])

# Fluctuations
flucts <- fluct.nma(modes, mode.inds=seq(7,9))

# Write to PDB files (see Figure 6.) --> J'ai du mal à voir les prot dans pymol, je vais essayer un autre logiciel 
write.pdb(pdb=NULL, xyz=modes$xyz, file="6MSM-defor.pdb", b=defsums)
write.pdb(pdb=NULL, xyz=modes$xyz, file="6MSM-fluct.pdb", b=flucts)

# Overlap analysis

# Closed state of the subunit
pdb.closed <- trim.pdb(pdb.full, atom.select(pdb.full, chain="B"))

# Align closed and open PDBs
aln <- struct.aln(pdb.open, pdb.closed, max.cycles=0, web.args=list(email="anmayroy@gmail.com")) ## Error in as.fasta(aln, id = id) : provide a sequence character matrix/vector
pdb.closed$xyz <- aln$xyz # Pas d'objet aln

# Caclulate a difference vector
xyz <- rbind(pdb.open$xyz[aln$a.inds$xyz], pdb.closed$xyz[aln$a.inds$xyz]) # Pas d'objet aln
diff <- difference.vector(xyz) # pas d'objet xyz

# Calculate overlap, ça marche pas :(
oa <- overlap(modes, diff) # pas d'objet diff
plot(oa$overlap, type='h', xlab="Mode index", ylab="Squared overlap", ylim=c(0,1))
points(oa$overlap, col=1)
lines(oa$overlap.cum, type='b', col=2, cex=0.5)
text(c(1,5)+.5, oa$overlap[c(1,5)], c("Mode 1", "Mode 5"), adj=0)

# Example 2: All-atom normal mode analysis (ENM)

# # Example 2A
# # keep only protein + methotrexate
# pdb <- trim(pdb, "notwater")
# 
# # calculate all-atom NMA with ENM, output calpha only
# m.aa <- aanma(pdb, outmodes="calpha")
# # write summary
# m.aa
# # compare with c-alpha modes ##plot(rac, xlab="PFANM", ylab="C-alpha FF") pour la visualisation 
# m.ca <- nma(pdb)
# rmsip(m.aa, m.ca)
# 
# # Example 2B: Reducing the computational load
# 
# # use rotation and translation of blocks
# m.rtb <- aanma(pdb, rtb=TRUE)
# rmsip(m.aa, m.rtb)
# 
# # use reduced-atom ENM
# m.red <- aanma(pdb, reduced=TRUE)
# rmsip(m.aa, m.red)
# 
# # use both 4-bead and RTB approach
# m.rr  <- aanma(pdb, reduced=TRUE, rtb=TRUE)

# Example 3: Ensemble normal mode analysis

## obtention des PDB 
aa <- get.seq("6MSM_A")
blast <- blast.pdb(aa)
hits <- plot.blast(blast, cutoff = 300) #DEfault cutoff = 368, hits = 12 
#head(hits$pdb.id)
files <- get.pdb(hits$pdb.id, path = "pdbs", split = TRUE, gzip = TRUE)
pdbs <- pdbaln(files, fit=TRUE, web.args=list(email="anmayroy@gmail.com"))
ids <- basename.pdb(pdbs$id)

# Sequence identity
summary( c(seqidentity(pdbs)) )

# NMA on all structures
modes <- nma(pdbs)

# Structure-based clustering
rd <- rmsd(pdbs)
hc <- hclust(as.dist(rd))
grps <- cutree(hc, k=4)

plot(modes, pdbs, col=grps)


