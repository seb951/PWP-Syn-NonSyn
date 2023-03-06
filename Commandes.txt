#####################################
#####################################
###  Version 2015-01-06
###
###
#####################################


#####################################
#####################################
###  Étape 1: Avoir la liste des exons avec les positions.

# "Exons.txt"

#####################################


#####################################
#####################################
###  Étape 2: Avoir la liste des SNPs avec leur position sur le B37/38.

# "SNPs.txt"

#####################################


#####################################
#####################################
###  Étape 3: Graphique

snps <- read.table("SNPs.txt", header = T, sep = "\t", stringsAsFactors = F)
exons <- read.table("Exons.txt", header = T, sep = "\t")

snps.order <- snps[order(snps$SNP.Position, decreasing = T) , ]

###
###
###

png("SERPINA1-PWP_2015-01-07.png", height = 3840, width = 5120, pointsize = 75)

par(mar = c(6,0,0,0))
plot(x=0,y=0, xlim = c(94857029, 94843083), ylim = c(0,15), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", cex.lab = 1.5)

box(col = "white", lwd = 10)

## Ligne de base
rect(94843083, 1, 94857029, 2, col = "black", border = NA)


## Axis 1
axis(1, at = c(94843000, 94844000, 94845000, 94846000, 94847000, 94848000, 94849000, 94850000, 94851000, 94852000, 94853000, 94854000, 94855000, 94856000, 94857000), labels = c("94,843", "94,844", "94,845", "94,846", "94,847", "94,848", "94,849", "94,850", "94,851", "94,852", "94,853", "94,854", "94,855", "94,856", "94,857"), cex.axis = 1.1, lwd = 15, line = 1)

mtext("Chromosome 14 (kb)", 1, cex = 1.5, line = 4)


## Exons #darkgrey = non codant
rect(exons$Exon.End[1], 3, exons$Exon.Start[1], 0, col = "darkgrey", border = NA)
rect(exons$Exon.End[2], 3, exons$Exon.Start[2], 0, col = "black", border = NA)
rect(exons$Exon.End[3], 3, exons$Exon.Start[3], 0, col = "black", border = NA)
rect(exons$Exon.End[4], 3, exons$Exon.Start[4], 0, col = "black", border = NA)
rect(exons$Exon.End[5], 3, exons$Exon.Start[5], 0, col = "black", border = NA)
rect(exons$Exon.End[6], 3, exons$Exon.Start[6], 0, col = "black", border = NA)
rect(exons$Exon.End[7], 3, exons$Exon.Start[7], 0, col = "black", border = NA)
rect(exons$Exon.End[8], 3, exons$Exon.Start[8], 0, col = "black", border = NA)

text((((exons$Exon.End[1]-exons$Exon.Start[1])/2)+exons$Exon.Start[1]), -0.3, "1a", cex = 1.5)
text((((exons$Exon.End[2]-exons$Exon.Start[2])/2)+exons$Exon.Start[2])+50, -0.3, "1b", cex = 1.5)
text((((exons$Exon.End[3]-exons$Exon.Start[3])/2)+exons$Exon.Start[3])-50, -0.3, "1c", cex = 1.5)
text((((exons$Exon.End[4]-exons$Exon.Start[4])/2)+exons$Exon.Start[4]), -0.3, "2", cex = 1.5)
text((((exons$Exon.End[5]-exons$Exon.Start[5])/2)+exons$Exon.Start[5]), -0.3, "3", cex = 1.5)
text((((exons$Exon.End[6]-exons$Exon.Start[6])/2)+exons$Exon.Start[6]), -0.3, "4", cex = 1.5)
text((((exons$Exon.End[8]-exons$Exon.Start[7])/2)+exons$Exon.Start[7]), -0.3, "5", cex = 1.5)


## SNPs-Position sur ligne de base
for (i in 1:nrow(snps.order))
{
	rect(snps.order$SNP.Position[i]-4, 2, snps.order$SNP.Position[i]+4, 4, col = "black", border = NA)
}


dif = ((94857029 - 94843083) / nrow(snps)) #103 + des peanuts. 

## SNP rs
for (i in 1:nrow(snps.order))
{
	segments(snps.order$SNP.Position[i], 4, (94857029 - (i * dif)), 6, lwd = 4)
	rect((94857029 - (i * dif)-6), 6, (94857029 - (i * dif)+6), 7, col = "black", border = NA)
	text(x = (94857029 + 30 - ((i-1) * dif)), y = 7.1, labels = snps.order$SNP[i], col = snps.order$SNP.Couleur[i], srt = 90, cex = 0.5, pos = 4)
}

dev.off()

#####################################
