#####################################
###  Version 2023-01-06
#three required arguments to the function
# "Exons.txt" Avoir la liste des exons avec les positions.
# "SNPs.txt" Avoir la liste des SNPs avec leur position sur le B37/38.
# figure  Le nom de la figure


#function 
pwp_syn_nonsyn = function(snps_file = "input/SNPs.txt", exons_file = "input/Exons.txt" , figure = paste0("output/SERPINA1-PWP_",Sys.Date(),".png")) {


snps <- read.table(snps_file, header = T, sep = "\t", stringsAsFactors = F)
exons <- read.table(exons_file, header = T, sep = "\t")

snps.order <- snps[order(snps$SNP.Position, decreasing = T) , ]


png(figure, height = 3840, width = 5120, pointsize = 75)

par(mar = c(6,0,0,0))
plot(x=0,y=0, xlim = c(94857029, 94843083), ylim = c(0,15), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", cex.lab = 1.5)

box(col = "white", lwd = 10)

## Ligne de base
rect(94843083, 1, 94857029, 2, col = "black", border = NA)


## Axis 1
axis(1, at = c(94843000, 94844000, 94845000, 94846000, 94847000, 94848000, 94849000, 94850000, 94851000, 94852000, 94853000, 94854000, 94855000, 94856000, 94857000), labels = c("94,843", "94,844", "94,845", "94,846", "94,847", "94,848", "94,849", "94,850", "94,851", "94,852", "94,853", "94,854", "94,855", "94,856", "94,857"), cex.axis = 1.1, lwd = 15, line = 1)

mtext("Chromosome 14 (kb)", 1, cex = 1.5, line = 4)


## Exons #darkgrey = non codant
for(i in seq_along(exons$Exon.Start))
{
   if(i ==1) rect(exons$Exon.End[i], 3, exons$Exon.Start[i], 0, col = "darkgrey", border = NA)
   if(i!= 1) rect(exons$Exon.End[i], 3, exons$Exon.Start[i], 0, col = "black", border = NA)

   if((i!=2) & (i!=3)) spacer = 0 
   if(i==2)  spacer = 50
   if(i==3)  spacer = -50 
  
  text((((exons$Exon.End[i]-exons$Exon.Start[i])/2)+exons$Exon.Start[i])+spacer, -0.3, exons$Exon.Name[i], cex = 1.5)
}

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

message(paste0('Done preparing figure ',figure,' --- Time is: ',Sys.time()))

}

#Run the function below...
pwp_syn_nonsyn()