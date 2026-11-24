# Script to run sNMF
#Leilton W Luna
#04-24-2024

#Installation: You must install the devtools package to install LEA

#if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("LEA")

#Create a work directory

library(LEA)

#Convert vcf to geno
vcf2geno("jacana__N27_filtered.vcf.gz")

#Create a project
project = NULL

project = snmf("jacana__N27_filtered.vcf.gz", K = 1:5, tolerance = 0.001, CPU = 6,
               percentage = 0.05, iterations = 1000,  entropy = TRUE, alpha = 100,
               repetitions = 200, ploidy = 2, seed = -1, project = "new")
show(project)

summary(project)

#Plot your cross-entropy result to see which is the best ancestral
#population coefficient
# Plot the data and save the plot in TIFF format
plot(project, col = "blue", pch = 19, cex = 1.2)

ce = cross.entropy(project, K = 2)

ce

#Choose the best run to visualize the population structure
best = which.min(cross.entropy(project, K = 3))

res = Q(project, K =3, run = best)

res

my.colors <- c("gold","red","black")
barplot(t(Q(project, K = 3, run = best)), col = my.colors)

#write Q values
write.csv(res, file = "sNMF_aestiva_ochrocephala_K2.csv")
