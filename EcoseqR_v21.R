# Get arguments
DIR=commandArgs(trailingOnly=TRUE)[1]
SAMP=commandArgs(trailingOnly=TRUE)[2]

# Set environment
setwd(DIR)
options(scipen=100) 
library(data.table)
library(tidyverse)
Sys.setenv('R_MAX_VSIZE'=32000000000)

# Load DCS depth data
DIR=paste(DIR,"/",sep="")
de <- fread(text = gsub(":", "\t", readLines(paste(DIR,SAMP,".f5an0.uniqb.dcs.fil.mpileup.rdc.bed",sep=""))), data.table = F, fill = TRUE)
colnames(de)=c("Chr", "Start", "End", "A.de", "T.de", "C.de", "G.de", "in.de", "del.de")
de$depth.de = de$A.de + de$T.de + de$C.de + de$G.de + de$in.de + de$del.de

# Load original sequence depth data
mpe <- fread(text = gsub(":", "\t", readLines(paste(DIR,SAMP,".f5an0.uniqb.dcs.fil.mpileup.pe.mpileup.txt",sep=""))), data.table = F, fill = TRUE)
colnames(mpe)=c("Chr", "End", "pe", "A.pe", "T.pe", "C.pe", "G.pe", "in.pe", "del.pe")
mpe$depth.pe = mpe$A.pe + mpe$T.pe + mpe$C.pe + mpe$G.pe + mpe$in.pe + mpe$del.pe

# Load common SNP data
snp <- fread(paste(DIR,"/source/dbSnp153Common3.bed",sep=""), data.table = F, fill = TRUE)
colnames(snp)=c("Chr", "Start", "dbSnp153Common")

# Merge data
de.mpe = merge(de, mpe, by.x=c("Chr", "End"), by.y=c("Chr", "End"), all.x=T)
all = merge(de.mpe, snp, by.x=c("Chr", "Start"), by.y=c("Chr", "Start"), all.x=T)

# Filter by depth
all.1.20 <- all[all$depth.pe>=20, ]

# Write output for different depths
all.2.20 <- all.1.20[all.1.20$depth.de>=2, ]
all.3.20 <- all.1.20[all.1.20$depth.de>=3, ] 
all.5.20 <- all.1.20[all.1.20$depth.de>=5, ] 
all.10.20 <- all.1.20[all.1.20$depth.de>=10, ] 
all.20.20 <- all.1.20[all.1.20$depth.de>=20, ]
fwrite(all.1.20,paste(SAMP,".uniqb.dcs.mpileup.1.20.bed",sep=""), quote=F, sep="\t",col.names=F, row.names=F, append=F)
fwrite(all.5.20,paste(SAMP,".uniqb.dcs.mpileup.5.20.bed",sep=""), quote=F, sep="\t",col.names=F, row.names=F, append=F)

# Count common SNP
all.1.20.n <- replace_na(all.1.20, replace = list(dbSnp153Common = 0))
all.1.20.d <- subset(all.1.20.n, dbSnp153Common==0)
all.2.20.d <- all.1.20.d[all.1.20.d$depth.de>=2, ]
all.3.20.d <- all.1.20.d[all.1.20.d$depth.de>=3, ] 
all.5.20.d <- all.1.20.d[all.1.20.d$depth.de>=5, ] 
all.10.20.d <- all.1.20.d[all.1.20.d$depth.de>=10, ] 
all.20.20.d <- all.1.20.d[all.1.20.d$depth.de>=20, ]

# Create summary table for sequencing depth
re=data.frame(DCS_depth=c("number_commonSNP","length_region","number_nucleotide_DCS","length_region_exc_commonSNP","number_nucleotide_DCS_exc_commonSNP"),D1=c(sum(table(all.1.20$dbSnp153Common)),nrow(all.1.20),sum(all.1.20$depth.de, na.rm=TRUE),nrow(all.1.20.d),sum(all.1.20.d$depth.de, na.rm=TRUE)),D2=c(sum(table(all.2.20$dbSnp153Common)),nrow(all.2.20),sum(all.2.20$depth.de, na.rm=TRUE),nrow(all.2.20.d),sum(all.2.20.d$depth.de, na.rm=TRUE)),D3=c(sum(table(all.3.20$dbSnp153Common)),nrow(all.3.20),sum(all.3.20$depth.de, na.rm=TRUE),nrow(all.3.20.d),sum(all.3.20.d$depth.de, na.rm=TRUE)),D5=c(sum(table(all.5.20$dbSnp153Common)),nrow(all.5.20),sum(all.5.20$depth.de, na.rm=TRUE),nrow(all.5.20.d),sum(all.5.20.d$depth.de, na.rm=TRUE)),D10=c(sum(table(all.10.20$dbSnp153Common)),nrow(all.10.20),sum(all.10.20$depth.de, na.rm=TRUE),nrow(all.10.20.d),sum(all.10.20.d$depth.de, na.rm=TRUE)),D20=c(sum(table(all.20.20$dbSnp153Common)),nrow(all.20.20),sum(all.20.20$depth.de, na.rm=TRUE),nrow(all.20.20.d),sum(all.20.20.d$depth.de, na.rm=TRUE)))
fwrite(re,paste(SAMP,".uniqb.depth.results.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)

# Load list of mutations
S.uniqb = read.table(paste(SAMP,".f5an0.uniq1.dcs.fil.mutcount.bed",sep=""), header=F, fill = TRUE)
colnames(S.uniqb) = c("Chr","Start.m","End", "Ref.m", "Alt.m", "depth.m", "vN.m", "bases.m", "A_C_G_T", "misRate","strandRatio")

# Merge list of mutations and regions analyzed
S.uniqb.1.20.mu = merge(S.uniqb, all.1.20, by.x=c("Chr", "End"), by.y=c("Chr", "End"))
fwrite(S.uniqb.1.20.mu,paste(SAMP,".uniqb.1.20.mu.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.1.20.mu1 = S.uniqb.1.20.mu[S.uniqb.1.20.mu$vN.m==1,]
fwrite(S.uniqb.1.20.mu1,paste(SAMP,".uniqb.1.20.mu1.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.2.20.mu = merge(S.uniqb, all.2.20, by.x=c("Chr", "End"), by.y=c("Chr", "End"))
fwrite(S.uniqb.2.20.mu,paste(SAMP,".uniqb.2.20.mu.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.2.20.mu1 = S.uniqb.2.20.mu[S.uniqb.2.20.mu$vN.m==1,]
fwrite(S.uniqb.2.20.mu1,paste(SAMP,".uniqb.2.20.mu1.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.3.20.mu = merge(S.uniqb, all.3.20, by.x=c("Chr", "End"), by.y=c("Chr", "End"))
fwrite(S.uniqb.3.20.mu,paste(SAMP,".uniqb.3.20.mu.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.3.20.mu1 = S.uniqb.3.20.mu[S.uniqb.3.20.mu$vN.m==1,]
fwrite(S.uniqb.3.20.mu1,paste(SAMP,".uniqb.3.20.mu1.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.5.20.mu = merge(S.uniqb, all.5.20, by.x=c("Chr", "End"), by.y=c("Chr", "End"))
fwrite(S.uniqb.5.20.mu,paste(SAMP,".uniqb.5.20.mu.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.5.20.mu1 = S.uniqb.5.20.mu[S.uniqb.5.20.mu$vN.m==1,]
fwrite(S.uniqb.5.20.mu1,paste(SAMP,".uniqb.5.20.mu1.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.10.20.mu = merge(S.uniqb, all.10.20, by.x=c("Chr", "End"), by.y=c("Chr", "End"))
fwrite(S.uniqb.10.20.mu,paste(SAMP,".uniqb.10.20.mu.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.10.20.mu1 = S.uniqb.10.20.mu[S.uniqb.10.20.mu$vN.m==1,]
fwrite(S.uniqb.10.20.mu1,paste(SAMP,".uniqb.10.20.mu1.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.20.20.mu = merge(S.uniqb, all.20.20, by.x=c("Chr", "End"), by.y=c("Chr", "End"))
fwrite(S.uniqb.20.20.mu,paste(SAMP,".uniqb.20.20.mu.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
S.uniqb.20.20.mu1 = S.uniqb.20.20.mu[S.uniqb.20.20.mu$vN.m==1,]
fwrite(S.uniqb.20.20.mu1,paste(SAMP,".uniqb.20.20.mu1.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)

# Create summary table for mutations
S.uniqb.mu=data.frame(DCS_depth_m=c("all_variant_DCS","all_variant_DCS_exc_commonSNP","mutation_DCS","mutation_DCS_exc_commonSNP"),D1=c(sum(table(S.uniqb.1.20.mu$Ref.m)),sum(table(S.uniqb.1.20.mu$dbSnp153Common)),sum(table(S.uniqb.1.20.mu1$Ref.m)),sum(table(S.uniqb.1.20.mu1$dbSnp153Common))), D2=c(sum(table(S.uniqb.2.20.mu$Ref.m)),sum(table(S.uniqb.2.20.mu$dbSnp153Common)),sum(table(S.uniqb.2.20.mu1$Ref.m)),sum(table(S.uniqb.2.20.mu1$dbSnp153Common))), D3=c(sum(table(S.uniqb.3.20.mu$Ref.m)),sum(table(S.uniqb.3.20.mu$dbSnp153Common)),sum(table(S.uniqb.3.20.mu1$Ref.m)),sum(table(S.uniqb.3.20.mu1$dbSnp153Common))), D5=c(sum(table(S.uniqb.5.20.mu$Ref.m)),sum(table(S.uniqb.5.20.mu$dbSnp153Common)),sum(table(S.uniqb.5.20.mu1$Ref.m)),sum(table(S.uniqb.5.20.mu1$dbSnp153Common))), D10=c(sum(table(S.uniqb.10.20.mu$Ref.m)),sum(table(S.uniqb.10.20.mu$dbSnp153Common)),sum(table(S.uniqb.10.20.mu1$Ref.m)),sum(table(S.uniqb.10.20.mu1$dbSnp153Common))), D20=c(sum(table(S.uniqb.20.20.mu$Ref.m)),sum(table(S.uniqb.20.20.mu$dbSnp153Common)),sum(table(S.uniqb.20.20.mu1$Ref.m)),sum(table(S.uniqb.20.20.mu1$dbSnp153Common))))
fwrite(S.uniqb.mu,paste(SAMP,".uniqb.mu.results.txt",sep=""), quote=F, sep="\t",col.names=T, row.names=F, append=F)
