###################
#
#   Count # Transcription Factor Binding Sites (TFBS) for a set of species and a set of transcription factors 
#   for sequences analyzed by JASPAR
#     KC June 2026
###################
#
# Inputs:
#   Need JASPAR file for each species and sequence region as analyzed by JASPAR.
#      Go to JASPAR and add TF matrices to cart. Then go to View Cart. Under SCAN, enter sequence and enter.
#      Display all, copy and paste results into text file. Do not need header of column labels.
#
# Things to choose here
#   choose which gene (SWS2, RH2, LWS) and region (promoter, lcr) to analyze. These should be part of JASPAR filenames
#      Usually names are Species nickname, gene and region, e.g. TilSWS2apro for Tliapia SWS2a gene promoter
#   choose which analysis:
#     Sens : analyzes all species for one TF vs a range of sensitivities (e.g 0.8, 0.82,...1.0)  
#     allTF : counts all TFs and all species for one sensitivity
#     TFdist for one sens and all TFs for #cichlid TFBS - #other fish TFBS
#   List the TFs that you want to analyze / count in the vector, TFnames
#   Add list of species nicknames that match beginning of species JASPAR filenames, e.g Til for Tilapia or Dre for zebrafish

## Function to fixup file with JASPAR output: replace Matrix name with TF name and make all caps
File_TFnames<-function(datafile) {
  ## Note if you paste data directly from JASPAR into text file you may get extra space after each
  ## TF name which will lead to zero counts below. Copy from JASPAR to Excel and then to text file
  ## avoids that problem
  size<-nrow(datafile); hold1<-matrix(data=NA,nrow=size, ncol=3)
  hold1<-sapply(datafile[2], strsplit,"\\.")
  hold2<-matrix(unlist(sapply(datafile[2], strsplit,"\\.")), nrow=size, byrow= TRUE)
  hold3<-toupper(hold2)
  datafile[,2]<- hold3[,3]
  return(datafile) }
########

###################
# Choose gene region to analyze"
##################
Anal2do<-"SWS2 promoter"  # these set the file name to be analyzed 
#     "SWS2 promoter"   SWS2pro
# other choices (no data in this Github version) :
#     "LWS promoter"    LWSpro
#     "LWS lcr"         LWSlcr
#     "RH2 promoter"    RH2pro
#     "RH2 lcr"         RH2lcr
#################
# Choose Plot ("allTF", "sens", or "TFdist")
#################
Plottype<-"allTF" 
# choices: 
#     "sens" pick one TF and adjust Y axis (Ymax)  
if (Plottype=="sens") {
  TFonly<-"RAX2"
  Ymax<-40
}

# "allTF" pick one sensitivity, plot # TFBS for all TFs for all species (TFs are numbered across bottom)
#     chose one sensitivity (sens2plot) and adjust Y axis (Ymax)
#     dotshift shifts the x axis value for each species so dots do not overlap
if (Plottype=="allTF") {
  sens2plot<-0.9
  dotshift<-0.05
  Xmin<-0; Ymax<-90
}
# "TFdist" for one TF plots # of cichlid TFBS minus # of other fish TFBS
#   choose one TF (TFonly) and adjust Y axis (Ymax)
#   set a displacement for each species so no dots overlap each other
if (Plottype=="TFdist") {
  TFonly<-"RAX2"
  Ymax<-20
}
# determine whether results are saved to a file or not
save2file<-"Y"

# List TFs to be analyzed. Note all TFs are capital. The function File_TFnames converts all TF names to the capitalized vertion
TFnames<-c("RAX2", "TBX2", "TBX4", "CREM", "THRB", "SOX5", "SOX6", "FOXB1", "NFIX", "RORC", "RARG", "RXRG","SIX3","RARA::RXRG")
#plot each TF in a different color
TFcolors<-c("blue","red","orange","grey80", "gold2", "pink", "hotpink", "lightblue","grey20","maroon1","violet","darkorchid1", "cyan", "purple")
TFmax<-length(TFnames)

# set sensitivity range
min_thresh<-0.8   # lowest score
max_thresh<-1     # highest score
step_thresh<-0.02 # step size from low to high
sens_number<-(sens2plot-min_thresh)/step_thresh+1

# JASPAR output is named SpeciesGeneRegion. Here are the Species names
species1<-"Til"     # Tilapia, Oreochromis niloticus. so e.g. JASPAR output file called TilSWS2pro.txt
species2<-"Abur"    # Astatotilapia burtoni
species3<-"Pnye"    # Pundamilia nyererei
species4<-"Cgun"    # Chromotilapia guntheri
species5<-"Acit"    # Amphilophus citrinellus
species6<-"Dre"     # Danio rerio
species7<-"Gac"     # Gasterosteus aculeatus, three spine stickleback
species8<-"Med"     # Medaka, Oryzias latipes
species9<-"Pam"     # Pseudopleuronectes americanus, winter flounder
species10<-"Pre"    # Poecilius reticulatus, guppy

# assign names to vectors and colors to species; here cichlids are blue to green and non cichlids are red to purple
species_names<-c(species1,species2,species3,species4,species5,species6,species7,species8,species9,species10)
species_colors<-c("blue","cyan","lightblue","darkgreen","green","purple","pink","orange","darkorange","red","black")
#
#Set the threshold of Jaspar score to include that TFBS 
#
# Need separate section for each regulatory region to be analyzed. Here we just include the SWS2 opsin promoter
# and set up for the data folders and save folders

if (Anal2do=="SWS2 promoter") {
directory<-"~/yourdirectory/Opsindata/"
savedirectory<-"~/yourdirectory/Opsinoutput/"
suffix1<-"S2"
suffix2<-"pro.txt"
# some species have S2 gene and some have S2a gene - so add the extra "a" here
file1<-paste(directory,species1,suffix1,"a",suffix2,sep="")
file2<-paste(directory,species2,suffix1,"a",suffix2,sep="")
file3<-paste(directory,species3,suffix1,"a",suffix2,sep="")
file4<-paste(directory,species4,suffix1,"a",suffix2,sep="")
file5<-paste(directory,species5,suffix1,"a",suffix2,sep="")
file6<-paste(directory,species6,suffix1,suffix2,sep="")
file7<-paste(directory,species7,suffix1,suffix2,sep="")
file8<-paste(directory,species8,suffix1,"a",suffix2,sep="")
file9<-paste(directory,species9,suffix1,suffix2,sep="")
file10<-paste(directory,species10,suffix1,"a",suffix2,sep="")
}

# Read in Jaspar data
data1<-read.table(file1, header=FALSE,sep="\t")
data2<-read.table(file2, header=FALSE,sep="\t")
data3<-read.table(file3, header=FALSE,sep="\t")
data4<-read.table(file4, header=FALSE,sep="\t")
data5<-read.table(file5, header=FALSE,sep="\t")
data6<-read.table(file6, header=FALSE,sep="\t")
data7<-read.table(file7, header=FALSE,sep="\t")
data8<-read.table(file8, header=FALSE,sep="\t")
data9<-read.table(file9, header=FALSE,sep="\t")
data10<-read.table(file10, header=FALSE,sep="\t")
#
# Add JASPAR column names
datacolnames<-c("MatrixID","Name","Score","Relative_score","SequenceID","Start","End","Strand","PredictedSequence")
colnames(data1)<-datacolnames
colnames(data2)<-datacolnames
colnames(data3)<-datacolnames
colnames(data4)<-datacolnames
colnames(data5)<-datacolnames
colnames(data6)<-datacolnames
colnames(data7)<-datacolnames
colnames(data8)<-datacolnames
colnames(data9)<-datacolnames
colnames(data10)<-datacolnames

# Use function File_TFnames to truncate matrix name to just include TF name and change to all uppercase 

TFdata1<-File_TFnames(data1)
TFdata2<-File_TFnames(data2)
TFdata3<-File_TFnames(data3)
TFdata4<-File_TFnames(data4)
TFdata5<-File_TFnames(data5)
TFdata6<-File_TFnames(data6)
TFdata7<-File_TFnames(data7)
TFdata8<-File_TFnames(data8)
TFdata9<-File_TFnames(data9)
TFdata10<-File_TFnames(data10)

# count all sites that have relative score above a threshold and then vary this threshold
threshold_set<-seq(from = min_thresh, to = max_thresh, by=step_thresh)
#colnames(threshold_set)<-"Threshold"
threshold_steps<-length(threshold_set)
#Set up Matrices for counts
TFcounts1<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFcounts2<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFcounts3<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFcounts4<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFcounts5<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFcounts6<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFcounts7<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFcounts8<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFcounts9<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFcounts10<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)

# Set up summary files
TFdist<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFWelch<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFGWelch<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)
TFLWelch<-as.data.frame(NA, nrow=threshold_steps,ncol=TFmax)

#Count number of TFBS that are above varying theshold for each of 10 species
for (i in (1:threshold_steps)) {
    score<-threshold_set[i]
    for (j in 1:TFmax) {
    TFcounts1[i,j]<-sum(TFdata1$Name==TFnames[j] & TFdata1$Relative_score>=score )
  }
}
colnames(TFcounts1)<-TFnames

for (i in (1:threshold_steps)) {
  score<-threshold_set[i]
  for (j in 1:TFmax) {
    TFcounts2[i,j]<-sum(TFdata2$Name==TFnames[j] & TFdata2$Relative_score>=score )
  }
}
colnames(TFcounts2)<-TFnames

for (i in (1:threshold_steps)) {
  score<-threshold_set[i]
  for (j in 1:TFmax) {
    TFcounts3[i,j]<-sum(TFdata3$Name==TFnames[j] & TFdata3$Relative_score>=score )
  }
}
colnames(TFcounts3)<-TFnames

for (i in (1:threshold_steps)) {
  score<-threshold_set[i]
  for (j in 1:TFmax) {
    TFcounts4[i,j]<-sum(TFdata4$Name==TFnames[j] & TFdata4$Relative_score>=score )
  }
}
colnames(TFcounts4)<-TFnames

for (i in (1:threshold_steps)) {
  score<-threshold_set[i]
  for (j in 1:TFmax) {
    TFcounts5[i,j]<-sum(TFdata5$Name==TFnames[j] & TFdata5$Relative_score>=score )
  }
}
colnames(TFcounts5)<-TFnames

for (i in (1:threshold_steps)) {
  score<-threshold_set[i]
  for (j in 1:TFmax) {
    TFcounts6[i,j]<-sum(TFdata6$Name==TFnames[j] & TFdata6$Relative_score>=score )
  }
}
colnames(TFcounts6)<-TFnames

for (i in (1:threshold_steps)) {
  score<-threshold_set[i]
  for (j in 1:TFmax) {
    TFcounts7[i,j]<-sum(TFdata7$Name==TFnames[j] & TFdata7$Relative_score>=score )
  }
}
colnames(TFcounts7)<-TFnames

for (i in (1:threshold_steps)) {
  score<-threshold_set[i]
  for (j in 1:TFmax) {
    TFcounts8[i,j]<-sum(TFdata8$Name==TFnames[j] & TFdata8$Relative_score>=score )
  }
}
colnames(TFcounts8)<-TFnames

for (i in (1:threshold_steps)) {
  score<-threshold_set[i]
  for (j in 1:TFmax) {
    TFcounts9[i,j]<-sum(TFdata9$Name==TFnames[j] & TFdata9$Relative_score>=score )
  }
}
colnames(TFcounts9)<-TFnames

for (i in (1:threshold_steps)) {
  score<-threshold_set[i]
  for (j in 1:TFmax) {
    TFcounts10[i,j]<-sum(TFdata10$Name==TFnames[j] & TFdata10$Relative_score>=score )
  }
}
colnames(TFcounts10)<-TFnames

## compare 5 cichlids to other 5 fish
for (j in 1:TFmax) {
    for (i in (1:threshold_steps)) {
    TFdist[i,j]<-((TFcounts1[i,j]+TFcounts2[i,j]+TFcounts3[i,j]+TFcounts4[i,j]+TFcounts5[i,j])/5-(TFcounts6[i,j]+TFcounts7[i,j]+TFcounts8[i,j]+TFcounts9[i,j]+TFcounts10[i,j])/5)
  }
}

# Save counts for Welch stats
TF4W1<-TFcounts1
TF4W2<-TFcounts2
TF4W3<-TFcounts3
TF4W4<-TFcounts4
TF4W5<-TFcounts5
TF4W6<-TFcounts6
TF4W7<-TFcounts7
TF4W8<-TFcounts8
TF4W9<-TFcounts9
TF4W10<-TFcounts10

# Add threshold scores for data output and plotting
TFcounts1<-cbind.data.frame(threshold_set,TFcounts1)
TFcounts2<-cbind.data.frame(threshold_set,TFcounts2)
TFcounts3<-cbind.data.frame(threshold_set,TFcounts3)
TFcounts4<-cbind.data.frame(threshold_set,TFcounts4)
TFcounts5<-cbind.data.frame(threshold_set,TFcounts5)
TFcounts6<-cbind.data.frame(threshold_set,TFcounts6)
TFcounts7<-cbind.data.frame(threshold_set,TFcounts7)
TFcounts8<-cbind.data.frame(threshold_set,TFcounts8)
TFcounts9<-cbind.data.frame(threshold_set,TFcounts9)
TFcounts10<-cbind.data.frame(threshold_set,TFcounts10)
TFdist<-cbind.data.frame(threshold_set,TFdist)

par(mar=c(4,4,3,6), xpd=FALSE) 
clip=c(4,4,3,6)

if (Plottype=="sens") {   # this would plot counts vs threshold score  

TFnumber<-which(TFnames==TFonly)
ytext<-paste("# of TFBS for ",TFonly,sep="")
plot(TFcounts1[,1],TFcounts1[,(TFnumber+1)], xlim=c(min_thresh,max_thresh),ylim=c(0,Ymax),xlab="Sensitivity threshold", ylab=ytext,col=species_colors[1], pch=19)
lines(TFcounts1[,1],TFcounts1[,(TFnumber+1)], type="b", col=species_colors[1], pch=19)
lines(TFcounts2[,1],TFcounts2[,(TFnumber+1)], type="b", col=species_colors[2], pch=19)
lines(TFcounts3[,1],TFcounts3[,(TFnumber+1)], type="b", col=species_colors[3], pch=19)
lines(TFcounts4[,1],TFcounts4[,(TFnumber+1)], type="b", col=species_colors[4], pch=19)
lines(TFcounts5[,1],TFcounts5[,(TFnumber+1)], type="b", col=species_colors[5], pch=19)
lines(TFcounts6[,1],TFcounts6[,(TFnumber+1)], type="b", col=species_colors[6], pch=19)
lines(TFcounts7[,1],TFcounts7[,(TFnumber+1)], type="b", col=species_colors[7], pch=19)
lines(TFcounts8[,1],TFcounts8[,(TFnumber+1)], type="b", col=species_colors[8], pch=19)
lines(TFcounts9[,1],TFcounts9[,(TFnumber+1)], type="b", col=species_colors[9], pch=19)
lines(TFcounts10[,1],TFcounts10[,(TFnumber+1)], type="b", col=species_colors[10], pch=19)
lines(TFdist[,1],TFcounts5[,(TFnumber+1)], type="b", col=species_colors[11], pch=19)

par(xpd=TRUE)
legend(1.01,Ymax,legend=species_names,col=species_colors, pch=19)
title(main=paste(Anal2do,": # of TFBS above sensitivity for TF:", TFonly))
}

if (Plottype=="allTF") {   #
#Plot # of TFBS for all TFs for 5 species
  xplot<-c(1:14)
  matrix1<-rbind(c(1:14),TFcounts1[sens_number, 2:15])
  ytext<-paste("# of TFBS for sensitivity of",sens2plot,sep="")
#plot(matrix1[1,], matrix1[2,], xlim=c(0,14),ylim=c(0,Ymax),xlab="holding this", ylab=ytext,col=species_colors[1], pch=17)
  plot(xplot, TFcounts1[sens_number,2:15], xlim=c(1,14),ylim=c(0,Ymax),xlab="Transcription factor", ylab=ytext,col=species_colors[1], pch=17)
  points(xplot+dotshift,TFcounts2[sens_number,2:15], col=species_colors[2], pch=19)
  points(xplot+2*dotshift,TFcounts3[sens_number,2:15], col=species_colors[3], pch=19)
  points(xplot+3*dotshift,TFcounts4[sens_number,2:15], col=species_colors[4], pch=19)
  points(xplot+4*dotshift,TFcounts5[sens_number,2:15], col=species_colors[5], pch=19)
  points(xplot+5*dotshift,TFcounts6[sens_number,2:15], col=species_colors[6], pch=19)
  points(xplot+6*dotshift,TFcounts7[sens_number,2:15], col=species_colors[7], pch=19)
  points(xplot+7*dotshift,TFcounts8[sens_number,2:15], col=species_colors[8], pch=19)
  points(xplot+8*dotshift,TFcounts9[sens_number,2:15], col=species_colors[9], pch=19)
  points(xplot+9*dotshift,TFcounts10[sens_number,2:15], col=species_colors[10], pch=19)
par(xpd=TRUE)
  legend(15,Ymax,legend=species_names,col=species_colors, pch=19)
  title(main=paste(Anal2do,": # of TFBS for TFs at sensitivity:",sens2plot))
  plot_name<-paste(savedirectory,"SWS2_10species_TFcntsallTF.png",sep="")
}

if (Plottype=="TFdist") {
  baseline<-matrix(c(0.8,1,0,0),nrow=2,ncol=2)
  TFBSthresh<-matrix(c(0.8,1,2,2),nrow=2,ncol=2)
  TFnumber<-which(TFnames==TFonly)
  ytext<-paste("# ",TFonly," TFBS for cichlids - fishes ",sep="")
  plot(TFdist[,1],TFdist[,(TFnumber+1)], xlim=c(min_thresh,max_thresh),ylim=c(-Ymax,Ymax),xlab="Sensitivity threshold", ylab=ytext,col=species_colors[6], pch=19)
  lines(TFdist[,1],TFdist[,(TFnumber+1)], type="b", col="black", pch=19)
  lines(baseline[,1],baseline[,2], col="gray")
  lines(TFBSthresh[,1],TFBSthresh[,2], col="gray",lty=2)
  
  #  par(xpd=TRUE)
  title(main=paste(Anal2do,": #TFBS cichlid-#others vs sensitivity for TF:", TFonly))
}

if (save2file=="Y") {
   Allresults<-rbind(TFcounts1,TFcounts2,TFcounts3,TFcounts4,TFcounts5,TFcounts6,TFcounts7,TFcounts8,TFcounts9,TFcounts10)
   outfile<-paste(savedirectory,"Allspecies",Anal2do,".txt",sep="")
   write.table(Allresults,file=outfile,sep="\t",quote=FALSE)
}
# Use Welch test to calc if sign diff btn cichlids and other fish - this can be two sided (default), or for whether group 1 is "greater" or "less" than group 2 (alternative string)
for (j in 1:TFmax) {
  for (i in (1:threshold_steps)) {
    group1<-c(TF4W1[i,j],TF4W2[i,j],TF4W3[i,j],TF4W4[i,j],TF4W5[i,j])
    group2<-c(TF4W6[i,j],TF4W7[i,j],TF4W8[i,j],TF4W9[i,j],TF4W10[i,j])
    if (((group1[1]==group1[2]) & (group1[1]==group1[3]) & (group1[1]==group1[4]) & (group1[1]==group1[5])) & ((group2[1]==group2[2]) & (group2[1]==group2[3]) & (group2[1]==group2[4]) & (group2[1]==group2[5]))) {
#    if ((group1[1]==group1[2]) & (group1[1]==group1[3]) & (group1[1]==group1[4]) & (group1[1]==group1[5])) {
      TFWelch[i,j]<-NA
      TFGWelch[i,j]<-NA
      TFLWelch[i,j]<-NA
    } else {
    Wresults<-t.test(group1,group2)
    WGresults<-t.test(group1,group2,alternative="greater")
    WLresults<-t.test(group1,group2,alternative="less")
    TFWelch[i,j]<-Wresults$p.value
    TFGWelch[i,j]<-WGresults$p.value
    TFLWelch[i,j]<-WLresults$p.value
    }
  }
}
colnames(TFWelch)<-TFnames
colnames(TFGWelch)<-TFnames
colnames(TFLWelch)<-TFnames

TFWelch<-cbind.data.frame(threshold_set,TFWelch)
TFGWelch<-cbind.data.frame(threshold_set,TFGWelch)
TFLWelch<-cbind.data.frame(threshold_set,TFLWelch)

if (save2file=="Y") {
  Allresults<-rbind(TFcounts1,TFcounts2,TFcounts3,TFcounts4,TFcounts5,TFcounts6,TFcounts7,TFcounts8,TFcounts9,TFcounts10)
  outfile2<-paste(savedirectory,"Allspecies",Anal2do,"Welch2side.txt",sep="")
  outfile3<-paste(savedirectory,"Allspecies",Anal2do,"WelchGreater.txt",sep="")
  outfile4<-paste(savedirectory,"Allspecies",Anal2do,"WelchLess.txt",sep="")
  write.table(TFWelch,file=outfile2,sep="\t",quote=FALSE)
  write.table(TFGWelch,file=outfile3,sep="\t",quote=FALSE)
  write.table(TFLWelch,file=outfile4,sep="\t",quote=FALSE)
}
