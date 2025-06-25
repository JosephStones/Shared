#### Getting fasta sequences for ORFs, and control proteins
#### and calculate the number of residues able to be targetted by enzymes
#Joseph Stones: 25-06-25


# load libraries
library(ggplot2)
library(Biostrings)

# Fasta sequences to check
orfFasta<-readAAStringSet("fastafile") #change "fastafile" to input fasta file name


# Enzyme list, examples from:
# https://www.promega.co.uk/products/mass-spectrometry/proteases-and-surfactants/#:~:text=Available%20proteases%20for%20optimal%20digestion,%2C%20rAsp%2DN%20and%20more.
# format = "PROTEASE-CLEAVAGE_TERMINAL-TARGET_RESIDUES"
enzymesMS<-c("trypsin-c-R|K",
             "lys-c-K",
             "arg-c-R",
             "chymotrypsin-c-Y|F|W",
             "asp-n-D",
             "proalanase-c-P|A",
             "glu-c-E",
             "pepsin-c-F|L|W|Y")



## get number of specific residues
getResidues<-function(aaseq,enzyme,output){
  
  #get enzyme info, pos1 = name, 2=C/N terminus, 3=AA code
  enzyme<-unlist(strsplit(enzyme,split="-")) 
  #positions of target residues
  residuePos<-gregexpr(enzyme[3],aaseq)
  
  if (residuePos[[1]][1] %in% (-1)) {return(list(0,0))} # if no cleavage sites return 0
  
  if (enzyme[1]=="trypsin") {
    # check for P, as cleavage won't occur on amino acid after this
    Psites<-gregexpr("P",aaseq)
    
    # if no Psites, skip this, otherwise remove cleavage positions from the list
    if (Psites[[1]][1] %in% (-1)) {
      
    } else {
      Psites[[1]]<-Psites[[1]]-1 # see if position before P is any cleavage residues
      residuePos[[1]]<-residuePos[[1]][which(!(residuePos[[1]] %in% Psites[[1]]))] #which of residue positions are not before a P
    }
    
    if (length(residuePos[[1]]) == 0) {return(list(0,0))} # check again for no cleavage sites
  }
  
  #if a target at the end of aa seq for C cleavage, ignore
  # or at start for N cleavage
  if (enzyme[2]=="c"){
    residuePos<-residuePos[[1]][grep(nchar(aaseq),residuePos[[1]],invert = T)] #remove aa if at end
    residuePos<-sort(residuePos,decreasing = F)
    
    peplengths<-lapply(1:(length(residuePos)+1),function(x){
      if (x==1){startpos<-0} else {startpos<-residuePos[x-1]} # set start position of fragment
      if (x==(length(residuePos)+1)){csite<-nchar(aaseq)} else {csite<-residuePos[x]}# end position of fragment
      peplength<-csite-startpos
      return(peplength)
      
    })
    
  } else if (enzyme[2]=="n") {
    residuePos<-residuePos[[1]][grep("1",residuePos[[1]],invert = T)] # remove aa if in position 1
    residuePos<-sort(residuePos,decreasing = F)
    
    peplengths<-lapply(1:(length(residuePos)+1),function(x){
      if (x==1){startpos<-1} else {startpos<-residuePos[x-1]} # set start position of fragment
      if (x==(length(residuePos)+1)){csite<-(nchar(aaseq)+1)} else {csite<-residuePos[x]}# end position of fragment
      peplength<-(csite)-startpos
      return(peplength)
    })
  }
  
  peplengths<-unlist(peplengths)
  
  return(list(peplengths,length(residuePos)))
}



### length of residues that may be cleaved by trypsin (R=arginine, K=lysine)
cleaveInfo<-lapply(enzymesMS,function(enzyme){
  print(enzyme)
  cleaveInfo<-lapply(orfFasta,function(x){
    getResidues(as.character(x),enzyme)
  })
  
  # get frequency of peptides by length
  peptideLengths<-lapply(cleaveInfo,`[[`,1)
  peptideLengths<-as.data.frame(table(unlist(peptideLengths)))
  colnames(peptideLengths)<-c("peptideLength","Frequency")
  peptideLengths$enzyme<-enzyme
  
  # number of cleavage sites
  cleaveSites<-lapply(cleaveInfo,`[[`,2)
  cleaveSites<-data.frame(table(unlist(cleaveSites)))
  colnames(cleaveSites)<-c("nCleaveSites","Frequency")
  cleaveSites$enzyme<-enzyme
  
  return(list(peptideLengths,cleaveSites))
})

# Combine peptide length info for all enzymes
peptideLengths<-lapply(cleaveInfo,`[[`,1)
peptideLengths<-do.call(rbind,peptideLengths)
peptideLengths$peptideLength<-as.integer(peptideLengths$peptideLength)-1 # from factor to integer

# Combine cleavage site info for all enzyme
cleaveSites<-lapply(cleaveInfo,`[[`,2)
cleaveSites<-do.call(rbind,cleaveSites)
cleaveSites$nCleaveSites<-as.integer(cleaveSites$nCleaveSites)-1 # from factor to integer


## Modify scale_x_continuous boundaries depending on data distribution
# plot pep lengths
ggplot(peptideLengths,aes(x=peptideLength,y=Frequency,group=enzyme))+geom_col(fill="steelblue4",colour="white")+
  geom_vline(xintercept=c(6.5,40.5),linetype="dashed")+
  coord_cartesian(expand=0)+
  scale_x_continuous(limits=c(-1,42),breaks=seq(0,42,3))+
  facet_wrap(as.formula(". ~enzyme"))+
  xlab("Peptide length (0=no cleavage)")+
  ylab("Number of peptides at specified length")+
  ggtitle("Length of peptides after cleavage for typical MS proteases")+
  theme(plot.background = element_blank(),axis.text=element_text(size=10),
        axis.title = element_text(size=12),
        plot.title = element_text(size=14))

# plot number of cleave sites
ggplot(cleaveSites,aes(x=nCleaveSites,y=Frequency,group=enzyme))+geom_col(fill="steelblue4",colour="white")+
  coord_cartesian(expand=0)+
  scale_x_continuous(limits=c(-1,30),breaks=seq(0,30,3))+
  facet_wrap(as.formula(". ~enzyme"))+
  xlab("Number of cleavage sites (0=no cleavage)")+
  ylab("Number of fasta sequences with specified cleavage sites")+
  ggtitle("Number of potential cleavage sites for typical MS proteases")+
  theme(plot.background = element_blank(),axis.text=element_text(size=10),
        axis.title = element_text(size=12),
        plot.title = element_text(size=14))

