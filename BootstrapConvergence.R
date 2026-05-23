# BOOTSTRAP CONVERGENCE

#load libraries
library(sleuth)
library(ggplot2)
library(data.table)

#read isoform RNA-seq data for a file
kallh5<-read_kallisto_h5("kallisto/abundance.h5") # bootstrap values - for estimating technical variance
kalltsv<-read_kallisto_tsv("kallisto/abundance.tsv") # deterministic feature counts for isoforms

##generic ggplot theme
gentheme<-theme(panel.background=element_blank(),
                panel.grid = element_blank(), 
                panel.border=element_rect(fill=NA,colour="black"))

## combine all bootstrap values, convert to table
# (each bootstrap estimate is separate in default data)
bsEsts<-lapply(kallh5$bootstrap,function(y){
  y$est_counts
})
bsEsts<-do.call(cbind,bsEsts)
rownames(bsEsts)<-kallh5$abundance$target_id # rename rows as transcript IDs
bsEsts<-as.data.frame(bsEsts)

### highest difference between tsv abundance + mean bootstrap value
### getting bootstraps that show highest standard deviation
bsm<-apply(bsEsts,1,function(x) {mean(x)}) # bootstrap means
kest<-kallh5$abundance$est_counts # kallisto estimated counts
names(kest)<-kallh5$abundance$target_id # name after transcripts IDs


all(names(bsm)==names(kest)) #check datasets match

# get 50 transcripts with highest scale difference between bootstrap estimated and deterministic counts
diff<-kest/bsm 
diff[is.nan(diff)]<-0 # replace NaN with 0
diff<-diff[order(diff,decreasing=T)]
topmd<-names(diff[1:50])
topmd<-bsEsts[topmd,]

#### estimated cumulative mean vs no .bootstrap
#for each transcript id
estmean<-lapply(1:nrow(topmd),function(x) {
  
  # From 1 to the total number of bootstraps, calculate the cumulative mean of estimated counts
  estmean<-unlist(lapply(1:length(kallh5$bootstrap),function(y) {
    bsv<-unlist(as.vector(topmd[x,1:y]))#estimated counts for 'y' number of bootstraps
    estmean<- mean(bsv) #mean of estimated counts
  }))  
  
  return(estmean)
})

#join + reformat
estmeans<-as.data.table(do.call(rbind,estmean))
estmeans$ID<-rownames(topmd)

# get the abundance values for the highest standard deviation isoforms
abundancecounts<-kallh5$abundance[kallh5$abundance$target_id %in% estmeans$ID,]
estmeans<-merge(x=estmeans,by.x="ID",y=abundancecounts[,c("target_id","est_counts")],by.y="target_id")
estmeans<-melt(estmeans,id.vars = c("ID","est_counts"))
levels(estmeans$variable)<-seq(1,1000) # 1000 was the total number of bootstraps I used - alter accordingly


## bootstrap cumulative mean and est_count used in abundance file
ggplot(estmeans,aes(x=variable,y=value,group=ID))+geom_line()+
  geom_line(data=estmeans,aes(x=variable,y=est_counts,group=ID),colour="steelblue4",linetype="dashed")+  
  facet_wrap(facets=~ID,scales="free",)+
  scale_x_discrete(breaks = c(0,100,250,500,750,1000))+ # set to values reasonable for the number of bootstraps you used
  xlab("Number bootstraps")+
  ylab("Mean estimated counts")+
  ggtitle("Genes with the greatest difference between kallisto estimated counts and bootstrap counts")+
  gentheme+theme(axis.text.y=element_blank(),axis.title.x = element_text(size=15),
                 axis.title.y = element_text(size=15),plot.title = element_text(size=17),
                 axis.text = element_text(size=11))


diffhist<-diff[diff!=0] # remove 0 values

## histogram of factor difference between est_counts in abundance + full bootstrap mean
ggplot(as.data.frame(diffhist),aes(x=diffhist))+geom_histogram(fill="steelblue4",colour="white")+
  coord_cartesian(expand=0)+
  xlab("Factor")+
  ylab("Frequency")+
  ggtitle("Difference between mean of bootstrap estimated counts \nand kallisto estimated count")+
  gentheme+theme(axis.title.x = element_text(size=14),axis.title.y = element_text(size=14),
                 plot.title = element_text(size=16),axis.text = element_text(size=10))
