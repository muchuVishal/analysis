
install.packages("car")
library(car)



# put data from all chr together, to filter the variant files, take only the variants annotated in ClinVar and the columns with the patho tools: 

######################
# for testing this can be also done in R only for chr21:
c<-read.delim("/external.data/Variants/dbNSFP/dbNSFPv3.1/dbNSFP3.1a/dbNSFP3.1a_variant.chr21", header=TRUE)

#head(c)

#select Columns with the variant possition, RA, AA, clinvar singnificance, and the prediction tools. 
#colnames(c)
pred<-c[,c(1,2,3,4,139,7,26,33,38,42,49,52,55,58,60,63,67,71,74,75)]

#select the prediction values for those variants present in clinvar

pred.clinvar<-pred[which(pred$clinvar_clnsig != "."),]
##########################

pred.clinvar<-read.delim("/external.data/Variants/dbNSFP/dbNSFPv3.1/dbNSFP3.1a/dbNSFP3.1a_variant.ClinVar.PathoTools.txt", sep="\t")

head(pred.clinvar)
pred.clinvar<-pred.clinvar[,c(1:4,20,5,6:19)]


name<-c("Chr", "Pos","R","A","clinvar_clnsig","rs_dbSNP144","SIFT", "Polyphen2_HVAR","LRT_pred","MutationTaster_pred","MutationAssessor_pred",
         "FATHMM", "PROVEAN","VEST3","CADD", "DANN", "fathmm.MKL","MetaSVM","MetaLR_pred","Reliability_index")
name
colnames(pred.clinvar)<-name

pred.clinvar[is.na(pred.clinvar)]<-"."

head(pred.clinvar,n = 10 )
summary(pred.clinvar)


#I am going to replace all the predictions by: disease causing = 1, tolerated=-1, unknown=0;

# SIFT:     Polyphen      LRT_pred  MT      MA      FATHMM  PROVEAN   FATHMM.MKL  MetaSVM   MetaLR
#   D = 1   > 0.447 = 1   D = 1     A = 1   H = 1   D = 1   D = 1     D = 1       D = 1     D = 1
#   T = -1  < 0.447 = -1  N = -1    D = 1   M = 1   T = -1  N = -1    N = -1       T = -1    T = -1
#   . = 0   . = 0         . = 0     N = -1  L= -1   . = 0   . = 0     . = 0       . = 0     . = 0
#                         U = 0     P= -1   N= -1
#                                   . = 0   . = 0
  
colnames(pred.clinvar)

#create a loop to replace letters by 0, -1 and 1 !!!!!!!!!!!!!!!!!!!!

list.pred<-list("")
list.pred
new.x<-list("")
z=1
for (i in c(7,9:13,17:19))
  { 
  
  pred.clinvar[,i]
  x<-strsplit(as.character(pred.clinvar[,i]),";")
  head(x, n=10)
  length(x)


#I use the recode function to select the most pathogenic score, when there are multiple predictions


        for (r in 1:length(x)){
          
          a<-x[r][[1]][
            which.max(recode(x[r][[1]],"c('T','N','P','L')='-1';c('D','A','H','M')='1';c('.','U')='0'"))
            ]
          #x[10][[1]][2]
          new.x[r]<-unlist(a)
          
          #head(a)
          #head(new.x)
          
                   
          #length(new.x)
          #head(new.x)
          #new.x[6]
          #name[7]
          
        }
  
  
  list.pred[z]<-list(new.x)
  z<-z+1
  }

#end of the for loop

#give a name to each list
names(list.pred)<-name[c(7,9:13,17:19)]
summary(list.pred)

#I replace the letters describing the pathogenicity predictions by numbers (-1, 0, 1):

list.pred$SIFT<-recode(list.pred$SIFT,"c('T','N','P','L')='-1';c('D','A','H','M')='1';c('.','U')='0'")

list.pred$LRT_pred<-recode(list.pred$LRT_pred,"c('T','N','P','L')='-1';c('D','A','H','M')='1';c('.','U')='0'")

list.pred$MutationTaster_pred<-recode(list.pred$MutationTaster_pred,"c('T','N','P','L')='-1';c('D','A','H','M')='1';c('.','U')='0'")

list.pred$MutationAssessor_pred<-recode(list.pred$MutationAssessor_pred,"c('T','N','P','L')='-1';c('D','A','H','M')='1';c('.','U')='0'")

list.pred$FATHMM<-recode(list.pred$FATHMM,"c('T','N','P','L')='-1';c('D','A','H','M')='1';c('.','U')='0'")

list.pred$PROVEAN<-recode(list.pred$PROVEAN,"c('T','N','P','L')='-1';c('D','A','H','M')='1';c('.','U')='0'")

list.pred$fathmm.MKL<-recode(list.pred$fathmm.MKL,"c('T','N','P','L')='-1';c('D','A','H','M')='1';c('.','U')='0'")

list.pred$MetaSVM<-recode(list.pred$MetaSVM,"c('T','N','P','L')='-1';c('D','A','H','M')='1';c('.','U')='0'")

list.pred$MetaLR_pred<- recode(list.pred$MetaLR_pred,"c('T','N','P','L')='-1';c('D','A','H','M')='1';c('.','U')='0'")

         


#convert the list to a dataframe
temp
pred.w<- data.frame(lapply(data.frame((sapply(list.pred, `[`))), unlist))

head(pred.w)

head(pred.clinvar)

summary(pred.w)

#now, for the patho tools with scores:

#I will take the maximum score for this tools, still we need to decide which score we will accept for a predicted pathogenic variant
list.pred2<-list("")
list.pred2
new.x<-list("")
z=1
for (i in c(8,14,15,16))
  { 
  i=8
  pred.clinvar[,i]
  x<-strsplit(as.character(pred.clinvar[,i]),";")
  x[x=="."]<- -9999999999999999999
 
  
  for (r in 1:length(x))
    {
      a<-max((x[r][[1]]))
      new.x[r]<-unlist(a)
    }
  
  new.x
  list.pred2[z]<-list(new.x)
  z<-z+1
}

#give a name to each list
names(list.pred2)<-name[c(8,14,15,16)]
list.pred2


  
#convert the list of lists in data.frame
pred.w2<- data.frame(lapply(data.frame((sapply(list.pred2, `[`))), unlist))

head(pred.w2)
summary(pred.w2)

#head(pred.clinvar)
#put together the new patho scores: 

pred.w.CV<-data.frame(pred.clinvar[,1:5], pred.w, pred.w2)
dim(pred.w.CV)

head(pred.w.CV)
summary(pred.w.CV)

write.table(pred.w.CV, file ="Transformed_prediction_scores.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep="\t")



pred.sum<-data.frame(pred.w.CV$clinvar_clnsig, "sum"=apply(pred.w.CV[,6:14], 1, sum))
head(pred.sum)
boxplot(pred.sum[pred.sum$pred.w.CV.clinvar_clnsig==5,2],
        pred.sum[pred.sum$pred.w.CV.clinvar_clnsig==4,2],
        pred.sum[pred.sum$pred.w.CV.clinvar_clnsig==3,2],
        pred.sum[pred.sum$pred.w.CV.clinvar_clnsig==2,2],
        pred.sum[pred.sum$pred.w.CV.clinvar_clnsig==6,2],
        pred.sum[pred.sum$pred.w.CV.clinvar_clnsig==7,2],
        names = c(5,4,3,2,6,7))







