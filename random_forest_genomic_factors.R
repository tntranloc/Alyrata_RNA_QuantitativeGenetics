library(ranger)
library(tidyverse)
library(patchwork)

theme_set(theme_minimal(base_size = 20)%+replace% theme(panel.grid.minor.y = element_blank()))


#################RUN the Random Forest###############################
names(ase)=make.names(names(ase))
motifsnames=make.names(motifsnames)


features=c("GenesDEG","GenesSweep","GenesSPsweep","GenesPLsweep","Fst","Density","transcript.length","Length","Exons","total") #


ase_rf=ranger(data=ase[,colnames(ase) %in% c("ASE",motifsnames,features)],formula = ASE~.,num.trees = 500,mtry=200,importance = "impurity_corrected")
ase_rf


save(ase_rf,file="ase_rf.RData")


cor.pred <- sum(ase$ASE==ase_rf$predictions) / length(ase$ASE)


barplot(sort(importance(ase_rf))[1:50],horiz = T)
save.image("random_forest_ase.RData")


importance=as.data.frame(importance_pvalues(ase_rf))
importance=as.data.frame(importance_pvalues(ase_rf,method = "altmann",data=ase[,colnames(ase) %in% c("ASE",motifsnames,features)],formula = ASE~.,num.threads=7,num.permutations = 100))
importance=importance[order(importance[,1],decreasing = T),]
importance$motif=rownames(importance)
importance$sig=importance$pvalue<0.05


save(importance,file="ase_rf_importance.RData")


summary(importance)


svg("ASE_importance.svg",height=12,width=10)
ggplot(data=importance[1:50,],aes(x=reorder(motif,importance),y=importance,fill=sig))+geom_bar(stat="identity")+coord_flip()+xlab("factor")+
  theme(panel.grid.major.y = element_blank(),legend.position = "top")+scale_fill_manual(name="significant?",values=c("grey25","#BADA55"))
dev.off()


#check the top 5 significant factors
ggplot(ase, aes(x=ASE, y=transcript.length)) + geom_boxplot()
tapply(ase$transcript.length, ase$ASE, mean)
tapply(ase$transcript.length, ase$ASE, median)


ggplot(ase, aes(x=ASE, y=Length)) + geom_boxplot()
tapply(ase$Length, ase$ASE, mean)
tapply(ase$Length, ase$ASE, median)


####Rerun only for TFs#####
ase2 <- ase[,c(1,2,17:438)]


names(ase2)=make.names(names(ase2))
motifsnames=make.names(motifsnames)


ase_rf2=ranger(data=ase2[,colnames(ase2) %in% c("ASE",motifsnames,"total")],formula = ASE~.,num.trees = 500,mtry=200,importance = "impurity_corrected")
ase_rf2


save(ase_rf2,file="ase_rf2.RData")


cor.pred <- sum(ase2$ASE==ase_rf2$predictions) / length(ase2$ASE)


barplot(sort(importance(ase_rf2))[1:50],horiz = T)
save.image("random_forest_ase2.RData")


importance2=as.data.frame(importance_pvalues(ase_rf2))
importance2=as.data.frame(importance_pvalues(ase_rf2,method = "altmann",data=ase2[,colnames(ase2) %in% c("ASE",motifsnames,"total")],formula = ASE~.,num.threads=7,num.permutations = 100))
importance2=importance2[order(importance2[,1],decreasing = T),]
importance2$motif=rownames(importance2)
importance2$sig=importance2$pvalue<0.05


save(importance2,file="ase_rf_importance2.RData")


summary(importance2)


svg("ASE_importance2.svg",height=12,width=10)
ggplot(data=importance2[1:50,],aes(x=reorder(motif,importance),y=importance,fill=sig))+geom_bar(stat="identity")+coord_flip()+xlab("factor")+
  theme(panel.grid.major.y = element_blank(),legend.position = "top")+scale_fill_manual(name="significant?",values=c("grey25","#BADA55"))
dev.off()


load("C:/Users/mzt5590/Downloads/ASEthaliana/ASEthaliana/deg_ASEthal_rf_importance.RData")
load("ase_rf_importance2.RData")


tmp1 <- importance[importance$sig == TRUE,]
tmp2 <- importance2[importance2$sig ==TRUE,]
intersect(tmp1[,3], tmp2[,3])

