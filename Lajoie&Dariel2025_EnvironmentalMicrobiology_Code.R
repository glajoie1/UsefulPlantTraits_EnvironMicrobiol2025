library(reshape2)
library(stringr)
library(readxl)
library(ggplot2)
library(forcats)
library(dplyr) 
library(rstatix)
library(Matrix)

## Importing data

# Extraction data
my_data_init <- read_xlsx("~/GitHub_SuppTab1.xlsx")
my_data_init<-my_data_init[,c(1:14)]
my_data<-my_data_init

# Merge categories
my_data[which(my_data$Microbial_taxo%in%c('Not specified')),"Microbial_taxo"]<-'All'
my_data[which(my_data$Microbial_taxo%in%c('Bacteria & Archaea')),"Microbial_taxo"]<-'Bacteria'
my_data[which(my_data$Microbial_taxo%in%c('Bacteria & Fungi & Nematodes')),"Microbial_taxo"]<-'Bacteria & Fungi'
my_data[which(my_data$Microbial_taxo%in%c('Bacteria & Fungi')),"Microbial_taxo"]<-'All'
my_data[which(my_data$Microbial_compartment=='NA'),'Microbial_compartment']<-'Mixed (Epiphytes & Endophytes)'
my_data[which(my_data$Microbial_compartment=='Mixed (Epiphytes & Endophytes)'),'Microbial_compartment']<-'Mixed'

## Basic statistics on no-taxo dataset
my_data_notaxo<-my_data[-which(my_data$Plant_trait_category=='Taxonomy'),]
table(my_data_notaxo$Microbial_fitness_measure)

# How many observations
nrow(my_data_notaxo) # 1115

# How many studies overall
length(unique(my_data_notaxo$Study_ID)) # 116

# Plant traits 
# Define categories to summarize?
sort(unique(my_data_notaxo$Plant_trait_category)) # 116
sort(unique(my_data_notaxo$Plant_trait)) # 231

# Number of plant species
plant.spp<-my_data_notaxo$Plant_species
plant.spp<-str_split(plant.spp, ', ')
plant.spp<-unlist(plant.spp)
sort(unique(plant.spp))

# Number of unique plant species
# Refine this to remove cultivars and varieties, or possibly add varieties documented in the table
plant.spp.uni<-unique(plant.spp)
plant.spp.uni<-sort(plant.spp.uni) # 
plant.spp.uni[grep(' +[?|x]+ ', plant.spp.uni)]<-gsub(' ', '.', plant.spp.uni[grep(' +[?|x]+ ', plant.spp.uni)]) 
plant.spp.uni[grep('(\\S+\\s+\\S+)', plant.spp.uni, perl=TRUE)]<-regmatches(plant.spp.uni, regexpr('(\\S+\\s+\\S+)', plant.spp.uni, perl=TRUE))
plant.spp.uni<-unique(plant.spp.uni)
length(plant.spp.uni)-2 # 611 plant species or hybrids

# Distribution across compartments and microorg.
# split experimental setting column
my_data_notaxo<-cbind(my_data_notaxo, colsplit(my_data_notaxo$Experimental_setting, ':', c('experiment','system')))

# Data description (overall)
# Microbial characterization
mc<-as.data.frame(table(my_data_notaxo$Microbial_characterization))
mc$variable<-'Microbial characterization'

# Microbial taxonomy
mt<-as.data.frame(table(my_data_notaxo$Microbial_taxo))
mt$variable<-'Microbial taxonomy'

# Microbial compartment
mcomp<-as.data.frame(table(my_data_notaxo$Microbial_compartment))
mcomp$variable<-'Compartment'

# Experiment
ex<-as.data.frame(table(my_data_notaxo$experiment))
ex$variable<-'Experiment type'

# System type
sys<-as.data.frame(table(my_data_notaxo$system))
sys$variable<-'System type'

# Within vs among plant species
spp<-as.data.frame(table(my_data_notaxo$Plant_scale_simp))
spp$variable<-'Plant scale'

# Pop vs community
sca<-as.data.frame(table(my_data_notaxo$Microbial_scale))
sca$variable<-'Microbial scale'

# Characterization
car<-rbind(mc,mt,mcomp,ex,sys, sca, spp)

threshold <- 50

out <- car %>% 
  group_by(variable) %>%
  count(Var1 = fct_collapse(as.character(Var1), Other = unique(as.character(Var1)[Freq < threshold])),  
        wt = Freq)
out
out$variable<-as.factor(out$variable)

out<-out[order(out$n, decreasing=T),]
out$Var1<-factor(out$Var1, levels=c(unique(out$Var1)[-which(unique(out$Var1)=='Other')],as.factor('Other')))
out$variable<-factor(out$variable, levels=c('Experiment type','System type','Compartment','Plant scale','Microbial taxonomy','Microbial scale','Microbial characterization'))

# Plot : Data distribution
ggplot(out, aes(x=Var1, y=n))+
  geom_bar(position='dodge',stat='identity')+
#  scale_fill_viridis(discrete = T, option = "E") +
  facet_wrap(~variable, scales='free', ncol=2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.text=element_text(size=10.5))

  
# How many unique traits*study per Microbial taxo x Compartment?
un.trt<-my_data_notaxo[!duplicated(my_data_notaxo[,c('Study_ID','Plant_trait_category','Microbial_taxo','Microbial_compartment')]),] 
un.trt$Line_ID<-c(1:nrow(un.trt))

# How often was each trait measured?
# By total number of observations
trt.study<-aggregate(Line_ID~Plant_trait_category, un.trt, length) # +Microbial_taxo+Microbial_compartment
trt.study<-trt.study[order(trt.study$Line_ID, decreasing=T),] 
# OR
# By number of studies
trt.study<-aggregate(Study_ID~Plant_trait_category, un.trt, FUN=function (x) length(unique(x))) # +Microbial_taxo+Microbial_compartment
trt.study<-trt.study[order(trt.study$Study_ID, decreasing=T),]
trt.study$Plant_trait_category<-factor(trt.study$Plant_trait_category, levels=trt.study[order(trt.study$Study_ID, decreasing=T),'Plant_trait_category'])
trt.study<-trt.study[order(trt.study$Study_ID, decreasing=T),] 

## Which traits are having impacts on microbial communities significantly more often than expected by chance (50/50) ?

# Only keep traits that have been measured in at least 5 studies & 10 observations
my_data_notaxo$Line_ID<-rownames(my_data_notaxo)
my_data_sub_spp<-my_data_notaxo

sum.study<-aggregate(Study_ID~Plant_trait_category, my_data_sub_spp, FUN=function(x) length(unique(x))) # change here to input all 3
sum.study<-sum.study[which(sum.study$Study_ID>=5),]
sum.line<-aggregate(Line_ID~Plant_trait_category, my_data_sub_spp, length) # change here to input all 3
sum.line<-sum.line[which(sum.line$Line_ID>=10),]
sum.line<-sum.line[order(sum.line$Line_ID),]

my_data_sub_spp<-my_data_sub_spp[which(my_data_sub_spp$Plant_trait_category%in%intersect(sum.study$Plant_trait_category, sum.line$Plant_trait_category)),]
my_data_sub_spp$Effect_YN<-ifelse(my_data_sub_spp$Effect_YN=="Yes",1,0)

length(unique(my_data_sub_spp$Plant_trait_category))
length(unique(my_data_sub_spp$Study_ID))

# Effect of taxonomy
tabtax<-rev(table(my_data[which(my_data$Plant_trait_category=='Taxonomy'),'Effect_YN']))

# Add Plant taxonomy effect
binom.test(tabtax)

# Plot all variables

# All data
yestab<-table(my_data_sub_spp$Plant_trait_category,my_data_sub_spp$Effect_YN)
yestab.df<-as.data.frame(yestab)

taxx<-as.data.frame(table(my_data[which(my_data$Plant_trait_category=='Taxonomy'),'Effect_YN']))
taxx<-cbind(c('Taxonomy','Taxonomy'),taxx)
colnames(taxx)<-c('Var1','Var2','Freq')
taxx$Var2<-ifelse(taxx$Var2=="Yes",1,0)

yestab.df<-rbind(yestab.df, taxx)

# Plot variables
ggplot(yestab.df, aes(x=Var2, y=Var1, label=Freq, color=Var2))+
  geom_point(aes(size=Freq))+
  scale_size(range = c(0, 12))+
  geom_text(col='white', size=3)+
  scale_y_discrete(limits=rev)+
  scale_colour_manual(values=c('black','darkgreen'))+
  theme_bw()



# Binomial tests for each trait (all data)

yestab2<-yestab[,c(2,1)]

dm<-NULL
for(i in 11:nrow(yestab2)){ # nrow(yestab2)
  dm<-rbind(dm,c(rownames(yestab2)[i],binom.test(yestab2[i,],sum(yestab2[i,]), alternative='greater'),yestab2[i,]))
  print(c(rownames(yestab2)[i],binom.test(yestab2[i,],sum(yestab2[i,]), alternative='greater'),yestab2[i,]))
}

dmm<-as.data.frame.matrix(dm)


# Binomial tests for each trait (median effect per trait per study)

# Median 
my_data_sub_med<-aggregate(Effect_YN~Study_ID+Plant_trait_category+Microbial_taxo+Microbial_compartment+Plant_scale_simp, my_data_sub_spp, FUN=function (x) round(median(x)+0.01)) # Median # min # max # moyenne
yestab.med<-table(my_data_sub_med$Plant_trait_category,my_data_sub_med$Effect_YN)


# Pairwise comparisons

yestab2.med<-yestab.med[,c(2,1)]

dm.med<-NULL
for(i in 11:nrow(yestab2.med)){ # 
  dm.med<-rbind(dm.med,c(rownames(yestab2.med)[i],binom.test(yestab2.med[i,],sum(yestab2.med[i,]), alternative='greater'),yestab2.med[i,]))
  print(c(rownames(yestab2.med)[i],binom.test(yestab2.med[i,],sum(yestab2.med[i,]), alternative='greater'),yestab2.med[i,]))
  }
dmm.med<-as.data.frame.matrix(dm.med)

###############################################################################

# Testing differences among factor levels


## General model
table(my_data_sub_spp$Microbial_taxo, my_data_sub_spp$Microbial_compartment, my_data_sub_spp$Plant_scale)


# AMONG TAXA
my_data_amongspp_sub<-my_data_sub_spp[which(my_data_sub_spp$Microbial_taxo%in%c('Fungi','Bacteria')),]

# Only traits that have at least 5 studies for each microbial taxo
tt0<-aggregate(my_data_amongspp_sub$Study_ID~my_data_amongspp_sub$Plant_trait_category+my_data_amongspp_sub$Microbial_taxo, FUN=function(x) length(unique(x)))
tt2<-tt0[which(tt0$`my_data_amongspp_sub$Study_ID`>=5),]
tt3<-table(tt2$`my_data_amongspp_sub$Plant_trait_category`)
tt3<-tt3[which(tt3==2)]

# Only traits that have at least 10 observations per taxa
tt<-as.data.frame.matrix(table(my_data_amongspp_sub$Plant_trait_category, my_data_amongspp_sub$Microbial_taxo)) #, my_data_amongspp_sub$Study_ID

# Select traits that fulfill both criteria
tt4<-tt[which(tt$Bacteria>=10&tt$Fungi>=10),]
tt4<-tt4[which(rownames(tt4)%in%names(tt3)),]

# All data
mds<-my_data_amongspp_sub[which(my_data_amongspp_sub$Plant_trait_category%in%rownames(tt4)),]


unique(mds$Plant_trait_category)

for (i in unique(mds$Plant_trait_category)){
  test<-mds[which(mds$Plant_trait_category==i),]
  test.t<-as.data.frame.matrix(table(test$Microbial_taxo,test$Effect_YN))
  test.t<-as.matrix(test.t[,c(2,1)])
  print(i)
  print(test.t)
  print(binom.test(colSums(test.t), alternative='greater'))
  print(fisher_test(test.t, detailed=T))
}

## Median

mds2<-aggregate(Effect_YN~Study_ID+Plant_trait_category+Microbial_taxo+Microbial_compartment+Plant_scale_simp, mds, FUN=function (x) round(median(x)+0.01))

unique(mds2$Plant_trait_category)

for (i in unique(mds2$Plant_trait_category)){
  test<-mds2[which(mds2$Plant_trait_category==i),]
  test.t<-as.data.frame.matrix(table(test$Microbial_taxo,test$Effect_YN))
  test.t<-as.matrix(test.t[,c(2,1)])
  print(i)
  print(test.t)
  print(binom.test(colSums(test.t), alternative='greater'))
  print(fisher_test(test.t, detailed=T))
}





# AMONG PLANT SCALES
tt0<-aggregate(my_data_sub_spp$Study_ID~my_data_sub_spp$Plant_trait_category+my_data_sub_spp$Plant_scale_simp, FUN=function(x) length(unique(x)))

tt2<-tt0[which(tt0$`my_data_sub_spp$Study_ID`>=5),]
tt3<-table(tt2$`my_data_sub_spp$Plant_trait_category`)
tt3<-tt3[which(tt3==2)]

tt<-as.data.frame.matrix(table(my_data_sub_spp$Plant_trait_category, my_data_sub_spp$Plant_scale_simp)) #, my_data_amongspp_sub$Study_ID

# Select traits that have at least 5 studies and 10 observations per taxa 
tt4<-tt[which(tt[,1]>=10&tt[,2]>=10),] # Four traits
tt4<-tt4[which(rownames(tt4)%in%names(tt3)),]

# Test impact of taxon for these traits
mds<-my_data_sub_spp[which(my_data_sub_spp$Plant_trait_category%in%rownames(tt4)),] 


# All data
unique(mds$Plant_trait_category)

for (i in unique(mds$Plant_trait_category)){
  test<-mds[which(mds$Plant_trait_category==i),]
  test.t<-as.data.frame.matrix(table(test$Plant_scale_simp,test$Effect_YN))
  test.t<-as.matrix(test.t[,c(2,1)])
  print(i)
  print(test.t)
  print(binom.test(colSums(test.t), alternative='greater'))
  print(fisher_test(test.t, detailed=T))
}

## Median

mds2<-aggregate(Effect_YN~Study_ID+Plant_trait_category+Microbial_taxo+Microbial_compartment+Plant_scale_simp, mds, FUN=function (x) round(median(x)+0.01))

unique(mds2$Plant_trait_category)

for (i in unique(mds2$Plant_trait_category)){
  test<-mds2[which(mds2$Plant_trait_category==i),]
  test.t<-as.data.frame.matrix(table(test$Plant_scale_simp,test$Effect_YN))
  test.t<-as.matrix(test.t[,c(2,1)])
  print(i)
  print(test.t)
  print(binom.test(colSums(test.t), alternative='greater'))
  print(fisher_test(test.t, detailed=T))
}











