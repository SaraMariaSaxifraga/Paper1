### beta diversity and null models (?) ####


# Script for performing null model analysis on functional diversity data
library(vegan)
install.packages("simba" )
library(FD)

source("~/Documents/Data/R_scripts/species_matrix.R")
source("~/Documents/Data/R_scripts/FD_calc.R")
### load tables from the species and the trait matrix ####
B
Tm

##### devide the dataset into two datasets: calcareous and acidic ####

calc<-myData %>% filter(Bedrock=="calcareous")
acidic<-myData %>% filter(Bedrock=="acidic")

B.subset<-B %>% filter(plotID %in% myData$plotID)
B.calc<- B.subset %>% filter(plotID %in% calc$plotID)
B.acidic<- B.subset %>% filter(plotID %in% acidic$plotID)

###########################################
###########################################
#the community matrix might have to be presence absence 
df.calc.pa <- B.calc[,-1]
x.pa.calc<-df.calc.pa %>% mutate_if(is.numeric, ~1 * (. != 0)) #x.pa means x.presence absence

df.acidic.pa <- B.acidic[,-1]
x.pa.acidic<-df.acidic.pa %>% mutate_if(is.numeric, ~1 * (. != 0)) #x.pa means x.presence absence

# delete all the zero columns 

#No I delete all species with zero columns
#put all zero columns in a vector 
b.calc<-x.pa.calc[,colSums(x.pa.calc!= 0) == 0]  #72 species 

bn.calc<-colnames(b.calc)

for (i in 1: length(bn.calc)) {
  
  x.pa.calc<-x.pa.calc[, !names(x.pa.calc) %in% bn.calc[i]]
  
}

#-------------------------------------------------------------#
b.acidic<-x.pa.acidic[,colSums(x.pa.acidic!= 0) == 0]  #58 species 

bn.acidic<-colnames(b.acidic)

for (i in 1: length(bn.acidic)) {
  
  x.pa.acidic<-x.pa.acidic[, !names(x.pa.acidic) %in% bn.acidic[i]]
  
}

#------------------------------------------------------------------------#
#now I have to adapt the trait matrix 
colnames.calc<-colnames(x.pa.calc)
traits.calc <- Y %>% filter (species %in% colnames.calc)
traits.calc.m <- traits.calc[,-1]
rownames(traits.calc.m) <- traits.calc$species

colnames.acidic<-colnames(x.pa.acidic)
traits.acidic <- Y %>% filter (species %in% colnames.acidic)
traits.acidic.m <- traits.acidic[,-1]
rownames(traits.acidic.m) <- traits.acidic$species

#### Perform null model analysis using vegan::oecosimu function and the quasiswap algorithm ####

fd_fric.calc<-function(x)
{f<-dbFD(Tm,x.pa,corr="cailliez",calc.FGR=F,clust.type = "ward")
print(f$FRic)
}
Fric<-fd_fric(x.pa)

cbind(Fric, FD$FRic)

FRic<-oecosimu(x, fd_fric, "quasiswap", nsimul = 20)
boxplot(FRic$oecosimu$z~EnvData$Alter[1:40])

# -------------------------------------------- # 

# -------------------------------------------- # 

# -------------------------------------------- # 
#### Fdis ####

# ------- ACIDIC -------------------#
fd_dis.acidic<-function(x)
{f.acidic <- dbFD(traits.acidic.m, x, w.abun=TRUE, stand.x=T, corr="lingoes", calc.FGR=F)
print(f.acidic$FDis)
}

Fdis.acidic<-fd_dis.acidic(x.pa.acidic)

FDis.acidic.oeco<-oecosimu(x.pa.acidic, fd_dis.acidic, "curveball", nsimul = 20)

#boxplot(FDis$oecosimu$simulated~all_plot_data$Bedrock)

# ------- CALC -------------------#

fd_dis.calc<-function(x)
{f.calc <- dbFD(traits.calc.m, x, w.abun=TRUE, stand.x=T, corr="lingoes", calc.FGR=F)
print(f.calc$FDis)
}

Fdis.calc<-fd_dis.calc(x.pa.calc)

FDis.calc.oeco<-oecosimu(x.pa.calc, fd_dis.calc, "quasiswap",nsimul = 20)


#20 simulated FDis are under FDis$oecosimu$simulated 

FDis.calc.oeco$oecosimu$simulated
FDis.acidic.oeco$oecosimu$simulated

FDis.calc.oeco$oecosimu$statistic
var(FDis.calc.oeco$oecosimu$statistic)   # 0.0003764175

FDis.acidic.oeco$oecosimu$statistic
var(FDis.acidic.oeco$oecosimu$statistic) # 0.0003484844

View(FDis.acidic.oeco$oecosimu$simulated)
View(FDis.calc.oeco$oecosimu$simulated)

###########################################
###########################################

###### Now the entire procedure not with presence absence data, but with abundance data ####

 
df.calc.ab <- B.calc[,-1]  #the ab stands for abundance


df.acidic.ab <- B.acidic[,-1]


# delete all the zero columns 

#No I delete all species with zero columns
#put all zero columns in a vector 
b.calc<-df.calc.ab[,colSums(df.calc.ab!= 0) == 0]  #72 species 

bn.calc<-colnames(b.calc)

for (i in 1: length(bn.calc)) {
  
  df.calc.ab<-df.calc.ab[, !names(df.calc.ab) %in% bn.calc[i]]
  
}
#-------------------------------------------------------------#
b.acidic<-df.acidic.ab[,colSums(df.acidic.ab!= 0) == 0]  #58 species 

bn.acidic<-colnames(b.acidic)

for (i in 1: length(bn.acidic)) {
  
  df.acidic.ab<-df.acidic.ab[, !names(df.acidic.ab) %in% bn.acidic[i]]
  
}

#------------------------------------------------------------------------#
#now I have to adapt the trait matrix 
colnames.calc<-colnames(df.calc.ab)
traits.calc <- Y %>% filter (species %in% colnames.calc)
traits.calc.m <- traits.calc[,-1]
rownames(traits.calc.m) <- traits.calc$species

colnames.acidic<-colnames(df.acidic.ab)
traits.acidic <- Y %>% filter (species %in% colnames.acidic)
traits.acidic.m <- traits.acidic[,-1]
rownames(traits.acidic.m) <- traits.acidic$species

#### Perform null model analysis using vegan::oecosimu function and the quasiswap algorithm ####

# -------------------------------------------- # 

# -------------------------------------------- # 

# -------------------------------------------- # 
#### Fdis ####

# ------- ACIDIC -------------------#
fd_dis.acidic<-function(x)
{f.acidic <- dbFD(traits.acidic.m , x, w.abun=TRUE, stand.x=T, corr="lingoes", calc.FGR=FALSE,)
print(f.acidic$FDis)
}

Fdis.acidic<-fd_dis.acidic(df.acidic.ab)

#torun
FDis.acidic.oeco<-oecosimu(df.acidic.ab, fd_dis.acidic, "swsh_samp", burnin=100, nsimul = 1000)

#boxplot(FDis$oecosimu$simulated~all_plot_data$Bedrock)

# ------- CALC -------------------#

fd_dis.calc<-function(x)
{f.calc <- dbFD(traits.calc.m , x, w.abun=TRUE, stand.x=T, corr="lingoes", calc.FGR=F)
print(f.calc$FDis)
}

Fdis.calc<-fd_dis.calc(df.calc.ab)

#torun
FDis.calc.oeco<-oecosimu(df.calc.ab, fd_dis.calc, "swsh_samp", burnin=100, nsimul = 1000)


#20 simulated FDis are under FDis$oecosimu$simulated 

FDis.calc.oeco$oecosimu$simulated
FDis.acidic.oeco$oecosimu$simulated

FDis.calc.oeco$oecosimu$statistic
var(FDis.calc.oeco$oecosimu$statistic)   # 0.001324794

FDis.acidic.oeco$oecosimu$statistic
var(FDis.acidic.oeco$oecosimu$statistic) # 0.001265174

FDis.acidic.oeco$statistic
acidic$FDis

View(FDis.acidic.oeco$oecosimu$simulated)
View(FDis.calc.oeco$oecosimu$simulated)

FDis.acidic.oeco$oecosimu$z

length(which(FDis.acidic.oeco$oecosimu$pval<0.05))
length(which(FDis.calc.oeco$oecosimu$pval<0.05))

boxplot(FDis.acidic.oeco$oecosimu$simulated~acidic$DLI_cv)
boxplot(FDis.calc.oeco$oecosimu$simulated~calc$DLI_cv)

boxplot(FDis.acidic.oeco$oecosimu$simulated~acidic$DLI_cv)
boxplot(FDis.calc.oeco$oecosimu$simulated~calc$DLI_cv)

z_acidic<-as.numeric(FDis.acidic.oeco$oecosimu$z)
z_calc<-as.numeric(FDis.calc.oeco$oecosimu$z)

z_calc<-as.data.frame(z_calc)
z_calc$Bedrock<-"calcareous"
colnames(z_calc)[1]<-"z"
z_acidic<-as.data.frame(z_acidic)
z_acidic$Bedrock<-"acidic"
colnames(z_acidic)[1]<-"z"


z_scores<-rbind(z_calc,z_acidic)

ggplot(data=z_scores)+
 # geom_boxplot(aes(y=z, x= Bedrock ))+
  geom_boxplot(aes(y= z, x= myData$Bedrock))

summary(aov(z~Bedrock,data=z_scores)) #significant!  0.00818 **



##### FDis null model plot #####
## make a summary 
z_scores_sum<-data_summary(data=z_scores, varname = "z",  # this is the summary of FDis values across plotID and Bedrock
                           groupnames= c("Bedrock"))

View(FDis.acidic.oeco$oecosimu$simulated)
simu_acidic<-as.data.frame(FDis.acidic.oeco$oecosimu$simulated)
simu_calc<-as.data.frame(FDis.calc.oeco$oecosimu$simulated)

simu_acidic$plotID<-acidic$plotID
simu_calc$plotID<-calc$plotID

simu_calc$Bedrock<-calc$Bedrock
simu_acidic$Bedrock<-acidic$Bedrock

simu<-rbind(simu_acidic,simu_calc)

simu_long<-data.table::melt(setDT(simu),id.vars=c("plotID","Bedrock"), variable.name="FDis_simu")

source("~/Documents/Data/Standard_protocoll_startup.R")

#alternative version to reduce the dataframe   (df.ds = df.data_summary)
df.ds<- data_summary(data=simu_long, varname = "value",    # this is the summary of bedrock independent of plotID
                            groupnames= c("Bedrock"))

 detach("package:dplyr", unload = TRUE)
simu_summary<- data_summary(data=simu_long, varname = "value",  # this is the summary of FDis values across plotID and Bedrock
                            groupnames= c("Bedrock","plotID"))

df.ds2<- data_summary(data=simu_summary, varname = "value",  #this is the summary of Bedrock, when plotID was already summarized
                     groupnames= c("Bedrock"))



FDis_simulated<-ggplot()+
  geom_jitter(data=simu_summary,aes(x=Bedrock, y=value),alpha=0.2,width=0.3)+
  geom_point(stat="summary", fun="mean",data=myData, aes(x=Bedrock, y=FDis,stat=mean),color="red")+
  #stat_summary()
  geom_errorbar(data=df.ds2, aes(ymin = value-2*se , ymax = value+2*se,x=Bedrock,color=Bedrock),width=.2)+
  geom_point(data=df.ds2,aes(x=Bedrock, y=value,color=Bedrock),size=2)+
  ylab("FDis(simulated and actual")+
  scale_color_manual(values=c("#1f78b4", "#4dac26"))+
  scale_fill_manual(values=c("#1f78b4", "#4dac26"))+
theme_set(
  theme_bw() +
    theme(
      # panel.background = element_rect(fill = NA, colour = "black"),
      # strip.background = element_rect(fill = "white", colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
      # panel.grid.major = element_blank(),  # remove major-grid
      # panel.grid.minor = element_blank(),  # remove minor-grid
      axis.ticks = element_line(color = "black", size = 0.75),
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size=15),
      axis.title.y = element_text(size=15),
      legend.position = "none",
      #legend.text = element_text(size = 12),
      #legend.title = element_text(size = 15),
      # legend.background = element_rect(fill = NA),
      strip.background = element_rect(fill = "#FFFFFF"))) # size=1.5, linetype="solid", face = "bold.italic" for facet background

pdf(file="FDis_simulated.pdf",width = 15, height = 6)
FDis_simulated
dev.off() 
#### z-scores plot #####  

#z-score :RaoSES (z) is observed Rao BA for each plot expressed relative to random expectation 
#using the standardized effect size (SES) 


z_score.bedrock<-ggplot()+
  geom_point(data=z_scores,aes(x=z,y=Bedrock),alpha=0.3)+
  geom_errorbar(data=z_scores_sum,aes(xmin = z-2*se , xmax = z+2*se,y=Bedrock,color=Bedrock),width=.2)+
  geom_point(data=z_scores_sum,aes(y=Bedrock, x= z, color=Bedrock),size=3)+
  scale_color_manual(values=c("#1f78b4", "#4dac26"))+
  scale_fill_manual(values=c("#1f78b4", "#4dac26"))+
  xlim(-4,4)+
  theme_set(
    theme_bw() +
      theme(
        # panel.background = element_rect(fill = NA, colour = "black"),
       # strip.background = element_rect(fill = "white", colour = "black"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.75),
        # panel.grid.major = element_blank(),  # remove major-grid
        # panel.grid.minor = element_blank(),  # remove minor-grid
        axis.ticks = element_line(color = "black", size = 0.75),
        axis.text = element_text(size = 15, color = "black"),
        axis.title = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.position = "none",
        #legend.text = element_text(size = 12),
        #legend.title = element_text(size = 15),
       # legend.background = element_rect(fill = NA),
        strip.background = element_rect(fill = "#FFFFFF"))) # size=1.5, linetype="solid", face = "bold.italic" for facet background

z_score.bedrock

pdf(file="z_scores_bedrock.pdf",width = 15, height = 6)
z_score.bedrock
dev.off() 





# nm<-simulate(nullmodel(df.calc.ab,"swsh_samp"),99)
# oecosimu(nm,fd_dis.calc)

#### OLD ####
# FDis_true<-FDis$oecosimu$statistic
# View(FDis_true)
# 
# FDis_simu<-FDis$oecosimu$simulated
# FDis$oecosimu$simulated[]
# 
# all_plot_data<-all_plot_data[order(all_plot_data$plotID),] #really important !!
# all_plot_data[,]
# 
# FDis_true.df<-as.data.frame(FDis_true)
# FDis_true.df$plotID <- all_plot_data$plotID #add the plotID column to the true Fdis values
# FDis_true.df<-FDis_true.df %>% filter(plotID %in% myData$plotID)
# 
# FDis_simu.df<-as.data.frame(FDis_simu)
# FDis_simu.df$plotID <- all_plot_data$plotID #add the plotID column to the true Fdis values
# FDis_simu.df<-FDis_simu.df %>% filter(plotID %in% myData$plotID)
# 
# boxplot(FDis_true.df$FDis_true~myData$Bedrock)
# 
# ggplot()+
#   geom_jitter()


# --------------------------------------------------- #

## Significance test using Null model communities.
## The current choice fixes both species and site totals.

# --------------------------------------------------- #

# --------------------------------------------------- #
