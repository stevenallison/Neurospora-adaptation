setwd(getwd())
library(ggplot2)
library(nlme)
library(reshape2)
library(gmodels)
library(multcomp)
library(lattice)
library(nlstools)
library(plyr)
library(dplyr)
library(tidyr)
library(xlsx)
library(pracma)
library(car)
library(gridGraphics)
library(gridExtra)
library(cowplot)

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
#################################################################################
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
  #  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  # datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
#################################################################################

#function to convert a 1-d list to a numeric vector, converting NULL to NA
#################################################################################
list.as.vectorNA <- function(dataset)
{
  if(is.null(dataset[[1]]))
  {
    my.vector <- NA
  } else
  {
    my.vector <- dataset[[1]]
  }
  i <- 1
  while(i<length(dataset))
  {
    i <- i+1
    if(is.null(dataset[[i]]))
    {
      my.vector <- c(my.vector,NA)
    } else
    {
      my.vector <- c(my.vector,dataset[[i]])
    }
  }
  my.vector
}
#################################################################################

# Read the CO2 data into a dataframe
CO2.0 <- read.csv("CO2_Raw_Data_Processed.csv")
# Copy dataframe for later merge with enzyme data
CO2 <- CO2.0

CO2$GrowthTemp <- factor(CO2$GrowthTemp)
CO2$EvolTemp <- factor(CO2$EvolTemp)
CO2$Substrate <- revalue(CO2$Substrate, c("L"="Sucrose", "R"="Lignin"))
CO2$Generation <- revalue(CO2$Generation, c("A"="0", "B"="500", "C" = "1000", "F" = "1500"))
write.csv(CO2,"CO2_SI.csv",row.names=F)

df <- CO2[CO2$Generation %in% c("0","1500"),]
df$top.label <- paste0(df$Generation,".",df$EvolTemp)
df$top.label <- revalue(df$top.label, c("0.16"="Ancestor", "1500.16"="Evolved 16°C", "1500.22" = "Evolved 22°C"))

pdf("FigS1.pdf",height=6,width=30)
p.1 <- ggplot(df, aes(Day, CO2, color=GrowthTemp, linetype=GrowthTemp)) + ylab(expression(paste(CO[2]," (mmol ",day^-1,")"))) + labs(color="Growth T", linetype="Growth T") +
  geom_point() + geom_line() +
  theme_linedraw() + theme(strip.background = element_blank(), strip.text = element_text(color="black"), panel.grid = element_blank(), legend.key.width = unit(1,"cm")) +
  scale_color_manual(values=c("blue3","firebrick1"), labels=c("16°C","22°C")) +
  scale_linetype_manual(values=c("dashed", "solid"), labels=c("16°C","22°C"))
p.1 + facet_grid(Substrate ~ Strain + top.label, scales="free")
dev.off()

CO2.a <- CO2 %>% group_by(Strain,Generation,GrowthTemp,EvolTemp,Substrate) %>%
  summarize(max.day = Day[which.max(CO2)],cum.CO2 = trapz(Day,CO2))
CO2.a$max.day <- as.factor(CO2.a$max.day)
CO2.b <- CO2.a[CO2.a$Strain %in% c("MM02", "MM18", "MM20", "MM24", "MM31"),]
CO2.c <- summarySE(CO2.b, "cum.CO2" ,c("Generation","GrowthTemp","EvolTemp","Substrate"))
CO2.d <- summarySE(CO2.a, "cum.CO2" ,c("Generation","GrowthTemp","EvolTemp","Substrate"))
summarySE(CO2.a[CO2.a$Substrate=="Sucrose",], "cum.CO2" ,c("Generation","EvolTemp","Substrate"))
summarySE(CO2.a[CO2.a$Substrate=="Lignin",], "cum.CO2" ,c("Generation","EvolTemp","Substrate"))
CO2.b$GrowthTemp <- revalue(CO2.b$GrowthTemp, c("16"="16°C Growth T", "22" = "22°C Growth T"))

options(contrasts = c("contr.sum", "contr.poly")) # Default is: options(contrasts = c("contr.treatment", "contr.poly"))

sink("CO2stats.txt")
# Strong substrate effect; no growth temp effect and no interaction
# Log-transform results in normal residuals
m.1 <- aov(log(unlist(CO2.a[CO2.a$Generation=="0","cum.CO2"]))~unlist(CO2.a[CO2.a$Generation=="0","Substrate"])*unlist(CO2.a[CO2.a$Generation=="0","GrowthTemp"]))
rs <- resid(m.1)
shapiro.test(rs)
hist(rs)
summary(m.1)

# Use Type III results if there is an interaction and using Type II results if there is not.
# Examine generation effect, which is highly significant
# Log-transform results in near-normal residuals
m.2 <- lme(log(cum.CO2)~Substrate*GrowthTemp*Generation,random=~1|Strain,data=CO2.b,na.action=na.omit, contrasts=list(Substrate=contr.sum, GrowthTemp=contr.sum, Generation=contr.sum))
rs <- resid(m.2)
shapiro.test(rs)
hist(rs)
Anova(m.2, type = 2)
Anova(m.2, type = 3)

# The only interaction with Generation is for EvolTemp
# Residuals are perfectly normal with log transform
m.3 <- lme(log(cum.CO2)~Substrate*GrowthTemp*EvolTemp*Generation,random=~1|Strain,data=CO2.b[CO2.b$Generation!="0",],na.action=na.omit)
m.4 <- lme(log(cum.CO2)~Substrate*GrowthTemp*EvolTemp,random=~1|Strain,data=CO2.b[CO2.b$Generation!="0",],na.action=na.omit)
anova.lme(m.3) # using this as the main model for results in paper; car::Anova does not work on this model for some reason (see next line)
Anova(m.3, type = 3) # This gives a subscript error
Anova(m.4, type = 2)
Anova(m.4, type = 3)
anova.lme(m.4,type="marginal") # This gives almost identical results to car::Anova type 2

sink()

# Means paired by evolution temp over generations by substrate and growth temp
pdf("CO2_Total.pdf",height=6,width=6)
p.2 <- ggplot(CO2.c, aes(Generation, mean, color=EvolTemp)) + ylab("Total CO2 (mmol)") +
  geom_point() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1, size=0.3) +
  scale_color_manual(values=c("blue3","firebrick1"))
p.2 + facet_grid(Substrate ~ GrowthTemp, scales="free_y")
dev.off()

# Means paired by evolution temp over generations by substrate and growth temp with strains as symbols
pdf("Fig4.pdf",height=6,width=8)
p.2.1 <- ggplot(CO2.b, aes(Generation, cum.CO2, color=EvolTemp)) + ylab(expression(paste("Cumulative respiration (mmol ",CO[2],")"))) + labs(color="Evolution T") +
  geom_point(aes(shape=Strain)) +
  stat_summary(fun.data = "mean_se",geom="crossbar",width=0.3,fatten=1) +
  theme_linedraw() + theme(strip.background = element_blank(), strip.text = element_text(color="black"),
                           panel.grid = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values=c("blue3","firebrick1"), labels=c("16°C","22°C"))
p.2.1 + facet_grid(Substrate ~ GrowthTemp, scales="free_y")
dev.off()

# Means paired by growth temp on both substrates for all ancestral strains
pdf("CO2_Total_A.pdf",height=6,width=6)
p.2.2 <- ggplot(CO2.d[CO2.d$Generation=="0",], aes(GrowthTemp, mean, color=GrowthTemp)) + ylab("Total CO2 (mmol)") +
  geom_point() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1, size=0.3) +
  scale_color_manual(values=c("blue3","firebrick1"))
p.2.2 + facet_grid(Substrate ~ ., scales="free_y")
dev.off()

df <- CO2.a[CO2.a$Generation %in% c("0","1500"),]
df$top.label <- paste0(df$Generation,".",df$EvolTemp)
df$top.label <- revalue(df$top.label, c("0.16"="Ancestor", "1500.16"="Evolved 16°C", "1500.22" = "Evolved 22°C"))

# Histogram of max CO2 production by day over generations, growth temps, and evolution temps
pdf("Fig2.pdf",height=5,width=6)
p3 <- ggplot(df, aes(max.day, fill=GrowthTemp, color=NA)) + 
  geom_bar() + xlab("Day of maximum respiration") + ylab("# of strains") + labs(fill="Growth T") +
  theme_linedraw() + theme(strip.background = element_blank(), strip.text = element_text(color="black")) +
  scale_fill_manual(values=c("blue3","firebrick1"), labels=c("16°C","22°C")) + 
  scale_color_manual(values=c("blue3","firebrick1"))
p3 + facet_grid(Substrate ~ top.label, scales="free_y")
dev.off()

# For enzyme analyses, if you want to quality filter the data, leave this as TRUE
rm.outliers <- TRUE

# Read in enzyme activity data
Enz3 <- read.csv("Enz3.csv")

enzymes <- c("AG","AP","BG","BX","CBH","LAP","NAG","PPO","OX","PER")
Enz3$Temp <- factor(Enz3$Temp)
Enz3$Sample <- gsub("W1$","W01",Enz3$Sample)
Enz3$Sample <- gsub("W7$","W07",Enz3$Sample)
Enz3$Sample <- gsub("C1$","C01",Enz3$Sample)
Enz3$Sample <- gsub("C7$","C07",Enz3$Sample)
Enz3$Sample <- gsub("MM2L","MM02L",Enz3$Sample)
Enz3$Sample <- gsub("MM2R","MM02R",Enz3$Sample)
Enz3$Sample <- gsub("MM6L","MM06L",Enz3$Sample)
Enz3$Sample <- gsub("MM6R","MM06R",Enz3$Sample)
Enz3$Sample <- gsub("MM2A","MM02A",Enz3$Sample)
Enz3$Strain <- gsub("MM2$","MM02",Enz3$Strain)
Enz3$Strain <- gsub("MM6$","MM06",Enz3$Strain)

if(rm.outliers==TRUE){
  remove <- data.frame(rbind(
    list("MM02RC14",34,19.5,"PER"),
    list("MM10LW14",22,7.8,"AG"),
    list("MM10RW14",22,7.8,"AG"),
    list("MM20LC14",22,15.6,"AG"),
    list("MM23LW28",34,500,"BG"),
    list("MM24RW28",4,500,"LAP"),
    list("MM30LW28",34,62.5,"NAG"),
    list("MM31LW07",4,62.5,"BX"),
    list("MM26LC07",22,500,"AG")
  ))
  attach(Enz3)
  for(i in 1:dim(remove)[1]){
    Enz3[Sample==remove[[i,1]] & Temp==remove[[i,2]] & Concen==remove[[i,3]] & Variate==remove[[i,4]],"Activity"]<-NA
  }
  detach(Enz3)
  
  # Convert negatives and very low numbers to small number
  Enz3$Activity[Enz3$Activity<0.0001] <- 0.00001
}

Enz3$Generation <- substr(Enz3$Sample,8,8)
Enz3$Generation[Enz3$Generation!="F"] <- "A"
Enz3$Day <- gsub(28,21,Enz3$Day)
Enz3$Day <- as.integer(Enz3$Day)
Enz3$EvolTemp[is.na(Enz3$EvolTemp)] <- 16

Enz3 <- merge(Enz3,CO2.0)

Enz3$Activity.C <- Enz3$Activity/Enz3$CO2

if(rm.outliers==TRUE){
  Enz3$Activity.C[Enz3$Activity.C<0.0001] <- 0.00001
}

# Switch activity variables normalized to CO2 or not
Enz3$dep <- Enz3$Activity

# Remove missing values that cause plotting errors
Enz3 <- Enz3[!is.na(Enz3$dep),]

attach(Enz3)
# Sort the dataset
Enz3 <- Enz3[order(Variate,Strain,Generation,Substrate,EvolTemp,GrowthTemp,Day,Temp,Concen),]
detach(Enz3)

# Fit the M-M equation
# Use try() to skip over fits that fail
# This step takes a while
# Skip to read.csv to get params.frame below to save time
attach(Enz3)
pdf("fits.pdf")
x.len <- seq(len=dim(Enz3)[1])
param <- tapply(x.len, list(Temp,Day,GrowthTemp,EvolTemp,Substrate,Generation,Strain,Variate), 
                function(i){
                  nls.m <- try(nls(y~Vmax*x/(Km+x),data=list(y=dep[i],x=Concen[i]),start=list(Vmax=max(dep[i],na.rm=T),Km=max(Concen[i]/4,na.rm=T))))
                  Vmax <- try(coef(nls.m)[[1]],silent=T)
                  Km <- try(coef(nls.m)[[2]],silent=T)
                  lowerV <- try(confint2(nls.m)[[1]],silent=T)
                  upperV <- try(confint2(nls.m)[[3]],silent=T)
                  lowerK <- try(confint2(nls.m)[[2]],silent=T)
                  upperK <- try(confint2(nls.m)[[4]],silent=T)
                  fit <- try(100*(upperV - lowerV)/Vmax,silent=T)
                  model.pts <- try(Vmax*Concen[i]/(Km+Concen[i]),silent=T)
                  plot(dep[i]~Concen[i],main=paste("Enz:",Variate[i][1]," Str:",Strain[i][1]," Gen:",Generation[i][1]," Sub:",Substrate[i][1],
                                                   " EvolT:",EvolTemp[i][1]," GrowthT:",GrowthTemp[i][1]," Day:",Day[i][1]," Temp:",Temp[i][1],sep=""))
                  try(lines(model.pts~Concen[i],col="firebrick1"),silent=T)
                  text(Concen[i],dep[i],Concen[i],pos=4)
                  mtext(paste("Vmax=",ifelse(is.numeric(Vmax),round(Vmax,3),NA)," Km=",ifelse(is.numeric(Km),round(Km,3),NA),
                              " Fit=",ifelse(is.numeric(fit),round(fit,3),NA),sep=""))
                  c(Vmax,Km,lowerV,upperV,lowerK,upperK)
                })
dev.off()

params.frame <- as.data.frame.table(tapply(x.len, list(Temp,Day,GrowthTemp,EvolTemp,Substrate,Generation,Strain,Variate),function(i)1))
names(params.frame) <- c("Temp","Day","GrowthTemp","EvolTemp","Substrate","Generation","Strain","Variate","Vmax")
detach(Enz3)

# Extract the Vmax and Km values
params.frame$Vmax <- as.numeric(list.as.vectorNA(sapply(sapply(param,'[',1),'[[',1)))
params.frame$Km <- as.numeric(list.as.vectorNA(sapply(sapply(param,'[',2),'[[',1)))
params.frame$lowerV <- as.numeric(list.as.vectorNA(sapply(sapply(param,'[',3),'[[',1)))
params.frame$upperV <- as.numeric(list.as.vectorNA(sapply(sapply(param,'[',4),'[[',1)))
params.frame$lowerK <- as.numeric(list.as.vectorNA(sapply(sapply(param,'[',5),'[[',1)))
params.frame$upperK <- as.numeric(list.as.vectorNA(sapply(sapply(param,'[',6),'[[',1)))
params.frame$fitV <- (params.frame$upperV - params.frame$lowerV)/params.frame$Vmax
params.frame$fitK <- (params.frame$upperK - params.frame$lowerK)/params.frame$Km

# Convert temperature to a numeric value
params.frame$Temp <- as.numeric(as.vector(params.frame$Temp))

# Remove missing observations
params.frame <- subset(params.frame,!is.na(Vmax))

# Insert NA for incorrect model fits (negative Km) or poor fits
params.frame[(!is.na(params.frame$Km) & params.frame$Km<0),5:dim(params.frame)[2]] <- NA
params.frame[(!is.na(params.frame$fitV) & params.frame$fitV>2),5:dim(params.frame)[2]] <- NA
params.frame[(!is.na(params.frame$fitK) & params.frame$fitK>4),"Km"] <- NA

# Remove missing observations again
params.frame <- params.frame[!is.na(params.frame$Vmax),]

params.frame$logVmax <- log(params.frame$Vmax)
params.frame$logKm <- log(params.frame$Km)

params.frame$Substrate <- revalue(params.frame$Substrate, c("L"="Sucrose", "R"="Lignin"))

write.csv(params.frame,"params.frame.csv",row.names=F)

# Skip to here if you just want the data frame with quality filtered, fitted parameters
params.frame <- read.csv("params.frame.csv")

if(rm.outliers==TRUE){
  # Remove questionable data based on plots of Vmax versus assay temperature
  remove <- data.frame(rbind(
    list("A","Lignin",21,"W1303","LAP",16,16,16),
    list("A","Lignin",21,"W1289","NAG",16,16,34),
    list("A","Lignin",21,"W1289","AG",16,16,34),
    list("F","Sucrose",14,"MM02","AG",16,16,34),
    list("F","Sucrose",14,"MM18","AG",16,16,4:34),
    list("F","Sucrose",14,"MM24","AG",16,16,34),
    list("F","Sucrose",21,"MM18","AG",16,16,4:34),
    list("F","Sucrose",21,"MM20","AG",16,16,34),
    list("F","Sucrose",21,"MM24","AG",16,16,4:34),
    list("F","Sucrose",21,"MM20","AG",22,16,4:34),
    list("F","Sucrose",14,"MM18","AP",16,16,28),
    list("F","Lignin",14,"MM18","BG",16,16,16:28),
    list("F","Lignin",14,"MM31","BG",16,22,4:28),
    list("F","Lignin",21,"MM02","BG",16,22,4:28),
    list("F","Lignin",21,"MM18","BG",16,22,28),
    list("F","Lignin",21,"MM24","BG",16,22,4:28),
    list("F","Sucrose",14,"MM18","BG",16,16,4:34),
    list("F","Sucrose",14,"MM24","BG",16,16,16),
    list("F","Sucrose",21,"MM02","BG",16,22,28),
    list("F","Sucrose",21,"MM18","BG",16,22,28),
    list("F","Sucrose",21,"MM20","BG",16,16,16),
    list("F","Sucrose",21,"MM24","BG",16,16,16),
    list("F","Sucrose",21,"MM24","BG",16,22,28),
    list("F","Lignin",14,"MM02","BG",22,22,22:28),
    list("F","Lignin",14,"MM24","BG",22,22,22:28),
    list("F","Lignin",21,"MM20","BG",22,16,22:28),
    list("F","Sucrose",14,"MM02","BG",22,16,16),
    list("F","Sucrose",14,"MM02","BG",22,22,16:34),
    list("F","Sucrose",14,"MM18","BG",22,16,4:34),
    list("F","Sucrose",14,"MM20","BG",22,16,4:34),
    list("F","Sucrose",14,"MM24","BG",22,22,4:34),
    list("F","Sucrose",14,"MM31","BG",22,16,16),
    list("F","Sucrose",21,"MM20","BG",22,16,4:34),
    list("F","Lignin",14,"MM18","BX",16,16,28),
    list("F","Lignin",14,"MM31","BX",16,22,4:28),
    list("F","Lignin",21,"MM18","BX",16,22,28),
    list("F","Lignin",21,"MM24","BX",16,22,10:28),
    list("F","Sucrose",14,"MM18","BX",16,16,4:34),
    list("F","Sucrose",14,"MM24","BX",16,16,4:34),
    list("F","Sucrose",21,"MM20","BX",16,16,4:34),
    list("F","Sucrose",21,"MM24","BX",16,16,4:34),
    list("F","Lignin",14,"MM02","BX",22,22,22:28),
    list("F","Lignin",14,"MM24","BX",22,22,22:28),
    list("F","Lignin",21,"MM02","BX",22,22,22),
    list("F","Lignin",21,"MM24","BX",22,22,22:28),
    list("F","Lignin",21,"MM20","BX",22,16,22:28),
    list("F","Sucrose",14,"MM02","BX",22,22,22:28),
    list("F","Sucrose",14,"MM24","BX",22,22,22:28),
    list("F","Sucrose",21,"MM20","BX",22,16,4:34),
    list("F","Lignin",14,"MM18","CBH",16,16,16:28),
    list("F","Lignin",14,"MM31","CBH",16,22,16:28),
    list("F","Lignin",21,"MM02","CBH",16,22,4:28),
    list("F","Lignin",21,"MM18","CBH",16,22,28),
    list("F","Lignin",21,"MM24","CBH",16,22,4:28),
    list("F","Sucrose",14,"MM18","CBH",16,16,4:34),
    list("F","Sucrose",14,"MM24","CBH",16,16,4:34),
    list("F","Sucrose",21,"MM20","CBH",16,16,16),
    list("F","Sucrose",21,"MM24","CBH",16,16,16),
    list("F","Lignin",14,"MM02","CBH",22,22,22),
    list("F","Lignin",14,"MM24","CBH",22,22,22),
    list("F","Lignin",21,"MM20","CBH",22,16,22:28),
    list("F","Sucrose",14,"MM02","CBH",22,16,16:34),
    list("F","Sucrose",14,"MM02","CBH",22,22,16:28),
    list("F","Sucrose",14,"MM18","CBH",22,16,4:34),
    list("F","Sucrose",14,"MM20","CBH",22,16,4:34),
    list("F","Sucrose",14,"MM24","CBH",22,22,4:28),
    list("F","Sucrose",14,"MM31","CBH",22,16,16),
    list("F","Sucrose",21,"MM20","CBH",22,16,4:28),
    list("F","Sucrose",14,"MM18","LAP",16,16,4:34),
    list("F","Sucrose",14,"MM31","LAP",16,22,22:28),
    list("F","Sucrose",14,"MM02","LAP",22,22,22:28),
    list("F","Sucrose",14,"MM24","LAP",22,22,22)
  ))
  attach(params.frame)
  for(i in 1:dim(remove)[1]){
    params.frame[Generation==remove[[i,1]] & Substrate==remove[[i,2]] & Day==remove[[i,3]] & Strain==remove[[i,4]] & Variate==remove[[i,5]] & 
                   GrowthTemp==remove[[i,6]] & EvolTemp==remove[[i,7]] & Temp %in% remove[[i,8]],c("Vmax","logVmax","Km","logKm")] <- NA
  }
  detach(params.frame)
}

params.frame$EvolTemp <- factor(params.frame$EvolTemp)

# Plotting Vmax and Km for various combinations of the independent variables
pdf("Vmax1.pdf",height=20,width=50)
d <- ggplot(params.frame[params.frame$Day %in% c(14,21),], aes(Temp, Vmax, color=EvolTemp)) + geom_point() + geom_line()
d + facet_grid(Variate + GrowthTemp ~ Generation + Substrate + Day + Strain, scales="free") + scale_color_manual(values=c("blue3","firebrick1")) + theme_gray()
dev.off()

pdf("VmaxAF.pdf",height=20,width=25)
d <- ggplot(params.frame[params.frame$Day %in% c(14,21) & params.frame$Strain %in% c("MM02","MM18","MM20","MM24","MM31"),], aes(Temp, Vmax, color=EvolTemp)) + geom_point() +geom_line()
d + facet_grid(Variate + GrowthTemp ~ Substrate + Day + Strain + Generation, scales="free") + scale_color_manual(values=c("blue3","firebrick1")) + theme_gray()
dev.off()

pdf("VmaxF.pdf",height=20,width=15)
d <- ggplot(params.frame[params.frame$Day %in% c(14,21) & params.frame$Generation=="F",], aes(Temp, Vmax, color=EvolTemp)) + geom_point() +geom_line()
d + facet_grid(Variate + GrowthTemp ~ Generation + Substrate + Day + Strain, scales="free") + scale_color_manual(values=c("blue3","firebrick1")) + theme_gray()
dev.off()

pdf("VmaxF22-14.pdf",height=6,width=10)
d <- ggplot(params.frame[params.frame$Variate %in% c("BG","BX","CBH") & params.frame$GrowthTemp==22 & params.frame$Day==14 & params.frame$Generation=="F",], 
            aes(Temp, Vmax, color=EvolTemp)) + geom_point() + geom_line()
d + facet_grid(Variate ~ Substrate + Strain, scales="free_y") + scale_color_manual(values=c("blue3","firebrick1")) + theme_gray()
dev.off()

pdf("VmaxF16-14.pdf",height=6,width=10)
d <- ggplot(params.frame[params.frame$Variate %in% c("BG","BX","CBH") & params.frame$GrowthTemp==16 & params.frame$Day==14 & params.frame$Generation=="F",], 
            aes(Temp, Vmax, color=EvolTemp)) + geom_point() +geom_line()
d + facet_grid(Variate ~ Substrate + Strain, scales="free_y") + scale_color_manual(values=c("blue3","firebrick1")) + theme_gray()
dev.off()

pdf("VmaxF22-21.pdf",height=6,width=10)
d <- ggplot(params.frame[params.frame$Variate %in% c("BG","BX","CBH") & params.frame$GrowthTemp==22 & params.frame$Day==21 & params.frame$Generation=="F",], 
            aes(Temp, Vmax, color=EvolTemp)) + geom_point() +geom_line()
d + facet_grid(Variate ~ Substrate + Strain, scales="free_y") + scale_color_manual(values=c("blue3","firebrick1")) + theme_gray()
dev.off()

pdf("Km1.pdf",height=20,width=50)
d <- ggplot(params.frame[params.frame$Day %in% c(14,21),], aes(Temp, Km, color=EvolTemp)) + geom_point() +geom_line()
d + facet_grid(Variate + GrowthTemp ~ Generation + Substrate + Day + Strain, scales="free") + scale_color_manual(values=c("blue3","firebrick1")) + theme_gray()
dev.off()

pdf("KmF.pdf",height=20,width=15)
d <- ggplot(params.frame[params.frame$Day %in% c(14,21) & params.frame$Generation=="F",], aes(Temp, Km, color=EvolTemp)) + geom_point() +geom_line()
d + facet_grid(Variate + GrowthTemp ~ Generation + Substrate + Day + Strain, scales="free") + scale_color_manual(values=c("blue3","firebrick1")) + theme_gray()
dev.off()

# Regressions of log(parameter) versus assay temperature
attach(params.frame)

Vmax.regs <- by(params.frame, list(Variate,Strain,Generation,Substrate,Day,GrowthTemp,EvolTemp),
                function(i){
                  sumlm <- try(summary(lm(logVmax~Temp, data=i)))
                  if(is.character(sumlm)){
                    c(Vmax.Int=NA, Vmax.Slope=NA, Vmax.SE=NA, Vmax.Rsq=NA)}
                  else{
                    c(Vmax.Int=sumlm[[4]][1], Vmax.Slope=sumlm[[4]][2], Vmax.SE=sumlm[[4]][4], Vmax.Rsq=sumlm[[8]])}
                }
)
Vmax.stats <- do.call("rbind",Vmax.regs)

Km.regs <- by(params.frame, list(Variate,Strain,Generation,Substrate,Day,GrowthTemp,EvolTemp),
              function(i){
                sumlm <- try(summary(lm(logKm~Temp, data=i)))
                if(is.character(sumlm)){
                  c(Km.Int=NA, Km.Slope=NA, Km.SE=NA, Km.Rsq=NA)}
                else{
                  c(Km.Int=sumlm[[4]][1], Km.Slope=sumlm[[4]][2], Km.SE=sumlm[[4]][4], Km.Rsq=sumlm[[8]])}
              }
)
Km.stats <- do.call("rbind",Km.regs)

# Fitting the macromolecular rate theory; not ideal because of temperature range, though
# Not used further
Vmax.MMRT <- by(params.frame, list(Variate,Strain,Generation,Substrate,Day,GrowthTemp,EvolTemp),
                function(i){
                  T0 = 289 #K
                  kB = 1.3806e-26 #kJ/K
                  h = 6.6261e-37 #kJ s
                  R = 0.008314 #kJ/(mol K)
                  nls.m <- try(nls(logVmax~log(kB*(Temp+273)/h) - (H+C*((Temp+273)-T0))/(R*(Temp+273)) + (S+C*(log((Temp+273))-log(T0)))/R,data=i,start=list(H=0,S=0,C=-1)))
                  H <- try(coef(nls.m)[[1]],silent=T)
                  S <- try(coef(nls.m)[[2]],silent=T)
                  C <- try(coef(nls.m)[[3]],silent=T)
                  if(is.character(nls.m)){
                    c(delH=NA,delS=NA,C=NA)}
                  else{
                    c(delH=H,delS=S,C=C)}
                }
)
MMRT.stats <- do.call("rbind",Vmax.MMRT)

slopes.frame <- as.data.frame.table(tapply(seq(len=dim(params.frame)[1]), list(Variate,Strain,Generation,Substrate,Day,GrowthTemp,EvolTemp),
                                           function(i)1))
names(slopes.frame) <- c("Variate","Strain","Generation","Substrate","Day","GrowthTemp","EvolTemp","dummy")
slopes.frame <- subset(slopes.frame,!is.na(slopes.frame$dummy))

detach(params.frame)

# Combine stats across all samples
stats.frame <- cbind(slopes.frame, Vmax.stats, Km.stats, MMRT.stats)
stats.frame$dummy <- NULL

# Remove linear fits to only 1 or 2 points
stats.frame[is.na(stats.frame$Vmax.SE),c("Vmax.Int","Vmax.Slope")] <- NA
stats.frame[is.na(stats.frame$Km.SE),c("Km.Int","Km.Slope")] <- NA

# Calculating Q10 values
stats.frame$Vmax.Q10 <- exp(10*stats.frame$Vmax.Slope)
stats.frame$Km.Q10 <- exp(10*stats.frame$Km.Slope)
stats.frame$Q10.MMRT <- exp(10*(stats.frame$delH-5*stats.frame$C)/(0.008314*289^2))

# Parameters normalized to 16 degrees based on the regression
stats.frame$Vmax.16 <- stats.frame$Vmax.Int + 16*stats.frame$Vmax.Slope
stats.frame$Km.16 <- stats.frame$Km.Int + 16*stats.frame$Km.Slope

if(rm.outliers==TRUE){
  # Insert minimum Vmax for each enzyme as the detection limit
  stats.frame <- stats.frame %>% group_by(Variate) %>% mutate(Vmax.Min=min(Vmax.16,na.rm=T))
  stats.frame$Vmax.16[is.na(stats.frame$Vmax.16)] <- stats.frame$Vmax.Min[is.na(stats.frame$Vmax.16)]
}

write.csv(stats.frame,"VmaxKmTS.csv",quote=F,row.names=FALSE)

stats.frame <- read.csv("VmaxKmTS.csv",header=T)
stats.frame$GrowthTemp <- factor(stats.frame$GrowthTemp)
stats.frame$EvolTemp <- factor(stats.frame$EvolTemp)
stats.frame$Day <- factor(stats.frame$Day)

i <- stats.frame$Generation=="A" & stats.frame$Day==21
SF.1 <- stats.frame[i,]
SF.1$GrowthTemp <- factor(SF.1$GrowthTemp)

# Use Type III results if there is an interaction and use Type II results if there is not.
# Some models fail due to missing data; using try() to skip them
enzymes <- c("AG","AP","BG","BX","CBH","LAP","NAG","OX","PER","PPO")
sink("GrowthTempSubstrateVmax.txt")
for (enz in enzymes) {
  m.5 <- try(lme(min_rank(Vmax.16)~Substrate*GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz,],na.action=na.omit))
  rs <- try(resid(m.5))
  print(enz)
  print(try(shapiro.test(rs)))
  try(hist(rs))
  print(try(Anova(m.5,type = 2)))
  print(try(anova.lme(m.5,type="marginal")))
}
for (enz in enzymes) {
  print(enz)
  m.5.1 <- try(lme(min_rank(Vmax.16)~Substrate,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$GrowthTemp=="16",],na.action=na.omit))
  print(try(anova.lme(m.5.1)))
  m.5.2 <- try(lme(min_rank(Vmax.16)~Substrate,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$GrowthTemp=="22",],na.action=na.omit))
  print(try(anova.lme(m.5.2)))
}
for (enz in enzymes) {
  print(enz)
  m.5.1 <- try(lme(min_rank(Vmax.16)~Substrate,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$GrowthTemp=="16",],na.action=na.omit))
  print(try(anova.lme(m.5.1)))
  m.5.2 <- try(lme(min_rank(Vmax.16)~Substrate,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$GrowthTemp=="22",],na.action=na.omit))
  print(try(anova.lme(m.5.2)))
}
for (enz in enzymes) {
  print(enz)
  m.5.3 <- try(lme(min_rank(Vmax.16)~GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$Substrate=="Lignin",],na.action=na.omit))
  print(try(anova.lme(m.5.3)))
  m.5.4 <- try(lme(min_rank(Vmax.16)~GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$Substrate=="Sucrose",],na.action=na.omit))
  print(try(anova.lme(m.5.4)))
}
sink()

sink("GrowthTempSubstrateVmaxTS.txt")
for (enz in enzymes) {
  m.5 <- try(lme(min_rank(Vmax.Slope)~Substrate*GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz,],na.action=na.omit))
  rs <- try(resid(m.5))
  print(enz)
  print(try(shapiro.test(rs)))
  try(hist(rs))
  print(try(Anova(m.5,type = 2)))
  print(try(anova.lme(m.5,type="marginal")))
}
for (enz in enzymes) {
  print(enz)
  m.5.3 <- try(lme(min_rank(Vmax.Slope)~GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$Substrate=="Lignin",],na.action=na.omit))
  print(try(anova.lme(m.5.3)))
  m.5.4 <- try(lme(min_rank(Vmax.Slope)~GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$Substrate=="Sucrose",],na.action=na.omit))
  print(try(anova.lme(m.5.4)))
}
sink()

sink("GrowthTempSubstrateKm.txt")
for (enz in enzymes) {
  m.5 <- try(lme(min_rank(Km.16)~Substrate*GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz,],na.action=na.omit))
  rs <- try(resid(m.5))
  print(enz)
  print(try(shapiro.test(rs)))
  try(hist(rs))
  print(try(Anova(m.5,type = 2)))
  print(try(Anova(m.5,type = 3)))
  print(try(anova.lme(m.5,type="marginal")))
}
for (enz in enzymes) {
  print(enz)
  m.5.3 <- try(lme(min_rank(Km.16)~GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$Substrate=="Lignin",],na.action=na.omit))
  print(try(anova.lme(m.5.3)))
  m.5.4 <- try(lme(min_rank(Km.16)~GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$Substrate=="Sucrose",],na.action=na.omit))
  print(try(anova.lme(m.5.4)))
}
sink()

sink("GrowthTempSubstrateKmTS.txt")
for (enz in enzymes) {
  m.5 <- try(lme(min_rank(Km.Slope)~Substrate*GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz,],na.action=na.omit))
  rs <- try(resid(m.5))
  print(enz)
  print(try(shapiro.test(rs)))
  try(hist(rs))
  print(try(Anova(m.5,type = 2)))
  print(try(anova.lme(m.5,type="marginal")))
}
for (enz in enzymes) {
  print(enz)
  m.5.3 <- try(lme(min_rank(Km.Slope)~GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$Substrate=="Lignin",],na.action=na.omit))
  print(try(anova.lme(m.5.3)))
  m.5.4 <- try(lme(min_rank(Km.Slope)~GrowthTemp,random=~1|Strain,data=SF.1[SF.1$Variate==enz & SF.1$Substrate=="Sucrose",],na.action=na.omit))
  print(try(anova.lme(m.5.4)))
}
sink()

pdf("Fig3.pdf",width = 7,height = 10)
# Boxplots paired by growth temp
sigs <- data.frame(Variate=enzymes,Substrate=c(rep("Lignin",10),rep("Sucrose",10)),GrowthTemp="16",Vmax.16=4.5,sig=NA)
if(rm.outliers==TRUE){
  sigs$sig[c(1,4,5,16)] <- "*"
}
if(rm.outliers==FALSE){
  sigs$sig[c(1,5,16)] <- "*"
}
p.6 <- ggplot(SF.1, aes_string("Variate", "Vmax.16", color="GrowthTemp", shape="GrowthTemp")) + xlab(NULL) + ylab("log(Vmax)") + labs(title = "a)", color="Growth T", shape="Growth T") + 
  geom_boxplot(position=position_dodge(width=0.8),width=0.7,outlier.size=0) + geom_point(position=position_jitterdodge(jitter.width=0, dodge.width=0.8)) +
  geom_vline(xintercept=seq(1.5,9.5,1),color="black",size=0.1) +
  geom_text(data=sigs,aes(label=sig),color="black",size=10) +
  theme_linedraw() + theme(strip.background = element_blank(), strip.text = element_text(color="black"),
                           panel.grid = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values=c("blue3","firebrick1"), labels=c("16°C","22°C")) + 
  scale_shape_manual(values=c(16,17), labels=c("16°C","22°C")) +
  facet_grid(Substrate~., scales="fixed")

sigs <- data.frame(Variate=enzymes,Substrate=c(rep("Lignin",10),rep("Sucrose",10)),GrowthTemp="16",Vmax.Slope=0.25,sig=NA)
if(rm.outliers==FALSE){
  sigs$sig[c(15)] <- "*"
}
# No significant growth T effects
p.6.1 <- ggplot(SF.1, aes_string("Variate", "Vmax.Slope", color="GrowthTemp", shape="GrowthTemp")) + xlab(NULL) + ylab("Vmax temperature sensitivity") + labs(title = "b)", color="Growth T", shape="Growth T") + 
  geom_boxplot(position=position_dodge(width=0.8),width=0.7,outlier.size=0) + geom_point(position=position_jitterdodge(jitter.width=0, dodge.width=0.8)) +
  geom_vline(xintercept=seq(1.5,9.5,1),color="black",size=0.1) +
  geom_text(data=sigs,aes(label=sig),color="black",size=10) +
  theme_linedraw() + theme(strip.background = element_blank(), strip.text = element_text(color="black"),
                           panel.grid = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values=c("blue3","firebrick1"), labels=c("16°C","22°C")) +
  scale_shape_manual(values=c(16,17), labels=c("16°C","22°C")) +
  facet_grid(Substrate~., scales="fixed")

plot_grid(p.6, p.6.1, ncol=1, align="v")

dev.off()

pdf("FigS2.pdf",width = 7,height = 10)
# Boxplots paired by growth temp
sigs <- data.frame(Variate=enzymes,Substrate=c(rep("Lignin",10),rep("Sucrose",10)),GrowthTemp="16",Km.16=8,sig=NA)
sigs$sig[13] <- "*"
p.6.2 <- ggplot(SF.1, aes_string("Variate", "Km.16", color="GrowthTemp", shape="GrowthTemp")) + xlab(NULL) + ylab("log(Km)") + labs(title = "a)", color="Growth T", shape="Growth T") + 
  geom_boxplot(position=position_dodge(width=0.8),width=0.7,outlier.size=0) + geom_point(position=position_jitterdodge(jitter.width=0, dodge.width=0.8)) +
  geom_vline(xintercept=seq(1.5,9.5,1),color="black",size=0.1) +
  geom_text(data=sigs,aes(label=sig),color="black",size=10) +
  theme_linedraw() + 
  theme(strip.background = element_blank(), strip.text = element_text(color="black"), panel.grid = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values=c("blue3","firebrick1"), labels=c("16°C","22°C")) + 
  scale_shape_manual(values=c(16,17), labels=c("16°C","22°C")) +
  facet_grid(Substrate~., scales="fixed")

# No significant growth T effects on Km TS
p.6.3 <- ggplot(SF.1, aes_string("Variate", "Km.Slope", color="GrowthTemp", shape="GrowthTemp")) + xlab(NULL) + ylab("Km temperature sensitivity") + labs(title = "b)", color="Growth T", shape="Growth T") + 
  geom_boxplot(position=position_dodge(width=0.8),width=0.7,outlier.size=0) + geom_point(position=position_jitterdodge(jitter.width=0, dodge.width=0.8)) +
  geom_vline(xintercept=seq(1.5,9.5,1),color="black",size=0.1) +
  theme_linedraw() + 
  theme(strip.background = element_blank(), strip.text = element_text(color="black"), panel.grid = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values=c("blue3","firebrick1"), labels=c("16°C","22°C")) +
  scale_shape_manual(values=c(16,17), labels=c("16°C","22°C")) +
  facet_grid(Substrate~., scales="fixed")

plot_grid(p.6.2, p.6.3, ncol=1, align="v")

dev.off()

SF.1$GrowthTemp <- revalue(SF.1$GrowthTemp, c("16"="16°C Growth T", "22" = "22°C Growth T"))
sigs <- data.frame(Variate=enzymes,GrowthTemp=c(rep("16°C Growth T",10),rep("22°C Growth T",10)),Substrate="Lignin",Vmax.16=5.2,sig="*")
sigs$sig[c(4,8,14,18,5,10,20)] <- NA
# Boxplot paired by substrate
pdf("Fig1.pdf",width = 7,height = 7)
p.7 <- ggplot(SF.1, aes_string("Variate", "Vmax.16", color="Substrate", shape="Substrate")) + xlab(NULL) + ylab("log(Vmax)") +
  geom_boxplot(position=position_dodge(width=0.8),width=0.7,outlier.size=0) + geom_point(position=position_jitterdodge(jitter.width=0, dodge.width=0.8)) +
  geom_vline(xintercept=seq(1.5,9.5,1),color="black",size=0.1) +
  geom_text(data=sigs,aes(label=sig),color="black",size=10) +
  theme_linedraw() + 
  theme(strip.background = element_blank(), strip.text = element_text(color="black"), panel.grid = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values=c("black","orange")) +
  scale_shape_manual(values=c(16,17))
print(p.7 + facet_grid(GrowthTemp~., scales="fixed"))
dev.off()

# Other variables by substrate
y.labs <- c("log(Vmax)", "Vmax temperature sensitivity", "log(Km)", "Km temperature sensitivity"); titles <- c("a)","b)","c)","d)"); k <- 1
for (j in c("Vmax.Slope", "Km.16", "Km.Slope")) {
  pdf(paste0("SubBoxplot.",j,".pdf"),width = 7,height = 7)
  p.7 <- ggplot(SF.1, aes_string("Variate", j, color="Substrate", shape="Substrate")) + xlab(NULL) + ylab(y.labs[k]) +
    geom_boxplot(position=position_dodge(width=0.8),width=0.7,outlier.size=0) + geom_point(position=position_jitterdodge(jitter.width=0, dodge.width=0.8)) +
    geom_vline(xintercept=seq(1.5,9.5,1),color="black",size=0.1) +
    theme_linedraw() + 
    theme(strip.background = element_blank(), strip.text = element_text(color="black"), panel.grid = element_blank(), text = element_text(size = 16)) +
    scale_color_manual(values=c("black","orange")) +
    scale_shape_manual(values=c(16,17))
    print(p.7 + facet_grid(GrowthTemp~., scales="fixed"))
  dev.off()
  k <- k+1
}

# All data from A and F generations for evolved strains
ALL.EVO <- stats.frame[stats.frame$Strain %in% c("MM02","MM18","MM20","MM24","MM31") & stats.frame$EvolTemp %in% c(16,22) & stats.frame$Day %in% c(14,21),]
# All data from A and F generations at 16C
ALL.16 <- stats.frame[stats.frame$Strain %in% c("MM02","MM18","MM20","MM24","MM31") & stats.frame$EvolTemp==16 & stats.frame$Day %in% c(14,21),]
# All data from generation
F.ALL <- stats.frame[stats.frame$Generation=="F",]
# F.ALL$Day <- revalue(F.ALL$Day, c("14"="14 Days", "21"="21 Days"))
F.ALL$GrowthTemp <- revalue(F.ALL$GrowthTemp, c("16"="16°C Growth T", "22" = "22°C Growth T"))

# NAG Vmax temperature sensitivity
df.1 <- summarySE(F.ALL[F.ALL$Variate=="NAG",], "Vmax.Slope" ,c("Substrate","GrowthTemp","Day","EvolTemp"))
pdf("Fig6.pdf", height = 5, width = 6)
p.8 <- ggplot(df.1, aes(Day, mean, color=EvolTemp, shape=EvolTemp)) + ylab("Vmax temperature sensitivity") + xlab("Day") + labs(color="Evolution T", shape="Evolution T") +
  geom_errorbar(data=df.1[df.1$EvolTemp=="22",], aes(ymin=mean-se, ymax=mean+se), width=0.1, size=0.3) +
  geom_point(data=df.1[df.1$EvolTemp=="22",],size=3) +
  geom_errorbar(data=df.1[df.1$EvolTemp=="16",], aes(ymin=mean-se, ymax=mean+se), width=0.1, size=0.3) +
  geom_point(data=df.1[df.1$EvolTemp=="16",],size=3) +
  theme_linedraw() + theme(strip.background = element_blank(), strip.text = element_text(color="black"),
                           panel.grid = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values=c("blue3","firebrick1"), labels=c("16°C","22°C")) +
  scale_shape_manual(values=c(16,17), labels=c("16°C","22°C"))
  p.8 + facet_grid(Substrate ~ GrowthTemp, scales="fixed")
dev.off()

# All enzymes; Vmax.16 by generation for 16C strains
df.1 <- summarySE(ALL.16, "Vmax.16" ,c("Variate","Substrate","GrowthTemp","Day","Generation"))
df.1$Day <- revalue(df.1$Day, c("14"="14 Days", "21" = "21 Days"))
df.1$Generation <- revalue(df.1$Generation, c("A"="0", "F" = "1500"))
pdf("Vmax.Gen.pdf", height = 10, width = 8)
p.8.1 <- ggplot(df.1, aes(Generation, mean, color=GrowthTemp, shape=GrowthTemp)) + ylab("log(Vmax)") + labs(color="Growth T", shape="Growth T") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1, size=0.3, data=df.1[df.1$GrowthTemp==22,]) +
  geom_point(data=df.1[df.1$GrowthTemp==22,]) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1, size=0.3, data=df.1[df.1$GrowthTemp==16,]) +
  geom_point(data=df.1[df.1$GrowthTemp==16,]) +
  theme_linedraw() + theme(strip.background = element_blank(), strip.text = element_text(color="black"),
                           panel.grid = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values=c("blue3","firebrick1"), labels=c("16°C","22°C")) +
  scale_shape_manual(values=c(16,17), labels=c("16°C","22°C"))
p.8.1 + facet_grid(Variate ~ Substrate + Day, scales="fixed")
dev.off()

# All enzymes; Vmax.16 by generation for 16 and 22 evolved strains averaged over growth T
df.2 <- summarySE(ALL.EVO, "Vmax.16" ,c("Variate","Substrate","EvolTemp","Day","Generation"))
df.2$Day <- revalue(df.2$Day, c("14"="14 Days", "21" = "21 Days"))
df.2$Generation <- revalue(df.2$Generation, c("A"="0", "F" = "1500"))
pdf("Fig5.pdf", height = 10, width = 8)
p.8.2 <- ggplot(df.2, aes(Generation, mean, color=EvolTemp, shape=EvolTemp)) + ylab("log(Vmax)") + labs(color="Evolution T", shape="Evolution T") +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1, size=0.3, data=df.2[df.2$EvolTemp==22,]) +
  geom_point(data=df.2[df.2$EvolTemp==22,]) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1, size=0.3, data=df.2[df.2$EvolTemp==16,]) +
  geom_point(data=df.2[df.2$EvolTemp==16,]) +
  theme_linedraw() + theme(strip.background = element_blank(), strip.text = element_text(color="black"),
                           panel.grid = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values=c("blue3","firebrick1"), labels=c("16°C","22°C")) +
  scale_shape_manual(values=c(16,17), labels=c("16°C","22°C"))
p.8.2 + facet_grid(Variate ~ Substrate + Day, scales="fixed")
dev.off()

# All enzymes; Vmax.16 by evolution temperature
for (enz in c("AG","AP","BG","BX","CBH","LAP","NAG","PPO","OX","PER")) {
  df.1 <- summarySE(F.ALL[F.ALL$Variate==enz,], "Vmax.16" ,c("Substrate","GrowthTemp","Day","EvolTemp"))
  pdf(paste0("Vmax.16.",enz,".pdf"))
  p.9 <- ggplot(df.1, aes(GrowthTemp, mean, color=EvolTemp)) + ylab(paste0(enz," Vmax.16")) +
    geom_point() +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1, size=0.3) +
    scale_color_manual(values=c("blue3","firebrick1"))
  print(p.9 + facet_grid(Substrate ~ Day, scales="fixed"))
  dev.off()
}

# All enzymes; Vmax.16 by generation
for (enz in c("AG","AP","BG","BX","CBH","LAP","NAG","PPO","OX","PER")) {
  df.1 <- summarySE(ALL.16[ALL.16$Variate==enz,], "Vmax.16" ,c("Substrate","GrowthTemp","Day","Generation"))
  df.1$GrowthTemp <- factor(df.1$GrowthTemp)
  pdf(paste0("Gen.Vmax.16.",enz,".pdf"))
  p.10 <- ggplot(df.1, aes(GrowthTemp, mean, color=Generation)) + ylab(paste0(enz," Vmax.16")) +
    geom_point() +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.1, size=0.3) +
    scale_color_manual(values=c("black","orange"))
  print(p.10 + facet_grid(Substrate ~ Day, scales="fixed"))
  dev.off()
}

# Testing NAG for effects on Vmax.Slope
# Corresponds to Fig. 6
sink("ANOVA.AF.txt")
m.7 <- lme(min_rank(Vmax.Slope)~Substrate*GrowthTemp*EvolTemp*Day,random=~1|Strain,data=F.ALL[F.ALL$Variate=="NAG",],na.action=na.omit)
rs <- resid(m.7)
print("NAG Vmax.Slope")
print(shapiro.test(rs))
hist(rs)
anova.lme(m.7,type="marginal")

for (enz in c("AG","BG","CBH")) {
  print(paste(enz,"Vmax.Slope"))
  m.7.1 <- lme(min_rank(Vmax.Slope)~GrowthTemp*EvolTemp*Day,random=~1|Strain,data=F.ALL[F.ALL$Variate==enz & F.ALL$Substrate=="Lignin",],na.action=na.omit)
  rs <- resid(m.7.1)
  print(shapiro.test(rs))
  hist(rs)
  print(anova.lme(m.7.1,type="marginal"))
}

# Testing all enzymes for evolution temperature effect on Vmax
# LAP and oxidase models fail due to missing data
for (enz in c("AG","AP","BG","BX","CBH","NAG")) {
  print(paste(enz,"Vmax.16"))
  m.8 <- lme(min_rank(Vmax.16)~Substrate*GrowthTemp*EvolTemp*Day,random=~1|Strain,data=F.ALL[F.ALL$Variate==enz,],na.action=na.omit)
  rs <- resid(m.8)
  print(shapiro.test(rs))
  hist(rs)
  print(anova.lme(m.8,type="marginal"))
}

# On lignin: EvolTemp significant for CBH
for (enz in c("AG","AP","BG","BX","CBH","NAG","PPO")) {
  print(paste(enz,"Vmax.16 Lignin"))
  m.8.1 <- lme(min_rank(Vmax.16)~GrowthTemp*EvolTemp*Day,random=~1|Strain,data=F.ALL[F.ALL$Variate==enz & F.ALL$Substrate=="Lignin",],na.action=na.omit)
  rs <- resid(m.8.1)
  print(shapiro.test(rs))
  hist(rs)
  print(anova.lme(m.8.1,type="marginal"))
}

# On sucrose: EvolTemp by Day interaction for AP, EvolTemp by GrowthTemp for NAG (but ns in m.8)
for (enz in c("AG","AP","BG","CBH","LAP","NAG")) {
  print(paste(enz,"Vmax.16 Sucrose"))
  m.8.2 <- lme(min_rank(Vmax.16)~GrowthTemp*EvolTemp*Day,random=~1|Strain,data=F.ALL[F.ALL$Variate==enz & F.ALL$Substrate=="Sucrose",],na.action=na.omit)
  rs <- resid(m.8.2)
  print(shapiro.test(rs))
  hist(rs)
  print(anova.lme(m.8.2,type="marginal"))
}

# Testing all enzymes for generation effect on Vmax
# OX and PPO models fail due to missing data
for (enz in c("AG","AP","BG","BX","CBH","LAP","NAG","PER")) {
  print(paste(enz,"Vmax.16"))
  m.9 <- lme(min_rank(Vmax.16)~Substrate*GrowthTemp*Generation*Day,random=~1|Strain,data=ALL.16[ALL.16$Variate==enz,],na.action=na.omit)
  rs <- resid(m.9)
  print(shapiro.test(rs))
  hist(rs)
  print(anova.lme(m.9,type="marginal"))
}

# On lignin: Generation significant for CBH and oxidases; interaction with Day for NAG
for (enz in c("AG","AP","BG","BX","CBH","LAP","NAG","PPO","OX","PER")) {
  print(paste(enz,"Vmax.16 Lignin"))
  m.9.1 <- lme(min_rank(Vmax.16)~GrowthTemp*Generation*Day,random=~1|Strain,data=ALL.16[ALL.16$Variate==enz & ALL.16$Substrate=="Lignin",],na.action=na.omit)
  rs <- resid(m.9.1)
  print(shapiro.test(rs))
  hist(rs)
  print(anova.lme(m.9.1,type="marginal"))
}

# On sucrose: Generation significant for AG, BX, CBH, PER; interactions for AP, LAP, NAG
for (enz in c("AG","AP","BG","BX","CBH","LAP","NAG","PER")) {
  print(paste(enz,"Vmax.16 Sucrose"))
  m.9.2 <- lme(min_rank(Vmax.16)~GrowthTemp*Generation*Day,random=~1|Strain,data=ALL.16[ALL.16$Variate==enz & ALL.16$Substrate=="Sucrose",],na.action=na.omit)
  rs <- resid(m.9.2)
  print(shapiro.test(rs))
  hist(rs)
  print(anova.lme(m.9.2,type="marginal"))
}
sink()

# Calculate SE from 95% CI limits
params.frame$Vmax.SE <- (params.frame$Vmax - params.frame$lowerV)/1.96
params.frame$Km.SE <- (params.frame$Km - params.frame$lowerK)/1.96

# Print out for SI Table
SI.table <- params.frame[c("Temp","Variate","Strain","Generation","Substrate","Day","GrowthTemp","EvolTemp","Vmax","Vmax.SE","Km","Km.SE")]
names(SI.table) <- c("Assay.Temperature","Enzyme","Strain","Generation","Substrate","Day","GrowthTemp","EvolTemp","Vmax","Vmax.SE","Km","Km.SE")
SI.table$Generation <- revalue(SI.table$Generation, c("A"="0", "B"="500", "C" = "1000", "F" = "1500"))
write.csv(SI.table,"Enzymes_SI.csv",quote=F,row.names = F)
