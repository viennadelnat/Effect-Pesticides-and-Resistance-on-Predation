#Resistance to a chemical pesticide increases vulnerability to a biopesticide: Effects on direct mortality and mortality by predation
#Vienna Delnat, Tam Tran, Lizanne Janssens and Robby Stoks
#Aquatic Toxicology (2019)
#R code tested on 15/03/2022


#####Packages#####

install.packages("doBy")     
install.packages("drc")     
install.packages("afex")     
install.packages("lme4")     
install.packages("lmerTest")     
install.packages("brglm")     
install.packages("car")     
install.packages("emmeans")     
install.packages("effects")     
install.packages("rptR")

library(doBy)
library(drc)
library(afex)
library(lme4)
library(lmerTest)
library(brglm)
library(car)     
library(emmeans)
library(effects)     
library(rptR)

sessionInfo()

##Set working directory to source file
#RStudio -> Session -> Set Working Directory...-> To Source File Location


#####Datasets#####

#Range finder data
dataRF=read.csv("Delnat-et-al_Range finding results.csv", sep=",", na.strings=c(""))
##Set correct data types 
#for drc:: drm, Concentration needs to be numeric (nominal concentrations are used and not measured concentrations!)
#for drc:: drm, Mortality/Survival needs to be in percentage
Factors <- c("Treatment", "Replicate", "Date", "Month","Strain", "Marked")
dataRF[Factors] <- do.call(cbind.data.frame, lapply(dataRF[Factors], as.factor))
Numerics <- c("Concentration", "Mortality48hPercent")
dataRF[Numerics] <- do.call(cbind.data.frame, lapply(dataRF[Numerics], as.numeric))
str(dataRF)
#Subsets
dataRFbti=subset(dataRF, Treatment!="CPF" & Strain=="S-Lab")
dataRFcpf=subset(dataRF, Treatment!="Bti" & Strain=="S-Lab")

#Pesticide data
dataPesticide50=read.csv("Delnat-et-al_Pesticide_After50percent_glm.csv", sep=",", na.strings=c(""))
##Set correct data types 
Factors <- c("StartDate", "Generation", "Treatment", "CPF", "Bti", "Strain", "Replicate")
dataPesticide50[Factors] <- do.call(cbind.data.frame, lapply(dataPesticide50[Factors], as.factor))
##Set levels in factor
dataPesticide50$Treatment=factor(dataPesticide50$Treatment, levels=c("Control", "Bti", "CPF", "Mixture"))
dataPesticide50$Strain=factor(dataPesticide50$Strain, levels=c("S-Lab", "Ace-1R"))
str(dataPesticide50)

#Delayed mortality data 
dataDelayed=read.csv("Delnat-et-al_DelayedMortality_glm.csv", sep=",", na.strings=c(""))
##Set correct data types 
Factors <- c("DelayedDate", "Strain", "Treatment", "CPF", "Bti", "Replicate")
dataDelayed[Factors] <- do.call(cbind.data.frame, lapply(dataDelayed[Factors], as.factor))
##Set levels in factor
dataDelayed$Treatment=factor(dataDelayed$Treatment, levels=c("control", "Bti", "CPF", "mixture"))
dataDelayed$Strain=factor(dataDelayed$Strain, levels=c("S-Lab", "Ace-1R"))
str(dataDelayed)

#Predation data
dataPredation50=read.csv("Delnat-et-al_Predation_After50percent_wilcox.csv", sep=",", na.strings=c(""))
##Set correct data types 
Factors <- c("Generation", "Day", "Treatment", "CPF", "Bti", "Replicate")
dataPredation50[Factors] <- do.call(cbind.data.frame, lapply(dataPredation50[Factors], as.factor))
dataPredation50$Treatment=factor(dataPredation50$Treatment, levels=c("Control", "Bti", "CPF", "Mixture"))
str(dataPredation50)
#Subsets
dataPredation50CTR=subset(dataPredation50, Treatment=="Control")
dataPredation50BTI=subset(dataPredation50, Treatment=="Bti")
dataPredation50CPF=subset(dataPredation50, Treatment=="CPF")
dataPredation50MIX=subset(dataPredation50, Treatment=="Mixture")

#Freezing behaviour
dataFreeze=read.csv("Delnat-et-al_FreezingBehaviour.csv", sep=",", na.strings=c(""))
##Set correct data types 
Factors <- c("BehaviorDate", "Strain", "Treatment", "CPF", "Bti", "Replicate", "Observer")
dataFreeze[Factors] <- do.call(cbind.data.frame, lapply(dataFreeze[Factors], as.factor))
##Set levels in factor
dataFreeze$Treatment=factor(dataFreeze$Treatment, levels=c("control", "Bti", "CPF", "mixture"))
dataFreeze$Strain=factor(dataFreeze$Strain, levels=c("S-Lab", "Ace-1R"))
str(dataFreeze)

#Swimming Speed
dataSwim=read.csv("Delnat-et-al_SwimmingSpeed.csv", sep=",", na.strings=c(""))
##Set correct data types 
Factors <- c("BehaviorDate", "Strain", "Treatment", "CPF", "Bti", "Replicate", "Tapped")
dataSwim[Factors] <- do.call(cbind.data.frame, lapply(dataSwim[Factors], as.factor))
##Set levels in factor
dataSwim$Treatment=factor(dataSwim$Treatment, levels=c("Control", "Bti", "CPF", "Mixture"))
dataSwim$Strain=factor(dataSwim$Strain, levels=c("S-Lab", "Ace-1R"))
str(dataSwim)

#Swimming Speed - repeated measures
dataSwimRepeat=read.csv("Delnat-et-al_SwimmingSpeedRepeat.csv", sep=",", na.strings=c(""))
##Set correct data types 
Factors <- c("BehaviorDate", "Strain", "Treatment", "CPF", "Bti", "Replicate", "Tapped")
dataSwimRepeat[Factors] <- do.call(cbind.data.frame, lapply(dataSwimRepeat[Factors], as.factor))
##Set levels in factor
dataSwimRepeat$Treatment=factor(dataSwimRepeat$Treatment, levels=c("Control", "Bti", "CPF", "Mixture"))
dataSwimRepeat$Strain=factor(dataSwimRepeat$Strain, levels=c("S-Lab", "Ace-1R"))
str(dataSwimRepeat)


#####Range finder Bti & CPF - drc#####

#drc and EC10,48h - Bti
EC48bti <- drc::drm(Mortality48hPercent ~ Concentration, data = dataRFbti, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "EC50")))
summary(EC48bti)
ED(EC48bti, c(10),interval='delta')
png("Figure_EC48bti.png", width = 9, height = 6, unit = "cm", res=300)
plot(EC48bti, type="obs", log="x", pch=16,cex=0.8, ylab="Mortality (%)", xlab="Bti concentration (log)")
plot(EC48bti, type="confidence", log="x", add=TRUE)
dev.off()

#drc and EC10,48h - CPF
EC48cpf <- drc::drm(Mortality48hPercent ~ Concentration, data = dataRFcpf, fct = LL.4(names=c("Slope","Lower Limit","Upper Limit", "EC50")))
summary(EC48cpf)
ED(EC48cpf, c(10),interval='delta')
png("Figure_EC48cpf.png", width = 9, height = 6, unit = "cm", res=300)
plot(EC48cpf, type="obs", log="x", pch=16,cex=0.8, ylab="Mortality (%)", xlab="Chlorpyrifos concentration (log)")
plot(EC48cpf, type="confidence", log="x", add=TRUE)
dev.off()


#####Pesticide#####

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#Generalized linear mixed models with a binomial error structure and the logit link
#Correct for pseudoreplication (Replicate = Vial/Container) and start date (StartDate) by adding them as random factor
glmerPesticide=glmer(Mortality48h ~ (Strain+Bti+CPF)^3 + (1|Replicate) + (1|StartDate), data=dataPesticide50, na.action=na.omit, family=binomial(link=logit),
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerPesticide, type="III") 

#Posthoc test - contrast analysis with fdr correction
interact1<-pairs(emmeans(glmerPesticide, ~Strain|Bti*CPF, adjust="none"))
interact2<-pairs(emmeans(glmerPesticide, ~Bti*CPF|Strain, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#Quick effects plot - not used in manuscript
plot(effect(mod=glmerPesticide, term="Strain*Bti*CPF"),type = "response")
#Emmeans and standard errors for figure
PesticidePlotData <- summary(emmeans(glmerPesticide, ~ Strain*Bti*CPF, type = "response"))

#Assumption - Dispersion parameter
glmPesticide=glm(Mortality48h~(Strain+Bti+CPF)^3, data=dataPesticide50, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmPesticide) 


#####Delayed Mortality#####

#When using set_sum_contrasts() all P-values > 0.99 and all Chisq 0e+00; normal values when using set_treatment_contrasts()
set_treatment_contrasts()

#Generalized linear mixed models with a binomial error structure and the logit link
#Correct for pseudoreplication (Replicate = Vial/Container) and start date (StartDate) by adding them as random factor
glmerDelayedMortality=glmer(Mortality48h ~ (Strain+Bti+CPF)^3 + (1|Replicate) + (1|DelayedDate), data=dataDelayed, na.action=na.omit, family=binomial(link=logit),
                            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerDelayedMortality, type="III") 

#Quick effects plot - not used in manuscript
plot(effect(mod=glmerDelayedMortality, term="Strain*Bti*CPF"),type = "response")
#Emmeans and standard errors for figure
DelayedMortalityPlotData <- summary(emmeans(glmerDelayedMortality, ~ Strain*Bti*CPF, type = "response"))

#Assumption - Dispersion parameter
glmDelayedMortality=glm(Mortality48h~(Strain+Bti+CPF)^3, data=dataDelayed, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmDelayedMortality) 


#####Predation#####

#non-parametric test (no interaction in model)
set_treatment_contrasts()

#Separate one-sided Wilcoxon signed rank tests for each pesticide treatment
wilcox.test(dataPredation50CTR$Proportion, mu=1, alternative = "less", conf.int=TRUE, exact=FALSE)
wilcox.test(dataPredation50BTI$Proportion, mu=1, alternative = "less", conf.int=TRUE, exact=FALSE)
wilcox.test(dataPredation50CPF$Proportion, mu=1, alternative = "greater", conf.int=TRUE, exact=FALSE)
wilcox.test(dataPredation50MIX$Proportion, mu=1, alternative = "greater", conf.int=TRUE, exact=FALSE)

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#General linear model with a normal error structure and the identity link
lmDuration=lm(HourCount ~ Bti*CPF, data=dataPredation50, na.action=na.omit) 
Anova(lmDuration, type="III") 

#Quick effects plot - not used in manuscript
plot(effect(mod=lmDuration, term="Bti*CPF"))
#Emmeans and standard errors for figure
DurationPlotData <- summary(emmeans(lmDuration, ~ Bti*CPF, type = "response"))

#Assumption - Normality of residuals
shapiro.test(resid(lmDuration))                  
hist(resid(lmDuration))    
#Assumption - Homogeneity of variance
leveneTest(HourCount ~ Bti*CPF, data = dataPredation50)

#Outliers and influential observations
outlierTest(lmDuration)
cd=cooks.distance(lmDuration); which(cd>1)
influenceIndexPlot(lmDuration, vars = c("studentized", "Bonf"))


#####Freezing behaviour#####

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#Normality not ok, try Box-Cox transformation
BoxCox=lm(DurationTotal ~ Strain*Bti*CPF, data=dataFreeze, na.action=na.omit) 
#Plot profile Log-likelihood: range determined by 'seq' --> change the range till you find the peak
car::boxCox(BoxCox, family="yjPower", plotit=TRUE, lambda = seq(0.08,0.18 , length = 10))
#Plot profile Log-likelihood: peak = 0.122 <-- power used in box cox transformation
bcFreeze <- car::yjPower(dataFreeze$Duration5min, 0.122)
#Add Box-Cox transformed Freeze to your dataset
dataFreeze$bcFreeze=bcFreeze
#Remark, by accident DurationTotal is used to determine the peak (power=0.122)
#However, not a problem as assumption of normality was met and if 0.20 as power would have been used, the same significance levels would have been achieved
# BoxCox=lm(Duration5min ~ Strain*Bti*CPF, data=dataFreeze, na.action=na.omit) 
# car::boxCox(BoxCox, family="yjPower", plotit=TRUE, lambda = seq(0.18,0.22 , length = 10))
# bcFreeze <- car::yjPower(dataFreeze$Duration5min, 0.20)

#General linear mixed model with a normal error structure and the identity link
lmerFreeze=lmer(bcFreeze ~ Strain*Bti*CPF + (1|BehaviorDate)+ (1|Replicate), data=dataFreeze, na.action=na.omit,
            control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5))) 
Anova(lmerFreeze, type="III") 

#Posthoc test - contrast analysis with fdr correction
interact1<-pairs(emmeans(lmerFreeze, ~Strain|Bti, adjust="none"))
interact2<-pairs(emmeans(lmerFreeze, ~Bti|Strain, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")
#A priori contrasts - specifically test for both strains whether the CPF-Bti mixture differed from the single pesticide exposures and solvent control
interact1<-pairs(emmeans(lmerFreeze, ~Strain|Bti*CPF, adjust="none"))
interact2<-pairs(emmeans(lmerFreeze, ~Bti*CPF|Strain, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#Quick effects plot - not used in manuscript
plot(effect(mod=lmerFreeze, term="Strain*Bti*CPF"))
plot(effect(mod=lmerFreeze, term="Strain*Bti"))
plot(effect(mod=lmerFreeze, term="CPF"))
#Emmeans and standard errors for figure
FreezePlotData <- summary(emmeans(lmerFreeze, ~ Strain*Bti*CPF, type = "response"))

#Assumption - Normality of residuals
shapiro.test(resid(lmerFreeze))                  
hist(resid(lmerFreeze))    
#Assumption - Homogeneity of variance
leveneTest(bcFreeze ~ Strain*Bti*CPF, data = dataFreeze)
#Thumb of rule - if minimum and maximum variance do not differ more than a factor 5 - assumption still met
aggregate(bcFreeze ~ Strain*Bti*CPF, data = dataFreeze, var)

#Outliers and influential observations
outlierTest(lmerFreeze)
cd=cooks.distance(lmerFreeze); which(cd>1)
influenceIndexPlot(lmerFreeze, vars = c("studentized", "Bonf"))


#####Swimming speed - 10 first frames#####

#Interaction in model --> use set_sum_contrasts() and type=3 in Anova
set_sum_contrasts()

#General linear mixed model with a normal error structure and the identity link
#Larvae started swimming after being tapped (~68%) or before being tapped (~32%), to correct for this tap it was added as random factor (this was by accident not mentioned in Mat&Meth in manuscript)
lmerSwimmingSpeed10fr=lmer(SwimmingSpeed10fr ~ Strain*Bti*CPF + (1|BehaviorDate) + (1|Tapped) + (1|Replicate), data=dataSwim, 
                           na.action=na.omit, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5))) 
Anova(lmerSwimmingSpeed10fr, type="III") 

#Posthoc test - contrast analysis with fdr correction
interact1<-pairs(emmeans(lmerSwimmingSpeed10fr, ~Strain|CPF, adjust="none"))
interact2<-pairs(emmeans(lmerSwimmingSpeed10fr, ~CPF|Strain, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")
#A priori contrasts - specifically test for both strains whether the CPF-Bti mixture differed from the single pesticide exposures and solvent control
interact1<-pairs(emmeans(lmerSwimmingSpeed10fr, ~Strain|Bti*CPF, adjust="none"))
interact2<-pairs(emmeans(lmerSwimmingSpeed10fr, ~Bti*CPF|Strain, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

#Quick effects plot - not used in manuscript
plot(effect(mod=lmerSwimmingSpeed10fr, term="Strain*Bti*CPF"))
plot(effect(mod=lmerSwimmingSpeed10fr, term="Strain*CPF"))
#Emmeans and standard errors for figure
SwimmingSpeedPlotData <- summary(emmeans(lmerSwimmingSpeed10fr, ~ Strain*Bti*CPF, type = "response"))

#Assumption - Normality of residuals
shapiro.test(resid(lmerSwimmingSpeed10fr))                  
hist(resid(lmerSwimmingSpeed10fr))    
#Assumption - Homogeneity of variance
leveneTest(SwimmingSpeed10fr ~ Strain*Bti*CPF, data = dataSwim)

#Outliers and influential observations
outlierTest(lmerSwimmingSpeed10fr)
cd=cooks.distance(lmerSwimmingSpeed10fr); which(cd>1)
influenceIndexPlot(lmerSwimmingSpeed10fr, vars = c("studentized", "Bonf"))

#repeatability - R = 0.43 en p<0.001 dus ok
REPSwimmingSpeed10fr <- rpt(SwimmingSpeed10fr ~   (1 |ID), grname = "ID", data = dataSwimRepeat, 
                            datatype = "Gaussian", nboot = 1000, npermut = 0)
print(REPSwimmingSpeed10fr) 


######Save Rdata######
save.image(file="Chapter2_Rdata_20220315_NotPublished.Rdata")
