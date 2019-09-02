library("ggplot2")
library("lmerTest") # calculates p-values for fixed effect
library("lattice")
library("cowplot")

### Preferences ###
options(digits.secs=6)

######################### 1. FILE I/O ################################
n = 31
ids <- 1:10
ids2 <- 11:31 # from second collection phase
gend <<- c("M","F","M","M","F","F","F","M","F","F","F", "F","M","M","F","F","M","M","M","M","M","F","M","M","M","F","F","F","F","M","F")
logs <- c("log00.csv", "log01.csv", "log02.csv", "log03.csv","log04.csv",
          "log06.csv", "log07.csv", "log19.csv", "log22.csv", "log23.csv")
mfiles <- paste("markers_",sep="",logs)
logs2 <- c("P1.csv", "P3.csv", "P4.csv", "P5.csv", "P6.csv", "P7.csv", "P8.csv","P9.csv",
           "P10.csv","P11.csv","P12.csv","P14.csv","P16.csv","P17.csv","P18.csv","P19.csv","P20.csv","P21.csv","P22.csv","P23.csv","P24.csv")
mfiles2 <- c("P1markers.csv", "P3markers.csv", "P4markers.csv", "P5markers.csv", "P6markers.csv", "P7markers.csv", "P8markers.csv","P9markers.csv",
             "P10markers.csv", "P11markers.csv", "P12markers.csv", "P14markers.csv", "P16markers.csv", "P17markers.csv", "P18markers.csv",
             "P19markers.csv", "P20markers.csv", "P21markers.csv", "P22markers.csv", "P23markers.csv", "P24markers.csv")
mapping <<- read.csv("stimuliOrder.csv", header=TRUE) # has metadata about trials based on ID
names(mapping)
mapping$trialID <- factor(mapping$trialID)
excludedTrial <- c(7,10,11,14,19)  # these trials, performed in the first study phase, were too similar to other body parts so excluded for simplicity
bp_order <- (c("Genitals","Butt","Breast","Heart","Inner thigh","Foot","Eye","Ear","Back","Neck","Forehead","Arm","Hand"))
bar_order <- c("lo","med","hi")

############## 1.1 Event times from a Q-sensor marker file ####################
keepVars <- c("pID", "gender", "condition", "bodyPart", "bar",  "phase", "time")
clean <- function(mfile, id, session=1){
  m <- read.csv(mfile, header=TRUE, skip=1)
  m <- m[!is.na(m$trial),]   # remove non-event markers
  m$time <- strptime(m$Start,"%H:%M:%OS") # convert to POSIXt format
  m$phase <- c("prompt", "action", "response")   # categorize markers
  m$pID <- factor(id)
  m$gender <- gend[id]
  if(session==1){
    m <- merge(m, mapping, by.x = "trial") 
  }else{
    m$trialID2 <- m$trial
    m <- merge(m, mapping, by = "trialID2")
  }
  m <- subset(m, !(m$trialID %in% excludedTrial), drop = T)
  m <- m[,keepVars]
  m$bodyPart <- factor(m$bodyPart, levels=bp_order, ordered = T)
  m$bar <- factor(m$bar, levels=bar_order, ordered = T)
  
  return(m)
}
m <- mapply(clean, mfiles, ids) # apply function to all files
m2 <- mapply(clean, mfiles2, ids2, 2) # apply function to all files 
m <- apply(m, MARGIN = 2 , FUN=as.data.frame) # make a list of 1 dataframe / id
m2 <- apply(m2, MARGIN = 2 , FUN=as.data.frame)
m <- do.call(rbind, m) # make single sheet, long form; can also call: Reduce(function(...) merge(...,all=T), m)
m2 <- do.call(rbind, m2)
m <- rbind(m, m2) # join m's from collection phases
nrow(m) == n * 3 * 26 # verify the correct number of trial

################# 1.2 EDA measurements from a Q-sensor log file ################
keepVars <- c("pID","time","EDA.uS.")
clean_log <- function(efile, id, sec=F){ # if sec=T, aggregate data to 1 sample per second
  e <- read.csv(efile, header=TRUE, skip=6)[-1,]
  e$time <- strptime(e$Time, "%H:%M:%OS") # convert to POSIXt format
  if(sec==T){
    e$cut <- as.POSIXct(cut.POSIXt(e$time,"sec", include.lowest = T))
    e <- aggregate( EDA.uS. ~ cut, data = e, mean )
    colnames(e)[colnames(e)=="cut"] <- "time"
  }
  e$pID <- factor(id)
  e <- e[,keepVars]
  return (e)
}

e <- mapply(clean_log, efile = logs, id = ids, sec=T)
e <- apply(e, MARGIN = 2 , FUN=as.data.frame)
e <- do.call(rbind, e)
e2 <- mapply(clean_log, efile = logs2, id = ids2, sec=T)
e2 <- apply(e2, MARGIN = 2 , FUN=as.data.frame)
e2 <- do.call(rbind, e2)
e <- rbind(e, e2) # join e's from collection phases

str(m)
str(e)

m$EDA.uS. <- NA
for(p in unique(m$pID) ){
  markers <- subset(m, pID==p)
  EDAs <- subset(e, pID==p)
  
  # determine the EDA at each of the marked times using look-up (ref: https://sapa-project.org/blog/2013/06/28/vlookup-in-R/)
  m[m$pID==p, "EDA.uS."] <- EDAs$EDA.uS.[match(markers$time, EDAs$time)]
  # if looking at action-response,
  # m[m$pID==p & m$phase=="response", "EDA.uS."] <- EDAs$EDA.uS.[match(markers$time[markers$phase=="action"]+3, EDAs$time)]
}

# calculate IVs
deltaEDA <- m[which(m$phase=="action"),"EDA.uS."] - m[which(m$phase=="prompt"),"EDA.uS."]
responseTime <- m[which(m$phase=="action"),"time"] - m[which(m$phase=="prompt"),"time"]
dat <- cbind(m[which(m$phase=="action"),!colnames(m) %in% c("phase", "EDA.uS.","time")], deltaEDA, responseTime)

# write summary data file
write.csv(dat, "summary.csv")
save(dat, file="dat.RData")

# analysis

# outliers - outlier value, sd's from mean, replace with mean (replacement value)
dat$deltaEDA[dat$bodyPart=="Genitals" & dat$condition=="touch" & dat$pID == 1]
summary(dat$deltaEDA[dat$bodyPart=="Genitals" & dat$condition=="touch" & dat$pID != 1])[4] - 2.5 * sd(dat[dat$bodyPart=="Genitals" & dat$condition=="touch"  & dat$pID != 1,"deltaEDA"], na.rm=T)
m[m$bodyPart == "Genitals" & m$condition=="touch" & m$pID == 1 & m$phase=="prompt", "EDA.uS."] + summary(dat$deltaEDA[dat$bodyPart=="Genitals" & dat$condition=="touch" & dat$pID != 1])[4]
# 5.309094
dat$deltaEDA[dat$bodyPart=="Heart" & dat$condition=="touch" & dat$pID == 1]
summary(dat$deltaEDA[dat$bodyPart=="Heart" & dat$condition=="touch" & dat$pID != 1])[4] - 5 * sd(dat[dat$bodyPart=="Heart" & dat$condition=="touch" & dat$pID != 1,"deltaEDA"], na.rm=T)
m[m$bodyPart == "Heart" & m$condition=="touch" & m$pID == 1 & m$phase=="prompt", "EDA.uS."] + summary(dat$deltaEDA[dat$bodyPart=="Heart" & dat$condition=="touch" & dat$pID != 1])[4]
dat$responseTime[dat$bodyPart=="Genitals" & dat$condition=="touch" & dat$pID == 5]


# aggregate data over body parts for each accesibility level
datbar <- aggregate(dat$deltaEDA, by = list(dat$pID, dat$condition, dat$bar), FUN = 'mean')
colnames(datbar) <- c("pID","condition","bar", "deltaEDA")

datbar2 <- aggregate(dat$responseTime, by = list(dat$pID, dat$condition, dat$bar), FUN = 'mean')
colnames(datbar2) <- c("pID","condition","bar", "responseTime")

# normality test
(normTest <- aggregate(cbind(P.value=deltaEDA) ~ condition, datbar,
                       FUN = function(x) shapiro.test(x)$p.value))

# RM-ANOVA
summary(aov(deltaEDA ~ condition * bar + Error(pID / (condition * bar)), data=datbar))

# linear mixed-effects model fit by maximum likelihood, with participant as a random effect. Followed by contrasts
lmer1 <- lmer(deltaEDA ~ condition * bar + (1| pID), data=datbar, REML=FALSE, contrasts="Simple") # .L .Q because ordinal

mat <- rbind(c(0,1,-1),c(1,0,-1)) # c(1,0,-1) is significant: low accessibilty vs high accessibility is different
require(MASS)
cMat <- ginv(mat)
lmer1.contrast <- lmer(deltaEDA ~ condition * bar + (1| pID), data=datbar, REML=FALSE, contrasts=list(bar=cMat))

mat <- rbind(c(0,1))
cMat <- ginv(mat)
lmer1.contrast2 <- lmer(deltaEDA ~ condition * bar + (1| pID), data=datbar, REML=FALSE, contrasts=list(condition=cMat))

summary(lmer1)
summary(lmer1.contrast)
confint(lmer1,method = "Wald")
dotplot(ranef(lmer1, condVar = TRUE))

summary(lmer(as.numeric(responseTime) ~ condition * bar + (1| pID), data=datbar2, REML=FALSE))

# Notes:

# 1. Statistical test
#    RM ANOVA using the aov function and either aggregated or unaggregated data gives similar results
#    sample code: aov(deltaEDA ~ condition * bar + Error(pID / (condition * bar)), data=dat)
# 2. Aggregation vs. non-aggregated data:
#    should we average the body parts into their respective accesibility groups?
#    LMER without aggregation, bodyPart is a RE. estimating the influence of all lo, med, hi accesibility body parts
#    when testing using lmer without aggregation, it estimates a different mean that's different than the actual mean of the 4 body parts because 
#    body part is a random effect; sample size for low curve fit based on the 4 points, so can't say that much about the curve
#    therefore, have fewer DoF b/c DoF is low when there's few observations
#    usu DoF = # of people - # of conditions
#    if millions of people, b0*lo + b1*med + b2*high + ai (P RE) + gi (body part RE); both ai and gi approach 0
#    report the estimates (for the population) from the MODEL, because it takes into account the variability, M(data), it's M(population)
#    report mean from data, M(population SE is lower) is simplier, M(population) matched the model
#    m(touchhi) = (Intercept) + conditiontouch + barhi + conditiontouch:barhi = M(touchhi)  (==contrast)
#    bar1 in 2nd model = barmed - barhi  (==contrast)
#    Aggregation test code: uses data = datbar instead of datbar
lmer.noagg <- lmer(deltaEDA ~ condition * bar + (1 | pID), data=dat, REML=FALSE)  # LMER without aggregation
lmer.aggregation <- summary(lmer(deltaEDA ~ condition * bar + (1| pID), data=datbar, REML=FALSE)) # LMER with aggregation

# table of values

means <- cbind(aggregate(dat[, c("deltaEDA","responseTime")], list(dat$bodyPart, dat$condition), mean), # could add dat$condition to list
               aggregate(dat[, c("deltaEDA","responseTime")], list(dat$bodyPart, dat$condition), sd)[,3:4])
(means <- means[,c(2,1,3,5,4,6)])
aggregate(dat[, c("deltaEDA","responseTime")], list(dat$bar, dat$condition), mean)


# plots

library(ggplot2)
require(scales)

science_theme = theme(panel.grid.major = element_line(size = 0.5, color = "grey"), 
                      axis.line.x = element_line(size = 0.7, color = "black"), axis.line.y = element_line(size = 0.7, color = "black"),
                      legend.position = c(0.87, 0.92), text = element_text(size = 14)) # reduce 0.92: move legend down, reduce 0.87: move right

#(plotdat <- summarySE(dat, measurevar="deltaEDA", groupvars=c("bar","condition"), na.rm=TRUE)) # if we want CIs instead of SE
png(filename="figure1.png", 
    type="cairo",
    units="in", 
    width=5, 
    height=5,  
    pointsize=10, 
    res=800)
fig1 <- ggplot(dat, aes(x = bar, y = deltaEDA, group = condition) ) +
  geom_errorbar(aes(color=condition, linetype=condition), size = 1.3, position = position_dodge(width=.9), data=dat,fun.data = mean_se,stat = 'summary') +
  geom_point(aes(color=condition, shape=condition, linetype=condition), size = 2.6, position = position_dodge(width=.9), data=dat,fun.data = mean_se,stat = 'summary') +
  scale_linetype_manual(values=c("solid", "solid")) +
  scale_color_manual(values=c('#999999','#E69F00')) +
  ylab(expression(paste("Change in physiological arousal (",mu, "S)"))) +
  xlab("Accessibility of body part on robot") +
  coord_flip() +
  theme_bw() +
  theme(axis.title.y = element_text(face="plain"),
        legend.title=element_blank(), legend.key = element_rect(color = "white", fill = NA, size=1), legend.key.width = unit(0.4, "in"),
        panel.grid.major.y = element_blank()) +
  science_theme +
  annotate("text", x=3.3, y=0.005, label= "A", size=16) +
  guides(color = guide_legend(reverse=T), linetype = guide_legend(reverse=T), shape = guide_legend(reverse=T))
#guides(shape = guide_legend(keywidth = 4, keyheight = 1))
dev.off()

png(filename="figure2.png", 
    type="cairo",
    units="in", 
    width=5, 
    height=10,  
    pointsize=10, 
    res=800)
fig2 <- ggplot(dat, aes(x = interaction(pID), y = deltaEDA, group = bar) ) +
  geom_errorbar(aes(color=bar), size = 0.6, position = position_dodge(width=.9), data=dat,fun.data = mean_se,stat = 'summary') +
  geom_point(aes(color=bar, shape=bar), size = 1.2, position = position_dodge(width=.9), data=dat,fun.data = mean_se,stat = 'summary') +
  facet_wrap(~condition) +
  scale_color_manual(values=c('#FFC500','#FF5627','#CB0035')) +
  ylab(expression(paste("Change in physiological arousal (",mu, "S)"))) +
  xlab("Participant") +
  coord_flip() +
  theme_bw() +
  science_theme +
  #annotate("text", x=31, y=0.005, label= "B", size=16) +
  theme(axis.title.y = element_text(face="plain"),
        legend.title=element_blank(), legend.key = element_rect(color = "white", fill = NA, size=1),
        legend.key.width = unit(0.4, "in"), legend.position = c(0.89, 0.955),
        panel.grid.major.y = element_blank()) +
  geom_text(data=data.frame(x=30.33, y=-0.5, label=c("B"," "), 
                            condition=c("point","touch")), 
            aes(x,y,label=label), inherit.aes=FALSE, size=16) +
  guides(color = guide_legend(reverse=T), shape = guide_legend(reverse=T))
dev.off()

require(gridExtra)
png(filename="figure1combined2.png", 
    type="cairo",
    units="in", 
    width=15, 
    height=5,  
    pointsize=10, 
    res=800)
grid.arrange(fig1, fig2, ncol=2, heights=c(2, 1))
dev.off()

# generate supplmental files - full EDA information
ebig <- mapply(clean_log, efile = logs, id = ids, sec=F)
ebig<- apply(ebig, MARGIN = 2 , FUN=as.data.frame)
ebig <- do.call(rbind, ebig)
e2big <- mapply(clean_log, efile = logs2, id = ids2, sec=F)
e2big <- apply(e2big, MARGIN = 2 , FUN=as.data.frame)
e2big <- do.call(rbind, e2big)
ecombinedbig <- rbind(ebig, e2big)
ebyp <- split(ecombinedbig, ecombinedbig$pID)
mbyp <- split(m, m$pID)
t_tail <- 10
ecombinedbig2 <- lapply(names(ebyp), # remove excess measurements from before and after session
            function(p)
              subset(ebyp[[p]],
                     ebyp[[p]]$time >= min(mbyp[[p]]$time) & ebyp[[p]]$time <= max(mbyp[[p]]$time)+t_tail) # and subset based on time
)
ecombinedbig2 <- do.call(rbind, ecombinedbig2)

# add marker information to e datafile
appendVars <- c("gender", "condition", "bodyPart", "bar", "phase")
ecombinedbig2[,c(appendVars)] <- m[1,c(appendVars)]
ecombinedbig2[,"mark"] <- NA
for(p in unique(m$pID) ){
  markers <- subset(m, pID==p)
  EDAs <- subset(ecombinedbig2, pID==p)
  # determine the EDA at each of the marked times (ref: https://sapa-project.org/blog/2013/06/28/vlookup-in-R/)
  ecombinedbig2[ecombinedbig2$pID==p, "mark"] <- (EDAs$time %in% markers$time) != FALSE
  ecombinedbig2[ecombinedbig2$pID==p, appendVars] <- markers[findInterval(EDAs[,"time"], markers[,"time"]), appendVars]
}

alldat <- ecombinedbig2[,c(1,4,5,6,7,8,9,2,3)]
write.csv(ecombinedbig2, "completemeasurements.csv")
save(alldat, file="alldat.Rdata")
