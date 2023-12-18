
setwd('D:/PhD/Undergrads/Bailey_Pyle')

masses <- read.csv('masses.csv')
cercs <- read.csv('b.cercs.csv')
dub.whammy <- read.csv('b.surv.prev.csv')

### Let's look at egg masses first; will need to subset to find true infections
hist(masses$Eggs)
control <- subset(masses, Group == 'Control')
schisto <- subset(masses, Group == 'S')
E1 <- subset(masses, Group == 'E1')
E2 <- subset(masses, Group == 'E2')
E3 <- subset(masses, Group == 'E3')
T1 <- subset(masses, Group == 'T1')
T2 <- subset(masses, Group == 'T2')
T3 <- subset(masses, Group == 'T3')

Es <- rbind(E1,E2,E3)
Esm <- subset(Es, Meta.inf == 1)

s1 <- subset(schisto, Schisto.Inf == 1)
Ts <- rbind(T1,T2,T3)
Ts1 <- subset(Ts, Total.fact == 'p')

ce123 <- rbind(control, Esm)


library(glmmTMB)
fit_hzip <- glmmTMB(Eggs ~ Group + Week + (1|ID), data=ce123, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)

#check on egg masses
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom1, family=list(family="truncated_nbinom1",link="log"))
library(bbmle)
AICtab(fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)

#the best one
summary(fit_hnbinom2)

library(emmeans)
multzip <- emmeans(fit_hnbinom2, pairwise ~ Group)
test(multzip, adjust = 'none')


st123 <- rbind(s1, Ts1)
fit_hzip <- glmmTMB(Eggs ~ Group + (1|ID), data=st123, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)

#check on egg masses
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom1, family=list(family="truncated_nbinom1",link="log"))
AICtab(fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)

#the best one
summary(fit_zinbinom2)

multzip <- emmeans(fit_zinbinom2, pairwise ~ Group)
test(multzip, adjust = 'none')


cs <- rbind(control, s1)
cs.3 <- subset(cs, Week < 5)
fit_hzip <- glmmTMB(Eggs ~ Group + (1|ID), data=cs.3, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)

#check on egg masses
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom1, family=list(family="truncated_nbinom1",link="log"))
AICtab(fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)

#the best one
summary(fit_hnbinom2)


E1s <- subset(Esm, Group == 'E1')
T1s <- subset(Ts1, Group == 'T1')
ET1 <- rbind(E1s, T1s)
ET1.3 <- subset(ET1, Week < 5)
fit_hzip <- glmmTMB(Eggs ~ Group + (1|ID), data=ET1.3, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)


E2s <- subset(Esm, Group == 'E2')
T2s <- subset(Ts1, Group == 'T2')
ET2 <- rbind(E2s, T2s)
ET2.3 <- subset(ET2, Week < 5)
fit_hzip <- glmmTMB(Eggs ~ Group + (1|ID), data=ET2.3, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)


E3s <- subset(Esm, Group == 'E3')
T3s <- subset(Ts1, Group == 'T3')
ET3 <- rbind(E3s, T3s)
ET3.3 <- subset(ET3, Week < 5)
fit_hzip <- glmmTMB(Eggs ~ Group + (1|ID), data=ET2.3, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)


library("ggpubr")
ggboxplot(ce123, x = 'Group', y = 'Eggs', color = 'Group')
ggboxplot(st123, x = 'Group', y = 'Eggs', color = 'Group')
ggboxplot(cs.3, x = 'Group', y = 'Eggs', color = 'Group')
ggboxplot(ET1.3, x = 'Group', y = 'Eggs', color = 'Group')
ggboxplot(ET2.3, x = 'Group', y = 'Eggs', color = 'Group')
ggboxplot(ET3.3, x = 'Group', y = 'Eggs', color = 'Group')

ce123.g <- ggplot(ce123, aes(x = Group, y = Eggs, fill = Group)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Number of egg masses laid per week per snail", hjust=10) +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle=90, hjust=1))

colabels <- c('CON', 'C.PRE', 'C.SIM', 'C.POST')

ce123.g + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 0.5, position=position_dodge(1)) +
  scale_fill_manual(values=c("green4","orangered" ,'blue4' , 'yellow3')) + theme_classic() +
  theme(text = element_text(size=10), legend.position='none') + scale_x_discrete(labels= colabels)

st123.g <- ggplot(st123, aes(x = Group, y = Eggs, fill = Group)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Number of egg masses laid per week per snail", hjust=10) +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle=90, hjust=1))

colabels <- c('S', 'PRE', 'SIM', 'POST')

st123.g + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 0.5, position=position_dodge(1)) +
  scale_fill_manual(values=c("green4","orangered" ,'blue4' , 'yellow3')) + theme_classic() +
  theme(text = element_text(size=10), legend.position='none') + scale_x_discrete(labels= colabels)


#### check out prevalence
pT1 <- subset(dub.whammy, Treatment == 'T1')
pT2 <- subset(dub.whammy, Treatment == 'T2')
pT3 <- subset(dub.whammy, Treatment == 'T3')
pS <- subset(dub.whammy, Treatment == 'S')
pE1 <- subset(dub.whammy, Treatment == 'E1')
pE2 <- subset(dub.whammy, Treatment == 'E2')
pE3 <- subset(dub.whammy, Treatment == 'E3')

#meta prevalence
pE1.3 <- rbind(pE1, pE2, pE3)
m.prev <- glm(Meta.inf. ~ Treatment, data = pE1.3, family = 'binomial')
summary(m.prev)
multzip <- emmeans(m.prev, pairwise ~ Treatment, adjust = 'none')
multzip

pET1.3 <- rbind(pE1, pT1)
m1.prev <- glm(Meta.inf. ~ Treatment, data = pET1.3, family = 'binomial')
summary(m1.prev)

H<- c(0.7, 0.712)
par(mar=c(8.1, 5, 4.1, 2.1))
barplot(H, xlab = "PRE-control                            PRE", ylab = "Proportion of hosts with metacercarie", ylim = c(0,1), cex.axis = 2, cex.lab= 2,
        col = c("#d1495b", "#edae49"))

pET2.3 <- rbind(pE2, pT2)
m2.prev <- glm(Meta.inf. ~ Treatment, data = pET2.3, family = 'binomial')
summary(m2.prev)

H<- c(0.722, 0.605)
par(mar=c(8.1, 5, 4.1, 2.1))
barplot(H, xlab = "SIM-control                            SIM", ylab = "Proportion of hosts with metacercarie", ylim = c(0,1), cex.axis = 2, cex.lab= 2,
        col = c("#d1495b", "#edae49"))

pET3.3 <- rbind(pE3, pT3)
m3.prev <- glm(Meta.inf. ~ Treatment, data = pET3.3, family = 'binomial')
summary(m3.prev)

H<- c(0.55, 0.789)
par(mar=c(8.1, 5, 4.1, 2.1))
barplot(H, xlab = "POST-control                            POST", ylab = "Proportion of hosts with metacercarie", ylim = c(0,1), cex.axis = 2, cex.lab= 2,
        col = c("#d1495b", "#edae49"))


#schisto prev
pST1.3 <- rbind(pT1, pT2, pT3, pS)
s.prev <- glm(Schisto.inf. ~ Treatment, data = pST1.3, family = 'binomial')
summary(s.prev)
multzip <- emmeans(s.prev, pairwise ~ Treatment, adjust = 'none')
multzip

H<- c(0.45, 0.282, 0.342, 0.526)
par(mar=c(8.1, 5, 4.1, 2.1))
barplot(H, xlab = "S             PRE             SIM          POST", ylab = 'Proportion of hosts with Schistosoma', ylim = c(0,1), cex.axis = 2, cex.lab= 2,
        col = c("green4","orangered" ,'blue4' , 'yellow3'))
text(0.7, 0.8, 'a ab b', cex = 2)

#coinfection prev
pT1.3 <- rbind(pT1, pT2, pT3)
t.prev <- glm(Total.inf. ~ Treatment, data = pT1.3, family = 'binomial')
summary(t.prev)
multzip <- emmeans(t.prev, pairwise ~ Treatment, adjust = 'none')
multzip

H<- c(0.231, 0.211, 0.395)
par(mar=c(8.1, 5, 4.1, 2.1))
barplot(H, xlab = "PRE                          SIM                          POST", ylab = "Proportion of hosts with both infections", ylim = c(0,1), cex.axis = 2, cex.lab= 2,
        col = c("orangered" ,'blue4' , 'yellow3'))


#meta intensity
posE1.3 <- subset(pE1.3, Metas > 0)
m.int <- glm(Metas ~ Treatment, data = posE1.3, family = 'poisson')
summary(m.int)
multzip <- emmeans(m.int, pairwise ~ Treatment, adjust = 'none')
multzip

emeta.g <- ggplot(posE1.3, aes(x = Treatment, y = Metas, fill = Treatment)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Number of metacerciae in infected snails", hjust=10) +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle=90, hjust=1))

colabels <- c('PRE-control', 'SIM-control', 'POST-control')

emeta.g + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 0.5, position=position_dodge(1)) +
  scale_fill_manual(values=c("orangered" ,'blue4' , 'yellow3')) + theme_classic() +
  theme(text = element_text(size=10), legend.position='none') + scale_x_discrete(labels= colabels)

posT1.3 <- subset(pT1.3, Metas > 0)
m.int <- glm(Metas ~ Treatment, data = posT1.3, family = 'poisson')
summary(m.int)
multzip <- emmeans(m.int, pairwise ~ Treatment, adjust = 'none')
multzip
tmeta.g <- ggplot(posT1.3, aes(x = Treatment, y = Metas, fill = Treatment)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Number of metacerciae in infected snails", hjust=10) +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle=90, hjust=1))

colabels <- c('PRE', 'SIM', 'POST')

tmeta.g + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 0.5, position=position_dodge(1)) +
  scale_fill_manual(values=c("orangered" ,'blue4' , 'yellow3')) + theme_classic() +
  theme(text = element_text(size=10), legend.position='none') + scale_x_discrete(labels= colabels)


#E vs T metas, all infection statuses
pET1.3 <- rbind(pE1, pT1)
posET1.3 <- subset(pET1.3, Metas > 0)
m1.prev <- glm(Metas ~ Treatment, data = posET1.3, family = 'poisson')
summary(m1.prev)

pET2.3 <- rbind(pE2, pT2)
posET2.3 <- subset(pET2.3, Metas > 0)
m2.prev <- glm(Metas ~ Treatment, data = posET2.3, family = 'poisson')
summary(m2.prev)

pET3.3 <- rbind(pE3, pT3)
posET3.3 <- subset(pET3.3, Metas > 0)
m3.prev <- glm(Metas ~ Treatment, data = posET3.3, family = 'poisson')
summary(m3.prev)

# E vs T for just coinfects
pTa1 <- subset(pT1.3, Total.fact == 'p')
posE1.3
newcolT <- rep('T', 32)
newcolE <- rep('E', 38)
pTaT <- cbind(pTa1, newcolT)
posEE <- cbind(posE1.3, newcolE)
names(posEE)[14] <- "newcolT"
allmet <- rbind(pTaT, posEE)
m4.prev <- glm(Metas ~ Treatment, data = allmet, family = 'poisson')
summary(m4.prev)
multzip <- emmeans(m4.prev, pairwise ~ Treatment, adjust = 'none')
multzip

m5.prev <- glm(Metas ~ newcolT, data = allmet, family = 'poisson')
summary(m5.prev)


#cercs shedding
hist(cercs1$Cercs)
cercs1 <- subset(cercs, Schisto.inf > 0)
Tco <- subset(cercs1, Total.fact == 'p')
scerc <- subset(cercs1, Treatment == 'S')
shed1 <- rbind(Tco, scerc)
shed <- subset(shed1, Week > 1)
fit_hzip <- glmmTMB(Cercs ~ Treatment + (1|ID), data=shed1, ziformula=~., family=poisson(link="log"))
summary(fit_hzip)
fit_zinbinom2 <- update(fit_hzip, family=nbinom2)
fit_zinbinom1 <- update(fit_hzip, family=nbinom1)
fit_hnbinom2 <- update(fit_zinbinom2, family=list(family="truncated_nbinom1",link="log"))
AICtab(fit_hzip, fit_zinbinom2, fit_zinbinom1, fit_hnbinom2)

#the best one
summary(fit_zinbinom1)

multzip <- emmeans(fit_zinbinom1, pairwise ~ Treatment)
test(multzip, adjust = 'none')

ggboxplot(shed, x = 'Week', y = 'Cercs', color = 'Treatment')
poshed <- subset(shed, Cercs > 0)
hist(poshed$Cercs)
library(lme4)
cercsnb <- glmer.nb(Cercs ~ Treatment + (1|ID), data=shed)
summary(cercsnb)
multzip <- emmeans(cercsnb, pairwise ~ Treatment, adjust = 'none')
multzip

cerc.g <- ggplot(shed1, aes(x = Treatment, y = Cercs, fill = Treatment)) + 
  geom_violin(trim = T, scale = 'width', draw_quantiles = c(0.25, 0.5, 0.75), adjust = 1.5) +
  labs(x="", y = "Number of cercariae per hour per snail", hjust=10) +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle=90, hjust=1))

colabels <- c('S', 'PRE', 'SIM', 'POST')

cerc.g + #geom_dotplot(binaxis='y', stackdir='centerwhole', stackratio = 0.5, dotsize = 0.5, position=position_dodge(1)) +
  scale_fill_manual(values=c("green4","orangered" ,'blue4' , 'yellow3')) + theme_classic() +
  theme(text = element_text(size=25), legend.position='none') + scale_x_discrete(labels= colabels)


##run as binomial
shed4 <- subset(shed1, Week == 1)
shed4$Cercs[shed4$Cercs > 0] <- 1 # replace all cerc values greater than 0 with 1
cerc.prev <- glm(Cercs ~ Treatment, data = shed4, family = 'binomial')
summary(cerc.prev)
multzip <- emmeans(cerc.prev, pairwise ~ Treatment, adjust = 'none')
multzip

cerc.prev <- glm(Cercs ~ Total.fact, data = shed4, family = 'binomial')
summary(cerc.prev)


library("survival")
library("survminer")

fit <- survfit(Surv(Time.from.start, Status) ~ Treatment, data = dub.whammy)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

surv_diff <- survdiff(Surv(Time.from.start, Status) ~ Treatment, data = dub.whammy)
surv_diff

cph <- coxph(Surv(Time.from.start, Status) ~ Treatment, data = dub.whammy)
cph
summary(cph)

psurv_diff <- pairwise_survdiff(Surv(Time.from.start, Status) ~ Treatment, data = dub.whammy, p.adjust.method = 'none')
psurv_diff


#from metas
survmetas <- rbind(pE1, pE2, pE3, pT1, pT2, pT3)
fit <- survfit(Surv(Time.from.meta.exposure, Status) ~ Treatment, data = survmetas)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

surv_diff <- survdiff(Surv(Time.from.meta.exposure, Status) ~ Treatment, data = survmetas)
surv_diff

cph <- coxph(Surv(Time.from.meta.exposure, Status) ~ Treatment, data = survmetas)
cph
summary(cph)

psurv_diff <- pairwise_survdiff(Surv(Time.from.meta.exposure, Status) ~ Treatment, data = survmetas, p.adjust.method = 'none')
psurv_diff


#from schisto
survschisto <- rbind(pS, pT1, pT2, pT3)
fit <- survfit(Surv(Time.from.schisto, Status) ~ Treatment, data = survschisto)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

surv_diff <- survdiff(Surv(Time.from.schisto, Status) ~ Treatment, data = survschisto)
surv_diff

cph <- coxph(Surv(Time.from.schisto, Status) ~ Treatment, data = survschisto)
cph
summary(cph)

psurv_diff <- pairwise_survdiff(Surv(Time.from.schisto, Status) ~ Treatment, data = survschisto, p.adjust.method = 'none')
psurv_diff

#groups
pC <- subset(dub.whammy, Treatment == 'C')
pT1 <- subset(dub.whammy, Treatment == 'T1')
pT2 <- subset(dub.whammy, Treatment == 'T2')
pT3 <- subset(dub.whammy, Treatment == 'T3')
pS <- subset(dub.whammy, Treatment == 'S')
pE1 <- subset(dub.whammy, Treatment == 'E1')
pE2 <- subset(dub.whammy, Treatment == 'E2')
pE3 <- subset(dub.whammy, Treatment == 'E3')


schisto <- rbind (pC, pS)
fit <- survfit(Surv(Time.from.start, Status) ~ Treatment, data = schisto)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

surv_diff <- survdiff(Surv(Time.from.start, Status) ~ Treatment, data = schisto)
surv_diff


cmeta <- rbind (pC, pE1, pE2, pE3)
fit <- survfit(Surv(Time.from.start, Status) ~ Treatment, data = cmeta)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

psurv_diff <- pairwise_survdiff(Surv(Time.from.start, Status) ~ Treatment, data = cmeta, p.adjust.method = 'none')
psurv_diff


ct <- rbind (pC, pT1, pT2, pT3)
fit <- survfit(Surv(Time.from.start, Status) ~ Treatment, data = ct)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = F, conf.int = TRUE,
           palette = c('#CC79A7','#FF9933','#0072B2' , '#F0E442'),
           legend.labs = c('CON', 'PRE', 'SIM', 'POST'),
           legend = c('right'),
           font.legend = 20,
           font.y = 20,
           font.tickslab = 20,
           font.x = 20,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

psurv_diff <- pairwise_survdiff(Surv(Time.from.start, Status) ~ Treatment, data = ct, p.adjust.method = 'none')
psurv_diff


t1 <- rbind (pE1, pT1)
fit <- survfit(Surv(Time.from.start, Status) ~ Treatment, data = t1)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

surv_diff <- survdiff(Surv(Time.from.start, Status) ~ Treatment, data = t1)
surv_diff


t2 <- rbind (pE2, pT2)
fit <- survfit(Surv(Time.from.start, Status) ~ Treatment, data = t2)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

surv_diff <- survdiff(Surv(Time.from.start, Status) ~ Treatment, data = t2)
surv_diff


t3 <- rbind (pE3, pT3)
fit <- survfit(Surv(Time.from.start, Status) ~ Treatment, data = t3)
print(fit)
summary(fit)

ggsurvplot(fit,
           pval = F, conf.int = TRUE,
           palette = c('#999999', '#F0E442'),
           legend.labs = c('C.POST', 'POST'),
           legend = c('right'),
           font.legend = 20,
           font.y = 20,
           font.tickslab = 20,
           font.x = 20,
           risk.table = F, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw() # Change ggplot2 theme
)

surv_diff <- survdiff(Surv(Time.from.start, Status) ~ Treatment, data = t3)
surv_diff










