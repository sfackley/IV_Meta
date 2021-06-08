rm(list=ls())
require('knitr')
require('tidyr')
require('dplyr')
require('gtools')
require('plotly')
require('gdata')
require('ggplot2')
require('magrittr')
require('reshape')
require('rlang')
require('pracma')
require('ggthemes')
require("binom")
require('drc')
z.crit <- qnorm(0.975)

# ADUCANUMAB Data
# ADAS COG 13

emerge.data <- data.frame(treatment=c("Placebo","Low","High"),
                       suvr=c(0.019,-0.165,-.272),
                       se.suvr=c(0.01,0.01,0.01),
                       cog=c(-3.3,-3.4,-2.7),
                       se.cog=c(0.25,0.25,0.25),
                       cog.total=c(288,293,298),
                       suvr.total=c(74,79,87),
                       a.cog=c(-5.162,-5.162+0.701,-5.162+1.4),
                       se.a.cog=c(0.3,0.3,0.3),
                       c.cog=c(-1.74,-1.74+0.26,-1.74+0.39),
                       se.c.cog=c(0.125,0.125,0.125)
)

# ADAS COG 13
engage.data <- data.frame(treatment=c("Placebo","Low","High"),
                          suvr=c(-0.005,-0.168,-.238),
                          se.suvr=c(0.01,0.01,0.01),
                          cog=c(-3.5,-3.3,-3.6),
                          se.cog=c(0.1,0.1,0.1),
                          cog.total=c(332,334,296),
                          suvr.total=c(104,116,97),
                          a.cog=c(-5.140,-5.140+0.583,-5.140+0.588),
                          se.a.cog=c(0.3,0.3,0.3),
                          c.cog=c(-1.56,-1.56+0.18,-1.56-0.03),
                          se.c.cog=c(0.125,0.125,0.125)
)

### BAN 2401 Data

dat <- data.frame(row.names= c("placebo","2.5.biweekly","5.monthly",
                               "5.biweekly",
                               "10.monthly","10.biweekly"))
dat$suvr.12 <- c(0.018,-0.061,-0.060,-0.128,-0.136,-0.217)
dat$suvr.18 <- c(0.029,-0.077,-0.061,-0.145,-0.156,-0.242)
dat$se.suvr.12 <- c(0.009,0.015,0.015,0.016,0.009,0.012)
dat$se.suvr.18 <- c(0.010,0.017,0.017,0.018,0.010,0.014)
apoe.e4.pos <- c(173,38,40,84,225,48)
apoe.e4.neg <- c(72,14,11,8,28,113)
dat$cog.12 <- -c(0.131,0.158,0.149,0.139,0.102,0.085)
dat$cog.18 <- -c(0.193,0.173,0.192,0.199,0.166,0.136)
dat$se.cog.12 <- c(0.013,0.027,0.027,0.021,0.014,0.017)
dat$se.cog.18 <- c(0.017,0.035,0.035,0.026,0.018,0.022)
dat$total <- apoe.e4.pos + apoe.e4.neg
dat$prop.e4.pos <- apoe.e4.pos/dat$total
dat %<>% mutate(se.suvr.12=se.suvr.12)
dat %<>% mutate(se.cog.12=se.cog.12)
dat %<>% mutate(se.suvr.18=se.suvr.18)
dat %<>% mutate(se.cog.18=se.cog.18)
dat$treatment <- c("Placebo","2.5mg Biweekly","5mg Monthly",
                   "5mg Biweekly",
                   "10mg Monthly","10mg Biweekly")

dat$suvr.total.12 <- c(96,27,27,25,88,43) 
dat$suvr.total.18 <- c(88,23,23,24,82,37)
dat$cog.total.12 <- c(187,38,42,67,165,93) 
dat$cog.total.18 <- c(160,33,35,61,146,79)

ban2401.dat <- dat

## Reformat Ban2401 Dat
select <- dplyr::select
extract <- tidyr::extract
rename <- dplyr::rename
ban.dat <- ban2401.dat %>% 
  select(c("treatment","total","prop.e4.pos"),everything())
ban.dat <- ban.dat %>% 
  gather(value,key,suvr.12:cog.total.18) %>% 
  extract("value",
          c("name","month"),
          "(.+).([:digit:]{2})") %>%
  spread(name,key) %>% select(-total)

# From https://jnnp.bmj.com/content/87/9/993.full
ban.covert.dat <- data.frame(
  ADCOMS.diff=c(0.027,0.057,0.082,0.002),
  MMSE.diff=c(0.51,1.58,1.24,0.12),
  num=c(sqrt(208^2+182^2),sqrt(71^2+74^2),sqrt(71^2+56^2),sqrt(208^2+208^2))
  #num=c(sqrt(182^2),sqrt(71^2),sqrt(56^2),sqrt(208^2))
  )
mod <- lm(MMSE.diff~ADCOMS.diff,weights=num,data=ban.covert.dat)
ban.dat.mmse <- ban.dat
ban.dat.mmse$cog <- predict(mod,data.frame(ADCOMS.diff=ban.dat.mmse$cog))
ban.dat.mmse$se.cog <- ban.dat.mmse$se.cog*mod$coefficients[[2]]

# ## ACC NCT01227564
# # Used MMSE bc was available at right time points. 
# 
# acc.data <- read.csv("ACC_data_NCT01227564.csv")
# acc.data$se.cog <- acc.data$half.ci.cog/qt(0.975,df=-1+acc.data$total.cog)
# acc.data$se.suvr <- acc.data$half.ci.suvr/qt(0.975,df=-1+acc.data$total.suvr)
# acc.data$cog.total <- acc.data$total.cog
# acc.data$suvr.total <-acc.data$total.suvr

## Bexarotene NCT01782742
# MMSE
# ADAS COG 0 (best) to 85 (worse)

bex.data <- data.frame(treatment=c("Bexarotene","Placebo"),
                       suvr=c(-0.028,0.023),
                       suvr.ci.length=c(0.064+0.008,0.049+0.096),
                       mmse.cog=c(0.750,1.750),
                       mmse.cog.ci.length=c(0.78+2.28,1.32+4.82),
                       adas.cog=c(-0.375,0.250),
                       adas.cog.ci.length=c(2.153+2.903,5.307+4.807),
                       total=c(16,4)
)
bex.data$se.suvr <- bex.data$suvr.ci.length/2/qt(0.975,-1+bex.data$total)
bex.data$se.cog <- bex.data$mmse.cog.ci.length/2/qt(0.975,-1+bex.data$total)
bex.data$cog <- bex.data$mmse.cog
bex.data$cog.total <- bex.data$total
bex.data$suvr.total <- bex.data$total
bex.data$a.cog <- bex.data$adas.cog
bex.data$se.a.cog <- bex.data$adas.cog.ci.length/2/qt(0.975,-1+bex.data$total)


## Solanezumab NCT00904683
# MMSE
# Expedition 1 & 2 
# Section 3.3 of Siemans 2016
# 6.21 (0.49)	4.08 (0.50)
# ADAS COG 11
sol.data <- data.frame(treatment=c("Solanezumab","Placebo"),
                       suvr=c(-0.01,0.02),
                       se.suvr=c(0.017,0.019),
                       cog=c(-1.83,-2.76),
                       se.cog=c(0.22,0.22),
                       cog.total=c(663,659),
                       suvr.total=c(97,98),
                       a.cog=c(-6.21,-4.08),
                       se.a.cog=c(0.49,0.50),
                       c.cog=c(-1.69,-1.53),
                       se.c.cog=c(0.13,0.13)
)
# Assumed about half in each group for SUVRs 
# C-check

## Solanezumab NCT01127633
# Expedition Ext
# MMSE
# ADAS COG 14
sol.data2 <- data.frame(treatment=c("Placebo","Solanezumab"),
                        suvr=c(0,-0.01),
                        se.suvr=c(0.131,0.222),
                        cog=c(-8.08,-7.60),
                        se.cog=c(0.270,0.263),
                        cog.total=c(412,446),
                        suvr.total=c(42,48),
                        a.cog=c(-17.58,-17.56),
                        se.a.cog=c(0.532,0.517),
                        adas.total=c(409,451),
                        c.cog=c(-5.59,-5.27),
                        se.c.cog=c(0.174,0.169),
                        c.total=c(416,452)
)
sol.data2$se.suvr <- sol.data2$se.suvr/sqrt(sol.data2$suvr.total)
# to convert from SD to SE
# C-check

## Solanezumab NCT01900665
# Expedition 3
#SUVr: Mean whole cerebellum corrected vs. subject specific white matter corrected. 
#MMSE
# ADAS COG 14

sol.data3 <- data.frame(treatment=c("Solanezumab","Placebo"),
                        suvr=c(-0.01,0.00),
                        se.suvr=c(0.005,0.005),
                        cog=c(-3.17,-3.66),
                        se.cog=c(0.154,0.156),
                        cog.total=c(893,876),
                        suvr.total=c(805,791),
                        a.cog=c(-6.65,-7.44),
                        se.a.cog=c(0.355,0.356),
                        adas.total=c(1053,1067),
                        c.cog=c(-1.91,-2.23),
                        se.c.cog=c(2.442,2.692),
                        c.total=c(903,895)
)
# C-check

## Conversion
dat <- read.csv(file="ADAS2MMSE.csv")
mod.adas.cog <- drm(MMSE ~ ADAS.Cog, data = dat, 
                    fct = l3u(upper=30), type = "continuous")
mod.cdr.sb <- drm(MMSE ~ CDR.SOB, data = dat, 
                  fct = l3u(upper=30), type = "continuous")
prediction.fun.adas.cog <- function(acog,model=mod.adas.cog){
  PR(mod.adas.cog,c(acog))
}
dat$pre.adas.cog <- prediction.fun.adas.cog(dat$ADAS.Cog)
prediction.fun.cdr.sb <- function(cdrcog,model=mod.cdr.sb){
  PR(mod.cdr.sb,c(cdrcog))
}
deriv.pred.fun.adas.cog <- function(acog,model=mod.adas.cog){
  numderiv(prediction.fun.adas.cog,acog,model=model)$df
}
deriv.pred.fun.cdr.sb <- function(acog,model=mod.cdr.sb){
  numderiv(prediction.fun.cdr.sb,acog,model=model)$df
}


## LY450139 NCT00762411

## "LY450139 dosing stopped due to evidence of dose-dependent cognitive/functional worsening."
## Used ADAS-Cog14. 
## Other measures were missing due to trial being stopped. 
## Higher score means worse cognitive function.

ly.data <- data.frame(treatment=c("LY450139","Placebo"),
                      suvr=c(-0.36,0.16),
                      se.suvr=c(0.23,0.11),
                      cog=c(-3.56,-3.35),
                      se.cog=c(0.38,0.35),
                      ac11.cog=c(7.37,6.77),
                      ac11.se.cog=c(0.79,0.72),
                      cog.total=c(555,553),
                      suvr.total=c(555,553),
                      starting=c(24.5,24.2),
                      a.cog=c(-9.23,-8.32),
                      se.a.cog=c(0.9,0.82),
                      c.cog=c(-3.05,-4.00),
                      se.c.cog=c(1.10,0.92),
                      c.total=c(22,34)
)

# treatment.start <- 24.5 
# placebo.start <- 24.2
# 
# treatment.change <- 7.37
# placebo.change <- 6.77
# 
# treatment.conversion =
#   deriv.pred.fun.adas.cog(treatment.change/2+treatment.start)
# placebo.conversion = 
#   deriv.pred.fun.adas.cog(placebo.change/2+placebo.start)
# 
# convert.vec <- c(treatment.conversion,placebo.conversion)
# 

ly.data.mmse <- ly.data

# ly.data.mmse$cog <- ly.data.mmse$cog*convert.vec
# ly.data.mmse$se.cog <- ly.data.mmse$se.cog*abs(convert.vec)
  

## LY450139 NCT00594568
## Used MMSE. 
## ADAS-COG 14

ly2.data <- data.frame(treatment=c("Placebo","LY450139 100mg","LY450139 140mg"),
                       suvr=c(0.08,0.06,0.09),
                       se.suvr=c(0.06,0.06,0.07),
                       cog=c(-2.95,-3.14,-3.71),
                       se.cog=c(0.26,0.27,0.27),
                       cog.total=c(376,289,274),
                       suvr.total=c(18,23,18),
                       a.cog=c(-7.42,-8.97,-9.48),
                       se.a.cog=c(0.63,0.67,0.69),
                       c.cog=c(-2.31,-2.73,-3.04),
                       se.c.cog=c(0.17,0.18,0.18),
                       c.total=c(267,216,195)
)

## Gantenerumab NCT01224106
# "The study was stopped early for futility, 
# but dose-dependent effects observed in 
# exploratory analyses on select clinical and 
# biomarker endpoints suggest that higher dosing 
# with gantenerumab may be necessary to achieve clinical efficacy."
# MMSE
# Ostrowitzki et al. (2017)
# 5.77 (4.54, 6.99) 	5.14 (3.91, 6.38) 	– 	5.54 (4.21, 6.87) 
# ADAS Cog 13

gan.data <- data.frame(treatment=c("Placebo","Gantenerumab 105mg","Gantenerumab 225mg"),
                       suvr=c(-0.02,-0.00,-0.09),
                       se.suvr=c(0.13,0.2,0.14),
                       cog=c(-2.93,-3.02,-2.73),
                       ci.length.cog=c(3.50-2.35,3.60-2.44,3.33-2.14),
                       cog.total=c(266,271,260),
                       suvr.total=c(21,15,19),
                       a.cog=c(-5.77,-5.14,-5.54),
                       ci.length.a.cog=c(-4.54+6.99,-3.91+6.38,-4.21+6.87),
                       c.cog=c(-1.60,-1.69,-1.73),
                       ci.length.c.cog=c(1.91-1.28,2.01-1.37,2.04-1.42)
)
#1.60 [1.28, 1.91], 1.69 [1.37, 2.01], and 1.73 [1.42, 2.04]
gan.data$se.cog <- gan.data$ci.length.cog/2/qt(0.975,-1+gan.data$cog.total)
gan.data$se.a.cog <- gan.data$ci.length.a.cog/2/qt(0.975,-1+gan.data$cog.total)
gan.data$se.c.cog <- gan.data$ci.length.c.cog/2/qt(0.975,-1+gan.data$cog.total)


## Bapineuzumab NCT00575055, NCT00676143, NCT00667810

# Note: only have change in cogntion relative to placebo group here (shouldn't matter for slope).
# MMSE
# –7·62 (–14·96 to –0·27)

bap.data1 <- data.frame(treatment=c("Bapineuzumab","Placebo"),
                       suvr=c(-0.09,0.20),
                       se.suvr=c(0.16,0.09),
                       cog=c(-3.02,0),
                       se.cog=c(7.41+1.38,0),
                       cog.total=c(19,7),
                       suvr.total=c(19,7),
                       a.cog=c(7.62,0),
                       se.a.cog=c(14.96-0.27,0)
                      )
bap.data1$se.cog <- bap.data1$se.cog/2/qt(0.975,sum(bap.data1$cog.total)-2)

# MMSE #ADAS COG 11
bap.data2 <- data.frame(treatment=c("Placebo","Bapineuzumab 0.5mg","Bapineuzumab 1mg"),
                        suvr_e4=c(0.102,0.001,NA),
                        se.suvr_e4=c(0.0264,0.0207,NA),
                        suvr_noe4=c(-0.046,0.039,-0.094),
                        se.suvr_noe4=c(0.0443,0.0452,0.0471),
                        cog_e4=c(-4.5,-4.7,NA),
                        se.cog_e4=c(0.2,0.2,NA),
                        cog_noe4=c(-3.9,-3.5,-3.7),
                        se.cog_noe4=c(0.2,0.3,0.3),
                        cog.total_e4=c(432,658,NA),
                        cog.total_noe4=c(493,314,307),
                        suvr.total_e4=c(40,75,NA),
                        suvr.total_noe4=c(15,12,12),
                        a.cog_e4=c(-8.7,-8.5,NA),
                        se.a.cog_e4=c(0.5,0.5,NA),
                        a.cog_noe4=c(-7.4,-7.1,-7.8),
                        se.a.cog_noe4=c(0.5,0.6,0.6),
                        c.cog_e4=c(-3.0,-3.3,NA),
                        se.c.cog_e4=c(0.2,0.1,NA),
                        c.cog_noe4=c(-2.6,-2.6,-2.8),
                        se.c.cog_noe4=c(0.2,0.2,0.2)
)
#c-check

bap.reg.data <- bap.data2 %>% 
  gather(key=att,value=meas,-treatment) %>% 
  separate(att,into=c("var","E4"),sep="_") %>%
  spread(var,meas)
bap.reg.data %<>% drop_na


## Verubecestat (MK-8931) NCT01739348 
# MMSE #ADAS COG11
ver.data <- data.frame(treatment=c("Verubecestat 12mg","Verubecestat 40mg","Placebo"),
                       suvr=c(-0.02,-0.04,0.00),
                       conf.int.len.suvr=c(0.04-0.01,0.06-0.03,0.01+0.01),
                       cog=c(-3.9,-3.6,-4.1),
                       conf.int.len.cog=c(4.3-3.6,4.0-3.3,4.5-3.8),
                       cog.total=c(610,600,628),
                       suvr.total=c(20,10,14),
                       a.cog.total=c(631,626,644),
                       a.cog=c(-7.9,-8.0,-7.7),
                       conf.int.len.a.cog=c(8.6-7.2,8.7-7.3,8.4-7.0),
                       c.cog=c(-2.1,-2.1,-2.1),
                       conf.int.len.c.cog=c(2.3-1.8,2.4-1.9,2.3-1.9),
                       c.cog.total=c(611,600,623)
)
ver.data$se.cog <- ver.data$conf.int.len.cog/2/qt(0.975,df=-1+ver.data$cog.total)
ver.data$se.suvr <-ver.data$conf.int.len.suvr/2/qt(0.975,df=-1+ver.data$suvr.total)
ver.data$se.a.cog <- ver.data$conf.int.len.a.cog/2/qt(0.975,df=-1+ver.data$a.cog.total)
ver.data$se.c.cog <- ver.data$conf.int.len.c.cog/2/qt(0.975,df=-1+ver.data$c.cog.total)


## Verubecestat (MK-8931) NCT01953601 
# CDR-SB
ver.data2 <- data.frame(treatment=c("Verubecestat 12 mg",
                                    "Verubecestat 40 mg",
                                    "placebo"),
                       suvr=c(-0.03,-0.04,0.02),
                       conf.int.len.suvr=c(0.04-0.03,0.05-0.04,0.03-0.02),
                       cog=c(1.6,2.0,1.6),
                       conf.int.len.cog=c(1.9-1.4,2.3-1.8,1.8-1.3),
                       cog.total=c(465,458,469),
                       suvr.total=c(63,59,65)
)
ver.data2$se.cog <- ver.data$conf.int.len.cog/2/qt(0.975,df=-1+ver.data2$cog.total)
ver.data2$se.suvr <-ver.data$conf.int.len.suvr/2/qt(0.975,df=-1+ver.data2$suvr.total)

treatment.start.1 <- 2.7  
treatment.start.2 <- 2.7
placebo.start <- 2.6

treatment.change.1 <- 1.6
treatment.change.2 <- 2.0
placebo.change <- 1.6

treatment.conversion.1 =
  deriv.pred.fun.cdr.sb(treatment.change.1/2+treatment.start.1)
treatment.conversion.2 =
  deriv.pred.fun.cdr.sb(treatment.change.2/2+treatment.start.2)
placebo.conversion = 
  deriv.pred.fun.cdr.sb(placebo.change/2+placebo.start)

convert.vec <- c(treatment.conversion.1,
                 treatment.conversion.2,
                 placebo.conversion)

ver.data2.mmse <- ver.data2
ver.data2.mmse$cog <- ver.data2$cog*convert.vec
ver.data2.mmse$se.cog <- ver.data2$se.cog*abs(convert.vec)

###
trial.table <- data.frame(
  clin.trial.num = c(
    "NCT01782742",
    "NCT00904683",
    "NCT01127633",
    "NCT02760602",
    "NCT00762411",
    "NCT00594568",
    "NCT01224106",
    "NCT00575055",
    "NCT00667810 & NCT00676143",
    "NCT01739348",
    "NCT01953601")
)
trial.table$number.arms <- c(
  2,
  2,
  2,
  2,
  2,
  3,
  3,
  2,
  3,
  3,
  3
)
trial.table$number.cog <- c(
  20,
  1222,
  412+446,
  893+876,
  555+553,
  270+216+199,
  266+271+260,
  19+7,
  432+658+493+314+307,
  610+600+628,
  465+458+469
)
trial.table$number.pet <- c(
  20,
  251,
  90,
  805+791,
  555+553,
  40+42+43,
  40+42+43,
  19+7,
  45+75+15+12+12,
  20+10+14,
  63+59+65 
)
trial.table$drug <- c(
  "Bexarotene",
  "Solanezumab",
  "Solanezumab",
  "Solanezumab",
  "LY450139",
  "LY450139",
  "Gantenerumab",
  "Bapineuzumab",
  "Bapineuzumab",
  "Verubecestat",
  "Verubecestat"
)
trial.table$population <- c(
  "Mild to Moderate Alzheimer's Disease",
  "Mild to Moderate Alzheimer's Disease",
  "Alzheimer's Disease",
  "Prodromal Alzheimer's Disease",
  "Mild to Moderate Alzheimer's Disease",
  "Mild to Moderate Alzheimer's Disease",
  "Prodromal Alzheimer's Disease",
  "Mild to Moderate Alzheimer's Disease",
  "Mild to Moderate Alzheimer's Disease in APOE-E4 carriers and non-carriers",
  "Mild to Moderate Alzheimer's Disease",
  "Prodromal Alzheimer's Disease"
)
colnames(trial.table) <- c(
  "Clinical Trial Number(s)",
  "Number of Treatment Arms",
  "Number of Participants with Cognitive Assessment",
  "Number of Participants with PET",
  "Drug Name",
  "Trial Population"
)

save.image("data.RData")