# run12
library("imager")

###
### NaOH 0.28
###

## run12
plot(load.image(file="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.28/run12/00000001.jpg"))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.28/run12/00000001.jpg"), x<560, x>160, y<800, y>400))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.28/run12/00000230.jpg"), x<560, x>160, y<800, y>400))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.28/run12/00000360.jpg"), x<560, x>160, y<800, y>400))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.28/run12/00000992.jpg"), x<560, x>160, y<800, y>400))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.28/run12/00004000.jpg"), x<560, x>160, y<800, y>400))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.28/run12/00008000.jpg"), x<560, x>160, y<800, y>400))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.28/run12/00012000.jpg"), x<560, x>160, y<800, y>400))


file_list <- list.files(path="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.3/run12/", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<560, x>160, y<800, y>400)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run12 <- frames.mean
save(run12,file = "teste2_run12.RData")     

load(file = "teste2_run12.RData")
nframes = length(run12)
run12.df = matrix(unlist(run12), ncol = 1, nrow = nframes, byrow = TRUE)
plot(run12.df)

###--- fig artigo
frame.n <- seq(1,length(run12.df))
reaction.profile.df <- data.frame(run12.df,frame.n) 
library(ggplot2)
ggplot(data = reaction.profile.df, aes(x=frame.n, y=run12.df)) +
  geom_point() +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x= "frame number",
       y=expression(italic(bar(S)[ROI])))

###


#S
frames.i = which.max(run12.df) #frames.i=1400 gera b=-0.0035 (praticamente nao muda)
frames.f = which.min(run12.df[frames.i:length(run12.df)])
S.run12 = run12.df[frames.i:frames.f]
tempo = seq(0,(frames.f-frames.i))/30
plot(tempo,S.run12, xlab = "tempo /s", ylab = "mean(S)")

sensor1 = data.frame(S.run12,tempo)
c.0 <- min(S.run12) * 0.5
model.0 <- lm(log(S.run12 - c.0) ~ tempo, data=sensor1)
start <- list(a=exp(coef(model.0)[1]), b=coef(model.0)[2], c=c.0)
model <- nls(S.run12 ~ a * exp(b * tempo) + c, data = sensor1, start = start,
             control = nls.control(maxiter = 1000)) # b = 0.06834
model
#a        b        c 
#0.764262 -0.003175  0.268833
plot(tempo,S.run12)
p <- coef(model)
curve(p["a"] * exp(p["b"] * x) + p["c"], lwd=2, col="Red", add=TRUE)

# model.3 ln S = −kt + ln So
model.3 <- lm(log(S.run12) ~ tempo, data=sensor1)
model.3
# (Intercept)        tempo  
# -0.014737    -0.001992  
## usar dados normalizados 0 a 1 (similar abs)
S.run12.norm1 = (S.run12-S.run12[length(S.run12)])/(S.run12[1]-S.run12[length(S.run12)])
plot(tempo,S.run12.norm1)
# usando model.3
model.3 <- lm(log(S.run12.norm1[1:6000]) ~ tempo[1:6000])
# (Intercept)  tempo[1:6000]  
# 0.126876      -0.006185  

## S/So
model.4 <- lm(log(S.run12/S.run12[1]) ~ tempo)
model.4
# (Intercept)        tempo  
# -0.0002038   -0.0018805  


## run11
plot(load.image(file="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.28/run11/00000001.jpg"))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.28/run11/00005750.jpg"), x<560, x>160, y<760, y>360))

file_list <- list.files(path="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.3/run11/", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<560, x>160, y<760, y>360)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run11 <- frames.mean
save(run11,file = "teste2_run11.RData")  

rm(list=ls())

## run10
#plot(load.image(file="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.3/run10/00000001.jpg"))
#plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.3/run10/00000001.jpg"), x<560, x>160, y<740, y>340))

file_list <- list.files(path="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.3/run10/", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<560, x>160, y<760, y>360)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run10 <- frames.mean
save(run10,file = "teste2_run10.RData")  

##
## usando o modelo padrao ln (Aoo - At) ~ t
## para ln(Soo - St) ~ t
##

load("teste2_run10.RData")
S10 = matrix(unlist(run10), ncol = 1, nrow = length(run10), byrow = TRUE)
load("teste2_run11.RData")
S11 = matrix(unlist(run11), ncol = 1, nrow = length(run11), byrow = TRUE)
load("teste2_run12.RData")
S12 = matrix(unlist(run12), ncol = 1, nrow = length(run12), byrow = TRUE)

# intervalos de medidas
fps = 30 # depende do video
S10.r = S10[which.max(S10):length(S10)]
t10.r = seq(0,(length(S10)-which.max(S10)))/fps
S11.r = S11[which.max(S11):length(S11)]
t11.r = seq(0,(length(S11)-which.max(S11)))/fps
S12.r = S12[which.max(S12):length(S12)]
t12.r = seq(0,(length(S12)-which.max(S12)))/fps

process <- caret::preProcess(as.data.frame(S10.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S10.r))
S10.n <- norm_scale$S10.r
process <- caret::preProcess(as.data.frame(S11.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S11.r))
S11.n <- norm_scale$S11.r
process <- caret::preProcess(as.data.frame(S12.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S12.r))
S12.n <- norm_scale$S12.r

lnS10 = log( S10.n - S10.n[length(S10.n)] )
lnS11 = log( S11.n - S11.n[length(S11.n)] )
lnS12 = log( S12.n - S12.n[length(S12.n)] )

runln10.df <- data.frame(S=S10.n[1:length(lnS10)-1], 
                         LnS = lnS10[1:length(lnS10)-1],
                         time=t10.r[1:length(lnS10)-1],
                         run=rep("10",  length(lnS10)-1))
runln11.df <- data.frame(S=S11.n[1:length(lnS11)-1], 
                         LnS = lnS11[1:length(lnS11)-1],
                         time=t11.r[1:length(lnS11)-1],
                         run=rep("11",  length(lnS11)-1))
runln12.df <- data.frame(S=S12.n[1:length(lnS12)-1], 
                         LnS = lnS12[1:length(lnS12)-1],
                         time=t12.r[1:length(lnS12)-1],
                         run=rep("12",  length(lnS12)-1))
runln10_12.df <- rbind(runln10.df,runln11.df,runln12.df)

## --- figura artigo
# S x t
library(ggplot2)
ggplot(data = runln10_12.df, aes(x=time, y=S, color = run)) +
  geom_point() +
  labs(x= expression(italic(t[r] /s)), 
       y=expression(italic(bar(S)[n]))) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        legend.text=element_text(size=10))
# ---------------------------------  

## --- figura artigo
# ln(Soo - St) x t
library(latex2exp)
library(ggplot2)
ggplot(data = runln10.df, aes(x=time, y=S, color = run)) +
  geom_point() +
  labs(x = "COLOR",
       y = "REACTION TIME")+
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10)
  
## Graphical TOC
ggplot(data = runln10_12.df, aes(x=time, y=S, color = run)) +
  geom_point() +
  labs(x= expression(italic(t[r] /s)), 
       y=expression(italic(bar(S)[n]))) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        legend.text=element_text(size=10))
# ---------------------------------  


#--------------------------------
# de acordo com 
# file:///home/elcana/grupo/isabella/docs/M9201_FadingPhenolphthalein.pdf
# o ajuste para NaOH 0.3 M é 180 s

t_ajustado <- 180*30
runln10.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln10.df)
runln11.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln11.df)
runln12.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln12.df)

plot(runln10.df$time[1:5400], runln10.df$LnS[1:5400], 
     main = "NaOH 0.3", frame = FALSE,
     xlim = c(0,200), ylim=c(-2,0), col = "blue",
     xlab = "Time /s", ylab = "Ln(Soo - St)")
abline(lm(LnS[1:5400] ~ time[1:5400], data = runln10.df), col = "blue")
points(runln11.df$time[1:5400], runln11.df$LnS[1:5400], col = "red")
abline(lm(LnS[1:5400] ~ time[1:5400], data = runln11.df), col = "red")
points(runln12.df$time[1:5400], runln12.df$LnS[1:5400])
abline(lm(LnS[1:5400] ~ time[1:5400], data = runln12.df))
legend(150,-0.5,legend = c("run10","run11","run12"),
       col = c("blue", "red", "black"), lty = c(1,1,1), cex = 0.8)

k10_12 <- mean(runln10.fit$coefficients[2],
              runln11.fit$coefficients[2],
              runln12.fit$coefficients[2])
sd_k10_12 <- sd(c(runln10.fit$coefficients[2],
              runln11.fit$coefficients[2],
              runln12.fit$coefficients[2]))
conc.NaOH <- 0.28
slope_conc_10_12 <- k10_12/conc.NaOH # artigo 0.021


###
### NaOH 0.17
###

## run7
plot(load.image(file="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.2/run7/00000001.jpg"))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.2/run7/00000001.jpg"), x<560, x>160, y<800, y>400))

plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.2/run9/00000001.jpg"), x<550, x>150, y<770, y>370))

file_list <- list.files(path="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.2/run7/", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<560, x>160, y<800, y>400)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run7 <- frames.mean
save(run7,file = "teste2_run7.RData")  

rm(list=ls())
file_list <- list.files(path="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.2/run8/", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<560, x>160, y<800, y>400)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run8 <- frames.mean
save(run8,file = "teste2_run8.RData")    

rm(list=ls())
file_list <- list.files(path="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.2/run9/", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<550, x>150, y<770, y>370)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run9 <- frames.mean
save(run9,file = "teste2_run9.RData")    

load("teste2_run7.RData")
S7 = matrix(unlist(run7), ncol = 1, nrow = length(run7), byrow = TRUE)
load("teste2_run8.RData")
S8 = matrix(unlist(run8), ncol = 1, nrow = length(run8), byrow = TRUE)
load("teste2_run9.RData")
S9 = matrix(unlist(run9), ncol = 1, nrow = length(run9), byrow = TRUE)

fps = 30
S7.r = S7[which.max(S7):length(S7)]
t7.r = seq(0,(length(S7)-which.max(S7)))/fps
S8.r = S8[which.max(S8):length(S8)]
t8.r = seq(0,(length(S8)-which.max(S8)))/fps
S9.r = S9[which.max(S9):length(S9)]
t9.r = seq(0,(length(S9)-which.max(S9)))/fps

process <- caret::preProcess(as.data.frame(S7.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S7.r))
S7.n <- norm_scale$S7.r
process <- caret::preProcess(as.data.frame(S8.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S8.r))
S8.n <- norm_scale$S8.r
process <- caret::preProcess(as.data.frame(S9.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S9.r))
S9.n <- norm_scale$S9.r

lnS7 = log( S7.n - S7.n[length(S7.n)] )
lnS8 = log( S8.n - S8.n[length(S8.n)] )
lnS9 = log( S9.n - S9.n[length(S9.n)] )

runln7.df <- data.frame(S=S7.n[1:length(lnS7)-1], 
                        LnS = lnS7[1:length(lnS7)-1],
                        time=t7.r[1:length(lnS7)-1],
                        run=rep("7",  length(lnS7)-1))
runln8.df <- data.frame(S=S8.n[1:length(lnS8)-1], 
                        LnS = lnS8[1:length(lnS8)-1],
                        time=t8.r[1:length(lnS8)-1],
                        run=rep("8",  length(lnS8)-1))
runln9.df <- data.frame(S=S9.n[1:length(lnS9)-1], 
                        LnS = lnS9[1:length(lnS9)-1],
                        time=t9.r[1:length(lnS9)-1],
                        run=rep("9",  length(lnS9)-1))
runln7_9.df <- rbind(runln7.df,runln8.df,runln9.df)

t_ajustado <- 360*30
runln7.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln7.df)
runln8.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln8.df)
runln9.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln9.df)

k7_9 <- -mean(runln7.fit$coefficients[2],
              runln8.fit$coefficients[2],
              runln9.fit$coefficients[2])
k7_9
sd_k7_9 <- sd(c(-runln7.fit$coefficients[2],
                -runln8.fit$coefficients[2],
                -runln9.fit$coefficients[2]))
sd_k7_9

conc.NaOH <- 0.17
k7_9/conc.NaOH

## --- figura artigo
# S x t
library(ggplot2)
ggplot(data = runln7_9.df, aes(x=time, y=S, color = run)) +
  geom_point() +
  labs(x= expression(italic(t[r] /s)), 
       y=expression(italic(bar(S)[n]))) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        legend.text=element_text(size=10))
# ---------------------------------  
## --- figura artigo
# ln(Soo - St) x t
library(latex2exp)
library(ggplot2)
ggplot(data = runln7_9.df, aes(x=time, y=LnS, color = run)) +
  ylab(TeX("$ ln \\left( \\bar{S}_n - \\bar{S}_\\infty \\right)  $"))+
  xlab(TeX("$ t_r / s $"))+
  theme(axis.title=element_text(size=14,face="italic"),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        legend.text=element_text(size=10)) +
  xlim(0,360)+
  ylim(-1.2,0.3)+
  geom_smooth(data = runln7_9.df, aes(time,LnS), method=lm, se=FALSE, linetype = "dashed") 


#### Obs NaOH 0.17
#### usando a região entre 150 e 350 s (sugestão vídeo https://youtu.be/KjHCmjVtJ2Y)
runln7.fit <- lm(LnS[4500:10500] ~ time[4500:10500], data = runln7.df)
runln7.fit$coefficients[2]/conc.NaOH # 0.0217

####
####  NaOH 0.08
####

plot(imsub(load.image("/home/elcana/scr/NaOH 0.08/run4/00000001.jpg"), x<770, x>370, y<550, y>150))
plot(imsub(load.image("/home/elcana/scr/NaOH 0.08/run5/00000001.jpg"), x<790, x>390, y<570, y>170))
plot(imsub(load.image("/home/elcana/scr/NaOH 0.08/run6/00000001.jpg"), x<820, x>420, y<540, y>140))
file_list <- list.files(path="/home/elcana/scr/NaOH 0.08/run4", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i),  x<770, x>370, y<550, y>150)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run4 <- frames.mean
save(run4,file = "teste2_run4.RData")  

rm(list=ls())
file_list <- list.files(path="/home/elcana/scr/NaOH 0.08/run5", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i),  x<790, x>390, y<570, y>170)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run5 <- frames.mean
save(run5,file = "teste2_run5.RData")    

rm(list=ls())
file_list <- list.files(path="/home/elcana/scr/NaOH 0.08/run6", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<820, x>420, y<540, y>140)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run6 <- frames.mean
save(run6,file = "teste2_run6.RData") 


load("teste2_run4.RData")
S4 = matrix(unlist(run4), ncol = 1, nrow = length(run4), byrow = TRUE)
load("teste2_run5.RData")
S5 = matrix(unlist(run5), ncol = 1, nrow = length(run5), byrow = TRUE)
load("teste2_run6.RData")
S6 = matrix(unlist(run6), ncol = 1, nrow = length(run6), byrow = TRUE)

# modelos

fps = 30
S4.r = S4[which.max(S4):length(S4)]
t4.r = seq(0,(length(S4)-which.max(S4)))/fps
S5.r = S5[which.max(S5):length(S5)]
t5.r = seq(0,(length(S5)-which.max(S5)))/fps
S6.r = S6[which.max(S6):length(S6)]
t6.r = seq(0,(length(S6)-which.max(S6)))/fps

process <- caret::preProcess(as.data.frame(S4.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S4.r))
S4.n <- norm_scale$S4.r
process <- caret::preProcess(as.data.frame(S5.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S5.r))
S5.n <- norm_scale$S5.r
process <- caret::preProcess(as.data.frame(S6.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S6.r))
S6.n <- norm_scale$S6.r

lnS4 = log( S4.n - S4.n[length(S4.n)] )
lnS5 = log( S5.n - S5.n[length(S5.n)] )
lnS6 = log( S6.n - S6.n[length(S6.n)] )

runln4.df <- data.frame(S=S4.n[1:length(lnS4)-1], 
                        LnS = lnS4[1:length(lnS4)-1],
                        time=t4.r[1:length(lnS4)-1],
                        run=rep("4",  length(lnS4)-1))
runln5.df <- data.frame(S=S5.n[1:length(lnS5)-1], 
                        LnS = lnS5[1:length(lnS5)-1],
                        time=t5.r[1:length(lnS5)-1],
                        run=rep("5",  length(lnS5)-1))
runln6.df <- data.frame(S=S6.n[1:length(lnS6)-1], 
                        LnS = lnS6[1:length(lnS6)-1],
                        time=t6.r[1:length(lnS6)-1],
                        run=rep("6",  length(lnS6)-1))
runln4_6.df <- rbind(runln4.df,runln5.df,runln6.df)

t_ajustado <- 720*30
runln4.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln4.df)
runln5.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln5.df)
runln6.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln6.df)

k4_6 <- -mean(runln4.fit$coefficients[2],
              runln5.fit$coefficients[2],
              runln6.fit$coefficients[2])
k4_6
sd_k4_6 <- sd(c(-runln4.fit$coefficients[2],
                -runln5.fit$coefficients[2],
                -runln6.fit$coefficients[2]))
sd_k4_6

conc.NaOH <- 0.08
k4_6/conc.NaOH

## --- figura artigo
# S x t
library(ggplot2)
ggplot(data = runln4_6.df, aes(x=time, y=S, color = run)) +
  geom_point() +
  labs(x= expression(italic(t[r] /s)), 
       y=expression(italic(bar(S)[n]))) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        legend.text=element_text(size=10))
# ---------------------------------  
## --- figura artigo
# ln(Soo - St) x t
library(latex2exp)
library(ggplot2)
ggplot(data = runln4_6.df, aes(x=time, y=LnS, color = run)) +
  ylab(TeX("$ ln \\left( \\bar{S}_n - \\bar{S}_\\infty \\right)  $"))+
  xlab(TeX("$ t_r / s $"))+
  theme(axis.title=element_text(size=14,face="italic"),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        legend.text=element_text(size=10)) +
  geom_smooth(data = runln4_6.df, aes(time,LnS), method=lm, se=FALSE, linetype = "dashed") +
xlim(0,720)+
  ylim(-1.4,0.3)

#### Obs NaOH 0.08
#### usando a região entre 300 e 500 s (sugestão vídeo https://youtu.be/KjHCmjVtJ2Y)
runln4.fit <- lm(LnS[9000:15000] ~ time[9000:15000], data = runln4.df)
runln4.fit$coefficients[2]/conc.NaOH # 0.026

# avaliar S/So
S4.k = S4.k/S4.k[1]; S5.k = S5.k/S5.k[1]; S6.k = S6.k/S6.k[1] #(k/[] = 0.0232 nao muda com S/So )

# avaliar normalizar S entre 0 e 1 #(k/[] = 0.0232 nao muda com norm 0 1 ) [Normalize Data with Min-Max Scaling ]
process <- caret::preProcess(as.data.frame(S4.k), method=c("range"))
norm_scale <- predict(process, as.data.frame(S4.k))
S4.k <- norm_scale$S4.k
process <- caret::preProcess(as.data.frame(S5.k), method=c("range"))
norm_scale <- predict(process, as.data.frame(S5.k))
S5.k <- norm_scale$S5.k
process <- caret::preProcess(as.data.frame(S6.k), method=c("range"))
norm_scale <- predict(process, as.data.frame(S6.k))
S6.k <- norm_scale$S6.k

plot(runln4.df$time[1:t_ajustado], runln4.df$LnS[1:t_ajustado], 
     main = "NaOH 0.08", frame = FALSE,
     xlim = c(0,800), ylim=c(-1.5,0.5), col = "blue",
     xlab = "Time /s", ylab = "Ln(Soo - St)")
legend(600,0.3,legend = c("run4","run5","run6"),
       col = c("blue", "red", "black"), lty = c(1,1,1), cex = 0.8)


####
####  NaOH 0.045
####


## run4
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.045/run1/00000001.jpg"), x<800, x>400, y<540, y>140))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.045/run2/00000001.jpg"), x<800, x>400, y<540, y>140))
plot(imsub(load.image("/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.045/run3/00000001.jpg"), x<800, x>400, y<540, y>140))
file_list <- list.files(path="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.045/run1", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<800, x>400, y<540, y>140)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run1 <- frames.mean
save(run1,file = "teste2_run1.RData")  

rm(list=ls())
file_list <- list.files(path="/media/elcana/ELCANA1TB/Fenolftaleina/NaOH 0.045/run2", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<800, x>400, y<540, y>140)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run2 <- frames.mean
save(run2,file = "teste2_run2.RData")    

rm(list=ls())
file_list <- list.files(path="/home/elcana/scr/run3", full.names = TRUE)
frames.mean  = list()
for (i in file_list) {
  frame.x = imsub(load.image(i), x<800, x>400, y<540, y>140)
  s.mean = mean((RGBtoHSV(frame.x))[,,2])
  frames.mean[[i]] = c(s.mean)
}
run3 <- frames.mean
save(run3,file = "teste2_run3.RData") 

# modelos

load("teste2_run1.RData")
S1 = matrix(unlist(run1), ncol = 1, nrow = length(run1), byrow = TRUE)
load("teste2_run2.RData")
S2 = matrix(unlist(run2), ncol = 1, nrow = length(run2), byrow = TRUE)
load("teste2_run3.RData")
S3 = matrix(unlist(run3), ncol = 1, nrow = length(run3), byrow = TRUE)

fps = 30
S1.r = S1[which.max(S1):length(S1)]
t1.r = seq(0,(length(S1)-which.max(S1)))/fps
S2.r = S2[which.max(S2):length(S2)]
t2.r = seq(0,(length(S2)-which.max(S2)))/fps
S3.r = S3[which.max(S3):length(S3)]
t3.r = seq(0,(length(S3)-which.max(S3)))/fps

process <- caret::preProcess(as.data.frame(S1.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S1.r))
S1.n <- norm_scale$S1.r
process <- caret::preProcess(as.data.frame(S2.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S2.r))
S2.n <- norm_scale$S2.r
process <- caret::preProcess(as.data.frame(S3.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S3.r))
S3.n <- norm_scale$S3.r

lnS1 = log( S1.n - S1.n[length(S1.n)] )
lnS2 = log( S2.n - S2.n[length(S2.n)] )
lnS3 = log( S3.n - S3.n[length(S3.n)] )

runln1.df <- data.frame(S=S1.n[1:length(lnS1)-1], 
                        LnS = lnS1[1:length(lnS1)-1],
                        time=t1.r[1:length(lnS1)-1],
                        run=rep("1",  length(lnS1)-1))
runln2.df <- data.frame(S=S2.n[1:length(lnS2)-1], 
                        LnS = lnS2[1:length(lnS2)-1],
                        time=t2.r[1:length(lnS2)-1],
                        run=rep("2",  length(lnS2)-1))
runln3.df <- data.frame(S=S3.n[1:length(lnS3)-1], 
                        LnS = lnS3[1:length(lnS3)-1],
                        time=t3.r[1:length(lnS3)-1],
                        run=rep("3",  length(lnS3)-1))
runln1_3.df <- rbind(runln1.df,runln2.df,runln3.df)

t_ajustado <- 720*30
runln1.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln1.df)
runln2.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln2.df)
runln3.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln3.df)

k1_3 <- -mean(runln1.fit$coefficients[2],
              runln2.fit$coefficients[2],
              runln3.fit$coefficients[2])
k1_3
sd_k1_3 <- sd(c(-runln1.fit$coefficients[2],
                -runln2.fit$coefficients[2],
                -runln3.fit$coefficients[2]))
sd_k1_3

conc.NaOH <- 0.045
k1_3/conc.NaOH


## --- figura artigo
# S x t
library(ggplot2)
ggplot(data = runln1_3.df, aes(x=time, y=S, color = run)) +
  geom_point() +
  labs(x= expression(italic(t[r] /s)), 
       y=expression(italic(bar(S)[n]))) +
  theme(axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        legend.text=element_text(size=10))
# ---------------------------------  
## --- figura artigo
# ln(Soo - St) x t
library(latex2exp)
library(ggplot2)
ggplot(data = runln1_3.df, aes(x=time, y=LnS, color = run)) +
  ylab(TeX("$ ln \\left( \\bar{S}_n - \\bar{S}_\\infty \\right)  $"))+
  xlab(TeX("$ t_r / s $"))+
  theme(axis.title=element_text(size=14,face="italic"),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        legend.text=element_text(size=10)) +
  geom_smooth(data = runln1_3.df, aes(time,LnS), method=lm, se=FALSE, linetype = "dashed") +
  xlim(0,720)+
  ylim(-0.72,0.2)

#----------------------

#### Obs NaOH 0.045
#### usando a região entre 300 e 500 s (sugestão vídeo https://youtu.be/KjHCmjVtJ2Y)
runln1.fit <- lm(LnS[9000:15000] ~ time[9000:15000], data = runln1.df)
runln1.fit$coefficients[2]/conc.NaOH # 0.021


