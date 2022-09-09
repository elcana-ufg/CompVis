# load RData file a calculate the rate of reaction

##
## usando o modelo padrao ln (Aoo - At) ~ t
## para ln(Soo - St) ~ t
##

load("teste2_run10.RData")
S10 = matrix(unlist(run10), ncol = 1, nrow = length(run10), byrow = TRUE)


# intervalos de medidas
fps = 30 # depende do video
S10.r = S10[which.max(S10):length(S10)]
t10.r = seq(0,(length(S10)-which.max(S10)))/fps

process <- caret::preProcess(as.data.frame(S10.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S10.r))
S10.n <- norm_scale$S10.r

lnS10 = log( S10.n - S10.n[length(S10.n)] )

#--------------------------------
# de acordo com 
# file:///home/elcana/grupo/isabella/docs/M9201_FadingPhenolphthalein.pdf
# o ajuste para NaOH 0.3 M é 180 s

t_ajustado <- 180*30
runln10.fit <- lm(LnS[1:t_ajustado] ~ time[1:t_ajustado], data = runln10.df)

k <- runln10.fit$coefficients[2]
              
conc.NaOH <- 0.28
slope_conc_10_12 <- k10_12/conc.NaOH # artigo 0.021



