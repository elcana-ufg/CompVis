# load RData file and calculate the rate of reaction, k

load("Smean.RData")
S = matrix(unlist(run), ncol = 1, nrow = length(run), byrow = TRUE)

fps = 30 # video frame per second
# S mean values during the time of reaction
S.r = S[which.max(S):length(S)] 
t.r = seq(0,(length(S)-which.max(S)))/fps
# S min-max normalization
process <- caret::preProcess(as.data.frame(S.r), method=c("range"))
norm_scale <- predict(process, as.data.frame(S.r))
S.n <- norm_scale$S.r
# Ln(S)
lnS = log( S.n - S.n[length(S.n)] )
t_adj <- 180*30 # time of 180 s
runln.fit <- lm(LnS[1:t_adj] ~ time[1:t_adj])

k <- runln.fit$coefficients[2]
              
cont.NaOH <- 0.28
slope_cont <- k/cont.NaOH 



