library(mlmRev)
options(show.signif.stars = FALSE)
Early$tos <- Early$age - 0.5
(fm1 <- lmer(cog ~ tos * trt + (tos|id), Early, con = list(EMv = 1, msV = 1)))
q("no")
