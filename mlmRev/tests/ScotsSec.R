library(mlmRev)
options(show.signif.stars = FALSE)
(fm1 <- lmer(attain ~ verbal * sex + (1|primary) + (1|second), ScotsSec,
             control = list(EMv = 1, msV = 1)))
(fm2 <- lmer(attain ~ verbal * sex + (1|primary) + (1|second), ScotsSec,
             control = list(EMv = 1, msV = 1, niterEM = 40)))
q("no")
