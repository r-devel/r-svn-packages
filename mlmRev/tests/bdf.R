library(mlmRev)
(fm5 <- lme(langPOST ~ IQ.ver.cen + avg.IQ.ver.cen,
                       data = bdf, random = ~ IQ.ver.cen | schoolNR))
fixef(fm5)
ranef(fm5)
(fm1 <- lmer(langPOST ~ IQ.ver.cen+avg.IQ.ver.cen+(IQ.ver.cen|schoolNR),
             bdf))
q("no")

