library(MEMSS)
options(show.signif.stars = FALSE)
lmer(pixel ~ day + I(day^2) + (1|DS) + (day|Dog), Pixel)
q("no")
