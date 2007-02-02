### $Id: BOD2.R,v 1.1 1999/10/21 22:55:32 bates Exp $
### second Biochemical Oxygen Demand (BOD) data set from Marske
"BOD2" <-
  data.frame(Time = c(1, 2, 3, 4, 5, 7, 9, 11),
             demand = c(0.47, 0.74, 1.17, 1.42, 1.6, 1.84, 2.19, 2.17))
attr(BOD2, "reference") <- "A4.1, p. 305"
