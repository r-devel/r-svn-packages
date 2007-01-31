### $Id: Lipo.R,v 1.1 1999/10/21 22:55:32 bates Exp $
### Lipoprotein concentrations versus time
Lipo <-
  data.frame(time = c(0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10),
                 conc = c(46.1, 25.9, 17, 12.1, 7.22, 4.51, 3.19, 2.4, 1.82, 
                   1.41, 1, 0.94))
attr(Lipo, "reference") <- "A1.16, p. 282"
     
