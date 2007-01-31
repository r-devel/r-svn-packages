### $Id: Saccharin.R,v 1.1 1999/10/21 22:55:32 bates Exp $
Saccharin <-
  data.frame(start = c(0, 5, 15, 30, 45, 60, 75, 90, 105), 
             length = c(5, 10, 15, 15, 15, 15, 15, 15, 15),
             sacch = c(7518, 6275, 4989, 2580, 1485, 861, 561, 363, 300))
attr(Saccharin, "reference") <- "A1.11, p. 277"

