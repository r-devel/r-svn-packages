### $Id: Rumford.R,v 1.1 1999/10/21 22:55:32 bates Exp $
Rumford <-
  data.frame(time = c(4, 5, 7, 12, 14, 16, 20, 24, 28, 31, 34, 37.5, 41),
             temp = c(126, 125, 123, 120, 119, 118, 116, 115, 114, 113, 112, 111, 110))
attr(Rumford, "reference") <- "A1.2, p. 268"
