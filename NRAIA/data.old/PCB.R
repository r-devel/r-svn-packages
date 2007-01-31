### $Id: PCB.R,v 1.1 1999/10/21 22:55:32 bates Exp $
### Concentration of PCB's in fish that were tagged by year of birth
"PCB" <-
  do.call("data.frame", list(age = as.integer(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7,
                   7, 7, 8, 8, 8, 9, 11, 12, 12, 12)),
                 conc = c(0.6, 1.6, 0.5, 1.2, 2, 1.3, 2.5, 2.2, 2.4, 1.2, 3.5,
                   4.1, 5.1, 5.7, 3.4, 9.7, 8.6, 4, 5.5, 10.5, 17.5, 13.4, 4.5,
                   30.4, 12.4, 13.4, 26.2, 7.4)))
attr(PCB, "reference")  = "A1.1, p. 267"

