library(foreign)
read.S("mySobj")
data.restore("dumpdata", print = TRUE)
print(myobj)
data.restore("tsdumpdata", print = TRUE)
sunspot
utils::str(carbon.dioxide)
carbon.dioxide[1:60]
q()
