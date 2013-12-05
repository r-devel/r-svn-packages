DATA LIST FILE= "datafile.dat"  free (",")
/ X1 X2 X3  .

VARIABLE LABELS
X1 "X1" 
 X2 "X2" 
 X3 "X3" 
 .

VALUE LABELS
/
X3  
1 "str_1" 
 2 "str_2" 
 3 "str_3" 
.

EXECUTE.
