��    c      4  �   L      p  �   q  7   �     1	     M	  '   e	  3   �	      �	  2   �	  L   
  #   b
  9   �
  1   �
  ?   �
  5   2  ?   h  $   �  M   �  &        B  )   S     }  C   �  K   �  @   ,  $   m  <   �  )   �  ;   �  !   5      W  8   x  J   �  -   �  �   *  4   �  V     >   t  <   �     �  L   p  2   �  .   �  +     ;   K  ?   �  T   �  G     E   d  K   �  2   �  >   )  v   h  /   �  ,     ;   <  1   x     �  4   �  8   �  7   4  0   l  8   �  -   �  5     )   :  M   d  6   �  0   �          9  $   M     r     �     �     �     �     �  %     B   .  "   q  '   �  7   �      �       /   "  1   R  \   �  5   �  F     #   ^  2   �  :   �     �       )     '   D     l  !   }    �  �   �  7   <     t     �  '   �  7   �  $     6   -  P   d  '   �  A   �  5      C   U   9   �   C   �   (   !  Q   @!  *   �!     �!  -   �!  "    "  G   #"  O   k"  D   �"  (    #  @   )#  -   j#  ?   �#  %   �#  $   �#  8   #$  J   \$  1   �$  �   �$  4   �%  V   �%  >   +&  @   j&     �&  L   +'  2   x'  .   �'  /   �'  C   
(  ?   N(  T   �(  G   �(  E   +)  K   q)  2   �)  B   �)  v   3*  7   �*  ,   �*  ;   +  1   K+     }+  4   �+  8   �+  7   ,  0   ?,  @   p,  -   �,  5   �,  )   -  Q   ?-  6   �-  0   �-     �-     .  (   0.     Y.     y.     �.     �.     �.     �.  %   �.  B   !/  "   d/  +   �/  7   �/      �/     0  /   0  1   I0  l   {0  A   �0  J   *1  #   u1  :   �1  :   �1     2     "2  )   92  '   c2     �2  !   �2     
   L       ,   K       '             _   &   D   c       -       T   [   B               1               F              W            "      @   2   4   =   ^       M                          b          5           6   ;       >      0   9   G       !      ]   Q       `       R   I   E   ?   	       7      C       *   Y      :   .   V   #   S              a   Z      O   U       )                 %   A      X   N          8   3   \              P       <       /              +         H      $       (   J        %d observation (%s) has *only* NAs --> omit them for clustering! %d observations (%s ...) have *only* NAs --> omit them for clustering! %s has constant columns %s; these are standardized to 0 %s has invalid column names %s must be in 1:ncol(x) %s must contain column names or numbers 'A' must be p x p  cov-matrix defining an ellipsoid 'B' has to be a positive integer 'col.clus' should have length 4 when color is TRUE 'diss = TRUE' is not implemented yet.  Please write to maintainer("cluster") 'diss' must be one of {TRUE, FALSE} 'dmatrix' is not a dissimilarity matrix compatible to 'x' 'full' must be FALSE, TRUE, or a number in [0, 1] 'iniMem.p' must be a nonnegative n * k matrix with rowSums == 1 'k' (number of clusters) must be in {1,2, .., n/2 -1} 'm', a membership matrix, must be nonnegative with rowSums == 1 'maxit' must be non-negative integer 'medoids' must be NULL or vector of %d distinct indices in {1,2, .., n}, n=%d 'memb.exp' must be a finite number > 1 'n' must be >= 2 'par.method' must be of length 1, 3, or 4 'samples' should be at least 1 'sampsize' = %d should not be larger than the number of objects, %d 'sampsize' must be a positive integer, not larger than .Machine$integer.max 'sampsize' should be at least %d = max(2, 1+ number of clusters) 'weights' must be of length p (or 1) 'x' is a "dist" object, but should be a data matrix or frame 'x' is not a numeric dataframe or matrix. 'x' is not and cannot be converted to class "dissimilarity" 'x' must be numeric  n x p matrix 'x' must only have integer codes >>>>> funny case in clusplot.default() -- please report! All variables must be binary (e.g., a factor with 2 levels, both present). Cannot keep data when 'x' is a dissimilarity! Distance computations with NAs: using correct instead of pre-2016 wrong formula.
Use  'correct.d=FALSE'  to get previous results or set 'correct.d=TRUE' explicitly
to suppress this warning. Distances must be result of dist or a square matrix. Each of the random samples contains objects between which no distance can be computed. Error in C routine for the spanning ellipsoid,
 rank problem?? FANNY algorithm has not converged in 'maxit' = %d iterations For each of the %d samples, at least one object was found which could not be assigned to a cluster (because of missing values). Missing values were displaced by the median of the corresponding variable(s) NA values in the dissimilarity matrix not allowed. NA-values are not allowed in clustering vector NA-values are not allowed in dist-like 'x'. Need either a dissimilarity 'dist' or diss.matrix 'dmatrix' No clustering performed, NA values in the dissimilarity matrix. No clustering performed, a variable was found with all non missing values identical. No clustering performed, all variables have at least one missing value. No clustering performed, an object was found with all values missing. No clustering performed, found variable with more than half values missing. No valid silhouette information (#{clusters} =? 1) Number of clusters 'k' must be in {1,2, .., n-1}; hence n >= 2 Observation %s has *only* NAs --> omit it for clustering Observations %s have *only* NAs --> omit them for clustering! Set either 'variant' or 'pamonce', but not both The clustering vector is of incorrect length The number of cluster should be at least 1 and at most n-1. algorithm possibly not converged in %d iterations ambiguous clustering method at least one binary variable has more than 2 levels. at least one binary variable has not 2 different levels. at least one binary variable has values not in {0,1,NA} binary variable(s) %s treated as interval scaled clustering 'x' and dissimilarity 'dist' are incompatible computed some negative or all 0 probabilities ellipsoidPoints() not yet implemented for p >= 3 dim. error from .C(cl_pam, *): invalid medID's full silhouette is only available for results of 'clara(*, keep.data = TRUE)' have %d observations, but not more than %d are allowed index has to be a function or a list of function invalid %s; must be named list invalid 'correct.d' invalid 'jstop' from .C(cl_clara,.): invalid 'silhouette' object invalid 'spaceH0': invalid 'twins' object invalid clustering method invalid partition object invalid silhouette structure invalid type %s for column numbers %s mona() needs at least p >= 2 variables (in current implementation) need at least 2 objects to cluster no diss nor data found for 'clusplot()' no diss nor data found, nor the original argument of %s no points without missing values omitting NAs one or more objects contain only missing values one or more variables contain only missing values setting 'logical' variable %s to type 'asymm' setting 'logical' variables %s to type 'asymm' specified both 'full' and 'subset'; will use 'subset' the memberships are all very close to 1/k. Maybe decrease 'memb.exp' ? the square matrix is not symmetric. when 'medoids.x' is FALSE, 'keep.data' must be too with mixed variables, metric "gower" is used automatically x has zero columns x is not a data matrix x is not a dataframe or a numeric matrix. x is not a numeric dataframe or matrix. x is not numeric x must be a matrix or data frame. Project-Id-Version: cluster 2.1.9
PO-Revision-Date: 2025-01-14 15:25
Last-Translator: Automatically generated
Language-Team: none
Language: en
MIME-Version: 1.0
Content-Type: text/plain; charset=UTF-8
Content-Transfer-Encoding: 8bit
Plural-Forms: nplurals=2; plural=(n != 1);
 %d observation (%s) has *only* NAs --> omit them for clustering! %d observations (%s ...) have *only* NAs --> omit them for clustering! %s has constant columns %s; these are standardized to 0 %s has invalid column names %s must be in 1:ncol(x) %s must contain column names or numbers ‘A’ must be p x p  cov-matrix defining an ellipsoid ‘B’ has to be a positive integer ‘col.clus’ should have length 4 when color is TRUE ‘diss = TRUE’ is not implemented yet.  Please write to maintainer("cluster") ‘diss’ must be one of {TRUE, FALSE} ‘dmatrix’ is not a dissimilarity matrix compatible to ‘x’ ‘full’ must be FALSE, TRUE, or a number in [0, 1] ‘iniMem.p’ must be a nonnegative n * k matrix with rowSums == 1 ‘k’ (number of clusters) must be in {1,2, .., n/2 -1} ‘m’, a membership matrix, must be nonnegative with rowSums == 1 ‘maxit’ must be non-negative integer ‘medoids’ must be NULL or vector of %d distinct indices in {1,2, .., n}, n=%d ‘memb.exp’ must be a finite number > 1 ‘n’ must be >= 2 ‘par.method’ must be of length 1, 3, or 4 ‘samples’ should be at least 1 ‘sampsize’ = %d should not be larger than the number of objects, %d ‘sampsize’ must be a positive integer, not larger than .Machine$integer.max ‘sampsize’ should be at least %d = max(2, 1+ number of clusters) ‘weights’ must be of length p (or 1) ‘x’ is a "dist" object, but should be a data matrix or frame ‘x’ is not a numeric dataframe or matrix. ‘x’ is not and cannot be converted to class "dissimilarity" ‘x’ must be numeric  n x p matrix ‘x’ must only have integer codes >>>>> funny case in clusplot.default() -- please report! All variables must be binary (e.g., a factor with 2 levels, both present). Cannot keep data when ‘x’ is a dissimilarity! Distance computations with NAs: using correct instead of pre-2016 wrong formula.
Use  ‘correct.d=FALSE’  to get previous results or set ‘correct.d=TRUE’ explicitly
to suppress this warning. Distances must be result of dist or a square matrix. Each of the random samples contains objects between which no distance can be computed. Error in C routine for the spanning ellipsoid,
 rank problem?? FANNY algorithm has not converged in ‘maxit’ = %d iterations For each of the %d samples, at least one object was found which could not be assigned to a cluster (because of missing values). Missing values were displaced by the median of the corresponding variable(s) NA values in the dissimilarity matrix not allowed. NA-values are not allowed in clustering vector NA-values are not allowed in dist-like ‘x’. Need either a dissimilarity ‘dist’ or diss.matrix ‘dmatrix’ No clustering performed, NA values in the dissimilarity matrix. No clustering performed, a variable was found with all non missing values identical. No clustering performed, all variables have at least one missing value. No clustering performed, an object was found with all values missing. No clustering performed, found variable with more than half values missing. No valid silhouette information (#{clusters} =? 1) Number of clusters ‘k’ must be in {1,2, .., n-1}; hence n >= 2 Observation %s has *only* NAs --> omit it for clustering Observations %s have *only* NAs --> omit them for clustering! Set either ‘variant’ or ‘pamonce’, but not both The clustering vector is of incorrect length The number of cluster should be at least 1 and at most n-1. algorithm possibly not converged in %d iterations ambiguous clustering method at least one binary variable has more than 2 levels. at least one binary variable has not 2 different levels. at least one binary variable has values not in {0,1,NA} binary variable(s) %s treated as interval scaled clustering ‘x’ and dissimilarity ‘dist’ are incompatible computed some negative or all 0 probabilities ellipsoidPoints() not yet implemented for p >= 3 dim. error from .C(cl_pam, *): invalid medID's full silhouette is only available for results of ‘clara(*, keep.data = TRUE)’ have %d observations, but not more than %d are allowed index has to be a function or a list of function invalid %s; must be named list invalid ‘correct.d’ invalid ‘jstop’ from .C(cl_clara,.): invalid ‘silhouette’ object invalid ‘spaceH0’: invalid ‘twins’ object invalid clustering method invalid partition object invalid silhouette structure invalid type %s for column numbers %s mona() needs at least p >= 2 variables (in current implementation) need at least 2 objects to cluster no diss nor data found for ‘clusplot()’ no diss nor data found, nor the original argument of %s no points without missing values omitting NAs one or more objects contain only missing values one or more variables contain only missing values setting ‘logical’ variable %s to type ‘asymm’ setting ‘logical’ variables %s to type ‘asymm’ specified both ‘full’ and ‘subset’; will use ‘subset’ the memberships are all very close to 1/k. Maybe decrease ‘memb.exp’ ? the square matrix is not symmetric. when ‘medoids.x’ is FALSE, ‘keep.data’ must be too with mixed variables, metric "gower" is used automatically x has zero columns x is not a data matrix x is not a dataframe or a numeric matrix. x is not a numeric dataframe or matrix. x is not numeric x must be a matrix or data frame. 