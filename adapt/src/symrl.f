      subroutine symrl(s, center, hwidth, minord, maxord, intvls,
     *     intcls, numsms, weghts, fulsms, fail)
c  multidimensional fully symmetric rule integration subroutine
c
c   this subroutine computes a sequence of fully symmetric rule
c   approximations to a fully symmetric multiple integral.
c   written by a. genz, mathematical institute, university of kent,
c   canterbury, kent ct2 7nf, england
c
c**************	 parameters for symrl  ********************************
c*****input parameters
c  s	   integer number of variables, must exceed 0 but not exceed 20
c  f	   externally declared user defined real function integrand.
c	   it must have parameters (s,x), where x is a real array
c	   with dimension s.
c  minord  integer minimum order parameter.  on entry minord specifies
c	   the current highest order approximation to the integral,
c	   available in the array intvls.  for the first call of symrl
c	   minord should be set to 0.  otherwise a previous call is
c	   assumed that computed intvls(1), ... , intvls(minord).
c	   on exit minord is set to maxord.
c  maxord  integer maximum order parameter, must be greater than minord
c	   and not exceed 20. the subroutine computes intvls(minord+1),
c	   intvls(minord+2),..., intvls(maxord).
c  g	   real array of dimension(maxord) of generators.
c	   all generators must be distinct and nonnegative.
c  numsms  integer length of array fulsms, must be at least the sum of
c	   the number of distinct partitions of length at most s
c	   of the integers 0,1,...,maxord-1.  an upper bound for numsms
c	   when s+maxord is less than 19 is 200
c******output parameters
c  intvls  real array of dimension(maxord).  upon successful exit
c	   intvls(1), intvls(2),..., intvls(maxord) are approximations
c	   to the integral.  intvls(d+1) will be an approximation of
c	   polynomial degree 2d+1.
c  intcls  integer total number of f values needed for intvls(maxord)
c  weghts  real working storage array with dimension (numsms). on exit
c	   weghts(j) contains the weight for fulsms(j).
c  fulsms  real working storage array with dimension (numsms). on exit
c	   fulsms(j) contains the fully symmetric basic rule sum
c	   indexed by the jth s-partition of the integers
c	   0,1,...,maxord-1.
c  fail	   integer failure output parameter
c	   fail=0 for successful termination of the subroutine
c	   fail=1 when numsms is too small for the subroutine to
c		   continue.  in this case weghts(1), weghts(2), ...,
c		   weghts(numsms), fulsms(1), fulsms(2), ...,
c		   fulsms(numsms) and intvls(1), intvls(2),...,
c		   intvls(j) are returned, where j is maximum value of
c		   maxord compatible with the given value of numsms.
c	   fail=2 when parameters s,minord, maxord or g are out of
c		   range
c***********************************************************************
cmmm	  external f
ctsl	real f
ctsl       double precision f
c***  for double precision change real to double precision
c      in the next statement
      integer d, i, fail, k(20), intcls, prtcnt, l, m(20), maxord,
     * minord, modofm, numsms, s, sumcls
      double precision intvls(maxord), center(s), hwidth(s), gisqrd,
     * glsqrd,
     * intmpa, intmpb, intval, one, fulsms(numsms), weghts(numsms),
     * two, momtol, momnkn, momprd(20,20), moment(20), zero, g(20)
      double precision flsm, wht
c	patterson generators
      data g(1), g(2) /0.0000000000000000,0.7745966692414833/
      data g(3), g(4) /0.9604912687080202,0.4342437493468025/
      data g(5), g(6) /0.9938319632127549,0.8884592328722569/
      data g(7), g(8) /0.6211029467372263,0.2233866864289668/
      data g(9), g(10), g(11), g(12) /0.1, 0.2, 0.3, 0.4/
c
c***  parameter checking and initialisation
      fail = 2
      maxrdm = 20
      maxs = 20
      if (s.gt.maxs .or. s.lt.1) return
      if (minord.lt.0 .or. minord.ge.maxord) return
      if (maxord.gt.maxrdm) return
      zero = 0
      one = 1
      two = 2
      momtol = one
   10 momtol = momtol/two
      if (momtol+one.gt.one) go to 10
      hundrd = 100
      momtol = hundrd*two*momtol
      d = minord
      if (d.eq.0) intcls = 0
c***  calculate moments and modified moments
      do 20 l=1,maxord
	floatl = l + l - 1
	moment(l) = two/floatl
   20 continue
      if (maxord.ne.1) then
         do 40 l=2,maxord
            intmpa = moment(l-1)
            glsqrd = g(l-1)**2
            do 30 i=l,maxord
               intmpb = moment(i)
               moment(i) = moment(i) - glsqrd*intmpa
               intmpa = intmpb
 30         continue
            if (moment(l)**2.lt.(momtol*moment(1))**2) moment(l) = zero
 40      continue
      endif
      do 70 l=1,maxord
	if (g(l).lt.zero) return
	momnkn = one
	momprd(l,1) = moment(1)
	if (maxord.eq.1) go to 70
	glsqrd = g(l)**2
	do 60 i=2,maxord
	  if (i.le.l) gisqrd = g(i-1)**2
	  if (i.gt.l) gisqrd = g(i)**2
	  if (glsqrd.eq.gisqrd) return
	  momnkn = momnkn/(glsqrd-gisqrd)
	  momprd(l,i) = momnkn*moment(i)
   60	continue
   70 continue
      fail = 1
c
c***  begin LOOP
c      for each d find all distinct partitions m with mod(m))=d
c
   80 prtcnt = 0
      intval = zero
      modofm = 0
      call nxprt(prtcnt, s, m)
   90 if (prtcnt.gt.numsms) return
c
c***  calculate weight for partition m and fully symmetric sums
c***	 when necessary
c
      if (d.eq.modofm) weghts(prtcnt) = zero
      if (d.eq.modofm) fulsms(prtcnt) = zero
      fulwgt = wht(s,moment,m,k,modofm,d,maxrdm,momprd)
      sumcls = 0
      if (weghts(prtcnt).eq.zero .and. fulwgt.ne.zero) fulsms(prtcnt) =
     * flsm(s, center, hwidth, moment, m, k, maxord, g, sumcls)
      intcls = intcls + sumcls
      intval = intval + fulwgt*fulsms(prtcnt)
      weghts(prtcnt) = weghts(prtcnt) + fulwgt
      call nxprt(prtcnt, s, m)
      if (m(1).gt.modofm) modofm = modofm + 1
      if (modofm.le.d) go to 90
c
c***  end loop for each d
      if (d.gt.0) intval = intvls(d) + intval
      intvls(d+1) = intval
      d = d + 1
      if (d.lt.maxord) go to 80
c
c***  set failure parameter and return
      fail = 0
      minord = maxord
      return
      end
