cMM this is the original adapt code with one modification.
cMM instead of calling the external function "functn", a fixed
cMM external routine adphlp is always called, and passed a pointer
cMM to the external s function.
cMM					Michael Meyer, October 1989.

      subroutine adapt(ndim,a,b,minpts,maxpts,eps,relerr,
     *	   lenwrk,wrkstr,finest,ifail)
c***begin prologue adapt
c  adaptive multidimensional integration subroutine
c	    author: A. C. Genz, Washington State University
c		     19 March 1984
c**************	 parameters for adapt  ********************************
c***** input parameters
c  ndim	   number of variables, must exceed 1, but not exceed 20
c  a	   real array of lower limits, with dimension ndim
c  b	   real array of upper limits, with dimension ndim
c  minpts  minimum number of function evaluations to be allowed.
c	   on the first call to adapt minpts should be set to a
c	   non negative value (caution... minpts is altered by adapt).
c	   It is possible to continue a calculation to greater accuracy
c	   by calling adapt again by decreasing eps (described below)
c	   and resetting minpts to any negative value.
c	   minpts must not exceed maxpts.
c  maxpts  maximum number of function evaluations to be allowed,
c	   which must be at least rulcls, where
c	   rulcls =  2**ndim+2*ndim**2+6*ndim+1
c
c	       for ndim =  2   3   4   5   6   7   8   9   10	12    15      20
c      maxpts >= rulcls = 25  45  73 113 173 269 433 729 1285 4457 33309 1049497
c
c	  a SUGGESTED value for maxpts is 100 times the above values.
c
c  functn  externally declared user defined function to be integrated.
c	   it must have parameters (ndim,z), where z is a real array
c	   of dimension ndim.
cTSL this function has been replaced by the fixed function adhlp
c  eps	   required relative accuracy
c  lenwrk  length of array wrkstr of working storage, the routine
c	   needs (2*ndim+3)*(1+maxpts/rulcls)/2 for lenwrk if
c	   maxpts function calls are used.
c	   for guidance, if you set maxpts to 100*rulcls (see table
c	   above) then acceptable values for lenwrk are
c
c	     for ndim = 2    3	  4    5    6	 7    8	    9
c	     lenwrk =  357  561	 1785 3417 6681 13209 26265 52377
c
c***** OUTPUT parameters
c
c  minpts  actual number of function evaluations used by adapt
c  wrkstr  real array of working storage of dimension (lenwrk).
c  relerr  estimated relative accuracy of finest
c  finest  estimated value of integral ["FINal ESTimate"]
c  ifail : return code
c
c	ifail=0 for normal exit, when estimated relative accuracy
c		relerr is less than eps with maxpts or less function
c		calls made.
c	ifail=1 if maxpts was too small for adapt to obtain the
c		required relative accuracy eps.
c		In this case adapt returns a value of finest
c		with estimated relative accuracy relerr.
c	ifail=2 if lenwrk too small for maxpts function calls.
c		In this case adapt returns a value of finest with
c		estimated accuracy relerr using the working storage
c		available, but relerr will be greater than eps.
c	ifail=3 if ndim < 2, ndim > 20,
c	ifail=4 if minpts > maxpts,
c	ifail=5 if maxpts < rulcls or other memory problems
c		      (which will only be found later)
c***********************************************************************
c***end prologue adapt

      implicit none

C-- Arguments:
C      double precision functn
C      external functn

      integer ndim, minpts,maxpts, lenwrk, ifail
      double precision a(ndim), b(ndim), eps, relerr, wrkstr(lenwrk),
     &	   finest

C-- Local Variables:
      double precision center(20), width(20)
      double precision errmin, rgnerr, rgnval, half, zero,one,two

      integer divaxo, divaxn, divflg, funcls, index1, index2,
     *     j, k,  maxcls, rgnstr, rulcls, sbrgns, sbtmpp, subrgn, subtmp

      data zero/0d0/, one/1d0/, two/2d0/

c Check arguments; fail w/ code '3' or '4'
      relerr=one
      funcls=0
      ifail=3
      if(ndim.lt.2.or.ndim.gt.20) goto 990
      ifail=4
      if(minpts.gt.maxpts) goto 990
      ifail=5
c
c*****	initialisation of subroutine
c
      half=one/two
      rgnstr =2*ndim+3
      errmin = zero
      maxcls = 2**ndim + 2*ndim**2 + 6*ndim+1
      maxcls = min0(maxcls,maxpts)
      divaxo=0
c
c*****	end subroutine initialisation
      if(minpts.lt.0) then
         sbrgns=wrkstr(lenwrk-1)
         goto 280
      endif

      do 30 j=1,ndim
	 width(j)=(b(j)-a(j))*half
 30	 center(j)=a(j)+width(j)
      finest=zero
      wrkstr(lenwrk)=zero
      divflg=1
      subrgn=rgnstr
      sbrgns=rgnstr

C-- REPEAT --- (outermost loop) -------
   40 call bsrl(ndim,center,width,maxcls,rulcls,
     *	   errmin,rgnerr,rgnval,divaxo,divaxn)
      finest=finest+rgnval
      wrkstr(lenwrk)=wrkstr(lenwrk)+rgnerr
      funcls = funcls + rulcls
c
c*****	place results of basic rule into partially ordered list
c*****	according to subregion error
      if(divflg .eq. 0) then
c
c*****	when divflg=0 start at top of list and move down list tree to
c     find correct position for results from first half of recently
c     divided subregion
 200	 subtmp=2*subrgn
	 if(subtmp.le.sbrgns) then
	    if(subtmp.ne.sbrgns) then
	       sbtmpp=subtmp+rgnstr
	       if(wrkstr(subtmp).lt.wrkstr(sbtmpp)) subtmp=sbtmpp
	    endif
 210	    if(rgnerr.lt.wrkstr(subtmp)) then
	       do 220 k=1,rgnstr
		  index1=subrgn-k+1
		  index2=subtmp-k+1
		  wrkstr(index1)=wrkstr(index2)
 220	       continue
	       subrgn=subtmp

	       goto 200
	    endif
	 endif

      else
c
c*****when divflg=1 start at bottom right branch and move up list
c     tree to find correct position for results from second half of
c     recently divided subregion
 230	 subtmp=(subrgn/(rgnstr*2))*rgnstr
	 if(subtmp.ge.rgnstr) then
	    if(rgnerr.gt.wrkstr(subtmp)) then
	       do 240 k=1,rgnstr
		  index1=subrgn-k+1
		  index2=subtmp-k+1
		  wrkstr(index1)=wrkstr(index2)
 240	       continue
	       subrgn=subtmp
	       goto 230
	    endif
	 endif
      endif

c*****	store results of basic rule in correct position in list
  250 wrkstr(subrgn)=rgnerr
      wrkstr(subrgn-1)=rgnval
      wrkstr(subrgn-2)=divaxn
      do 260 j=1,ndim
	subtmp=subrgn-2*(j+1)
	wrkstr(subtmp+1)=center(j)
	wrkstr(subtmp)=width(j)
 260  continue
      if(divflg .eq. 0) then
c***  when divflg=0 prepare for second application of basic rule
	 center(divaxo)=center(divaxo)+two*width(divaxo)
	 sbrgns=sbrgns+rgnstr
	 subrgn=sbrgns
	 divflg=1
c***  loop back to apply basic rule to other half of subregion
	 go to 40
      endif
c
c*****	end ordering and storage of basic rule results
c*****	make checks for possible termination of routine
c
  270 relerr=one
      if(wrkstr(lenwrk).le.zero) wrkstr(lenwrk)=zero
      if(dabs(finest).ne.zero) relerr=wrkstr(lenwrk)/dabs(finest)
      if(relerr.gt.one) relerr=one

      if(sbrgns+rgnstr.gt.lenwrk-2) ifail=2
      if(funcls+funcls*rgnstr/sbrgns.gt.maxpts) ifail=1
      if(relerr.lt.eps.and.funcls.ge.minpts) ifail=0
      if(ifail.lt.3) goto 990
c
c*****	prepare to use basic rule on each half of subregion with largest
c	error
  280 divflg=0
      subrgn=rgnstr
      subtmp = 2*sbrgns/rgnstr
      maxcls = maxpts/subtmp
      errmin = dabs(finest)*eps/dfloat(subtmp)
      wrkstr(lenwrk)=wrkstr(lenwrk)-wrkstr(subrgn)
      finest=finest-wrkstr(subrgn-1)
      divaxo=wrkstr(subrgn-2)
      do 290 j=1,ndim
	subtmp=subrgn-2*(j+1)
	center(j)=wrkstr(subtmp+1)
  290	width(j)=wrkstr(subtmp)
      width(divaxo)=width(divaxo)*half
      center(divaxo)=center(divaxo)-width(divaxo)
c
c*****	loop back to apply basic rule
c
      goto 40
c
c*****	termination point
c
  990 minpts=funcls
      wrkstr(lenwrk-1)=sbrgns
      return
      end
