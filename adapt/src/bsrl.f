      subroutine bsrl(s, center,hwidth, maxvls,funcls,
     *                  errmin,errest,basest,divaxo,divaxn)

      implicit none
C-- Arguments:
      integer s
      double precision center(s), hwidth(s)
      integer maxvls,funcls, divaxo,divaxn
      double precision errmin, errest, basest
C
      EXTERNAL adphlp
      double precision adphlp
C-- Local Variables:
      double precision  intvls(20), z(20), fulsms(200), weghts(200)

      integer intcls, i, mindeg, maxdeg, maxord, minord
      integer ifail

      double precision zero, one, two, ten, dif, errorm,
     *     sum0, sum1, sum2, difmax, x1, x2

      data zero/0d0/, one/1d0/, two/2d0/, ten/10d0/

      maxdeg = 12
      mindeg = 4
      minord = 0
      do 10 maxord = mindeg,maxdeg
         call symrl(s, center, hwidth, minord, maxord, intvls,
     *        intcls, 200, weghts, fulsms, ifail)
         if (ifail.eq.2) goto 20
         errest = dabs(intvls(maxord)  -intvls(maxord-1))
         errorm = dabs(intvls(maxord-1)-intvls(maxord-2))
         if (errest.ne.zero)
     &        errest = errest*
     &        dmax1(one/ten,errest/dmax1(errest/two,errorm))
         if (errorm.le. 5.*errest) goto 20
         if (2*intcls.gt.maxvls) goto 20
         if (errest.lt.errmin) goto 20
 10   continue
 20   difmax = -1
      x1 = one/two**2
      x2 = 3.*x1
      do 30 i = 1,s
         z(i) = center(i)
 30   continue
cmmm
      sum0 = adphlp(s,z)
      do 40 i = 1,s
         z(i) = center(i) - x1*hwidth(i)
cmmm
         sum1 = adphlp(s,z)
         z(i) = center(i) + x1*hwidth(i)

         sum1 = sum1 + adphlp(s,z)
         z(i) = center(i) - x2*hwidth(i)

         sum2 = adphlp(s,z)
         z(i) = center(i) + x2*hwidth(i)

         sum2 = sum2 + adphlp(s,z)
         z(i) = center(i)

         dif = dabs((sum1-two*sum0) - (x1/x2)**2*(sum2-two*sum0))
         if (dif.ge.difmax) then
            difmax = dif
            divaxn = i
         endif
 40   continue
      if (sum0.eq.sum0+difmax/two) divaxn = mod(divaxo,s) + 1
      basest = intvls(minord)
      funcls = intcls + 4*s
      return
      end
