      double precision function flsm(s,center,hwidth,x,m,mp,maxord,
     * g,sumcls)
c
c***  function to compute fully symmetric basic rule sum
c
      integer s, m(s), mp(s), maxord, sumcls, ixchng, lxchng, i, l,
     * ihalf, mpi, mpl
      double precision g(maxord), x(s), intwgt, zero, one,two, intsum,
     * center(s), hwidth(s)
      double precision adphlp

      zero = 0
      one = 1
      two = 2

      intwgt = one
      do 10 i=1,s
        mp(i) = m(i)
        if (m(i).ne.0) intwgt = intwgt/two
        intwgt = intwgt*hwidth(i)
   10 continue
      sumcls = 0
      flsm = zero
c
c*******  compute centrally symmetric sum for permutation mp
   20 intsum = zero
      do 30 i=1,s
        mpi = mp(i) + 1
        x(i) = center(i) + g(mpi)*hwidth(i)
   30 continue
   40 sumcls = sumcls + 1
cmmm
      intsum = intsum + adphlp(s,x)
      do 50 i=1,s
        mpi = mp(i) + 1
        if(g(mpi).ne.zero) hwidth(i) = -hwidth(i)
        x(i) = center(i) + g(mpi)*hwidth(i)
        if (x(i).lt.center(i)) go to 40
   50 continue
c*******  end integration loop for mp
c
      flsm = flsm + intwgt*intsum
      if (s.eq.1) return
c
c*******  find next distinct permutation of m and loop back
c          to compute next centrally symmetric sum
      do 80 i=2,s
        if (mp(i-1).le.mp(i)) go to 80
        mpi = mp(i)
        ixchng = i - 1
        if (i.eq.2) go to 70
        ihalf = ixchng/2
        do 60 l=1,ihalf
          mpl = mp(l)
          imnusl = i - l
          mp(l) = mp(imnusl)
          mp(imnusl) = mpl
          if (mpl.le.mpi) ixchng = ixchng - 1
          if (mp(l).gt.mpi) lxchng = l
   60   continue
        if (mp(ixchng).le.mpi) ixchng = lxchng
   70   mp(i) = mp(ixchng)
        mp(ixchng) = mpi
        go to 20
   80 continue
c*****  end loop for permutations of m and associated sums
c
      return
      end
