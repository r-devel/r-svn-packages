      double precision function wht(s, intrps, m, k, modofm, d,
     *     maxrdm, momprd)
c***  subroutine to calculate weight for partition m
c
      integer s, m(s), k(s), d, maxrdm, mi, ki, m1, k1, modofm
      double precision intrps(s), zero, momprd(maxrdm,maxrdm)

      zero = 0
      do 10 i=1,s
        intrps(i) = zero
        k(i) = 0
   10 continue
      m1 = m(1) + 1
      k1 = d - modofm + m1
   20 intrps(1) = momprd(m1,k1)
      if (s.eq.1) go to 40
      do 30 i=2,s
        mi = m(i) + 1
        ki = k(i) + mi
        intrps(i) = intrps(i) + momprd(mi,ki)*intrps(i-1)
        intrps(i-1) = zero
        k1 = k1 - 1
        k(i) = k(i) + 1
        if (k1.ge.m1) go to 20
        k1 = k1 + k(i)
        k(i) = 0
   30 continue
   40 wht = intrps(s)
      return
      end
