      subroutine nxprt(prtcnt, s, m)
c
c***  subroutine to compute the next s partition
c
      implicit none
      integer s, m(s), prtcnt

      integer i,k, msum

      if (prtcnt.gt.0) go to 20
      do 10 i=1,s
        m(i) = 0
   10 continue
      prtcnt = 1
      return
   20 prtcnt = prtcnt + 1
      msum = m(1)
      if (s.eq.1) go to 60
      do 50 i=2,s
        msum = msum + m(i)
        if (m(1).le.m(i)+1) go to 40
        m(1) = msum - (i-1)*(m(i)+1)
        do 30 k=2,i
          m(k) = m(i) + 1
   30   continue
        return
   40   m(i) = 0
   50 continue
   60 m(1) = msum + 1
      return
      end
