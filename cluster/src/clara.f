c     Clustering LARge Applications
c     ~          ~~~   ~
c     Clustering program based upon the k-medoid approach,
c     and suitable for data sets of at least 100 objects.
c     (for smaller data sets, please use program pam.)
c     
      subroutine clara(nn,jpp,kk,x,nran,nsam,dys,mdata,valmd,jtmd,ndyst,
     f     nrepr,nsel,nbest,nr,nrx,radus,ttd,ratt,ttbes,rdbes,rabes,
     f     mtt,azba,avsyl,ttsyl,sylinf,jstop,
     f     tmp1,tmp2,tmp3,ntmp1,ntmp2,ntmp3,ntmp4,ntmp5,ntmp6)

      implicit none
      integer nn, jpp, kk, nran, nsam, mdata, ndyst, jstop
c     nn   = number of objects
c     jpp  = number of variables
c     kk   = number of clusters
c     nran = ??
c     nsam = number of objects drawn from data set
      double precision x(nn*jpp), dys(1 + nsam*(nsam-1)/2),
     1     valmd(jpp), radus(kk),ttd(kk),ratt(kk), 
     2     ttbes(kk),rdbes(kk),rabes(kk), azba, avsyl(kk), ttsyl,
     3     sylinf(nsam,4), tmp1(nsam),tmp2(nsam),tmp3(nsam)
      integer jtmd(jpp), nrepr(nsam),nsel(nsam),nbest(nsam),
     1     nr(kk),nrx(kk), mtt(kk),
     2     ntmp1(nsam),ntmp2(nsam),ntmp3(nsam),
     3     ntmp4(nsam),ntmp5(nsam),ntmp6(nsam)
c Var      
      integer j,jjb,jk,jkk,jn,js,jsm,jran,jhalt, kall,kans,kkm,kkp,kran,
     1     l,less, nad,nadv,nadvp,nafs,nexap,nexbp, nhalf, nneq,nnpp,
     2     nrun,nsamb,nsub, nunfs,nsm,ntt
      double precision rnn, ran, zba, s,sx,z,zb

      jstop=0
      rnn=nn
      if(nn.eq.nsam) then
         nneq=1
      else
         nneq=0
      endif
      nhalf= nsam*(nsam-1)/2 + 1
      nsamb=2*nsam
      if(nn.lt.nsamb) then
         less=nn-nsam
      else
         less=nsam
      endif
      nnpp=nn*jpp

      nunfs=0
      kall=0
      nrun=0
c     
c     Loop :  random subsamples are drawn and partitioned into kk clusters
      do 400 jran=1,nran
         jhalt=0
         if(nneq.eq.0)go to 140
         if(nneq.eq.2)go to 400
c     else nneq == 1 when above nn == nsam :
         nneq=2
         do 130 j=1,nsam
            nsel(j)=j
 130     continue
         go to 320

c     nneq = 0 :
 140     ntt=0
         if(jran.ne.1 .and. nunfs.ne.jran .and. nn.ge.nsamb) then
            do 150 jk=1,kk
               nsel(jk)=nrx(jk)
 150        continue
            kkm=kk-1
            do 170 jk=1,kkm
               nsm=nsel(jk)
               kkp=jk+1
               jsm=jk
               do 160 jkk=kkp,kk
                  if(nsel(jkk).ge.nsm)go to 160
                  nsm=nsel(jkk)
                  jsm=jkk
 160           continue
               nsel(jsm)=nsel(jk)
               nsel(jk)=nsm
 170        continue
            ntt=kk

         else

 180        call randm(nrun,ran)
            kran=rnn*ran+1.
            if(kran.gt.nn)kran=nn
            if(jran.ne.1) then
               do 190 jk=1,kk
                  if(kran.eq.nrx(jk))go to 180
 190           continue
            endif
 200        ntt=ntt+1
            nsel(ntt)=kran
            if (less.eq.ntt)go to 290
         endif
C     Loop
 210     call randm(nrun,ran)
         kran=rnn*ran+1.
         if(kran.gt.nn)kran=nn
         if(jran.eq.1)go to 230
         if(nn.ge.nsamb)go to 230
         do 220 jk=1,kk
            if(kran.eq.nrx(jk))go to 210
 220     continue
 230     do 260 kans=1,ntt
            if(nsel(kans).lt.kran)go to 260
            if(nsel(kans).eq.kran)go to 210
            go to 270
 260     continue

         ntt=ntt+1
         nsel(ntt)=kran
         go to 290

 270     do 280 nad=kans,ntt
            nadv=ntt-nad+kans
            nadvp=nadv+1
            nsel(nadvp)=nsel(nadv)
 280     continue
         ntt=ntt+1
         nsel(kans)=kran

 290     if(ntt.lt.less)go to 210
         if(nn.ge.nsamb)go to 320
         nexap=1
         nexbp=1
c     do 305  jn=1, nn
         jn=0
 300     jn=jn+1
         if(nsel(nexap).eq.jn)then
            nexap=nexap+1
         else
            nrepr(nexbp)=jn
            nexbp=nexbp+1
         endif
         if(jn.lt.nn)go to 300
c     305  continue

         do 310 nsub=1,nsam
            nsel(nsub)=nrepr(nsub)
 310     continue

 320     call dysta2(nsam,jpp,nsel,x,nn,dys,ndyst,jtmd,valmd, jhalt)
         if(jhalt.eq.1)go to 400
         kall=1
         s=0.0

         l=1
 340     l=l+1
         if(dys(l).gt.s)s=dys(l)
         if(l.lt.nhalf)go to 340

         call bswap2(kk,nsam,nrepr,dys,z,s,tmp1,tmp2,tmp3)
         call selec(kk,nn,jpp,ndyst,zb,nsam,mdata,
     f        jtmd,valmd,nrepr,nsel,dys,x,nr,nafs,ttd,radus,ratt,
     f        ntmp1,ntmp2,ntmp3,ntmp4,ntmp5,ntmp6, tmp1,tmp2)
         nunfs=nunfs+nafs
         if(nafs.eq.1)go to 400
         if(jran.ne.1) then
            if(zb.ge.zba) go to 400
         endif
         zba=zb
         do 345 jjb=1,kk
            ttbes(jjb)=ttd(jjb)
            rdbes(jjb)=radus(jjb)
            rabes(jjb)=ratt(jjb)
 345     continue
         do 360 jk=1,kk
            nrx(jk)=nr(jk)
 360     continue
         do 370 js=1,nsam
            nbest(js)=nsel(js)
 370     continue
         sx=s
 400  continue
c--- end random sampling loop 

      if(nunfs.ge.nran) then
         jstop=1
         return
      endif
c     
c     for the best subsample, the objects of the entire data set
c     are assigned to their clusters 
c     
 450  if(kall.ne.1) then
         jstop=2
         return
      endif
 460  azba=zba/rnn
 470  call dysta2(nsam,jpp,nbest,x,nn,dys,ndyst,jtmd,valmd, jhalt)
      call resul(kk,nn,jpp,ndyst,mdata,jtmd,valmd, x,nrx,mtt)
      if(kk.gt.1) then
 480     call black(kk,jpp,nn,nsam,nbest,dys,sx,x,avsyl,ttsyl,sylinf,
     +        ntmp1,ntmp2,ntmp3,ntmp4,tmp1,tmp2)
      endif
 500  end
c     --- of clara() -----------------------------------------------------
c     
c     
      subroutine dysta2(nsam,jpp,nsel,x,nn,dys,ndyst,jtmd, valmd,jhalt)

      implicit none

      integer nsam,jpp, nn, nsel(nsam), ndyst, jtmd(jpp), jhalt
      double precision x(nn*jpp), dys(1+nsam*(nsam-1)/2), valmd(jpp) 
c
      double precision pp, rpres, clk
      integer j,k,l, nlk, lsubt,lsel, ksel,numlj,numkj, npres
c
      pp=jpp
      nlk=1
      dys(1)=0.0
      do 100 l=2,nsam
         lsubt=l-1
         lsel=nsel(l)
         do 20 k=1,lsubt
            ksel=nsel(k)
            clk=0.0
            nlk=nlk+1
            npres=0
            do 30 j=1,jpp
               numlj=(lsel-1)*jpp+j
               numkj=(ksel-1)*jpp+j
               if(jtmd(j).lt.0) then
                  if(x(numlj).eq.valmd(j))go to 30
                  if(x(numkj).eq.valmd(j))go to 30
               endif
 40            npres=npres+1
               if(ndyst.eq.1)then
                  clk=clk+(x(numlj)-x(numkj))*(x(numlj)-x(numkj))
               else
                  clk=clk+dabs(x(numlj)-x(numkj))
               endif
 30         continue
            rpres=npres
            if(npres.ne.0)go to 60
            jhalt=1
            dys(nlk)=-1.0
            go to 20
 60         if(ndyst.ne.1)go to 70
            dys(nlk)=dsqrt(clk*(pp/rpres))
            go to 20
 70         dys(nlk)=clk*(pp/rpres)
 20      continue
 100  continue
      end
c     -----------------------------------------------------------
c
      subroutine randm(nrun,ran)
c   we programmed this generator ourselves because we wanted it
c   to be machine independent. it should run on most computers
c   because the largest integer used is less than 2**30 . the period
c   is 2**16=65536, which is good enough for our purposes.

      implicit none
      integer nrun
      double precision ran

      integer k
      double precision ry

      nrun=nrun*5761+999
      k=nrun/65536
      nrun=nrun-k*65536
      ry=nrun
      ran=ry/65536.0
      return
      end
c     -----------------------------------------------------------
c     
      subroutine bswap2(kk,nsam,nrepr,dys,sky,s,dysma,dysmb,beter)

      implicit none
      integer kk,nsam, nrepr(nsam)
      double precision dys(1+nsam*(nsam-1)/2), sky,s,
     1     dysma(nsam),dysmb(nsam),beter(nsam)
c
      integer meet
      external meet
c Var      
      logical nafs
      integer j,jk,jn,jkabc, newf, na,nb,nrjk
      double precision zb,pp,dsum,dnull, tra, pres, abc

c     
c     nafs = .true. if a distance cannot be calculated
c FIXME: nafs is not really used; maybe should really be returned! (MM)
      nafs=.false.
c     
c     identification of representative objects, and initializations
c     
      jk=0
      do 10 j=1,nsam
         if(nrepr(j).eq.0)go to 10
         jk=jk+1
         nr(jk)=nsel(j)
         ns(jk)=0
         ttd(jk)=0.
         radus(jk)=-1.
         np(jk)=j
 10   continue
c     
c     assignment of the objects of the entire data set to a cluster,
c     computation of some statistics, determination of the
c     new ordering of the clusters
c     
      zb=0.
      pp=jpp
      newf=0
      jn=0
 15   jn=jn+1

      if(mdata.eq.0) then
         do 30 jk=1,kk
            dsum=0.
            nrjk=nr(jk)
            if (nrjk.ne.jn) then
               do 20 jp=1,jpp
                  na=(nrjk-1)*jpp+jp
                  nb=(jn-1)*jpp+jp
                  tra=dabs(x(na)-x(nb))
                  if(ndyst.eq.1)tra=tra*tra
                  dsum=dsum+tra
 20            continue
               if(jk.ne.1) then
                  if(dsum.ge.dnull)go to 30
               endif
            endif
 25         dnull=dsum
            jkabc=jk
 30      continue

      else

 40      pres=0.
         do 70 jk=1,kk
            dsum=0.
            nrjk=nr(jk)
            if (nrjk.ne.jn) then
               abc=0.
               do 50 jp=1,jpp
                  na=(nrjk-1)*jpp+jp
                  nb=(jn-1)*jpp+jp
                  if(jtmd(jp).lt.0)then
                     if(x(na).eq.valmd(jp))go to 50
                     if(x(nb).eq.valmd(jp))go to 50
                  endif
 45               abc=abc+1.
                  tra=dabs(x(na)-x(nb))
                  if(ndyst.eq.1)tra=tra*tra
                  dsum=dsum+tra
 50            continue
               if(abc.lt.0.5)go to 70
               dsum=dsum*abc/pp
            endif
 64         if(pres.le.0.5) then
               pres=1.
            else
               if(dsum.ge.dnull)go to 70
            endif
 65         dnull=dsum
            jkabc=jk
 70      continue
         if(pres.gt.0.5)go to 80
         nafs=.true.
         return
      endif

 80   if(ndyst.eq.1)dnull=dsqrt(dnull)
      zb=zb+dnull
      ttd(jkabc)=ttd(jkabc)+dnull
      if(dnull.gt.radus(jkabc))radus(jkabc)=dnull
      ns(jkabc)=ns(jkabc)+1
      if(newf.lt.kk) then
         if(newf.ne.0) then
            do 82 jnew=1,newf
               if(jkabc.eq.new(jnew))go to 90
 82         continue
         endif
         newf=newf+1
         new(newf)=jkabc
      endif
 90   if(jn.lt.nn)go to 15

c     
c     a permutation is carried out on vectors nr,ns,np,ttd,radus
c     using the information in vector new.
c     
      do 92 jk=1,kk
         njk=new(jk)
         nrnew(jk)=nr(njk)
         nsnew(jk)=ns(njk)
         npnew(jk)=np(njk)
         ttnew(jk)=ttd(njk)
         rdnew(jk)=radus(njk)
 92   continue
      do 94 jk=1,kk
         nr(jk)=nrnew(jk)
         ns(jk)=nsnew(jk)
         np(jk)=npnew(jk)
         ttd(jk)=ttnew(jk)
         radus(jk)=rdnew(jk)
 94   continue

      do 101 j=1,kk
         rns=ns(j)
         ttd(j)=ttd(j)/rns
 101  continue
      if(kk.eq.1)go to 150
c     
c     computation of minimal distance of medoid ka to any
c     other medoid for comparison with the radius of cluster ka.
c     
      do 120 ka=1,kk
         nstrt=0
         npa=np(ka)
         do 110 kb=1,kk
            if(kb.eq.ka)go to 110
            npb=np(kb)
            npab=meet(npa,npb)
            if(nstrt.eq.0) then
               nstrt=1
            else
               if(dys(npab).ge.ratt(ka))go to 110
            endif
 104        ratt(ka)=dys(npab)
            if(ratt(ka).ne.0.)go to 110
            ratt(ka)=-1.
 110     continue
         if(ratt(ka).gt. -0.5) ratt(ka)=radus(ka)/ratt(ka)
 120  continue

 150  return
      end
c     -----------------------------------------------------------
c     
      subroutine resul(kk,nn,jpp,ndyst,mdata,jtmd,
     f     valmd,x,nrx,mtt)

      implicit double precision (a-h,o-z)
      dimension x(nn*jpp),nrx(kk),jtmd(jpp),valmd(jpp),mtt(kk)
      pp=jpp
c     
c     clustering vector is incorporated into x, and printed.
c     
      jn=0
 100  jn=jn+1
      njnb=(jn-1)*jpp
      do 145 jk=1,kk
         if(nrx(jk).eq.jn)go to 220
 145  continue
      jna=(jn-1)*jpp+1
      if(mdata.ne.0)go to 170
      do 160 jk=1,kk
         dsum=0.
         nrjk=(nrx(jk)-1)*jpp
         do 150 j=1,jpp
            na=nrjk+j
            nb=njnb+j
            tra=dabs(x(na)-x(nb))
            if(ndyst.eq.1)tra=tra*tra
            dsum=dsum+tra
 150     continue
         if(ndyst.eq.1)dsum=dsqrt(dsum)
         if(jk.eq.1)dnull=dsum+0.1
         if(dsum.ge.dnull)go to 160
         dnull=dsum
         jksky=jk
 160  continue
      go to 200
 170  do 190 jk=1,kk
         dsum=0.
         nrjk=(nrx(jk)-1)*jpp
         abc=0.
         do 180 j=1,jpp
            na=nrjk+j
            nb=njnb+j
            if(jtmd(j).ge.0)go to 185
            if(x(na).eq.valmd(j))go to 180
            if(x(nb).eq.valmd(j))go to 180
 185        abc=abc+1.
            tra=dabs(x(na)-x(nb))
            if(ndyst.eq.1)tra=tra*tra
            dsum=dsum+tra
 180     continue
         if(ndyst.eq.1)dsum=dsqrt(dsum)
         dsum=dsum*abc/pp
         if(jk.eq.1)dnull=dsum+0.1
         if(dsum.ge.dnull)go to 190
         dnull=dsum
         jksky=jk
 190  continue
 200  x(jna)=jksky
 220  if(jn.lt.nn)go to 100
      do 230 jk=1,kk
         nrjk=nrx(jk)
         nrjka=(nrjk-1)*jpp+1
         x(nrjka)=jk
 230  continue

 300  do 320 ka=1,kk
         mtt(ka)=0
         j=0
 325     j=j+1
         ja=(j-1)*jpp+1
         nxja=idint(x(ja)+0.1)
         if(nxja.eq.ka)mtt(ka)=mtt(ka)+1
         if(j.lt.nn)go to 325
 320  continue
 330  end
c     -----------------------------------------------------------
c     
      subroutine black(kk,jpp,nn,nsam,nbest,dys,sx,x,avsyl,ttsyl,sylinf,
     f     ncluv,nsend,nelem,negbr,syl,srank)

      implicit double precision (a-h,o-z)

      dimension ncluv(nsam),nsend(nsam),nelem(nsam),negbr(nsam)
      dimension syl(nsam),srank(nsam),avsyl(kk),nbest(nsam)
      dimension x(nn*jpp),dys(1+nsam*(nsam-1)/2),sylinf(nsam,4)
c     
c     construction of clustering vector (ncluv)
c     of selected sample (nbest).
c     
      do 12 l=1,nsam
         ncase=nbest(l)
         jna=(ncase-1)*jpp+1
         ncluv(l)=idint(x(jna)+0.1)
 12   continue
c     
c     drawing of the silhouettes
c     
      nsylr=0
      ttsyl=0.0
      do 100 numcl=1,kk
         ntt=0
         do 30 j=1,nsam
            if(ncluv(j).ne.numcl)go to 30
            ntt=ntt+1
            nelem(ntt)=j
 30      continue

         do 40 j=1,ntt
            nj=nelem(j)
            dysb=1.1*sx+1.0
            negbr(j)=-1
            do 41 nclu=1,kk
               if(nclu.eq.numcl)go to 41
               nbb=0
               db=0.0
               do 43 l=1,nsam
                  if(ncluv(l).ne.nclu)go to 43
                  nbb=nbb+1
                  db=db + dys(meet(nj,l))
 43            continue
               btt=nbb
               db=db/btt
               if(db.ge.dysb)go to 41
               dysb=db
               negbr(j)=nclu
 41         continue
            if(ntt.eq.1)go to 50
            dysa=0.0
            do 45 l=1,ntt
               nl=nelem(l)
               dysa=dysa + dys(meet(nj,nl))
 45         continue
            att=ntt-1
            dysa=dysa/att
            if(dysa.gt.0.0)go to 51
            if(dysb.gt.0.0)go to 52
 50         syl(j)=0.0
            go to 40
 52         syl(j)=1.0
            go to 40
 51         if(dysb.gt.0.)then
               if(dysb.gt.dysa)syl(j)=1.0-dysa/dysb
               if(dysb.lt.dysa)syl(j)=dysb/dysa-1.0
               if(dysb.eq.dysa)syl(j)=0.0
            else
 53            syl(j)=-1.0
            endif
 54         if(syl(j).le.(-1.0))syl(j)=-1.0
            if(syl(j).ge.1.0)syl(j)=1.0
 40      continue

         avsyl(numcl)=0.0
         do 60 j=1,ntt
            symax= -2.
            do 70 l=1,ntt
               if(syl(l).le.symax)go to 70
               symax=syl(l)
               lang=l
 70         continue
            nsend(j)=lang
            srank(j)=syl(lang)
            avsyl(numcl)=avsyl(numcl)+srank(j)
            syl(lang)= -3.
 60      continue
         ttsyl=ttsyl+avsyl(numcl)
         rtt=ntt
         avsyl(numcl)=avsyl(numcl)/rtt

         if(ntt.ge.2)goto 75
         ncase=nelem(1)
         nsylr=nsylr+1
         sylinf(nsylr,1)=numcl
         sylinf(nsylr,2)=negbr(1)
         sylinf(nsylr,3)=0.0
         sylinf(nsylr,4)=nbest(ncase)
         goto 100
 75      do 80 l=1,ntt
            lplac=nsend(l)
            ncase=nelem(lplac)
            nsylr=nsylr+1
            sylinf(nsylr,1)=numcl
            sylinf(nsylr,2)=negbr(lplac)
            sylinf(nsylr,3)=srank(l)
            sylinf(nsylr,4)=nbest(ncase)
 80      continue
 100  continue
      rsam=nsam
      ttsyl=ttsyl/rsam
      end
