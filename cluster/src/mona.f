      subroutine mona(nn,pp, x,jerr, nban,ner,kwan,lava, jlack)
cc
cc   MONothetic Analysis
cc
cc   Program for divisive hierarchical clustering of binary data,
cc   using association analysis.
cc
cc   list of functions and subroutines:
cc       function kab

c Args
      integer nn, pp, jerr
cc	nn   = number of objects
cc      pp   = number of variables
cc      jerr : error return code in {1,2,3,4}
      integer x(nn,pp), jlack(pp), nban(nn),ner(nn),kwan(nn),lava(nn)

c Function called:
      integer kab
c VARs
      integer j, ja, jb, jnat, jma, jtel, jtelz, jtel2, jres
      integer j0,j1, jptwe
      integer jva, jvb, jvc, jvd
      integer k, ka, kb, kal, kalf, kva, kvb, kvc, kvd, km
      integer l, lbb, laa, lcc, ldd, lee, lack,lama, lams
      integer nel, nelbb, nzf, nhalf, npass, nclu,nsyn, myst, mysca

      lack=0
      nhalf=(nn+1)/2
      jptwe=(pp+4)/5
      myst=0
      do 70 l=1,nn
         mysca=0
         do 60 j=1,pp
            if(x(l,j).eq.0)go to 60
            if(x(l,j).eq.1)go to 60
c	    else: missing
	    mysca=mysca+1
 60	 continue
	 if(mysca .eq. pp) then
c     all variables missing for this object
	    jerr=1
	    return
	 endif
	 myst=myst+mysca
 70   continue
      if(myst.eq.0)go to 290
      do 100 j=1,pp
         j0=0
         j1=0
         do 80 l=1,nn
            if(x(l,j).eq.0)j0=j0+1
            if(x(l,j).eq.1)j1=j1+1
 80      continue
         jlack(j)=nn-j0-j1
         if(jlack(j).ne.0)lack=lack+1
         if(jlack(j).ge.nhalf) then
c     at least 50% of the objects have missing values for this variable
	    jerr=2
	    return
	 endif
 90	 if(j0.eq.0 .or. j1.eq.0) then
c	    all non missing values are identical for this variable
	    jerr=3
	    return
	 endif
 100  continue

      if(lack .eq. pp) then
c     all variables have missing values
         jerr=4
         return
      endif
cc
cc   filling in missing values
cc
      do 260 j=1,pp
         if(jlack(j).eq.0)go to 260
         lama=-1
         nsyn=1
         do 240 ja=1,pp
            if(jlack(ja).ne.0)go to 240
            jva=0
            jvb=0
            jvc=0
            jvd=0
            do 230 k=1,nn
               if(x(k,j).eq.1)go to 220
               if(x(k,ja).eq.0)jva=jva+1
               if(x(k,ja).eq.1)jvb=jvb+1
               go to 230
 220           if(x(k,ja).eq.0)jvc=jvc+1
               if(x(k,ja).eq.1)jvd=jvd+1
 230        continue
            kal=jva*jvd-jvb*jvc
            kalf=kab(kal)
            if(kalf.ge.lama)then
               lama=kalf
               jma=ja
               if(kal.lt.0) nsyn=-1
            endif
 240     continue
         do 250 l=1,nn
            if(x(l,j).eq.0)go to 250
            if(x(l,j).eq.1)go to 250
            if(nsyn.eq.1)then
               x(l,j)=x(l,jma)
            else
               if(x(l,jma).eq.1)x(l,j)=0
               if(x(l,jma).eq.0)x(l,j)=1
            endif
 250     continue
 260  continue
c--- end of treating missing values ----

cc
cc   initialization
cc
 290  do 300 k=1,nn
         kwan(k)=0
         ner(k)=k
         lava(k)=0
 300  continue
      npass=1
      kwan(1)=nn
cc
cc   algorithm
cc
      nclu=1
      ka=1
C --- Loop ---
 310  kb=ka+kwan(ka)-1
      lama=-1
      jnat=pp
      do 370 j=1,pp
         if(nclu.eq.1)go to 330
         j0=0
         j1=0
         do 325 l=ka,kb
            nel=ner(l)
            if(x(nel,j).eq.0)j0=j0+1
            if(x(nel,j).eq.1)j1=j1+1
 325     continue
         if(j1.eq.0)go to 370
         if(j0.eq.0)go to 370
 330     jnat=jnat-1
         lams=0
         do 360 jb=1,pp
            if(jb.eq.j)go to 360
            kva=0
            kvb=0
            kvc=0
            kvd=0
            do 350 l=ka,kb
               nel=ner(l)
               if(x(nel,j).eq.1)go to 340
               if(x(nel,jb).eq.0)kva=kva+1
               if(x(nel,jb).eq.1)kvb=kvb+1
               go to 350
 340           if(x(nel,jb).eq.0)kvc=kvc+1
               if(x(nel,jb).eq.1)kvd=kvd+1
 350        continue
            lams=lams+kab(kva*kvd-kvb*kvc)
 360     continue
         if(lams.le.lama)go to 370
         jtel=kvc+kvd
         jtelz=kva+kvb
         lama=lams
         jma=j
 370  continue
      if(jnat.lt.pp)go to 375
      kwan(ka)=-kwan(ka)
      go to 400
cc
cc    splitting
cc
 375  nel=ner(ka)
      if(x(nel,jma).eq.1)then
         nzf=0
         jtel2=jtel
      else
         nzf=1
         jtel2=jtelz
      endif
      jres=kb-ka+1-jtel2
      km=ka+jtel2
      l=ka
c  -- inner loop --
 378  nel=ner(l)
      if(x(nel,jma).eq.nzf)go to 380
      l=l+1
      if(l.lt.km)go to 378
      go to 390

 380  do 381 lbb=l,kb
         nelbb=ner(lbb)
         if(x(nelbb,jma).eq.nzf)go to 381
         lcc=lbb-1
         go to 382
 381  continue
 382  do 383 laa=l,lcc
         ldd=lcc+l-laa
         lee=ldd+1
         ner(lee)=ner(ldd)
 383  continue
      ner(l)=nelbb
      go to 378

 390  nclu=nclu+1
      nban(km)=npass
      kwan(ka)=jtel2
      kwan(km)=jres
      lava(km)=jma
      ka=ka+kwan(ka)
 400  if(kb.eq.nn)go to 500
 410  ka=ka+kab(kwan(ka))
      if(ka.gt.nn)go to 500
      if(kwan(ka).lt.2)go to 410
      go to 310

 500  npass=npass+1
      do 510 ka=1,nn
         if(kwan(ka).ge.2)go to 310
 510  continue
      end
cc
cc kab(j) = |j|
cc
      integer function kab(j)
      integer j
      kab=j
      if(j.lt.0) kab=-j
      return
      end
