      SUBROUTINE TWINS(NN,JPP,X,DYS,DYS2,JDYSS,VALMD,JTMD,NDYST,JALG,
     F METHOD,KWAN,NER,BAN,COEF,MERGE)
CC
CC   THIS PROGRAM PERFORMS AGGLOMERATIVE NESTING (AGNES) USING THE
CC   GROUP AVERAGE METHOD OF SOKAL AND MICHENER (1958), AS WELL AS
CC   DIVISIVE ANALYSIS (DIANA) USING THE METHOD OF MCNAUGHTON-SMITH,
CC   WILLIAMS, DALE, AND MOCKETT (1964).
CC
CC   LIST OF FUNCTIONS AND SUBROUTINES: 
CC       MAIN UNIT
CC       FUNCTION MEET
CC       SUBROUTINE DYSTA4
CC       SUBROUTINE AVERL 
CC       SUBROUTINE SPLYT 
CC       SUBROUTINE SUPCL 
CC
CC   THE FOLLOWING VECTORS AND MATRICES MUST BE DIMENSIONED IN THE
CC   MAIN PROGRAM ONLY:
CC       KWAN(NN),NER(NN),BAN(NN)
CC       X(NN,JPP),JTMD(JPP),VALMD(JPP),DYS((NN*(NN-1))/2 + 1)
CC   WHERE: 
CC       NN = MAXIMAL NUMBER OF OBJECTS
CC       JPP = MAXIMAL NUMBER OF VARIABLES USED IN THE ANALYSIS
CC
CC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION KWAN(NN),NER(NN),BAN(NN)
      DIMENSION X(NN,JPP),JTMD(JPP),VALMD(JPP),DYS((NN*(NN-1))/2 + 1)
      DIMENSION DYS2((NN*(NN-1))/2 + 1)
      DIMENSION MERGE(NN-1,2)

      IF(JDYSS.EQ.0)GO TO 70
      JPP=1
      GO TO 100
   70 JHALT=0
      CALL DYSTA4(NN,JPP,X,DYS,NDYST,JTMD,VALMD,JHALT)
      IF(JHALT.EQ.0)GO TO 100
      JDYSS=-1
      RETURN

  100 DO 110 I=1,(NN*(NN-1))/2 + 1
         DYS2(I)=DYS(I)
  110 CONTINUE
      IF(JALG.EQ.2)GO TO 200
      CALL AVERL(NN,KWAN,NER,BAN,DYS,METHOD,MERGE)
      CALL BANAG(NN,BAN,NER,COEF)
      GO TO 300
  200 CALL SPLYT(NN,KWAN,NER,BAN,DYS,MERGE)
      CALL BANDY(NN,BAN,NER,COEF)
  300 END 
CC
CC
      SUBROUTINE DYSTA4(NN,JPP,X,DYS,NDYST,JTMD,VALMD,JHALT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(NN,JPP),DYS(1+NN*(NN-1)/2),JTMD(JPP),VALMD(JPP)
      PP=JPP
      NLK=1 
      DYS(1)=0.0
      DO 100 L=2,NN 
      LSUBT=L-1 
      DO 20 K=1,LSUBT 
      CLK=0.0 
      NLK=NLK+1 
      NPRES=0 
      DO 30 J=1,JPP 
      IF(JTMD(J).GE.0)GOTO 40 
      IF(X(L,J).EQ.VALMD(J))GOTO 30 
      IF(X(K,J).EQ.VALMD(J))GOTO 30 
   40 NPRES=NPRES+1 
      IF(NDYST.NE.1)GOTO 50 
      CLK=CLK+(X(L,J)-X(K,J))*(X(L,J)-X(K,J)) 
      GOTO 30 
   50 CLK=CLK+DABS(X(L,J)-X(K,J))
   30 CONTINUE
      RPRES=NPRES 
      IF(NPRES.NE.0)GOTO 60 
      JHALT=1 
      DYS(NLK)=-1.0
      GOTO 20 
   60 IF(NDYST.NE.1)GOTO 70 
      DYS(NLK)=DSQRT(CLK*(PP/RPRES)) 
      GOTO 20 
   70 DYS(NLK)=CLK*(PP/RPRES) 
   20 CONTINUE
  100 CONTINUE
      END 
CC
CC
      SUBROUTINE AVERL(NN,KWAN,NER,BAN,DYS,METHOD,MERGE) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION KWAN(NN),DYS(1+NN*(NN-1)/2),NER(NN),BAN(NN)
      DIMENSION MERGE(NN-1,2)
      NCLU=NN-1 
CC      INITIALIZATION
      DO 10 L=1,NN
      KWAN(L)=1 
      NER(L)=L
   10 CONTINUE
CC
CC    FIND CLOSEST CLUSTERS
CC
      NMERGE=1
  100 J=1 
   80 J=J+1 
      IF(KWAN(J).EQ.0)GOTO 80 
      NEJ=MEET(1,J) 
      SMALD=DYS(NEJ)*1.1+1.0
      NNS=NN-1
      DO 120 L=1,NNS  
      IF(KWAN(L).EQ.0)GO TO 120 
      LMUCH=L+1 
      DO 110 J=LMUCH,NN 
      IF(KWAN(J).EQ.0)GO TO 110 
      NLJ=MEET(L,J) 
      IF(DYS(NLJ).GT.SMALD)GO TO 110
      SMALD=DYS(NLJ)
      LA=L
      LB=J  
  110 CONTINUE
  120 CONTINUE
CC
CC    MERGE-STRUCTURE FOR PLOTTING TREE IN S-PLUS
CC
      L1=-LA
      L2=-LB
      IF(NMERGE.EQ.1)GO TO 121
      DO 122 J=1,(NMERGE-1)
      IF((MERGE(J,1).EQ.L1).OR.(MERGE(J,2).EQ.L1))L1=J
      IF((MERGE(J,1).EQ.L2).OR.(MERGE(J,2).EQ.L2))L2=J
  122 CONTINUE
  121 MERGE(NMERGE,1)=L1
      MERGE(NMERGE,2)=L2
      NMERGE=NMERGE+1
CC
CC    DETERMINE LFYRS AND LLAST
CC
      DO 200 L=1,NN 
      IF(NER(L).EQ.LA)LFYRS=L 
      IF(NER(L).EQ.LB)LLAST=L 
  200 CONTINUE
      BAN(LLAST)=SMALD
CC
CC    IF THE TWO CLUSTERS ARE NEXT TO EACH OTHER,
CC    NER MUST NOT BE CHANGED
CC
      LNEXT=LFYRS+KWAN(LA)
      IF(LNEXT.EQ.LLAST)GOTO 230
CC
CC    UPDATING NER AND BAN
CC
      LPUT=LFYRS+KWAN(LA) 
      LNUM=LLAST-LPUT 
      DO 220 L=1,LNUM 
      LKA=NER(LPUT) 
      AKB=BAN(LPUT) 
      LENDA=LLAST+KWAN(LB)-2
      LENDB=LENDA+1 
      DO 210 J=LPUT,LENDA 
      NER(J)=NER(J+1) 
      BAN(J)=BAN(J+1) 
  210 CONTINUE
      NER(LENDB)=LKA
      BAN(LENDB)=AKB
  220 CONTINUE
CC
CC    CALCULATE NEW DISSIMILARITIES
CC
  230 DO 240 LQ=1,NN
      IF(LQ.EQ.LA.OR.LQ.EQ.LB)GO TO 240 
      IF(KWAN(LQ).EQ.0)GO TO 240
      NAQ=MEET(LA,LQ) 
      NBQ=MEET(LB,LQ) 
      IF(METHOD.EQ.2)GO TO 300
      IF(METHOD.EQ.3)GO TO 310
      IF(METHOD.EQ.4)GO TO 320
      IF(METHOD.EQ.5)GO TO 330
CC   GROUP AVERAGE METHOD
      TA=KWAN(LA) 
      TB=KWAN(LB) 
      FA=TA/(TA+TB) 
      FB=TB/(TA+TB) 
      DYS(NAQ)=FA*DYS(NAQ)+FB*DYS(NBQ)
      GO TO 240
CC   SINGLE LINKAGE
  300 DNEW=DYS(NAQ)
      IF(DYS(NBQ).LT.DNEW)DNEW=DYS(NBQ)
      DYS(NAQ)=DNEW
      GO TO 240
CC   COMPLETE LINKAGE
  310 DNEW=DYS(NAQ)
      IF(DNEW.LT.DYS(NBQ))DNEW=DYS(NBQ)
      DYS(NAQ)=DNEW
      GO TO 240
CC   WARD'S METHOD
  320 TA=KWAN(LA) 
      TB=KWAN(LB) 
      TQ=KWAN(LQ)
      FA=(TA+TQ)/(TA+TB+TQ)
      FB=(TB+TQ)/(TA+TB+TQ)
      FC=-TQ/(TA+TB+TQ)
      NAB=MEET(LA,LB)
      D=FA*DYS(NAQ)*DYS(NAQ)+FB*DYS(NBQ)*DYS(NBQ)
      D=D+FC*DYS(NAB)*DYS(NAB)
      DYS(NAQ)=SQRT(D)
      GO TO 240
CC   WEIGHTED AVERAGE LINKAGE
  330 DYS(NAQ)=(DYS(NAQ)+DYS(NBQ))/2.D0
  240 CONTINUE
  250 KWAN(LA)=KWAN(LA)+KWAN(LB)
      KWAN(LB)=0
      NCLU=NCLU-1 
      IF(NCLU.GT.0)GOTO 100 
      END 
CC
CC
      SUBROUTINE BANAG(NN,BAN,NER,AC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BAN(NN),NER(NN)
      SUP=0.0 
      DO 70 K=2,NN
      IF(BAN(K).GT.SUP)SUP=BAN(K) 
   70 CONTINUE
      AC=0.0
      DO 80 K=1,NN
      KEARL=K 
      IF(K.EQ.1)KEARL=2 
      KAFTE=K+1 
      IF(K.EQ.NN)KAFTE=NN 
      SYZE=BAN(KEARL) 
      IF(BAN(KAFTE).LT.SYZE)SYZE=BAN(KAFTE) 
      AC=AC+1.0-(SYZE/SUP)
   80 CONTINUE
      RNN=NN
      AC=AC/RNN 
      END 

CC
CC
      SUBROUTINE SPLYT(NN,KWAN,NER,BAN,DYS,MERGE) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION KWAN(NN),DYS(1+NN*(NN-1)/2),NER(NN),BAN(NN)
      DIMENSION MERGE(NN-1,2)
CC
CC    INITIALIZATION
CC
      NCLU=1
      NHALF=NN*(NN-1)/2+1 
      DO 10 L=1,NN
      KWAN(L)=0 
      BAN(L)=0. 
      NER(L)=L
   10 CONTINUE
      KWAN(1)=NN
      JA=1
CC
CC    COMPUTATION OF DIAMETER OF DATA SET
CC
      CS=0.0
      K=0 
   20 K=K+1 
      IF(DYS(K).GT.CS)CS=DYS(K) 
      IF(K.LT.NHALF)GO TO 20
CC
CC    PREPARE FOR SPLITTING
CC
   30 JB=JA+KWAN(JA)-1
      JMA=JB
CC
CC    SPECIAL CASE OF A PAIR OF OBJECTS
CC
      IF(KWAN(JA).NE.2)GO TO 50 
      KWAN(JA)=1
      KWAN(JB)=1
      JAN=NER(JA) 
      JBN=NER(JB) 
      JAB=MEET(JAN,JBN) 
      BAN(JB)=DYS(JAB)
      GO TO 400 
CC
CC    FINDING FIRST OBJECT TO BE SHIFTED
CC
   50 BYGSD=-1. 
      DO 110 L=JA,JB
      LNER=NER(L) 
      SD=0. 
      DO 100 J=JA,JB
      JNER=NER(J) 
      NLJ=MEET(LNER,JNER) 
      SD=SD+DYS(NLJ)
  100 CONTINUE
      IF(SD.LE.BYGSD)GO TO 110
      BYGSD=SD
      LNDSD=L 
  110 CONTINUE
CC
CC    SHIFTING THE FIRST OBJECT
CC
      KWAN(JA)=KWAN(JA)-1 
      KWAN(JB)=1
      IF(JB.EQ.LNDSD)GO TO 115
      LCHAN=NER(LNDSD)
      LMM=JB-1
      DO 112 LMMA=LNDSD,LMM 
      LMMB=LMMA+1 
      NER(LMMA)=NER(LMMB) 
  112 CONTINUE
      NER(JB)=LCHAN 
  115 SPLYN=0.
      JMA=JB-1
CC
CC    FINDING THE NEXT OBJECT TO BE SHIFTED
CC
  120 SPLYN=SPLYN+1.
      REST=JMA-JA 
      BDYFF=-1. 
      DO 150 L=JA,JMA 
      LNER=NER(L) 
      DA=0. 
      DO 130 J=JA,JMA 
      JNER=NER(J) 
      NLJ=MEET(LNER,JNER) 
      DA=DA+DYS(NLJ)
  130 CONTINUE
      DA=DA/REST
      DB=0. 
      JMB=JMA+1 
      DO 140 J=JMB,JB 
      JNER=NER(J) 
      NLJ=MEET(LNER,JNER) 
      DB=DB+DYS(NLJ)
  140 CONTINUE
      DB=DB/SPLYN 
      DYFF=DA-DB
      IF(DYFF.LE.BDYFF)GO TO 150
      BDYFF=DYFF
      JAWAY=L 
  150 CONTINUE
      JMB=JMA+1
CC
CC    SHIFTING THE NEXT OBJECT WHEN NECESSARY
CC
      IF(BDYFF.LE.0.)GO TO 200
      IF(JMA.EQ.JAWAY)GO TO 165 
      LCHAN=NER(JAWAY)
      LMZ=JMA-1 
      DO 160 LXX=JAWAY,LMZ
      LXXP=LXX+1
      NER(LXX)=NER(LXXP)
  160 CONTINUE
      NER(JMA)=LCHAN
  165 DO 170 LXX=JMB,JB 
      LXY=LXX-1 
      IF(NER(LXY).LT.NER(LXX))GO TO 180 
      LCHAN=NER(LXY)
      NER(LXY)=NER(LXX) 
      NER(LXX)=LCHAN
  170 CONTINUE
  180 KWAN(JA)=KWAN(JA)-1 
      KWAN(JMA)=KWAN(JMB)+1 
      KWAN(JMB)=0 
      JMA=JMA-1 
      JMB=JMA+1 
      IF(JMA.NE.JA)GO TO 120
CC
CC    SWITCH THE TWO PARTS WHEN NECESSARY
CC
  200 IF(NER(JA).LT.NER(JMB))GO TO 300
      LXXA=JA 
      DO 220 LGRB=JMB,JB
      LXXA=LXXA+1 
      LCHAN=NER(LGRB) 
      DO 210 LXY=LXXA,LGRB
      LXF=LGRB-LXY+LXXA 
      LXG=LXF-1 
      NER(LXF)=NER(LXG) 
  210 CONTINUE
      NER(LXG)=LCHAN
  220 CONTINUE
      LLQ=KWAN(JMB) 
      KWAN(JMB)=0 
      JMA=JA+JB-JMA-1 
      JMB=JMA+1 
      KWAN(JMB)=KWAN(JA)
      KWAN(JA)=LLQ
CC
CC    COMPUTE LEVEL FOR BANNER
CC
  300 IF(NCLU.EQ.1)BAN(JMB)=CS
      IF(NCLU.EQ.1)GO TO 400
      CALL SUPCL(DYS,JA,JB,AREST,NN,NER) 
      BAN(JMB)=AREST
  400 NCLU=NCLU+1 
      IF(NCLU.EQ.NN)GOTO 500
CC
CC    CONTINUE SPLITTING UNTIL ALL OBJECTS ARE SEPARATED
CC
      IF(JB.EQ.NN)GO TO 430 
  420 JA=JA+KWAN(JA)
      IF(JA.GT.NN)GO TO 430 
      IF(KWAN(JA).LE.1)GO TO 420
      GO TO 30
  430 JA=1
      IF(KWAN(JA).EQ.1)GO TO 420
      GO TO 30
CC
CC    MERGE-STRUCTURE FOR PLOTTING TREE IN S-PLUS
CC
  500 DO 550 NMERGE=1,(NN-1)
      DMIN=CS
      DO 560 J=2,NN
      IF ((KWAN(J).LT.0).OR.(BAN(J).GT.DMIN))GOTO 560
      DMIN=BAN(J)
      NJ=J
  560 CONTINUE
      KWAN(NJ)=-1
      L1=-NER(NJ-1)
      L2=-NER(NJ)
      IF(NMERGE.EQ.1)GO TO 570
      DO 580 J=1,(NMERGE-1)
      IF((MERGE(J,1).EQ.L1).OR.(MERGE(J,2).EQ.L1))L1=J
      IF((MERGE(J,1).EQ.L2).OR.(MERGE(J,2).EQ.L2))L2=J
  580 CONTINUE
  570 MERGE(NMERGE,1)=L1
      MERGE(NMERGE,2)=L2
  550 CONTINUE
      END 
CC
CC  
      SUBROUTINE SUPCL(DYS,KKA,KKB,AREST,NN,NER) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION DYS(1+NN*(NN-1)/2),NER(NN) 
      KKC=KKB-1 
      AREST=0.
      DO 20 L=KKA,KKC 
      LNER=NER(L) 
      KKD=L+1 
      DO 10 J=KKD,KKB 
      JNER=NER(J) 
      MLJ=MEET(LNER,JNER) 
      IF(DYS(MLJ).GT.AREST)AREST=DYS(MLJ) 
   10 CONTINUE
   20 CONTINUE
      RETURN
      END 
CC
CC
      SUBROUTINE BANDY(NN,BAN,NER,DC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION BAN(NN),NER(NN)
      SUP=0.0 
      DO 70 K=2,NN
      IF(BAN(K).GT.SUP)SUP=BAN(K) 
   70 CONTINUE
      DC=0.0
      DO 80 K=1,NN
      KEARL=K 
      IF(K.EQ.1)KEARL=2 
      KAFTE=K+1 
      IF(K.EQ.NN)KAFTE=NN 
      SYZE=BAN(KEARL) 
      IF(BAN(KAFTE).LT.SYZE)SYZE=BAN(KAFTE) 
      DC=DC+1.0-(SYZE/SUP)
   80 CONTINUE
      RNN=NN
      DC=DC/RNN 
      END

