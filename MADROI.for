      SUBROUTINE MADROI       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION PPP(5),VVV(5) 
      COMMON/S0/R(5,16000),T(5,16000),R1(5,16000),T1(5,16000),DIST,K,L  
      COMMON/S3/ATT,VV0,VVR,VVT,INCL    
      COMMON/S5/DIA,GRAD,POLAR
      COMMON VMAX,GAP,SIV,ER0(4),ET0(4),ER1(4),ET1(4),
     1DERDR(4),DETDT(4),BS0(4),DBSDR(4),DBSDT(4),BBM,DBM,       
     2RO0,TETA0,RRM,AAM,AMAX,WIDTH,POTINT,POTEX,EN,   
     3DERDT(4),DETDR(4),EP,ITEST       
      IS=0
11    IS=IS+1       
      GO TO (1,2,3,4,5,10),IS 
1     RR=R(K,L)     
      TT=T(K,L)     
      GO TO 6       
2     RR=R(K,L)+.001
      GO TO 6       
3     RR=R(K,L)-.001
      GO TO 6       
4     RR=R(K,L)     
      TT=T(K,L)+.001
      GO TO 6       
5     TT=T(K,L)-.001
6     ECA=RR*DCOS(TT-AAM)-RO0*DCOS(TETA0-AAM)     
      DI=RR*DSIN(TT-AAM)-RO0*DSIN(TETA0-AAM)

      IF(K.NE.2) GO TO 7
     
      IF(IS.NE.1) GO TO 7     
      IF(ABS(ECA).GT.EP/2.) GO TO 7     
      DIST=DI       
      IF(DI.GT.DIA/2.) ITEST=1
      IF(DI.LT.-DIA/2.) ITEST=2         
7     GGG=(2.*ECA+EP)/DIA/2.+1.         
      HHH=(2.*ECA-EP)/DIA/2.-1.         
      IF(GGG.GT.30.) GGG=30.  
      IF(GGG.LT.-30.) GGG=-30.
      IF(HHH.GT.30.) HHH=30.  
      IF(HHH.LT.-30.) HHH=-30.
      FFF=(DTANH(GGG)-DTANH(HHH))*.5    
      VVV(IS)=1.-FFF*ATT      
      EEE=GRAD*DI*POLAR+BBM   
      PPP(IS)=EEE*FFF         
      GO TO 11      
10    BS0(K)=PPP(1) 
      DBSDR(K)=(PPP(2)-PPP(3))*500.     
      DBSDT(K)=(PPP(4)-PPP(5))*500.     
      VV0=VVV(1)    
      VVR=(VVV(2)-VVV(3))*500.
      VVT=(VVV(4)-VVV(5))*500.
      RETURN        
      END 