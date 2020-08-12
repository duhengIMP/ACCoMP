c-----curving magnet element   
      SUBROUTINE MACOUR       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION PPP(5),VVV(5) 
      COMMON/S0/R(5,16000),T(5,16000),R1(5,16000),T1(5,16000),DIST,K,L  
      COMMON/S3/ATT,VV0,VVR,VVT,INCL    
      COMMON VMAX,GAP,SIV,ER0(4),ET0(4),ER1(4),ET1(4),
     1DERDR(4),DETDT(4),BS0(4),DBSDR(4),DBSDT(4),BBM,DBM,       
     2RO0,TETA0,RRM,AAM,AMAX,WIDTH,POTINT,POTEX,EN,   
     3DERDT(4),DETDR(4),EP,ITEST        
      TINCL=DBLE(FLOAT(INCL))*3.141593/1800.      
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
6     AAA=RO0*DSIN(AAM-TETA0)-RR*DSIN(AAM-TT)     
      BBB=-RO0*DCOS(AAM-TETA0)+RR*DCOS(AAM-TT)    
      ECA=DATAN(AAA/BBB)      
      RRMP2=RO0*RO0+RR*RR-2.*RO0*RR*DCOS(TT-TETA0)
      DDD=DSQRT(RRMP2)        
      DI=DDD-RRM    

      IF(K.NE.2) GO TO 7  
    
      IF(IS.NE.1) GO TO 7     
      IF(ABS(ECA).GT.AMAX/2.) GO TO 7   
      DIST=DI-WIDTH/2.        
      IF(DI.GT.WIDTH) ITEST=1 
      IF(DI.LT.0.) ITEST=2    
7     AECA=DABS(ECA)
      PM=DDD*DSIN(AMAX/2.-AECA)         
      IF(INCL.EQ.0) GO TO 8   
      CP=DDD*DCOS(AMAX/2.-AECA)-RRM-WIDTH/2.      
      ANG=DATAN2(CP,PM)       
      CM=DSQRT(PM*PM+CP*CP)   
      PM=CM*DCOS(ANG+TINCL)   
8     FFF=(DTANH(PM/GAP+.5)+1.)*.5      
      VVV(IS)=1.-FFF*ATT      
      EEE=BBM+(DI-WIDTH/2.)/WIDTH*DBM   
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