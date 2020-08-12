c-----electric static deflact panle      
      SUBROUTINE ELDROI       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RR(3),POT(3,3)
      COMMON/S0/R(5,16000),T(5,16000),R1(5,16000),T1(5,16000),DIST,K,L  
      COMMON/S4/APOT,BPOT     
      COMMON VMAX,GAP,SIV,ER0(4),ET0(4),ER1(4),ET1(4),
     1DERDR(4),DETDT(4),BS0(4),DBSDR(4),DBSDT(4),BBM,DBM,       
     2RO0,TETA0,RRM,AAM,AMAX,WIDTH,POTINT,POTEX,EN,   
     3DERDT(4),DETDR(4),EP,ITEST    
      DELTAR=.001   
      DELTAT=.001   
      DO 2 I=1,3    
      RI=DBLE(FLOAT(I-2))     
      TT=T(K,L)+RI*DELTAT     
      DO 2 J=1,3    
      RJ=DBLE(FLOAT(J-2))     
      RR(J)=R(K,L)+RJ*DELTAR  
      ECA=RR(J)*DCOS(TT-AAM)-RO0*DCOS(TETA0-AAM)  
      DI=RR(J)*DSIN(TT-AAM)-RO0*DSIN(TETA0-AAM) 
   
      IF(K.NE.2) GO TO 1  
      
      IF(I.NE.2) GO TO 1      
      IF(J.NE.2) GO TO 1      
      IF(ABS(ECA).GT.EP/2.) GO TO 1     
      IF(DI.GT.GAP/2.) ITEST=1
      IF(DI.LT.-GAP/2.) ITEST=2         
      DIST=DI       
1     CONTINUE      
      GGG=(2.*ECA+EP)/GAP+1.  
      HHH=(2.*ECA-EP)/GAP-1.  
      IF(GGG.GT.30.) GGG=30.  
      IF(GGG.LT.-30.) GGG=-30.
      IF(HHH.GT.30.) HHH=30.  
      IF(HHH.LT.-30.) HHH=-30.
      FFF=(DTANH(GGG)-DTANH(HHH))*.5    
      AX=DI+GAP/2.  
      POT(I,J)=(POTINT+APOT*AX+BPOT*AX*AX/2.)*FFF 
2     CONTINUE      
      ER0(K)=-(POT(2,3)-POT(2,1))/DELTAR/2.       
      ET0(K)=-(POT(3,2)-POT(1,2))/DELTAT/R(K,L)/2.
      DERDR(K)=-(POT(2,3)-2.*POT(2,2)+POT(2,1))/DELTAR/DELTAR         
      DERDT(K)=-(POT(3,3)-POT(3,1)-POT(1,3)+POT(1,1))/DELTAR/DELTAT/4.
      DETDT(K)=-(POT(3,2)-2.*POT(2,2)+POT(1,2))/DELTAT/R(K,L)/DELTAT  
      DETDR(K)=-((POT(3,3)-POT(1,3))/RR(3)-(POT(3,1)-POT(1,1))/RR(1))/
     1DELTAT/DELTAR/4.        
      ER1(K)=DERDR(K)*R1(K,L)+DERDT(K)*T1(K,L)    
      ET1(K)=DETDR(K)*R1(K,L)+DETDT(K)*T1(K,L)    
      RETURN        
      END 