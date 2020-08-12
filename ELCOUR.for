      SUBROUTINE ELCOUR       
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RR(3),POT(3,3)
      COMMON/S0/R(5,16000),T(5,16000),R1(5,16000),T1(5,16000),DIST,K,L  
      COMMON VMAX,GAP,SIV,ER0(4),ET0(4),ER1(4),ET1(4),
     1DERDR(4),DETDT(4),BS0(4),DBSDR(4),DBSDT(4),BBM,DBM,     
     2RO0,TETA0,RRM,AAM,AMAX,WIDTH,POTINT,POTEX,EN,   
     3DERDT(4),DETDR(4),EP,ITEST     
      BPOT=-EN*(POTEX-POTINT)/GAP/(RRM+GAP/2.)    
      APOT=(POTEX-POTINT)/GAP-BPOT*GAP/2.         
      DELTAR=.001   
      DELTAT=.001   
      DO 2 I=1,3    
      RI=DBLE(FLOAT(I-2))     
      TT=T(K,L)+RI*DELTAT     
      DO 2 J=1,3    
      RJ=DBLE(FLOAT(J-2))     
      RR(J)=R(K,L)+RJ*DELTAR  
      AAA=RO0*DSIN(AAM-TETA0)-RR(J)*DSIN(AAM-TT)  
      BBB=-RO0*DCOS(AAM-TETA0)+RR(J)*DCOS(AAM-TT) 
      ECA=DATAN(AAA/BBB)      
      RRMP2=RO0*RO0+RR(J)*RR(J)-2.*RO0*RR(J)*DCOS(TT-TETA0) 
      DDD=DSQRT(RRMP2)        
      DI=DDD-RRM   
      
      IF(K.NE.2) GO TO 1   

      IF(I.NE.2) GO TO 1      
      IF(J.NE.2) GO TO 1      
      IF(ABS(ECA).GT.AMAX/2.) GO TO 1   
      IF(DI.GT.GAP) ITEST=1   
      IF(DI.LT.0.) ITEST=2    
      DIST=DI-GAP/2.
1     CONTINUE      
      GGG=(2.*ECA+AMAX)*RRM/GAP+1.      
      HHH=(2.*ECA-AMAX)*RRM/GAP+1.      
      IF(GGG.GT.30.) GGG=30.  
      IF(GGG.LT.-30.) GGG=-30.
      IF(HHH.GT.30.) HHH=30.  
      IF(HHH.LT.-30.) HHH=-30.
      FFF=(DTANH(GGG)-DTANH(HHH))*.5    
      POT(I,J)=(POTINT+APOT*DI+BPOT*DI*DI/2.)*FFF 
2     CONTINUE      
      ER0(K)=-(POT(2,3)-POT(2,1))/DELTAR/2.       
      ET0(K)=-(POT(3,2)-POT(1,2))/DELTAT/R(K,L)/2.
      DERDR(K)=-(POT(2,3)-2.*POT(2,2)+POT(2,1))/DELTAR/DELTAR         
      DETDT(K)=-(POT(3,2)-2.*POT(2,2)+POT(1,2))/DELTAT/R(K,L)/DELTAT  
      DERDT(K)=-(POT(3,3)-POT(3,1)-POT(1,3)+POT(1,1))/DELTAR/DELTAT/4.
      DETDR(K)=-((POT(3,3)-POT(1,3))/RR(3)-(POT(3,1)-POT(1,1))/RR(1))/
     1DELTAT/DELTAR/4.        
      ER1(K)=DERDR(K)*R1(K,L)+DERDT(K)*T1(K,L)    
      ET1(K)=DETDR(K)*R1(K,L)+DETDT(K)*T1(K,L)    
      RETURN        
      END 