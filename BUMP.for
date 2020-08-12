
      SUBROUTINE BUMP         
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON/S0/R(5,16000),T(5,16000),R1(5,16000),T1(5,16000),DIST,K,L  
      COMMON/S6/BUMP0(4),BUMPR(4),BUMPT(4),BRO0,BTETA0,BRRM,BAMAX,        
     1BWIDTH,BGAP,BBBM,IBUMP  
c      COMMON VMAX,GAP,SIV,ER0(4),ET0(4),ER1(4),ET1(4),DERDR(4),DETDT(4),BS0(4),DBSDR(4),DBSDT(4),BBM,DBM,       
c     2RO0,TETA0,RRM,AAM,AMAX,WIDTH,POTINT,POTEX,EN,   DERDT(4),DETDR(4),EP,ITEST        
      PI=3.141593   
      BAAM=112.5*PI/180.      
      IF(IBUMP.EQ.-1) BAA=BAAM+PI       
      IF(IBUMP.EQ.-1) BTE=BTETA0+PI     
      SINA=DSIN(BTE-T(K,L))   
      COSA=DCOS(BTE-T(K,L))   
      SINO=DSIN(BAA-BTE)      
      COSO=DCOS(BAA-BTE)      
      AAA=BRO0*SINO-R(K,L)*DSIN(BAA-T(K,L))       
      BBB=-BRO0*COSO+R(K,L)*DCOS(BAA-T(K,L))       
      ECA=DATAN(AAA/BBB)      
      RRMP2=AAA*AAA+BBB*BBB   
      DDD=DSQRT(RRMP2)        
      BDIST=DDD-(BRRM+BWIDTH*.5)        
      GGG=(2.*ECA+BAMAX)*BRRM/BGAP+1.   
      HHH=(2.*ECA-BAMAX)*BRRM/BGAP-1.   
      FFF=(DTANH(GGG)-DTANH(HHH))*.5    
      HH=BDIST*40.-1.   
      GG=BDIST*40.+1.
      PP=DCOSH(BDIST*10.)     
      QQ=DSINH(BDIST*10.)     
      FF=(.6*(DTANH(GG)-DTANH(HH))-.15/PP)/.764   
      SIBUMP=DBLE(FLOAT(IBUMP))         
      BUMP0(K)=BBBM*FF*FFF*SIBUMP       
      DECADR=BRO0*SINA/RRMP2  
      DECADT=R(K,L)*(R(K,L)-BRO0*COSA)/RRMP2      
      U1=DCOSH(GGG) 
      U2=DCOSH(HHH) 
      UU=BRRM*(1./U1/U1-1./U2/U2)/BGAP  
      DFFFDR=DECADR*UU        
      DFFFDT=DECADT*UU        
      DIDR=(R(K,L)-BRO0*COSA)/DDD       
      DIDT=-R(K,L)*BRO0*SINA/DDD        
      V1=DCOSH(GG)  
      V2=DCOSH(HH)  
      VV=((1./V1/V1-1./V2/V2)*24.+QQ/PP/PP*1.5)/.764        
      DFFDR=DIDR*VV 
      DFFDT=DIDT*VV 
      BUMPR(K)=BBBM*(FF*DFFFDR+FFF*DFFDR)*SIBUMP  
      BUMPT(K)=BBBM*(FF*DFFFDT+FFF*DFFDT)*SIBUMP  
      RETURN        
      END 