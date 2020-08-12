       PROGRAM ACCMP
C     20200730--DUHENG ---MulitParticle extend 
C     THE PROGRAM 'ACC' IS USED TO CALCULATE THE BEAM TRAJECTORIES FROM A POINT
C     BEFORE INJECTION SYSTEM TO A POINT AFTER EXTRECTION SYSTEM
c      initial version N particle max =16
c      current version N particle max =16000      

      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION BField(720,250),BBS(720)
      DIMENSION TIN(16000),BRho(16000),V00(16000),V0(16000)
      DIMENSION VF(16000),GAMMA0(16000) 
      DIMENSION R2(4),T2(4),R3(4),T3(4)
      DIMENSION R0(4),T0(4),CB0(4),DBDR(4),DBDT(4),CB1(4)  
      DIMENSION QSM(4),VI2(4),C8(4),C9(4)        
      DIMENSION A11(4),A12(4),A13(4),A14(4),A31(4),A32(4),A33(4),A34(4)
      DIMENSION VSR(4),VSRP(4)         
      DIMENSION DTETA(4),DCOEF(4) 
      DIMENSION DT(4),D(4,24),Z(5,24)
      DIMENSION RAND_DATA(6,16000)
        
      COMMON/S0/R(5,16000),T(5,16000),R1(5,16000),T1(5,16000),DIST,K,L
      COMMON/S1/TOUV,PHI(5,16000),PHI0(16000),CHF,PHIH,CHARM,DER0(4),DET0(4),DRDR(4),DRDT(4),DTDR(4),DTDT(4),DER1(4),
     &DET1(4),DGAP,DWIDTH,IHARM      
      COMMON/S1A/OME         
      COMMON/S3/ATT,VV0,VVR,VVT,INCL   
      COMMON/S4/APOT,BPOT    
      COMMON/S5/DIA,GRAD,POLAR         
      COMMON/S6/BUMP0(4),BUMPR(4),BUMPT(4),BRO0,BTETA0,BRRM,BAMAX,      BWIDTH,BGAP,BBBM,IBUMP 
      COMMON/S9/FIPRI,FIPUN,PHID,IPUN,NC,V000,VFF,VSR0,VSRF,VG0,VGF,VG02,VGF2,G0,G02,GF,GF2,PMatrix(6,6),S(2,2),IMAT        
      COMMON VMAX,GAP,SIV,ER0(4),ET0(4),ER1(4),ET1(4),         DERDR(4),DETDT(4),BS0(4),DBSDR(4),DBSDT(4),BBM,DBM,     
     *RO0,TETA0,RRM,AAM,AMAX,WIDTH,POTINT,POTEX,EN,  DERDT(4),DETDR(4),EP,ITEST        
C -----------------输入文件!!!!!!!!!!@@@@@@@@@@@
C------------------in.dat
      OPEN(UNIT=13,NAME='in.dat',TYPE='OLD',ERR=8302) 
      OPEN(UNIT=11,NAME='fort_mat.dat',TYPE='UNKNOWN',ERR=8304)

      READ(13,*) OR,OT,PASR,NPAST,NPASR
      READ(13,*) NIMP,IXY,NParticle,TFF,EPSPAS,DTTT,ICAR,IPUN,IMAT
c---- initial version      
c      READ(13,*) (R(1,I),T(1,I),TIN(I),BRho(I),PHI0(I),I=1,NParticle)
c----current version 20200731   
      I =1
      CALL RANDOM_SEED()
      CALL RANDOM_NUMBER(RAND_DATA)
      READ(13,*) R(1,I),T(1,I),TIN(I),BRho(I),PHI0(I)
      READ(13,*)Rsigma,Tsigma,TINsigma,BRhoSigma,PHI0sigma
      DO I=2,Nparticle
c----uniform distribution     
      R(1,I)=R(1,1)+Rsigma/float(Nparticle)*float(I-Nparticle/2)
      T(1,I)=T(1,1)+Tsigma/float(Nparticle)*float(I-Nparticle/2)
      TIN(I)=TIN(1)+TINsigma/float(Nparticle)*float(I-Nparticle/2)
      BRho(I)=BRho(1)+BRhoSigma/float(Nparticle)*float(I-Nparticle/2)
      PHI0(I)=PHI0(1)+PHI0sigma/float(Nparticle)*float(I-Nparticle/2)
c----Gaussian distribution
C          R(1,I)=R(1,1)+Rsigma*(RAND_DATA(2,I)-0.5)*2.0
C          T(1,I)=T(1,1)+Tsigma*(RAND_DATA(2,I)-0.5)*2.0
C          TIN(I)=TIN(1)+TINsigma*(RAND_DATA(3,I)-0.5)*2.0
C          BRho(I)=BRho(1)+BRhosigma*(RAND_DATA(4,I)-0.5)*2.0
C          PHI0(I)=PHI0(1)+PHI0sigma*(RAND_DATA(5,I)-0.5)*2.0
          
      END DO
C----------------------------------------20200731      
      
      
      
      READ(13,*) RF,TF,TEXF,PHIF,ITER
      READ(13,*) NQ,AM,IV,COEF
      READ(13,*) BRO0,BTETA0,BRRM,BGAP,BAMAX,BBBM,BWIDTH,NBUMP
      READ(13,*) FR,CHF,DWIDTH,DGAP,OUV
      READ(13,*) (DTETA(I),I=1,4)
      READ(13,*) (DCOEF(I),I=1,4)
      READ(13,*) DSTART,DSTOP,NDEE,FIPRI,FIPUN
      
      OPEN(UNIT=9,NAME='out.dat',TYPE='UNKNOWN')
      OPEN(UNIT=29,NAME='Trajector.dat',TYPE='UNKNOWN')
      OPEN(UNIT=28,NAME='1Particle.dat',TYPE='UNKNOWN')
      
      
      write(9,4000)   
4000  FORMAT(1X,'INTEGRATION BY RK'//)  
      write(9,10) OR,OT,PASR,NPAST,NPASR  
10    FORMAT(5X,3HOR=E13.6,5X,3HOT=E13.6,3X,5HPASR=E13.6,2X,6HNPAST=I4, 
     * 2X,6HNPASR=I4)  
c-----NIMP-output data step size
      
      write(9,20) NIMP,IXY,NParticle,TFF,EPSPAS,DTTT,ICAR,IPUN,IMAT 
20    FORMAT(3X,5HNIMP=I5,10X,4HIXY=I5,10X,10HNParticle=I5/4X,4HTFF=  
     * E13.6,1X,7HEPSPAS=E13.6,3X,5HDTTT=E13.6,3X,5HICAR=I5,9X,   
     * 5HIPUN=I5,9X,5HIMAT=I5)    
cccc      
      WRITE(9,*)'SIG=',Rsigma,Tsigma,TINsigma,BRhosigma,PHI0sigma
cccc      
      write(9,30) (I,R(1,I),T(1,I),TIN(I),BRho(I),PHI0(I),I=1,NParticle)
30    FORMAT(1X,'NUMBER OF PARTICLES=',I5,3X,2HR=E13.6,6X,2HT=E13.6,  
     * 4X,4HTIN=E13.6,5X,5HBRho=E13.6,3X,5HPHI0=E13.6)    
      write(9,40) RF,TF,TEXF,PHIF,ITER    
40    FORMAT(5X,3HRF=E13.6,5X,3HTF=E13.6,3X,5HTEXF=E13.6,  
     * 3X,5HPHIF=E13.6,3X,5HITER=I4)  
      write(9,50) NQ,AM,IV,COEF 
50    FORMAT(5X,3HNQ=I5,11X,3HAM=E13.6,5X,3HIV=I5,9X,5HCOEF=E13.6)   
      write(9,60) BRO0,BTETA0,BRRM,BGAP,BAMAX,BBBM,BWIDTH,NBUMP         
60    FORMAT(3X,5HBRO0=E13.6,1X,7HBTETA0=E13.6,3X,5HBRRM=E13.6,      
     * 3X,5HBGAP=E13.6,2X,6HBAMAX=E13.6/3X,5HBBBM=E13.6,    
     * 1X,7HBWIDTH=E13.6,2X,6HNBUMP=I4)   
      write(9,80) FR,CHF,DWIDTH,DGAP,OUV,DTETA(1),DTETA(2),DTETA(3), 
     * DTETA(4),DCOEF(1),DCOEF(2),DCOEF(3),DCOEF(4),DSTART,DSTOP,NDEE
80    FORMAT(2X,13HFREQUENCE FR=E13.6,4X,'VOLTAGE AMPLITUDE HF',       
     & 5H CHF=E13.6/4X,'HORIZONTAL GAP=',E13.6,4X,'VERTICAL GAP=',E13.6,
     & 4X,16HANGLE DE CALCUL=E13.6 /        
     & 5X,'DEE AXIS 1=',E13.6,5X,'DEE AXIS 2=',E13.6,5X,'DEE AXIS 3=',
     &E13.6,5X,'DEE AXIS 4=',E13.6/4X,'DEE COEF 1=',E13.6,4X,      
     &'DEE COEF 2=',E13.6,4X,'DEE COEF 3=',E13.6,4X,'DEE COEF 4=',
     &E13.6/9X,6HSTART=E13.6,9X,5HSTOP=E13.6/4X,5HNDEE=I4)   
      write(9,90) FIPRI,FIPUN    
90    FORMAT(2X,6HFIPRI=E13.6,2X,6HFIPUN=E13.6)   
      IHARM=1 
      CHARM=0.
      PHIH=0.
      LL=0
      NPAST1=720       
      PI=3.141593   
      CONVER=PI/180.
      FIPRI=FIPRI*CONVER      
      FIPUN=FIPUN*CONVER      
      RNPAST=DBLE(FLOAT(NPAST))         
      PAST=2.*PI/RNPAST

c      OPEN(UNIT=2,NAME='Ar22M',STATUS='OLD',RECORDSIZE=6000,
c     &FORM='UNFORMATTED',ACCESS='DIRECT',ERR=8301)
      OPEN(UNIT=12,NAME='B_FIELD.DAT',STATUS='OLD',ERR=8301)
      DO 810 J=1,NPASR      
      READ(12,*) (BBS(I),I=1,NPAST1)
      DO 810 I=1,NPAST1    
810   BField(I,J)=DBLE(BBS(I))*COEF
  
      CLOSE(UNIT=12)
      NPAS1=40
      NPAS2=153
      write(9,*) NPAS1,NPAS2     
      write(9,811) NPAS1         
811   FORMAT(5X,6HFOR J=I5)  
      write(9,1830) (BField(I,NPAS1),I=1,NPAST1) 
      write(9,812) NPAS2         
812   FORMAT(5X,6HFOR J=I5)  
      write(9,1830) (BField(I,NPAS2),I=1,NPAST1) 
      write(9,818)     
818   FORMAT(5X,'SECTOR AXIS'/)         
      write(9,1830) (BField(90,J),J=1,NPASR)      
1830  FORMAT(T1,720F8.4)
1831  FORMAT(T1,250F8.4)    

      IQ=IABS(NPAST)
      SIV=DBLE(FLOAT(IV))     
      IPT=ISIGN(1,NPAST)      
      SIPT=DBLE(FLOAT(IPT))   
      EPS=1.D-9     
      OT=OT*CONVER  
      DO 7040 I=1,NParticle         
      PHI0(I)=PHI0(I)*CONVER  
      PHI(1,I)=PHI0(I)        
7040  T(1,I)=T(1,I)*CONVER    
      TF=TF*CONVER  
      TFF=TFF*CONVER
      PHIH=PHIH*CONVER        
      PHIF=PHIF*CONVER        
      BTETA0=BTETA0*CONVER    
      BAMAX=BAMAX*CONVER      
      DSTART=DSTART*CONVER    
      DSTOP=DSTOP*CONVER      
      OUV=OUV*CONVER
      DTETA(1)=DTETA(1)*CONVER
      DTETA(2)=DTETA(2)*CONVER
      DTETA(3)=DTETA(3)*CONVER
      DTETA(4)=DTETA(4)*CONVER
      OME=2.*PI*FR  
      EXF=DATAN(TEXF)         
      CL=.299792458      
      AMU=.93149432   
      QN=DBLE(FLOAT(NQ))      
      QSM0=CL*CL*QN/AM/AMU    
      AA=AM*AMU/QN  
      E0=AMU*1000.  
      DO 7050 I=1,NParticle         
      AB=BRho(I)*CL   
      QSM1=CL*CL/DSQRT(AA*AA+AB*AB)     
c-----V00-- particle velosity      
      V00(I)=BRho(I)*QSM1       
      GAMMA2=1./(1.-V00(I)*V00(I)/(CL*CL))  
      IF(I.EQ.1) G02=GAMMA2   
      GAMMA0(I)=DSQRT(GAMMA2)     
      R1(1,I)=V00(I)*TIN(I)/DSQRT(1.+TIN(I)*TIN(I))         
      T1(1,I)=SIV*V00(I)/R(1,I)/DSQRT(1.+TIN(I)*TIN(I))     
      V0(I)=V00(I)  
7050  CONTINUE
      V000=V00(1)   
      G0=GAMMA0(1)  
      write(9,8500)    
8500  FORMAT(/12X,'RADIAL MATRIX  A  AP  B  BP'/)      
      DO 851 J=1,4  
      READ(13,*) Z(1,8+J),Z(1,12+J),Z(1,16+J),Z(1,20+J)      
      write(9,8511) Z(1,8+J),Z(1,12+J),Z(1,16+J),Z(1,20+J)     
851   CONTINUE      
8511  FORMAT(4E14.6)
      write(9,8512)    
8512  FORMAT(/16X,'AXIAL MATRIX  Z  ZP'/)     
      READ(13,*) Z(1,5),Z(1,7),Z(1,6),Z(1,8)       
      write(9,8513) Z(1,5),Z(1,7),Z(1,6),Z(1,8)      
8513  FORMAT(14X,2E14.6)      
      write(9,9000)    
      CLOSE(UNIT=13)
      NC=0
      ITEST=0       
      DVT=DTTT/5.       
      IBUMP=0 
      IDEE=0        
      DO 871 K=1,4  
        BUMP0(K)=0.       
        BUMPR(K)=0.       
        BUMPT(K)=0.     
        DER0(K)=0.    
        DET0(K)=0.    
        DRDR(K)=0.    
        DRDT(K)=0.    
        DTDR(K)=0.    
        DTDT(K)=0.    
        DER1(K)=0.    
        DET1(K)=0.    
871   CONTINUE      
      BORNEA=(270.+19.)*CONVER
      BORNEB=(270.+45.)*CONVER
C------------------------- other input file  
C------------------------may be Ar-2.dat inject element
      OPEN(UNIT=14,NAME='INOUT_ELEMENT.dat',TYPE='OLD',ERR=8303)
      N=0 
101   N=N+1         
      DTT=DVT       
      NPR=0         
      DIST=99.      
      VV0=1.        
      VVR=0.        
      VVT=0.        
      IF(N.NE.1) THEN   
        DO L=1,NParticle         
          V0(L)=VF(L)   
          R(1,L)=R(5,L) 
          T(1,L)=T(5,L) 
          R1(1,L)=R1(5,L)         
          T1(1,L)=T1(5,L)         
          PHI(1,L)=PHI(5,L)       
        END DO    
        DO I=5,24
          Z(1,I)=Z(5,I) 
        END DO
        IF(NBFIN.EQ.1) TEST1=SIV*(T(1,1)-BFIN)      
        IF(NBFIN.EQ.2) TEST1=R(1,1)-BFIN  
        IF(NBFIN.EQ.4) TEST1=BFIN-R(1,1)  
        IF(TEST1.LT.0.) GO TO 120 
      END IF
      READ(14,*) BFIN,NBFIN,IP        
      IF(ICAR.EQ.0) GO TO 1003
      write(9,1002) BFIN,NBFIN,IP
1002  FORMAT(//2X,'CHARACTERISTICS OF THE COMPONENT'/4X,5HBFIN=E13.6,   
     &4X,6HNBFIN=I3,5X,3HIP=I3)         
1003  CONTINUE      
      IF(NBFIN.EQ.1) BFIN=BFIN*CONVER   
      DO K=1,4  
        ER0(K)=0.     
        ER1(K)=0.     
        ET0(K)=0.     
        ET1(K)=0.     
        DERDR(K)=0.   
        DETDT(K)=0.   
        DERDT(K)=0.   
        DETDR(K)=0.   
        BS0(K)=0.     
        DBSDR(K)=0.   
        DBSDT(K)=0.   
      END DO      
      GO TO (111,112,113,114,115,116,117),IP      
111   GO TO 116         
112   READ(14,*) RO0,TETA0,RRM,AAM,GAP,AMAX        
      IF(ICAR.EQ.0) GO TO 1124
      write(9,1121) RO0,TETA0,RRM,AAM,GAP,AMAX       
1121  FORMAT(2X,'CURVED ELECTROSTATIC DEVICE'/4X,4HRO0=E13.6,2X,     
     &6HTETA0=E13.6,4X,4HRRM=E13.6,4X,4HAAM=E13.6,4X,4HGAP=E13.6,3X,  
     &5HAMAX=E13.6) 
1124  CONTINUE      
      READ(14,*) POTINT,POTEX,EN         
      IF(ICAR.EQ.0) GO TO 1125
      write(9,1123) POTINT,POTEX,EN        
1123  FORMAT(1X,7HPOTINT=E13.6,2X,6HPOTEX=E13.6,5X,3HEN=E13.6)        
1125  CONTINUE      
      GO TO 1181    
113   READ(14,*) RO0,TETA0,RRM,AAM,GAP,AMAX      
      IF(ICAR.EQ.0) GO TO 1134
      write(9,1131) RO0,TETA0,RRM,AAM,GAP,AMAX       
1131  FORMAT(2X,'CURVED MAGNETIC DEVICE'/4X,4HRO0=E13.6,2X,     
     &6HTETA0=E13.6,4X,4HRRM=E13.6,4X,4HAAM=E13.6,4X,4HGAP=E13.6,3X,  
     &5HAMAX=E13.6) 
1134  CONTINUE      
      READ(14,*) BBM,DBM,WIDTH,ATT,INCL         
      IF(ICAR.EQ.0) GO TO 1135
      write(9,1133) BBM,DBM,WIDTH,ATT,INCL 
1133  FORMAT(4X,4HBBM=E13.6,4X,4HDBM=E13.6,2X,6HWIDTH=E13.6,4X,       
     &4HATT=E13.6,3X,5HINCL=I5)         
1135  CONTINUE      
      GO TO 1181    
114   READ(14,*) RO0,TETA0,AAM,GAP,EP    
      READ(14,*) POTINT,APOT,BPOT
      IF(ICAR.EQ.0) GO TO 1143
      write(9,1142) RO0,TETA0,AAM,GAP,EP,POTINT,APOT,BPOT      
1142  FORMAT(2X,'STRAIGHT ELECTROSTATIC DEVICE'/4X,4HRO0=E13.6,         
     &2X,6HTETA0=E13.6,4X,4HAAM=E13.6,4X,4HGAP=E13.6,5X,3HEP=E13.6/   
     &1X,7HPOTINT=E13.6,3X,5HAPOT=E13.6,3X,5HBPOT=E13.6)    
1143  CONTINUE      
      GO TO 118     
115   READ(14,*) RO0,TETA0,AAM,DIA,EP    
      READ(14,*) GRAD,POLAR,BBM,ATT
      IF(ICAR.EQ.0) GO TO 1153
      write(9,1152) RO0,TETA0,AAM,DIA,EP,GRAD,POLAR,BBM,ATT    
1152  FORMAT(2X,'STRAIGHT MAGNET WITH GRADIENT'/4X,4HRO0=E13.6,2X,
     &6HTETA0=E13.6,4X,4HAAM=E13.6,4X,4HDIA=E13.6,5X,3HEP=E13.6/5X,
     &5HGRAD=E13.6,2X,6HPOLAR=E13.6,4X,4HBBM=E13.6,4X,4HATT=E13.6)       
1153  CONTINUE      
      GO TO 118     
117   READ(14,*) RO0,TETA0,RRM,AAM,GAP,AMAX        
      READ(14,*) BBM,WIDTH
      IF(ICAR.EQ.0) GO TO 1173
      write(9,1172) RO0,TETA0,RRM,AAM,GAP,AMAX,BBM,WIDTH       
1172  FORMAT(2X,'PERTURBATION FIELD'/4X,4HRO0=E13.6,2X,6HTETA0=E13.6,    
     &4X,4HRRM=E13.6,4X,4HAAM=E13.6,4X,4HGAP=E13.6,3X,5HAMAX=E13.6/   
     &1X,4HBBM=E13.6,2X,6HWIDTH=E13.6)  
1173  CONTINUE      
      GO TO 1181    
116   IF(ICAR.EQ.0) GO TO 119 
      write(9,1161)    
1161  FORMAT(2X,'ONLY HAVING PRINCIPAL FIELD AND ACCELERATING FIELD')       
      GO TO 119     
1181  AMAX=AMAX*CONVER        
118   TETA0=TETA0*CONVER      
      AAM=AAM*CONVER
119   CONTINUE      
      IF(N.NE.1)  GO TO 120   
      write(9,1195)  
      WRITE(28,1195)
1195  FORMAT(/6X,4HR(N),12X,4HT(N),10X,3HTEX,12X,4HDIST,10X,1HN,12X,
     14HPHID,9X,2HEC,10X,6HS(1,1),10X,6HS(1,2),6X,6HS(2,1),8X,6HS(2,2),    
     212X,4HBUMP,10X,3HDTT)        
      TD=T(1,1)*180./PI
      TEX=R1(1,1)/R(1,1)/T1(1,1)*DSIGN(1D0,T1(1,1))
      PHID=PHI(1,1)*180./PI   
      NN=N-1        
      GAMMA=1./DSQRT(1.-V0(1)*V0(1)/CL/CL)        
      EC=E0*(GAMMA-1)         
      write(9,1198) R(1,1),TD,TEX,NN,PHID,EC   
C      WRITE(28,1198) R(1,1),TD,TEX,NN,PHID,EC
1198  FORMAT(1X,F13.6,1X,F13.6,1X,E13.6,2X,I6,2X,F13.6,2X,E13.6/)
C1198  FORMAT(1X,F7.5,1X,F10.3,1X,E11.5,9X,I5,F10.3,2X,E13.6/)  
120   TN=T(1,1)+DTT*T1(1,1)*.5
      RN=R(1,1)+DTT*R1(1,1)*.5
      DIF=TN
      SDIF=DSIGN(1D0,DIF)      
      I=DIF/PAST+(1.+SDIF*SIPT)/2.      
      TI=DBLE(FLOAT(I-1))     
      TTI=TI*PAST
      ECT=(TN-TTI)/PAST       
      IF(ECT.GT..5) I=I+1     
      TI=DBLE(FLOAT(I-1))     
      J=(RN-OR)/PASR+1.       
      RJ=DBLE(FLOAT(J-1))     
      RRJ=RJ*PASR+OR
      ECR=(RN-RRJ)/PASR       
      IF(ECR.GT..5) J=J+1     
      RJ=DBLE(FLOAT(J-1))     
      TT0=TI*PAST
      RR0=RJ*PASR+OR
      IF(I.GE.(IQ+1)) I=MOD(I,IQ)       
      IF(I.LT.1)   I=IQ+MOD(I,IQ)       
      IF(I.EQ.0)   I=IQ       
      IPL=I+1       
      IMN=I-1       
      IF(I.EQ.1) IMN=2        
      IF(I.EQ.(IQ+1)) IPL=I-1
1281  B0=BField(I,J)
      B1=(BField(I,J+1)-BField(I,J-1))/PASR*.5    
      B2=(BField(IPL,J)-BField(IMN,J))/PAST*.5    
      B3=(BField(I,J+1)+BField(I,J-1)-2.*BField(I,J))/PASR/PASR  
      B4=(BField(IPL,J)+BField(IMN,J)-2.*BField(I,J))/PAST/PAST  
      B5=(BField(IPL,J+1)-BField(IPL,J-1)-BField(IMN,J+1)+
     &BField(IMN,J-1))/PASR/PAST*.25
      IF(NBUMP.EQ.0) GO TO 129
      IBUMPP=IBUMP
      IBUMP=0 
      IBUMP=0       
      DIF=DMOD(T(1,1),2.*PI)  
      IF(DIF.LT.0.) DIF=DIF+2.*PI       
      IF(DIF.GT.BORNEA.AND.DIF.LT.BORNEB) IBUMP=-1
129   CONTINUE      
      IF(NDEE.EQ.0) GO TO 1283
      IDEEE=IDEE    
      IDEE=0        
      TS1=(T(1,1)-DSTART)*SIV 
      TS2=(T(1,1)-DSTOP)*SIV  
      IF(TS1.LT.0.) GO TO 1283
      IF(TS2.GT.0.) GO TO 1283
      ARG2=2.*PI    
      DIF=DMOD(T(1,1),ARG2)   
      IF(DIF.LT.0.) DIF=DIF+ARG2        
      IF(ABS(DIF-DTETA(1)).LT.OUV) IDEE=1         
      IF(ABS(DIF-DTETA(2)).LT.OUV) IDEE=2         
      IF(ABS(DIF-DTETA(3)).LT.OUV) IDEE=3         
      IF(ABS(DIF-DTETA(4)).LT.OUV) IDEE=4         
      IF(IDEE.EQ.0) GO TO 1283
      TOUV=DTETA(IDEE)        
      VMAX=DCOEF(IDEE)        
1283  CONTINUE      
      DO 350 L=1,NParticle
        IF(L.NE.1) DTT=DVT      
        MMM=4         
        IF(L.EQ.1.AND.IMAT.EQ.0) MMM=8    
        IF(L.EQ.1.AND.IMAT.NE.0) MMM=24   
        Z(1,1)=R(1,L) 
        Z(1,2)=R1(1,L)
        Z(1,3)=T(1,L) 
        Z(1,4)=T1(1,L)
130     K=0 
131     K=K+1         
133     R0(K)=R(K,L)-RR0        
        T0(K)=T(K,L)-TT0        
        IF(NBUMP.EQ.0) GO TO 135    
        IF(IBUMP.NE.0) GO TO 134    
        IF(IBUMPP.EQ.0) GO TO 135   
        BUMP0(K)=0.       
        BUMPR(K)=0.       
        BUMPT(K)=0.       
        GO TO 135         
134     CALL BUMP         
135     CONTINUE
        IF(IDEE.EQ.0.AND.IDEEE.EQ.0) GO TO 140
        ER0(K)=0.         
        ET0(K)=0.        
        DERDR(K)=0.       
        DERDT(K)=0.       
        DETDR(K)=0.       
        DETDT(K)=0.       
        ER1(K)=0.         
        ET1(K)=0.         
        IF(IDEE.NE.0) GO TO 139 
        DER0(K)=0.    
        DET0(K)=0.    
        DRDR(K)=0.    
        DRDT(K)=0.    
        DTDR(K)=0.    
        DTDT(K)=0.    
        DER1(K)=0.    
        DET1(K)=0.    
        GO TO 140     
139   CALL DEE      
140   CONTINUE      
      GO TO (150,142,143,144,145,150,146),IP      
142   CALL ELCOUR   
      GO TO 147     
143   CALL MACOUR   
      GO TO 147     
144   CALL ELDROI   
      GO TO 147     
145   CALL MADROI   
      GO TO 147     
146   CALL PERTUR   
      GO TO 150     
147   IF(K.EQ.2) GO TO 150    
      IF(ITEST.EQ.0) GO TO 150
      TD=T(1,L)*180./PI       
      TEX=SIV*R1(1,L)/R(1,L)/T1(1,L)*DSIGN(1D0,T1(1,L))      
      write(9,1470) R(1,L),TD,TEX,N,L,DIST 
1470  FORMAT(2X,'FOR THE VALUES   ',5HR(N)=E13.6,3X,5HT(N)=E13.6,   
     &4X,4HTEX=E13.6,5X,2HN=I5,3X,'PARTICLE NO. 'I3,3X,5HDIST=E13.6)    
      IF(ITEST.EQ.2) GO TO 148
      write(9,1471)    
1471  FORMAT(2X,'INTERNAL LIMIT IS HITTED')       
      GO TO 999     
148   write(9,1481)    
1481  FORMAT(2X,'EXTERNAL LIMIT IS HITTED')       
      GO TO 999     
150   CONTINUE      
      IF(NDEE.EQ.0.OR.IDEE.EQ.0) GO TO 1505       
      ER0(K)=ER0(K)+DER0(K)   
      ET0(K)=ET0(K)+DET0(K)   
      DERDR(K)=DERDR(K)+DRDR(K)         
      DERDT(K)=DERDT(K)+DRDT(K)         
      DETDR(K)=DETDR(K)+DTDR(K)         
      DETDT(K)=DETDT(K)+DTDT(K)         
      ER1(K)=ER1(K)+DER1(K)   
      ET1(K)=ET1(K)+DET1(K)   
1505  CONTINUE      
      CBB0=B0+B1*R0(K)+B2*T0(K)+B3*R0(K)*R0(K)*.5+
     1B4*T0(K)*T0(K)*.5+B5*R0(K)*T0(K)  
      CB0(K)=-SIV*(CBB0*VV0+BS0(K)+BUMP0(K))      
      VI2(K)=R1(K,L)*R1(K,L)+R(K,L)*R(K,L)*T1(K,L)*T1(K,L)  
      GAMMA2=1./(1.-VI2(K)/CL/CL)       
      GAMMA=DSQRT(GAMMA2)     
      QSM(K)=QSM0/GAMMA       
      WT=ER0(K)*R1(K,L)+ET0(K)*R(K,L)*T1(K,L)     
      W=WT/CL/CL    
      GR=QSM(K)*(CB0(K)*R(K,L)*T1(K,L)+ER0(K)-R1(K,L)*W)    
      GT=QSM(K)*(-CB0(K)*R1(K,L)+ET0(K)-R(K,L)*T1(K,L)*W)   
      R2(K)=R(K,L)*T1(K,L)*T1(K,L)+GR   
      T2(K)=(-2.*R1(K,L)*T1(K,L)+GT)/R(K,L)
      IF(L.NE.1) GO TO 1520   
      DW=(ER1(K)*R1(K,L)+ER0(K)*R2(K)+ET1(K)*R(K,L)*T1(K,L)+
     1ET0(K)*(R1(K,L)*T1(K,L)+R(K,L)*T2(K)))/CL/CL
      CBDR=B1+B3*R0(K)+B5*T0(K)         
      CBDT=B2+B4*T0(K)+B5*R0(K)         
      DBDR(K)=CBB0*VVR+CBDR*VV0+DBSDR(K)+BUMPR(K) 
      DBDT(K)=CBB0*VVT+CBDT*VV0+DBSDT(K)+BUMPT(K) 
      CB1(K)=-SIV*(DBDR(K)*R1(K,L)+DBDT(K)*T1(K,L))         
      GPR=QSM(K)*(-GR*W+CB1(K)*R(K,L)*T1(K,L)+CB0(K)*(R1(K,L)*T1(K,L)+
     1R(K,L)*T2(K))+ER1(K)-R2(K)*W-R1(K,L)*DW)-T1(K,L)*GT   
      GPT=QSM(K)*(-GT*W-CB1(K)*R1(K,L)-CB0(K)*R2(K)+ET1(K)-(R1(K,L)*  
     1T1(K,L)+R(K,L)*T2(K))*W-R(K,L)*T1(K,L)*DW)+T1(K,L)*GR 
      R3(K)=GPR+3.*T1(K,L)*(R(K,L)*T2(K)+R1(K,L)*T1(K,L))   
      T3(K)=(GPT-3.*R1(K,L)*T2(K)-3.*R2(K)*T1(K,L)+R(K,L)*  
     1T1(K,L)*T1(K,L)*T1(K,L))/R(K,L)   
      VDBZDB=-SIV*(-DBDR(K)*R(K,L)*T1(K,L)+DBDT(K)*R1(K,L)/R(K,L))        
      DEZDZ=-(ER0(K)/R(K,L)+DERDR(K)+DETDT(K)/R(K,L)) 
      C8(K)=QSM(K)*(VDBZDB+DEZDZ) 
      C9(K)=-W*QSM(K)         
      IF(K.NE.1) GO TO 1510   
      IF(NPR.EQ.1) GO TO 1510 
      DELTAT=DABS(R(1,1)*T3(1)/6./V0(1))
      IF(NC.EQ.1) GO TO 1510  
      DELTAR=DABS(R3(1)/6./V0(1))       
      U=DMAX1(DELTAR,DELTAT,EPSPAS)/EPSPAS        
      FAC=DSQRT(U)  
      DVT=DTTT/FAC  
      IF(N.EQ.1) DVT=DTTT/5.  
      DTT=DVT       
1510  IF(IMAT.EQ.0) GO TO 1520
      VSR(K)=(R1(K,L)*GT-R(K,L)*T1(K,L)*GR)/VI2(K)
      IF(N.NE.1) GO TO 151        
      IF(K.NE.1) GO TO 152        
      VSR0=VSR(1)   
      VG0=(GAMMA2+GAMMA)/V000 
      VG02=GAMMA2/V000  
      GO TO 152         
151   IF(K.NE.1) GO TO 152        
      IF(NC.EQ.1.OR.NPR.EQ.1) GO TO 152 
      DD=VSRA-VSR(1)    
      Z(1,10)=Z(5,10)-DD*Z(5,11)  
      Z(1,14)=Z(5,14)-DD*Z(5,15) 
      Z(1,18)=Z(5,18)-DD*Z(5,19)  
      Z(1,22)=Z(5,22)-DD*Z(5,23)  
152   CONTINUE      
      WN=-R(K,L)*T1(K,L)*ER0(K)+R1(K,L)*ET0(K)    
      AA=R1(K,L)*DERDR(K)+T1(K,L)*DERDT(K)        
      BB=R1(K,L)*DETDR(K)+T1(K,L)*DETDT(K)        
      DETDA=(R1(K,L)*AA+R(K,L)*T1(K,L)*BB-T1(K,L)*WN)/VI2(K)    
      DETDB=(-R(K,L)*T1(K,L)*AA+R1(K,L)*BB+T1(K,L)*WT)/VI2(K)   
      VDBZDA=CB1(K) 
      GTSV=(R1(K,L)*GR+R(K,L)*T1(K,L)*GT)/VI2(K)  
      GTSR=GTSV*VSR(K)        
      VSRP(K)=(R1(K,L)*GPT-R(K,L)*T1(K,L)*GPR)/VI2(K)-2.*GTSR         
      WR=WN/VI2(K)      
      A11(K)=QSM(K)*(DETDA+VSR(K)*WR)/GAMMA2
      A12(K)=-QSM(K)*W*3.         
      A13(K)=QSM(K)*DETDB/GAMMA2+VSRP(K)-A12(K)*VSR(K)
      A14(K)=QSM(K)*WR/GAMMA2+VSR(K)       
      A31(K)=QSM(K)*(-VDBZDA+DETDB-VSR(K)*W)-VSRP(K)  
      A32(K)=-QSM(K)*WR-VSR(K)*GAMMA2      
      A33(K)=QSM(K)*(-DETDA+VSR(K)*WR)-C8(K)+VSR(K)*VSR(K)*(GAMMA2-1)     
      A34(K)=-QSM(K)*W  
1520  CONTINUE      
      RK=DBLE(FLOAT(K/3))     
      DT(K)=DTT*(1.+RK)*.5    
      D(K,1)=R1(K,L)
      D(K,2)=R2(K)  
      D(K,3)=T1(K,L)
      D(K,4)=T2(K)  
      IF(L.NE.1) GO TO 153    
      D(K,5)=Z(K,6) 
      D(K,6)=C8(K)*Z(K,5)+C9(K)*Z(K,6)  
      D(K,7)=Z(K,8) 
      D(K,8)=C8(K)*Z(K,7)+C9(K)*Z(K,8)      
      IF(IMAT.EQ.0) GO TO 153 
      D(K,9)=Z(K,10)
      D(K,10)=A11(K)*Z(K,9)+A12(K)*Z(K,10)+A13(K)*Z(K,11)+  
     1A14(K)*Z(K,12)
      D(K,11)=Z(K,12)         
      D(K,12)=A31(K)*Z(K,9)+A32(K)*Z(K,10)+A33(K)*Z(K,11)+  
     1A34(K)*Z(K,12)
      D(K,13)=Z(K,14)         
      D(K,14)=A11(K)*Z(K,13)+A12(K)*Z(K,14)+A13(K)*Z(K,15)+ 
     1A14(K)*Z(K,16)
      D(K,15)=Z(K,16)         
      D(K,16)=A31(K)*Z(K,13)+A32(K)*Z(K,14)+A33(K)*Z(K,15)+ 
     1A34(K)*Z(K,16)
      D(K,17)=Z(K,18)         
      D(K,18)=A11(K)*Z(K,17)+A12(K)*Z(K,18)+A13(K)*Z(K,19)+ 
     1A14(K)*Z(K,20)
      D(K,19)=Z(K,20)         
      D(K,20)=A31(K)*Z(K,17)+A32(K)*Z(K,18)+A33(K)*Z(K,19)+ 
     1A34(K)*Z(K,20)
      D(K,21)=Z(K,22)         
      D(K,22)=A11(K)*Z(K,21)+A12(K)*Z(K,22)+A13(K)*Z(K,23)+ 
     1A14(K)*Z(K,24)
      D(K,23)=Z(K,24)         
      D(K,24)=A31(K)*Z(K,21)+A32(K)*Z(K,22)+A33(K)*Z(K,23)+ 
     1A34(K)*Z(K,24)
153   PHI(K+1,L)=PHI(1,L)+DT(K)*OME     
      DO 155 MM=1,MMM         
      IF(K.EQ.4) GO TO 154    
      Z(K+1,MM)=Z(1,MM)+D(K,MM)*DT(K)   
      GO TO 155     
154   Z(5,MM)=Z(1,MM)+DT(4)*(D(1,MM)+2.*D(2,MM)+2.*D(3,MM)+D(4,MM))/6.
155   CONTINUE      
      R(K+1,L)=Z(K+1,1)       
      R1(K+1,L)=Z(K+1,2)      
      T(K+1,L)=Z(K+1,3)       
      T1(K+1,L)=Z(K+1,4)      
      IF(K.NE.4) GO TO 131    
      IF(NC.EQ.1) GO TO 157   
      IF(L.NE.1) GO TO 1560   
      IF(TFF.GT.(361.*CONVER)) GO TO 1560  
      IF(NPR.EQ.1) GO TO 157  
c       write(9,*)'     n=',n
c       write(9,*)'r(5,1)=',r(5,1)
      ITF=T(5,L)/TFF
      IF(T(5,L).GT.0.) ITF=ITF+1        
      IF(N.EQ.1) GO TO 156    
      IF(ITF.NE.ITFPRE) NPR=1 
      RIT=DBLE(FLOAT(ITFPRE))+(SIV-1.)/2.         
156   ITFPRE=ITF    
1560  IF(N.EQ.1) GO TO 163    
157   CONTINUE      
      IF(ITER.NE.3) GO TO 158 
      IF(DABS(PHI(5,L)-PHIF).GT.1.) GO TO 159     
      SVV=DSIGN(1D0,T1(5,L))   
      EXFF=EXF-SVV*(TF-T(5,L))
      DELTA=R(5,L)*DSIN(EXFF)-RF*DSIN(EXF)        
158   IF(ITER.EQ.1) DELTA=SIV*(T(5,L)-TF)         
      IF(ITER.EQ.2) DELTA=R(5,L)-RF     
      IF(ITER.EQ.4)  DELTA=RF-R(5,L)        
      IF(NC.EQ.1) GO TO 160   
      IF(L.EQ.1.AND.DELTA.GT.0.) NC=1   
      IF(NC.EQ.1) GO TO 1580      
159   TEX=R1(5,L)/R(5,L)/T1(5,L)*DSIGN(1D0,T1(5,L))
      IF(ABS(TEX).GT..6) NPR=0
      IF(NPR.EQ.0) GO TO 163  
      DELTA=SIV*(T(5,L)-TFF*RIT)        
      IF(DABS(DELTA).LT.EPS) GO TO 163  
      DELTAP=SIV*T1(5,L)      
      GO TO 162     
1580  write(9,1581)        
1581  FORMAT(/2X,'FINISH OF TRAJECTORIES'/)
c      write(9,1195)        
160   IF(DABS(DELTA).LT.EPS) GO TO 163  
      IF(ITER.EQ.1) DELTAP=SIV*T1(5,L)  
      IF(ITER.EQ.2) DELTAP=R1(5,L)      
      IF(ITER.EQ.4)  DELTAP=-R1(5,L)        
      IF(ITER.EQ.3) DELTAP=R1(5,L)*DSIN(EXFF)+R(5,L)*SVV*T1(5,L)*     
     1DCOS(EXFF)    
162   DTT=DTT-DELTA/DELTAP    
c      write(9,*)'npr=',npr
      GO TO 130     
163   CONTINUE      
      VF2=R1(5,L)*R1(5,L)+R(5,L)*R(5,L)*T1(5,L)*T1(5,L)     
      VF(L)=DSQRT(VF2)        
      IF(L.NE.1) GO TO 164       
      IF(IMAT.EQ.0) GO TO 164     
      VSRF=VSR(1)+DT(4)*(VSRP(1)+2.*VSRP(2)+2.*VSRP(3)+VSRP(4))/6.    
      VSRA=VSRF         
164   CONTINUE      
      IF(NC.EQ.1) GO TO 201   
      IF(PHI(5,L).LT.FIPRI) GO TO 261   
      IF(NPR.EQ.1) GO TO 201  
      IPRINT=MOD(N,NIMP)      
      IF(IPRINT.NE.0) GO TO 350         
201   TD=T(5,L)*180./PI       
      TEX=R1(5,L)/R(5,L)/T1(5,L)*DSIGN(1D0,T1(5,L))
c      write(12,*) R1(5,L),T1(5,l),R(5,L),TEX,TD,PHID
      NN=N
      PHID=PHI(5,L)*180./PI   
      GAMMA2=1./(1.-VF2/CL/CL)
      GAMMA=DSQRT(GAMMA2)     
      EC=E0*(GAMMA-1)         
      IF(L.NE.1) GO TO 23     
      PMatrix(5,5)=Z(5,5) 
      PMatrix(5,6)=Z(5,7) 
      PMatrix(6,5)=Z(5,6) 
      PMatrix(6,6)=Z(5,8) 
      S(1,1)=PMatrix(5,5) 
      S(1,2)=PMatrix(5,6)*V000      
      S(2,1)=PMatrix(6,5)/VF(1)     
      S(2,2)=PMatrix(6,6)*V000/VF(1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc-clear      
c      ITD=-INT((TD+1.0)/360.)
c      write(11,*) ITD,R(5,L) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      write(9,2500) R(5,1),TD,TEX,DIST,NN,PHID,EC,   
     &((S(I,J),J=1,2),I=1,2),BUMP0(2),DVT  
      WRITE(28,2500) R(5,1),TD,TEX,DIST,NN,PHID,EC,   
     &((S(I,J),J=1,2),I=1,2),BUMP0(2),DVT
C2500  FORMAT(1X,F7.5,1X,F10.3,1X,E13.6,1X,F6.4,1X,I6,F10.3,2X,E13.6,
2500  FORMAT(1X,F13.6,1X,F15.6,1X,E15.6,1X,F13.6,1X,I8,2X,F15.6,2X,
     &E13.6,4(1X,E13.6),1X,F13.6,1X,F13.6) 
      IF(IPUN.EQ.0.OR.IPUN.EQ.3) GO TO 24
      IF(PHI(5,1).LT.FIPUN) GO TO 24
      WRITE(11,"(5(2x,F13.6))") S(1,1),S(1,2),S(2,1),S(2,2),PHID
      GO TO 24      
23    write(9,2300) R(5,L),TD,TEX,DIST,NN,PHID,EC   
      WRITE(29,2300)R(5,L),TD,TEX,DIST,NN,PHID,EC
2300  FORMAT(1X,F13.6,1X,F15.6,1X,E15.6,1X,F13.6,1X,I8,F15.6,1X,E13.6)
C2300  FORMAT(1X,F7.5,1X,F10.3,1X,E11.5,1X,F6.4,1X,I6,F10.3,1X,E13.6)    
24    CONTINUE      
      IF(IXY.NE.1) GO TO 261  
      X=R(5,L)*DCOS(T(5,L))   
      Y=R(5,L)*DSIN(T(5,L))   
      write(9,2600) X,Y,GAMMA,CB0(2),BS0(2),DBDR(2),DBDT(2)    
2600  FORMAT(10X,2HX=F13.6,2X,2HY=F13.6,2X,6HGAMMA=E13.6,     
     &2X,4HCB0=E13.6,2X,4HBS0=E13.6,1X,5HDBDR=E13.6,1X,5HDBDT=E13.6)  
      X=X*1000.0                                                                                                                                                                            !2003.10.15
      Y=Y*1000.0                                                                                                                                                                            !2003.10.15
      WRITE(10,2605) X,Y                                                                                                                                                                    !2003.10.15 
2605      FORMAT(10X,F8.2,5X,F8.2)                                                                                                                                                          !2003.10.15
 
261   IF(L.NE.1.OR.IMAT.EQ.0) GO TO 27  
      VFF=VF(1)     
      VGF=(GAMMA2+GAMMA)/VF(1)
      GF2=GAMMA2    
      GF=GAMMA      
      VGF2=GAMMA2/VF(1) 
      DO 263 I=1,4  
      DO 263 J=1,4  
      JJ=I+4*(J+1)  
      PMatrix(I,J)=Z(5,JJ)
263   CONTINUE      
      CALL MAT
      LL=LL+1
27    IF(NC.EQ.1) GO TO 300   
      GO TO 350     
300   CONTINUE      
      RIG=GAMMA*VF(L)/QSM0
C      write(12,*) RIG    
      write(9,3040) V00(L),VF(L),GAMMA0(L),GAMMA,EC,RIG,LL  
3040  FORMAT(5X,4HV00=E13.6,3X,3HVF=E13.6,3X,7HGAMMA0=E13.6,
     &3X,7HGAMMAF=E13.6,3X,3HEC=E13.6,5X,5HBRho=E13.6,3X,3HLL=I5/)     
350   CONTINUE      
      
      if (nc.eq.1) go to 990      
      IF(NPR.NE.1) GO TO 101  
      write(9,9000)    
9000  FORMAT(1X)    
      GO TO 101     
990   IF(IMAT.EQ.0) GO TO 999 
      IMAT=9        
      CALL MAT
      CLOSE(UNIT=14)
      CLOSE(UNIT=11)
      close(9)
      GO TO 999
8301  WRITE(8311)
8311  FORMAT(T1,'UNIT=12 FILE OPEN ERROR!')
      CALL EXIT
8302  WRITE(8312)
8312  FORMAT(T1,'UNIT=13 FILE OPEN ERROR!')
      CALL EXIT
8303  WRITE(8313)
8313  FORMAT(T1,'UNIT=14 FILE OPEN ERROR!')
      CALL EXIT
8304  WRITE(8314)
8314  FORMAT(T1,'UNIT=11 FILE OPEN ERROR!')
      CALL EXIT      
999   STOP
      END 
      

