      SUBROUTINE MAT
c-------calculate matrix 
c-------O----matrix in R R' Theta Theta'      
c-------S----matrix in Z Z'                    
c
c--------------------------------------------      
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION O(4,4),Q(4,9),U(2,2)    
      COMMON/S1A/OME
      COMMON/S9/FIPRI,FIPUN,PHID,IPUN,NC,V000,VFF,VSR0,VSRF,     
     1VG0,VGF,VG02,VGF2,G0,G02,GF,GF2,PMatrix(6,6),S(2,2),IMAT         
      CONVER=3.141593/180.    
      PHI=PHID*CONVER         
      O(1,1)=PMatrix(3,3)+PMatrix(3,2)*VSR0         
      O(1,2)=PMatrix(3,4)*V000      
      O(1,3)=PMatrix(3,1)/OME*V000*CONVER*1000.         
      O(1,4)=PMatrix(3,2)/VG0*1000. 
      O(2,1)=(PMatrix(4,3)+PMatrix(4,2)*VSR0)/VFF   
      O(2,2)=PMatrix(4,4)*V000/VFF  
      O(2,3)=PMatrix(4,1)/OME*V000/VFF*CONVER*1000.     
      O(2,4)=PMatrix(4,2)/VG0/VFF*1000.       
      O(3,1)=(PMatrix(1,3)+PMatrix(1,2)*VSR0)*OME/VFF/CONVER/1000.      
      O(3,2)=PMatrix(1,4)*OME*V000/VFF/CONVER/1000.     
      O(3,3)=PMatrix(1,1)*V000/VFF  
      O(3,4)=PMatrix(1,2)*OME/VG0/VFF/CONVER  
      O(4,1)=(PMatrix(2,3)+PMatrix(2,2)*VSR0-(PMatrix(3,3)+PMatrix(3,2)*VSR0)*VSRF)*VGF/1000. 
      O(4,2)=(PMatrix(2,4)-PMatrix(3,4)*VSRF)*V000*VGF/1000.  
      O(4,3)=(PMatrix(2,1)-PMatrix(3,1)*VSRF)/OME*V000*VGF*CONVER       
      O(4,4)=(PMatrix(2,2)-PMatrix(3,2)*VSRF)*VGF/VG0         
      IF(IMAT.EQ.9) GO TO 1   
      IF(IPUN.LT.2) GO TO 1   
      IF(PHI.LT.FIPUN) GO TO 1
C      WRITE(11) O(1,1),O(1,2),O(2,1),O(2,2),PHID
1     IF(PHI.LT.FIPRI) GO TO 7
      IF(IMAT.GE.8) GO TO 2   
      CALL DETINV(Q,O,DETRAD,U,S,DETAXI,IMAT)     
      write(9,1004) ((O(I,J),J=1,4),I=1,2),DETAXI    
1004  FORMAT(1X,8E12.5,5X,7HDETAXI=E12.5)         
      write(9,1005) ((O(I,J),J=1,4),I=3,4),DETRAD    
1005  FORMAT(1X,8E12.5,5X,7HDETRAD=E12.5)         
      write(9,1006)    
1006  FORMAT(/)     
      GO TO 7       
2     NP=0
      IF(IMAT.EQ.9) write(9,1007)
1007  FORMAT(1H1)   
      CALL DETINV(Q,O,DETRAD,U,S,DETAXI,IMAT)     
      AXINOR=DETAXI*VFF*GF/(V000*G0)    
      RADNOR=DETRAD*VFF*(GF2-GF)/(V000*(G02-G0))  
      write(9,2000)    
2000  FORMAT(/4X,44HRADIAL  MATRIX   X  XP   DELTA PHI  DELTAW/W,     
     18X,22HAXIAL  MATRIX   Z  ZP ,12X,12HDETERMINANTS/)    
      GO TO 6       
3     O(1,3)=PMatrix(3,1) 
      O(1,4)=PMatrix(3,2)/VG02*1000.
      O(2,3)=PMatrix(4,1)/VFF       
      O(2,4)=PMatrix(4,2)/VG02/VFF*1000.      
      O(3,1)=PMatrix(1,3)+PMatrix(1,2)*VSR0         
      O(3,2)=PMatrix(1,4)*V000      
      O(3,3)=PMatrix(1,1) 
      O(3,4)=PMatrix(1,2)/VG02*1000.
      O(4,1)=(PMatrix(2,3)+PMatrix(2,2)*VSR0-(PMatrix(3,3)+PMatrix(3,2)*VSR0)*VSRF)*VGF2/1000.
      O(4,2)=(PMatrix(2,4)-PMatrix(3,4)*VSRF)*VGF2*V000/1000. 
      O(4,3)=(PMatrix(2,1)-PMatrix(3,1)*VSRF)*VGF2/1000.      
      O(4,4)=(PMatrix(2,2)-PMatrix(3,2)*VSRF)*VGF2/VG02       
      CALL DETINV(Q,O,DETRAD,U,S,DETAXI,IMAT)     
      RADNOR=DETRAD*(VGF2*VFF-1.)/(VG02*V000-1.)  
      write(9,2001)    
2001  FORMAT(/4X,44HRADIAL  MATRIX   X  XPR   DELTA L   DELTAP/P,     
     18X,22HAXIAL  MATRIX   Z  ZP ,12X,12HDETERMINANTS/)    
      GO TO 6       
4     AMU=.931478   
      CL=.2998      
      AM0=AMU/(CL*CL)         
      G03=G02*G0*V000         
      GF3=GF2*GF*VFF
      S(1,1)=PMatrix(5,5) 
      S(1,2)=PMatrix(5,6)/(AM0*G0)  
      S(2,1)=PMatrix(6,5)*AM0*GF    
      S(2,2)=PMatrix(6,6)*GF/G0     
      O(1,1)=PMatrix(3,3)+VSR0*PMatrix(3,2)         
      O(1,2)=PMatrix(3,4)/(AM0*G0)  
      O(1,3)=PMatrix(3,1)*V000/OME  
      O(1,4)=PMatrix(3,2)/(AM0*G03) 
      O(2,1)=(PMatrix(4,3)+PMatrix(4,2)*VSR0)*AM0*GF
      O(2,2)=PMatrix(4,4)*GF/G0     
      O(2,3)=PMatrix(4,1)*AM0*GF*V000/OME     
      O(2,4)=PMatrix(4,2)*GF/G03    
      O(3,1)=(PMatrix(1,3)+VSR0*PMatrix(1,2))*OME/VFF         
      O(3,2)=PMatrix(1,4)*OME/(VFF*AM0*G0)    
      O(3,3)=PMatrix(1,1)*V000/VFF  
      O(3,4)=PMatrix(1,2)*OME/(VFF*AM0*G03)   
      O(4,1)=(PMatrix(2,3)+VSR0*PMatrix(2,2)-VSRF*(PMatrix(3,3)+VSR0*PMatrix(3,2)))*AM0*GF3   
      O(4,2)=(PMatrix(2,4)-VSRF*PMatrix(3,4))*GF3/G0
      O(4,3)=(PMatrix(2,1)-VSRF*PMatrix(3,1))*AM0*GF3*V000/OME
      O(4,4)=(PMatrix(2,2)-VSRF*PMatrix(3,2))*GF3/G03         
      CALL DETINV(Q,O,DETRAD,U,S,DETAXI,IMAT)     
      AXINOR=DETAXI 
      RADNOR=DETRAD 
      write(9,2002)    
2002  FORMAT(/4X,44HRADIAL  MATRIX   X  PX   DELTA PHI   DELTA W,     
     18X,22HAXTIAL MATRIX   Z   PZ,12X,12HDETERMINANTS/)    
      GO TO 6       
5     S(1,1)=PMatrix(5,5) 
      S(1,2)=PMatrix(5,6) 
      S(2,1)=PMatrix(6,5) 
      S(2,2)=PMatrix(6,6) 
      DO 50 I=1,4   
      DO 50 J=1,4   
50    O(I,J)=PMatrix(I,J) 
      CALL DETINV(Q,O,DETRAD,U,S,DETAXI,IMAT)     
      AXINOR=DETAXI*GF/G0     
      RADNOR=DETRAD*GF2*GF2/(G02*G02)   
      write(9,2003)    
2003  FORMAT(/4X,35HRADIAL  MATRIX   A  AP     B    BP ,    
     117X,22HMATRICE AXIALE  Z  ZPT,12X,12HDETERMINANTS/)   
6     write(9,3000) (O(1,J),J=1,4),S(1,1),S(1,2),DETRAD        
      write(9,3001) (O(2,J),J=1,4),S(2,1),S(2,2),DETAXI        
      write(9,3002) (O(3,J),J=1,4),RADNOR  
      write(9,3003) (O(4,J),J=1,4),AXINOR  
      write(9,3005)    
      write(9,3006) (Q(1,J),J=6,9),U(1,1),U(1,2)     
      write(9,3006) (Q(2,J),J=6,9),U(2,1),U(2,2)     
      write(9,3007) ((Q(I,J),J=6,9),I=3,4) 
3000  FORMAT(1X,4E12.5,5X,2E12.5,12X,7HRADIAL=E12.5)        
3001  FORMAT(1X,4E12.5,5X,2E12.5,13X,6HAXIAL=E12.5)         
3002  FORMAT(1X,4E12.5,31X,17HRADIAL NORMALISE=E12.5)       
3003  FORMAT(1X,4E12.5,32X,16HAXIAL NORMALISE=E12.5)        
3005  FORMAT(/14X,23HRADIAL  INVERSE  MATRIX,19X, 
     122HAXTIAL  INVERSE  MATRIX/)       
3006  FORMAT(1X,4E12.5,5X,2E12.5)       
3007  FORMAT(1X,4E12.5)       
      write(9,1006)    
      NP=NP+1       
      GO TO (3,4,5,7),NP      
7     RETURN        
      END 