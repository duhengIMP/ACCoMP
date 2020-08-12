c-------calculate matrix value and inv matrix      
c-------O----matrix in R R' Theta Theta'      input
c-------S----matrix in Z Z'                   input 
c
c-------Q----inv matrix of O                  output
c-------U----inv matrix of S                  output
c--------------------------------------------
      SUBROUTINE DETINV(Q,O,DETRAD,U,S,DETAXI,IMAT)         
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Q(4,9),O(4,4),U(2,2),S(2,2)       
      DIMENSION D(4)
c------ hanglie shi      
      DETAXI=S(1,1)*S(2,2)-S(1,2)*S(2,1)
      N=0 
1     N=N+1         
      JB=MOD(N,4)+1 
      JC=MOD(N+1,4)+1         
      JD=MOD(N+2,4)+1         
      D(N)=O(2,JB)*(O(3,JC)*O(4,JD)-O(3,JD)*O(4,JC))+       
     1O(2,JC)*(O(3,JD)*O(4,JB)-O(3,JB)*O(4,JD))+  
     2O(2,JD)*(O(3,JB)*O(4,JC)-O(3,JC)*O(4,JB))   
      IF(N.LE.3) GO TO 1      
      DETRAD=O(1,1)*D(1)-O(1,2)*D(2)+O(1,3)*D(3)-O(1,4)*D(4)
      IF(IMAT.LT.8) GO TO 35  
      U(1,1)=S(2,2)/DETAXI    
      U(1,2)=-S(1,2)/DETAXI   
      U(2,1)=-S(2,1)/DETAXI   
      U(2,2)=S(1,1)/DETAXI  
c------ inv matrix      
      DO 20 I=1,4   
      DO 20 J=1,4   
20    Q(I,J)=O(I,J) 
      NH=4
      NA=NH+1       
      NE=NH+NH+1    
      NP1=NA+1      
      DO 51 I=1,NH  
      DO 51 J=NP1,NE
51    Q(I,J)=0.     
      DO 50 I=1,NH  
      Q(I,NA)=1.    
      IPNA=I+NA     
50    Q(I,IPNA)=1.  
      I=1 
32    N=I+1         
      DO 30 J=N,NE  
30    Q(I,J)=Q(I,J)/Q(I,I)    
      IF(I-NH)33,34,34        
33    DO 31 K=N,NH  
      DO 31 J=N,NE  
31    Q(K,J)=Q(K,J)-Q(K,I)*Q(I,J)       
      I=I+1         
      GO TO 32     
34    N=I-1        
      DO 40 K=1,N  
      DO 40 J=NA,NE
40    Q(K,J)=Q(K,J)-Q(K,I)*Q(I,J)      
      I=I-1        
      IF(I-1)34,35,34        
35    RETURN       
      END