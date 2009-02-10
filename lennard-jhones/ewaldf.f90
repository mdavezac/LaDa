subroutine EWALDF(IPR,EEWALD,FEWA,FEWAC,STRESS, &
     NATOT,RC,ZZ,RCUT,AVEC,MXDNAT)
!
!!C     COMPUTES THE EWALD ENERGY, FORCES AND STRESSES.
!!C     ADAPTED FROM JLM PLANE-WAVE PROGRAM
!!C     MARCH 12, 1996. NB
!!C     THE G-SUM IS ORDER N^2
!C
!!C     INPUT:
!!C     IPR         CONTROLS PRINTING
!!C     NATOT       NUMBER OF ATOMS 
!!C     RC(K,N)     K-TH COMPONENT (IN LATTICE COORDINATES) OF THE
!!C                 POSITION OF THE N-TH ATOM
!!C     ZZ(N)       CHARGE OF ATOM N 
!!C     RCUT        REAL SPACE CUT-OFF FOR EWALD SUM (AU)
!!C     AVEC(I,J)   I-TH COMPONENT OF J-TH DIRECT LATTICE VECTOR
!!C     MXDNAT      ARRAY DIMENSION FOR NUMBER OF ATOMS
!C
!!C     OUTPUT:
!!C     EEWALD      EWALD ENERGY IN RYDBERG
!!C     FEWA(K,N)   K-TH COMPONENT (IN LATTICE COORDINATES) OF THE
!!C                 FORCE OF THE N-TH ATOM (RYD/AU)
!!C     FEWAC(K,N)  K-TH COMPONENT (IN CARTESIAN COORDINATES) OF THE
!!C                 FORCE OF THE N-TH ATOM (RYD/AU)
!!C     STRESS(I)   STRESS TENSOR (RYD) I=XX,YY,ZZ,XY,YZ,ZX
!C
!!C     NOTE: THE SUPERCELL MUST BE NEUTRAL !
!C

  use ep_param    ! XZ

  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  
  PARAMETER (PI = 3.141592653589793E0, TWOPI = 2.0E0*PI)
  PARAMETER (ZERO = 0.0E0, UM = 1.0E0 , DOIS = 2.0E0)
  PARAMETER (TRES = 3.0E0, QUATRO = 4.0E0)
  
!!C     SMALL HAS BEEN CHANGED FROM 1.E-12 TO 1.E-8 FORM MD

  PARAMETER (SMALL = 1.0E-12)
 
  DIMENSION AVEC(3,3),ZZ(MXDNAT),RC(3,MXDNAT)
  DIMENSION FEWA(3,MXDNAT),FEWAC(3,MXDNAT),STRESS(6)
  
  DIMENSION ADOT(3,3),BDOT(3,3),BVEC(3,3) 
  DIMENSION TAU(3),IR(3),IG(3),RP(3)
  DIMENSION SSUMG(6),SSUMR(6),FSUB(3),SSUB(6),SSUM0(6)  
!     WORK ARRAYS FSUMG, FSUMR
!  DIMENSION FSUMR(3,MXDNAT),FSUMG(9,MXDNAT)

  REAL(kind=dbl) :: FSUMR(3,MXDNAT),FSUMG(9,MXDNAT)
  real(kind=dbl), external :: boost_erfc

 
!C     COMPUTE VARIOUS PARAMETERS
 
       TOL = -LOG(SMALL)   
       EPS = TOL/RCUT**2
       SEPS = SQRT(EPS)
       SEPI = DOIS*SEPS/SQRT(PI)
       GCUT2 = 4.*EPS*TOL      
       GCUT = SQRT(GCUT2)
 
!C     BDOT(I,J)   METRIC IN RECIPROCAL SPACE
 
!C     COMPUTE THE LATTICE WAVE-VECTORS, THE CELL VOLUME AND
!C     THE INNER PRODUCTS.
 
       BVEC(1,1) = AVEC(2,2)*AVEC(3,3) - AVEC(3,2)*AVEC(2,3)
       BVEC(2,1) = AVEC(3,2)*AVEC(1,3) - AVEC(1,2)*AVEC(3,3)
       BVEC(3,1) = AVEC(1,2)*AVEC(2,3) - AVEC(2,2)*AVEC(1,3)
       BVEC(1,2) = AVEC(2,3)*AVEC(3,1) - AVEC(3,3)*AVEC(2,1)
       BVEC(2,2) = AVEC(3,3)*AVEC(1,1) - AVEC(1,3)*AVEC(3,1)
       BVEC(3,2) = AVEC(1,3)*AVEC(2,1) - AVEC(2,3)*AVEC(1,1)
       BVEC(1,3) = AVEC(2,1)*AVEC(3,2) - AVEC(3,1)*AVEC(2,2)
       BVEC(2,3) = AVEC(3,1)*AVEC(1,2) - AVEC(1,1)*AVEC(3,2)
       BVEC(3,3) = AVEC(1,1)*AVEC(2,2) - AVEC(2,1)*AVEC(1,2)
 
!C     CELL VOLUME
 
       VCELL = BVEC(1,1)*AVEC(1,1) + BVEC(2,1)*AVEC(2,1) + &
               BVEC(3,1)*AVEC(3,1)

!       write(6,*)' PUPPPPPPPPPPPAAAAAAAAAAAAaa         ', VCELL

       IF(ABS(VCELL) .LT. SMALL) CALL FATAL(50,VCELL,IDUM)
 
       DO 10 J=1,3
         BVEC(1,J) = TWOPI*BVEC(1,J)/VCELL
         BVEC(2,J) = TWOPI*BVEC(2,J)/VCELL
         BVEC(3,J) = TWOPI*BVEC(3,J)/VCELL
 10    CONTINUE
       VCELL = ABS(VCELL) 
       QPV = QUATRO*PI/VCELL
 
!C     COMPUTE METRIC BDOT(I,J)
 
       BDOT(1,1) = BVEC(1,1)*BVEC(1,1) + BVEC(2,1)*BVEC(2,1) + &
                   BVEC(3,1)*BVEC(3,1)
       BDOT(2,2) = BVEC(1,2)*BVEC(1,2) + BVEC(2,2)*BVEC(2,2) + &
                   BVEC(3,2)*BVEC(3,2)
       BDOT(3,3) = BVEC(1,3)*BVEC(1,3) + BVEC(2,3)*BVEC(2,3) + &
                   BVEC(3,3)*BVEC(3,3)
       BDOT(1,2) = BVEC(1,1)*BVEC(1,2) + BVEC(2,1)*BVEC(2,2) + &
                   BVEC(3,1)*BVEC(3,2)
       BDOT(1,3) = BVEC(1,1)*BVEC(1,3) + BVEC(2,1)*BVEC(2,3) + &
                   BVEC(3,1)*BVEC(3,3)
       BDOT(2,3) = BVEC(1,2)*BVEC(1,3) + BVEC(2,2)*BVEC(2,3) + &
                   BVEC(3,2)*BVEC(3,3)
       BDOT(2,1) = BDOT(1,2)
       BDOT(3,1) = BDOT(1,3)
       BDOT(3,2) = BDOT(2,3)
 
       IMX = INT(RCUT*SQRT(BDOT(1,1))/TWOPI) + 1
       JMX = INT(RCUT*SQRT(BDOT(2,2))/TWOPI) + 1
       KMX = INT(RCUT*SQRT(BDOT(3,3))/TWOPI) + 1
 
!C     COMPUTE METRIC IN REAL SPACE
 
       FACTOR = (VCELL/(TWOPI*TWOPI))*(VCELL/(TWOPI*TWOPI))
       ADOT(1,1) = FACTOR*(BDOT(2,2)*BDOT(3,3) - BDOT(2,3)*BDOT(2,3))
       ADOT(2,2) = FACTOR*(BDOT(3,3)*BDOT(1,1) - BDOT(3,1)*BDOT(3,1))
       ADOT(3,3) = FACTOR*(BDOT(1,1)*BDOT(2,2) - BDOT(1,2)*BDOT(1,2))
       ADOT(1,2) = FACTOR*(BDOT(1,3)*BDOT(3,2) - BDOT(1,2)*BDOT(3,3))
       ADOT(1,3) = FACTOR*(BDOT(1,2)*BDOT(2,3) - BDOT(1,3)*BDOT(2,2))
       ADOT(2,3) = FACTOR*(BDOT(2,1)*BDOT(1,3) - BDOT(2,3)*BDOT(1,1))
       ADOT(2,1) = ADOT(1,2)
       ADOT(3,1) = ADOT(1,3)
       ADOT(3,2) = ADOT(2,3)
 
       IMK = INT(GCUT*SQRT(ADOT(1,1))/TWOPI) + 1
       JMK = INT(GCUT*SQRT(ADOT(2,2))/TWOPI) + 1
       KMK = INT(GCUT*SQRT(ADOT(3,3))/TWOPI) + 1
 
!C     INITIALIZE SUMS : ESUM(G,R)   -   ENERGY
!C                       FSUM(G,R)   -   FORCE
!C                       SSUM(G,R)   -   STRESS
 
!C     NOTE THAT THE SUMS ARE IN UNITS OF BASIS VECTORS
!C     (G) IN UNITS OF BI'S, AND (R) IN UNITS OF AI'S
 
       ESUMG = ZERO
       ESUMR = ZERO
       DO 16 I=1,3
         SSUMG(I+3) = ZERO
         SSUMR(I+3) = ZERO
         SSUMG(I)   = ZERO
         SSUMR(I)   = ZERO
         DO 14 J=1,NATOT
           FSUMG(I,J) = ZERO
           FSUMR(I,J) = ZERO
 14      CONTINUE
 16    CONTINUE
 
!C     START SUM IN G SPACE
 
       IM2 = 2*IMK+1
       JM2 = 2*JMK+1
       KM2 = 2*KMK+1
       DO 30 I=1,IM2
       IG(1) = I-IMK-1
       DO 30 J=1,JM2
       IG(2) = J-JMK-1
       DO 30 K=1,KM2
       IG(3) = K-KMK-1
         GMOD2 = ZERO
         DO 18 L=1,3
         DO 18 M=1,3
           GMOD2 = GMOD2 + REAL(IG(L))*BDOT(L,M)*REAL(IG(M))
 18      CONTINUE
         IF (GMOD2 .NE. ZERO) THEN
           ARG = GMOD2/(QUATRO*EPS)       
           EXPG = EXP(-ARG)/GMOD2                      
           SFACR = ZERO     
           SFACI = ZERO 
           DO 20 N1= 1,NATOT        
             GDT = TWOPI*(REAL(IG(1))*RC(1,N1) +  &
                          REAL(IG(2))*RC(2,N1) +  &
                          REAL(IG(3))*RC(3,N1))
             COSG = ZZ(N1)*COS(GDT)  
             SING = ZZ(N1)*SIN(GDT)
             SFACR = SFACR + COSG  
             SFACI = SFACI + SING
             FSUMG(4,N1) = - REAL(IG(1))*COSG 
             FSUMG(5,N1) = - REAL(IG(2))*COSG 
             FSUMG(6,N1) = - REAL(IG(3))*COSG 
             FSUMG(7,N1) = REAL(IG(1))*SING
             FSUMG(8,N1) = REAL(IG(2))*SING
             FSUMG(9,N1) = REAL(IG(3))*SING
 20        CONTINUE
           SFAC2 = SFACR*SFACR + SFACI*SFACI   
           EXP1 = SFAC2*EXPG      
           ESUMG = ESUMG + EXP1
           EXP2 = - (UM/(EPS*DOIS) + DOIS/GMOD2) * EXP1
           EXPGR = SFACR*EXPG
           EXPGI = SFACI*EXPG
           DO 26 N1=1,NATOT   
             FSUMG(1,N1) = FSUMG(1,N1) + FSUMG(7,N1)*EXPGR +  &
                           FSUMG(4,N1)*EXPGI
             FSUMG(2,N1) = FSUMG(2,N1) + FSUMG(8,N1)*EXPGR +  &
                           FSUMG(5,N1)*EXPGI
             FSUMG(3,N1) = FSUMG(3,N1) + FSUMG(9,N1)*EXPGR +  &
                           FSUMG(6,N1)*EXPGI
 26        CONTINUE 
           SSUMG(1) = SSUMG(1) + EXP2 * REAL(IG(1)*IG(1))
           SSUMG(2) = SSUMG(2) + EXP2 * REAL(IG(2)*IG(2))
           SSUMG(3) = SSUMG(3) + EXP2 * REAL(IG(3)*IG(3))
           SSUMG(4) = SSUMG(4) + EXP2 * REAL(IG(1)*IG(2))
           SSUMG(5) = SSUMG(5) + EXP2 * REAL(IG(2)*IG(3))
           SSUMG(6) = SSUMG(6) + EXP2 * REAL(IG(3)*IG(1))
         ENDIF
 30    CONTINUE
!C                                     
       ESUMG = QPV*ESUMG
       DO 32 L=1,6
         SSUMG(L) = QPV*SSUMG(L)
 32    CONTINUE
       DO 36 L=1,3
         DO 34 N1=1,NATOT
           FSUMG(L,N1) = DOIS*QPV*FSUMG(L,N1)
 34      CONTINUE
 36    CONTINUE
 
!C     START SUM IN R SPACE
 
       IM2 = 2*IMX+1
       JM2 = 2*JMX+1
       KM2 = 2*KMX+1
       ESUM0 = ZERO
       DO 38 I = 1,6
         SSUM0(I) = ZERO
 38    CONTINUE
       DO 40 I=1,IM2
       IR(1) = I-IMX-1
       DO 40 J=1,JM2
       IR(2) = J-JMX-1
       DO 40 K=1,KM2
       IR(3) = K-KMX-1
         RMOD = ZERO
         DO 39 L=1,3
         DO 39 M=1,3
           RMOD = RMOD + REAL(IR(L))*ADOT(L,M)*REAL(IR(M))
 39      CONTINUE
         IF (RMOD .NE. ZERO) THEN
           RMOD = SQRT(RMOD)
           ARG = SEPS*RMOD
           IF (ARG .LT. 25.0) THEN
             EXP1 = boost_erfc(ARG) / RMOD
             EXP2 = (EXP1 + SEPI*EXP(-ARG*ARG))/(RMOD*RMOD)
             ESUM0 = ESUM0 + EXP1
             SSUM0(1) = SSUM0(1) + EXP2 * REAL(IR(1)*IR(1))
             SSUM0(2) = SSUM0(2) + EXP2 * REAL(IR(2)*IR(2))
             SSUM0(3) = SSUM0(3) + EXP2 * REAL(IR(3)*IR(3))
             SSUM0(4) = SSUM0(4) + EXP2 * REAL(IR(1)*IR(2))
             SSUM0(5) = SSUM0(5) + EXP2 * REAL(IR(2)*IR(3))
             SSUM0(6) = SSUM0(6) + EXP2 * REAL(IR(3)*IR(1))
           ENDIF
         ENDIF
 40    CONTINUE
       ESUM0 = ESUM0 - SEPI
 
!C     START LOOP OVER ATOMS IN CELL
 
       DO 52 I=1,NATOT
 
!C       TERM WITH A=B
 
         ESUMR = ESUMR + ZZ(I)*ZZ(I)*ESUM0
         DO 42 J=1,6
           SSUMR(J) = SSUMR(J) + ZZ(I)*ZZ(I) * SSUM0(J)
 42      CONTINUE
         IM = I-1
         IF (IM .NE. 0) THEN
 
!C         TERMS WITH A#B
 
           DO 50 J=1,IM
 
!C           LOOP OVER LATTICE POINTS
!C                                   
!C           ATOMS NOT IN THE UNIT CELL ARE SENT BACK
 
             R1CC = RC(1,I) - RC(1,J)    
             R2CC = RC(2,I) - RC(2,J)
             R3CC = RC(3,I) - RC(3,J)
             IR1CC = INT(R1CC)
             IF(R1CC .LT. ZERO)  IR1CC = IR1CC - 1
             IR2CC = INT(R2CC)
             IF(R2CC .LT. ZERO) IR2CC = IR2CC - 1
             IR3CC = INT(R3CC)
             IF(R3CC .LT. ZERO) IR3CC = IR3CC - 1
             R1CC = R1CC - REAL(IR1CC)
             R2CC = R2CC - REAL(IR2CC)
             R3CC = R3CC - REAL(IR3CC)
 
             ESUB = ZERO
             DO 43 K=1,3
               FSUB(K) = ZERO
               SSUB(K+3) = ZERO
               SSUB(K) = ZERO
 43          CONTINUE
             DO 46 II=1,IM2
             IR(1) = II-IMX-1
             DO 46 JJ=1,JM2
             IR(2) = JJ-JMX-1
             DO 46 KK=1,KM2
             IR(3) = KK-KMX-1
               RP(1) = REAL(IR(1)) + R1CC
               RP(2) = REAL(IR(2)) + R2CC
               RP(3) = REAL(IR(3)) + R3CC
               RMOD = ZERO
               DO 44 L=1,3
               DO 44 M=1,3
                 RMOD = RMOD + RP(L)*ADOT(L,M)*RP(M)
 44            CONTINUE
               RMOD = SQRT(RMOD)
               ARG = SEPS*RMOD
               IF (ARG .LT. 25.0) THEN
                 EXP1 = boost_erfc(ARG) / RMOD
                 EXP2 = (EXP1 + SEPI*EXP(-ARG*ARG))/(RMOD*RMOD)
                 ESUB = ESUB + EXP1
                 FSUB(1) = FSUB(1) + RP(1) * EXP2
                 FSUB(2) = FSUB(2) + RP(2) * EXP2
                 FSUB(3) = FSUB(3) + RP(3) * EXP2
                 SSUB(1) = SSUB(1) + RP(1) * EXP2 * RP(1)
                 SSUB(2) = SSUB(2) + RP(2) * EXP2 * RP(2)
                 SSUB(3) = SSUB(3) + RP(3) * EXP2 * RP(3)
                 SSUB(4) = SSUB(4) + RP(1) * EXP2 * RP(2)
                 SSUB(5) = SSUB(5) + RP(2) * EXP2 * RP(3)
                 SSUB(6) = SSUB(6) + RP(3) * EXP2 * RP(1)
               ENDIF
 46          CONTINUE
             ESUMR = ESUMR + DOIS*ZZ(I)*ZZ(J)*ESUB
             DO 48 K=1,6
               SSUMR(K) = SSUMR(K) + DOIS*ZZ(I)*ZZ(J)*SSUB(K)
 48          CONTINUE
             DO 49 K=1,3
               FSUMR(K,I) = FSUMR(K,I) + DOIS*ZZ(I)*ZZ(J)*FSUB(K)
               FSUMR(K,J) = FSUMR(K,J) - DOIS*ZZ(I)*ZZ(J)*FSUB(K)
 49          CONTINUE
 50        CONTINUE
         ENDIF
 52    CONTINUE
 
!C     END R SUM
 
       EEWALD = ESUMG + ESUMR
 
!C     FORCE
!C     NOTE - RETURNED FORCE IN UNITS OF LATTICE VECTORS (FEWA)
!C            AND CARTESIAN COORDINATES (FEWAC)
!C            PRINTED FORCE IN CARTESIAN COORDINATES
 
       DO 62 I=1,NATOT
         DO 60 K=1,3
           FEWA(K,I) = FSUMR(K,I)
           DO 58 L=1,3
             FEWA(K,I) = FEWA(K,I) + BDOT(K,L)*FSUMG(L,I)/TWOPI
 58        CONTINUE
 60      CONTINUE
 62    CONTINUE
 
!C     STRESS
!C     NOTE - BOTH RETURNED AND PRINTED STRESS ARE
!C            IN CARTESIAN COORDINATES
 
       DO 70 I=1,6
         J = I
         K = I
         IF (I .GT. 3) J = I - 3
         IF (I .GT. 3) K = J + 1
         IF (K .GT. 3) K = 1
         STRESS(I) = &
           BVEC(J,1)*SSUMG(1)*BVEC(K,1) + AVEC(J,1)*SSUMR(1)*AVEC(K,1) &
         + BVEC(J,2)*SSUMG(2)*BVEC(K,2) + AVEC(J,2)*SSUMR(2)*AVEC(K,2) &
         + BVEC(J,3)*SSUMG(3)*BVEC(K,3) + AVEC(J,3)*SSUMR(3)*AVEC(K,3) &
         + BVEC(J,1)*SSUMG(4)*BVEC(K,2) + AVEC(J,1)*SSUMR(4)*AVEC(K,2) &
         + BVEC(J,2)*SSUMG(5)*BVEC(K,3) + AVEC(J,2)*SSUMR(5)*AVEC(K,3) &
         + BVEC(J,3)*SSUMG(6)*BVEC(K,1) + AVEC(J,3)*SSUMR(6)*AVEC(K,1) &
         + BVEC(K,1)*SSUMG(4)*BVEC(J,2) + AVEC(K,1)*SSUMR(4)*AVEC(J,2) &
         + BVEC(K,2)*SSUMG(5)*BVEC(J,3) + AVEC(K,2)*SSUMR(5)*AVEC(J,3) &
         + BVEC(K,3)*SSUMG(6)*BVEC(J,1) + AVEC(K,3)*SSUMR(6)*AVEC(J,1)
         IF (I .LE. 3) THEN
           STRESS(I) = STRESS(I) + ESUMG 
         ENDIF
 70    CONTINUE
!C          
!C     FORCES IN CARTESIAN COORDINATES
 
       DO 82 I=1,NATOT
         DO 80 K=1,3
           FEWAC(K,I) = AVEC(K,1)*FEWA(1,I) &
                    + AVEC(K,2)*FEWA(2,I) &
                    + AVEC(K,3)*FEWA(3,I)
 80      CONTINUE
 82    CONTINUE
!C      
!C     PRINTOUT
       ENORM = EEWALD*VCELL**(UM/TRES)
!C      WRITE(6,100) ENORM
       IF(IPR .NE. 0) THEN 
         WRITE(6,101) EEWALD,ENORM 
          DO 84 N=1,NATOT
            WRITE(6,102) N,(RC(K,N),K=1,3),(FEWA(K,N),K=1,3), (FEWAC(K,N),K=1,3)
 84       CONTINUE
 
         WRITE(6,103) (STRESS(I),I=1,6)
         WRITE(6,103) (STRESS(I)/VCELL,I=1,6)
       ENDIF
! change stress
       RETURN
 100   FORMAT(/,26H NORMALIZED EWALD ENERGY =,F16.10)
 101   FORMAT(/,' EWALD ANALYSIS',/,1X,14('*'),/,50X,'1/3',/, &
       ' ENERGY :',12X,'ENERGY (RY)',16X, &
       '*V    (RY*A.U.)',/,16X,F16.10,14X,F16.9,// &
       ' FORCES :',3X,'N',11X,'COORD',19X,'FORCE (RY/A.U.)',/ &
       19X,'A1',4X,'A2',4X,'A3',11X,'-X-',7X,'-Y-',7X,'-Z-')
 102   FORMAT(10X,I3,3X,3F6.3,5X,3F14.6,5X,3F14.6)
 103   FORMAT(/,' STRESS :',25X,'SIGMA * V (RY)',/, &
       14X,'-XX-',6X,'-YY-',6X,'-ZZ-',6X,'-XY-',6X,'-YZ-',6X,'-ZX-',/, &
       9X,6F10.6)

end subroutine EWALDF

subroutine fatal(i,xarg,iarg)
!c
!c      handles the fatal errors
!c      and stops the execution of the program
!c
       implicit double precision (a-h,o-z)
!c
       write(6,1000)
       i10 = i/10
       ir = i - i10*10
       if(i10 .eq. 3) then
!c      rdpp
       if (ir .eq. 1) then
          write(6,1031) iarg
       endif
       endif
!c
       if(i10 .eq. 4) then
!c      vpair
       if (ir .eq. 1) then
          write(6,1041) xarg
       endif
       endif
       if(i10 .eq. 5) then
!c        crstl or updg or forts
         if(ir .eq. 0) then
           write(6,1050) xarg
         else if(ir .eq. 1) then
           write(6,1051) iarg
         else if(ir .eq. 2) then
           write(6,1052) iarg
         else if(ir .eq. 3) then
           write(6,1053) iarg
         endif
       endif
       if (i10 .eq. 6) then 
!c        forts
         if (ir .eq. 1) then 
           write(6,1061) iarg
         endif
       endif
       stop

 1000  format('  ***fatal error***')
 1031  format('  number of mesh points too large = ',i5)
 1041  format('  interatomic distance too small = ',d12.5)
 1050  format('  cell volume = ',d12.5)
 1051  format('  types of atoms = ',i5)
 1052  format('  number of atoms = ',i5)
 1053  format('  number of steps = ',i5)
 1061  format('  potential index = ',i5)

end subroutine fatal 
