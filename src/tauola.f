!     ==================================================================
!     This is a C-library compliant reformating of TAUOLA
!
!     The original FORTRAN code is available from the tauolapp website:
!     https://tauolapp.web.cern.ch/tauolapp (v1.1.8, LHC).
!
!     The main library function is the `tauola_decay` routine, defined
!     below. It wraps TAUOLA internal routines as closures. In addition,
!     common blocks used by TAUOLA are defined with explicit C binding,
!     using the `tauola_` prefix.
!     ==================================================================
      SUBROUTINE TAUOLA_DECAY(KTO,HX) BIND(C)
      INTEGER KTO,ITMP
      REAL*8 HX(4)
      COMMON /IPChT/ IVER
      INTEGER        IVER
      DATA           IVER/1/
      BIND(C,NAME='tauola_ipcht') /IPChT/
      COMMON /TAUOLA_VERSION/ VERSION
      CHARACTER(len=16)::     VERSION='1.1.8-LHC'//CHAR(0)
      BIND(C) /TAUOLA_VERSION/
!     ==================================================================
!     External routines needed by TAUOLA
!     ==================================================================
      INTERFACE
        SUBROUTINE TAUOLA_RANDOM(RVEC,LENV) BIND(C)
          DIMENSION RVEC(*)
        END SUBROUTINE TAUOLA_RANDOM
        SUBROUTINE FILHEP(N,STATUS,PID,MF,ML,DF,DL,P,AM,PFLAG)
     &  BIND(C,NAME='tauola_filhep')
          LOGICAL PFLAG
          INTEGER N,STATUS,PID,MF,DF,DL
          REAL P(4),AM
        END SUBROUTINE FILHEP
        SUBROUTINE TAUOLA_PRINT(S) BIND(C)
          USE ISO_C_BINDING
          IMPLICIT NONE
          CHARACTER(kind=C_CHAR), dimension(*), intent(in) :: S
        END SUBROUTINE TAUOLA_PRINT
        SUBROUTINE TAUOLA_STOP() BIND(C)
        END SUBROUTINE TAUOLA_STOP
      END INTERFACE
!     ==================================================================
!     Call the TAUOLA.DEKAY routine
!     ==================================================================
      IF(KTO.EQ.-1) THEN
        ITMP=IVER ! Backup version
        IVER=1    ! Enable new currents
        CALL RCHL_PARAMETERS(IVER)
        CALL INIMAS()
        CALL INITDK()
        CALL INIPHY(1.0D-1) ! XK0=0.1 is not used
        CALL DEKAY(KTO,HX)
        IVER=ITMP ! Restore version
      ELSE
        IF (IVER.EQ.1) THEN
          CALL RCHL_PARAMETERS(IVER)
        ENDIF
        CALL DEKAY(KTO,HX)
      ENDIF
      CONTAINS
      SUBROUTINE CURR_CLEO(MNUM,PIM1,PIM2,PIM3,PIM4,HADCUR)
C     ==================================================================
C AJW, 11/97 - based on original CURR from TAUOLA:
C     hadronic current for 4 pi final state
C     R. Fisher, J. Wess and F. Wagner Z. Phys C3 (1980) 313
C     R. Decker Z. Phys C36 (1987) 487.
C     M. Gell-Mann, D. Sharp, W. Wagner Phys. Rev. Lett 8 (1962) 261.
C BUT, rewritten to be more general and less "theoretical",
C  using parameters tuned by Vasia and DSC.
C     ==================================================================
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL  PIM1(4),PIM2(4),PIM3(4),PIM4(4)
      COMPLEX HADCUR(4)
      INTEGER K,L,MNUM,K1,K2,IRO,I,J,KK
      REAL PA(4),PB(4),PAA(4)
      REAL AA(4,4),PP(4,4)
      REAL A,XM,XG,G1,G2,G,AMRO2,GAMRO2,AMRO3,GAMRO3,AMOM,GAMOM
      REAL FRO,COEF1,FPI,COEF2,QQ,SK,DENOM,SIG,QQA,SS23,SS24,SS34,QP1P2
      REAL QP1P3,QP1P4,P1P2,P1P3,P1P4,SIGN
      REAL AMPA
      COMPLEX ALF0,ALF1,ALF2,ALF3
      COMPLEX LAM0,LAM1,LAM2,LAM3
      COMPLEX BET1,BET2,BET3
      COMPLEX FORM1,FORM2,FORM3,FORM4,FORM2PI
      COMPLEX AMPL(7),AMPR
      COMPLEX BWIGN
C
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      BWIGN(A,XM,XG)=1.0/CMPLX(A-XM**2,XM*XG)
C*******************************************************************************
C
C --- masses and constants
      IF (G1.NE.12.924) THEN
      G1=12.924
      G2=1475.98
      FPI=93.3E-3
      G =G1*G2
      FRO=0.266*AMRO**2
      COEF1=2.0*SQRT(3.0)/FPI**2
      COEF2=FRO*G ! overall constant for the omega current
      COEF2= COEF2*0.56  ! factor 0.56 reduces contribution of omega from 68% to 40 %
C masses and widths for for rho-prim and rho-bis:
      AMRO2 = 1.465
      GAMRO2= 0.310
      AMRO3=1.700
      GAMRO3=0.235
C
      AMOM  = PKORB(1,14)
      GAMOM = PKORB(2,14)
      AMRO2 = PKORB(1,21)
      GAMRO2= PKORB(2,21)
      AMRO3 = PKORB(1,22)
      GAMRO3= PKORB(2,22)
C
C Amplitudes for (pi-pi-pi0pi+) -> PS, rho0, rho-, rho+, omega.
      AMPL(1) = CMPLX(PKORB(3,31)*COEF1,0.)
      AMPL(2) = CMPLX(PKORB(3,32)*COEF1,0.)*CEXP(CMPLX(0.,PKORB(3,42)))
      AMPL(3) = CMPLX(PKORB(3,33)*COEF1,0.)*CEXP(CMPLX(0.,PKORB(3,43)))
      AMPL(4) = CMPLX(PKORB(3,34)*COEF1,0.)*CEXP(CMPLX(0.,PKORB(3,44)))
      AMPL(5) = CMPLX(PKORB(3,35)*COEF2,0.)*CEXP(CMPLX(0.,PKORB(3,45)))
C Amplitudes for (pi0pi0pi0pi-) -> PS, rho-.
      AMPL(6) = CMPLX(PKORB(3,36)*COEF1)
      AMPL(7) = CMPLX(PKORB(3,37)*COEF1)
C
C rho' contributions to rho' -> pi-omega:
      ALF0 = CMPLX(PKORB(3,51),0.0)
      ALF1 = CMPLX(PKORB(3,52)*AMRO**2,0.0)
      ALF2 = CMPLX(PKORB(3,53)*AMRO2**2,0.0)
      ALF3 = CMPLX(PKORB(3,54)*AMRO3**2,0.0)
C rho' contribtions to rho' -> rhopipi:
      LAM0 = CMPLX(PKORB(3,55),0.0)
      LAM1 = CMPLX(PKORB(3,56)*AMRO**2,0.0)
      LAM2 = CMPLX(PKORB(3,57)*AMRO2**2,0.0)
      LAM3 = CMPLX(PKORB(3,58)*AMRO3**2,0.0)
C rho contributions to rhopipi, rho -> 2pi:
      BET1 = CMPLX(PKORB(3,59)*AMRO**2,0.0)
      BET2 = CMPLX(PKORB(3,60)*AMRO2**2,0.0)
      BET3 = CMPLX(PKORB(3,61)*AMRO3**2,0.0)
C
      END IF
C**************************************************
C
C --- initialization of four vectors
      DO 7 K=1,4
      DO 8 L=1,4
 8    AA(K,L)=0.0
      HADCUR(K)=CMPLX(0.0)
      PAA(K)=PIM1(K)+PIM2(K)+PIM3(K)+PIM4(K)
      PP(1,K)=PIM1(K)
      PP(2,K)=PIM2(K)
      PP(3,K)=PIM3(K)
 7    PP(4,K)=PIM4(K)
C
      IF (MNUM.EQ.1) THEN
C ===================================================================
C pi- pi- p0 pi+ case                                            ====
C ===================================================================
       QQ=PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2
C  Add M(4pi)-dependence to rhopipi channels:
       FORM4= LAM0+LAM1*BWIGN(QQ,AMRO,GAMRO)
     *            +LAM2*BWIGN(QQ,AMRO2,GAMRO2)
     *            +LAM3*BWIGN(QQ,AMRO3,GAMRO3)
C --- loop over five contributions of the rho-pi-pi
       DO 201 K1=1,3
       DO 201 K2=3,4
C
         IF (K2.EQ.K1) THEN
           GOTO 201
         ELSEIF (K2.EQ.3) THEN
C rho-
            AMPR = AMPL(3)
            AMPA = AMPIZ
         ELSEIF (K1.EQ.3) THEN
C rho+
            AMPR = AMPL(4)
            AMPA = AMPIZ
         ELSE
C rho0
            AMPR = AMPL(2)
            AMPA = AMPI
         END IF
C
         SK=(PP(K1,4)+PP(K2,4))**2-(PP(K1,3)+PP(K2,3))**2
     $     -(PP(K1,2)+PP(K2,2))**2-(PP(K1,1)+PP(K2,1))**2
C -- definition of AA matrix
C -- cronecker delta
        DO 202 I=1,4
         DO 203 J=1,4
 203     AA(I,J)=0.0
 202    AA(I,I)=1.0
C ... and the rest ...
        DO 204 L=1,4
         IF (L.NE.K1.AND.L.NE.K2) THEN
          DENOM=(PAA(4)-PP(L,4))**2-(PAA(3)-PP(L,3))**2
     $         -(PAA(2)-PP(L,2))**2-(PAA(1)-PP(L,1))**2
          DO 205 I=1,4
          DO 205 J=1,4
                      SIG= 1.0
           IF(J.NE.4) SIG=-SIG
           AA(I,J)=AA(I,J)
     $            -SIG*(PAA(I)-2.0*PP(L,I))*(PAA(J)-PP(L,J))/DENOM
 205      CONTINUE
         ENDIF
 204    CONTINUE
C
C --- lets add something to HADCURR
C        FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKM(SQRT(QQ),AMPI,AMPI)
C        FORM1= AMPL(1)+AMPR*FPIKM(SQRT(SK),AMPI,AMPI)
        FORM2PI= BET1*BWIGM(SK,AMRO,GAMRO,AMPA,AMPI)
     1          +BET2*BWIGM(SK,AMRO2,GAMRO2,AMPA,AMPI)
     2          +BET3*BWIGM(SK,AMRO3,GAMRO3,AMPA,AMPI)
        FORM1= AMPL(1)+AMPR*FORM2PI
C
       DO 206 I=1,4
       DO 206 J=1,4
        HADCUR(I)=HADCUR(I)+FORM1*FORM4*AA(I,J)*(PP(K1,J)-PP(K2,J))
 206   CONTINUE
C --- end of the rho-pi-pi current (5 possibilities)
 201   CONTINUE
C
C ===================================================================
C Now modify the coefficient for the omega-pi current:              =
C ===================================================================
       IF (AMPL(5).EQ.CMPLX(0.,0.)) GOTO 311
C Overall rho+rhoprime for the 4pi system:
C       FORM2=AMPL(5)*(BWIGN(QQ,AMRO,GAMRO)+ELPHA*BWIGN(QQ,AMROP,GAMROP))
C Modified M(4pi)-dependence:
       FORM2=AMPL(5)*(ALF0+ALF1*BWIGN(QQ,AMRO,GAMRO)
     *                    +ALF2*BWIGN(QQ,AMRO2,GAMRO2)
     *                    +ALF3*BWIGN(QQ,AMRO3,GAMRO3))
C
C --- there are two possibilities for omega current
C --- PA PB are corresponding first and second pi-s
       DO 301 KK=1,2
        DO 302 I=1,4
         PA(I)=PP(KK,I)
         PB(I)=PP(3-KK,I)
 302    CONTINUE
C --- lorentz invariants
         QQA=0.0
         SS23=0.0
         SS24=0.0
         SS34=0.0
         QP1P2=0.0
         QP1P3=0.0
         QP1P4=0.0
         P1P2 =0.0
         P1P3 =0.0
         P1P4 =0.0
        DO 303 K=1,4
                     SIGN=-1.0
         IF (K.EQ.4) SIGN= 1.0
         QQA=QQA+SIGN*(PAA(K)-PA(K))**2
         SS23=SS23+SIGN*(PB(K)  +PIM3(K))**2
         SS24=SS24+SIGN*(PB(K)  +PIM4(K))**2
         SS34=SS34+SIGN*(PIM3(K)+PIM4(K))**2
         QP1P2=QP1P2+SIGN*(PAA(K)-PA(K))*PB(K)
         QP1P3=QP1P3+SIGN*(PAA(K)-PA(K))*PIM3(K)
         QP1P4=QP1P4+SIGN*(PAA(K)-PA(K))*PIM4(K)
         P1P2=P1P2+SIGN*PA(K)*PB(K)
         P1P3=P1P3+SIGN*PA(K)*PIM3(K)
         P1P4=P1P4+SIGN*PA(K)*PIM4(K)
 303    CONTINUE
C
C omega -> rho pi for the 3pi system:
C       FORM3=BWIGN(QQA,AMOM,GAMOM)*(BWIGN(SS23,AMRO,GAMRO)+
C     $        BWIGN(SS24,AMRO,GAMRO)+BWIGN(SS34,AMRO,GAMRO))
C No omega -> rho pi; just straight omega:
        FORM3=BWIGN(QQA,AMOM,GAMOM)
C
        DO 304 K=1,4
         HADCUR(K)=HADCUR(K)+FORM2*FORM3*(
     $             PB  (K)*(QP1P3*P1P4-QP1P4*P1P3)
     $            +PIM3(K)*(QP1P4*P1P2-QP1P2*P1P4)
     $            +PIM4(K)*(QP1P2*P1P3-QP1P3*P1P2) )
 304    CONTINUE
 301   CONTINUE
 311   CONTINUE
C
      ELSE
C ===================================================================
C pi0 pi0 p0 pi- case                                            ====
C ===================================================================
       QQ=PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2
C --- loop over three contribution of the non-omega current
       DO 101 K=1,3
        SK=(PP(K,4)+PIM4(4))**2-(PP(K,3)+PIM4(3))**2
     $    -(PP(K,2)+PIM4(2))**2-(PP(K,1)+PIM4(1))**2
C -- definition of AA matrix
C -- cronecker delta
        DO 102 I=1,4
         DO 103 J=1,4
 103     AA(I,J)=0.0
 102    AA(I,I)=1.0
C
C ... and the rest ...
        DO 104 L=1,3
         IF (L.NE.K) THEN
          DENOM=(PAA(4)-PP(L,4))**2-(PAA(3)-PP(L,3))**2
     $         -(PAA(2)-PP(L,2))**2-(PAA(1)-PP(L,1))**2
          DO 105 I=1,4
          DO 105 J=1,4
                      SIG=1.0
           IF(J.NE.4) SIG=-SIG
           AA(I,J)=AA(I,J)
     $            -SIG*(PAA(I)-2.0*PP(L,I))*(PAA(J)-PP(L,J))/DENOM
 105      CONTINUE
         ENDIF
 104    CONTINUE
C --- lets add something to HADCURR
C       FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKMD(SQRT(QQ),AMPI,AMPI)
CCCCCCCCCCCCC       FORM1=WIGFOR(SK,AMRO,GAMRO)        (tests)
C       FORM1= FPIKM(SQRT(SK),AMPI,AMPI) *FPIKM(SQRT(QQ),AMPI,AMPI)
       FORM1 = AMPL(6)+AMPL(7)*FPIKM(SQRT(SK),AMPI,AMPI)
        DO 106 I=1,4
        DO 106 J=1,4
         HADCUR(I)=HADCUR(I)+FORM1*AA(I,J)*(PP(K,J)-PP(4,J))
 106    CONTINUE
C --- end of the non omega current (3 possibilities)
 101   CONTINUE
      ENDIF
      END
      FUNCTION FORMOM(XMAA,XMOM)
C     ==================================================================
C     formfactorfor pi-pi0 gamma final state
C      R. Decker, Z. Phys C36 (1987) 487.
C     ==================================================================
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON /TESTA1/ KEYA1
      COMPLEX BWIGN,FORMOM
      DATA ICONT /1/
* THIS INLINE FUNCT. CALCULATES THE SCALAR PART OF THE PROPAGATOR
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_testa1') /TESTA1/
!
      BWIGN(XM,AM,GAMMA)=1./CMPLX(XM**2-AM**2,GAMMA*AM)
* HADRON CURRENT
      FRO  =0.266*AMRO**2
      ELPHA=- 0.1
      AMROP = 1.7
      GAMROP= 0.26
      AMOM  =0.782
      GAMOM =0.0085
      AROMEG= 1.0
      GCOUP=12.924
      GCOUP=GCOUP*AROMEG
      FQED  =SQRT(4.0*3.1415926535/137.03604)
      FORMOM=FQED*FRO**2/SQRT(2.0)*GCOUP**2*BWIGN(XMOM,AMOM,GAMOM)
     $     *(BWIGN(XMAA,AMRO,GAMRO)+ELPHA*BWIGN(XMAA,AMROP,GAMROP))
     $     *(BWIGN( 0.0,AMRO,GAMRO)+ELPHA*BWIGN( 0.0,AMROP,GAMROP))
      END
      COMPLEX FUNCTION FK1AB(XMSQ,INDX)
C     ==================================================================
C     complex form-factor for a1+a1prime.                       AJW 1/98
C     ==================================================================
      COMPLEX F1,F2,AMPA,AMPB
      INTEGER IFIRST,INDX
      DATA IFIRST/0/
      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
        XM1 = PKORB(1,19)
        XG1 = PKORB(2,19)
        XM2 = PKORB(1,20)
        XG2 = PKORB(2,20)
        XM1SQ = XM1*XM1
        GF1 = GFUN(XM1SQ)
        GG1 = XM1*XG1/GF1
        XM2SQ = XM2*XM2
        GF2 = GFUN(XM2SQ)
        GG2 = XM2*XG2/GF2
      END IF
      IF (INDX.EQ.1) THEN
        AMPA = CMPLX(PKORB(3,81),0.)
        AMPB = CMPLX(PKORB(3,82),0.)
      ELSE IF (INDX.EQ.2) THEN
        AMPA = CMPLX(PKORB(3,83),0.)
        AMPB = CMPLX(PKORB(3,84),0.)
      ELSEIF (INDX.EQ.3) THEN
        AMPA = CMPLX(PKORB(3,85),0.)
        AMPB = CMPLX(PKORB(3,86),0.)
      ELSEIF (INDX.EQ.4) THEN
        AMPA = CMPLX(PKORB(3,87),0.)
        AMPB = CMPLX(PKORB(3,88),0.)
      END IF
      GF = GFUN(XMSQ)
      FG1 = GG1*GF
      FG2 = GG2*GF
      F1 = CMPLX(-XM1SQ,0.0)/CMPLX(XMSQ-XM1SQ,FG1)
      F2 = CMPLX(-XM2SQ,0.0)/CMPLX(XMSQ-XM2SQ,FG2)
      FK1AB = AMPA*F1+AMPB*F2
      RETURN
      END
      FUNCTION FORM1(MNUM,QQ,S1,SDWA)
C     ==================================================================
C     formfactorfor F1 for 3 scalar final state
C     R. Fisher, J. Wess and F. Wagner Z. Phys C3 (1980) 313
C     H. Georgi, Weak interactions and modern particle theory,
C     The Benjamin/Cummings Pub. Co., Inc. 1984.
C     R. Decker, E. Mirkes, R. Sauer, Z. Was Karlsruhe preprint TTP92-25
C     and erratum !!!!!!
C     ==================================================================
C
      COMPLEX FORM1,WIGNER,WIGFOR
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /IPChT/ IVER
      INTEGER        IVER
      COMPLEX FORMA1,FORMK1,FORMRO,FORMKS
C
!
      BIND(C,NAME='tauola_ipcht') /IPChT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      IF     (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
C       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
C       FORMRO = F3PI(1,QQ,S1,SDWA)
C       FORMA1 = FA1A1P(QQ)
C       FORM1 = FORMA1*FORMRO
      IF (IVER.EQ.0) THEN
        FORM1 = F3PI(1,QQ,S1,SDWA)
      ELSE
        FORM1 = F3PI_RCHT(1,QQ,S1,SDWA)
      ENDIF
      ELSEIF (MNUM.EQ.1) THEN
C ------------ K- pi- K+ (K*0 K-)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM1 = FORMA1*FORMKS
      ELSEIF (MNUM.EQ.2) THEN
C ------------ K0 pi- K0B (K*- K0)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM1 = FORMA1*FORMKS
      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- pi0 K0 (K*0 K-)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM1 = FORMA1*FORMKS
      ELSEIF (MNUM.EQ.4) THEN
C ------------ pi0 pi0 K-  (K*-pi0)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMK1 = FK1AB(QQ,3)
       FORM1 = FORMK1*FORMKS
      ELSEIF (MNUM.EQ.5) THEN
C ------------ K- pi- pi+ (rho0 K-)
       FORMK1 = FK1AB(QQ,4)
       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
       FORM1 = FORMK1*FORMRO
      ELSEIF (MNUM.EQ.6) THEN
C ------------ pi- K0B pi0 (pi- K*0B)
       FORMK1 = FK1AB(QQ,1)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMK,AMPI)
       FORM1 = FORMK1*FORMKS
      ELSEIF (MNUM.EQ.7) THEN
C -------------- eta pi- pi0 final state
       FORM1=0.0
      ENDIF
      END
      FUNCTION FORM2(MNUM,QQ,S1,SDWA)
C     ==================================================================
C     formfactorfor F2 for 3 scalar final state
C     R. Fisher, J. Wess and F. Wagner Z. Phys C3 (1980) 313
C     H. Georgi, Weak interactions and modern particle theory,
C     The Benjamin/Cummings Pub. Co., Inc. 1984.
C     R. Decker, E. Mirkes, R. Sauer, Z. Was Karlsruhe preprint TTP92-25
C     and erratum !!!!!!
C     ==================================================================
C
      COMPLEX FORM2,WIGNER,WIGFOR
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /IPChT/ IVER
      INTEGER        IVER
      COMPLEX FORMA1,FORMK1,FORMRO,FORMKS
!
      BIND(C,NAME='tauola_ipcht') /IPChT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      IF     (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
C       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
C       FORMRO = F3PI(2,QQ,S1,SDWA)
C       FORMA1 = FA1A1P(QQ)
C       FORM2 = FORMA1*FORMRO
C       FORM2 = F3PI(2,QQ,S1,SDWA)
      IF (IVER.EQ.0) THEN
       FORM2 = F3PI(2,QQ,S1,SDWA)
      ELSE
       FORM2 = F3PI_RCHT(2,QQ,S1,SDWA)
      ENDIF
      ELSEIF (MNUM.EQ.1) THEN
C ------------ K- pi- K+ (rho0 pi-)
       FORMRO = FPIKM(SQRT(S1),AMK,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM2 = FORMA1*FORMRO
      ELSEIF (MNUM.EQ.2) THEN
C ------------ K0 pi- K0B (rho0 pi-)
       FORMRO = FPIKM(SQRT(S1),AMK,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM2 = FORMA1*FORMRO
      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- pi0 K0 (rho- pi0)
       FORMRO = FPIKM(SQRT(S1),AMK,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM2 = FORMA1*FORMRO
      ELSEIF (MNUM.EQ.4) THEN
C ------------ pi0 pi0 K-  (K*-pi0)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMK1 = FK1AB(QQ,3)
       FORM2 = FORMK1*FORMKS
      ELSEIF (MNUM.EQ.5) THEN
C ------------ K- pi- pi+  (K*0B pi-)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPI,AMK)
       FORMK1 = FK1AB(QQ,1)
       FORM2 = FORMK1*FORMKS
C
      ELSEIF (MNUM.EQ.6) THEN
C ------------ pi- K0B pi0 (rho- K0B)
       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
       FORMK1 = FK1AB(QQ,2)
       FORM2 = FORMK1*FORMRO
C
      ELSEIF (MNUM.EQ.7) THEN
C -------------- eta pi- pi0 final state
       FORM2=0.0
      ENDIF
C
      END
      COMPLEX FUNCTION BWIGM(S,M,G,XM1,XM2)
C **********************************************************
C     P-WAVE BREIT-WIGNER  FOR RHO
C **********************************************************
      REAL S,M,G,XM1,XM2
      REAL PI,QS,QM,W,GS
      DATA INIT /0/
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0) THEN
      INIT=1
      PI=3.141592654
C -------  BREIT-WIGNER -----------------------
         ENDIF
       IF (S.GT.(XM1+XM2)**2) THEN
         QS=SQRT(ABS((S   -(XM1+XM2)**2)*(S   -(XM1-XM2)**2)))/SQRT(S)
         QM=SQRT(ABS((M**2-(XM1+XM2)**2)*(M**2-(XM1-XM2)**2)))/M
         W=SQRT(S)
         GS=G*(M/W)**2*(QS/QM)**3
       ELSE
         GS=0.0
       ENDIF
         BWIGM=M**2/CMPLX(M**2-S,-SQRT(S)*GS)
      RETURN
      END
      COMPLEX FUNCTION FPIKM(W,XM1,XM2)
C **********************************************************
C     PION FORM FACTOR
C **********************************************************
      REAL ROM,ROG,ROM1,ROG1,BETA1,PI,PIM,S,W
      EXTERNAL BWIG
      DATA  INIT /0/
C
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0 ) THEN
      INIT=1
      PI=3.141592654
      PIM=.140
      ROM=0.773
      ROG=0.145
      ROM1=1.370
      ROG1=0.510
      BETA1=-0.145
      ENDIF
C -----------------------------------------------
      S=W**2
      FPIKM=(BWIGM(S,ROM,ROG,XM1,XM2)+BETA1*BWIGM(S,ROM1,ROG1,XM1,XM2))
     & /(1+BETA1)
      RETURN
      END
      COMPLEX FUNCTION FPIKMD(W,XM1,XM2)
C **********************************************************
C     PION FORM FACTOR
C **********************************************************
      REAL ROM,ROG,ROM1,ROG1,PI,PIM,S,W
      EXTERNAL BWIG
      DATA  INIT /0/
C
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0 ) THEN
      INIT=1
      PI=3.141592654
      PIM=.140
      ROM=0.773
      ROG=0.145
      ROM1=1.500
      ROG1=0.220
      ROM2=1.750
      ROG2=0.120
      BETA=6.5
      DELTA=-26.0
      ENDIF
C -----------------------------------------------
      S=W**2
      FPIKMD=(DELTA*BWIGM(S,ROM,ROG,XM1,XM2)
     $      +BETA*BWIGM(S,ROM1,ROG1,XM1,XM2)
     $      +     BWIGM(S,ROM2,ROG2,XM1,XM2))
     & /(1+BETA+DELTA)
      RETURN
      END
      FUNCTION FORM3(MNUM,QQ,S1,SDWA)
C     ==================================================================
C     formfactorfor F3 for 3 scalar final state
C     R. Fisher, J. Wess and F. Wagner Z. Phys C3 (1980) 313
C     H. Georgi, Weak interactions and modern particle theory,
C     The Benjamin/Cummings Pub. Co., Inc. 1984.
C     R. Decker, E. Mirkes, R. Sauer, Z. Was Karlsruhe preprint TTP92-25
C     and erratum !!!!!!
C     ==================================================================
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /IPChT/ IVER
      INTEGER        IVER
      COMPLEX FORM3
      COMPLEX FORMA1,FORMK1,FORMRO,FORMKS
      COMPLEX F3PI_RCHT
C
!
      BIND(C,NAME='tauola_ipcht') /IPChT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      IF (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
C       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
C       FORMRO = F3PI(3,QQ,S1,SDWA)
C       FORMA1 = FA1A1P(QQ)
C       FORM3 = FORMA1*FORMRO
      IF (IVER.EQ.0) THEN
        FORM3 = F3PI(3,QQ,S1,SDWA)
      ELSE
        FORM3 = (0.,0.)
      ENDIF
      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- pi0 K0  (K*- K0)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMPIZ,AMK)
       FORMA1 = FA1A1P(QQ)
       FORM3 = FORMA1*FORMKS
      ELSEIF (MNUM.EQ.6) THEN
C ------------ pi- K0B pi0 (K*- pi0)
       FORMKS = BWIGM(S1,AMKST,GAMKST,AMK,AMPI)
       FORMK1 = FK1AB(QQ,3)
       FORM3 = FORMK1*FORMKS
      ELSE
       FORM3=CMPLX(0.,0.)
      ENDIF
      END
      FUNCTION FORM4(MNUM,QQ,S1,S2,S3)
C     ==================================================================
C     formfactorfor F4 for 3 scalar final state
C     R. Decker, in preparation
C     R. Decker, E. Mirkes, R. Sauer, Z. Was Karlsruhe preprint TTP92-25
C     and erratum !!!!!!
C     ==================================================================
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /IPChT/ IVER
      INTEGER        IVER
      COMPLEX FORM4,WIGNER
      REAL*4 M
C ---- this formfactor is switched off for cleo version
!
      BIND(C,NAME='tauola_ipcht') /IPChT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
       FORM4=CMPLX(0.0,0.0)
      IF (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
C       FORMRO = FPIKM(SQRT(S1),AMPI,AMPI)
C       FORMRO = F3PI(3,QQ,S1,SDWA)
C       FORMA1 = FA1A1P(QQ)
C       FORM3 = FORMA1*FORMRO
C        FORM4 = CMPLX(-1., 0.)
       IF (IVER.EQ.1) FORM4 = F3PI_RCHT(4,QQ,S1,S2)* CMPLX(0., 1.)
      ENDIF
      END
      FUNCTION FORM5(MNUM,QQ,S1,S2)
C     ==================================================================
C     formfactorfor F5 for 3 scalar final state
C     G. Kramer, W. Palmer, S. Pinsky, Phys. Rev. D30 (1984) 89.
C     G. Kramer, W. Palmer             Z. Phys. C25 (1984) 195.
C     R. Decker, E. Mirkes, R. Sauer, Z. Was Karlsruhe preprint TTP92-25
C     and erratum !!!!!!
C     ==================================================================
C
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMPLEX FORM5,WIGNER
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      IF     (MNUM.EQ.0) THEN
C ------------  3 pi hadronic state (a1)
        FORM5=0.0
      ELSEIF (MNUM.EQ.1) THEN
C ------------ K- pi- K+
         ELPHA=-0.2
         FORM5=FPIKMD(SQRT(QQ),AMPI,AMPI)/(1+ELPHA)
     $        *(       FPIKM(SQRT(S2),AMPI,AMPI)
     $          +ELPHA*BWIGM(S1,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.2) THEN
C ------------ K0 pi- K0B
         ELPHA=-0.2
         FORM5=FPIKMD(SQRT(QQ),AMPI,AMPI)/(1+ELPHA)
     $        *(       FPIKM(SQRT(S2),AMPI,AMPI)
     $          +ELPHA*BWIGM(S1,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.3) THEN
C ------------ K- K0 pi0
        FORM5=0.0
      ELSEIF (MNUM.EQ.4) THEN
C ------------ pi0 pi0 K-
        FORM5=0.0
      ELSEIF (MNUM.EQ.5) THEN
C ------------ K- pi- pi+
        ELPHA=-0.2
        FORM5=BWIGM(QQ,AMKST,GAMKST,AMPI,AMK)/(1+ELPHA)
     $       *(       FPIKM(SQRT(S1),AMPI,AMPI)
     $         +ELPHA*BWIGM(S2,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.6) THEN
C ------------ pi- K0B pi0
        ELPHA=-0.2
        FORM5=BWIGM(QQ,AMKST,GAMKST,AMPI,AMKZ)/(1+ELPHA)
     $       *(       FPIKM(SQRT(S2),AMPI,AMPI)
     $         +ELPHA*BWIGM(S1,AMKST,GAMKST,AMPI,AMK))
      ELSEIF (MNUM.EQ.7) THEN
C -------------- eta pi- pi0 final state
       FORM5=FPIKMD(SQRT(QQ),AMPI,AMPI)*FPIKM(SQRT(S1),AMPI,AMPI)
      ENDIF
C
      END
      FUNCTION WIGFOR(S,XM,XGAM)
      COMPLEX WIGFOR,WIGNOR
      WIGNOR=CMPLX(-XM**2,XM*XGAM)
      WIGFOR=WIGNOR/CMPLX(S-XM**2,XM*XGAM)
      END
      SUBROUTINE JAKER(JAK)
C     *********************
C
C **********************************************************************
C                                                                      *
C           *********TAUOLA LIBRARY: VERSION 2.9 ********              *
C           **************October  2011******************              *
C           **      AUTHORS: S.JADACH, Z.WAS        *****              *
C           **  R. DECKER, M. JEZABEK, J.H.KUEHN,   *****              *
C           ********AVAILABLE FROM: www.cern.ch/wasm ****              *
C           *******PUBLISHED IN COMP. PHYS. COMM.********              *
C           *** 76 (1993) 361                        ****              *
C           *** 64 (1990) 275                        ****              *
C           *** 70 (1992) 69                         ****              *
C           *** CLEO initialization:                 ****              *
C           *** Alain Weinstein www home page:       ****              *
C           *** http://www.cithep.caltech.edu/~ajw/  ****              *
C           *** RChL initialization:                 ****              *
C           *** O. Shekhovtsova, T. Przedzinski,     ****              *
C           *** P. Roig and Z. Was                   ****              *
C           *** IFJPAN-2013-5, UAB-FT-731            ****              *
C **********************************************************************
C
C ----------------------------------------------------------------------
c SUBROUTINE JAKER,
C CHOOSES DECAY MODE ACCORDING TO LIST OF BRANCHING RATIOS
C JAK=1 ELECTRON MODE
C JAK=2 MUON MODE
C JAK=3 PION MODE
C JAK=4 RHO  MODE
C JAK=5 A1   MODE
C JAK=6 K    MODE
C JAK=7 K*   MODE
C JAK=8 nPI  MODE
C
C     called by : DEXAY
C ----------------------------------------------------------------------
      COMMON / TAUBRA / GAMPRT(30),JLIST(30),NCHAN
C      REAL   CUMUL(20)
      REAL   CUMUL(30),RRR(1)
C
!
      BIND(C,NAME='tauola_taubra') /TAUBRA/
!
      IF(NCHAN.LE.0.OR.NCHAN.GT.30) GOTO 902
      CALL RANMAR(RRR,1)
      SUM=0
      DO 20 I=1,NCHAN
      SUM=SUM+GAMPRT(I)
  20  CUMUL(I)=SUM
      DO 25 I=NCHAN,1,-1
      IF(RRR(1).LT.CUMUL(I)/CUMUL(NCHAN)) JI=I
  25  CONTINUE
      JAK=JLIST(JI)
      RETURN
 902  CALL TAUOLA_PRINT('----- JAKER: WRONG NCHAN'//CHAR(0))
!902  PRINT 9020
!9020 FORMAT(' ----- JAKER: WRONG NCHAN')
      CALL TAUOLA_STOP()
      END
      SUBROUTINE DEKAY(KTO,HX)
C     ***********************
C THIS DEKAY IS IN SPIRIT OF THE 'DECAY' WHICH
C WAS INCLUDED IN KORAL-B PROGRAM, COMP. PHYS. COMMUN.
C VOL. 36 (1985) 191, SEE COMMENTS  ON GENERAL PHILOSOPHY THERE.
C KTO=0 INITIALISATION (OBLIGATORY)
C KTO=1,11 DENOTES TAU+ AND KTO=2,12 TAU-
C DEKAY(1,H) AND DEKAY(2,H) IS CALLED INTERNALLY BY MC GENERATOR.
C H DENOTES THE POLARIMETRIC VECTOR, USED BY THE HOST PROGRAM FOR
C CALCULATION OF THE SPIN WEIGHT.
C USER MAY OPTIONALLY CALL DEKAY(11,H) DEKAY(12,H) IN ORDER
C TO TRANSFORM DECAY PRODUCTS TO CMS AND WRITE LUND RECORD IN /LUJETS/.
C KTO=100, PRINT FINAL REPORT  (OPTIONAL).
C DECAY MODES:
C JAK=1 ELECTRON DECAY
C JAK=2 MU  DECAY
C JAK=3 PI  DECAY
C JAK=4 RHO DECAY
C JAK=5 A1  DECAY
C JAK=6 K   DECAY
C JAK=7 K*  DECAY
C JAK=8 NPI DECAY
C JAK=0 INCLUSIVE:  JAK=1,2,3,4,5,6,7,8
      REAL  H(4)
      REAL*8 HX(4)
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / IDFC  / IDF
      COMMON /TAUPOS/ NP1,NP2
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / TAUDCD /IDFFIN(9,NMODE),MULPIK(NMODE)
     &                ,NAMES
      CHARACTER NAMES(NMODE)*31
      COMMON / INOUT / INUT,IOUT
      COMMON /IPChT/ IVER
      INTEGER        IVER
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4),PDUM5(4),HDUM(4)
      REAL  PDUMX(4,9)
      DATA IWARM/0/
!
      BIND(C,NAME='tauola_idfc') /IDFC/
      BIND(C,NAME='tauola_inout') /INOUT/
      BIND(C,NAME='tauola_ipcht') /IPChT/
      BIND(C,NAME='tauola_jaki') /JAKI/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
      BIND(C,NAME='tauola_taudcd') /TAUDCD/
      BIND(C,NAME='tauola_taupos') /TAUPOS/
!
      KTOM=KTO
      IF(KTO.EQ.-1) THEN
C     ==================
C       INITIALISATION OR REINITIALISATION
C       first or second tau positions in HEPEVT as in KORALB/Z
        NP1=3
        NP2=4
        KTOM=1
        IF (IWARM.EQ.1) X=5/(IWARM-1)
        IWARM=1
        CALL TAUOLA_PRINT('tauola.f:909: (suppressed)'//CHAR(0))
!       WRITE(IOUT,7001) JAK1,JAK
        CALL TAUOLA_PRINT('*IVER   = hadronic current version'//CHAR(0))
!       WRITE(IOUT,7002) IVER
        NEVTOT=0
        NEV1=0
        NEV2=0
        IF(JAK1.NE.-1.OR.JAK2.NE.-1) THEN
          CALL DADMEL(-1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DADMMU(-1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DADMPI(-1,IDUM,HDUM,PDUM1,PDUM2)
          CALL DADMRO(-1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4)
          CALL DADMAA(-1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5,JDUM)
          CALL DADMKK(-1,IDUM,HDUM,PDUM1,PDUM2)
          CALL DADMKS(-1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,JDUM)
          CALL DADNEW(-1,IDUM,HDUM,PDUM1,PDUM2,PDUMX,JDUM)
        ENDIF
        DO 21 I=1,30
        NEVDEC(I)=0
        GAMPMC(I)=0
 21     GAMPER(I)=0
      ELSEIF(KTO.EQ.1) THEN
C     =====================
C DECAY OF TAU+ IN THE TAU REST FRAME
        NEVTOT=NEVTOT+1
        IF(IWARM.EQ.0) GOTO 902
        ISGN= IDF/IABS(IDF)
C AJWMOD to change BRs depending on sign:
        CALL TAURDF(KTO)
        CALL DEKAY1(0,H,ISGN)
      ELSEIF(KTO.EQ.2) THEN
C     =================================
C DECAY OF TAU- IN THE TAU REST FRAME
        NEVTOT=NEVTOT+1
        IF(IWARM.EQ.0) GOTO 902
        ISGN=-IDF/IABS(IDF)
C AJWMOD to change BRs depending on sign:
        CALL TAURDF(KTO)
        CALL DEKAY2(0,H,ISGN)
      ELSEIF(KTO.EQ.11) THEN
C     ======================
C REST OF DECAY PROCEDURE FOR ACCEPTED TAU+ DECAY
        NEV1=NEV1+1
        ISGN= IDF/IABS(IDF)
        CALL DEKAY1(1,H,ISGN)
      ELSEIF(KTO.EQ.12) THEN
C     ======================
C REST OF DECAY PROCEDURE FOR ACCEPTED TAU- DECAY
        NEV2=NEV2+1
        ISGN=-IDF/IABS(IDF)
        CALL DEKAY2(1,H,ISGN)
      ELSEIF(KTO.EQ.100) THEN
C     =======================
        IF(JAK1.NE.-1.OR.JAK2.NE.-1) THEN
          CALL DADMEL( 1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DADMMU( 1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
          CALL DADMPI( 1,IDUM,HDUM,PDUM1,PDUM2)
          CALL DADMRO( 1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4)
          CALL DADMAA( 1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5,JDUM)
          CALL DADMKK( 1,IDUM,HDUM,PDUM1,PDUM2)
          CALL DADMKS( 1,IDUM,HDUM,PDUM1,PDUM2,PDUM3,PDUM4,JDUM)
          CALL DADNEW( 1,IDUM,HDUM,PDUM1,PDUM2,PDUMX,JDUM)
          CALL TAUOLA_PRINT('tauola.f:971: (suppressed)'//CHAR(0))
!         WRITE(IOUT,7010) NEV1,NEV2,NEVTOT
          CALL TAUOLA_PRINT('tauola.f:973: (suppressed)'//CHAR(0))
!         WRITE(IOUT,7011) (NEVDEC(I),GAMPMC(I),GAMPER(I),I= 1,7)
          CALL TAUOLA_PRINT('tauola.f:975: (suppressed)'//CHAR(0))
!         WRITE(IOUT,7012)
!    $         (NEVDEC(I),GAMPMC(I),GAMPER(I),NAMES(I-7),I=8,7+NMODE)
          CALL TAUOLA_PRINT('tauola.f:978: (suppressed)'//CHAR(0))
!         WRITE(IOUT,7013)
        ENDIF
      ELSE
C     ====
        GOTO 910
      ENDIF
C     =====
        DO 78 K=1,4
 78     HX(K)=H(K)
      RETURN
!7001 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'*****TAUOLA LIBRARY: VERSION 2.9 ******',9X,1H*,
!    $ /,' *',     25X,'***********October 2011 ***************',9X,1H*,
!    $ /,' *',     25X,'**AUTHORS: S.JADACH, Z.WAS*************',9X,1H*,
!    $ /,' *',     25X,'**R. DECKER, M. JEZABEK, J.H.KUEHN*****',9X,1H*,
!    $ /,' *',     25X,'**AVAILABLE FROM: www.cern.ch/wasm**** ',9X,1H*,
!    $ /,' *',     25X,'***** PUBLISHED IN COMP. PHYS. COMM.***',9X,1H*,
!    $ /,' *',     25X,'0: Physics initialization  CLEO collab ',9X,1H*,
!    $ /,' *',     25X,' see Alain Weinstein www home page:    ',9X,1H*,
!    $ /,' *',     25X,'http://www.cithep.caltech.edu/~ajw/    ',9X,1H*,
!    $ /,' *',     25X,'/korb_doc.html#files                   ',9X,1H*,
!    $ /,' *',     25X,'1: Physics initialization RChL of:     ',9X,1H*,
!    $ /,' *',     25X,' O. Shekhovtsova, T. Przedzinski,      ',9X,1H*,
!    $ /,' *',     25X,' P. Roig and Z. Was                    ',9X,1H*,
!    $ /,' *',     25X,' IFJPAN-2013-5, UAB-FT-731             ',9X,1H*,
!    $ /,' *',     25X,'*******CPC 76 (1993) 361          *****',9X,1H*,
!    $ /,' *',     25X,'**5 or more pi dec.: precision limited ',9X,1H*,
!    $ /,' *',     25X,'****DEKAY ROUTINE: INITIALIZATION******',9X,1H*,
!    $ /,' *',I20  ,5X,'JAK1   = DECAY MODE TAU+               ',9X,1H*,
!    $ /,' *',I20  ,5X,'JAK2   = DECAY MODE TAU-               ',9X,1H*,
!    $  /,1X,15(5H*****)/)
!7002 FORMAT(' *',I20  ,5X,'IVER   = hadronic current version  ',9X,1H*)
!7010 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'*****TAUOLA LIBRARY: VERSION 2.9 ******',9X,1H*,
!    $ /,' *',     25X,'***********October 2011 ***************',9X,1H*,
!    $ /,' *',     25X,'**AUTHORS: S.JADACH, Z.WAS*************',9X,1H*,
!    $ /,' *',     25X,'**R. DECKER, M. JEZABEK, J.H.KUEHN*****',9X,1H*,
!    $ /,' *',     25X,'**AVAILABLE FROM: www.cern.ch/wasm ****',9X,1H*,
!    $ /,' *',     25X,'***** PUBLISHED IN COMP. PHYS. COMM.***',9X,1H*,
!    $ /,' *',     25X,'******* 64 (1990) 275             *****',9X,1H*,
!    $ /,' *',     25X,'******* 70 (1992) 69              *****',9X,1H*,
!    $ /,' *',     25X,'******* 76 (1993) 361             *****',9X,1H*,
!    $ /,' *',     25X,'******* IFJPAN-2013-5, UAB-FT-731    **',9X,1H*,
!    $ /,' *',     25X,'*****DEKAY ROUTINE: FINAL REPORT*******',9X,1H*,
!    $ /,' *',I20  ,5X,'NEV1   = NO. OF TAU+ DECS. ACCEPTED    ',9X,1H*,
!    $ /,' *',I20  ,5X,'NEV2   = NO. OF TAU- DECS. ACCEPTED    ',9X,1H*,
!    $ /,' *',I20  ,5X,'NEVTOT = SUM                           ',9X,1H*,
!    $ /,' *','    NOEVTS ',
!    $   ' PART.WIDTH     ERROR       ROUTINE    DECAY MODE    ',9X,1H*)
!7011 FORMAT(1X,'*'
!    $       ,I10,2F12.7       ,'     DADMEL     ELECTRON      ',9X,1H*
!    $ /,' *',I10,2F12.7       ,'     DADMMU     MUON          ',9X,1H*
!    $ /,' *',I10,2F12.7       ,'     DADMPI     PION          ',9X,1H*
!    $ /,' *',I10,2F12.7,       '     DADMRO     RHO (->2PI)   ',9X,1H*
!    $ /,' *',I10,2F12.7,       '     DADMAA     A1  (->3PI)   ',9X,1H*
!    $ /,' *',I10,2F12.7,       '     DADMKK     KAON          ',9X,1H*
!    $ /,' *',I10,2F12.7,       '     DADMKS     K*            ',9X,1H*)
!7012 FORMAT(1X,'*'
!    $       ,I10,2F12.7,A31                                    ,8X,1H*)
!7013 FORMAT(1X,'*'
!    $       ,20X,'THE ERROR IS RELATIVE AND  PART.WIDTH      ',10X,1H*
!    $ /,' *',20X,'IN UNITS GFERMI**2*MASS**5/192/PI**3       ',10X,1H*
!    $  /,1X,15(5H*****)/)
 902  CALL TAUOLA_PRINT('----- DEKAY: LACK OF INITIALISATION'//CHAR(0))
!902  PRINT 9020
!9020 FORMAT(' ----- DEKAY: LACK OF INITIALISATION')
      CALL TAUOLA_STOP()
 910  CALL TAUOLA_PRINT('----- DEKAY: WRONG VALUE OF KTO'//CHAR(0))
!910  PRINT 9100
!9100 FORMAT(' ----- DEKAY: WRONG VALUE OF KTO ')
      CALL TAUOLA_STOP()
      END
      SUBROUTINE DEKAY1(IMOD,HH,ISGN)
C     *******************************
C THIS ROUTINE  SIMULATES TAU+  DECAY
      COMMON / DECP4 / PP1(4),PP2(4),KF1,KF2
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      REAL  HH(4)
      REAL  HV(4),PNU(4),PPI(4)
      REAL  PWB(4),PMU(4),PNM(4)
      REAL  PRHO(4),PIC(4),PIZ(4)
      REAL  PAA(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PKK(4),PKS(4)
      REAL  PNPI(4,9)
      REAL  PHOT(4)
      REAL  PDUM(4)
      DATA NEV,NPRIN/0,10/
!
      BIND(C,NAME='tauola_decp4') /DECP4/
      BIND(C,NAME='tauola_jaki') /JAKI/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
!
      KTO=1
      IF(JAK1.EQ.-1) RETURN
      IMD=IMOD
      IF(IMD.EQ.0) THEN
C     =================
      JAK=JAK1
      IF(JAK1.EQ.0) CALL JAKER(JAK)
      IF(JAK.EQ.1) THEN
        CALL DADMEL(0, ISGN,HV,PNU,PWB,PMU,PNM,PHOT)
      ELSEIF(JAK.EQ.2) THEN
        CALL DADMMU(0, ISGN,HV,PNU,PWB,PMU,PNM,PHOT)
      ELSEIF(JAK.EQ.3) THEN
        CALL DADMPI(0, ISGN,HV,PPI,PNU)
      ELSEIF(JAK.EQ.4) THEN
        CALL DADMRO(0, ISGN,HV,PNU,PRHO,PIC,PIZ)
      ELSEIF(JAK.EQ.5) THEN
        CALL DADMAA(0, ISGN,HV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
      ELSEIF(JAK.EQ.6) THEN
        CALL DADMKK(0, ISGN,HV,PKK,PNU)
      ELSEIF(JAK.EQ.7) THEN
        CALL DADMKS(0, ISGN,HV,PNU,PKS ,PKK,PPI,JKST)
      ELSE
        CALL DADNEW(0, ISGN,HV,PNU,PWB,PNPI,JAK-7)
      ENDIF
      DO 33 I=1,3
 33   HH(I)=HV(I)
      HH(4)=1.0
      ELSEIF(IMD.EQ.1) THEN
C     =====================
      NEV=NEV+1
        IF (JAK.LT.31) THEN
           NEVDEC(JAK)=NEVDEC(JAK)+1
         ENDIF
      DO 34 I=1,4
 34   PDUM(I)=.0
      IF(JAK.EQ.1) THEN
        CALL DWLUEL(1,ISGN,PNU,PWB,PMU,PNM)
        CALL DWRPH(KTOM,PHOT)
        DO 10 I=1,4
 10     PP1(I)=PMU(I)
      ELSEIF(JAK.EQ.2) THEN
        CALL DWLUMU(1,ISGN,PNU,PWB,PMU,PNM)
        CALL DWRPH(KTOM,PHOT)
        DO 20 I=1,4
 20     PP1(I)=PMU(I)
      ELSEIF(JAK.EQ.3) THEN
        CALL DWLUPI(1,ISGN,PPI,PNU)
        DO 30 I=1,4
 30     PP1(I)=PPI(I)
      ELSEIF(JAK.EQ.4) THEN
        CALL DWLURO(1,ISGN,PNU,PRHO,PIC,PIZ)
        DO 40 I=1,4
 40     PP1(I)=PRHO(I)
      ELSEIF(JAK.EQ.5) THEN
        CALL DWLUAA(1,ISGN,PNU,PAA,PIM1,PIM2,PIPL,JAA)
        DO 50 I=1,4
 50     PP1(I)=PAA(I)
      ELSEIF(JAK.EQ.6) THEN
        CALL DWLUKK(1,ISGN,PKK,PNU)
        DO 60 I=1,4
 60     PP1(I)=PKK(I)
      ELSEIF(JAK.EQ.7) THEN
        CALL DWLUKS(1,ISGN,PNU,PKS,PKK,PPI,JKST)
        DO 70 I=1,4
 70     PP1(I)=PKS(I)
      ELSE
CAM     MULTIPION DECAY
        CALL DWLNEW(1,ISGN,PNU,PWB,PNPI,JAK)
        DO 80 I=1,4
 80     PP1(I)=PWB(I)
      ENDIF
      ENDIF
C     =====
      END
      SUBROUTINE DEKAY2(IMOD,HH,ISGN)
C     *******************************
C THIS ROUTINE  SIMULATES TAU-  DECAY
      COMMON / DECP4 / PP1(4),PP2(4),KF1,KF2
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      REAL  HH(4)
      REAL  HV(4),PNU(4),PPI(4)
      REAL  PWB(4),PMU(4),PNM(4)
      REAL  PRHO(4),PIC(4),PIZ(4)
      REAL  PAA(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PKK(4),PKS(4)
      REAL  PNPI(4,9)
      REAL  PHOT(4)
      REAL  PDUM(4)
      DATA NEV,NPRIN/0,10/
!
      BIND(C,NAME='tauola_decp4') /DECP4/
      BIND(C,NAME='tauola_jaki') /JAKI/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
!
      KTO=2
      IF(JAK2.EQ.-1) RETURN
      IMD=IMOD
      IF(IMD.EQ.0) THEN
C     =================
      JAK=JAK2
      IF(JAK2.EQ.0) CALL JAKER(JAK)
      IF(JAK.EQ.1) THEN
        CALL DADMEL(0, ISGN,HV,PNU,PWB,PMU,PNM,PHOT)
      ELSEIF(JAK.EQ.2) THEN
        CALL DADMMU(0, ISGN,HV,PNU,PWB,PMU,PNM,PHOT)
      ELSEIF(JAK.EQ.3) THEN
        CALL DADMPI(0, ISGN,HV,PPI,PNU)
      ELSEIF(JAK.EQ.4) THEN
        CALL DADMRO(0, ISGN,HV,PNU,PRHO,PIC,PIZ)
      ELSEIF(JAK.EQ.5) THEN
        CALL DADMAA(0, ISGN,HV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
      ELSEIF(JAK.EQ.6) THEN
        CALL DADMKK(0, ISGN,HV,PKK,PNU)
      ELSEIF(JAK.EQ.7) THEN
        CALL DADMKS(0, ISGN,HV,PNU,PKS ,PKK,PPI,JKST)
      ELSE
        CALL DADNEW(0, ISGN,HV,PNU,PWB,PNPI,JAK-7)
      ENDIF
      DO 33 I=1,3
 33   HH(I)=HV(I)
      HH(4)=1.0
      ELSEIF(IMD.EQ.1) THEN
C     =====================
      NEV=NEV+1
        IF (JAK.LT.31) THEN
           NEVDEC(JAK)=NEVDEC(JAK)+1
         ENDIF
      DO 34 I=1,4
 34   PDUM(I)=.0
      IF(JAK.EQ.1) THEN
        CALL DWLUEL(2,ISGN,PNU,PWB,PMU,PNM)
        CALL DWRPH(KTOM,PHOT)
        DO 10 I=1,4
 10     PP2(I)=PMU(I)
      ELSEIF(JAK.EQ.2) THEN
        CALL DWLUMU(2,ISGN,PNU,PWB,PMU,PNM)
        CALL DWRPH(KTOM,PHOT)
        DO 20 I=1,4
 20     PP2(I)=PMU(I)
      ELSEIF(JAK.EQ.3) THEN
        CALL DWLUPI(2,ISGN,PPI,PNU)
        DO 30 I=1,4
 30     PP2(I)=PPI(I)
      ELSEIF(JAK.EQ.4) THEN
        CALL DWLURO(2,ISGN,PNU,PRHO,PIC,PIZ)
        DO 40 I=1,4
 40     PP2(I)=PRHO(I)
      ELSEIF(JAK.EQ.5) THEN
        CALL DWLUAA(2,ISGN,PNU,PAA,PIM1,PIM2,PIPL,JAA)
        DO 50 I=1,4
 50     PP2(I)=PAA(I)
      ELSEIF(JAK.EQ.6) THEN
        CALL DWLUKK(2,ISGN,PKK,PNU)
        DO 60 I=1,4
 60     PP1(I)=PKK(I)
      ELSEIF(JAK.EQ.7) THEN
        CALL DWLUKS(2,ISGN,PNU,PKS,PKK,PPI,JKST)
        DO 70 I=1,4
 70     PP1(I)=PKS(I)
      ELSE
CAM     MULTIPION DECAY
        CALL DWLNEW(2,ISGN,PNU,PWB,PNPI,JAK)
        DO 80 I=1,4
 80     PP1(I)=PWB(I)
      ENDIF
C
      ENDIF
C     =====
      END
      SUBROUTINE DADMEL(MODE,ISGN,HHV,PNU,PWB,Q1,Q2,PHX)
C ----------------------------------------------------------------------
C
C     called by : DEXEL,(DEKAY,DEKAY1)
C ----------------------------------------------------------------------
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      REAL*4         PHX(4)
      COMMON / INOUT / INUT,IOUT
      REAL  HHV(4),HV(4),PWB(4),PNU(4),Q1(4),Q2(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4),PDUM5(4)
      REAL*4 RRR(3)
      REAL*8 SWT, SSWT
      DATA PI /3.141592653589793238462643/
      DATA IWARM/0/
C
      COMMON /TAUOLA_WEIGHT_DADMEL/ WTMAX
      REAL*4 WTMAX
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_inout') /INOUT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
      BIND(C,NAME='tauola_weight_dadmel') /TAUOLA_WEIGHT_DADMEL/
!
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        NEVRAW=0
        NEVACC=0
        NEVOVR=0
        SWT=0
        SSWT=0
        WTMAX=1E-20
        DO 15 I=1,500
        CALL DPHSEL(WT,HV,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
        IF(WT.GT.WTMAX/1.2) WTMAX=WT*1.2
15      CONTINUE
CC      CALL HBOOK1(803,'WEIGHT DISTRIBUTION  DADMEL    $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
300     CONTINUE
        IF(IWARM.EQ.0) GOTO 902
        NEVRAW=NEVRAW+1
        CALL DPHSEL(WT,HV,PNU,PWB,Q1,Q2,PHX)
CC      CALL HFILL(803,WT/WTMAX)
        SWT=SWT+WT
        SSWT=SSWT+WT**2
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX) NEVOVR=NEVOVR+1
        IF(RN*WTMAX.GT.WT) GOTO 300
C ROTATIONS TO BASIC TAU REST FRAME
        RR2=RRR(2)
        COSTHE=-1.+2.*RR2
        THET=ACOS(COSTHE)
        RR3=RRR(3)
        PHI =2*PI*RR3
        CALL ROTOR2(THET,PNU,PNU)
        CALL ROTOR3( PHI,PNU,PNU)
        CALL ROTOR2(THET,PWB,PWB)
        CALL ROTOR3( PHI,PWB,PWB)
        CALL ROTOR2(THET,Q1,Q1)
        CALL ROTOR3( PHI,Q1,Q1)
        CALL ROTOR2(THET,Q2,Q2)
        CALL ROTOR3( PHI,Q2,Q2)
        CALL ROTOR2(THET,HV,HV)
        CALL ROTOR3( PHI,HV,HV)
        CALL ROTOR2(THET,PHX,PHX)
        CALL ROTOR3( PHI,PHX,PHX)
        DO 44,I=1,3
 44     HHV(I)=-ISGN*HV(I)
        NEVACC=NEVACC+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVRAW.EQ.0) RETURN
        PARGAM=SWT/FLOAT(NEVRAW+1)
        ERROR=0
        IF(NEVRAW.NE.0) ERROR=SQRT(SSWT/SWT**2-1./FLOAT(NEVRAW))
        RAT=PARGAM/GAMEL
        CALL TAUOLA_PRINT('tauola.f:1335: (suppressed)'//CHAR(0))
!       WRITE(IOUT, 7010) NEVRAW,NEVACC,NEVOVR,PARGAM,RAT,ERROR
CC      CALL HPRINT(803)
        GAMPMC(1)=RAT
        GAMPER(1)=ERROR
CAM     NEVDEC(1)=NEVACC
      ENDIF
C     =====
      RETURN
!7010 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADMEL FINAL REPORT  ******** ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVRAW = NO. OF EL  DECAYS TOTAL       ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVACC = NO. OF EL   DECS. ACCEPTED    ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
!    $ /,' *',E20.5,5X,'PARTIAL WTDTH ( ELECTRON) IN GEV UNITS ',9X,1H*
!    $ /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
!    $ /,' *',F20.9,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
!    $ /,' *',25X,     'COMPLETE QED CORRECTIONS INCLUDED      ',9X,1H*
!    $ /,' *',25X,     'BUT ONLY V-A CUPLINGS                  ',9X,1H*
!    $  /,1X,15(5H*****)/)
 902  CALL TAUOLA_PRINT('----- DADMEL: LACK OF INITIALISATION'//CHAR(0))
!902  WRITE(IOUT, 9020)
!9020 FORMAT(' ----- DADMEL: LACK OF INITIALISATION')
      CALL TAUOLA_STOP()
      END
      SUBROUTINE DADMMU(MODE,ISGN,HHV,PNU,PWB,Q1,Q2,PHX)
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL*4         PHX(4)
      REAL  HHV(4),HV(4),PNU(4),PWB(4),Q1(4),Q2(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4),PDUM5(4)
      REAL*4 RRR(3)
      REAL*8 SWT, SSWT
      DATA PI /3.141592653589793238462643/
      DATA IWARM /0/
C
      COMMON /TAUOLA_WEIGHT_DADMMU/ WTMAX
      REAL*4 WTMAX
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_inout') /INOUT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
      BIND(C,NAME='tauola_weight_dadmmu') /TAUOLA_WEIGHT_DADMMU/
!
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        NEVRAW=0
        NEVACC=0
        NEVOVR=0
        SWT=0
        SSWT=0
        WTMAX=1E-20
        DO 15 I=1,500
        CALL DPHSMU(WT,HV,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5)
        IF(WT.GT.WTMAX/1.2) WTMAX=WT*1.2
15      CONTINUE
CC      CALL HBOOK1(802,'WEIGHT DISTRIBUTION  DADMMU    $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
300     CONTINUE
        IF(IWARM.EQ.0) GOTO 902
        NEVRAW=NEVRAW+1
        CALL DPHSMU(WT,HV,PNU,PWB,Q1,Q2,PHX)
CC      CALL HFILL(802,WT/WTMAX)
        SWT=SWT+WT
        SSWT=SSWT+WT**2
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX) NEVOVR=NEVOVR+1
        IF(RN*WTMAX.GT.WT) GOTO 300
C ROTATIONS TO BASIC TAU REST FRAME
        COSTHE=-1.+2.*RRR(2)
        THET=ACOS(COSTHE)
        PHI =2*PI*RRR(3)
        CALL ROTOR2(THET,PNU,PNU)
        CALL ROTOR3( PHI,PNU,PNU)
        CALL ROTOR2(THET,PWB,PWB)
        CALL ROTOR3( PHI,PWB,PWB)
        CALL ROTOR2(THET,Q1,Q1)
        CALL ROTOR3( PHI,Q1,Q1)
        CALL ROTOR2(THET,Q2,Q2)
        CALL ROTOR3( PHI,Q2,Q2)
        CALL ROTOR2(THET,HV,HV)
        CALL ROTOR3( PHI,HV,HV)
        CALL ROTOR2(THET,PHX,PHX)
        CALL ROTOR3( PHI,PHX,PHX)
        DO 44,I=1,3
 44     HHV(I)=-ISGN*HV(I)
        NEVACC=NEVACC+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVRAW.EQ.0) RETURN
        PARGAM=SWT/FLOAT(NEVRAW+1)
        ERROR=0
        IF(NEVRAW.NE.0) ERROR=SQRT(SSWT/SWT**2-1./FLOAT(NEVRAW))
        RAT=PARGAM/GAMEL
        CALL TAUOLA_PRINT('tauola.f:1446: (suppressed)'//CHAR(0))
!       WRITE(IOUT, 7010) NEVRAW,NEVACC,NEVOVR,PARGAM,RAT,ERROR
CC      CALL HPRINT(802)
        GAMPMC(2)=RAT
        GAMPER(2)=ERROR
CAM     NEVDEC(2)=NEVACC
      ENDIF
C     =====
      RETURN
!7010 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADMMU FINAL REPORT  ******** ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVRAW = NO. OF MU  DECAYS TOTAL       ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVACC = NO. OF MU   DECS. ACCEPTED    ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
!    $ /,' *',E20.5,5X,'PARTIAL WTDTH (MU  DECAY) IN GEV UNITS ',9X,1H*
!    $ /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
!    $ /,' *',F20.9,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
!    $ /,' *',25X,     'COMPLETE QED CORRECTIONS INCLUDED      ',9X,1H*
!    $ /,' *',25X,     'BUT ONLY V-A CUPLINGS                  ',9X,1H*
!    $  /,1X,15(5H*****)/)
 902  CALL TAUOLA_PRINT('----- DADMMU: LACK OF INITIALISATION'//CHAR(0))
!902  WRITE(IOUT, 9020)
!9020 FORMAT(' ----- DADMMU: LACK OF INITIALISATION')
      CALL TAUOLA_STOP()
      END
      SUBROUTINE DPHSEL(DGAMX,HVX,XNX,PAAX,QPX,XAX,PHX)
C XNX,XNA was flipped in parameters of dphsel and dphsmu
C *********************************************************************
C *   ELECTRON DECAY MODE                                             *
C *********************************************************************
      REAL*4         PHX(4)
      REAL*4  HVX(4),PAAX(4),XAX(4),QPX(4),XNX(4)
      REAL*8  HV(4),PH(4),PAA(4),XA(4),QP(4),XN(4)
      REAL*8  DGAMT
      IELMU=1
      CALL DRCMU(DGAMT,HV,PH,PAA,XA,QP,XN,IELMU)
      DO 7 K=1,4
        HVX(K)=HV(K)
        PHX(K)=PH(K)
        PAAX(K)=PAA(K)
        XAX(K)=XA(K)
        QPX(K)=QP(K)
        XNX(K)=XN(K)
  7   CONTINUE
      DGAMX=DGAMT
      END
      SUBROUTINE DPHSMU(DGAMX,HVX,XNX,PAAX,QPX,XAX,PHX)
C XNX,XNA was flipped in parameters of dphsel and dphsmu
C *********************************************************************
C *   MUON     DECAY MODE                                             *
C *********************************************************************
      REAL*4         PHX(4)
      REAL*4  HVX(4),PAAX(4),XAX(4),QPX(4),XNX(4)
      REAL*8  HV(4),PH(4),PAA(4),XA(4),QP(4),XN(4)
      REAL*8  DGAMT
      IELMU=2
      CALL DRCMU(DGAMT,HV,PH,PAA,XA,QP,XN,IELMU)
      DO 7 K=1,4
        HVX(K)=HV(K)
        PHX(K)=PH(K)
        PAAX(K)=PAA(K)
        XAX(K)=XA(K)
        QPX(K)=QP(K)
        XNX(K)=XN(K)
  7   CONTINUE
      DGAMX=DGAMT
      END
      SUBROUTINE DRCMU(DGAMT,HV,PH,PAA,XA,QP,XN,IELMU)
      IMPLICIT REAL*8 (A-H,O-Z)
C ----------------------------------------------------------------------
* IT SIMULATES E,MU CHANNELS OF TAU  DECAY IN ITS REST FRAME WITH
* QED ORDER ALPHA CORRECTIONS
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / INOUT / INUT,IOUT
      COMMON / TAURAD / XK0DEC,ITDKRC
      REAL*8            XK0DEC
      REAL*8  HV(4),PT(4),PH(4),PAA(4),XA(4),QP(4),XN(4)
      REAL*8  PR(4)
      REAL*4 RRR(6)
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_inout') /INOUT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taurad') /TAURAD/
!
      LOGICAL IHARD
      DATA PI /3.141592653589793238462643D0/
C AJWMOD to satisfy compiler, comment out this unused function.
C AMRO, GAMRO IS ONLY A PARAMETER FOR GETING HIGHT EFFICIENCY
C
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1./2**17/PI**8
      AMTAX=AMTAU
C TAU MOMENTUM
      PT(1)=0.D0
      PT(2)=0.D0
      PT(3)=0.D0
      PT(4)=AMTAX
C
      CALL RANMAR(RRR,6)
C
        IF (IELMU.EQ.1) THEN
          AMU=AMEL
        ELSE
          AMU=AMMU
        ENDIF
C
        PRHARD=0.30D0
        IF (  ITDKRC.EQ.0) PRHARD=0D0
        PRSOFT=1.-PRHARD
         IF(PRSOFT.LT.0.1) THEN
           CALL TAUOLA_PRINT('ERROR IN DRCMU; PRSOFT='//CHAR(0))
!          PRINT *, 'ERROR IN DRCMU; PRSOFT=',PRSOFT
           CALL TAUOLA_STOP()
         ENDIF
C
        RR5=RRR(5)
        IHARD=(RR5.GT.PRSOFT)
       IF (IHARD) THEN
C                     TAU DECAY TO 'TAU+photon'
          RR1=RRR(1)
          AMS1=(AMU+AMNUTA)**2
          AMS2=(AMTAX)**2
          XK1=1-AMS1/AMS2
          XL1=LOG(XK1/2/XK0DEC)
          XL0=LOG(2*XK0DEC)
          XK=EXP(XL1*RR1+XL0)
          AM3SQ=(1-XK)*AMS2
          AM3 =SQRT(AM3SQ)
          PHSPAC=PHSPAC*AMS2*XL1*XK
          PHSPAC=PHSPAC/PRHARD
        ELSE
          AM3=AMTAX
          PHSPAC=PHSPAC*2**6*PI**3
          PHSPAC=PHSPAC/PRSOFT
        ENDIF
C MASS OF NEUTRINA SYSTEM
        RR2=RRR(2)
        AMS1=(AMNUTA)**2
        AMS2=(AM3-AMU)**2
CAM
CAM
* FLAT PHASE SPACE;
      AM2SQ=AMS1+   RR2*(AMS2-AMS1)
      AM2 =SQRT(AM2SQ)
      PHSPAC=PHSPAC*(AMS2-AMS1)
* NEUTRINA REST FRAME, DEFINE XN AND XA
        ENQ1=(AM2SQ+AMNUTA**2)/(2*AM2)
        ENQ2=(AM2SQ-AMNUTA**2)/(2*AM2)
        PPI=         ENQ1**2-AMNUTA**2
        PPPI=SQRT(ABS(ENQ1**2-AMNUTA**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* NU TAU IN NUNU REST FRAME
        CALL SPHERD(PPPI,XN)
        XN(4)=ENQ1
* NU LIGHT IN NUNU REST FRAME
        DO 30 I=1,3
 30     XA(I)=-XN(I)
        XA(4)=ENQ2
* TAU-prim REST FRAME, DEFINE QP (muon
*       NUNU  MOMENTUM
        PR(1)=0
        PR(2)=0
        PR(4)=1.D0/(2*AM3)*(AM3**2+AM2**2-AMU**2)
        PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
        PPI  =          PR(4)**2-AM2**2
*       MUON MOMENTUM
        QP(1)=0
        QP(2)=0
        QP(4)=1.D0/(2*AM3)*(AM3**2-AM2**2+AMU**2)
        QP(3)=-PR(3)
      PHSPAC=PHSPAC*(4*PI)*(2*PR(3)/AM3)
* NEUTRINA BOOSTED FROM THEIR FRAME TO TAU-prim REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTD3(EXE,XN,XN)
      CALL BOSTD3(EXE,XA,XA)
      RR3=RRR(3)
      RR4=RRR(4)
      IF (IHARD) THEN
        EPS=4*(AMU/AMTAX)**2
        XL1=LOG((2+EPS)/EPS)
        XL0=LOG(EPS)
        ETA  =EXP(XL1*RR3+XL0)
        CTHET=1+EPS-ETA
        THET =ACOS(CTHET)
        PHSPAC=PHSPAC*XL1/2*ETA
        PHI = 2*PI*RR4
        CALL ROTPOX(THET,PHI,XN)
        CALL ROTPOX(THET,PHI,XA)
        CALL ROTPOX(THET,PHI,QP)
        CALL ROTPOX(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE TAU-prim AND GAMMA MOMENTA
* tau-prim  MOMENTUM
        PAA(1)=0
        PAA(2)=0
        PAA(4)=1/(2*AMTAX)*(AMTAX**2+AM3**2)
        PAA(3)= SQRT(ABS(PAA(4)**2-AM3**2))
        PPI   =          PAA(4)**2-AM3**2
        PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAX)
* GAMMA MOMENTUM
        PH(1)=0
        PH(2)=0
        PH(4)=PAA(3)
        PH(3)=-PAA(3)
* ALL MOMENTA BOOSTED FROM TAU-prim REST FRAME TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO PHOTON MOMENTUM
        EXE=(PAA(4)+PAA(3))/AM3
        CALL BOSTD3(EXE,XN,XN)
        CALL BOSTD3(EXE,XA,XA)
        CALL BOSTD3(EXE,QP,QP)
        CALL BOSTD3(EXE,PR,PR)
      ELSE
        THET =ACOS(-1.+2*RR3)
        PHI = 2*PI*RR4
        CALL ROTPOX(THET,PHI,XN)
        CALL ROTPOX(THET,PHI,XA)
        CALL ROTPOX(THET,PHI,QP)
        CALL ROTPOX(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE TAU-prim AND GAMMA MOMENTA
* tau-prim  MOMENTUM
        PAA(1)=0
        PAA(2)=0
        PAA(4)=AMTAX
        PAA(3)=0
* GAMMA MOMENTUM
        PH(1)=0
        PH(2)=0
        PH(4)=0
        PH(3)=0
      ENDIF
C PARTIAL WIDTH CONSISTS OF PHASE SPACE AND AMPLITUDE
      CALL DAMPRY(ITDKRC,XK0DEC,PH,XA,QP,XN,AMPLIT,HV)
      DGAMT=1/(2.*AMTAX)*AMPLIT*PHSPAC
      END
      SUBROUTINE DAMPRY(ITDKRC,XK0DEC,XK,XA,QP,XN,AMPLIT,HV)
      IMPLICIT REAL*8 (A-H,O-Z)
C ----------------------------------------------------------------------
C IT CALCULATES MATRIX ELEMENT FOR THE
C TAU --> MU(E) NU NUBAR DECAY MODE
C INCLUDING COMPLETE ORDER ALPHA QED CORRECTIONS.
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      REAL*8  HV(4),QP(4),XN(4),XA(4),XK(4)
C
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      HV(4)=1.D0
      AK0=XK0DEC*AMTAU
      IF(XK(4).LT.0.1D0*AK0) THEN
        AMPLIT=THB(ITDKRC,QP,XN,XA,AK0,HV)
      ELSE
        AMPLIT=SQM2(ITDKRC,QP,XN,XA,XK,AK0,HV)
      ENDIF
      RETURN
      END
      FUNCTION SQM2(ITDKRC,QP,XN,XA,XK,AK0,HV)
C
C **********************************************************************
C     REAL PHOTON MATRIX ELEMENT SQUARED                               *
C     PARAMETERS:                                                      *
C     HV- POLARIMETRIC FOUR-VECTOR OF TAU                              *
C     QP,XN,XA,XK - 4-momenta of electron (muon), NU, NUBAR and PHOTON *
C                   All four-vectors in TAU rest frame (in GeV)        *
C     AK0 - INFRARED CUTOFF, MINIMAL ENERGY OF HARD PHOTONS (GEV)      *
C     SQM2 - value for S=0                                             *
C     see Eqs. (2.9)-(2.10) from CJK ( Nucl.Phys.B(1991) )             *
C **********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
      REAL*8    QP(4),XN(4),XA(4),XK(4)
      REAL*8    R(4)
      REAL*8   HV(4)
      REAL*8 S0(3),RXA(3),RXK(3),RQP(3)
      DATA PI /3.141592653589793238462643D0/
C
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_qedprm') /QEDPRM/
!
      TMASS=AMTAU
      GF=GFERMI
      ALPHAI=ALFINV
      TMASS2=TMASS**2
      EMASS2=QP(4)**2-QP(1)**2-QP(2)**2-QP(3)**2
      R(4)=TMASS
C     SCALAR PRODUCTS OF FOUR-MOMENTA
      DO 7 I=1,3
        R(1)=0.D0
        R(2)=0.D0
        R(3)=0.D0
        R(I)=TMASS
        RXA(I)=R(4)*XA(4)-R(1)*XA(1)-R(2)*XA(2)-R(3)*XA(3)
C       RXN(I)=R(4)*XN(4)-R(1)*XN(1)-R(2)*XN(2)-R(3)*XN(3)
        RXK(I)=R(4)*XK(4)-R(1)*XK(1)-R(2)*XK(2)-R(3)*XK(3)
        RQP(I)=R(4)*QP(4)-R(1)*QP(1)-R(2)*QP(2)-R(3)*QP(3)
  7   CONTINUE
      QPXN=QP(4)*XN(4)-QP(1)*XN(1)-QP(2)*XN(2)-QP(3)*XN(3)
      QPXA=QP(4)*XA(4)-QP(1)*XA(1)-QP(2)*XA(2)-QP(3)*XA(3)
      QPXK=QP(4)*XK(4)-QP(1)*XK(1)-QP(2)*XK(2)-QP(3)*XK(3)
c     XNXA=XN(4)*XA(4)-XN(1)*XA(1)-XN(2)*XA(2)-XN(3)*XA(3)
      XNXK=XN(4)*XK(4)-XN(1)*XK(1)-XN(2)*XK(2)-XN(3)*XK(3)
      XAXK=XA(4)*XK(4)-XA(1)*XK(1)-XA(2)*XK(2)-XA(3)*XK(3)
      TXN=TMASS*XN(4)
      TXA=TMASS*XA(4)
      TQP=TMASS*QP(4)
      TXK=TMASS*XK(4)
C
      X= XNXK/QPXN
      Z= TXK/TQP
      A= 1+X
      B= 1+ X*(1+Z)/2+Z/2
      S1= QPXN*TXA*( -EMASS2/QPXK**2*A + 2*TQP/(QPXK*TXK)*B-
     $TMASS2/TXK**2)  +
     $QPXN/TXK**2* ( TMASS2*XAXK - TXA*TXK+ XAXK*TXK) -
     $TXA*TXN/TXK - QPXN/(QPXK*TXK)* (TQP*XAXK-TXK*QPXA)
      CONST4=256*PI/ALPHAI*GF**2
      IF (ITDKRC.EQ.0) CONST4=0D0
      SQM2=S1*CONST4
      DO 5 I=1,3
        S0(I) = QPXN*RXA(I)*(-EMASS2/QPXK**2*A + 2*TQP/(QPXK*TXK)*B-
     $  TMASS2/TXK**2) +
     $  QPXN/TXK**2* (TMASS2*XAXK - TXA*RXK(I)+ XAXK*RXK(I))-
     $  RXA(I)*TXN/TXK - QPXN/(QPXK*TXK)*(RQP(I)*XAXK- RXK(I)*QPXA)
  5     HV(I)=S0(I)/S1-1.D0
      RETURN
      END
      FUNCTION THB(ITDKRC,QP,XN,XA,AK0,HV)
C
C **********************************************************************
C     BORN +VIRTUAL+SOFT PHOTON MATRIX ELEMENT**2  O(ALPHA)            *
C     PARAMETERS:                                                      *
C     HV- POLARIMETRIC FOUR-VECTOR OF TAU                              *
C     QP,XN,XA - FOUR-MOMENTA OF ELECTRON (MUON), NU AND NUBAR IN GEV  *
C     ALL FOUR-VECTORS IN TAU REST FRAME                               *
C     AK0 - INFRARED CUTOFF, MINIMAL ENERGY OF HARD PHOTONS            *
C     THB - VALUE FOR S=0                                              *
C     SEE EQS. (2.2),(2.4)-(2.5) FROM CJK (NUCL.PHYS.B351(1991)70      *
C     AND (C.2) FROM JK (NUCL.PHYS.B320(1991)20 )                      *
C **********************************************************************
C
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_qedprm') /QEDPRM/
!
      DIMENSION QP(4),XN(4),XA(4)
      REAL*8 HV(4)
      DIMENSION R(4)
      REAL*8 RXA(3),RXN(3),RQP(3)
      REAL*8 BORNPL(3),AM3POL(3),XM3POL(3)
      DATA PI /3.141592653589793238462643D0/
C
      TMASS=AMTAU
      GF=GFERMI
      ALPHAI=ALFINV
C
      TMASS2=TMASS**2
      R(4)=TMASS
      DO 7 I=1,3
        R(1)=0.D0
        R(2)=0.D0
        R(3)=0.D0
        R(I)=TMASS
        RXA(I)=R(4)*XA(4)-R(1)*XA(1)-R(2)*XA(2)-R(3)*XA(3)
        RXN(I)=R(4)*XN(4)-R(1)*XN(1)-R(2)*XN(2)-R(3)*XN(3)
C       RXK(I)=R(4)*XK(4)-R(1)*XK(1)-R(2)*XK(2)-R(3)*XK(3)
        RQP(I)=R(4)*QP(4)-R(1)*QP(1)-R(2)*QP(2)-R(3)*QP(3)
  7   CONTINUE
C     QUASI TWO-BODY VARIABLES
      U0=QP(4)/TMASS
      U3=SQRT(QP(1)**2+QP(2)**2+QP(3)**2)/TMASS
      W3=U3
      W0=(XN(4)+XA(4))/TMASS
      UP=U0+U3
      UM=U0-U3
      WP=W0+W3
      WM=W0-W3
      YU=LOG(UP/UM)/2
      YW=LOG(WP/WM)/2
      EPS2=U0**2-U3**2
      EPS=SQRT(EPS2)
      Y=W0**2-W3**2
      AL=AK0/TMASS
C     FORMFACTORS
      F0=2*U0/U3*(  DILOGT(1-(UM*WM/(UP*WP)))- DILOGT(1-WM/WP) +
     $DILOGT(1-UM/UP) -2*YU+ 2*LOG(UP)*(YW+YU) ) +
     $1/Y* ( 2*U3*YU + (1-EPS2- 2*Y)*LOG(EPS) ) +
     $ 2 - 4*(U0/U3*YU -1)* LOG(2*AL)
      FP= YU/(2*U3)*(1 + (1-EPS2)/Y ) + LOG(EPS)/Y
      FM= YU/(2*U3)*(1 - (1-EPS2)/Y ) - LOG(EPS)/Y
      F3= EPS2*(FP+FM)/2
C     SCALAR PRODUCTS OF FOUR-MOMENTA
      QPXN=QP(4)*XN(4)-QP(1)*XN(1)-QP(2)*XN(2)-QP(3)*XN(3)
      QPXA=QP(4)*XA(4)-QP(1)*XA(1)-QP(2)*XA(2)-QP(3)*XA(3)
      XNXA=XN(4)*XA(4)-XN(1)*XA(1)-XN(2)*XA(2)-XN(3)*XA(3)
      TXN=TMASS*XN(4)
      TXA=TMASS*XA(4)
      TQP=TMASS*QP(4)
C     DECAY DIFFERENTIAL WIDTH WITHOUT AND WITH POLARIZATION
      CONST3=1/(2*ALPHAI*PI)*64*GF**2
      IF (ITDKRC.EQ.0) CONST3=0D0
      XM3= -( F0* QPXN*TXA +  FP*EPS2* TXN*TXA +
     $FM* QPXN*QPXA + F3* TMASS2*XNXA )
      AM3=XM3*CONST3
C V-A  AND  V+A COUPLINGS, BUT IN THE BORN PART ONLY
      BRAK= (GV+GA)**2*TQP*XNXA+(GV-GA)**2*TXA*QPXN
     &     -(GV**2-GA**2)*TMASS*AMNUTA*QPXA
      BORN= 32*(GFERMI**2/2.)*BRAK
      DO 5 I=1,3
        XM3POL(I)= -( F0* QPXN*RXA(I) +  FP*EPS2* TXN*RXA(I) +
     $  FM* QPXN* (QPXA + (RXA(I)*TQP-TXA*RQP(I))/TMASS2 ) +
     $  F3* (TMASS2*XNXA +TXN*RXA(I) -RXN(I)*TXA)  )
        AM3POL(I)=XM3POL(I)*CONST3
C V-A  AND  V+A COUPLINGS, BUT IN THE BORN PART ONLY
        BORNPL(I)=BORN+(
     &            (GV+GA)**2*TMASS*XNXA*QP(I)
     &           -(GV-GA)**2*TMASS*QPXN*XA(I)
     &           +(GV**2-GA**2)*AMNUTA*TXA*QP(I)
     &           -(GV**2-GA**2)*AMNUTA*TQP*XA(I) )*
     &                                             32*(GFERMI**2/2.)
  5     HV(I)=(BORNPL(I)+AM3POL(I))/(BORN+AM3)-1.D0
      THB=BORN+AM3
      IF (THB/BORN.LT.0.1D0) THEN
        CALL TAUOLA_PRINT('ERROR IN THB, THB/BORN='//CHAR(0))
!       PRINT *, 'ERROR IN THB, THB/BORN=',THB/BORN
        THB=0.D0
      ENDIF
      RETURN
      END
      SUBROUTINE DADMPI(MODE,ISGN,HV,PPI,PNU)
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL  PPI(4),PNU(4),HV(4)
      DATA PI /3.141592653589793238462643/
C
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_inout') /INOUT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
!
      IF(MODE.EQ.-1) THEN
C     ===================
        NEVTOT=0
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
        NEVTOT=NEVTOT+1
        EPI= (AMTAU**2+AMPI**2-AMNUTA**2)/(2*AMTAU)
        ENU= (AMTAU**2-AMPI**2+AMNUTA**2)/(2*AMTAU)
        XPI= SQRT(EPI**2-AMPI**2)
C PI MOMENTUM
        CALL SPHERA(XPI,PPI)
        PPI(4)=EPI
C TAU-NEUTRINO MOMENTUM
        DO 30 I=1,3
30      PNU(I)=-PPI(I)
        PNU(4)=ENU
        PXQ=AMTAU*EPI
        PXN=AMTAU*ENU
        QXN=PPI(4)*PNU(4)-PPI(1)*PNU(1)-PPI(2)*PNU(2)-PPI(3)*PNU(3)
        BRAK=(GV**2+GA**2)*(2*PXQ*QXN-AMPI**2*PXN)
     &      +(GV**2-GA**2)*AMTAU*AMNUTA*AMPI**2
        DO 40 I=1,3
40      HV(I)=-ISGN*2*GA*GV*AMTAU*(2*PPI(I)*QXN-PNU(I)*AMPI**2)/BRAK
        HV(4)=1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVTOT.EQ.0) RETURN
        FPI=0.1284
C        GAMM=(GFERMI*FPI)**2/(16.*PI)*AMTAU**3*
C     *       (BRAK/AMTAU**4)**2
CZW 7.02.93 here was an error affecting non standard model
C       configurations only
        GAMM=(GFERMI*FPI)**2/(16.*PI)*AMTAU**3*
     $       (BRAK/AMTAU**4)*
     $       SQRT((AMTAU**2-AMPI**2-AMNUTA**2)**2
     $            -4*AMPI**2*AMNUTA**2           )/AMTAU**2
        ERROR=0
        RAT=GAMM/GAMEL
        CALL TAUOLA_PRINT('tauola.f:1980: (suppressed)'//CHAR(0))
!       WRITE(IOUT, 7010) NEVTOT,GAMM,RAT,ERROR
        GAMPMC(3)=RAT
        GAMPER(3)=ERROR
CAM     NEVDEC(3)=NEVTOT
      ENDIF
C     =====
      RETURN
!7010 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADMPI FINAL REPORT  ******** ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVTOT = NO. OF PI  DECAYS TOTAL       ',9X,1H*
!    $ /,' *',E20.5,5X,'PARTIAL WTDTH ( PI DECAY) IN GEV UNITS ',9X,1H*
!    $ /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
!    $ /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH (STAT.)',9X,1H*
!    $  /,1X,15(5H*****)/)
      END
      SUBROUTINE DADMRO(MODE,ISGN,HHV,PNU,PRO,PIC,PIZ)
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL  HHV(4)
      REAL  HV(4),PRO(4),PNU(4),PIC(4),PIZ(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4)
      REAL*4 RRR(3)
      REAL*8 SWT, SSWT
      DATA PI /3.141592653589793238462643/
      DATA IWARM/0/
C
      COMMON /TAUOLA_WEIGHT_DADMRO/ WTMAX
      REAL*4 WTMAX
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_inout') /INOUT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
      BIND(C,NAME='tauola_weight_dadmro') /TAUOLA_WEIGHT_DADMRO/
!
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        NEVRAW=0
        NEVACC=0
        NEVOVR=0
        SWT=0
        SSWT=0
        WTMAX=1E-20
        DO 15 I=1,500
        CALL DPHSRO(WT,HV,PDUM1,PDUM2,PDUM3,PDUM4)
        IF(WT.GT.WTMAX/1.2) WTMAX=WT*1.2
15      CONTINUE
CC      CALL HBOOK1(801,'WEIGHT DISTRIBUTION  DADMRO    $',100,0,2)
CC      PRINT 7003,WTMAX
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
300     CONTINUE
        IF(IWARM.EQ.0) GOTO 902
        CALL DPHSRO(WT,HV,PNU,PRO,PIC,PIZ)
CC      CALL HFILL(801,WT/WTMAX)
        NEVRAW=NEVRAW+1
        SWT=SWT+WT
        SSWT=SSWT+WT**2
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX) NEVOVR=NEVOVR+1
        IF(RN*WTMAX.GT.WT) GOTO 300
C ROTATIONS TO BASIC TAU REST FRAME
        COSTHE=-1.+2.*RRR(2)
        THET=ACOS(COSTHE)
        PHI =2*PI*RRR(3)
        CALL ROTOR2(THET,PNU,PNU)
        CALL ROTOR3( PHI,PNU,PNU)
        CALL ROTOR2(THET,PRO,PRO)
        CALL ROTOR3( PHI,PRO,PRO)
        CALL ROTOR2(THET,PIC,PIC)
        CALL ROTOR3( PHI,PIC,PIC)
        CALL ROTOR2(THET,PIZ,PIZ)
        CALL ROTOR3( PHI,PIZ,PIZ)
        CALL ROTOR2(THET,HV,HV)
        CALL ROTOR3( PHI,HV,HV)
        DO 44 I=1,3
 44     HHV(I)=-ISGN*HV(I)
        NEVACC=NEVACC+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVRAW.EQ.0) RETURN
        PARGAM=SWT/FLOAT(NEVRAW+1)
        ERROR=0
        IF(NEVRAW.NE.0) ERROR=SQRT(SSWT/SWT**2-1./FLOAT(NEVRAW))
        RAT=PARGAM/GAMEL
        CALL TAUOLA_PRINT('tauola.f:2081: (suppressed)'//CHAR(0))
!       WRITE(IOUT, 7010) NEVRAW,NEVACC,NEVOVR,PARGAM,RAT,ERROR
CC      CALL HPRINT(801)
        GAMPMC(4)=RAT
        GAMPER(4)=ERROR
CAM     NEVDEC(4)=NEVACC
      ENDIF
C     =====
      RETURN
!7003 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADMRO INITIALISATION ********',9X,1H*
!    $ /,' *',E20.5,5X,'WTMAX  = MAXIMUM WEIGHT                ',9X,1H*
!    $  /,1X,15(5H*****)/)
!7010 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADMRO FINAL REPORT  ******** ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVRAW = NO. OF RHO DECAYS TOTAL       ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVACC = NO. OF RHO  DECS. ACCEPTED    ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
!    $ /,' *',E20.5,5X,'PARTIAL WTDTH (RHO DECAY) IN GEV UNITS ',9X,1H*
!    $ /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
!    $ /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
!    $  /,1X,15(5H*****)/)
 902  CALL TAUOLA_PRINT('----- DADMRO: LACK OF INITIALISATION'//CHAR(0))
!902  WRITE(IOUT, 9020)
!9020 FORMAT(' ----- DADMRO: LACK OF INITIALISATION')
      CALL TAUOLA_STOP()
      END
      SUBROUTINE DPHSRO(DGAMT,HV,PN,PR,PIC,PIZ)
C ----------------------------------------------------------------------
C IT SIMULATES RHO DECAY IN TAU REST FRAME WITH
C Z-AXIS ALONG RHO MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      REAL  HV(4),PT(4),PN(4),PR(4),PIC(4),PIZ(4),QQ(4),RR1(1)
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
C
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      PHSPAC=1./2**11/PI**5
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C MASS OF (REAL/VIRTUAL) RHO
      AMS1=(AMPI+AMPIZ)**2
      AMS2=(AMTAU-AMNUTA)**2
C FLAT PHASE SPACE
C     AMX2=AMS1+   RR1*(AMS2-AMS1)
C     AMX=SQRT(AMX2)
C     PHSPAC=PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH SAMPLING FOR RHO RESONANCE
      ALP1=ATAN((AMS1-AMRO**2)/AMRO/GAMRO)
      ALP2=ATAN((AMS2-AMRO**2)/AMRO/GAMRO)
CAM
 100  CONTINUE
      CALL RANMAR(RR1,1)
      ALP=ALP1+RR1(1)*(ALP2-ALP1)
      AMX2=AMRO**2+AMRO*GAMRO*TAN(ALP)
      AMX=SQRT(AMX2)
      IF(AMX.LT.2.*AMPI) GO TO 100
CAM
      PHSPAC=PHSPAC*((AMX2-AMRO**2)**2+(AMRO*GAMRO)**2)/(AMRO*GAMRO)
      PHSPAC=PHSPAC*(ALP2-ALP1)
C
C TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX**2)
      PN(3)=-SQRT(ABS((PN(4)-AMNUTA)*(PN(4)+AMNUTA)))
C RHO MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX**2)
      PR(3)=-PN(3)
      PHSPAC=PHSPAC*(4*PI)*(2*PR(3)/AMTAU)
C
CAM
      ENQ1=(AMX2+AMPI**2-AMPIZ**2)/(2.*AMX)
      ENQ2=(AMX2-AMPI**2+AMPIZ**2)/(2.*AMX)
      PPPI=SQRT((ENQ1-AMPI)*(ENQ1+AMPI))
      PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMX)
C CHARGED PI MOMENTUM IN RHO REST FRAME
      CALL SPHERA(PPPI,PIC)
      PIC(4)=ENQ1
C NEUTRAL PI MOMENTUM IN RHO REST FRAME
      DO 20 I=1,3
20    PIZ(I)=-PIC(I)
      PIZ(4)=ENQ2
      EXE=(PR(4)+PR(3))/AMX
C PIONS BOOSTED FROM RHO REST FRAME TO TAU REST FRAME
      CALL BOSTR3(EXE,PIC,PIC)
      CALL BOSTR3(EXE,PIZ,PIZ)
      CALL DAM2PI(0,PT,PN,PIC,PIZ,AMPLIT,HV)
      DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
      RETURN
      END
      SUBROUTINE DAM2PI(MNUM,PT,PN,PIM1,PIM2,AMPLIT,HV)
C ----------------------------------------------------------------------
* CALCULATES DIFFERENTIAL CROSS SECTION AND POLARIMETER VECTOR
* FOR TAU DECAY INTO 2 scalar MODES
* ALL SPIN EFFECTS IN THE FULL DECAY CHAIN ARE TAKEN INTO ACCOUNT.
* CALCULATIONS DONE IN TAU REST FRAME WITH Z-AXIS ALONG NEUTRINO MOMENT
C MNUM DECAY MODE IDENTIFIER.
C
C     called by : DPHSAA
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  HV(4),PT(4),PN(4),PIM1(4),PIM2(4)
      REAL  PIVEC(4),PIAKS(4),HVM(4)
      COMPLEX HADCUR(4)
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
C
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      IF     (MNUM.EQ.0) THEN
       CALL CURR_PIPI0(PIM1,PIM2,HADCUR)
      ELSEIF (MNUM.EQ.1) THEN
        CALL CURR_PIK0(PIM1,PIM2,HADCUR)
      ELSEIF (MNUM.EQ.2) THEN
        CALL CURR_KPI0(PIM1,PIM2,HADCUR)
      ELSEIF (MNUM.EQ.3) THEN
        CALL CURR_KK0(PIM1,PIM2,HADCUR)
      ELSE
       CALL TAUOLA_PRINT('DAM2PI: wrong MNUM='//CHAR(0))
!      write(*,*) 'DAM2PI: wrong MNUM= ',MNUM
        CALL TAUOLA_STOP()
      ENDIF
C
* CALCULATE PI-VECTORS: VECTOR AND AXIAL
      CALL CLVEC(HADCUR,PN,PIVEC)
      CALL CLAXI(HADCUR,PN,PIAKS)
      CALL CLNUT(HADCUR,BRAKM,HVM)
* SPIN INDEPENDENT PART OF DECAY DIFF-CROSS-SECT. IN TAU REST  FRAME
      BRAK= (GV**2+GA**2)*PT(4)*PIVEC(4) +2.*GV*GA*PT(4)*PIAKS(4)
     &     +2.*(GV**2-GA**2)*AMNUTA*AMTAU*BRAKM
      IF (MNUM.EQ.0.OR.MNUM.EQ.3) THEN
        AMPLIT=(CCABIB*GFERMI)**2*BRAK
      ELSE
        AMPLIT=(SCABIB*GFERMI)**2*BRAK
      ENDIF
C POLARIMETER VECTOR IN TAU REST FRAME
      DO 90 I=1,3
      HV(I)=-(AMTAU*((GV**2+GA**2)*PIAKS(I)+2.*GV*GA*PIVEC(I)))
     &      +(GV**2-GA**2)*AMNUTA*AMTAU*HVM(I)
C HV IS DEFINED FOR TAU-    WITH GAMMA=B+HV*POL
      HV(I)=-HV(I)/BRAK
 90   CONTINUE
      END
      SUBROUTINE CURR_PIPI0(PC,PN,HADCUR)
C standard TAUOLA current for tau to pi pi0 nu decay
C now it has universal form eg. it is straighforward to add
C scalar part
C NOTE:
C       PC 4-momentum of pi
C       PN 4-momentum of pi0
C       06.08.2011
      IMPLICIT NONE
      COMPLEX BWIGS,HADCUR(4),FKPIPL,FRHO_PI
      REAL  PC(4),PN(4),QQ(4),PKS(4)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON /IPChT/ IVER
      INTEGER        IVER
      INTEGER FF2PIRHO
!
      BIND(C,NAME='tauola_ipcht') /IPChT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      REAL PKSD,QQPKS
      INTEGER IK,K
        DO IK=1,4
         PKS(IK)=PC(IK)+ PN(IK)
          QQ(IK)=PC(IK)- PN(IK)
        ENDDO
C QQ transverse to PKS
        PKSD =PKS(4)*PKS(4)-PKS(3)*PKS(3)-PKS(2)*PKS(2)-PKS(1)*PKS(1)
        QQPKS=PKS(4)* QQ(4)-PKS(3)* QQ(3)-PKS(2)* QQ(2)-PKS(1)* QQ(1)
        DO 31 IK=1,4
 31      QQ(IK)=QQ(IK)-PKS(IK)*QQPKS/PKSD
      IF (IVER.EQ.1) THEN
       CALL GETFF2PIRHO(FF2PIRHO)
       IF (FF2PIRHO.EQ.2) THEN ! Belle,
C                                  ! all fit parameters, par(1...11), are free
        DO K=1,4
         HADCUR(K)=QQ(k)* fpibel(sqrt(pksd),0)
        ENDDO
       ELSEIF (FF2PIRHO.EQ.3) THEN ! Belle
c                             ! all fit parameter free except for
c                             !  par(1)=F_pi(0)=1-fixed
        DO K=1,4
         HADCUR(K)=QQ(k)* fpibel(sqrt(pksd),1)
        ENDDO
       ELSE
        CALL TAUOLA_PRINT('problem in 2-scalars current FF2PIRHO='//CHAR(0))
!       write(*,*) 'problem in 2-scalars current FF2PIRHO=',FF2PIRHO
        CALL TAUOLA_STOP()
       ENDIF
      ELSEIF (IVER.EQ.0) THEN ! cleo
        DO K=1,4
         HADCUR(K)=QQ(k)* sqrt(fpirho(sqrt(pksd)))
        ENDDO
      ELSE
        CALL TAUOLA_PRINT('problem in 2-scalars current IVER='//CHAR(0))
!       write(*,*) 'problem in 2-scalars current IVER=',IVER
        CALL TAUOLA_STOP()
      ENDIF
      END
      SUBROUTINE CURR_PIK0(PC,PN,HADCUR)
C standard TAUOLA current for tau to pi K0 nu decay
C now it has universal form eg. it is straighforward to add
C scalar part
C NOTE:
C       PC 4-momentum of pi
C       PN 4-momentum of K0
C       06.08.2011
      implicit none
      COMPLEX HADCUR(4),FKPIPL
      REAL  PC(4),PN(4),QQ(4),PKS(4),FKPISC,PKSD,QQPKS
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     &                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     &                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     &                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     &                 ,AMK,AMKZ,AMKST,GAMKST,FACT_K0PI
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      Integer           I,K
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taukle') /TAUKLE/
!
        DO I=1,4
         PKS(I)=PC(I)+ PN(I)
          QQ(I)=PC(I)- PN(I)
        ENDDO
C QQ transverse to PKS
        PKSD =PKS(4)*PKS(4)-PKS(3)*PKS(3)-PKS(2)*PKS(2)-PKS(1)*PKS(1)
        QQPKS=PKS(4)* QQ(4)-PKS(3)* QQ(3)-PKS(2)* QQ(2)-PKS(1)* QQ(1)
        DO 31 I=1,4
 31      QQ(I)=QQ(I)-PKS(I)*QQPKS/PKSD
        DO K=1,4
         HADCUR(K)=QQ(k)*BWIGS(pksd,AMKST,GAMKST)
        ENDDO
      END
      SUBROUTINE CURR_KPI0(PC,PN,HADCUR)
C standard TAUOLA current for tau to pi pi0 nu decay
C now it has universal form eg. it is straighforward to add
C scalar part
C NOTE:
C       PC 4-momentum of K
C       PN 4-momentum of pi0
C       06.08.2011
      implicit none
      COMPLEX HADCUR(4),FKPIPL
      REAL  PC(4),PN(4),QQ(4),PKS(4),FKPISC,PKSD,QQPKS
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST,FACT_KPI0
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      INTEGER        I,K
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taukle') /TAUKLE/
!
        DO 30 I=1,4
         PKS(I)=PC(I)+ PN(I)
 30       QQ(I)=PC(I)- PN(I)
C QQ transverse to PKS
        PKSD =PKS(4)*PKS(4)-PKS(3)*PKS(3)-PKS(2)*PKS(2)-PKS(1)*PKS(1)
        QQPKS=PKS(4)* QQ(4)-PKS(3)* QQ(3)-PKS(2)* QQ(2)-PKS(1)* QQ(1)
        DO 31 I=1,4
 31      QQ(I)=QQ(I)-PKS(I)*QQPKS/PKSD
        DO K=1,4
         HADCUR(K)=QQ(k)*BWIGS(pksd,AMKST,GAMKST)
        ENDDO
      END
      SUBROUTINE CURR_KK0(PC,PN,HADCUR)
C standard TAUOLA current for tau to K K0 nu decay
C now it has universal form eg. it is straighforward to add
C scalar part
C NOTE:
C       PC 4-momentum of K
C       PN 4-momentum of K0
C       06.08.2011
      IMPLICIT NONE
      COMPLEX BWIGS,HADCUR(4),FKK0_RCHT
      REAL  PC(4),PN(4),QQ(4),PKS(4),PKSD,QQPKS
      INTEGER I,K
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
        DO I=1,4
         PKS(I)=PC(I)+ PN(I)
          QQ(I)=PC(I)- PN(I)
        ENDDO
C QQ transverse to PKS
        PKSD =PKS(4)*PKS(4)-PKS(3)*PKS(3)-PKS(2)*PKS(2)-PKS(1)*PKS(1)
        QQPKS=PKS(4)* QQ(4)-PKS(3)* QQ(3)-PKS(2)* QQ(2)-PKS(1)* QQ(1)
        DO 31 I=1,4
 31      QQ(I)=QQ(I)-PKS(I)*QQPKS/PKSD
        DO K=1,4
         HADCUR(K)=QQ(k)*sqrt(fpirk(sqrt(pksd)))
        ENDDO
      END
      FUNCTION COEF(I,J)
C clebsh gordan (or so ...)  coefs for 3 scalar final states
      implicit none
C IVER=0  TAUOLA cleo COEF(I,J) =  COEFc(I,J)
C IVER=1  TAUOLA RChL COEF(I,J) =  COEFr(I,J)
      COMMON /IPChT/ IVER
      INTEGER        IVER
      REAL COEFc(1:5,0:7)
      REAL COEFr(1:5,0:7)
      REAL COEF
      DATA PI /3.141592653589793238462643/
      REAL PI
      DATA ICONT /0/
      INTEGER ICONT
      INTEGER I,J
      REAL FPIc,FPIr
C initialization of FPI matrix defined in ...
C FPIc is to be used with cleo initialization
C actual choice is made in ???
      DATA  FPIc /93.3E-3/
C initialization of COEF matrix defined in ...
C COEFc is to be used with cleo initialization
!
      BIND(C,NAME='tauola_ipcht') /IPChT/
!
      IF (ICONT.EQ.0) THEN
       ICONT=1
C
C*****COEFc(I,J)
       COEFc(1,0)= 2.0*SQRT(2.)/3.0
       COEFc(2,0)=-2.0*SQRT(2.)/3.0
C AJW 2/98: Add in the D-wave and I=0 3pi substructure:
       COEFc(3,0)= 2.0*SQRT(2.)/3.0
       COEFc(4,0)= FPIc
       COEFc(5,0)= 0.0
C
       COEFc(1,1)=-SQRT(2.)/3.0
       COEFc(2,1)= SQRT(2.)/3.0
       COEFc(3,1)= 0.0
       COEFc(4,1)= FPIc
       COEFc(5,1)= SQRT(2.)
C
       COEFc(1,2)=-SQRT(2.)/3.0
       COEFc(2,2)= SQRT(2.)/3.0
       COEFc(3,2)= 0.0
       COEFc(4,2)= 0.0
       COEFc(5,2)=-SQRT(2.)
C AJW 11/97: Add in the K*-prim-s, ala Finkemeier&Mirkes
       COEFc(1,3)= 1./3.
       COEFc(2,3)=-2./3.
       COEFc(3,3)= 2./3.
       COEFc(4,3)= 0.0
       COEFc(5,3)= 0.0
C
       COEFc(1,4)= 1.0/SQRT(2.)/3.0
       COEFc(2,4)=-1.0/SQRT(2.)/3.0
       COEFc(3,4)= 0.0
       COEFc(4,4)= 0.0
       COEFc(5,4)= 0.0
C
       COEFc(1,5)=-SQRT(2.)/3.0
       COEFc(2,5)= SQRT(2.)/3.0
       COEFc(3,5)= 0.0
       COEFc(4,5)= 0.0
       COEFc(5,5)=-SQRT(2.)
C
C AJW 11/97: Add in the K*-prim-s, ala Finkemeier&Mirkes
       COEFc(1,6)= 1./3.
       COEFc(2,6)=-2./3.
       COEFc(3,6)= 2./3.
       COEFc(4,6)= 0.0
       COEFc(5,6)=-2.0
C
       COEFc(1,7)= 0.0
       COEFc(2,7)= 0.0
       COEFc(3,7)= 0.0
       COEFc(4,7)= 0.0
       COEFc(5,7)=-SQRT(2.0/3.0)
      ENDIF
      IF (IVER.EQ.0.OR.J.NE.0) THEN   ! so far rchl only for 3pi modes
       COEF=COEFc(I,J)
      ELSEIF (IVER.EQ.1) THEN
       COEF=COEFrr(I,J)
      ELSE
       CALL TAUOLA_PRINT('wrong IVER='//CHAR(0))
!      write(*,*) 'wrong IVER=',IVER
       CALL TAUOLA_STOP()
      ENDIF
      END
      SUBROUTINE INIRChLget(I)
C routine to get version no for the currents physics initialization
C IVER=0  TAUOLA cleo
C IVER=1  TAUOLA RChL
      COMMON /IPChT/ IVER
      INTEGER        IVER
!
      BIND(C,NAME='tauola_ipcht') /IPChT/
!
      I=IVER
      end
      SUBROUTINE DADMAA(MODE,ISGN,HHV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
C ----------------------------------------------------------------------
* A1 DECAY UNWEIGHTED EVENTS
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL  HHV(4)
      REAL  HV(4),PAA(4),PNU(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4),PDUM5(4)
      REAL*4 RRR(3)
      REAL*8 SWT, SSWT
      DATA PI /3.141592653589793238462643/
      DATA IWARM/0/
C
      COMMON /TAUOLA_WEIGHT_DADMAA/ WTMAX
      REAL*4 WTMAX
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_inout') /INOUT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
      BIND(C,NAME='tauola_weight_dadmaa') /TAUOLA_WEIGHT_DADMAA/
!
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        NEVRAW=0
        NEVACC=0
        NEVOVR=0
        SWT=0
        SSWT=0
        WTMAX=1E-20
        DO 15 I=1,500
        CALL DPHSAA(WT,HV,PDUM1,PDUM2,PDUM3,PDUM4,PDUM5,JAA)
        IF(WT.GT.WTMAX/1.2) WTMAX=WT*1.2
15      CONTINUE
CC      CALL HBOOK1(801,'WEIGHT DISTRIBUTION  DADMAA    $',100,0,2)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
300     CONTINUE
        IF(IWARM.EQ.0) GOTO 902
        CALL DPHSAA(WT,HV,PNU,PAA,PIM1,PIM2,PIPL,JAA)
CC      CALL HFILL(801,WT/WTMAX)
        NEVRAW=NEVRAW+1
        SWT=SWT+WT
ccM.S.>>>>>>
cc        SSWT=SSWT+WT**2
        SSWT=SSWT+dble(WT)**2
ccM.S.<<<<<<
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX) NEVOVR=NEVOVR+1
        IF(RN*WTMAX.GT.WT) GOTO 300
C ROTATIONS TO BASIC TAU REST FRAME
        COSTHE=-1.+2.*RRR(2)
        THET=ACOS(COSTHE)
        PHI =2*PI*RRR(3)
        CALL ROTPOL(THET,PHI,PNU)
        CALL ROTPOL(THET,PHI,PAA)
        CALL ROTPOL(THET,PHI,PIM1)
        CALL ROTPOL(THET,PHI,PIM2)
        CALL ROTPOL(THET,PHI,PIPL)
        CALL ROTPOL(THET,PHI,HV)
        DO 44 I=1,3
 44     HHV(I)=-ISGN*HV(I)
        NEVACC=NEVACC+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVRAW.EQ.0) RETURN
        PARGAM=SWT/FLOAT(NEVRAW+1)
        ERROR=0
        IF(NEVRAW.NE.0) ERROR=SQRT(SSWT/SWT**2-1./FLOAT(NEVRAW))
        RAT=PARGAM/GAMEL
        CALL TAUOLA_PRINT('tauola.f:2609: (suppressed)'//CHAR(0))
!       WRITE(IOUT, 7010) NEVRAW,NEVACC,NEVOVR,PARGAM,RAT,ERROR
CC      CALL HPRINT(801)
        GAMPMC(5)=RAT
        GAMPER(5)=ERROR
CAM     NEVDEC(5)=NEVACC
      ENDIF
C     =====
      RETURN
!7003 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADMAA INITIALISATION ********',9X,1H*
!    $ /,' *',E20.5,5X,'WTMAX  = MAXIMUM WEIGHT                ',9X,1H*
!    $  /,1X,15(5H*****)/)
!7010 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADMAA FINAL REPORT  ******** ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVRAW = NO. OF A1  DECAYS TOTAL       ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVACC = NO. OF A1   DECS. ACCEPTED    ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
!    $ /,' *',E20.5,5X,'PARTIAL WTDTH (A1  DECAY) IN GEV UNITS ',9X,1H*
!    $ /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
!    $ /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
!    $  /,1X,15(5H*****)/)
 902  CALL TAUOLA_PRINT('----- DADMAA: LACK OF INITIALISATION'//CHAR(0))
!902  WRITE(IOUT, 9020)
!9020 FORMAT(' ----- DADMAA: LACK OF INITIALISATION')
      CALL TAUOLA_STOP()
      END
      SUBROUTINE DPHSAA(DGAMT,HV,PN,PAA,PIM1,PIM2,PIPL,JAA)
C ----------------------------------------------------------------------
* IT SIMULATES A1  DECAY IN TAU REST FRAME WITH
* Z-AXIS ALONG A1  MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      REAL  HV(4),PN(4),PAA(4),PIM1(4),PIM2(4),PIPL(4)
      REAL*4 RRR(1)
C MATRIX ELEMENT NUMBER:
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taukle') /TAUKLE/
!
      MNUM=0
C TYPE OF THE GENERATION:
      KEYT=1
      CALL RANMAR(RRR,1)
      RMOD=RRR(1)
      IF (RMOD.LT.BRA1) THEN
       JAA=1
       AMP1=AMPI
       AMP2=AMPI
       AMP3=AMPI
      ELSE
       JAA=2
       AMP1=AMPIZ
       AMP2=AMPIZ
       AMP3=AMPI
      ENDIF
       CALL CH3PISET(JAA)  ! information on sub-chanel passed for further use.
      CALL
     $   DPHTRE(DGAMT,HV,PN,PAA,PIM1,AMP1,PIM2,AMP2,PIPL,AMP3,KEYT,MNUM)
      END
      SUBROUTINE DADMKK(MODE,ISGN,HV,PKK,PNU)
C ----------------------------------------------------------------------
C FZ
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      REAL  PKK(4),PNU(4),HV(4)
      DATA PI /3.141592653589793238462643/
C
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_inout') /INOUT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
!
      IF(MODE.EQ.-1) THEN
C     ===================
        NEVTOT=0
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
        NEVTOT=NEVTOT+1
        EKK= (AMTAU**2+AMK**2-AMNUTA**2)/(2*AMTAU)
        ENU= (AMTAU**2-AMK**2+AMNUTA**2)/(2*AMTAU)
        XKK= SQRT(EKK**2-AMK**2)
C K MOMENTUM
        CALL SPHERA(XKK,PKK)
        PKK(4)=EKK
C TAU-NEUTRINO MOMENTUM
        DO 30 I=1,3
30      PNU(I)=-PKK(I)
        PNU(4)=ENU
        PXQ=AMTAU*EKK
        PXN=AMTAU*ENU
        QXN=PKK(4)*PNU(4)-PKK(1)*PNU(1)-PKK(2)*PNU(2)-PKK(3)*PNU(3)
        BRAK=(GV**2+GA**2)*(2*PXQ*QXN-AMK**2*PXN)
     &      +(GV**2-GA**2)*AMTAU*AMNUTA*AMK**2
        DO 40 I=1,3
40      HV(I)=-ISGN*2*GA*GV*AMTAU*(2*PKK(I)*QXN-PNU(I)*AMK**2)/BRAK
        HV(4)=1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVTOT.EQ.0) RETURN
        FKK=0.0354
CFZ THERE WAS BRAK/AMTAU**4 BEFORE
C        GAMM=(GFERMI*FKK)**2/(16.*PI)*AMTAU**3*
C     *       (BRAK/AMTAU**4)**2
CZW 7.02.93 here was an error affecting non standard model
C       configurations only
        GAMM=(GFERMI*FKK)**2/(16.*PI)*AMTAU**3*
     $       (BRAK/AMTAU**4)*
     $       SQRT((AMTAU**2-AMK**2-AMNUTA**2)**2
     $            -4*AMK**2*AMNUTA**2           )/AMTAU**2
        ERROR=0
        ERROR=0
        RAT=GAMM/GAMEL
        CALL TAUOLA_PRINT('tauola.f:2742: (suppressed)'//CHAR(0))
!       WRITE(IOUT, 7010) NEVTOT,GAMM,RAT,ERROR
        GAMPMC(6)=RAT
        GAMPER(6)=ERROR
CAM     NEVDEC(6)=NEVTOT
      ENDIF
C     =====
      RETURN
!7010 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADMKK FINAL REPORT   ********',9X,1H*
!    $ /,' *',I20  ,5X,'NEVTOT = NO. OF K  DECAYS TOTAL        ',9X,1H*,
!    $ /,' *',E20.5,5X,'PARTIAL WTDTH ( K DECAY) IN GEV UNITS  ',9X,1H*,
!    $ /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
!    $ /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH (STAT.)',9X,1H*
!    $  /,1X,15(5H*****)/)
      END
      SUBROUTINE DADMKS(MODE,ISGN,HHV,PNU,PKS,PKK,PPI,JKST)
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      COMMON / INOUT / INUT,IOUT
      REAL  HHV(4)
      REAL  HV(4),PKS(4),PNU(4),PKK(4),PPI(4)
      REAL  PDUM1(4),PDUM2(4),PDUM3(4),PDUM4(4)
      REAL*4 RRR(3),RMOD(1)
      REAL*8 SWT, SSWT
      DATA PI /3.141592653589793238462643/
      DATA IWARM/0/
C
      COMMON /TAUOLA_WEIGHT_DADMKS/ WTMAX
      REAL*4 WTMAX
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_inout') /INOUT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
      BIND(C,NAME='tauola_taukle') /TAUKLE/
      BIND(C,NAME='tauola_weight_dadmks') /TAUOLA_WEIGHT_DADMKS/
!
      IF(MODE.EQ.-1) THEN
C     ===================
        IWARM=1
        NEVRAW=0
        NEVACC=0
        NEVOVR=0
        SWT=0
        SSWT=0
        WTMAX=1E-20
        DO 15 I=1,5000
C THE INITIALISATION IS DONE WITH THE 66.7% MODE
        JKST=10
        CALL DPHSKS(WT,HV,PDUM1,PDUM2,PDUM3,PDUM4,JKST)
        IF(WT.GT.WTMAX/1.2) WTMAX=WT*1.2
15      CONTINUE
CC      CALL HBOOK1(801,'WEIGHT DISTRIBUTION  DADMKS    $',100,0,2)
CC      PRINT 7003,WTMAX
CC      CALL HBOOK1(112,'-------- K* MASS -------- $',100,0.,2.)
      ELSEIF(MODE.EQ. 0) THEN
C     =====================================
        IF(IWARM.EQ.0) GOTO 902
C  HERE WE CHOOSE RANDOMLY BETWEEN K0 PI+_ (66.7%)
C  AND K+_ PI0 (33.3%)
        DEC1=BRKS
400     CONTINUE
        CALL RANMAR(RMOD,1)
        IF(RMOD(1).LT.DEC1) THEN
          JKST=10
        ELSE
          JKST=20
        ENDIF
        CALL DPHSKS(WT,HV,PNU,PKS,PKK,PPI,JKST)
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX) NEVOVR=NEVOVR+1
        NEVRAW=NEVRAW+1
        SWT=SWT+WT
        SSWT=SSWT+WT**2
        IF(RN*WTMAX.GT.WT) GOTO 400
C ROTATIONS TO BASIC TAU REST FRAME
        COSTHE=-1.+2.*RRR(2)
        THET=ACOS(COSTHE)
        PHI =2*PI*RRR(3)
        CALL ROTOR2(THET,PNU,PNU)
        CALL ROTOR3( PHI,PNU,PNU)
        CALL ROTOR2(THET,PKS,PKS)
        CALL ROTOR3( PHI,PKS,PKS)
        CALL ROTOR2(THET,PKK,PKK)
        CALL ROTOR3(PHI,PKK,PKK)
        CALL ROTOR2(THET,PPI,PPI)
        CALL ROTOR3( PHI,PPI,PPI)
        CALL ROTOR2(THET,HV,HV)
        CALL ROTOR3( PHI,HV,HV)
        DO 44 I=1,3
 44     HHV(I)=-ISGN*HV(I)
        NEVACC=NEVACC+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        IF(NEVRAW.EQ.0) RETURN
        PARGAM=SWT/FLOAT(NEVRAW+1)
        ERROR=0
        IF(NEVRAW.NE.0) ERROR=SQRT(SSWT/SWT**2-1./FLOAT(NEVRAW))
        RAT=PARGAM/GAMEL
        CALL TAUOLA_PRINT('tauola.f:2856: (suppressed)'//CHAR(0))
!       WRITE(IOUT, 7010) NEVRAW,NEVACC,NEVOVR,PARGAM,RAT,ERROR
CC      CALL HPRINT(801)
        GAMPMC(7)=RAT
        GAMPER(7)=ERROR
CAM     NEVDEC(7)=NEVACC
      ENDIF
C     =====
      RETURN
!7003 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADMKS INITIALISATION ********',9X,1H*
!    $ /,' *',E20.5,5X,'WTMAX  = MAXIMUM WEIGHT                ',9X,1H*
!    $  /,1X,15(5H*****)/)
!7010 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADMKS FINAL REPORT   ********',9X,1H*
!    $ /,' *',I20  ,5X,'NEVRAW = NO. OF K* DECAYS TOTAL        ',9X,1H*,
!    $ /,' *',I20  ,5X,'NEVACC = NO. OF K*  DECS. ACCEPTED     ',9X,1H*,
!    $ /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
!    $ /,' *',E20.5,5X,'PARTIAL WTDTH (K* DECAY) IN GEV UNITS  ',9X,1H*,
!    $ /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
!    $ /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
!    $  /,1X,15(5H*****)/)
 902  CALL TAUOLA_PRINT('----- DADMKS: LACK OF INITIALISATION'//CHAR(0))
!902  WRITE(IOUT, 9020)
!9020 FORMAT(' ----- DADMKS: LACK OF INITIALISATION')
      CALL TAUOLA_STOP()
      END
      SUBROUTINE DPHSKS(DGAMT,HV,PN,PKS,PKK,PPI,JKST)
C ----------------------------------------------------------------------
C IT SIMULATES KAON* DECAY IN TAU REST FRAME WITH
C Z-AXIS ALONG KAON* MOMENTUM
C     JKST=10 FOR K* --->K0 + PI+-
C     JKST=20 FOR K* --->K+- + PI0
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      REAL  HV(4),PT(4),PN(4),PKS(4),PKK(4),PPI(4),QQ(4),RR1(1),RR2(1)
      DATA PI /3.141592653589793238462643/
C
      DATA ICONT /0/
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      PHSPAC=1./2**11/PI**5
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
      CALL RANMAR(RR1,1)
C HERE BEGIN THE K0,PI+_ DECAY
      IF(JKST.EQ.10)THEN
C     ==================
C MASS OF (REAL/VIRTUAL) K*
        AMS1=(AMPI+AMKZ)**2
        AMS2=(AMTAU-AMNUTA)**2
        ALP1=ATAN((AMS1-AMKST**2)/AMKST/GAMKST)
        ALP2=ATAN((AMS2-AMKST**2)/AMKST/GAMKST)
        CALL RANMAR(RR2,1)
        PROB1=0.2
        IF (RR2(1).LT.PROB1) THEN
C FLAT PHASE SPACE
         AMX2=AMS1+   RR1(1)*(AMS2-AMS1)
         AMX=SQRT(AMX2)
        ELSE
C PHASE SPACE WITH SAMPLING FOR K* RESONANCE
         ALP=ALP1+RR1(1)*(ALP2-ALP1)
         AMX2=AMKST**2+AMKST*GAMKST*TAN(ALP)
         AMX=SQRT(AMX2)
        ENDIF
C merging of the two channels
        PHSPAC1=(AMS2-AMS1)
        PHSPAC2=((AMX2-AMKST**2)**2+(AMKST*GAMKST)**2)
     &                /(AMKST*GAMKST)
        PHSPAC2=PHSPAC2*(ALP2-ALP1)
        A1=0.0
        A2=0.0
        IF (PHSPAC1.NE.0.0) A1=PROB1    /PHSPAC1
        IF (PHSPAC2.NE.0.0) A2=(1-PROB1)/PHSPAC2
        IF (A1+A2.NE.0.0) THEN
         PHSPAC=PHSPAC/(A1+A2)
        ELSE
         PHSPAC=0
        ENDIF
C
C TAU-NEUTRINO MOMENTUM
        PN(1)=0
        PN(2)=0
        PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX**2)
        PN(3)=-SQRT(ABS((PN(4)-AMNUTA)*(PN(4)+AMNUTA)))
C
C K* MOMENTUM
        PKS(1)=0
        PKS(2)=0
        PKS(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX**2)
        PKS(3)=-PN(3)
        PHSPAC=PHSPAC*(4*PI)*(2*PKS(3)/AMTAU)
C
CAM
        ENPI=( AMX**2+AMPI**2-AMKZ**2 ) / ( 2*AMX )
        PPPI=SQRT(abs(ENPI-AMPI)*(ENPI+AMPI))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMX)
C CHARGED PI MOMENTUM IN KAON* REST FRAME
        CALL SPHERA(PPPI,PPI)
        PPI(4)=ENPI
C NEUTRAL KAON MOMENTUM IN K* REST FRAME
        DO 20 I=1,3
20      PKK(I)=-PPI(I)
        PKK(4)=( AMX**2+AMKZ**2-AMPI**2 ) / ( 2*AMX )
        EXE=(PKS(4)+PKS(3))/AMX
C PION AND K  BOOSTED FROM K* REST FRAME TO TAU REST FRAME
        CALL BOSTR3(EXE,PPI,PPI)
        CALL BOSTR3(EXE,PKK,PKK)
        CALL DAM2PI(1,PT,PN,PPI,PKK,AMPLIT,HV)
        DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
C
C HERE BEGIN THE K+-,PI0 DECAY
      ELSEIF(JKST.EQ.20)THEN
C     ======================
C MASS OF (REAL/VIRTUAL) K*
        AMS1=(AMPIZ+AMK)**2
        AMS2=(AMTAU-AMNUTA)**2
        ALP1=ATAN((AMS1-AMKST**2)/AMKST/GAMKST)
        ALP2=ATAN((AMS2-AMKST**2)/AMKST/GAMKST)
        CALL RANMAR(RR2,1)
        PROB1=0.2
        IF (RR2(1).LT.PROB1) THEN
C FLAT PHASE SPACE
         AMX2=AMS1+   RR1(1)*(AMS2-AMS1)
         AMX=SQRT(AMX2)
        ELSE
C PHASE SPACE WITH SAMPLING FOR K* RESONANCE
         ALP=ALP1+RR1(1)*(ALP2-ALP1)
         AMX2=AMKST**2+AMKST*GAMKST*TAN(ALP)
         AMX=SQRT(AMX2)
        ENDIF
C merging of the two channels
        PHSPAC1=(AMS2-AMS1)
        PHSPAC2=((AMX2-AMKST**2)**2+(AMKST*GAMKST)**2)
     &                /(AMKST*GAMKST)
        PHSPAC2=PHSPAC2*(ALP2-ALP1)
        A1=0.0
        A2=0.0
        IF (PHSPAC1.NE.0.0) A1=PROB1    /PHSPAC1
        IF (PHSPAC2.NE.0.0) A2=(1-PROB1)/PHSPAC2
        IF (A1+A2.NE.0.0) THEN
         PHSPAC=PHSPAC/(A1+A2)
        ELSE
         PHSPAC=0
        ENDIF
C
C TAU-NEUTRINO MOMENTUM
        PN(1)=0
        PN(2)=0
        PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX**2)
        PN(3)=-SQRT(ABS((PN(4)-AMNUTA)*(PN(4)+AMNUTA)))
C KAON* MOMENTUM
        PKS(1)=0
        PKS(2)=0
        PKS(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX**2)
        PKS(3)=-PN(3)
        PHSPAC=PHSPAC*(4*PI)*(2*PKS(3)/AMTAU)
C
CAM
        ENPI=( AMX**2+AMPIZ**2-AMK**2 ) / ( 2*AMX )
        PPPI=SQRT(abs(ENPI-AMPIZ)*(ENPI+AMPIZ))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMX)
C NEUTRAL PI MOMENTUM IN K* REST FRAME
        CALL SPHERA(PPPI,PPI)
        PPI(4)=ENPI
C CHARGED KAON MOMENTUM IN K* REST FRAME
        DO 50 I=1,3
50      PKK(I)=-PPI(I)
        PKK(4)=( AMX**2+AMK**2-AMPIZ**2 ) / ( 2*AMX )
        EXE=(PKS(4)+PKS(3))/AMX
C PION AND K  BOOSTED FROM K* REST FRAME TO TAU REST FRAME
        CALL BOSTR3(EXE,PPI,PPI)
        CALL BOSTR3(EXE,PKK,PKK)
        CALL DAM2PI(2,PT,PN,PKK,PPI,AMPLIT,HV)
        DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
      ENDIF
      RETURN
      END
      SUBROUTINE DPHNPI(DGAMT,HVX,PNX,PRX,PPIX,JNPI)
C ----------------------------------------------------------------------
C IT SIMULATES MULTIPI DECAY IN TAU REST FRAME WITH
C Z-AXIS OPPOSITE TO NEUTRINO MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / TAUDCD /IDFFIN(9,NMODE),MULPIK(NMODE)
     &                ,NAMES
      CHARACTER NAMES(NMODE)*31
      REAL*8 WETMAX(20)
C
      REAL*8  PN(4),PR(4),PPI(4,9),HV(4)
      REAL*4  PNX(4),PRX(4),PPIX(4,9),HVX(4)
      REAL*8  PV(5,9),PT(4),UE(3),BE(3)
      REAL*8  PAWT,AMX,AMS1,AMS2,PA,PHS,PHSMAX,PMIN,PMAX
!!! M.S. to fix underflow >>>
      REAL*8  PHSPAC
!!! M.S. to fix underflow <<<
      REAL*8  GAM,BEP,PHI,A,B,C
      REAL*8  AMPIK
      REAL*4 RRR(9),RRX(2),RN(1),RR2(1)
C
      DATA PI /3.141592653589793238462643/
      DATA WETMAX /20*1D-15/
C
CC--      PAWT(A,B,C)=SQRT((A**2-(B+C)**2)*(A**2-(B-C)**2))/(2.*A)
C
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taudcd') /TAUDCD/
!
      PAWT(A,B,C)=
     $  SQRT(MAX(0.D0,(A**2-(B+C)**2)*(A**2-(B-C)**2)))/(2.D0*A)
C
      AMPIK(I,J)=DCDMAS(IDFFIN(I,J))
C
C
      IF ((JNPI.LE.0).OR.JNPI.GT.20) THEN
       CALL TAUOLA_PRINT(
     $   'JNPI OUTSIDE RANGE DEFINED BY WETMAX; JNPI='//CHAR(0))
!      WRITE(6,*) 'JNPI OUTSIDE RANGE DEFINED BY WETMAX; JNPI=',JNPI
       CALL TAUOLA_STOP()
      ENDIF
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C
 500  CONTINUE
C MASS OF VIRTUAL W
      ND=MULPIK(JNPI)
      PS=0.
      PHSPAC = 1./2.**5 /PI**2
      DO 4 I=1,ND
4     PS  =PS+AMPIK(I,JNPI)
      CALL RANMAR(RR2,1)
      AMS1=PS**2
      AMS2=(AMTAU-AMNUTA)**2
C
C
      AMX2=AMS1+   RR2(1)*(AMS2-AMS1)
      AMX =SQRT(AMX2)
      AMW =AMX
      PHSPAC=PHSPAC * (AMS2-AMS1)
C
C TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX2)
      PN(3)=-SQRT(ABS((PN(4)-AMNUTA)*(PN(4)+AMNUTA)))
C W MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX2)
      PR(3)=-PN(3)
      PHSPAC=PHSPAC * (4.*PI) * (2.*PR(3)/AMTAU)
C
C AMPLITUDE  (cf YS.Tsai Phys.Rev.D4,2821(1971)
C    or F.Gilman SH.Rhie Phys.Rev.D31,1066(1985)
C
        PXQ=AMTAU*PR(4)
        PXN=AMTAU*PN(4)
        QXN=PR(4)*PN(4)-PR(1)*PN(1)-PR(2)*PN(2)-PR(3)*PN(3)
C HERE WAS AN ERROR. 20.10.91 (ZW)
C       BRAK=2*(GV**2+GA**2)*(2*PXQ*PXN+AMX2*QXN)
        BRAK=2*(GV**2+GA**2)*(2*PXQ*QXN+AMX2*PXN)
     &      -6*(GV**2-GA**2)*AMTAU*AMNUTA*AMX2
CAM     Assume neutrino mass=0. and sum over final polarisation
C     BRAK= 2*(AMTAU**2-AMX2) * (AMTAU**2+2.*AMX2)
      AMPLIT=CCABIB**2*GFERMI**2/2. * BRAK * AMX2*SIGEE(AMX2,JNPI)
      DGAMT=1./(2.*AMTAU)*AMPLIT*PHSPAC
C
C   ISOTROPIC W DECAY IN W REST FRAME
      PHSMAX = 1.
      DO 200 I=1,4
  200 PV(I,1)=PR(I)
      PV(5,1)=AMW
      PV(5,ND)=AMPIK(ND,JNPI)
C    COMPUTE MAX. PHASE SPACE FACTOR
      PMAX=AMW-PS+AMPIK(ND,JNPI)
      PMIN=.0
      DO 220 IL=ND-1,1,-1
      PMAX=PMAX+AMPIK(IL,JNPI)
      PMIN=PMIN+AMPIK(IL+1,JNPI)
  220 PHSMAX=PHSMAX*PAWT(PMAX,PMIN,AMPIK(IL,JNPI))/PMAX
C --- 2.02.94 ZW  9 lines
      AMX=AMW
      DO 222 IL=1,ND-2
      AMS1=.0
      DO 223 JL=IL+1,ND
 223  AMS1=AMS1+AMPIK(JL,JNPI)
      AMS1=AMS1**2
      AMX =(AMX-AMPIK(IL,JNPI))
      AMS2=(AMX)**2
      PHSMAX=PHSMAX * (AMS2-AMS1)
 222  CONTINUE
      NCONT=0
  100 CONTINUE
      NCONT=NCONT+1
CAM  GENERATE ND-2 EFFECTIVE MASSES
      PHS=1.D0
      PHSPAC = 1./2.**(6*ND-7) /PI**(3*ND-4)
      AMX=AMW
      CALL RANMAR(RRR,ND-2)
      DO 230 IL=1,ND-2
      AMS1=.0D0
      DO 231 JL=IL+1,ND
  231 AMS1=AMS1+AMPIK(JL,JNPI)
      AMS1=AMS1**2
      AMS2=(AMX-AMPIK(IL,JNPI))**2
      RR1=RRR(IL)
      AMX2=AMS1+  RR1*(AMS2-AMS1)
      AMX=SQRT(AMX2)
      PV(5,IL+1)=AMX
      PHSPAC=PHSPAC * (AMS2-AMS1)
C ---  2.02.94 ZW 1 line
      PHS=PHS* (AMS2-AMS1)
      PA=PAWT(PV(5,IL),PV(5,IL+1),AMPIK(IL,JNPI))
      PHS   =PHS    *PA/PV(5,IL)
  230 CONTINUE
      PA=PAWT(PV(5,ND-1),AMPIK(ND-1,JNPI),AMPIK(ND,JNPI))
      PHS   =PHS    *PA/PV(5,ND-1)
      CALL RANMAR(RN,1)
      IF(PHSMAX.NE.0.0) THEN ! TP 5.10.2011 due to rounding errs.
                             ! PHSMAX may be zero, protect div. by it
        WETMAX(JNPI)=1.2D0*MAX(WETMAX(JNPI)/1.2D0,PHS/PHSMAX)
      ELSE
        WETMAX(JNPI)=1.2D0*WETMAX(JNPI)/1.2D0
      ENDIF
      IF (NCONT.EQ.500 000) THEN
          XNPI=0.0
          DO KK=1,ND
            XNPI=XNPI+AMPIK(KK,JNPI)
          ENDDO
       CALL TAUOLA_PRINT('ROUNDING INSTABILITY IN DPHNPI ?'//CHAR(0))
!      WRITE(6,*) 'ROUNDING INSTABILITY IN DPHNPI ?'
       CALL TAUOLA_PRINT('AMW=XNPI='//CHAR(0))
!      WRITE(6,*) 'AMW=',AMW,'XNPI=',XNPI
       CALL TAUOLA_PRINT(
     $   'IF =AMW= IS NEARLY EQUAL =XNPI= THAT IS IT'//CHAR(0))
!      WRITE(6,*) 'IF =AMW= IS NEARLY EQUAL =XNPI= THAT IS IT'
       CALL TAUOLA_PRINT('PHS=PHSMAX='//CHAR(0))
!      WRITE(6,*) 'PHS=',PHS,'PHSMAX=',PHSMAX
       GOTO 500
      ENDIF
      IF(RN(1)*PHSMAX*WETMAX(JNPI).GT.PHS) GO TO 100
C...PERFORM SUCCESSIVE TWO-PARTICLE DECAYS IN RESPECTIVE CM FRAME
  280 DO 300 IL=1,ND-1
      PA=PAWT(PV(5,IL),PV(5,IL+1),AMPIK(IL,JNPI))
      CALL RANMAR(RRX,2)
      UE(3)=2.*RRX(1)-1.
      PHI=2.*PI*RRX(2)
      UE(1)=SQRT(1.D0-UE(3)**2)*COS(PHI)
      UE(2)=SQRT(1.D0-UE(3)**2)*SIN(PHI)
      DO 290 J=1,3
      PPI(J,IL)=PA*UE(J)
  290 PV(J,IL+1)=-PA*UE(J)
      PPI(4,IL)=SQRT(PA**2+AMPIK(IL,JNPI)**2)
      PV(4,IL+1)=SQRT(PA**2+PV(5,IL+1)**2)
      PHSPAC=PHSPAC *(4.*PI)*(2.*PA/PV(5,IL))
  300 CONTINUE
C...LORENTZ TRANSFORM DECAY PRODUCTS TO TAU FRAME
      DO 310 J=1,4
  310 PPI(J,ND)=PV(J,ND)
      DO 340 IL=ND-1,1,-1
      DO 320 J=1,3
  320 BE(J)=PV(J,IL)/PV(4,IL)
      GAM=PV(4,IL)/PV(5,IL)
      DO 340 I=IL,ND
      BEP=BE(1)*PPI(1,I)+BE(2)*PPI(2,I)+BE(3)*PPI(3,I)
      DO 330 J=1,3
  330 PPI(J,I)=PPI(J,I)+GAM*(GAM*BEP/(1.D0+GAM)+PPI(4,I))*BE(J)
      PPI(4,I)=GAM*(PPI(4,I)+BEP)
  340 CONTINUE
C
            HV(4)=1.
            HV(3)=0.
            HV(2)=0.
            HV(1)=0.
      DO K=1,4
        PNX(K)=PN(K)
        PRX(K)=PR(K)
        HVX(K)=HV(K)
        DO L=1,ND
          PPIX(K,L)=PPI(K,L)
        ENDDO
      ENDDO
      RETURN
      END
      FUNCTION SIGEE(Q2,JNP)
C ----------------------------------------------------------------------
C  e+e- cross section in the (1.GEV2,AMTAU**2) region
C  normalised to sig0 = 4/3 pi alfa2
C  used in matrix element for multipion tau decays
C  cf YS.Tsai        Phys.Rev D4 ,2821(1971)
C     F.Gilman et al Phys.Rev D17,1846(1978)
C     C.Kiesling, to be pub. in High Energy e+e- Physics (1988)
C  DATSIG(*,1) = e+e- -> pi+pi-2pi0
C  DATSIG(*,2) = e+e- -> 2pi+2pi-
C  DATSIG(*,3) = 5-pion contribution (a la TN.Pham et al)
C                (Phys Lett 78B,623(1978)
C  DATSIG(*,5) = e+e- -> 6pi
C
C  4- and 6-pion cross sections from data
C  5-pion contribution related to 4-pion cross section
C
C     Called by DPHNPI
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
        REAL*4 DATSIG(17,6)
C
      DATA DATSIG/
     1  7.40,12.00,16.15,21.25,24.90,29.55,34.15,37.40,37.85,37.40,
     2 36.00,33.25,30.50,27.70,24.50,21.25,18.90,
     3  1.24, 2.50, 3.70, 5.40, 7.45,10.75,14.50,18.20,22.30,28.90,
     4 29.35,25.60,22.30,18.60,14.05,11.60, 9.10,
     5 17*.0,
     6 17*.0,
     7 9*.0,.65,1.25,2.20,3.15,5.00,5.75,7.80,8.25,
     8 17*.0/
      DATA SIG0 / 86.8 /
      DATA PI /3.141592653589793238462643/
      DATA INIT / 0 /
C
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
        JNPI=JNP
        IF(JNP.EQ.4) JNPI=3
        IF(JNP.EQ.3) JNPI=4
      IF(INIT.EQ.0) THEN
        INIT=1
C AJWMOD: initialize if called from outside QQ:
C        IF (AMPI.LT.0.139) AMPI = 0.1395675
        AMPI2=AMPI**2
        FPI = .943*AMPI
        DO 100 I=1,17
        DATSIG(I,2) = DATSIG(I,2)/2.
        DATSIG(I,1) = DATSIG(I,1) + DATSIG(I,2)
        S = 1.025+(I-1)*.05
        FACT=0.
        S2=S**2
        DO 200 J=1,17
        T= 1.025+(J-1)*.05
        IF(T . GT. S-AMPI ) GO TO 201
        T2=T**2
        FACT=(T2/S2)**2*SQRT((S2-T2-AMPI2)**2-4.*T2*AMPI2)/S2 *2.*T*.05
        FACT = FACT * (DATSIG(J,1)+DATSIG(J+1,1))
 200    DATSIG(I,3) = DATSIG(I,3) + FACT
 201    DATSIG(I,3) = DATSIG(I,3) /(2*PI*FPI)**2
        DATSIG(I,4) = DATSIG(I,3)
        DATSIG(I,6) = DATSIG(I,5)
 100    CONTINUE
C       WRITE(6,1000) DATSIG
!1000   FORMAT(///1X,' EE SIGMA USED IN MULTIPI DECAYS'/
!    %        (17F7.2/))
      ENDIF
      Q=SQRT(Q2)
      QMIN=1.
      IF(Q.LT.QMIN) THEN
        SIGEE=DATSIG(1,JNPI)+
     &       (DATSIG(2,JNPI)-DATSIG(1,JNPI))*(Q-1.)/.05
      ELSEIF(Q.LT.1.8) THEN
        DO 1 I=1,16
        QMAX = QMIN + .05
        IF(Q.LT.QMAX) GO TO 2
        QMIN = QMIN + .05
 1      CONTINUE
 2      SIGEE=DATSIG(I,JNPI)+
     &       (DATSIG(I+1,JNPI)-DATSIG(I,JNPI)) * (Q-QMIN)/.05
      ELSEIF(Q.GT.1.8) THEN
        SIGEE=DATSIG(17,JNPI)+
     &       (DATSIG(17,JNPI)-DATSIG(16,JNPI)) * (Q-1.8)/.05
      ENDIF
      IF(SIGEE.LT..0) SIGEE=0.
C
      SIGEE = SIGEE/(6.*PI**2*SIG0)
C
      RETURN
      END
      SUBROUTINE DPHSPK(DGAMT,HV,PN,PAA,PNPI,JAA)
C ----------------------------------------------------------------------
* IT SIMULATES THREE PI (K) DECAY IN THE TAU REST FRAME
* Z-AXIS ALONG HADRONIC SYSTEM
C ----------------------------------------------------------------------
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / TAUDCD /IDFFIN(9,NMODE),MULPIK(NMODE)
     &                ,NAMES
      CHARACTER NAMES(NMODE)*31
      REAL  HV(4),PN(4),PAA(4),PIM1(4),PIM2(4),PIPL(4),PNPI(4,9)
C MATRIX ELEMENT NUMBER:
!
      BIND(C,NAME='tauola_taudcd') /TAUDCD/
!
      MNUM=JAA
C TYPE OF THE GENERATION:
      KEYT=4
      IF(JAA.EQ.7) KEYT=3
C --- MASSES OF THE DECAY PRODUCTS
       AMP1=DCDMAS(IDFFIN(1,JAA+NM4+NM5+NM6))
       AMP2=DCDMAS(IDFFIN(2,JAA+NM4+NM5+NM6))
       AMP3=DCDMAS(IDFFIN(3,JAA+NM4+NM5+NM6))
      CALL
     $   DPHTRE(DGAMT,HV,PN,PAA,PIM1,AMP1,PIM2,AMP2,PIPL,AMP3,KEYT,MNUM)
            DO I=1,4
              PNPI(I,1)=PIM1(I)
              PNPI(I,2)=PIM2(I)
              PNPI(I,3)=PIPL(I)
            ENDDO
      END
      SUBROUTINE
     $   DPHTRE(DGAMT,HV,PN,PAA,PIM1,AMPA,PIM2,AMPB,PIPL,AMP3,KEYT,MNUM)
C ----------------------------------------------------------------------
* IT SIMULATES A1  DECAY IN TAU REST FRAME WITH
* Z-AXIS ALONG A1  MOMENTUM
* it can be also used to generate K K pi and K pi pi tau decays.
* INPUT PARAMETERS
* KEYT - algorithm controlling switch
*  2   - flat phase space PIM1 PIM2 symmetrized statistical factor 1/2
*  1   - like 1 but peaked around a1 and rho (two channels) masses.
*  3   - peaked around omega, all particles different
* other- flat phase space, all particles different
* AMP1 - mass of first pi, etc. (1-3)
* MNUM - matrix element type
*  0   - a1 matrix element
* 1-6  - matrix element for K pi pi, K K pi decay modes
*  7   - pi- pi0 gamma matrix element
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  HV(4),PT(4),PN(4),PAA(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PR(4)
      REAL*4 RRR(5)
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      XLAM(X,Y,Z)=SQRT(ABS((X-Y-Z)**2-4.0*Y*Z))
C AMRO, GAMRO IS ONLY A PARAMETER FOR GETING HIGHT EFFICIENCY
C
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1./2**17/PI**8
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C
      CALL RANMAR(RRR,5)
      RR=RRR(5)
C
      CALL CHOICE(MNUM,RR,ICHAN,PROB1,PROB2,PROB3,
     $            AMRX,GAMRX,AMRA,GAMRA,AMRB,GAMRB)
      IF     (ICHAN.EQ.1) THEN
        AMP1=AMPB
        AMP2=AMPA
      ELSEIF (ICHAN.EQ.2) THEN
        AMP1=AMPA
        AMP2=AMPB
      ELSE
        AMP1=AMPB
        AMP2=AMPA
      ENDIF
CAM
        RR1=RRR(1)
        AMS1=(AMP1+AMP2+AMP3)**2
        AMS2=(AMTAU-AMNUTA)**2
* PHASE SPACE WITH SAMPLING FOR A1  RESONANCE
        ALP1=ATAN((AMS1-AMRX**2)/AMRX/GAMRX)
        ALP2=ATAN((AMS2-AMRX**2)/AMRX/GAMRX)
        ALP=ALP1+RR1*(ALP2-ALP1)
        AM3SQ =AMRX**2+AMRX*GAMRX*TAN(ALP)
        AM3 =SQRT(AM3SQ)
        PHSPAC=PHSPAC*((AM3SQ-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)
        PHSPAC=PHSPAC*(ALP2-ALP1)
C MASS OF (REAL/VIRTUAL) RHO -
        RR2=RRR(2)
        AMS1=(AMP2+AMP3)**2
        AMS2=(AM3-AMP1)**2
      IF (ICHAN.LE.2) THEN
* PHASE SPACE WITH SAMPLING FOR RHO RESONANCE,
        ALP1=ATAN((AMS1-AMRA**2)/AMRA/GAMRA)
        ALP2=ATAN((AMS2-AMRA**2)/AMRA/GAMRA)
        ALP=ALP1+RR2*(ALP2-ALP1)
        AM2SQ =AMRA**2+AMRA*GAMRA*TAN(ALP)
        AM2 =SQRT(AM2SQ)
C --- THIS PART OF THE JACOBIAN WILL BE RECOVERED LATER ---------------
C     PHSPAC=PHSPAC*(ALP2-ALP1)
C     PHSPAC=PHSPAC*((AM2SQ-AMRA**2)**2+(AMRA*GAMRA)**2)/(AMRA*GAMRA)
C----------------------------------------------------------------------
      ELSE
* FLAT PHASE SPACE;
        AM2SQ=AMS1+   RR2*(AMS2-AMS1)
        AM2 =SQRT(AM2SQ)
        PHF0=(AMS2-AMS1)
      ENDIF
* RHO RESTFRAME, DEFINE PIPL AND PIM1
        ENQ1=(AM2SQ-AMP2**2+AMP3**2)/(2*AM2)
        ENQ2=(AM2SQ+AMP2**2-AMP3**2)/(2*AM2)
        PPI=         ENQ1**2-AMP3**2
        PPPI=SQRT(ABS(ENQ1**2-AMP3**2))
C --- this part of jacobian will be recovered later
        PHF1=(4*PI)*(2*PPPI/AM2)
* PI MINUS MOMENTUM IN RHO REST FRAME
        CALL SPHERA(PPPI,PIPL)
        PIPL(4)=ENQ1
* PI0 1 MOMENTUM IN RHO REST FRAME
        DO 30 I=1,3
 30     PIM1(I)=-PIPL(I)
        PIM1(4)=ENQ2
* A1 REST FRAME, DEFINE PIM2
*       RHO  MOMENTUM
        PR(1)=0
        PR(2)=0
        PR(4)=1./(2*AM3)*(AM3**2+AM2**2-AMP1**2)
        PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
        PPI  =          PR(4)**2-AM2**2
*       PI0 2 MOMENTUM
        PIM2(1)=0
        PIM2(2)=0
        PIM2(4)=1./(2*AM3)*(AM3**2-AM2**2+AMP1**2)
        PIM2(3)=-PR(3)
      PHF2=(4*PI)*(2*PR(3)/AM3)
* OLD PIONS BOOSTED FROM RHO REST FRAME TO A1 REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTR3(EXE,PIPL,PIPL)
      CALL BOSTR3(EXE,PIM1,PIM1)
      RR3=RRR(3)
      RR4=RRR(4)
CAM   THET =PI*RR3
      THET =ACOS(-1.+2*RR3)
      PHI = 2*PI*RR4
      CALL ROTPOL(THET,PHI,PIPL)
      CALL ROTPOL(THET,PHI,PIM1)
      CALL ROTPOL(THET,PHI,PIM2)
      CALL ROTPOL(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE A1 AND NEUTRINO MOMENTA
* A1  MOMENTUM
      PAA(1)=0
      PAA(2)=0
      PAA(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AM3**2)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM3**2))
      PPI   =          PAA(4)**2-AM3**2
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
* TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AM3**2)
      PN(3)=-PAA(3)
C HERE WE CORRECT FOR THE JACOBIANS OF THE TWO CHAINS
C ---FIRST CHANNEL ------- PIM1+PIPL
        AMS1=(AMP2+AMP3)**2
        AMS2=(AM3-AMP1)**2
        ALP1=ATAN((AMS1-AMRA**2)/AMRA/GAMRA)
        ALP2=ATAN((AMS2-AMRA**2)/AMRA/GAMRA)
       XPRO =      (PIM1(3)+PIPL(3))**2
     $            +(PIM1(2)+PIPL(2))**2+(PIM1(1)+PIPL(1))**2
       AM2SQ=-XPRO+(PIM1(4)+PIPL(4))**2
C JACOBIAN OF SPEEDING
       FF1   =       ((AM2SQ-AMRA**2)**2+(AMRA*GAMRA)**2)/(AMRA*GAMRA)
       FF1   =FF1     *(ALP2-ALP1)
C LAMBDA OF RHO DECAY
       GG1   =       (4*PI)*(XLAM(AM2SQ,AMP2**2,AMP3**2)/AM2SQ)
C LAMBDA OF A1 DECAY
       GG1   =GG1   *(4*PI)*SQRT(4*XPRO/AM3SQ)
       XJAJE=GG1*(AMS2-AMS1)
C ---SECOND CHANNEL ------ PIM2+PIPL
       AMS1=(AMP1+AMP3)**2
       AMS2=(AM3-AMP2)**2
        ALP1=ATAN((AMS1-AMRB**2)/AMRB/GAMRB)
        ALP2=ATAN((AMS2-AMRB**2)/AMRB/GAMRB)
       XPRO =      (PIM2(3)+PIPL(3))**2
     $            +(PIM2(2)+PIPL(2))**2+(PIM2(1)+PIPL(1))**2
       AM2SQ=-XPRO+(PIM2(4)+PIPL(4))**2
       FF2   =       ((AM2SQ-AMRB**2)**2+(AMRB*GAMRB)**2)/(AMRB*GAMRB)
       FF2   =FF2     *(ALP2-ALP1)
       GG2   =       (4*PI)*(XLAM(AM2SQ,AMP1**2,AMP3**2)/AM2SQ)
       GG2   =GG2   *(4*PI)*SQRT(4*XPRO/AM3SQ)
       XJADW=GG2*(AMS2-AMS1)
C
       A1=0.0
       A2=0.0
       A3=0.0
       XJAC1=FF1*GG1
       XJAC2=FF2*GG2
       IF (ICHAN.EQ.2) THEN
         XJAC3=XJADW
       ELSE
         XJAC3=XJAJE
       ENDIF
       IF (XJAC1.NE.0.0) A1=PROB1/XJAC1
       IF (XJAC2.NE.0.0) A2=PROB2/XJAC2
       IF (XJAC3.NE.0.0) A3=PROB3/XJAC3
C
       IF (A1+A2+A3.NE.0.0) THEN
         PHSPAC=PHSPAC/(A1+A2+A3)
       ELSE
         PHSPAC=0.0
       ENDIF
       IF(ICHAN.EQ.2) THEN
        DO 70 I=1,4
        X=PIM1(I)
        PIM1(I)=PIM2(I)
 70     PIM2(I)=X
       ENDIF
* ALL PIONS BOOSTED FROM A1  REST FRAME TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO NEUTRINO MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM3
      CALL BOSTR3(EXE,PIPL,PIPL)
      CALL BOSTR3(EXE,PIM1,PIM1)
      CALL BOSTR3(EXE,PIM2,PIM2)
      CALL BOSTR3(EXE,PR,PR)
C PARTIAL WIDTH CONSISTS OF PHASE SPACE AND AMPLITUDE
      IF (MNUM.EQ.8) THEN
        CALL DAMPOG(PT,PN,PIM1,PIM2,PIPL,AMPLIT,HV)
C      ELSEIF (MNUM.EQ.0) THEN
C        CALL DAMPAA(PT,PN,PIM1,PIM2,PIPL,AMPLIT,HV)
      ELSE
        CALL DAMPPK(MNUM,PT,PN,PIM1,PIM2,PIPL,AMPLIT,HV)
      ENDIF
      IF (KEYT.EQ.1.OR.KEYT.EQ.2) THEN
C THE STATISTICAL FACTOR FOR IDENTICAL PI-S IS CANCELLED WITH
C TWO, FOR TWO MODES OF A1 DECAY NAMELLY PI+PI-PI- AND PI-PI0PI0
        PHSPAC=PHSPAC*2.0
        PHSPAC=PHSPAC/2.
      ENDIF
      DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
      END
      FUNCTION GFUN(QKWA)
C ****************************************************************
C     G-FUNCTION USED TO INRODUCE ENERGY DEPENDENCE IN A1 WIDTH
C ****************************************************************
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
       IF (QKWA.LT.(AMRO+AMPI)**2) THEN
          GFUN=4.1*(QKWA-9*AMPIZ**2)**3
     $        *(1.-3.3*(QKWA-9*AMPIZ**2)+5.8*(QKWA-9*AMPIZ**2)**2)
       ELSE
          GFUN=QKWA*(1.623+10.38/QKWA-9.32/QKWA**2+0.65/QKWA**3)
       ENDIF
      END
      COMPLEX FUNCTION BWIGS(S,M,G)
C **********************************************************
C     P-WAVE BREIT-WIGNER  FOR K*
C **********************************************************
      REAL S,M,G
      REAL PI,PIM,QS,QM,W,GS,MK
C AJW: add K*-prim possibility:
      REAL PM, PG, PBETA
      COMPLEX BW,BWP
      DATA INIT /0/
      P(A,B,C)=SQRT(ABS(ABS(((A+B-C)**2-4.*A*B)/4./A)
     $                    +(((A+B-C)**2-4.*A*B)/4./A))/2.0)
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0) THEN
      INIT=1
      PI=3.141592654
      PIM=.139
      MK=.493667
C AJW: add K*-prim possibility:
      PM = PKORB(1,16)
      PG = PKORB(2,16)
      PBETA = PKORB(3,16)
C -------  BREIT-WIGNER -----------------------
         ENDIF
         QS=P(S,PIM**2,MK**2)
         QM=P(M**2,PIM**2,MK**2)
         W=SQRT(S)
         GS=G*(M/W)*(QS/QM)**3
         BW=M**2/CMPLX(M**2-S,-M*GS)
         QPM=P(PM**2,PIM**2,MK**2)
         G1=PG*(PM/W)*(QS/QPM)**3
         BWP=PM**2/CMPLX(PM**2-S,-PM*G1)
         BWIGS= (BW+PBETA*BWP)/(1+PBETA)
      RETURN
      END
      COMPLEX FUNCTION BWIG(S,M,G)
C **********************************************************
C     P-WAVE BREIT-WIGNER  FOR RHO
C **********************************************************
      REAL S,M,G
      REAL PI,PIM,QS,QM,W,GS
      DATA INIT /0/
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0) THEN
      INIT=1
      PI=3.141592654
      PIM=.139
C -------  BREIT-WIGNER -----------------------
         ENDIF
       IF (S.GT.4.*PIM**2) THEN
         QS=SQRT(ABS(ABS(S/4.-PIM**2)+(S/4.-PIM**2))/2.0)
         QM=SQRT(M**2/4.-PIM**2)
         W=SQRT(S)
         GS=G*(M/W)*(QS/QM)**3
       ELSE
         GS=0.0
       ENDIF
         BWIG=M**2/CMPLX(M**2-S,-M*GS)
      RETURN
      END
      COMPLEX FUNCTION FPIK(W)
C **********************************************************
C     PION FORM FACTOR
C **********************************************************
      REAL ROM,ROG,ROM1,ROG1,BETA1,PI,PIM,S,W
      DATA  INIT /0/
C
C ------------ PARAMETERS --------------------
      IF (INIT.EQ.0 ) THEN
      INIT=1
      PI=3.141592654
      PIM=.140
      ROM=PKORB(1,9)
      ROG=PKORB(2,9)
      ROM1=PKORB(1,15)
      ROG1=PKORB(2,15)
      BETA1=PKORB(3,15)
      ENDIF
C -----------------------------------------------
      S=W**2
      FPIK= (BWIG(S,ROM,ROG)+BETA1*BWIG(S,ROM1,ROG1))
     & /(1+BETA1)
      RETURN
      END
      FUNCTION FPIRHO(W)
C **********************************************************
C     SQUARE OF PION FORM FACTOR
C **********************************************************
      FPIRHO=CABS(FPIK(W))**2
      END
      SUBROUTINE CLVEC(HJ,PN,PIV)
C ----------------------------------------------------------------------
* CALCULATES THE "VECTOR TYPE"  PI-VECTOR  PIV
* NOTE THAT THE NEUTRINO MOM. PN IS ASSUMED TO BE ALONG Z-AXIS
C
C     called by : DAMPAA
C ----------------------------------------------------------------------
      REAL PIV(4),PN(4)
      COMPLEX HJ(4),HN
C
      HN= HJ(4)*CMPLX(PN(4))-HJ(3)*CMPLX(PN(3))
      HH= REAL(HJ(4)*CONJG(HJ(4))-HJ(3)*CONJG(HJ(3))
     $        -HJ(2)*CONJG(HJ(2))-HJ(1)*CONJG(HJ(1)))
      DO 10 I=1,4
   10 PIV(I)=4.*REAL(HN*CONJG(HJ(I)))-2.*HH*PN(I)
      RETURN
      END
      SUBROUTINE CLAXI(HJ,PN,PIA)
C ----------------------------------------------------------------------
* CALCULATES THE "AXIAL TYPE"  PI-VECTOR  PIA
* NOTE THAT THE NEUTRINO MOM. PN IS ASSUMED TO BE ALONG Z-AXIS
C SIGN is chosen +/- for decay of TAU +/- respectively
C     called by : DAMPAA, CLNUT
C ----------------------------------------------------------------------
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / IDFC  / IDFF
      REAL PIA(4),PN(4)
      COMPLEX HJ(4),HJC(4)
C     DET2(I,J)=AIMAG(HJ(I)*HJC(J)-HJ(J)*HJC(I))
C -- here was an error (ZW, 21.11.1991)
!
      BIND(C,NAME='tauola_idfc') /IDFC/
      BIND(C,NAME='tauola_jaki') /JAKI/
!
      DET2(I,J)=AIMAG(HJC(I)*HJ(J)-HJC(J)*HJ(I))
C -- it was affecting sign of A_LR asymmetry in a1 decay.
C -- note also collision of notation of gamma_va as defined in
C -- TAUOLA paper and J.H. Kuhn and Santamaria Z. Phys C 48 (1990) 445
* -----------------------------------
      IF     (KTOM.EQ.1.OR.KTOM.EQ.-1) THEN
        SIGN= IDFF/ABS(IDFF)
      ELSEIF (KTOM.EQ.2) THEN
        SIGN=-IDFF/ABS(IDFF)
      ELSE
        CALL TAUOLA_PRINT('STOP IN CLAXI: KTOM='//CHAR(0))
!       PRINT *, 'STOP IN CLAXI: KTOM=',KTOM
        CALL TAUOLA_STOP()
      ENDIF
C
      DO 10 I=1,4
 10   HJC(I)=CONJG(HJ(I))
      PIA(1)= -2.*PN(3)*DET2(2,4)+2.*PN(4)*DET2(2,3)
      PIA(2)= -2.*PN(4)*DET2(1,3)+2.*PN(3)*DET2(1,4)
      PIA(3)=  2.*PN(4)*DET2(1,2)
      PIA(4)=  2.*PN(3)*DET2(1,2)
C ALL FOUR INDICES ARE UP SO  PIA(3) AND PIA(4) HAVE SAME SIGN
      DO 20 I=1,4
  20  PIA(I)=PIA(I)*SIGN
      END
      SUBROUTINE CLNUT(HJ,B,HV)
C ----------------------------------------------------------------------
* CALCULATES THE CONTRIBUTION BY NEUTRINO MASS
* NOTE THE TAU IS ASSUMED TO BE AT REST
C
C     called by : DAMPAA
C ----------------------------------------------------------------------
      COMPLEX HJ(4)
      REAL HV(4),P(4)
      DATA P /3*0.,1.0/
C
      CALL CLAXI(HJ,P,HV)
      B=REAL( HJ(4)*AIMAG(HJ(4)) - HJ(3)*AIMAG(HJ(3))
     &      - HJ(2)*AIMAG(HJ(2)) - HJ(1)*AIMAG(HJ(1))  )
      RETURN
      END
      SUBROUTINE DAMPOG(PT,PN,PIM1,PIM2,PIPL,AMPLIT,HV)
C ----------------------------------------------------------------------
* CALCULATES DIFFERENTIAL CROSS SECTION AND POLARIMETER VECTOR
* FOR TAU DECAY INTO A1, A1 DECAYS NEXT INTO RHO+PI AND RHO INTO PI+PI.
* ALL SPIN EFFECTS IN THE FULL DECAY CHAIN ARE TAKEN INTO ACCOUNT.
* CALCULATIONS DONE IN TAU REST FRAME WITH Z-AXIS ALONG NEUTRINO MOMENT
* THE ROUTINE IS WRITEN FOR ZERO NEUTRINO MASS.
C
C     called by : DPHSAA
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON /TESTA1/ KEYA1
      REAL  HV(4),PT(4),PN(4),PIM1(4),PIM2(4),PIPL(4)
      REAL  PAA(4),VEC1(4),VEC2(4)
      REAL  PIVEC(4),PIAKS(4),HVM(4)
      COMPLEX BWIGN,HADCUR(4),FNORM
      DATA ICONT /1/
* THIS INLINE FUNCT. CALCULATES THE SCALAR PART OF THE PROPAGATOR
C AJWMOD to satisfy compiler, comment out this unused function.
C
* FOUR MOMENTUM OF A1
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_testa1') /TESTA1/
!
      DO 10 I=1,4
      VEC1(I)=0.0
      VEC2(I)=0.0
      HV(I)  =0.0
   10 PAA(I)=PIM1(I)+PIM2(I)+PIPL(I)
      VEC1(1)=1.0
* MASSES OF A1, AND OF TWO PI-PAIRS WHICH MAY FORM RHO
      XMAA   =SQRT(ABS(PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2))
      XMOM   =SQRT(ABS( (PIM2(4)+PIPL(4))**2-(PIM2(3)+PIPL(3))**2
     $                 -(PIM2(2)+PIPL(2))**2-(PIM2(1)+PIPL(1))**2   ))
      XMRO2  =(PIPL(1))**2 +(PIPL(2))**2 +(PIPL(3))**2
* ELEMENTS OF HADRON CURRENT
      PROD1  =VEC1(1)*PIPL(1)
      PROD2  =VEC2(2)*PIPL(2)
      P12    =PIM1(4)*PIM2(4)-PIM1(1)*PIM2(1)
     $       -PIM1(2)*PIM2(2)-PIM1(3)*PIM2(3)
      P1PL   =PIM1(4)*PIPL(4)-PIM1(1)*PIPL(1)
     $       -PIM1(2)*PIPL(2)-PIM1(3)*PIPL(3)
      P2PL   =PIPL(4)*PIM2(4)-PIPL(1)*PIM2(1)
     $       -PIPL(2)*PIM2(2)-PIPL(3)*PIM2(3)
      DO 40 I=1,3
        VEC1(I)= (VEC1(I)-PROD1/XMRO2*PIPL(I))
 40   CONTINUE
        GNORM=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
      DO 41 I=1,3
        VEC1(I)= VEC1(I)/GNORM
 41   CONTINUE
      VEC2(1)=(VEC1(2)*PIPL(3)-VEC1(3)*PIPL(2))/SQRT(XMRO2)
      VEC2(2)=(VEC1(3)*PIPL(1)-VEC1(1)*PIPL(3))/SQRT(XMRO2)
      VEC2(3)=(VEC1(1)*PIPL(2)-VEC1(2)*PIPL(1))/SQRT(XMRO2)
      P1VEC1   =PIM1(4)*VEC1(4)-PIM1(1)*VEC1(1)
     $         -PIM1(2)*VEC1(2)-PIM1(3)*VEC1(3)
      P2VEC1   =VEC1(4)*PIM2(4)-VEC1(1)*PIM2(1)
     $         -VEC1(2)*PIM2(2)-VEC1(3)*PIM2(3)
      P1VEC2   =PIM1(4)*VEC2(4)-PIM1(1)*VEC2(1)
     $         -PIM1(2)*VEC2(2)-PIM1(3)*VEC2(3)
      P2VEC2   =VEC2(4)*PIM2(4)-VEC2(1)*PIM2(1)
     $         -VEC2(2)*PIM2(2)-VEC2(3)*PIM2(3)
* HADRON CURRENT
      FNORM=FORMOM(XMAA,XMOM)
      BRAK=0.0
      DO 120 JJ=1,2
        DO 45 I=1,4
       IF (JJ.EQ.1) THEN
        HADCUR(I) = FNORM *(
     $             VEC1(I)*(AMPI**2*P1PL-P2PL*(P12-P1PL))
     $            -PIM2(I)*(P2VEC1*P1PL-P1VEC1*P2PL)
     $            +PIPL(I)*(P2VEC1*P12 -P1VEC1*(AMPI**2+P2PL))  )
       ELSE
        HADCUR(I) = FNORM *(
     $             VEC2(I)*(AMPI**2*P1PL-P2PL*(P12-P1PL))
     $            -PIM2(I)*(P2VEC2*P1PL-P1VEC2*P2PL)
     $            +PIPL(I)*(P2VEC2*P12 -P1VEC2*(AMPI**2+P2PL))  )
       ENDIF
 45     CONTINUE
C
* CALCULATE PI-VECTORS: VECTOR AND AXIAL
      CALL CLVEC(HADCUR,PN,PIVEC)
      CALL CLAXI(HADCUR,PN,PIAKS)
      CALL CLNUT(HADCUR,BRAKM,HVM)
* SPIN INDEPENDENT PART OF DECAY DIFF-CROSS-SECT. IN TAU REST  FRAME
      BRAK=BRAK+(GV**2+GA**2)*PT(4)*PIVEC(4) +2.*GV*GA*PT(4)*PIAKS(4)
     &         +2.*(GV**2-GA**2)*AMNUTA*AMTAU*BRAKM
      DO 90 I=1,3
      HV(I)=HV(I)-(AMTAU*((GV**2+GA**2)*PIAKS(I)+2.*GV*GA*PIVEC(I)))
     &      +(GV**2-GA**2)*AMNUTA*AMTAU*HVM(I)
  90  CONTINUE
C HV IS DEFINED FOR TAU-    WITH GAMMA=B+HV*POL
 120  CONTINUE
      AMPLIT=(GFERMI*CCABIB)**2*BRAK/2.
C THE STATISTICAL FACTOR FOR IDENTICAL PI-S WAS CANCELLED WITH
C TWO, FOR TWO MODES OF A1 DECAY NAMELLY PI+PI-PI- AND PI-PI0PI0
C POLARIMETER VECTOR IN TAU REST FRAME
      DO 91 I=1,3
      HV(I)=-HV(I)/BRAK
 91   CONTINUE
      END
      SUBROUTINE DAMPPK(MNUM,PT,PN,PIM1,PIM2,PIM3,AMPLIT,HV)
C ----------------------------------------------------------------------
* CALCULATES DIFFERENTIAL CROSS SECTION AND POLARIMETER VECTOR
* FOR TAU DECAY INTO K K pi, K pi pi.
* ALL SPIN EFFECTS IN THE FULL DECAY CHAIN ARE TAKEN INTO ACCOUNT.
* CALCULATIONS DONE IN TAU REST FRAME WITH Z-AXIS ALONG NEUTRINO MOMENT
C MNUM DECAY MODE IDENTIFIER.
C
C     called by : DPHSAA
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON /IPChT/ IVER
      INTEGER        IVER
      REAL  HV(4),PT(4),PN(4),PIM1(4),PIM2(4),PIM3(4)
      REAL  PAA(4),VEC1(4),VEC2(4),VEC3(4),VEC4(4),VEC5(4)
      REAL  PIVEC(4),PIAKS(4),HVM(4)
      REAL FNORM(0:7)
      COMPLEX HADCUR(4),UROJ
      COMPLEX F1,F2,F3,F4,F5
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
CC
      DATA  FPIc /93.3E-3/
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_ipcht') /IPChT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      IF (ICONT.EQ.0) THEN
       ICONT=1
       UROJ=CMPLX(0.0,1.0)
       DWAPI0=SQRT(2.0)
      ENDIF
       IF (IVER.EQ.0.OR.MNUM.NE.0) THEN ! so far rchl only for 3pi modes
         FPI=FPIc
        ELSEIF (IVER.EQ.1) THEN
         FPI=GETFPIRPT(1) ! GET  defined in in ffwid3pi.f  of RChL-currents
        ELSE
         CALL TAUOLA_PRINT('wrong IVER='//CHAR(0))
!        write(*,*) 'wrong IVER=',IVER
         CALL TAUOLA_STOP()
        ENDIF
       FNORM(0)=CCABIB/FPI
       FNORM(1)=CCABIB/FPI
       FNORM(2)=CCABIB/FPI
       FNORM(3)=CCABIB/FPI
       FNORM(4)=SCABIB/FPI/DWAPI0
       FNORM(5)=SCABIB/FPI
       FNORM(6)=SCABIB/FPI
       FNORM(7)=CCABIB/FPI
C
      DO 10 I=1,4
   10 PAA(I)=PIM1(I)+PIM2(I)+PIM3(I)
      XMAA   =SQRT(ABS(PAA(4)**2-PAA(3)**2-PAA(2)**2-PAA(1)**2))
      XMRO1  =SQRT(ABS((PIM3(4)+PIM2(4))**2-(PIM3(1)+PIM2(1))**2
     $                -(PIM3(2)+PIM2(2))**2-(PIM3(3)+PIM2(3))**2))
      XMRO2  =SQRT(ABS((PIM3(4)+PIM1(4))**2-(PIM3(1)+PIM1(1))**2
     $                -(PIM3(2)+PIM1(2))**2-(PIM3(3)+PIM1(3))**2))
      XMRO3  =SQRT(ABS((PIM1(4)+PIM2(4))**2-(PIM1(1)+PIM2(1))**2
     $                -(PIM1(2)+PIM2(2))**2-(PIM1(3)+PIM2(3))**2))
* ELEMENTS OF HADRON CURRENT
      PROD1  =PAA(4)*(PIM2(4)-PIM3(4))-PAA(1)*(PIM2(1)-PIM3(1))
     $       -PAA(2)*(PIM2(2)-PIM3(2))-PAA(3)*(PIM2(3)-PIM3(3))
      PROD2  =PAA(4)*(PIM3(4)-PIM1(4))-PAA(1)*(PIM3(1)-PIM1(1))
     $       -PAA(2)*(PIM3(2)-PIM1(2))-PAA(3)*(PIM3(3)-PIM1(3))
      PROD3  =PAA(4)*(PIM1(4)-PIM2(4))-PAA(1)*(PIM1(1)-PIM2(1))
     $       -PAA(2)*(PIM1(2)-PIM2(2))-PAA(3)*(PIM1(3)-PIM2(3))
      DO 40 I=1,4
      VEC1(I)= PIM2(I)-PIM3(I) -PAA(I)*PROD1/XMAA**2
      VEC2(I)= PIM3(I)-PIM1(I) -PAA(I)*PROD2/XMAA**2
      VEC3(I)= PIM1(I)-PIM2(I) -PAA(I)*PROD3/XMAA**2
 40   VEC4(I)= PIM1(I)+PIM2(I)+PIM3(I)
      CALL PROD5(PIM1,PIM2,PIM3,VEC5)
* HADRON CURRENT
C be aware that sign of vec2 is opposite to sign of vec1 in a1 case
C Rationalize this code:
      F1 = CMPLX(COEF(1,MNUM))*FORM1(MNUM,XMAA**2,XMRO1**2,XMRO2**2)
      F2 = CMPLX(COEF(2,MNUM))*FORM2(MNUM,XMAA**2,XMRO2**2,XMRO1**2)
      F3 = CMPLX(COEF(3,MNUM))*FORM3(MNUM,XMAA**2,XMRO3**2,XMRO1**2)
      F4 = (-1.0*UROJ)*
     $CMPLX(COEF(4,MNUM))*FORM4(MNUM,XMAA**2,XMRO1**2,XMRO2**2,XMRO3**2)
      F5 = (-1.0)*UROJ/4.0/PI**2/FPI**2*
     $     CMPLX(COEF(5,MNUM))*FORM5(MNUM,XMAA**2,XMRO1**2,XMRO2**2)
      DO 45 I=1,4
      HADCUR(I)= CMPLX(FNORM(MNUM)) * (
     $  CMPLX(VEC1(I))*F1+CMPLX(VEC2(I))*F2+CMPLX(VEC3(I))*F3+
     $  CMPLX(VEC4(I))*F4+CMPLX(VEC5(I))*F5)
 45   CONTINUE
C
* CALCULATE PI-VECTORS: VECTOR AND AXIAL
      CALL CLVEC(HADCUR,PN,PIVEC)
      CALL CLAXI(HADCUR,PN,PIAKS)
      CALL CLNUT(HADCUR,BRAKM,HVM)
* SPIN INDEPENDENT PART OF DECAY DIFF-CROSS-SECT. IN TAU REST  FRAME
      BRAK= (GV**2+GA**2)*PT(4)*PIVEC(4) +2.*GV*GA*PT(4)*PIAKS(4)
     &     +2.*(GV**2-GA**2)*AMNUTA*AMTAU*BRAKM
      AMPLIT=(GFERMI)**2*BRAK/2.
      IF (MNUM.GE.9) THEN
        CALL TAUOLA_PRINT('MNUM='//CHAR(0))
!       PRINT *, 'MNUM=',MNUM
        ZNAK=-1.0
        XM1=0.0
        XM2=0.0
        XM3=0.0
        DO 77 K=1,4
        IF (K.EQ.4) ZNAK=1.0
        XM1=ZNAK*PIM1(K)**2+XM1
        XM2=ZNAK*PIM2(K)**2+XM2
        XM3=ZNAK*PIM3(K)**2+XM3
 77     CALL TAUOLA_PRINT('PIM1=PIM2=PIM3='//CHAR(0))
!77     PRINT *, 'PIM1=',PIM1(K),'PIM2=',PIM2(K),'PIM3=',PIM3(K)
        CALL TAUOLA_PRINT('XM1=XM2=XM3='//CHAR(0))
!       PRINT *, 'XM1=',SQRT(XM1),'XM2=',SQRT(XM2),'XM3=',SQRT(XM3)
        CALL TAUOLA_PRINT(
     $    '************************************************'//CHAR(0))
!       PRINT *, '************************************************'
      ENDIF
C POLARIMETER VECTOR IN TAU REST FRAME
      DO 90 I=1,3
      HV(I)=-(AMTAU*((GV**2+GA**2)*PIAKS(I)+2.*GV*GA*PIVEC(I)))
     &      +(GV**2-GA**2)*AMNUTA*AMTAU*HVM(I)
C HV IS DEFINED FOR TAU-    WITH GAMMA=B+HV*POL
      HV(I)=-HV(I)/BRAK
 90   CONTINUE
      END
      SUBROUTINE PROD5(P1,P2,P3,PIA)
C ----------------------------------------------------------------------
C external product of P1, P2, P3 4-momenta.
C SIGN is chosen +/- for decay of TAU +/- respectively
C     called by : DAMPAA, CLNUT
C ----------------------------------------------------------------------
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / IDFC  / IDFF
      REAL PIA(4),P1(4),P2(4),P3(4)
!
      BIND(C,NAME='tauola_idfc') /IDFC/
      BIND(C,NAME='tauola_jaki') /JAKI/
!
      DET2(I,J)=P1(I)*P2(J)-P2(I)*P1(J)
* -----------------------------------
      IF     (KTOM.EQ.1.OR.KTOM.EQ.-1) THEN
        SIGN= IDFF/ABS(IDFF)
      ELSEIF (KTOM.EQ.2) THEN
        SIGN=-IDFF/ABS(IDFF)
      ELSE
        CALL TAUOLA_PRINT('STOP IN PROD5: KTOM='//CHAR(0))
!       PRINT *, 'STOP IN PROD5: KTOM=',KTOM
        CALL TAUOLA_STOP()
      ENDIF
C
C EPSILON( p1(1), p2(2), p3(3), (4) ) = 1
C
      PIA(1)= -P3(3)*DET2(2,4)+P3(4)*DET2(2,3)+P3(2)*DET2(3,4)
      PIA(2)= -P3(4)*DET2(1,3)+P3(3)*DET2(1,4)-P3(1)*DET2(3,4)
      PIA(3)=  P3(4)*DET2(1,2)-P3(2)*DET2(1,4)+P3(1)*DET2(2,4)
      PIA(4)=  P3(3)*DET2(1,2)-P3(2)*DET2(1,3)+P3(1)*DET2(2,3)
C ALL FOUR INDICES ARE UP SO  PIA(3) AND PIA(4) HAVE SAME SIGN
      DO 20 I=1,4
  20  PIA(I)=PIA(I)*SIGN
      END
      SUBROUTINE DADNEW(MODE,ISGN,HV,PNU,PWB,PNPI,JNPI)
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUBMC / GAMPMC(30),GAMPER(30),NEVDEC(30)
      REAL*4            GAMPMC    ,GAMPER
      COMMON / INOUT / INUT,IOUT
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / TAUDCD /IDFFIN(9,NMODE),MULPIK(NMODE)
     &                ,NAMES
      CHARACTER NAMES(NMODE)*31
      REAL*4 PNU(4),PWB(4),PNPI(4,9),HV(4),HHV(4)
      REAL*4 PDUM1(4),PDUM2(4),PDUMI(4,9)
      REAL*4 RRR(3)
      REAL*4 WTMAX(NMODE)
      REAL*8              SWT(NMODE),SSWT(NMODE)
      COMMON /TAUOLA_WEIGHT_DADNEW/ WTMAX
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_inout') /INOUT/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taubmc') /TAUBMC/
      BIND(C,NAME='tauola_taudcd') /TAUDCD/
      BIND(C,NAME='tauola_weight_dadnew') /TAUOLA_WEIGHT_DADNEW/
!
      DIMENSION NEVRAW(NMODE),NEVOVR(NMODE),NEVACC(NMODE)
C
      DATA PI /3.141592653589793238462643/
      DATA IWARM/0/
C
      IF(MODE.EQ.-1) THEN
C     ===================
C -- AT THE MOMENT ONLY TWO DECAY MODES OF MULTIPIONS HAVE M. ELEM
        NMOD=NMODE
        IWARM=1
C       PRINT 7003
        DO 1 JNPI=1,NMOD
        NEVRAW(JNPI)=0
        NEVACC(JNPI)=0
        NEVOVR(JNPI)=0
        SWT(JNPI)=0
        SSWT(JNPI)=0
        WTMAX(JNPI)=-1.
C for 4pi phase space, need lots more trials at initialization,
C or use the WTMAX determined with many trials for default model:
        NTRIALS = 5000
        IF (JNPI.LE.NM4) THEN
C         11.Oct.11: fix for BINP and KARLSRUHE currents added
          WTMAX(JNPI) = PKORB(3,37+JNPI)
          NTRIALS = 20000
        END IF
        DO  I=1,NTRIALS
          IF    (JNPI.LE.0) THEN
            GOTO 903
          ELSEIF(JNPI.LE.NM4) THEN
            CALL DPH4PI(WT,HV,PDUM1,PDUM2,PDUMI,JNPI)
          ELSEIF(JNPI.LE.NM4+NM5) THEN
             CALL DPH5PI(WT,HV,PDUM1,PDUM2,PDUMI,JNPI)
          ELSEIF(JNPI.LE.NM4+NM5+NM6) THEN
            CALL DPHNPI(WT,HV,PDUM1,PDUM2,PDUMI,JNPI)
          ELSEIF(JNPI.LE.NM4+NM5+NM6+NM3) THEN
            INUM=JNPI-NM4-NM5-NM6
            CALL DPHSPK(WT,HV,PDUM1,PDUM2,PDUMI,INUM)
          ELSEIF(JNPI.LE.NM4+NM5+NM6+NM3+NM2) THEN
            INUM=JNPI-NM4-NM5-NM6-NM3
            CALL DPHSRK(WT,HV,PDUM1,PDUM2,PDUMI,INUM)
          ELSE
           GOTO 903
          ENDIF
        IF(WT.GT.WTMAX(JNPI)/1.2) WTMAX(JNPI)=WT*1.2
        ENDDO
C       PRINT *,' DADNEW JNPI,NTRIALS,WTMAX =',JNPI,NTRIALS,WTMAX(JNPI)
C       CALL HBOOK1(801,'WEIGHT DISTRIBUTION  DADNPI    $',100,0.,2.,.0)
C       PRINT 7004,WTMAX(JNPI)
1       CONTINUE
        CALL TAUOLA_PRINT('tauola.f:4169: (suppressed)'//CHAR(0))
!       WRITE(IOUT,7005)
C
      ELSEIF(MODE.EQ. 0) THEN
C     =======================
        IF(IWARM.EQ.0) GOTO 902
C
300     CONTINUE
          IF    (JNPI.LE.0) THEN
            GOTO 903
          ELSEIF(JNPI.LE.NM4) THEN
             CALL DPH4PI(WT,HHV,PNU,PWB,PNPI,JNPI)
          ELSEIF(JNPI.LE.NM4+NM5) THEN
             CALL DPH5PI(WT,HHV,PNU,PWB,PNPI,JNPI)
          ELSEIF(JNPI.LE.NM4+NM5+NM6) THEN
            CALL DPHNPI(WT,HHV,PNU,PWB,PNPI,JNPI)
          ELSEIF(JNPI.LE.NM4+NM5+NM6+NM3) THEN
            INUM=JNPI-NM4-NM5-NM6
            CALL DPHSPK(WT,HHV,PNU,PWB,PNPI,INUM)
          ELSEIF(JNPI.LE.NM4+NM5+NM6+NM3+NM2) THEN
            INUM=JNPI-NM4-NM5-NM6-NM3
            CALL DPHSRK(WT,HHV,PNU,PWB,PNPI,INUM)
          ELSE
           GOTO 903
          ENDIF
            DO I=1,4
              HV(I)=-ISGN*HHV(I)
            ENDDO
C       CALL HFILL(801,WT/WTMAX(JNPI))
        NEVRAW(JNPI)=NEVRAW(JNPI)+1
        SWT(JNPI)=SWT(JNPI)+WT
cccM.S.>>>>>>
cc        SSWT(JNPI)=SSWT(JNPI)+WT**2
        SSWT(JNPI)=SSWT(JNPI)+dble(WT)**2
cccM.S.<<<<<<
        CALL RANMAR(RRR,3)
        RN=RRR(1)
        IF(WT.GT.WTMAX(JNPI)) NEVOVR(JNPI)=NEVOVR(JNPI)+1
        IF(RN*WTMAX(JNPI).GT.WT) GOTO 300
C ROTATIONS TO BASIC TAU REST FRAME
        COSTHE=-1.+2.*RRR(2)
        THET=ACOS(COSTHE)
        PHI =2*PI*RRR(3)
        CALL ROTOR2(THET,PNU,PNU)
        CALL ROTOR3( PHI,PNU,PNU)
        CALL ROTOR2(THET,PWB,PWB)
        CALL ROTOR3( PHI,PWB,PWB)
        CALL ROTOR2(THET,HV,HV)
        CALL ROTOR3( PHI,HV,HV)
        ND=MULPIK(JNPI)
        DO 301 I=1,ND
        CALL ROTOR2(THET,PNPI(1,I),PNPI(1,I))
        CALL ROTOR3( PHI,PNPI(1,I),PNPI(1,I))
301     CONTINUE
        NEVACC(JNPI)=NEVACC(JNPI)+1
C
      ELSEIF(MODE.EQ. 1) THEN
C     =======================
        DO 500 JNPI=1,NMOD
          IF(NEVRAW(JNPI).EQ.0) GOTO 500
          PARGAM=SWT(JNPI)/FLOAT(NEVRAW(JNPI)+1)
          ERROR=0
          IF(NEVRAW(JNPI).NE.0)
     &    ERROR=SQRT(SSWT(JNPI)/SWT(JNPI)**2-1./FLOAT(NEVRAW(JNPI)))
          RAT=PARGAM/GAMEL
          CALL TAUOLA_PRINT('tauola.f:4234: (suppressed)'//CHAR(0))
!         WRITE(IOUT, 7010) NAMES(JNPI),
!    &     NEVRAW(JNPI),NEVACC(JNPI),NEVOVR(JNPI),PARGAM,RAT,ERROR
CC        CALL HPRINT(801)
          GAMPMC(8+JNPI-1)=RAT
          GAMPER(8+JNPI-1)=ERROR
CAM       NEVDEC(8+JNPI-1)=NEVACC(JNPI)
  500     CONTINUE
      ENDIF
C     =====
      RETURN
!7003 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADNEW INITIALISATION ********',9X,1H*
!    $ )
!7004 FORMAT(' *',E20.5,5X,'WTMAX  = MAXIMUM WEIGHT  ',9X,1H*/)
!7005 FORMAT(
!    $  /,1X,15(5H*****)/)
!7010 FORMAT(///1X,15(5H*****)
!    $ /,' *',     25X,'******** DADNEW FINAL REPORT  ******** ',9X,1H*
!    $ /,' *',     25X,'CHANNEL:',A31                           ,9X,1H*
!    $ /,' *',I20  ,5X,'NEVRAW = NO. OF DECAYS TOTAL           ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVACC = NO. OF DECAYS ACCEPTED        ',9X,1H*
!    $ /,' *',I20  ,5X,'NEVOVR = NO. OF OVERWEIGHTED EVENTS    ',9X,1H*
!    $ /,' *',E20.5,5X,'PARTIAL WTDTH IN GEV UNITS             ',9X,1H*
!    $ /,' *',F20.9,5X,'IN UNITS GFERMI**2*MASS**5/192/PI**3   ',9X,1H*
!    $ /,' *',F20.8,5X,'RELATIVE ERROR OF PARTIAL WIDTH        ',9X,1H*
!    $  /,1X,15(5H*****)/)
 902  CALL TAUOLA_PRINT('----- DADNEW: LACK OF INITIALISATION'//CHAR(0))
!902  WRITE(IOUT, 9020)
!9020 FORMAT(' ----- DADNEW: LACK OF INITIALISATION')
      CALL TAUOLA_STOP()
 903  CALL TAUOLA_PRINT('----- DADNEW: WRONG JNPI'//CHAR(0))
!903  WRITE(IOUT, 9030) JNPI,MODE
!9030 FORMAT(' ----- DADNEW: WRONG JNPI',2I5)
      CALL TAUOLA_STOP()
      END
      SUBROUTINE DPH4PI(DGAMT,HV,PN,PAA,PMULT,JNPI)
C ----------------------------------------------------------------------
* IT SIMULATES A1  DECAY IN TAU REST FRAME WITH
* Z-AXIS ALONG A1  MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  HV(4),PT(4),PN(4),PAA(4),PIM1(4),PIM2(4),PIPL(4),PMULT(4,9)
      REAL  PR(4),PIZ(4)
      REAL*4 RRR(9)
      REAL*8 UU,FF,FF1,FF2,FF3,FF4,GG1,GG2,GG3,GG4,RR
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      XLAM(X,Y,Z)=SQRT(ABS((X-Y-Z)**2-4.0*Y*Z))
C AMRO, GAMRO IS ONLY A PARAMETER FOR GETING HIGHT EFFICIENCY
C
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1./2**23/PI**11
      PHSP=1./2**5/PI**2
      IF (JNPI.EQ.1) THEN
       PREZ=0.7
       AMP1=AMPI
       AMP2=AMPI
       AMP3=AMPI
       AMP4=AMPIZ
       AMRX=PKORB(1,14)
       GAMRX=PKORB(2,14)
C AJW: cant simply change AMROP, etc, here!
C CHOICE is a by-hand tuning/optimization, no simple relationship
C to actual resonance masses (accd to Z.Was).
C What matters in the end is what you put in formf/curr .
        AMROP =1.2
        GAMROP=.46
      ELSE
       PREZ=0.0
       AMP1=AMPIZ
       AMP2=AMPIZ
       AMP3=AMPIZ
       AMP4=AMPI
       AMRX=1.4
       GAMRX=.6
        AMROP =AMRX
        GAMROP=GAMRX
      ENDIF
      RRB=0.3
      CALL CHOICE(100+JNPI,RRB,ICHAN,PROB1,PROB2,PROB3,
     $            AMROP,GAMROP,AMRX,GAMRX,AMRB,GAMRB)
      PREZ=PROB1+PROB2
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C
      CALL RANMAR(RRR,9)
C
* MASSES OF 4, 3 AND 2 PI SYSTEMS
C 3 PI WITH SAMPLING FOR RESONANCE
CAM
        RR1=RRR(6)
        AMS1=(AMP1+AMP2+AMP3+AMP4)**2
        AMS2=(AMTAU-AMNUTA)**2
        ALP1=ATAN((AMS1-AMROP**2)/AMROP/GAMROP)
        ALP2=ATAN((AMS2-AMROP**2)/AMROP/GAMROP)
        ALP=ALP1+RR1*(ALP2-ALP1)
        AM4SQ =AMROP**2+AMROP*GAMROP*TAN(ALP)
        AM4 =SQRT(AM4SQ)
        PHSPAC=PHSPAC*
     $         ((AM4SQ-AMROP**2)**2+(AMROP*GAMROP)**2)/(AMROP*GAMROP)
        PHSPAC=PHSPAC*(ALP2-ALP1)
C
        RR1=RRR(1)
        AMS1=(AMP2+AMP3+AMP4)**2
        AMS2=(AM4-AMP1)**2
        IF (RRR(9).GT.PREZ) THEN
          AM3SQ=AMS1+   RR1*(AMS2-AMS1)
          AM3 =SQRT(AM3SQ)
C --- this part of jacobian will be recovered later
          FF1=AMS2-AMS1
        ELSE
* PHASE SPACE WITH SAMPLING FOR OMEGA RESONANCE,
        ALP1=ATAN((AMS1-AMRX**2)/AMRX/GAMRX)
        ALP2=ATAN((AMS2-AMRX**2)/AMRX/GAMRX)
        ALP=ALP1+RR1*(ALP2-ALP1)
        AM3SQ =AMRX**2+AMRX*GAMRX*TAN(ALP)
        AM3 =SQRT(AM3SQ)
C --- THIS PART OF THE JACOBIAN WILL BE RECOVERED LATER ---------------
        FF1=((AM3SQ-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)
        FF1=FF1*(ALP2-ALP1)
        ENDIF
C MASS OF 2
        RR2=RRR(2)
        AMS1=(AMP3+AMP4)**2
        AMS2=(AM3-AMP2)**2
* FLAT PHASE SPACE;
        AM2SQ=AMS1+   RR2*(AMS2-AMS1)
        AM2 =SQRT(AM2SQ)
C --- this part of jacobian will be recovered later
        FF2=(AMS2-AMS1)
*  2 RESTFRAME, DEFINE PIZ AND PIPL
        ENQ1=(AM2SQ-AMP3**2+AMP4**2)/(2*AM2)
        ENQ2=(AM2SQ+AMP3**2-AMP4**2)/(2*AM2)
        PPI=         ENQ1**2-AMP4**2
        PPPI=SQRT(ABS(ENQ1**2-AMP4**2))
        PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AM2)
* PIZ   MOMENTUM IN 2 REST FRAME
        CALL SPHERA(PPPI,PIZ)
        PIZ(4)=ENQ1
* PIPL  MOMENTUM IN 2 REST FRAME
        DO 30 I=1,3
 30     PIPL(I)=-PIZ(I)
        PIPL(4)=ENQ2
* 3 REST FRAME, DEFINE PIM1
*       PR   MOMENTUM
        PR(1)=0
        PR(2)=0
        PR(4)=1./(2*AM3)*(AM3**2+AM2**2-AMP2**2)
        PR(3)= SQRT(ABS(PR(4)**2-AM2**2))
        PPI  =          PR(4)**2-AM2**2
*       PIM1  MOMENTUM
        PIM1(1)=0
        PIM1(2)=0
        PIM1(4)=1./(2*AM3)*(AM3**2-AM2**2+AMP2**2)
        PIM1(3)=-PR(3)
C --- this part of jacobian will be recovered later
        FF3=(4*PI)*(2*PR(3)/AM3)
* OLD PIONS BOOSTED FROM 2 REST FRAME TO 3 REST FRAME
      EXE=(PR(4)+PR(3))/AM2
      CALL BOSTR3(EXE,PIZ,PIZ)
      CALL BOSTR3(EXE,PIPL,PIPL)
      RR3=RRR(3)
      RR4=RRR(4)
      THET =ACOS(-1.+2*RR3)
      PHI = 2*PI*RR4
      CALL ROTPOL(THET,PHI,PIPL)
      CALL ROTPOL(THET,PHI,PIM1)
      CALL ROTPOL(THET,PHI,PIZ)
      CALL ROTPOL(THET,PHI,PR)
* 4  REST FRAME, DEFINE PIM2
*       PR   MOMENTUM
        PR(1)=0
        PR(2)=0
        PR(4)=1./(2*AM4)*(AM4**2+AM3**2-AMP1**2)
        PR(3)= SQRT(ABS(PR(4)**2-AM3**2))
        PPI  =          PR(4)**2-AM3**2
*       PIM2 MOMENTUM
        PIM2(1)=0
        PIM2(2)=0
        PIM2(4)=1./(2*AM4)*(AM4**2-AM3**2+AMP1**2)
        PIM2(3)=-PR(3)
C --- this part of jacobian will be recovered later
        FF4=(4*PI)*(2*PR(3)/AM4)
* OLD PIONS BOOSTED FROM 3 REST FRAME TO 4 REST FRAME
      EXE=(PR(4)+PR(3))/AM3
      CALL BOSTR3(EXE,PIZ,PIZ)
      CALL BOSTR3(EXE,PIPL,PIPL)
      CALL BOSTR3(EXE,PIM1,PIM1)
      RR3=RRR(7)
      RR4=RRR(8)
      THET =ACOS(-1.+2*RR3)
      PHI = 2*PI*RR4
      CALL ROTPOL(THET,PHI,PIPL)
      CALL ROTPOL(THET,PHI,PIM1)
      CALL ROTPOL(THET,PHI,PIM2)
      CALL ROTPOL(THET,PHI,PIZ)
      CALL ROTPOL(THET,PHI,PR)
C
* NOW TO THE TAU REST FRAME, DEFINE PAA AND NEUTRINO MOMENTA
* PAA  MOMENTUM
      PAA(1)=0
      PAA(2)=0
      PAA(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AM4**2)
      PAA(3)= SQRT(ABS(PAA(4)**2-AM4**2))
      PPI   =          PAA(4)**2-AM4**2
      PHSPAC=PHSPAC*(4*PI)*(2*PAA(3)/AMTAU)
      PHSP=PHSP*(4*PI)*(2*PAA(3)/AMTAU)
* TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AM4**2)
      PN(3)=-PAA(3)
C ZBW 20.12.2002 bug fix
        IF(RRR(9).LE.0.5*PREZ) THEN
         DO 72 I=1,4
         X=PIM1(I)
         PIM1(I)=PIM2(I)
 72      PIM2(I)=X
        ENDIF
C end of bug fix
C WE INCLUDE REMAINING PART OF THE JACOBIAN
C --- FLAT CHANNEL
        AM3SQ=(PIM1(4)+PIZ(4)+PIPL(4))**2-(PIM1(3)+PIZ(3)+PIPL(3))**2
     $       -(PIM1(2)+PIZ(2)+PIPL(2))**2-(PIM1(1)+PIZ(1)+PIPL(1))**2
        AMS2=(AM4-AMP2)**2
        AMS1=(AMP1+AMP3+AMP4)**2
        FF1=(AMS2-AMS1)
        AMS1=(AMP3+AMP4)**2
        AMS2=(SQRT(AM3SQ)-AMP1)**2
        FF2=AMS2-AMS1
        FF3=(4*PI)*(XLAM(AM2**2,AMP1**2,AM3SQ)/AM3SQ)
        FF4=(4*PI)*(XLAM(AM3SQ,AMP2**2,AM4**2)/AM4**2)
        UU=FF1*FF2*FF3*FF4
C --- FIRST CHANNEL
        AM3SQ=(PIM1(4)+PIZ(4)+PIPL(4))**2-(PIM1(3)+PIZ(3)+PIPL(3))**2
     $       -(PIM1(2)+PIZ(2)+PIPL(2))**2-(PIM1(1)+PIZ(1)+PIPL(1))**2
        AMS2=(AM4-AMP2)**2
        AMS1=(AMP1+AMP3+AMP4)**2
        ALP1=ATAN((AMS1-AMRX**2)/AMRX/GAMRX)
        ALP2=ATAN((AMS2-AMRX**2)/AMRX/GAMRX)
        FF1=((AM3SQ-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)
        FF1=FF1*(ALP2-ALP1)
        AMS1=(AMP3+AMP4)**2
        AMS2=(SQRT(AM3SQ)-AMP1)**2
        FF2=AMS2-AMS1
        FF3=(4*PI)*(XLAM(AM2**2,AMP1**2,AM3SQ)/AM3SQ)
        FF4=(4*PI)*(XLAM(AM3SQ,AMP2**2,AM4**2)/AM4**2)
        FF=FF1*FF2*FF3*FF4
C --- SECOND CHANNEL
        AM3SQ=(PIM2(4)+PIZ(4)+PIPL(4))**2-(PIM2(3)+PIZ(3)+PIPL(3))**2
     $       -(PIM2(2)+PIZ(2)+PIPL(2))**2-(PIM2(1)+PIZ(1)+PIPL(1))**2
        AMS2=(AM4-AMP1)**2
        AMS1=(AMP2+AMP3+AMP4)**2
        ALP1=ATAN((AMS1-AMRX**2)/AMRX/GAMRX)
        ALP2=ATAN((AMS2-AMRX**2)/AMRX/GAMRX)
        GG1=((AM3SQ-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)
        GG1=GG1*(ALP2-ALP1)
        AMS1=(AMP3+AMP4)**2
        AMS2=(SQRT(AM3SQ)-AMP2)**2
        GG2=AMS2-AMS1
        GG3=(4*PI)*(XLAM(AM2**2,AMP2**2,AM3SQ)/AM3SQ)
        GG4=(4*PI)*(XLAM(AM3SQ,AMP1**2,AM4**2)/AM4**2)
        GG=GG1*GG2*GG3*GG4
C --- JACOBIAN AVERAGED OVER THE TWO
        ! 05.10.2011 missing factor in IF( (0.5*PREZ* ... (1-PREZ)* ...) added
        IF ( (0.5*PREZ*(FF+GG)*UU+(1.0-PREZ)*FF*GG).GT.0.0D0) THEN
          RR=FF*GG*UU/(0.5*PREZ*(FF+GG)*UU+(1.0-PREZ)*FF*GG)
          PHSPAC=PHSPAC*RR
        ELSE
          PHSPAC=0.0
        ENDIF
* MOMENTA OF THE TWO PI-MINUS ARE RANDOMLY SYMMETRISED
       IF (JNPI.EQ.1) THEN
        RR5= RRR(5)
        IF(RR5.LE.0.5) THEN
         DO 70 I=1,4
         X=PIM1(I)
         PIM1(I)=PIM2(I)
 70      PIM2(I)=X
        ENDIF
        PHSPAC=PHSPAC/2.
       ELSE
C MOMENTA OF PI0-S ARE GENERATED UNIFORMLY ONLY IF PREZ=0.0
        RR5= RRR(5)
        IF(RR5.LE.0.5) THEN
         DO 71 I=1,4
         X=PIM1(I)
         PIM1(I)=PIM2(I)
 71      PIM2(I)=X
        ENDIF
        PHSPAC=PHSPAC/6.
       ENDIF
* ALL PIONS BOOSTED FROM  4  REST FRAME TO TAU REST FRAME
* Z-AXIS ANTIPARALLEL TO NEUTRINO MOMENTUM
      EXE=(PAA(4)+PAA(3))/AM4
      CALL BOSTR3(EXE,PIZ,PIZ)
      CALL BOSTR3(EXE,PIPL,PIPL)
      CALL BOSTR3(EXE,PIM1,PIM1)
      CALL BOSTR3(EXE,PIM2,PIM2)
      CALL BOSTR3(EXE,PR,PR)
C PARTIAL WIDTH CONSISTS OF PHASE SPACE AND AMPLITUDE
C CHECK ON CONSISTENCY WITH DADNPI, THEN, CODE BREAKES UNIFORM PION
C DISTRIBUTION IN HADRONIC SYSTEM
CAM     Assume neutrino mass=0. and sum over final polarisation
C      AMX2=AM4**2
C      BRAK= 2*(AMTAU**2-AMX2) * (AMTAU**2+2.*AMX2)
C      AMPLIT=CCABIB**2*GFERMI**2/2. * BRAK * AMX2*SIGEE(AMX2,1)
      IF     (JNPI.EQ.1) THEN
        CALL DAM4PI(JNPI,PT,PN,PIM1,PIM2,PIZ,PIPL,AMPLIT,HV)
      ELSEIF (JNPI.EQ.2) THEN
        CALL DAM4PI(JNPI,PT,PN,PIM1,PIM2,PIPL,PIZ,AMPLIT,HV)
      ENDIF
      DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
C PHASE SPACE CHECK
C      DGAMT=PHSPAC
      DO 77 K=1,4
        PMULT(K,1)=PIM1(K)
        PMULT(K,2)=PIM2(K)
        PMULT(K,3)=PIPL(K)
        PMULT(K,4)=PIZ (K)
 77   CONTINUE
      END
      SUBROUTINE DAM4PI(MNUM,PT,PN,PIM1,PIM2,PIM3,PIM4,AMPLIT,HV)
C ----------------------------------------------------------------------
* CALCULATES DIFFERENTIAL CROSS SECTION AND POLARIMETER VECTOR
* FOR TAU DECAY INTO 4 PI MODES
* ALL SPIN EFFECTS IN THE FULL DECAY CHAIN ARE TAKEN INTO ACCOUNT.
* CALCULATIONS DONE IN TAU REST FRAME WITH Z-AXIS ALONG NEUTRINO MOMENT
C MNUM DECAY MODE IDENTIFIER.
C
C     called by : DPHSAA
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL  HV(4),PT(4),PN(4),PIM1(4),PIM2(4),PIM3(4),PIM4(4)
      REAL  PIVEC(4),PIAKS(4),HVM(4)
      COMPLEX HADCUR(4),FORM1,FORM2,FORM3,FORM4,FORM5
      EXTERNAL FORM1,FORM2,FORM3,FORM4,FORM5
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
C
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      CALL CURR_CLEO(MNUM,PIM1,PIM2,PIM3,PIM4,HADCUR)
C
* CALCULATE PI-VECTORS: VECTOR AND AXIAL
      CALL CLVEC(HADCUR,PN,PIVEC)
      CALL CLAXI(HADCUR,PN,PIAKS)
      CALL CLNUT(HADCUR,BRAKM,HVM)
* SPIN INDEPENDENT PART OF DECAY DIFF-CROSS-SECT. IN TAU REST  FRAME
      BRAK= (GV**2+GA**2)*PT(4)*PIVEC(4) +2.*GV*GA*PT(4)*PIAKS(4)
     &     +2.*(GV**2-GA**2)*AMNUTA*AMTAU*BRAKM
      AMPLIT=(CCABIB*GFERMI)**2*BRAK/2.
C POLARIMETER VECTOR IN TAU REST FRAME
      DO 90 I=1,3
      HV(I)=-(AMTAU*((GV**2+GA**2)*PIAKS(I)+2.*GV*GA*PIVEC(I)))
     &      +(GV**2-GA**2)*AMNUTA*AMTAU*HVM(I)
C HV IS DEFINED FOR TAU-    WITH GAMMA=B+HV*POL
      HV(I)=-HV(I)/BRAK
 90   CONTINUE
      END
      SUBROUTINE DPH5PI(DGAMT,HV,PN,PAA,PMULT,JNPI)
C ----------------------------------------------------------------------
* IT SIMULATES 5pi DECAY IN TAU REST FRAME WITH
* Z-AXIS ALONG 5pi MOMENTUM
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / TAUDCD /IDFFIN(9,NMODE),MULPIK(NMODE)
     &                ,NAMES
      CHARACTER NAMES(NMODE)*31
      REAL  HV(4),PT(4),PN(4),PAA(4),PMULT(4,9)
      REAL*4 PR(4),PI1(4),PI2(4),PI3(4),PI4(4),PI5(4)
      REAL*8 AMP1,AMP2,AMP3,AMP4,AMP5,ams1,ams2,amom,gamom
      REAL*8 AM5SQ,AM4SQ,AM3SQ,AM2SQ,AM5,AM4,AM3
      REAL*4 RRR(10)
      REAL*8 gg1,gg2,gg3,ff1,ff2,ff3,ff4,alp,alp1,alp2
      REAL*8 XM,AM,GAMMA
ccM.S.>>>>>>
      real*8 phspac
ccM.S.<<<<<<
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
      data fpi /93.3e-3/
c
      COMPLEX BWIGN
C
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taudcd') /TAUDCD/
!
      BWIGN(XM,AM,GAMMA)=XM**2/CMPLX(XM**2-AM**2,GAMMA*AM)
C
      AMOM=.782
      GAMOM=0.0085
c
C 6 BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
C D**3 P /2E/(2PI)**3 (2PI)**4 DELTA4(SUM P)
      PHSPAC=1./2**29/PI**14
c     PHSPAC=1./2**5/PI**2
C init 5pi decay mode (JNPI)
      AMP1=DCDMAS(IDFFIN(1,JNPI))
      AMP2=DCDMAS(IDFFIN(2,JNPI))
      AMP3=DCDMAS(IDFFIN(3,JNPI))
      AMP4=DCDMAS(IDFFIN(4,JNPI))
      AMP5=DCDMAS(IDFFIN(5,JNPI))
c
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C
      CALL RANMAR(RRR,10)
C
c masses of 5, 4, 3 and 2 pi systems
c 3 pi with sampling for omega resonance
cam
c mass of 5   (12345)
      rr1=rrr(10)
      ams1=(amp1+amp2+amp3+amp4+amp5)**2
      ams2=(amtau-amnuta)**2
      am5sq=ams1+   rr1*(ams2-ams1)
      am5 =sqrt(am5sq)
      phspac=phspac*(ams2-ams1)
c
c mass of 4   (2345)
c flat phase space
      rr1=rrr(9)
      ams1=(amp2+amp3+amp4+amp5)**2
      ams2=(am5-amp1)**2
      am4sq=ams1+   rr1*(ams2-ams1)
      am4 =sqrt(am4sq)
      gg1=ams2-ams1
c
c mass of 3   (234)
C phase space with sampling for omega resonance
      rr1=rrr(1)
      ams1=(amp2+amp3+amp4)**2
      ams2=(am4-amp5)**2
      alp1=atan((ams1-amom**2)/amom/gamom)
      alp2=atan((ams2-amom**2)/amom/gamom)
      alp=alp1+rr1*(alp2-alp1)
      am3sq =amom**2+amom*gamom*tan(alp)
      am3 =sqrt(am3sq)
c --- this part of the jacobian will be recovered later ---------------
      gg2=((am3sq-amom**2)**2+(amom*gamom)**2)/(amom*gamom)
      gg2=gg2*(alp2-alp1)
c flat phase space;
C      am3sq=ams1+   rr1*(ams2-ams1)
C      am3 =sqrt(am3sq)
c --- this part of jacobian will be recovered later
C      gg2=ams2-ams1
c
C mass of 2  (34)
      rr2=rrr(2)
      ams1=(amp3+amp4)**2
      ams2=(am3-amp2)**2
c flat phase space;
      am2sq=ams1+   rr2*(ams2-ams1)
      am2 =sqrt(am2sq)
c --- this part of jacobian will be recovered later
      gg3=ams2-ams1
c
c (34) restframe, define pi3 and pi4
      enq1=(am2sq+amp3**2-amp4**2)/(2*am2)
      enq2=(am2sq-amp3**2+amp4**2)/(2*am2)
      ppi=          enq1**2-amp3**2
      pppi=sqrt(abs(enq1**2-amp3**2))
      ff1=(4*pi)*(2*pppi/am2)
c pi3   momentum in (34) rest frame
      call sphera(pppi,pi3)
      pi3(4)=enq1
c pi4   momentum in (34) rest frame
      do 30 i=1,3
 30   pi4(i)=-pi3(i)
      pi4(4)=enq2
c
c (234) rest frame, define pi2
c pr   momentum
      pr(1)=0
      pr(2)=0
      pr(4)=1./(2*am3)*(am3**2+am2**2-amp2**2)
      pr(3)= sqrt(abs(pr(4)**2-am2**2))
      ppi  =          pr(4)**2-am2**2
c pi2   momentum
      pi2(1)=0
      pi2(2)=0
      pi2(4)=1./(2*am3)*(am3**2-am2**2+amp2**2)
      pi2(3)=-pr(3)
c --- this part of jacobian will be recovered later
      ff2=(4*pi)*(2*pr(3)/am3)
c old pions boosted from 2 rest frame to 3 rest frame
      exe=(pr(4)+pr(3))/am2
      call bostr3(exe,pi3,pi3)
      call bostr3(exe,pi4,pi4)
      rr3=rrr(3)
      rr4=rrr(4)
      thet =acos(-1.+2*rr3)
      phi = 2*pi*rr4
      call rotpol(thet,phi,pi2)
      call rotpol(thet,phi,pi3)
      call rotpol(thet,phi,pi4)
C
C (2345)  rest frame, define pi5
c pr   momentum
      pr(1)=0
      pr(2)=0
      pr(4)=1./(2*am4)*(am4**2+am3**2-amp5**2)
      pr(3)= sqrt(abs(pr(4)**2-am3**2))
      ppi  =          pr(4)**2-am3**2
c pi5  momentum
      pi5(1)=0
      pi5(2)=0
      pi5(4)=1./(2*am4)*(am4**2-am3**2+amp5**2)
      pi5(3)=-pr(3)
c --- this part of jacobian will be recovered later
      ff3=(4*pi)*(2*pr(3)/am4)
c old pions boosted from 3 rest frame to 4 rest frame
      exe=(pr(4)+pr(3))/am3
      call bostr3(exe,pi2,pi2)
      call bostr3(exe,pi3,pi3)
      call bostr3(exe,pi4,pi4)
      rr3=rrr(5)
      rr4=rrr(6)
      thet =acos(-1.+2*rr3)
      phi = 2*pi*rr4
      call rotpol(thet,phi,pi2)
      call rotpol(thet,phi,pi3)
      call rotpol(thet,phi,pi4)
      call rotpol(thet,phi,pi5)
C
C (12345)  rest frame, define pi1
c pr   momentum
      pr(1)=0
      pr(2)=0
      pr(4)=1./(2*am5)*(am5**2+am4**2-amp1**2)
      pr(3)= sqrt(abs(pr(4)**2-am4**2))
      ppi  =          pr(4)**2-am4**2
c pi1  momentum
      pi1(1)=0
      pi1(2)=0
      pi1(4)=1./(2*am5)*(am5**2-am4**2+amp1**2)
      pi1(3)=-pr(3)
c --- this part of jacobian will be recovered later
      ff4=(4*pi)*(2*pr(3)/am5)
c old pions boosted from 4 rest frame to 5 rest frame
      exe=(pr(4)+pr(3))/am4
      call bostr3(exe,pi2,pi2)
      call bostr3(exe,pi3,pi3)
      call bostr3(exe,pi4,pi4)
      call bostr3(exe,pi5,pi5)
      rr3=rrr(7)
      rr4=rrr(8)
      thet =acos(-1.+2*rr3)
      phi = 2*pi*rr4
      call rotpol(thet,phi,pi1)
      call rotpol(thet,phi,pi2)
      call rotpol(thet,phi,pi3)
      call rotpol(thet,phi,pi4)
      call rotpol(thet,phi,pi5)
c
* now to the tau rest frame, define paa and neutrino momenta
* paa  momentum
      paa(1)=0
      paa(2)=0
c     paa(4)=1./(2*amtau)*(amtau**2-amnuta**2+am5**2)
c     paa(3)= sqrt(abs(paa(4)**2-am5**2))
c     ppi   =          paa(4)**2-am5**2
      paa(4)=1./(2*amtau)*(amtau**2-amnuta**2+am5sq)
      paa(3)= sqrt(abs(paa(4)**2-am5sq))
      ppi   =          paa(4)**2-am5sq
      phspac=phspac*(4*pi)*(2*paa(3)/amtau)
* tau-neutrino momentum
      pn(1)=0
      pn(2)=0
      pn(4)=1./(2*amtau)*(amtau**2+amnuta**2-am5**2)
      pn(3)=-paa(3)
c
      phspac=phspac * gg1*gg2*gg3*ff1*ff2*ff3*ff4
c
C all pions boosted from  5  rest frame to tau rest frame
C z-axis antiparallel to neutrino momentum
      exe=(paa(4)+paa(3))/am5
      call bostr3(exe,pi1,pi1)
      call bostr3(exe,pi2,pi2)
      call bostr3(exe,pi3,pi3)
      call bostr3(exe,pi4,pi4)
      call bostr3(exe,pi5,pi5)
c
C partial width consists of phase space and amplitude
C AMPLITUDE  (cf YS.Tsai Phys.Rev.D4,2821(1971)
C    or F.Gilman SH.Rhie Phys.Rev.D31,1066(1985)
C
      PXQ=AMTAU*PAA(4)
      PXN=AMTAU*PN(4)
      QXN=PAA(4)*PN(4)-PAA(1)*PN(1)-PAA(2)*PN(2)-PAA(3)*PN(3)
      BRAK=2*(GV**2+GA**2)*(2*PXQ*QXN+AM5SQ*PXN)
     &    -6*(GV**2-GA**2)*AMTAU*AMNUTA*AM5SQ
      fompp = cabs(bwign(am3,amom,gamom))**2
c normalisation factor (to some numerical undimensioned factor;
c cf R.Fischer et al ZPhys C3, 313 (1980))
      fnorm = 1/fpi**6
c     AMPLIT=CCABIB**2*GFERMI**2/2. * BRAK * AM5SQ*SIGEE(AM5SQ,JNPI)
      AMPLIT=CCABIB**2*GFERMI**2/2. * BRAK
      amplit = amplit * fompp * fnorm
c phase space test
c     amplit = amplit * fnorm
      DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
c ignore spin terms
      DO 40 I=1,3
 40   HV(I)=0.
c
      do 77 k=1,4
        pmult(k,1)=pi1(k)
        pmult(k,2)=pi2(k)
        pmult(k,3)=pi3(k)
        pmult(k,4)=pi4(k)
        pmult(k,5)=pi5(k)
 77   continue
      return
C missing: transposition of identical particles, startistical factors
C for identical matrices, polarimetric vector. Matrix element rather naive.
C flat phase space in pion system + with breit wigner for omega
C anyway it is better than nothing, and code is improvable.
      end
      SUBROUTINE DPHSRK(DGAMT,HV,PN,PR,PMULT,INUM)
C ----------------------------------------------------------------------
C IT SIMULATES RHO DECAY IN TAU REST FRAME WITH
C Z-AXIS ALONG RHO MOMENTUM
C Rho decays to K Kbar
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      REAL  HV(4),PT(4),PN(4),PR(4),PKC(4),PKZ(4),QQ(4),PMULT(4,9)
      REAL RR1(1)
      DATA PI /3.141592653589793238462643/
      DATA ICONT /0/
C
C THREE BODY PHASE SPACE NORMALISED AS IN BJORKEN-DRELL
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      PHSPAC=1./2**11/PI**5
C TAU MOMENTUM
      PT(1)=0.
      PT(2)=0.
      PT(3)=0.
      PT(4)=AMTAU
C MASS OF (REAL/VIRTUAL) RHO
      AMS1=(AMK+AMKZ)**2
      AMS2=(AMTAU-AMNUTA)**2
C FLAT PHASE SPACE
      CALL RANMAR(RR1,1)
      AMX2=AMS1+   RR1(1)*(AMS2-AMS1)
      AMX=SQRT(AMX2)
      PHSPAC=PHSPAC*(AMS2-AMS1)
C PHASE SPACE WITH SAMPLING FOR RHO RESONANCE
c     ALP1=ATAN((AMS1-AMRO**2)/AMRO/GAMRO)
c     ALP2=ATAN((AMS2-AMRO**2)/AMRO/GAMRO)
CAM
 100  CONTINUE
c     CALL RANMAR(RR1,1)
c     ALP=ALP1+RR1(1)*(ALP2-ALP1)
c     AMX2=AMRO**2+AMRO*GAMRO*TAN(ALP)
c     AMX=SQRT(AMX2)
c     IF(AMX.LT.(AMK+AMKZ)) GO TO 100
CAM
c     PHSPAC=PHSPAC*((AMX2-AMRO**2)**2+(AMRO*GAMRO)**2)/(AMRO*GAMRO)
c     PHSPAC=PHSPAC*(ALP2-ALP1)
C
C TAU-NEUTRINO MOMENTUM
      PN(1)=0
      PN(2)=0
      PN(4)=1./(2*AMTAU)*(AMTAU**2+AMNUTA**2-AMX**2)
      PN(3)=-SQRT(ABS((PN(4)-AMNUTA)*(PN(4)+AMNUTA)))
C RHO MOMENTUM
      PR(1)=0
      PR(2)=0
      PR(4)=1./(2*AMTAU)*(AMTAU**2-AMNUTA**2+AMX**2)
      PR(3)=-PN(3)
      PHSPAC=PHSPAC*(4*PI)*(2*PR(3)/AMTAU)
C
CAM
      ENQ1=(AMX2+AMK**2-AMKZ**2)/(2.*AMX)
      ENQ2=(AMX2-AMK**2+AMKZ**2)/(2.*AMX)
      PPPI=SQRT(abs(ENQ1-AMK)*(ENQ1+AMK))
      PHSPAC=PHSPAC*(4*PI)*(2*PPPI/AMX)
C CHARGED PI MOMENTUM IN RHO REST FRAME
      CALL SPHERA(PPPI,PKC)
      PKC(4)=ENQ1
C NEUTRAL PI MOMENTUM IN RHO REST FRAME
      DO 20 I=1,3
20    PKZ(I)=-PKC(I)
      PKZ(4)=ENQ2
      EXE=(PR(4)+PR(3))/AMX
C PIONS BOOSTED FROM RHO REST FRAME TO TAU REST FRAME
      CALL BOSTR3(EXE,PKC,PKC)
      CALL BOSTR3(EXE,PKZ,PKZ)
      CALL DAM2PI(3,PT,PN,PKC,PKZ,AMPLIT,HV)
      DGAMT=1/(2.*AMTAU)*AMPLIT*PHSPAC
      do 77 k=1,4
        pmult(k,1)=pkc(k)
        pmult(k,2)=pkz(k)
 77   continue
      RETURN
      END
      FUNCTION FPIRK(W)
C ----------------------------------------------------------
c     square of pion form factor
C ----------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
c     COMPLEX FPIKMK
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      FPIRK=CABS(FPIKM(W,AMK,AMKZ))**2
c     FPIRK=CABS(FPIKMK(W,AMK,AMKZ))**2
      END
      SUBROUTINE DWRPH(KTO,PHX)
C
C -------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*4         PHX(4)
      REAL*4 QHOT(4)
C
      DO  9 K=1,4
      QHOT(K)  =0.0
  9   CONTINUE
C CASE OF TAU RADIATIVE DECAYS.
C FILLING OF THE LUND COMMON BLOCK.
        DO 1002 I=1,4
 1002   QHOT(I)=PHX(I)
        IF (QHOT(4).GT.1.E-5) CALL DWLUPH(KTO,QHOT)
        RETURN
      END
      SUBROUTINE DWLUPH(KTO,PHOT)
C---------------------------------------------------------------------
C Lorentz transformation to CMsystem and
C Updating of HEPEVT record
C
C     called by : DEXAY1,(DEKAY1,DEKAY2)
C
C used when radiative corrections in decays are generated
C---------------------------------------------------------------------
C
      REAL  PHOT(4)
      COMMON /TAUPOS/ NP1,NP2
C
C check energy
!
      BIND(C,NAME='tauola_taupos') /TAUPOS/
!
      IF (PHOT(4).LE.0.0) RETURN
C
C position of decaying particle:
      IF((KTO.EQ. 1).OR.(KTO.EQ.11)) THEN
        NPS=NP1
      ELSE
        NPS=NP2
      ENDIF
C
      KTOS=KTO
      IF(KTOS.GT.10) KTOS=KTOS-10
C boost and append photon (gamma is 22)
      CALL TRALO4(KTOS,PHOT,PHOT,AM)
      CALL FILHEP(0,1,22,NPS,NPS,0,0,PHOT,0.0,.TRUE.)
C
      RETURN
      END
      SUBROUTINE DWLUEL(KTO,ISGN,PNU,PWB,PEL,PNE)
C ----------------------------------------------------------------------
C Lorentz transformation to CMsystem and
C Updating of HEPEVT record
C
C ISGN = 1/-1 for tau-/tau+
C
C     called by : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PWB(4),PEL(4),PNE(4)
      COMMON /TAUPOS/ NP1,NP2
C
C position of decaying particle:
!
      BIND(C,NAME='tauola_taupos') /TAUPOS/
!
      IF(KTO.EQ. 1) THEN
        NPS=NP1
      ELSE
        NPS=NP2
      ENDIF
C
C tau neutrino (nu_tau is 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C W boson (W+ is 24)
      CALL TRALO4(KTO,PWB,PWB,AM)
C     CALL FILHEP(0,2,-24*ISGN,NPS,NPS,0,0,PWB,AM,.TRUE.)
C
C electron (e- is 11)
      CALL TRALO4(KTO,PEL,PEL,AM)
      CALL FILHEP(0,1,11*ISGN,NPS,NPS,0,0,PEL,AM,.FALSE.)
C
C anti electron neutrino (nu_e is 12)
      CALL TRALO4(KTO,PNE,PNE,AM)
      CALL FILHEP(0,1,-12*ISGN,NPS,NPS,0,0,PNE,AM,.TRUE.)
C
      RETURN
      END
      SUBROUTINE DWLUMU(KTO,ISGN,PNU,PWB,PMU,PNM)
C ----------------------------------------------------------------------
C Lorentz transformation to CMsystem and
C Updating of HEPEVT record
C
C ISGN = 1/-1 for tau-/tau+
C
C     called by : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PWB(4),PMU(4),PNM(4)
      COMMON /TAUPOS/ NP1,NP2
C
C position of decaying particle:
!
      BIND(C,NAME='tauola_taupos') /TAUPOS/
!
      IF(KTO.EQ. 1) THEN
        NPS=NP1
      ELSE
        NPS=NP2
      ENDIF
C
C tau neutrino (nu_tau is 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C W boson (W+ is 24)
      CALL TRALO4(KTO,PWB,PWB,AM)
C     CALL FILHEP(0,2,-24*ISGN,NPS,NPS,0,0,PWB,AM,.TRUE.)
C
C muon (mu- is 13)
      CALL TRALO4(KTO,PMU,PMU,AM)
      CALL FILHEP(0,1,13*ISGN,NPS,NPS,0,0,PMU,AM,.FALSE.)
C
C anti muon neutrino (nu_mu is 14)
      CALL TRALO4(KTO,PNM,PNM,AM)
      CALL FILHEP(0,1,-14*ISGN,NPS,NPS,0,0,PNM,AM,.TRUE.)
C
      RETURN
      END
      SUBROUTINE DWLUPI(KTO,ISGN,PPI,PNU)
C ----------------------------------------------------------------------
C Lorentz transformation to CMsystem and
C Updating of HEPEVT record
C
C ISGN = 1/-1 for tau-/tau+
C
C     called by : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PPI(4)
      COMMON /TAUPOS/ NP1,NP2
C
C position of decaying particle:
!
      BIND(C,NAME='tauola_taupos') /TAUPOS/
!
      IF(KTO.EQ. 1) THEN
        NPS=NP1
      ELSE
        NPS=NP2
      ENDIF
C
C tau neutrino (nu_tau is 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C charged pi meson (pi+ is 211)
      CALL TRALO4(KTO,PPI,PPI,AM)
      CALL FILHEP(0,1,-211*ISGN,NPS,NPS,0,0,PPI,AM,.TRUE.)
C
      RETURN
      END
      SUBROUTINE DWLURO(KTO,ISGN,PNU,PRHO,PIC,PIZ)
C ----------------------------------------------------------------------
C Lorentz transformation to CMsystem and
C Updating of HEPEVT record
C
C ISGN = 1/-1 for tau-/tau+
C
C     called by : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PRHO(4),PIC(4),PIZ(4)
      COMMON /TAUPOS/ NP1,NP2
C
C position of decaying particle:
!
      BIND(C,NAME='tauola_taupos') /TAUPOS/
!
      IF(KTO.EQ. 1) THEN
        NPS=NP1
      ELSE
        NPS=NP2
      ENDIF
C
C tau neutrino (nu_tau is 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C charged rho meson (rho+ is 213)
      CALL TRALO4(KTO,PRHO,PRHO,AM)
      CALL FILHEP(0,2,-213*ISGN,NPS,NPS,0,0,PRHO,AM,.TRUE.)
C
C charged pi meson (pi+ is 211)
      CALL TRALO4(KTO,PIC,PIC,AM)
      CALL FILHEP(0,1,-211*ISGN,-1,-1,0,0,PIC,AM,.TRUE.)
C
C pi0 meson (pi0 is 111)
      CALL TRALO4(KTO,PIZ,PIZ,AM)
      CALL FILHEP(0,1,111,-2,-2,0,0,PIZ,AM,.TRUE.)
C
      RETURN
      END
      SUBROUTINE DWLUAA(KTO,ISGN,PNU,PAA,PIM1,PIM2,PIPL,JAA)
C ----------------------------------------------------------------------
C Lorentz transformation to CMsystem and
C Updating of HEPEVT record
C
C ISGN = 1/-1 for tau-/tau+
C JAA  = 1 (2) FOR A_1- DECAY TO PI+ 2PI- (PI- 2PI0)
C
C     called by : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PAA(4),PIM1(4),PIM2(4),PIPL(4)
      COMMON /TAUPOS/ NP1,NP2
C
C position of decaying particle:
!
      BIND(C,NAME='tauola_taupos') /TAUPOS/
!
      IF(KTO.EQ. 1) THEN
        NPS=NP1
      ELSE
        NPS=NP2
      ENDIF
C
C tau neutrino (nu_tau is 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C charged a_1 meson (a_1+ is 20213)
      CALL TRALO4(KTO,PAA,PAA,AM)
      CALL FILHEP(0,1,-20213*ISGN,NPS,NPS,0,0,PAA,AM,.TRUE.)
C
C two possible decays of the charged a1 meson
      IF(JAA.EQ.1) THEN
C
C A1  --> PI+ PI-  PI- (or charged conjugate)
C
C pi minus (or c.c.) (pi+ is 211)
        CALL TRALO4(KTO,PIM2,PIM2,AM)
        CALL FILHEP(0,1,-211*ISGN,-1,-1,0,0,PIM2,AM,.TRUE.)
C
C pi minus (or c.c.) (pi+ is 211)
        CALL TRALO4(KTO,PIM1,PIM1,AM)
        CALL FILHEP(0,1,-211*ISGN,-2,-2,0,0,PIM1,AM,.TRUE.)
C
C pi plus (or c.c.) (pi+ is 211)
        CALL TRALO4(KTO,PIPL,PIPL,AM)
        CALL FILHEP(0,1, 211*ISGN,-3,-3,0,0,PIPL,AM,.TRUE.)
C
      ELSE IF (JAA.EQ.2) THEN
C
C A1  --> PI- PI0  PI0 (or charged conjugate)
C
C pi zero (pi0 is 111)
        CALL TRALO4(KTO,PIM2,PIM2,AM)
        CALL FILHEP(0,1,111,-1,-1,0,0,PIM2,AM,.TRUE.)
C
C pi zero (pi0 is 111)
        CALL TRALO4(KTO,PIM1,PIM1,AM)
        CALL FILHEP(0,1,111,-2,-2,0,0,PIM1,AM,.TRUE.)
C
C pi minus (or c.c.) (pi+ is 211)
        CALL TRALO4(KTO,PIPL,PIPL,AM)
        CALL FILHEP(0,1,-211*ISGN,-3,-3,0,0,PIPL,AM,.TRUE.)
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE DWLUKK (KTO,ISGN,PKK,PNU)
C ----------------------------------------------------------------------
C Lorentz transformation to CMsystem and
C Updating of HEPEVT record
C
C ISGN = 1/-1 for tau-/tau+
C
C ----------------------------------------------------------------------
C
      REAL PKK(4),PNU(4)
      COMMON /TAUPOS/ NP1,NP2
C
C position of decaying particle
!
      BIND(C,NAME='tauola_taupos') /TAUPOS/
!
      IF (KTO.EQ.1) THEN
        NPS=NP1
      ELSE
        NPS=NP2
      ENDIF
C
C tau neutrino (nu_tau is 16)
      CALL TRALO4 (KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C K meson (K+ is 321)
      CALL TRALO4 (KTO,PKK,PKK,AM)
      CALL FILHEP(0,1,-321*ISGN,NPS,NPS,0,0,PKK,AM,.TRUE.)
C
      RETURN
      END
      SUBROUTINE DWLUKS(KTO,ISGN,PNU,PKS,PKK,PPI,JKST)
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
C ----------------------------------------------------------------------
C Lorentz transformation to CMsystem and
C Updating of HEPEVT record
C
C ISGN = 1/-1 for tau-/tau+
C JKST=10 (20) corresponds to K0B pi- (K- pi0) decay
C
C ----------------------------------------------------------------------
C
      REAL  PNU(4),PKS(4),PKK(4),PPI(4),XIO(1)
      COMMON /TAUPOS/ NP1,NP2
C
C position of decaying particle
!
      BIND(C,NAME='tauola_taukle') /TAUKLE/
      BIND(C,NAME='tauola_taupos') /TAUPOS/
!
      IF(KTO.EQ. 1) THEN
        NPS=NP1
      ELSE
        NPS=NP2
      ENDIF
C
C tau neutrino (nu_tau is 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C charged K* meson (K*+ is 323)
      CALL TRALO4(KTO,PKS,PKS,AM)
      CALL FILHEP(0,1,-323*ISGN,NPS,NPS,0,0,PKS,AM,.TRUE.)
C
C two possible decay modes of charged K*
      IF(JKST.EQ.10) THEN
C
C K*- --> pi- K0B (or charged conjugate)
C
C charged pi meson  (pi+ is 211)
        CALL TRALO4(KTO,PPI,PPI,AM)
        CALL FILHEP(0,1,-211*ISGN,-1,-1,0,0,PPI,AM,.TRUE.)
C
        BRAN=BRK0B
        IF (ISGN.EQ.-1) BRAN=BRK0
C K0 --> K0_long (is 130) / K0_short (is 310) = 1/1
        CALL RANMAR(XIO,1)
        IF(XIO(1).GT.BRAN) THEN
          K0TYPE = 130
        ELSE
          K0TYPE = 310
        ENDIF
C
        CALL TRALO4(KTO,PKK,PKK,AM)
        CALL FILHEP(0,1,K0TYPE,-2,-2,0,0,PKK,AM,.TRUE.)
C
      ELSE IF(JKST.EQ.20) THEN
C
C K*- --> pi0 K-
C
C pi zero (pi0 is 111)
        CALL TRALO4(KTO,PPI,PPI,AM)
        CALL FILHEP(0,1,111,-1,-1,0,0,PPI,AM,.TRUE.)
C
C charged K meson (K+ is 321)
        CALL TRALO4(KTO,PKK,PKK,AM)
        CALL FILHEP(0,1,-321*ISGN,-2,-2,0,0,PKK,AM,.TRUE.)
C
      ENDIF
C
      RETURN
      END
      SUBROUTINE DWLNEW(KTO,ISGN,PNU,PWB,PNPI,MODE)
C ----------------------------------------------------------------------
C Lorentz transformation to CMsystem and
C Updating of HEPEVT record
C
C ISGN = 1/-1 for tau-/tau+
C
C     called by : DEXAY,(DEKAY1,DEKAY2)
C ----------------------------------------------------------------------
C
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / TAUDCD /IDFFIN(9,NMODE),MULPIK(NMODE)
     &                ,NAMES
      COMMON /TAUPOS/ NP1,NP2
      CHARACTER NAMES(NMODE)*31
      REAL  PNU(4),PWB(4),PNPI(4,9)
      REAL  PPI(4)
C
!
      BIND(C,NAME='tauola_taudcd') /TAUDCD/
      BIND(C,NAME='tauola_taupos') /TAUPOS/
!
      JNPI=MODE-7
C position of decaying particle
      IF(KTO.EQ. 1) THEN
        NPS=NP1
      ELSE
        NPS=NP2
      ENDIF
C
C tau neutrino (nu_tau is 16)
      CALL TRALO4(KTO,PNU,PNU,AM)
      CALL FILHEP(0,1,16*ISGN,NPS,NPS,0,0,PNU,AM,.TRUE.)
C
C W boson (W+ is 24)
      CALL TRALO4(KTO,PWB,PWB,AM)
      CALL FILHEP(0,1,-24*ISGN,NPS,NPS,0,0,PWB,AM,.TRUE.)
C
C multi pi mode JNPI
C
C get multiplicity of mode JNPI
      ND=MULPIK(JNPI)
      DO I=1,ND
        KFPI=LUNPIK(IDFFIN(I,JNPI),-ISGN)
C for charged conjugate case, change charged pions only
C        IF(KFPI.NE.111)KFPI=KFPI*ISGN
        DO J=1,4
          PPI(J)=PNPI(J,I)
        END DO
        CALL TRALO4(KTO,PPI,PPI,AM)
        CALL FILHEP(0,1,KFPI,-I,-I,0,0,PPI,AM,.TRUE.)
      END DO
C
      RETURN
      END
      SUBROUTINE ROTOD2(PH1,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PVEC(4),QVEC(4),RVEC(4)
C
      PHI=PH1
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)= CS*RVEC(1)+SN*RVEC(3)
      QVEC(2)=RVEC(2)
      QVEC(3)=-SN*RVEC(1)+CS*RVEC(3)
      QVEC(4)=RVEC(4)
      RETURN
      END
      SUBROUTINE ROTOD3(PH1,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION PVEC(4),QVEC(4),RVEC(4)
      PHI=PH1
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)= CS*RVEC(1)-SN*RVEC(2)
      QVEC(2)= SN*RVEC(1)+CS*RVEC(2)
      QVEC(3)=RVEC(3)
      QVEC(4)=RVEC(4)
      END
      SUBROUTINE BOSTR3(EXE,PVEC,QVEC)
C ----------------------------------------------------------------------
C BOOST ALONG Z AXIS, EXE=EXP(ETA), ETA= HIPERBOLIC VELOCITY.
C
C     USED BY : TAUOLA KORALZ (?)
C ----------------------------------------------------------------------
      REAL*4 PVEC(4),QVEC(4),RVEC(4)
C
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      RPL=RVEC(4)+RVEC(3)
      RMI=RVEC(4)-RVEC(3)
      QPL=RPL*EXE
      QMI=RMI/EXE
      QVEC(1)=RVEC(1)
      QVEC(2)=RVEC(2)
      QVEC(3)=(QPL-QMI)/2
      QVEC(4)=(QPL+QMI)/2
      END
      SUBROUTINE BOSTD3(EXE,PVEC,QVEC)
C ----------------------------------------------------------------------
C BOOST ALONG Z AXIS, EXE=EXP(ETA), ETA= HIPERBOLIC VELOCITY.
C
C     USED BY : KORALZ RADKOR
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PVEC(4),QVEC(4),RVEC(4)
C
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      RPL=RVEC(4)+RVEC(3)
      RMI=RVEC(4)-RVEC(3)
      QPL=RPL*EXE
      QMI=RMI/EXE
      QVEC(1)=RVEC(1)
      QVEC(2)=RVEC(2)
      QVEC(3)=(QPL-QMI)/2
      QVEC(4)=(QPL+QMI)/2
      RETURN
      END
      SUBROUTINE ROTOR2(PH1,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : TAUOLA
C ----------------------------------------------------------------------
      IMPLICIT REAL*4(A-H,O-Z)
      REAL*4 PVEC(4),QVEC(4),RVEC(4)
C
      PHI=PH1
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)= CS*RVEC(1)+SN*RVEC(3)
      QVEC(2)=RVEC(2)
      QVEC(3)=-SN*RVEC(1)+CS*RVEC(3)
      QVEC(4)=RVEC(4)
      END
      SUBROUTINE ROTOR3(PHI,PVEC,QVEC)
C ----------------------------------------------------------------------
C
C     USED BY : TAUOLA
C ----------------------------------------------------------------------
      REAL*4 PVEC(4),QVEC(4),RVEC(4)
C
      CS=COS(PHI)
      SN=SIN(PHI)
      DO 10 I=1,4
  10  RVEC(I)=PVEC(I)
      QVEC(1)= CS*RVEC(1)-SN*RVEC(2)
      QVEC(2)= SN*RVEC(1)+CS*RVEC(2)
      QVEC(3)=RVEC(3)
      QVEC(4)=RVEC(4)
      END
      SUBROUTINE SPHERD(R,X)
C ----------------------------------------------------------------------
C GENERATES UNIFORMLY THREE-VECTOR X ON SPHERE  OF RADIUS R
C DOUBLE PRECISON VERSION OF SPHERA
C ----------------------------------------------------------------------
      REAL*8  R,X(4),PI,COSTH,SINTH
      REAL*4 RRR(2)
      DATA PI /3.141592653589793238462643D0/
C
      CALL RANMAR(RRR,2)
      COSTH=-1+2*RRR(1)
      SINTH=SQRT(1 -COSTH**2)
      X(1)=R*SINTH*COS(2*PI*RRR(2))
      X(2)=R*SINTH*SIN(2*PI*RRR(2))
      X(3)=R*COSTH
      RETURN
      END
      SUBROUTINE ROTPOX(THET,PHI,PP)
      IMPLICIT REAL*8 (A-H,O-Z)
C ----------------------------------------------------------------------
C
C ----------------------------------------------------------------------
      DIMENSION PP(4)
C
      CALL ROTOD2(THET,PP,PP)
      CALL ROTOD3( PHI,PP,PP)
      RETURN
      END
      SUBROUTINE SPHERA(R,X)
C ----------------------------------------------------------------------
C GENERATES UNIFORMLY THREE-VECTOR X ON SPHERE  OF RADIUS R
C
C     called by : DPHSxx,DADMPI,DADMKK
C ----------------------------------------------------------------------
      REAL  X(4)
      REAL*4 RRR(2)
      DATA PI /3.141592653589793238462643/
C
      CALL RANMAR(RRR,2)
      COSTH=-1.+2.*RRR(1)
      SINTH=SQRT(1.-COSTH**2)
      X(1)=R*SINTH*COS(2*PI*RRR(2))
      X(2)=R*SINTH*SIN(2*PI*RRR(2))
      X(3)=R*COSTH
      RETURN
      END
      SUBROUTINE ROTPOL(THET,PHI,PP)
C ----------------------------------------------------------------------
C
C     called by : DADMAA,DPHSAA
C ----------------------------------------------------------------------
      REAL  PP(4)
C
      CALL ROTOR2(THET,PP,PP)
      CALL ROTOR3( PHI,PP,PP)
      RETURN
      END
      SUBROUTINE RANMAR(RVEC,LENV)
      DIMENSION RVEC(*)
      CALL TAUOLA_RANDOM(RVEC,LENV)
      END SUBROUTINE RANMAR
      FUNCTION DILOGT(X)
C     *****************
      IMPLICIT REAL*8(A-H,O-Z)
CERN      C304      VERSION    29/07/71 DILOG        59                C
      Z=-1.64493406684822
      IF(X .LT.-1.0) GO TO 1
      IF(X .LE. 0.5) GO TO 2
      IF(X .EQ. 1.0) GO TO 3
      IF(X .LE. 2.0) GO TO 4
      Z=3.2898681336964
    1 T=1.0/X
      S=-0.5
      Z=Z-0.5* LOG(ABS(X))**2
      GO TO 5
    2 T=X
      S=0.5
      Z=0.
      GO TO 5
    3 DILOGT=1.64493406684822
      RETURN
    4 T=1.0-X
      S=-0.5
      Z=1.64493406684822 - LOG(X)* LOG(ABS(T))
    5 Y=2.66666666666666 *T+0.66666666666666
      B=      0.00000 00000 00001
      A=Y*B  +0.00000 00000 00004
      B=Y*A-B+0.00000 00000 00011
      A=Y*B-A+0.00000 00000 00037
      B=Y*A-B+0.00000 00000 00121
      A=Y*B-A+0.00000 00000 00398
      B=Y*A-B+0.00000 00000 01312
      A=Y*B-A+0.00000 00000 04342
      B=Y*A-B+0.00000 00000 14437
      A=Y*B-A+0.00000 00000 48274
      B=Y*A-B+0.00000 00001 62421
      A=Y*B-A+0.00000 00005 50291
      B=Y*A-B+0.00000 00018 79117
      A=Y*B-A+0.00000 00064 74338
      B=Y*A-B+0.00000 00225 36705
      A=Y*B-A+0.00000 00793 87055
      B=Y*A-B+0.00000 02835 75385
      A=Y*B-A+0.00000 10299 04264
      B=Y*A-B+0.00000 38163 29463
      A=Y*B-A+0.00001 44963 00557
      B=Y*A-B+0.00005 68178 22718
      A=Y*B-A+0.00023 20021 96094
      B=Y*A-B+0.00100 16274 96164
      A=Y*B-A+0.00468 63619 59447
      B=Y*A-B+0.02487 93229 24228
      A=Y*B-A+0.16607 30329 27855
      A=Y*A-B+1.93506 43008 6996
      DILOGT=S*T*(A-B)+Z
      RETURN
C=======================================================================
C===================END OF CPC PART ====================================
C=======================================================================
      END
      COMPLEX FUNCTION F3PI(IFORM,QQ,SA,SB)
C.......................................................................
C.
C. F3PI - 1 version of a1 form factor, used in TAUOLA
C.
C. Inputs    : None
C.           :
C. Outputs   : None
C.
C. COMMON    : None
C.
C. Calls     :
C. Called    : by FORM1-FORM3 in $C_CVSSRC/korb/koralb/formf.F
C. Author    : Alan Weinstein 2/98
C.
C. Detailed description
C.   First determine whether we are doing pi-2pi0 or 3pi.
C.   Then implement full form-factor from fit:
C.   [(rho-pi S-wave) + (rho-prim-pi S-wave) +
C.    (rho-pi D-wave) + (rho-prim-pi D-wave) +
C.    (f2 pi D-wave) + (sigmapi S-wave) + (f0pi S-wave)]
C.   based on fit to pi-2pi0 by M. Schmidler, CBX 97-64-Update (4/22/98)
C.   All the parameters in this routine are hard-coded!!
C.
C.......................................................................
* -------------------- Argument declarations ---------------
      INTEGER IFORM
      REAL QQ,SA,SB
* -------------------- EXTERNAL declarations ---------------
*
* -------------------- SEQUENCE declarations ---------------
*
* -------------------- Local    declarations ---------------
*
      CHARACTER*(*) CRNAME
      PARAMETER(    CRNAME = 'F3PI' )
*
      INTEGER IFIRST,IDK
      REAL MRO,GRO,MRP,GRP,MF2,GF2,MF0,GF0,MSG,GSG
      REAL M1,M2,M3,M1SQ,M2SQ,M3SQ,MPIZ,MPIC
      REAL S1,S2,S3,R,PI
      REAL F134,F150,F15A,F15B,F167
      REAL F34A,F34B,F35,F35A,F35B,F36A,F36B
      COMPLEX BT1,BT2,BT3,BT4,BT5,BT6,BT7
      COMPLEX FRO1,FRO2,FRP1,FRP2
      COMPLEX FF21,FF22,FF23,FSG1,FSG2,FSG3,FF01,FF02,FF03
      COMPLEX FORMA1
* -------------------- SAVE     declarations ---------------
*
* -------------------- DATA  initializations ---------------
*
      DATA IFIRST/0/
* ----------------- Executable code starts here ------------
*
C. Hard-code the fit parameters:
      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
C rho, rhoprime, f2(1275), f0(1186), sigma(made up!)
        MRO = 0.7743
        GRO = 0.1491
        MRP = 1.370
        GRP = 0.386
        MF2 = 1.275
        GF2 = 0.185
        MF0 = 1.186
        GF0 = 0.350
        MSG = 0.860
        GSG = 0.880
        MPIZ = PKORB(1,7)
        MPIC = PKORB(1,8)
C Fit coefficients for each of the contributions:
        PI = 3.14159
        BT1 = CMPLX(1.,0.)
        BT2 = CMPLX(0.12,0.)*CEXP(CMPLX(0., 0.99*PI))
        BT3 = CMPLX(0.37,0.)*CEXP(CMPLX(0.,-0.15*PI))
        BT4 = CMPLX(0.87,0.)*CEXP(CMPLX(0., 0.53*PI))
        BT5 = CMPLX(0.71,0.)*CEXP(CMPLX(0., 0.56*PI))
        BT6 = CMPLX(2.10,0.)*CEXP(CMPLX(0., 0.23*PI))
        BT7 = CMPLX(0.77,0.)*CEXP(CMPLX(0.,-0.54*PI))
        CALL TAUOLA_PRINT(
     $    'In F3pi: add (rho-pi S-wave) + (rhop-pi S-wave) +'//CHAR(0))
!       PRINT *,' In F3pi: add (rho-pi S-wave) + (rhop-pi S-wave) +'
        CALL TAUOLA_PRINT('(rho-pi D-wave) + (rhop-pi D-wave) +'//CHAR(0))
!       PRINT *,'              (rho-pi D-wave) + (rhop-pi D-wave) +'
        CALL TAUOLA_PRINT(
     $    '(f2 pi D-wave) + (sigmapi S-wave) + (f0pi S-wave)'//CHAR(0))
!       PRINT *,'   (f2 pi D-wave) + (sigmapi S-wave) + (f0pi S-wave)'
      END IF
C Initialize to 0:
      F3PI = CMPLX(0.,0.)
C.   First determine whether we are doing pi-2pi0 or 3pi.
C     PKORB is set up to remember what flavor of 3pi it gave to KORALB,
C     since KORALB doesnt bother to remember!!
      R = PKORB(4,11)
      IF (R.EQ.0.) THEN
C it is 2pi0pi-
        IDK = 1
        M1 = MPIZ
        M2 = MPIZ
        M3 = MPIC
      ELSE
C it is 3pi
        IDK = 2
        M1 = MPIC
        M2 = MPIC
        M3 = MPIC
      END IF
      M1SQ = M1*M1
      M2SQ = M2*M2
      M3SQ = M3*M3
C.   Then implement full form-factor from fit:
C.   [(rho-pi S-wave) + (rho-prim-pi S-wave) +
C.    (rho-pi D-wave) + (rho-prim-pi D-wave) +
C.    (f2 pi D-wave) + (sigmapi S-wave) + (f0pi S-wave)]
C.   based on fit to pi-2pi0 by M. Schmidler, CBX 97-64-Update (4/22/98)
C Note that for FORM1, the arguments are S1, S2;
C           for FORM2, the arguments are S2, S1;
C           for FORM3, the arguments are S3, S1.
C Here, we implement FORM1 and FORM2 at the same time,
C so the above switch is just what we need!
      IF (IFORM.EQ.1.OR.IFORM.EQ.2) THEN
        S1 = SA
        S2 = SB
        S3 = QQ-SA-SB+M1SQ+M2SQ+M3SQ
        IF (S3.LE.0..OR.S2.LE.0.) RETURN
        IF (IDK.EQ.1) THEN
C it is 2pi0pi-
C Lorentz invariants for all the contributions:
          F134 = -(1./3.)*((S3-M3SQ)-(S1-M1SQ))
          F150 =  (1./18.)*(QQ-M3SQ+S3)*(2.*M1SQ+2.*M2SQ-S3)/S3
          F167 =  (2./3.)
C Breit Wigners for all the contributions:
          FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1)
          FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1)
          FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1)
          FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1)
          FF23 = BWIGML(S3,MF2,GF2,M1,M2,2)
          FSG3 = BWIGML(S3,MSG,GSG,M1,M2,0)
          FF03 = BWIGML(S3,MF0,GF0,M1,M2,0)
          F3PI = BT1*FRO1+BT2*FRP1+
     1       BT3*CMPLX(F134,0.)*FRO2+BT4*CMPLX(F134,0.)*FRP2+
     1       BT5*CMPLX(F150,0.)*FF23+
     1       BT6*CMPLX(F167,0.)*FSG3+BT7*CMPLX(F167,0.)*FF03
C         F3PI = FPIKM(SQRT(S1),M2,M3)
        ELSEIF (IDK.EQ.2) THEN
C it is 3pi
C Lorentz invariants for all the contributions:
          F134 = -(1./3.)*((S3-M3SQ)-(S1-M1SQ))
          F15A = -(1./2.)*((S2-M2SQ)-(S3-M3SQ))
          F15B = -(1./18.)*(QQ-M2SQ+S2)*(2.*M1SQ+2.*M3SQ-S2)/S2
          F167 = -(2./3.)
C Breit Wigners for all the contributions:
          FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1)
          FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1)
          FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1)
          FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1)
          FF21 = BWIGML(S1,MF2,GF2,M2,M3,2)
          FF22 = BWIGML(S2,MF2,GF2,M3,M1,2)
          FSG2 = BWIGML(S2,MSG,GSG,M3,M1,0)
          FF02 = BWIGML(S2,MF0,GF0,M3,M1,0)
          F3PI = BT1*FRO1+BT2*FRP1+
     1       BT3*CMPLX(F134,0.)*FRO2+BT4*CMPLX(F134,0.)*FRP2
     1      -BT5*CMPLX(F15A,0.)*FF21-BT5*CMPLX(F15B,0.)*FF22
     1      -BT6*CMPLX(F167,0.)*FSG2-BT7*CMPLX(F167,0.)*FF02
C         F3PI = FPIKM(SQRT(S1),M2,M3)
        END IF
      ELSE IF (IFORM.EQ.3) THEN
        S3 = SA
        S1 = SB
        S2 = QQ-SA-SB+M1SQ+M2SQ+M3SQ
        IF (S1.LE.0..OR.S2.LE.0.) RETURN
        IF (IDK.EQ.1) THEN
C it is 2pi0pi-
C Lorentz invariants for all the contributions:
          F34A = (1./3.)*((S2-M2SQ)-(S3-M3SQ))
          F34B = (1./3.)*((S3-M3SQ)-(S1-M1SQ))
          F35  =-(1./2.)*((S1-M1SQ)-(S2-M2SQ))
C Breit Wigners for all the contributions:
          FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1)
          FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1)
          FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1)
          FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1)
          FF23 = BWIGML(S3,MF2,GF2,M1,M2,2)
          F3PI =
     1       BT3*(CMPLX(F34A,0.)*FRO1+CMPLX(F34B,0.)*FRO2)+
     1       BT4*(CMPLX(F34A,0.)*FRP1+CMPLX(F34B,0.)*FRP2)+
     1       BT5*CMPLX(F35,0.)*FF23
C         F3PI = CMPLX(0.,0.)
        ELSEIF (IDK.EQ.2) THEN
C it is 3pi
C Lorentz invariants for all the contributions:
          F34A = (1./3.)*((S2-M2SQ)-(S3-M3SQ))
          F34B = (1./3.)*((S3-M3SQ)-(S1-M1SQ))
          F35A = -(1./18.)*(QQ-M1SQ+S1)*(2.*M2SQ+2.*M3SQ-S1)/S1
          F35B =  (1./18.)*(QQ-M2SQ+S2)*(2.*M3SQ+2.*M1SQ-S2)/S2
          F36A = -(2./3.)
          F36B =  (2./3.)
C Breit Wigners for all the contributions:
          FRO1 = BWIGML(S1,MRO,GRO,M2,M3,1)
          FRP1 = BWIGML(S1,MRP,GRP,M2,M3,1)
          FRO2 = BWIGML(S2,MRO,GRO,M3,M1,1)
          FRP2 = BWIGML(S2,MRP,GRP,M3,M1,1)
          FF21 = BWIGML(S1,MF2,GF2,M2,M3,2)
          FF22 = BWIGML(S2,MF2,GF2,M3,M1,2)
          FSG1 = BWIGML(S1,MSG,GSG,M2,M3,0)
          FSG2 = BWIGML(S2,MSG,GSG,M3,M1,0)
          FF01 = BWIGML(S1,MF0,GF0,M2,M3,0)
          FF02 = BWIGML(S2,MF0,GF0,M3,M1,0)
          F3PI =
     1       BT3*(CMPLX(F34A,0.)*FRO1+CMPLX(F34B,0.)*FRO2)+
     1       BT4*(CMPLX(F34A,0.)*FRP1+CMPLX(F34B,0.)*FRP2)
     1      -BT5*(CMPLX(F35A,0.)*FF21+CMPLX(F35B,0.)*FF22)
     1      -BT6*(CMPLX(F36A,0.)*FSG1+CMPLX(F36B,0.)*FSG2)
     1      -BT7*(CMPLX(F36A,0.)*FF01+CMPLX(F36B,0.)*FF02)
C         F3PI = CMPLX(0.,0.)
        END IF
      END IF
C Add overall a1/a1prime:
      FORMA1 = FA1A1P(QQ)
      F3PI = F3PI*FORMA1
      RETURN
      END
      COMPLEX FUNCTION BWIGML(S,M,G,M1,M2,L)
C **********************************************************
C     L-WAVE BREIT-WIGNER
C **********************************************************
      REAL S,M,G,M1,M2
      INTEGER L,IPOW
      REAL MSQ,W,WGS,MP,MM,QS,QM
      MP = (M1+M2)**2
      MM = (M1-M2)**2
      MSQ = M*M
      W = SQRT(S)
      WGS = 0.0
      IF (W.GT.(M1+M2)) THEN
        QS=SQRT(ABS((S   -MP)*(S   -MM)))/W
        QM=SQRT(ABS((MSQ -MP)*(MSQ -MM)))/M
        IPOW = 2*L+1
        WGS=G*(MSQ/W)*(QS/QM)**IPOW
      ENDIF
      BWIGML=CMPLX(MSQ,0.)/CMPLX(MSQ-S,-WGS)
      RETURN
      END
      COMPLEX FUNCTION FA1A1P(XMSQ)
C     ==================================================================
C     complex form-factor for a1+a1prime.                       AJW 1/98
C     ==================================================================
      REAL XMSQ
      REAL XM1,XG1,XM2,XG2,XM1SQ,XM2SQ,GG1,GG2,GF,FG1,FG2
      COMPLEX BET,F1,F2
      INTEGER IFIRST/0/
      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
C The user may choose masses and widths that differ from nominal:
        XM1 = PKORB(1,10)
        XG1 = PKORB(2,10)
        XM2 = PKORB(1,17)
        XG2 = PKORB(2,17)
        BET = CMPLX(PKORB(3,17),0.)
C scale factors relative to nominal:
        GG1 = XM1*XG1/(1.3281*0.806)
        GG2 = XM2*XG2/(1.3281*0.806)
        XM1SQ = XM1*XM1
        XM2SQ = XM2*XM2
      END IF
      GF = WGA1(XMSQ)
      FG1 = GG1*GF
      FG2 = GG2*GF
      F1 = CMPLX(-XM1SQ,0.0)/CMPLX(XMSQ-XM1SQ,FG1)
      F2 = CMPLX(-XM2SQ,0.0)/CMPLX(XMSQ-XM2SQ,FG2)
      FA1A1P = F1+BET*F2
      RETURN
      END
      FUNCTION WGA1(QQ)
C mass-dependent M*Gamma of a1 through its decays to
C.   [(rho-pi S-wave) + (rho-pi D-wave) +
C.    (f2 pi D-wave) + (f0pi S-wave)]
C.  AND simple K*K S-wave
      REAL QQ,WGA1
      DOUBLE PRECISION MKST,MK,MK1SQ,MK2SQ,C3PI,CKST
      DOUBLE PRECISION S,WG3PIC,WG3PIN,GKST
      INTEGER IFIRST
C-----------------------------------------------------------------------
C
      IF (IFIRST.NE.987) THEN
        IFIRST = 987
C
C Contribution to M*Gamma(m(3pi)^2) from S-wave K*K:
        MKST = 0.894D0
        MK = 0.496D0
        MK1SQ = (MKST+MK)**2
        MK2SQ = (MKST-MK)**2
C coupling constants squared:
        C3PI = 0.2384D0**2
        CKST = 4.7621D0**2*C3PI
      END IF
C-----------------------------------------------------------------------
C Parameterization of numerical integral of total width of a1 to 3pi.
C From M. Schmidtler, CBX-97-64-Update.
      S = DBLE(QQ)
      WG3PIC = WGA1C(S)
      WG3PIN = WGA1N(S)
C Contribution to M*Gamma(m(3pi)^2) from S-wave K*K, if above threshold
      GKST = 0.D0
      IF (S.GT.MK1SQ) GKST = SQRT((S-MK1SQ)*(S-MK2SQ))/(2.*S)
      WGA1 = SNGL(C3PI*(WG3PIC+WG3PIN)+CKST*GKST)
      RETURN
      END
      DOUBLE PRECISION FUNCTION WGA1C(S)
C
C parameterization of m*Gamma(m^2) for pi-2pi0 system
C
      DOUBLE PRECISION S,STH,Q0,Q1,Q2,P0,P1,P2,P3,P4,G1_IM
C
      PARAMETER(Q0 =   5.80900D0,Q1 =  -3.00980D0,Q2 =   4.57920D0,
     1          P0 = -13.91400D0,P1 =  27.67900D0,P2 = -13.39300D0,
     2          P3 =   3.19240D0,P4 =  -0.10487D0)
C
      PARAMETER (STH   = 0.1753D0)
C---------------------------------------------------------------------
      IF(S.LT.STH) THEN
       G1_IM = 0.D0
      ELSEIF((S.GT.STH).AND.(S.LT.0.823D0)) THEN
       G1_IM = Q0*(S-STH)**3*(1. + Q1*(S-STH) + Q2*(S-STH)**2)
      ELSE
       G1_IM = P0 + P1*S + P2*S**2+ P3*S**3 + P4*S**4
      ENDIF
      WGA1C = G1_IM
      RETURN
      END
      DOUBLE PRECISION FUNCTION WGA1N(S)
C
C parameterization of m*Gamma(m^2) for pi-pi+pi- system
C
      DOUBLE PRECISION S,STH,Q0,Q1,Q2,P0,P1,P2,P3,P4,G1_IM
C
      PARAMETER(Q0 =   6.28450D0,Q1 =  -2.95950D0,Q2 =   4.33550D0,
     1          P0 = -15.41100D0,P1 =  32.08800D0,P2 = -17.66600D0,
     2          P3 =   4.93550D0,P4 =  -0.37498D0)
C
      PARAMETER (STH   = 0.1676D0)
C---------------------------------------------------------------------
      IF(S.LT.STH) THEN
       G1_IM = 0.D0
      ELSEIF((S.GT.STH).AND.(S.LT.0.823D0)) THEN
       G1_IM = Q0*(S-STH)**3*(1. + Q1*(S-STH) + Q2*(S-STH)**2)
      ELSE
       G1_IM = P0 + P1*S + P2*S**2+ P3*S**3 + P4*S**4
      ENDIF
      WGA1N = G1_IM
      RETURN
      END
      REAL FUNCTION PKORB(IF1,IF2)
**********************************************************************
*
* This function returns a real value
* needed in the 1 version of KORALB/TAUOLA
* corresponding to a mass, width, mixing amplitude, or branching fraction
* depending on whether IF1 = 1, 2, 3, 4 respectively.
* The idea is to make minimal mods to the 3-rd party KORALB/TAUOLA code,
* so this function supplies all the 1-specific parameters.
*
*  Alan Weinstein, ajw, 11/97
**********************************************************************
* Arguments:
      INTEGER IF1   ! input, flag for type of data required
      INTEGER IF2   ! input, flag for type of data required
* MC info
*#include "seq/clinc/qqpars.inc"
*#include "seq/clinc/qqprop.inc"
*#include "qqlib/seq/qqbrat.inc"
      INTEGER            JAK1,JAK2,JAKP,JAKM,KTOM
      COMMON / JAKI   /  JAK1,JAK2,JAKP,JAKM,KTOM
      REAL*4 RRR(1)
      REAL PARM(4,100)
      integer imixpp(300)
      INTEGER INIT,I,J
      REAL C1270,C1402,A1270_KSPI,A1270_KRHO,A1402_KSPI,A1402_KRHO
      REAL CG1,CG2,R,BRA1,BRKS
      SAVE INIT,PARM
      DATA INIT/0/
**********************************************************************
* Initialize return variable:
!
      BIND(C,NAME='tauola_jaki') /JAKI/
!
      PKORB = 0.
**********************************************************************
* Initialize:
      IF (INIT.EQ.0) THEN
        INIT = 1
C WARNING: Isospin symmetry enforced, cleo or babar were not enforcing it.
C          Simplification to be used for precision tau decay simulations.
        BRA1=0.0
        BRKS=0.0
C        CALL VZERO(PARM,400)
        DO I=1,4
        DO J=1,100
          PARM(I,J) = 0
        END DO
        END DO
C Youd better be using korb.dec, NOT decay.dec!!!!
C masses (needed in dist/inimas, formf/form*, etc)
        PARM(1, 1) = 1.777000   ! TAU
        PARM(1, 2) = 0.         ! NUTA
        PARM(1, 3) = 0.000511   ! EL
        PARM(1, 4) = 0.         ! NUEL
        PARM(1, 5) = 0.105658   ! MU
        PARM(1, 6) = 0.         ! NUMU
        PARM(1, 7) = 0.134976   ! PIZ
        PARM(1, 8) = 0.139570   ! PI+
        PARM(1, 9) = 0.769900   ! RHO+
        PARM(1,10) = 1.275000   ! A1+
        PARM(1,11) = 0.493677   ! K+
        PARM(1,12) = 0.497670   ! KZ
        PARM(1,13) = 0.891590   ! K*+
        PARM(1,14) = 0.781940   ! OMEG
        PARM(1,15) = 1.370000   ! RHOP+
        PARM(1,16) = 1.700000   ! K*P+
        PARM(1,17) = 1.461000   ! A1P+
        PARM(1,18) = 1.300000   ! PIP+
        PARM(1,19) = 1.270000   ! K1A+
        PARM(1,20) = 1.402000   ! K1B+
        PARM(1,21) = 1.465000   ! RHOPP+
        PARM(1,22) = 1.700000   ! RHOPPP+
C widths (needed in dist/inimas, formf/form*, etc)
        PARM(2, 1) = 0.          ! TAU
        PARM(2, 2) = 0.          ! NUTA
        PARM(2, 3) = 0.          ! EL
        PARM(2, 4) = 0.          ! NUEL
        PARM(2, 5) = 0.          ! MU
        PARM(2, 6) = 0.          ! NUMU
        PARM(2, 7) = 0.          ! PIZ
        PARM(2, 8) = 0.          ! PI+
        PARM(2, 9) = 0.1512      ! RHO+
        PARM(2,10) = 0.700       ! A1+
        PARM(2,11) = 0.          ! K+
        PARM(2,12) = 0.          ! KZ
        PARM(2,13) = 0.0498      ! K*+
        PARM(2,14) = 0.00843     ! OMEG
        PARM(2,15) = 0.510       ! RHOP+
        PARM(2,16) = 0.235       ! K*P+
        PARM(2,17) = 0.250       ! A1P+
        PARM(2,18) = 0.400       ! PIP+
        PARM(2,19) = 0.090       ! K1A+
        PARM(2,20) = 0.174       ! K1B+
        PARM(2,21) = 0.310       ! RHOPP+
        PARM(2,22) = 0.235       ! RHOPPP+
C Now store mixing parameters for 2pi and 4pi FFs
C needed in tauola/fpik, tauola/bwigs, formf/form* , formf/curr :
        PARM(3,15) = -0.145
        IMIXPP(205)=1
        IMIXPP(207)=1
        IMIXPP(209)=1
        IMIXPP(211)=1
        IMIXPP(201)=1
        IMIXPP(203)=1
        IMIXPP(213)=1
        IMIXPP(215)=1
        IF (IMIXPP(205).NE.0) PARM(3,15) = -0.110
        IF (IMIXPP(207).NE.0) PARM(3,16) = -0.038
        IF (IMIXPP(209).NE.0) PARM(3,17) = 0.00
        IF (IMIXPP(211).NE.0) PARM(3,18) = 0.00
        IF (IMIXPP(201).NE.0) PARM(3,19) = 1.0
        IF (IMIXPP(203).NE.0) PARM(3,20) = 0.8
        IF (IMIXPP(213).NE.0) PARM(3,21) = -0.110
        IF (IMIXPP(215).NE.0) PARM(3,22) = -0.110
        CALL TAUOLA_PRINT('KORB: rho/rhop -> pi-pi0 mixing:'//CHAR(0))
!       PRINT *,' KORB: rho/rhop -> pi-pi0 mixing:'
        CALL TAUOLA_PRINT('KORB: rho   ='//CHAR(0))
!       PRINT *,' KORB: rho   =',PARM(1,9) ,PARM(2,9)
        CALL TAUOLA_PRINT('KORB: rhop  ='//CHAR(0))
!       PRINT *,' KORB: rhop  =',PARM(1,15),PARM(2,15),PARM(3,15)
        CALL TAUOLA_PRINT('KORB: K*/K*prime -> Kpi mixing:'//CHAR(0))
!       PRINT *,' KORB: K*/K*prime -> Kpi mixing:'
        CALL TAUOLA_PRINT('KORB: kstp  ='//CHAR(0))
!       PRINT *,' KORB: kstp  =',PARM(1,16),PARM(2,16),PARM(3,16)
        CALL TAUOLA_PRINT('KORB: a1/a1prime -> 3pi, KKpi mixing:'//CHAR(0))
!       PRINT *,' KORB: a1/a1prime -> 3pi, KKpi mixing:'
        CALL TAUOLA_PRINT('KORB: a1    ='//CHAR(0))
!       PRINT *,' KORB: a1    =',PARM(1,10),PARM(2,10)
        CALL TAUOLA_PRINT('KORB: a1prim='//CHAR(0))
!       PRINT *,' KORB: a1prim=',PARM(1,17),PARM(2,17),PARM(3,17)
        CALL TAUOLA_PRINT('KORB: K1A/K1B -> Kpipi mixing:'//CHAR(0))
!       PRINT *,' KORB: K1A/K1B -> Kpipi mixing:'
        CALL TAUOLA_PRINT('KORB: K1A   ='//CHAR(0))
!       PRINT *,' KORB: K1A   =',PARM(1,19),PARM(2,19),PARM(3,19)
        CALL TAUOLA_PRINT('KORB: K1B   ='//CHAR(0))
!       PRINT *,' KORB: K1B   =',PARM(1,20),PARM(2,20),PARM(3,20)
        CALL TAUOLA_PRINT('KORB: rho/rhop/rhopp -> 4pi mixing:'//CHAR(0))
!       PRINT *,' KORB: rho/rhop/rhopp -> 4pi mixing:'
        CALL TAUOLA_PRINT('KORB: rho   ='//CHAR(0))
!       PRINT *,' KORB: rho   =',PARM(1,9) ,PARM(2,9)
        CALL TAUOLA_PRINT('KORB: rhopp ='//CHAR(0))
!       PRINT *,' KORB: rhopp =',PARM(1,21),PARM(2,21),PARM(3,21)
        CALL TAUOLA_PRINT('KORB: rhoppp='//CHAR(0))
!       PRINT *,' KORB: rhoppp=',PARM(1,22),PARM(2,22),PARM(3,22)
C amplitudes for curr_cleo.F:
C for (3pi)-pi0: 4pi phase space; rho0pi-pi0; rho-pi+pi-; rho+pi-pi-; pi-omega
        PARM(3,31) = 0.
        PARM(3,32) = 0.1242
        PARM(3,33) = 0.1604
        PARM(3,34) = 0.2711
        PARM(3,35) = 0.4443
C for pi-3pi0: 4pi phase space; rho-pi0pi0
        PARM(3,36) = 0.
        PARM(3,37) = 1.0
C Modify amplitudes for 4pi form-factor in formf/curr, from korb.dec:
CCC        IF (IPLIST(2,282).EQ.5) THEN
        IPLIST=0
        IF (IPLIST.EQ.5) THEN
        PARM(3,31) = 0.0000
        PARM(3,32) = 0.1242
        PARM(3,33) = 0.1604
        PARM(3,34) = 0.2711
        PARM(3,35) = 0.4443
        PARM(3,36) = 0.0000
        PARM(3,37) = 1.0000
        END IF
        CALL TAUOLA_PRINT('KORB: 3PI-PI0 PARAMS:'//CHAR(0))
!       PRINT *,' KORB: 3PI-PI0 PARAMS:',(PARM(3,I),I=31,35)
        CALL TAUOLA_PRINT('KORB: PI-3PI0 PARAMS:'//CHAR(0))
!       PRINT *,' KORB: PI-3PI0 PARAMS:',(PARM(3,I),I=36,37)
C The 4pi models are the most complicated in TAUOLA.
C If the user has not modified any parameters of the 4pi model,
C we can use the WTMAX determined with many trials.
        IF (ABS(PARM(3,31)-0.0000).GT.0.0001 .OR.
     1      ABS(PARM(3,32)-0.1242).GT.0.0001 .OR.
     1      ABS(PARM(3,33)-0.1604).GT.0.0001 .OR.
     1      ABS(PARM(3,34)-0.2711).GT.0.0001 .OR.
     1      ABS(PARM(3,35)-0.4443).GT.0.0001 ) THEN
           PARM(3,38) = -1.
        ELSE
           PARM(3,38) = 6.9673671E-14
        END IF
        IF (ABS(PARM(3,36)-0.0000).GT.0.0001 .OR.
     1      ABS(PARM(3,37)-1.0000).GT.0.0001 ) THEN
           PARM(3,39) = -1.
        ELSE
           PARM(3,39) = 3.5374880E-13
        END IF
C phases for curr_cleo.F:
        PARM(3,42) = -0.40
        PARM(3,43) =  0.00
        PARM(3,44) = -0.20+3.1416
        PARM(3,45) = -1.50
C rho' contributions to rho' -> pi-omega:
        PARM(3,51) = -0.10
        PARM(3,52) =  1.00
        PARM(3,53) = -0.10
        PARM(3,54) = -0.04
C rho' contribtions to rho' -> rhopipi:
        PARM(3,55) =  1.00
        PARM(3,56) =  0.14
        PARM(3,57) = -0.05
        PARM(3,58) = -0.05
C rho contributions to rhopipi, rho -> 2pi:
        PARM(3,59) =  1.000
        PARM(3,60) = -0.145
        PARM(3,61) =  0.000
C Set the BRs for (A1+ -> rho+ pi0) and (K*+ -> K0 pi+)
C needed in dist/taurdf:
        PARM(4,1) = 0.4920                 ! BRA1+
        PARM(4,2) = 0.4920                 ! BRA1-
        PARM(4,3) = 0.6660                 ! BRKS+
        PARM(4,4) = 0.6660                 ! BRKS-
        PARM(4,5) = 0.5                    ! BRK0
        PARM(4,6) = 0.5                    ! BRK0B
C amplitude coefficients for tau -> K1(1270) / K1(1402)
        C1270 = PARM(3,19)
        C1402 = PARM(3,20)
        IF (C1270.EQ.0.AND.C1402.EQ.0.) THEN
           C1270 = 1.
           C1402 = 0.6
        END IF
C From PDG96, square roots of branching fractions:
        A1270_KSPI = SQRT(0.16)
        A1270_KRHO = SQRT(0.42)
        A1402_KSPI = SQRT(0.94)
        A1402_KRHO = SQRT(0.03)
C C-G coefficients for K1- -> CG1 * |K- pi0> + CG2 * |K0bar pi->
        CG1 = -SQRT(2./3.)
        CG2 =  SQRT(1./3.)
C and the resulting amplitudes (times normalized FF):
        PARM(3,81) = C1270*A1270_KSPI*CG1  ! K1270 -> K*0B pi-
        PARM(3,82) = C1402*A1402_KSPI*CG1  ! K1402 -> K*0B pi-
        PARM(3,83) = C1270*A1270_KRHO*CG1  ! K1270 -> K0B rho-
        PARM(3,84) = C1402*A1402_KRHO*CG1  ! K1402 -> K0B rho-
        PARM(3,85) = C1270*A1270_KSPI*CG2  ! K1270 -> K*- pi0
        PARM(3,86) = C1402*A1402_KSPI*CG2  ! K1402 -> K*- pi0
        PARM(3,87) = C1270*A1270_KRHO*CG2  ! K1270 -> K- rho0
        PARM(3,88) = C1402*A1402_KRHO*CG2  ! K1402 -> K- rho0
      END IF
**********************************************************************
      R = 0.
      IF (IF1.GE.1 .AND. IF1.LE.4 .AND. IF2.GE.1 .AND. IF2.LE.100) THEN
         R = PARM(IF1,IF2)
CAJW 4/4/94  Better to decide on A1 br now, avoid DADMAA/DPHSAA problem.
        IF (IF1.EQ.4.AND.JAK1.EQ.5) THEN
          IF (IF2.EQ.11) THEN
C Return the BR used in the last call:
             R = BRA1
          ELSE IF (IF2.EQ.1) THEN
            BRA1 = R
            CALL RANMAR(RRR,1)
            IF (RRR(1).LT.BRA1) THEN
              R = 1.     ! 3pi
            ELSE
              R = 0.     ! pi-2pi0
            END IF
            BRA1 = R
          END IF
        ELSEIF (IF1.EQ.4.AND.JAK1.EQ.7) THEN
          IF (IF2.EQ.13) THEN
C Return the BR used in the last call:
             R = BRKS
          ELSE IF (IF2.EQ.3) THEN
            BRKS = R
            CALL RANMAR(RRR,1)
            IF (RRR(1).LT.BRKS) THEN
              R = 1.     ! K0 pi-
            ELSE
              R = 0.     ! K- pi0
            END IF
            BRKS = R
          END IF
        END IF
      END IF
      PKORB = R
      RETURN
      END
      Complex*16 FUNCTION FPIBEL(W,flpar)
C****************************************************************************
C                       PION FORM FACTOR FROM BELLE
C
C            Rho(770)+rho(1450)+rho(1700) model is used as in:
C             M.Fujikawa et al., Phys.Rev.D 78 (2008) 072006
C                          -------------------
C       * flpar: = 0 - (all free) fit result for par(1...11)
C                = 1 - (par(1)=F_pi(0)=1-fixed) fit result for par(2...11)
C
C       * parameter : par(1)=overall norm. factor (=F_pi(0))
C                   : par(2)=rho(770) mass
C                   : par(3)=rho(770) width
C                   : par(4)=rho(1450) mass
C                   : par(5)=rho(1450) width
C                   : par(6)=rho(1450) admixture abs. value (|beta|)
C                   : par(7)=rho(1450) admixture phase (phi=arg(beta))
C                   : par(8)=rho(1700) mass
C                   : par(9)=rho(1700) width
C                   : par(10)=rho(1700) admixture abs. value (|gamma|)
C                   : par(11)=rho(1700) admixture phase (phi3=arg(gamma))
C
C       * x:=s=(M_pipi0)^2=W^2
C
C       * Notice: the following code was extracted (and checked) directly
C                 from the Belle fit function with minor cosmetic changes
C                 (this is internal comment of Belle)
C
C       * Notice: the following code was extracted (and checked) directly
C                 from  tauola_with_belle_fpi.F sent by D. Epifanov and
C                 it is an identical copy of
C                 Complex*16 FUNCTION FPIBEL(W,flpar)
C                 used in tauola_with_belle_fpi.F
C
C       * Called: CURR_PIPI0  if (FF2PIRHO = 2  or FF2PIRHO =3)
C                 in ../RChL-currents/value_parameter.f
C                 FF2PIRHO = 2 is for flpar = 0
C                 FF2PIRHO = 3 is for flpar = 1
C
C****************************************************************************
      implicit none
      integer flpar
      real W
      real*8 x,par(11),pi,pimas,pi0mas,taumas
      real*8 beta,berho1,berho2,berho3
      real*8 ps,prho1,prho2,prho3
      real*8 gamma1,gamma2,gamma3
      real*8 hs,h1,h2,h3,dhds1,dhds2,dhds3
      real*8 d1,d2,d3,fs1,fs2,fs3
      real*8 a1,b1,c1,bw_re1,bw_im1
      real*8 a2,b2,c2,bw_re2,bw_im2
      real*8 a3,b3,c3,bw_re3,bw_im3
      real*8 sinphi,cosphi,sinphi3,cosphi3
      complex*16 mbw1,mbw2,mbw3,mff,mephi,mephi3
      parameter(pi=3.141592653589)
      x=dble(W**2)
C------------- PARAMETERS --------------------
C     Particle mass (unit: GeV)
      pimas  = 0.1395702
      pi0mas = 0.13498
      taumas = 1.77699
C---------------------------------------------
C Opt. par. from Table VII of PRD78 (2008) 072006
      if(flpar.eq.0) then
         par(1) = 1.02
         par(2) = 0.7749
         par(3) = 0.1486
         par(4) = 1.428
         par(5) = 0.413
         par(6) = 0.13
         par(7) = 197
         par(8) = 1.694
         par(9) = 0.135
         par(10)= 0.028
         par(11)= -3
      else
         par(1) = 1.0
         par(2) = 0.7746
         par(3) = 0.1481
         par(4) = 1.446
         par(5) = 0.434
         par(6) = 0.15
         par(7) = 202
         par(8) = 1.728
         par(9) = 0.164
         par(10)= 0.037
         par(11)= 24
      endif
C========= Beta (s, rho(770), rho(1450), rho(1700)) ==================
       beta  =sqrt( (1.0 - (pimas - pi0mas)**2/x)
     +             *(1.0 - (pimas + pi0mas)**2/x) )
c      ----- rho(770) -----
       berho1=sqrt( ( 1.0 - (pimas - pi0mas)**2/(par(2)*par(2)) )
     +             *( 1.0 - (pimas + pi0mas)**2/(par(2)*par(2)) ) )
c      ----- rho(1450) -----
       berho2=sqrt( ( 1.0 - (pimas - pi0mas)**2/(par(4)*par(4)) )
     +             *( 1.0 - (pimas + pi0mas)**2/(par(4)*par(4)) ) )
c      ----- rho(1700) -----
       berho3=sqrt( ( 1.0 - (pimas - pi0mas)**2/(par(8)*par(8)) )
     +             *( 1.0 - (pimas + pi0mas)**2/(par(8)*par(8)) ) )
C========= Momentum (s, rho(770), rho(1450), rho(1700)) ==============
C
C       ps    : momentum of pi at s
C       prho1 : momentum of pi at s=rho(770)**2
C       prho2 : momentum of pi at s=rho(1450)**2
C       prho3 : momentum of pi at s=rho(1700)**2
       ps  = 0.5 * sqrt(x) * beta
c      ----- rho(770) -----
       prho1 = 0.5* par(2) * berho1
c      ----- rho(1450) -----
       prho2 = 0.5* par(4) * berho2
c      ----- rho(1700) -----
       prho3 = 0.5* par(8) * berho3
C========= Width (rho(770), rho(1450), rho(1700)) ====================
C      ----- rho(770) -----
       gamma1 = par(3) * ( par(2)*par(2)/x ) * (ps/prho1)**3
C      ----- rho(1450) -----
       gamma2 = par(5) * ( par(4)*par(4)/x ) * (ps/prho2)**3
C      ----- rho(1700) -----
       gamma3 = par(9) * ( par(8)*par(8)/x ) * (ps/prho3)**3
C============== h(rho(770), rho(1450), rho(1700)) ====================
C
C       hs : h(s)...s(s=x)
C       h1 : h(s) for s=rho(770)**2
C       h2 : h(s) for s=rho(1450)**2
C       h3 : h(s) for s=rho(1700)**2
        hs = (2.0/pi)*(ps/sqrt(x))*
     +       log( (sqrt(x) + 2.0*ps)/(2.0*pimas) )
C      ----- rho(770) -----
        h1 = (2.0/pi)*(prho1/par(2))*
     +       log( (par(2) + 2.0*prho1)/(2.0*pimas) )
C      ----- rho(1450) -----
        h2 = (2.0/pi)*(prho2/par(4))*
     +       log( (par(4) + 2.0*prho2)/(2.0*pimas) )
C      ----- rho(1700) -----
        h3 = (2.0/pi)*(prho3/par(8))*
     +       log( (par(8) + 2.0*prho3)/(2.0*pimas) )
C============ dhds(rho(770), rho(1450), rho(1700)) ===================
C
C       dhds = dh/ds = h(M_rho)[1/(8*p(M_rho)^2) - 1/(2M_rho^2)
C                                              + 1/(2pi*M_rho^2)
C       dhds1 : dh/ds for M_rho=rho(770)
C       dhds2 : dh/ds for M_rho=rho(1450)
C       dhds3 : dh/ds for M_rho=rho(1700)
C      ----- rho(770) -----
        dhds1 = h1*( 1.0/(8.0*prho1**2) - 1.0/(2.0*par(2)**2) )
     +                       +  1.0/(2.0*pi*par(2)**2)
C      ----- rho(1450) -----
        dhds2 = h2*( 1.0/(8.0*prho2**2) - 1.0/(2.0*par(4)**2) )
     +                       +  1.0/(2.0*pi*par(4)**2)
C      ----- rho(1700) -----
        dhds3 = h3*( 1.0/(8.0*prho3**2) - 1.0/(2.0*par(8)**2) )
     +                       +  1.0/(2.0*pi*par(8)**2)
C============ d(rho(770), rho(1450), rho(1700)) ======================
C
C       d = 3/pi * (m_pi^2/p(M_rho)^2) *
C            ln {(M_rho+2*p(M_rho))/2M_pi} + M_rho/(2pi*p(M_rho))
C                               - (m_pi^2 * M_rho)/(pi*p(M_rho)^3)
C       d1 : d for M_rho=rho(770)
C       d2 : d for M_rho=rho(1450)
C       d2 : d for M_rho=rho(1700)
C      ----- rho(770) -----
        d1 = (3.0/pi)*(pimas**2/prho1**2)
     +         * log((par(2) + 2.0*prho1)/(2.0*pimas))
     +         + par(2)/(2.0*pi*prho1)
     +         - ((pimas**2)*par(2))/(pi*prho1**3)
C      ----- rho(1450) -----
        d2 = (3.0/pi)*(pimas**2/prho2**2)
     +         * log((par(4) + 2.0*prho2)/(2.0*pimas))
     +         + par(4)/(2.0*pi*prho2)
     +         - ((pimas**2)*par(4))/(pi*prho2**3)
C      ----- rho(1700) -----
        d3 = (3.0/pi)*(pimas**2/prho3**2)
     +         * log((par(8) + 2.0*prho3)/(2.0*pimas))
     +         + par(8)/(2.0*pi*prho3)
     +         - ((pimas**2)*par(8))/(pi*prho3**3)
C================ f(s) (rho(770), rho(1450), rho(1700) ) ==========
C
C       f(s) = gamma-rho * M_rhp**2/p(M_rho)[ p(s)^2(h(s)-h(M_rho)
C                                 + (M_rho^2 -s)p(M_rho)^2 * dh/ds]
C       fs1 : f(s) for M_rho=rho(770)
C       fs2 : f(s) for M_rho=rho(1450)
C       fs3 : f(s) for M_rho=rho(1700)
C
C      ----- rho(770) -----
        fs1 = par(3) * (par(2)**2/prho1**3) *
     +        ( ps**2 *(hs - h1) + (par(2)**2 -x)*prho1**2 * dhds1 )
C      ----- rho(1450) -----
        fs2 = par(5) * (par(4)**2/prho2**3) *
     +        ( ps**2 *(hs - h2) + (par(4)**2 -x)*prho2**2 * dhds2 )
C      ----- rho(1700) -----
        fs3 = par(9) * (par(8)**2/prho3**3) *
     +        ( ps**2 *(hs - h3) + (par(8)**2 -x)*prho3**2 * dhds3 )
C====== BW form G&S model (rho(770) + rho(1450) + rho(1700)) ======
C    bw_re - real part ; bw_im - imaginary part
C      ----- rho(770) -----
       a1=par(2)*par(2) - x + fs1
       b1=sqrt(x)*gamma1
       c1=par(2)*par(2) + d1*par(2)*par(3)
       bw_re1= (a1*c1)/(a1**2 + b1**2)
       bw_im1= (b1*c1)/(a1**2 + b1**2)
       mbw1=dcmplx(bw_re1,bw_im1)
c      ----- rho(1450) -----
       a2=par(4)*par(4) - x + fs2
       b2=sqrt(x)*gamma2
       c2=par(4)*par(4)+ d2*par(4)*par(5)
       bw_re2= (a2*c2)/(a2**2 + b2**2)
       bw_im2= (b2*c2)/(a2**2 + b2**2)
       mbw2=dcmplx(bw_re2,bw_im2)
c      ----- rho(1700) -----
       a3=par(8)*par(8) - x + fs3
       b3=sqrt(x)*gamma3
       c3=par(8)*par(8)+ d3*par(8)*par(9)
       bw_re3= (a3*c3)/(a3**2 + b3**2)
       bw_im3= (b3*c3)/(a3**2 + b3**2)
       mbw3=dcmplx(bw_re3,bw_im3)
C====== admixtures of rho(1450) and rho(1700) ======
C --------- exp(i*arg(beta)) ----------
       sinphi=sin(par(7)*(pi/180.0))  ! use degree
       cosphi=cos(par(7)*(pi/180.0))  ! use degree
       mephi=dcmplx(cosphi,sinphi)
C --------- exp(i*arg(gamma)) ---------
       sinphi3=sin(par(11)*(pi/180.0))  ! use degree
       cosphi3=cos(par(11)*(pi/180.0))  ! use degree
       mephi3=dcmplx(cosphi3,sinphi3)
C================= Form factor =====================
       mff = 1.0/(1.0+par(6)*mephi+par(10)*mephi3)
       FPIBEL=par(1)*mff*(mbw1+par(6)*mephi*mbw2+par(10)*mephi3*mbw3)
      RETURN
      END
      COMPLEX FUNCTION F3PI_RCHT(IFORM,QQ,SA,SB)
      IMPLICIT NONE
      INTEGER                    IFORM
      REAL                             QQ,SA,SB
C.......................................................................
C.
C.    F3PI - RchT version of the hadronic curent used in TAUOLA
C.    References: [1] arXiv:0911.4436 (hep-ph) P. Roig et al.
C.                eqs (5)-(9), (32) gives the main part of the current
C                 (the part without the sigma contribution)
C.                [2] arXiv:1203.3955,
C.                the manual for the 3pi current without the sigma contribution
C.                [3] http://annapurna.ifj.edu.pl/~wasm/RChL/RChLa1.pdf
C                     eq (3) the sigma meson contribution to the 3 pi current
C.    Inputs    : QQ,SA,SB - invariant masses**2  [GeV**2]
C.              : IFORM formfactor no.
C.    Outputs   : F3PI_RCHT formfactor value
C.
C.    COMMON    : RCHT_3PI content is defined in this routine
C.
C.    Calls     : functions from file ./funct_rpt.f,  ./fa1rchl.f
C.    Called    : from file ../formf.f
C.    Author    : O.S
C.    Created   :
C.    Modified  : 1. Feb 2011
C                 a part with scalars added on the 1st May 2012
C.......................................................................
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      INTEGER          J3PI
      REAL wid_a1_fit
*
      CHARACTER*(*) CRNAME
      PARAMETER(    CRNAME = 'F3PI' )
*
      INTEGER IDK
      DOUBLE PRECISION M1,M2,M3,M1SQ,M2SQ,M3SQ
      REAL GGMA1
      REAL S1,S2,S3,FACT_RPT
*
      REAL R_RPT,KAP_RPT,FACT_ADD_RPT
      COMPLEX ALP21_RPT,ALP22_RPT,ALP11_RPT,ALP12_RPT
     $       ,BETA1_RPT,BETA2_RPT,BETA1_RPT_RHO1,ALP11_RPT_RHO1
     $       ,ALP21_RPT_RHO1,ALP22_RPT_RHO1
      INTEGER IFIRST
      DATA IFIRST/0/
C. GENERAL INITIALIZATION
C. ======================
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_taukle') /TAUKLE/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
        CALL TAUOLA_PRINT(
     $    'In F3pi: chpt + 1 resonance + 2 resonance (RchT)'//CHAR(0))
!       PRINT *,' In F3pi: chpt + 1 resonance + 2 resonance (RchT)'
      END IF
C******************************************
C    Initilisation of the mass of the particles
C*****************************************
        call rchl_parameters(5)
C. We impose isospin symmetry requesting that charged and neutral pion mass
C. are equal. This may need to be changed MMPI_AV = (2.*MPIC+MPIZ)/3.
        M1 = MMPI_AV
        M2 = MMPI_AV
        M3 = MMPI_AV
C. normalization factor to compensate on different
C. convention for normalization constant  used in  TAUOLA and
C. TAUOLA documentation on one side and paper [1] on other.
        FACT_ADD_RPT = 1./FPI_RPT
C. FUNCTION VALUE, GENERAL CASE
C. ============================
C     Function value set to 0
      F3PI_RCHT = CMPLX(0.,0.)
C.    First determine whether we are doing pi-2pi0 or 3pi.
C. we  CH3PI GET information (eg from phase space geneator of tauola.f)
C. whether it is 3 prong J3pi=1 or 1 prong J3pi=2 final state of 3 pion.
      CALL CH3PIGET(J3pi)
      IF (J3pi.EQ.2) THEN  ! it is pi 2pi0
          IDK     = 1
          R_RPT   = -1.
          KAP_RPT = 0.5
      ELSE IF(J3pi.EQ.1) THEN              !  it is 3pi
          IDK     = 2
          R_RPT   = 1.
          KAP_RPT = 1.
      END IF
      M1SQ = M1*M1
      M2SQ = M2*M2
      M3SQ = M3*M3
C. Calculation, IFORM = 1 or 2.
C.   VECTOR 3 PION FORM FACTORS
C. ============================
      IF (IFORM.EQ.1.OR.IFORM.EQ.2) THEN
        S1 = SA     ! t variable in  [2] (!!! vec1 = v2 !!!)
        S2 = SB    ! s variable in [2]
        S3 = QQ-SA-SB+M1SQ+M2SQ+M3SQ
        IF (S3.LE.0..OR.S2.LE.0.) RETURN
C.    FUNCTIONS BETA_RPT, ALP1_RPT are coded in ./funct_rpt.f
C.    they are defined in Eq (6) of [2]
        BETA1_RPT = BETA_RPT(QQ,S1,S2,M1SQ,M2SQ,M3SQ,MRO,MRHO1,GRHO1)
        ALP11_RPT = ALP1_RPT(QQ,S1,S2,M1SQ,M2SQ,M3SQ,MRO,MRHO1,GRHO1)
        F3pi_RCHT = - 2.*SQRT(2.)/(3.*FPI_RPT) -
     $     SQRT(2.)*FV_RPT*GV_RPT/(3.*FPI_RPT**3)*ALP11_RPT +
     $     4.*FA_RPT*GV_RPT/(3.*FPI_RPT**3)*BETA1_RPT*QQ
     $     *FA1RCHL(QQ)
C.    FA1RCHL(QQ) is the a1 propagator, is it coded in ./fa1rchl.f
        F3pi_RCHT =  F3pi_RCHT * R_RPT
C.
C.      The contribution from the scalar (sigma) resonance
C.
      IF(FF3PISCAL.EQ.1) THEN
C.      This parametrization does not fit
         IF(J3PI.EQ.1) THEN
        F3pi_RCHT = F3pi_RCHT + (
     &            sqrt(2.)*(R0scal_3pi(QQ,S1)+R0scal_3pi(QQ,S2)) +
     &            (R2scal_3pi(QQ,S1)+R2scal_3pi(QQ,S2))
     &            )
         ELSE IF(J3PI.EQ.2) THEN
        F3pi_RCHT = F3pi_RCHT - (
     &            (R0scal_3pi(QQ,S1)+R0scal_3pi(QQ,S2)) -
     &            sqrt(2.)*(R2scal_3pi(QQ,S1)+R2scal_3pi(QQ,S2))
     &            )
         ENDIF
       ELSE IF(FF3PISCAL.EQ.2) THEN
C.    A new parametrization 10.02.2013 for the scalar contribution
C.       analytical formulae in [3], eqs (3)-(6),
C.       functions BWsig(Mm,Gg,Qx), FFsig(QQ,Qx) coded in ./funct_rpt.f
C.
         IF(J3PI.EQ.1) THEN
        F3pi_RCHT = F3pi_RCHT +
     &              SQRT(2.)*FV_RPT*GV_RPT/(3.*FPI_RPT**3)*
     &             (alpsig*BWsig(MSIG,GSIG,S1)*FFsig(QQ,S1) +
     &              betasig*BWsig(MSIG,GSIG,S2)*FFsig(QQ,S2))
     &             +4.*FA_RPT*GV_RPT/(3.*FPI_RPT**3)*
     &             (gamsig*BWsig(MSIG,GSIG,S1)*FFsig(QQ,S1) +
     &              delsig*BWsig(MSIG,GSIG,S2)*FFsig(QQ,S2))
     $     *FA1RCHL(QQ)
        ELSE IF(J3PI.EQ.2) THEN
        F3pi_RCHT = F3pi_RCHT - (
     &              SQRT(2.)*FV_RPT*GV_RPT/(3.*FPI_RPT**3)*
     &             alpsig*BWsig(MSIG,GSIG,S3)*FFsig(QQ,S3)
     &             +4.*FA_RPT*GV_RPT/(3.*FPI_RPT**3)*
     &             gamsig*BWsig(MSIG,GSIG,S3)*FFsig(QQ,S3)
     $     *FA1RCHL(QQ) )
         ENDIF
C.
C.     The Coulomb interaction effects for the final pions
C.          only for pi-pi+pi- should be included
C.     Functions fattcoul(m1,m2,s),frepcoul(m1,m2,s) are coded in ./funct_rpt.f
C.
      IF (FCOUL.EQ.1.AND.(J3pi.EQ.1)) THEN !  it is 3pi
        F3PI_RCHT = F3PI_RCHT*sqrt(fattcoul(m2,m3,s1)
     &      *fattcoul(m1,m3,s2)*frepcoul(m1,m2,s3))
      END IF
       ENDIF
C.
C. The factor 1/FACT_ADD_RPT = FPI_RPT comes to compensate an additional factor
C.  in the hadronic current F3pi_rcht  above compare with eq (6) in [2]
       F3pi_RCHT = F3pi_RCHT/FACT_ADD_RPT
C. Calculation, for  IFORM = 3 is not needed
C. =======================
C.   F3PI_RCHT is set to zero in ../formf.f.
C. Calculation, IFORM = 4.
C. PSEUDOSCALAR 3 PION FORM FACTOR
C. =======================
      ELSEIF (IFORM.EQ.4) THEN
        S1 = SA
        S2 = SB
        S3 = QQ-SA-SB+M1SQ+M2SQ+M3SQ
        IF (S3.LE.0..OR.S2.LE.0.) RETURN
C. Functions  ALP21_RPT,  ALP22_RPT are Eq(10) in [2]
C. GRHO_RCHT, GRHO1_RCHT s1 or s2 dependent widths of rho or rho1
C. coded in ./funct_rpt.f
        ALP21_RPT = 3.*GV_RPT*S1*M1SQ*(S3-S2)/
     &    (FV_RPT*QQ*(QQ-M1SQ)*(1.+BETA_RHO))*
     &    (1./(S1-MRO**2+i*MRO*GRHO_RCHT(S1,MRO))+
     &     BETA_RHO/(S1-MRHO1**2+i*MRHO1*GRHO1_RCHT(S1,MRHO1,GRHO1)))
        ALP22_RPT = 3.*GV_RPT*S2*M1SQ*(S3-S1)/
     &    (FV_RPT*QQ*(QQ-M1SQ)*(1.+BETA_RHO))*
     &    (1./(S2-MRO**2+i*MRO*GRHO_RCHT(S2,MRO))+
     &     BETA_RHO/(S2-MRHO1**2+i*MRHO1*GRHO1_RCHT(S2,MRHO1,GRHO1)))
C.  PSEUDOSCALAR 3 PION FORM FACTOR. Eqs (9) in [2]
        F3PI_RCHT = SQRT(2.)/(3.*FPI_RPT*QQ*(QQ-M1SQ))*
     &      M1SQ*(3.*(S3-M3SQ)-QQ*(1.+2.*KAP_RPT*R_RPT))
        F3PI_RCHT = F3PI_RCHT -  sqrt(2.)*FV_RPT*GV_RPT/(3.*FPI_RPT**3)*
     &    ( ALP21_RPT + ALP22_RPT)
C.
C. The factor 1/FACT_ADD_RPT = FPI_RPT comes to compensate an additional factor
C.  in the hadronic current F3pi_rcht  above compare with eq (9) in [2]
        F3PI_RCHT = R_RPT * F3PI_RCHT/FACT_ADD_RPT
C.
C.     The Coulomb interaction effects for the final pions
C.          only for pi-pi+pi- should be included
C.     Functions fattcoul(m1,m2,s),frepcoul(m1,m2,s) are coded in ./funct_rpt.f
C.
      IF (FCOUL.EQ.1.AND.(J3pi.EQ.1)) THEN !  it is 3pi
        F3PI_RCHT = F3PI_RCHT*sqrt(fattcoul(m2,m3,s1)
     &      *fattcoul(m1,m3,s2)*frepcoul(m1,m2,s3))
      END IF
      END IF
      RETURN
      END
      FUNCTION ALP1_RPT(Q,S,T,M1SQ,M2SQ,M3SQ,
     $                           MMM,MMM2,GGG)
      IMPLICIT NONE
      COMPLEX ALP1_RPT
      REAL                      Q,S,T
      DOUBLE PRECISION                M1SQ,M2SQ,M3SQ,
     &                           MMM,MMM2,GGG
C **********************************************************
C     COMPLEX FUNCTION ALP1_RPT;
C     one-resonance contribution to
C     F3pi_rcht(iform=1,2) at f3pi_rcht.f;
C     corrected  formula (8) of REF [1] including rho1,
C     the factor sqrt(2)*FV_RPT*GV_RPT/(3.*FPI_RPT**3) is included in
C     F3pi_rcht(iform1=1,2) at f3pi_rcht.f
C **********************************************************
      REAL U
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
       U= Q-S-T+M1SQ+M2SQ+M3SQ
       ALP1_RPT = - 3.*S/(1.+BETA_RHO)*
     $     (1./(S-MMM**2+i*MMM*GRHO_RCHT(S,MMM))+
     $      BETA_RHO/(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG)))
     $      +(2.*GV_RPT/FV_RPT-1.)*
     $      ((2.*Q-2.*S-U)/(1.+BETA_RHO)*
     $      (1./(S-MMM**2+i*MMM*GRHO_RCHT(S,MMM))+
     $      BETA_RHO/(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG)))+
     $      (U-S)/(1.+BETA_RHO)*
     $      (1./(T-MMM**2+i*MMM*GRHO_RCHT(T,MMM))+
     $       BETA_RHO/(T-MMM2**2+i*MMM2*GRHO1_RCHT(T,MMM2,GGG))))
      RETURN
      END
      FUNCTION BETA_RPT(Q,S,T,M1SQ,M2SQ,M3SQ,MMM,MMM2,GGG)
      IMPLICIT NONE
      COMPLEX BETA_RPT
      REAL                      Q,S,T
      DOUBLE PRECISION                M1SQ,M2SQ,M3SQ,MMM,MMM2,GGG
C **********************************************************
C     COMPLEX FUNCTION BETA_RPT ; two-resonance contribution to
C     F3pi_rcht(iform=1,2) at f3pi_rcht.f;
C     corrected  formula (8) of REF [1] including rho1,
C     the factor 4*FA_RPT*GV_RPT/(3.*FPI_RPT**3)*QQ/D_a1(QQ)
C     is included in F3pi_rcht(iform1=1,2) at f3pi_rcht.f
C **********************************************************
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      REAL U,LAM0_RPT,LAM1_RPT,LAM2_RPT,FF1_RPT,FF2_RPT,FF_REL
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      FF_REL = FPI_RPT*FPI_RPT/(FV_RPT*FV_RPT)
       LAM1_RPT = FPI_RPT*FPI_RPT/(2.D0*SQRT(2.D0)*FA_RPT*GV_RPT)
       LAM2_RPT = -(1.-2.*FF_REL)*LAM1_RPT
       LAM0_RPT = (LAM1_RPT + LAM2_RPT)/4.
        U= Q-S-T+M1SQ+M2SQ+M3SQ
      FF1_RPT = -LAM0_RPT*M1SQ/Q +LAM1_RPT*S/Q+LAM2_RPT
      FF2_RPT = -LAM0_RPT*M1SQ/Q +LAM1_RPT*T/Q+LAM2_RPT
      BETA_RPT = -3.*(LAM1_RPT+LAM2_RPT)*S/(1.+BETA_RHO)*
     $     (1./(S-MMM**2+i*MMM*GRHO_RCHT(S,MMM))+
     $      BETA_RHO/(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG)))+
     $         FF1_RPT*(2.*Q+S-U)/(1.+BETA_RHO)*
     $     (1./(S-MMM**2+i*MMM*GRHO_RCHT(S,MMM))+
     $      BETA_RHO/(S-MMM2**2+i*MMM2*GRHO1_RCHT(S,MMM2,GGG)))+
     $         FF2_RPT*(U-S)/(1.+BETA_RHO)*
     $     (1./(T-MMM**2+i*MMM*GRHO_RCHT(T,MMM))+
     $      BETA_RHO/(T-MMM2**2+i*MMM2*GRHO1_RCHT(T,MMM2,GGG)))
      RETURN
      END
      COMPLEX FUNCTION FA1RCHL(QQ)
      IMPLICIT NONE
      REAL                     QQ
      DOUBLE PRECISION M1,M2,M3
      REAL GGMA1
C.......................................................................
C.
C.    FA1CHL - RchT version of  the A1 propagator
C.
C.    Inputs    : QQ - invariant masses**2  [GeV**2]
C.    Outputs   : FA1RCHL formfactor value at QQ
C.
C.    COMMON    : RCHT_3PI content is defined in this routine
C.
C.    Calls     : functions from file ./wid_a1_fit.f
C.    Called    : from file f3pi_rcht.f, fkkpi.f, fkk0pi0.f
C************************************************************************
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
C******************************************
C    Initilisation of the mass of the particles
C*****************************************
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
        call rchl_parameters(5)
c$$$C we impose isospin symmetry requesting that charged and neutral pion mass
c$$$C are equal. This may need to be changed
c$$$        MMPI_AV = (2.*MPIC+MPIZ)/3.
        M1 = MMPI_AV
        M2 = MMPI_AV
        M3 = MMPI_AV
c
C   Function wid_a1_fit.f calculates the energy dependence of
C   the a1 meson width
C
      IF(QQ.GE.(M1+M2+M3)**2) THEN
        GGMA1 = wid_a1_fit(QQ)
      ELSE
        GGMA1 = 0.
      ENDIF
      FA1RCHL = 1./(QQ-MMA1*MMA1+i*MMA1*GGMA1)
      RETURN
      END
      subroutine getFF3PISCAL(INUM)
      INTEGER INUM
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      INUM=FF3PISCAL
      return
      end
      subroutine setFF3PISCAL(INUM)
      INTEGER INUM
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      FF3PISCAL=INUM
      return
      end
      DOUBLE PRECISION FUNCTION FFWID3PI(QQ,S1,S3)
      IMPLICIT NONE
      DOUBLE PRECISION                   QQ,S1,S3
C **************************************************************
C     Input:  QQ S1 S3   ! mpi-pi-pi+**2   mpi-pi+**2  mpi-pi-**2
C     Calls: functions FORM1, FORM2, FORM4
C     Uses constants: tau mass, pi mass, normalization constant.
C     Load: Initialized tauola library,
C     Output:  d\Gamma(tau --> 3pi nu)/(dQQ dS1 dS3)
C     Remark: If QQ S1 S3 are outside of the phase space
C             function FFWID3PI returns zero.
C **************************************************************
      COMPLEX F1,F2,F4
      DOUBLE PRECISION V11,V12,V22,GGF2,VUD2,ABS1,QQMIN,
     &                 QQMAX,S3MAX,S3MIN,S1MIN,S1MAX
      REAL XQQ,XS1,XS3, XS2,RQQ
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      DOUBLE PRECISION  XLAM,X,Y,Z
      DOUBLE PRECISION  XAMPI2
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      INTEGER  IMODE,IDUM,IFRCHL
      REAL RRQ
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      XLAM(X,Y,Z)= sqrt(abs((x-y-z)**2 - 4.*y*z))
      ABS1 = 1.d-5
      GGF2 = GFERMI**2
      VUD2 = CCABIB**2
C     TO CHANGE THE VARIABLES TO SINGLE PRECISION
C     INPUT FOR FORM1,FORM2,FORM4 IS SINGLE PRECISION
      XQQ = QQ
      XS1 = S1
      XS2 = QQ -S1-S3 + 3.*AMPI**2
      XS3 = S3
      XAMPI2 = AMPI**2
C     Limits for PHASE SPACE
C     Limits for QQ
      QQMIN = 9.D0*AMPI**2
      QQMAX  = (AMTAU-AMNUTA)**2
C     Limits for S1
      S1MAX=(DSQRT(QQ) - AMPI)**2 -ABS1
      S1MIN=4.D0*AMPI**2 +ABS1
C    LIMIT FOR XS3
      S3MAX = (QQ - AMPI**2)**2 -
     &    ( XLAM(QQ,S1,XAMPI2)
     &         - XLAM(S1,XAMPI2,XAMPI2) )**2
      S3MIN = (QQ - AMPI**2)**2 -
     &    (XLAM(QQ,S1,XAMPI2)
     &          + XLAM(S1,XAMPI2,XAMPI2) )**2
      S3MAX = S3MAX/4./S1
      S3MIN = S3MIN/4./S1
C     Check on PHASE SPACE
C
      IF((XS2.LE.0.) .OR.(S3MAX.LE.S3MIN)
     &   .OR.(XS1.LE.S1MIN).OR.(XS1.GE.S1MAX)
     &   .OR.(XS3.LE.S3MIN).OR.(XS3.GE.S3MAX)
     &   .OR.(QQ.LE.QQMIN).OR.(QQ.GE.QQMAX)
     &    )  THEN
        FFWID3PI = 0.D0
        RETURN
      ENDIF
C
      V11 = -XS1+4.D0*AMPI**2 -(XS2-XS3)**2/(4.D0*XQQ)
      V22 = -XS2+4.D0*AMPI**2 - (XS3-XS1)**2/(4.D0*XQQ)
      V12 = 0.5D0*(XS3-XS1-XS2+4.D0*AMPI**2)-0.25D0*(XS3-XS2)*(XS3-XS1)/XQQ
      F1 = FORM1(0,XQQ,XS1,XS2)
      F2 = FORM2(0,XQQ,XS2,XS1)
      F4 = FORM4(0,XQQ,XS2,XS1,XS3)
C formula 3.46 of [3]
      FFWID3PI = ABS(F1*CONJG(F1))*V11+ABS(F2*CONJG(F2))*V22+
     $       2.D0*REAL(F1*CONJG(F2))*V12
      CALL IFGFACT(2,IMODE,IDUM)
      CALL INIRChLget(IFRCHL)
      IF (IMODE.EQ.0) THEN
C      VERSION A: The 3 pion contribution to the a1 width
C      factor of a1 phase space and zeroing a1 propagator etc.is
C      done in RCHLWIDA1PI
       RQQ=QQ
       IF (IFRCHL.EQ.1) THEN
C       CASE OF RCHL
        FFWID3PI   = RCHLWIDA1PI(RQQ,FFWID3PI)
       ELSE
C       NOT READY YET
        CALL TAUOLA_PRINT(
     $    'FFWID3PI is not ready for non rchl currents'//CHAR(0))
!       WRITE(*,*) 'FFWID3PI is not ready for non rchl currents'
        CALL TAUOLA_STOP()
       ENDIF
C      to get the total a1 width the contribution from (KKPI)- and K-K0pi0
C      channels have to be added
      ELSE
C      VERSION B: calculation of 3 pion spectra in tau to 3pi nu channel.
C      factor for phase space of tau to XQQ nu decay and contribution from F4
C      (formula 3.21 of [3])
      FFWID3PI = (- FFWID3PI/3.D0*(1.D0+2.D0*XQQ/(AMTAU**2))+
     $       XQQ*ABS(F4*CONJG(F4)))*(AMTAU**2/XQQ-1.D0)**2
C Flux factor and normalization const.
        FFWID3PI =
     &       GGF2*VUD2/(128.D0*(2.D0*PI)**5*AMTAU)/2.d0
     &            *FFWID3PI
       IF (IFRCHL.EQ.1) THEN
C       CASE OF RCHL
C RChL normalization constant
        FFWID3PI =FFWID3PI/GETFPIRPT(1)**2
       ELSE
C       NOT READY YET
        CALL TAUOLA_PRINT(
     $    'FFWID3PI is not ready for non rchl currents'//CHAR(0))
!       WRITE(*,*) 'FFWID3PI is not ready for non rchl currents'
        CALL TAUOLA_STOP()
       ENDIF
      ENDIF
      RETURN
      END
      REAL FUNCTION RCHLWIDA1PI(RQQ,FFWID3PI)
C The 3 pion contribution to the a1 width
C (in [3] simple pretabulation is used through formula 3.48)
C for calculation of g(QQ) of 3.45 3.46 of [3] in RChL style
C a1 propagator has to be taken with the zero width.
      IMPLICIT NONE
      REAL RQQ
      DOUBLE PRECISION FFWID3PI
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      DOUBLE PRECISION  XLAM,X,Y,Z
      DOUBLE PRECISION  XAMPI2
      DOUBLE PRECISION  GETFPIRPT
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
C     AMA1 should be replaced by variable from the rchl namespace.
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
!
      RCHLWIDA1PI=- 1.0/REAL(FA1RCHL(RQQ)*CONJG(FA1RCHL(RQQ)))/RQQ**2
     $            /(96.D0*8.D0*PI**3*AMA1)/(FA_RPT**2*FPI_RPT**2)
     $            *FFWID3PI/2.d0
      END
      DOUBLE PRECISION FUNCTION  GETFPIRPT(I)
      IMPLICIT NONE
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      INTEGER I
!
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
!
      GETFPIRPT=FPI_RPT
      END
      FUNCTION GRHO_RCHT(XS,XMMM)
      IMPLICIT NONE
      REAL GRHO_RCHT
      REAL               XS
      DOUBLE PRECISION      XMMM
C **********************************************************
C     REAL FUNCTION Gamma Rho ; energy-dependent width of  rho meson
c     in SU(2) limit mpi=mpi0, mk=mk0;
C     formula (14) of REF [1]
C **********************************************************
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
       REAL MMPI_AV2,MMK_2
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      MMPI_AV2 = MMPI_AV**2
      MMK_2 = MMK**2
      IF(XS.GE.(4.*MMK_2)) THEN
           GRHO_RCHT=XMMM*XS*((1.-4.*MMPI_AV2/XS)**1.5
     $              +0.5*(1.-4.*MMK_2/XS)**1.5)
     $                  /(96.*PI*FPI_RPT**2)
      ELSE IF((XS.GE.(4.*MMPI_AV2)).AND.(XS.LE.(4.*MMK_2))) THEN
        GRHO_RCHT=XMMM*XS*(1.-4.*MMPI_AV2/XS)**1.5
     $            /(96.*PI*FPI_RPT**2)
      ELSE
        GRHO_RCHT = 0.
      ENDIF
      RETURN
      END
      FUNCTION GRHO1_RCHT(XS,XMMM,XGGG)
      IMPLICIT NONE
      REAL GRHO1_RCHT
      REAL                XS
      DOUBLE PRECISION       XMMM,XGGG
C **********************************************************
C     REAL FUNCTION Gamma Rho1 ; energy-dependent width of  rho1 meson;
c     in SU(2) limit mpi=mpi0; only rho' -> 2pi loop;
C     formula (33) of REF [1]
C **********************************************************
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
       REAL MMPI_AV2
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      MMPI_AV2 = MMPI_AV**2
      IF (XS.GE.(4.*MMPI_AV2)) THEN
        GRHO1_RCHT=XGGG*SQRT(XMMM**2/XS)*
     &      ((XS-4.*MMPI_AV2)/(XMMM**2-4.*MMPI_AV2))**1.5
      ELSE
        GRHO1_RCHT = 0.
      ENDIF
      RETURN
      END
      FUNCTION SIGP(SS)
      IMPLICIT NONE
      DOUBLE PRECISION SIGP
      DOUBLE PRECISION               SS
C***********************************************************
C     DOUBLRE PRECISION FUNCTION
c     of two body phase space threshold: equal mass scalars
C***********************************************************
      REAL TT
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      TT = 1. - 4.*MMPI_AV**2/SS
      IF (TT.GE.0) THEN
      SIGP = SQRT(TT)
      ELSE
      SIGP = 0.
      ENDIF
      RETURN
      END
      FUNCTION LAMB_RCHT(X1,X2,X3)
      IMPLICIT NONE
      REAL LAMB_RCHT
      REAL               X1,X2,X3,ARG_RCHT
C***********************************************************
C     REAL FUNCTION LAMDA of three body phase space
C***********************************************************
      ARG_RCHT = (X1-X2-X3)**2 - 4.*X3*X2
      IF(ARG_RCHT.GE.0.) THEN
      LAMB_RCHT = ARG_RCHT
      ELSE
      LAMB_RCHT = 0.
      ENDIF
      RETURN
      END
      FUNCTION R0SCAL_3PI(QX,SX)
C****************************************************************
C     Complex Function: R_0 function (Pablo notes)
C     for the scalar contribution for three pion modes
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      COMPLEX R0SCAL_3PI
      REAL            QX,SX,delta0_3piscal,xsx,xst
      DOUBLE PRECISION    XM1,XM2,DSX
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
c$$$      a00_3piscal = 0.220
c$$$      b00_3piscal = 0.268/mmpi_av**2
c$$$      c00_3piscal = -0.0139/mmpi_av**4
c$$$      d00_3piscal = -0.00139/mmpi_av**6
c$$$      x00_3piscal = 36.77*mmpi_av**2
c$$$c      MMF0 = 0.98
c$$$
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      dsx = sx
      xsx = sx/4.*sigp(dsx)**2
      xst = sqrt(sx)
      if(sx.le.0.7) then
      delta0_3piscal = sigp(dsx)*(a00_3piscal + b00_3piscal*xsx +
     &               c00_3piscal*xsx**2 + d00_3piscal*xsx**3)
      delta0_3piscal = delta0_3piscal*
     &         (4*mmpi_av**2 - x00_3piscal)/(sx-x00_3piscal)
      else if (xst.le.1.21) then
      delta0_3piscal = -10572.0+50658.0*xst-87903.0*xst**2+66886.0*xst**3
     & -18699.0*xst**4
      delta0_3piscal = delta0_3piscal*pi/180.0
      else
      delta0_3piscal = 255.0*pi/180.0
      endif
      delta0_3piscal = atan(delta0_3piscal)
      R0SCAL_3PI = ALPHA0_3PI/QX +
     &          ALPHA1_3PI/QX**2*(SX - MMF0**2)
      R0SCAL_3PI = R0SCAL_3PI*(cos(delta0_3piscal) +
     &              i*sin(delta0_3piscal))
      RETURN
      END
      FUNCTION R2SCAL_3PI(QX,SX)
C****************************************************************
C     Complex Function: R_2 function (Pablo notes)
C     for the scalar contribution for three pion modes
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      COMPLEX R2SCAL_3PI
      REAL            QX,SX,delta2_3piscal,xsx,xst
      DOUBLE PRECISION    XM1,XM2,DSX
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
c$$$      a02_3piscal = -0.0444
c$$$      b02_3piscal = -0.0857/mmpi_av**2
c$$$      c02_3piscal = -0.00221/mmpi_av**4
c$$$      d02_3piscal = -0.000129/mmpi_av**6
c$$$      x02_3piscal = -21.62*mmpi_av**2
c$$$c      MMF0 = 0.98
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      dsx = sx
      xsx = sx/4.*sigp(dsx)**2
      xst = sqrt(sx)
      if(sx.le.0.7) then
      delta2_3piscal = sigp(dsx)*(a02_3piscal + b02_3piscal*xsx +
     &               c02_3piscal*xsx**2 + d02_3piscal*xsx**3)
      delta2_3piscal = delta2_3piscal*
     &         (4*mmpi_av**2 - x02_3piscal)/(sx-x02_3piscal)
      else if(xst.le.1.21) then
      delta2_3piscal = 282.9-1314.9*xst+2153.4*xst**2-1574.5d0*xst**3+
     & 428.06d0*xst**4
      delta2_3piscal = delta2_3piscal*pi/180.0
      else
      delta2_3piscal = -27.0*pi/180.0
      endif
      delta2_3piscal = atan(delta2_3piscal)
      R2SCAL_3PI = GAMMA0_3PI/QX +
     &          GAMMA1_3PI/QX**2*(SX - MMF0**2)
      R2SCAL_3PI = R2SCAL_3PI*(cos(delta2_3piscal) +
     &              i*sin(delta2_3piscal))
      RETURN
      END
      FUNCTION FFsig(QX,XX)
C **********************************************************
C     Complex Function:
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      COMPLEX FFSIG
      REAL            QX,XX,mm2
      DOUBLE PRECISION    XM1,XM2,xphi
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      mm2 = MMPI_AV**2
      xphi = - rsigma**2* LAMB_RCHT(QX,XX,mm2)/(8.*QX)
      FFsig = dexp(xphi)
       RETURN
       END
      FUNCTION BWsig(XM,XG,XQ)
C **********************************************************
C     Complex Function: S-wave Breit-Wigner
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      COMPLEX BWSIG
      REAL            XQ
      DOUBLE PRECISION    XM,XG,XM2,XXQ,GAMMA
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      XXQ = XQ
      XM2 = XM**2
      GAMMA = XG*SIGP(XXQ)/SIGP(XM2)
      BWsig = XM*XM/CMPLX(XM*XM-XQ, -XM*GAMMA)
       RETURN
       END
      FUNCTION fattcoul(mm1,mm3,ss2)
C **********************************************************
C     Real Function: Coulomb attraction of two particles
C                    with mm1 and mm3
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      REAL FATTCOUL
      REAL            ss2
      DOUBLE PRECISION    mm1,mm3,betam1m3
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
C*******************************************
C   COMMON  block fixed in SUBROUTINE INIPHY
c*******************************************
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_qedprm') /QEDPRM/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      if(ss2.gt.(mm1+mm3)**2) then
      betam1m3 = 2.*dsqrt(1.d0 - (mm1 +mm3)**2/ss2)
     &            /(1.+ (1.d0 - (mm1 +mm3)**2/ss2))
         fattcoul = 2.*pi/betam1m3/ALFINV
     &               /(1.-exp(-2.*pi/betam1m3/ALFINV))
       else
         fattcoul = 1
          endif
       RETURN
       END
      FUNCTION frepcoul(mm1,mm3,ss2)
C **********************************************************
C     Real Function: Coulomb repuslcion of two particles
C                    with mm1 and mm3
C     Called by f3pi_rcht.f
C **********************************************************
      IMPLICIT NONE
      REAL FREPCOUL
      REAL            ss2
      DOUBLE PRECISION    mm1,mm3,betam1m3
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
C*******************************************
C   COMMON  block fixed in SUBROUTINE INIPHY
c*******************************************
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_qedprm') /QEDPRM/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
        if(ss2.gt.(mm1+mm3)**2) then
      betam1m3 = 2.*dsqrt(1.d0 - (mm1 +mm3)**2/ss2)
     &           /(1.+ (1.d0 - (mm1 +mm3)**2/ss2))
         frepcoul = 2.*pi/betam1m3/ALFINV
     &               /(-1.+exp(2.*pi/betam1m3/ALFINV))
       else
         frepcoul = 1
          endif
       RETURN
       END
      DOUBLE PRECISION FUNCTION GAUS(F,A,B,EEPS)
*     *************************
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      DIMENSION W(12),X(12)
      DATA CONST /1.0D-19/
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
C-----------------------------------------------------------------------------
C NEW INTEGRATION METHOD - constant step size (provided by fourth parameter)
C comment out following 5 lines to return to old solutions of advantages too
C-----------------------------------------------------------------------------
      RESULT=0
      EPS=0
c      CALL STEPGAUSS(F,A,B,6,RESULT,EPS)
      CALL CHANGEGAUSS(F,A,B,2,RESULT,EPS)
      GAUS=RESULT
      RETURN
C-----------------------------------------------------------------------------
C PREVIOUS INTEGRATION METHOD - step size depends on EPSSQ
C-----------------------------------------------------------------------------
      EPS=dABS(EEPS)
      DELTA=CONST*dABS(A-B)
      GAUS=0.d0
      AA=A
    5 Y=B-AA
      IF(dABS(Y) .LE. DELTA) RETURN
    2 BB=AA+Y
      C1=0.5d0*(AA+BB)
      C2=C1-AA
      S8=0.d0
      S16=0.d0
      DO 1 I=1,4
      U=X(I)*C2
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 3 I=5,12
      U=X(I)*C2
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF(EEPS.LT.0D0) THEN
        IF(dABS(S16-S8) .GT. EPS*dABS(S16)) GO TO 4
      ELSE
        IF(dABS(S16-S8) .GT. EPS) GO TO 4
      ENDIF
      GAUS=GAUS+S16
      AA=BB
      GO TO 5
    4 Y=0.5d0*Y
      IF(dABS(Y) .GT. DELTA) GOTO 2
      CALL TAUOLA_PRINT('tauola.f:8292: (suppressed)'//CHAR(0))
!     PRINT 7
      GAUS=0.d0
      RETURN
!   7 FORMAT(1X,36HGAUS  ... TOO HIGH ACCURACY REQUIRED)
      END
      DOUBLE PRECISION FUNCTION GAUS2(F,A,B,EEPS)
*     *************************
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL F
      DIMENSION W(12),X(12)
      DATA CONST /1.0D-19/
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
C-----------------------------------------------------------------------------
C NEW INTEGRATION METHOD - constant step size (provided by fourth parameter)
C comment out following 5 lines to return to old solutions of advantages too
C-----------------------------------------------------------------------------
      RESULT=0
      EPS=0
c      CALL STEPGAUSS2(F,A,B,6,RESULT,EPS)
      CALL CHANGEGAUSS2(F,A,B,2,RESULT,EPS)
      GAUS2=RESULT
      RETURN
C-----------------------------------------------------------------------------
C PREVIOUS INTEGRATION METHOD - step size depends on EPSSQ
C-----------------------------------------------------------------------------
      EPS=dABS(EEPS)
      DELTA=CONST*dABS(A-B)
      GAUS2=0.D0
      AA=A
    5 Y=B-AA
      IF(dABS(Y) .LE. DELTA) RETURN
    2 BB=AA+Y
      C1=0.5d0*(AA+BB)
      C2=C1-AA
      S8=0.d0
      S16=0.d0
      DO 1 I=1,4
      U=X(I)*C2
    1 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      DO 3 I=5,12
      U=X(I)*C2
    3 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S8=S8*C2
      S16=S16*C2
      IF(EEPS.LT.0D0) THEN
        IF(dABS(S16-S8) .GT. EPS*dABS(S16)) GO TO 4
      ELSE
        IF(dABS(S16-S8) .GT. EPS) GO TO 4
      ENDIF
      GAUS2=GAUS2+S16
      AA=BB
      GO TO 5
    4 Y=0.5d0*Y
      IF(dABS(Y) .GT. DELTA) GOTO 2
      CALL TAUOLA_PRINT('tauola.f:8356: (suppressed)'//CHAR(0))
!     PRINT 7
      GAUS2=0.D0
      RETURN
!   7 FORMAT(1X,36HGAUS2 ... TOO HIGH ACCURACY REQUIRED)
      END
      DOUBLE PRECISION FUNCTION F_CHANGE(X,F,AMS1,AMS2)
      IMPLICIT NONE
      EXTERNAL F
      DOUBLE PRECISION F
      DOUBLE PRECISION X, NEW_X
      DOUBLE PRECISION AMS1,AMS2
      DOUBLE PRECISION ALP1,ALP2,ALP
      DOUBLE PRECISION AMRX,GAMRX
      DATA AMRX/0.77/,GAMRX/1.8/
      ALP1  = ATAN((AMS1-AMRX**2)/AMRX/GAMRX) ! integration range for ALP
      ALP2  = ATAN((AMS2-AMRX**2)/AMRX/GAMRX) ! integration range for ALP
      ALP   = ALP1 + X*(ALP2-ALP1)            ! change of variables
      NEW_X = AMRX**2+AMRX*GAMRX*TAN(ALP)     ! second change of variables
      F_CHANGE = F(NEW_X)
! Jacobian for change of variables
      F_CHANGE = F_CHANGE * (ALP2-ALP1) * ((NEW_X-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)
      RETURN
      END
      DOUBLE PRECISION FUNCTION F_CHANGE2(X,F,AMS1,AMS2)
      IMPLICIT NONE
      EXTERNAL F
      DOUBLE PRECISION F
      DOUBLE PRECISION X, NEW_X
      DOUBLE PRECISION AMS1,AMS2
      DOUBLE PRECISION ALP1,ALP2,ALP
      DOUBLE PRECISION AMRX,GAMRX
      DATA AMRX/0.77/,GAMRX/1.8/
      ALP1  = ATAN((AMS1-AMRX**2)/AMRX/GAMRX) ! integration range for ALP
      ALP2  = ATAN((AMS2-AMRX**2)/AMRX/GAMRX) ! integration range for ALP
      ALP   = ALP1 + X*(ALP2-ALP1)            ! change of variables
      NEW_X = AMRX**2+AMRX*GAMRX*TAN(ALP)     ! second change of variables
      F_CHANGE2 = F(NEW_X)
! Jacobian for change of variables
      F_CHANGE2 = F_CHANGE2 * (ALP2-ALP1) * ((NEW_X-AMRX**2)**2+(AMRX*GAMRX)**2)/(AMRX*GAMRX)
      RETURN
      END
      SUBROUTINE CHANGEGAUSS(F,AMS1,AMS2,LAMBDA,RESULT,EPS)
      implicit real*8 (a-h,o-z)
      EXTERNAL F,F2
      DIMENSION W(12),X(12)
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      RESULT=0.
      EPS=0.
! Range is changed from [AMS1,AMS2] to [0,1]
! F, AMS1 and AMS2 are parameters of F_CHANGE
      A=0
      B=1
      AU=A
      C=(B-A)/FLOAT(LAMBDA)
      DO 1 I = 1,LAMBDA
      FI=I
      AO=A+FI*C
      C1=0.5*(AU+AO)
      C2=C1-AU
      S8=0.
      S16=0.
      DO 2 J = 1,4
      U=X(J)*C2
    2 S8=S8+W(J)*(F_CHANGE(C1+U,F,AMS1,AMS2)+F_CHANGE(C1-U,F,AMS1,AMS2))
      DO 3 J = 5,12
      U=X(J)*C2
    3 S16=S16+W(J)*(F_CHANGE(C1+U,F,AMS1,AMS2)+F_CHANGE(C1-U,F,AMS1,AMS2))
      RESULT=RESULT+C2*S16
      EPS=EPS+ABS(C2*(S16-S8))
    1 AU=AO
      RETURN
      END
      SUBROUTINE CHANGEGAUSS2(F,AMS1,AMS2,LAMBDA,RESULT,EPS)
      implicit real*8 (a-h,o-z)
      EXTERNAL F,F2
      DIMENSION W(12),X(12)
      DATA W
     1/0.10122 85362 90376, 0.22238 10344 53374, 0.31370 66458 77887,
     2 0.36268 37833 78362, 0.02715 24594 11754, 0.06225 35239 38648,
     3 0.09515 85116 82493, 0.12462 89712 55534, 0.14959 59888 16577,
     4 0.16915 65193 95003, 0.18260 34150 44924, 0.18945 06104 55069/
      DATA X
     1/0.96028 98564 97536, 0.79666 64774 13627, 0.52553 24099 16329,
     2 0.18343 46424 95650, 0.98940 09349 91650, 0.94457 50230 73233,
     3 0.86563 12023 87832, 0.75540 44083 55003, 0.61787 62444 02644,
     4 0.45801 67776 57227, 0.28160 35507 79259, 0.09501 25098 37637/
      RESULT=0.
      EPS=0.
! Range is changed from [AMS1,AMS2] to [0,1]
! F, AMS1 and AMS2 are parameters of F_CHANGE
      A=0
      B=1
      AU=A
      C=(B-A)/FLOAT(LAMBDA)
      DO 1 I = 1,LAMBDA
      FI=I
      AO=A+FI*C
      C1=0.5*(AU+AO)
      C2=C1-AU
      S8=0.
      S16=0.
      DO 2 J = 1,4
      U=X(J)*C2
    2 S8=S8+W(J)*(F_CHANGE2(C1+U,F,AMS1,AMS2)+F_CHANGE2(C1-U,F,AMS1,AMS2))
      DO 3 J = 5,12
      U=X(J)*C2
    3 S16=S16+W(J)*(F_CHANGE2(C1+U,F,AMS1,AMS2)+F_CHANGE2(C1-U,F,AMS1,AMS2))
      RESULT=RESULT+C2*S16
      EPS=EPS+ABS(C2*(S16-S8))
    1 AU=AO
      RETURN
      END
      DOUBLE PRECISION FUNCTION GFACT(QQ)
      IMPLICIT NONE
C factor G to be used as inteligent retabulation as in paper
C Kuhn Santamaria
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      DOUBLE PRECISION QQ,THR,THU,TH0,LAM
      INTEGER II,INIT,INUM
      DOUBLE PRECISION AF4,AF5,AA1,AA2,AA4
      DOUBLE PRECISION F1,F2,F3,F4,F5
      DOUBLE PRECISION AA,BB,CC,DD,A,B,C,D,X
      SAVE A,B,C,D,AF4,AF5,AA,BB,CC,DD,THR,THU,TH0
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      CALL IFGFACT(1,II,INIT)
      IF(INIT.EQ.0) THEN
       THR=(AMPI+AMRO)**2
       TH0=9*AMPI**2
       THU=THR
       CALL getFF3PISCAL(INUM)  ! we switch off part of the contr to 3 pi mode
       CALL setFF3PISCAL(0)
C below THR , we calculate at distance 1/4 1/2 and 1 between
C minimum and maximum for this range
       LAM=(THR-TH0)/4
       AA1=DGAMQQ(TH0+  LAM)
       AA2=DGAMQQ(TH0+2*LAM)
       AA4=DGAMQQ(TH0+4*LAM)
C above THR we calculate at THR times 1 2 3 and 4 for higher range
       F1=DGAMQQ(THU)
       F2=DGAMQQ(1.5*THU)
       F3=DGAMQQ(2.*THU)
       F4=DGAMQQ(3*THU)
       F5=DGAMQQ(3.5*THU)
       CALL setFF3PISCAL(INUM) ! we switch back part of the contr to 3 pi mode
C we calculate coefs for expansion  ( polynomial of order 2 once
C                                     X**2 factorized out, X=QQ-TH0 )
       AA1 = AA1/LAM**2            ! RESCALING due to factorized X**2
       AA2 = AA2/(2.*LAM)**2       !
       AA4 = AA4/(4.*LAM)**2       !
       BB=1./8.*(10.*AA2-AA4-16.*AA1)
       AA=(8.*AA1-AA2-4.*BB)/6.
       CC=AA1-AA-BB
       AA=AA/LAM
       BB=BB/LAM**2
       CC=CC/LAM**3
C calulate coefs, assuming it is polynomial order 3 note negative pwrs
       D=-9*(-4*F3-F1+4*F2+F4)
       C=3*(F4+F1-2*F3-11./18.*D)
       A=F3-F1+0.5*C+0.75*D
       B=F1-A-C-D
       A=A/THU
       B=B
       C=C*THU
       D=D*THU**2
       AF4=F4
       AF5=F5
c       write(*,*) "A=",AA,"B=",-BB,"C=",CC,"D=",A,"E=",-B,"F=",C
c       write(*,*) "G=",-D,"H=",AF4,"P=",AF5-AF4
      ENDIF
      IF (QQ.GT.3*THU) THEN
        GFACT=AF4+(AF5-AF4)*(QQ-3*THU)*2/THU
      ELSEIF (QQ.GT.THR) THEN
        GFACT=A*QQ+B+C/QQ+D/QQ**2
      ELSEIF(QQ.LE.TH0) THEN
        GFACT=0.0
      ELSE
        X=QQ-TH0
        GFACT=AA*X+BB*X**2+CC*X**3
        GFACT=X**2*GFACT
      ENDIF
      END
      DOUBLE PRECISION FUNCTION DGAMQQ(XQQB)
C **************************************************************
C     calculates \tau^- -> pi^- pi^- pi^+ nu width as function of QQ (XQQB)
C     formulas (19) of ref [2a] integration over S1
C     limit of integration (21) of ref [2a]  see also [4]
C     called from main function
C **************************************************************
      IMPLICIT NONE
      COMMON/PRECINT/          EPSSQ,ABS1
      DOUBLE PRECISION         EPSSQ,ABS1
      COMMON /EXTERNAL/ XQQA
      DOUBLE PRECISION  XQQA
      DOUBLE PRECISION  XQQB,EPS,UPS1,DOWNS1
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
!
      BIND(C,NAME='tauola_external') /EXTERNAL/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_precint') /PRECINT/
!
      XQQA = XQQB
      EPS = EPSSQ/3.D0
      UPS1=(DSQRT(XQQB) - AMPI)**2-ABS1              ! limits on S1
      DOWNS1=4.D0*AMPI**2+ABS1
      DGAMQQ = GAUS(DGAMQQS1,DOWNS1,UPS1,EPS)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DGAMQQS1(S1)
      IMPLICIT NONE
      DOUBLE PRECISION                   S1
C **************************************************************
C     calculates \tau^- -> pi^- pi^- pi^+ nu width
C           as function of QQ,S1
C     GAUS integrant in DGAMQQ (XQQA- hidden argument)
C     calculates tau^- -> pi^- pi^- pi^+ nu spectrum as function of S1
C     formulas (19) of ref [2a] see also [4]
C     limit of integration (21) of ref [2a]
C **************************************************************
      EXTERNAL          FFWID3PI
      DOUBLE PRECISION  FFWID3PI
      COMMON /INTERNAL/ S1A
      DOUBLE PRECISION  S1A
      COMMON /EXTERNAL/   XQQA
      DOUBLE PRECISION    XQQA
      COMMON/PRECINT/          EPSSQ,ABS1
      DOUBLE PRECISION         EPSSQ,ABS1
      DOUBLE PRECISION EPS,UPS3,DOWNS3
      DOUBLE PRECISION XLAM,X,Y,Z
      DOUBLE PRECISION XAMPI2
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
!
      BIND(C,NAME='tauola_external') /EXTERNAL/
      BIND(C,NAME='tauola_internal') /INTERNAL/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_precint') /PRECINT/
!
      Xlam(x,y,z) = sqrt(abs((x-y-z)**2 - 4.*y*z))
      S1A = S1
      EPS = EPSSQ/9.D0
      XAMPI2=AMPI**2
      UPS3 = (XQQA - AMPI**2)**2 -                              ! limits on S3
     &    ( XLAM(XQQA,S1,XAMPI2)
     &         - XLAM(S1,XAMPI2,XAMPI2) )**2
      DOWNS3 = (XQQA - AMPI**2)**2 -
     &    (XLAM(XQQA,S1,XAMPI2)
     &          + XLAM(S1,XAMPI2,XAMPI2) )**2
      UPS3 = UPS3/4./S1
      DOWNS3 = DOWNS3/4./S1
      DGAMQQS1 = GAUS2(DGAMQQS1S3,DOWNS3,UPS3,EPS)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DGAMQQS1S3(XS3)
      IMPLICIT NONE
C **************************************************************
C     calculates \tau^- -> pi^- pi^- pi^+ nu width
C           as function of QQ,S1,S3
C **************************************************************
C      EXTERNAL          GAUS2,FFWID3PI
      DOUBLE PRECISION  XS3
      COMMON /INTERNAL/ XS1A
      DOUBLE PRECISION  XS1A
      COMMON /EXTERNAL/   XQQA
      DOUBLE PRECISION    XQQA
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
!
      BIND(C,NAME='tauola_external') /EXTERNAL/
      BIND(C,NAME='tauola_internal') /INTERNAL/
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      DGAMQQS1S3 = FFWID3PI(XQQA,XS1A,XS3)
      RETURN
      END
       SUBROUTINE InitA1TAB(qk_min    ,w_qmin    )
       PARAMETER (Nq=1001)
       real*8               qk_min(Nq),w_qmin(Nq)
 ! this routine was automatically created by  a1-tabler
       qk_min(   1)=0.171492857E+00
       w_qmin(   1)=0.000000000E+00
       qk_min(   2)=0.174479093E+00
       w_qmin(   2)=0.629755564E-07
       qk_min(   3)=0.177465328E+00
       w_qmin(   3)=0.492358834E-06
       qk_min(   4)=0.180451564E+00
       w_qmin(   4)=0.162470857E-05
       qk_min(   5)=0.183437800E+00
       w_qmin(   5)=0.376709625E-05
       qk_min(   6)=0.186424036E+00
       w_qmin(   6)=0.720012925E-05
       qk_min(   7)=0.189410272E+00
       w_qmin(   7)=0.121806627E-04
       qk_min(   8)=0.192396508E+00
       w_qmin(   8)=0.189441916E-04
       qk_min(   9)=0.195382744E+00
       w_qmin(   9)=0.277070270E-04
       qk_min(  10)=0.198368980E+00
       w_qmin(  10)=0.386681867E-04
       qk_min(  11)=0.201355216E+00
       w_qmin(  11)=0.520111478E-04
       qk_min(  12)=0.204341452E+00
       w_qmin(  12)=0.679053591E-04
       qk_min(  13)=0.207327688E+00
       w_qmin(  13)=0.865076882E-04
       qk_min(  14)=0.210313924E+00
       w_qmin(  14)=0.107963592E-03
       qk_min(  15)=0.213300160E+00
       w_qmin(  15)=0.132408251E-03
       qk_min(  16)=0.216286396E+00
       w_qmin(  16)=0.159967766E-03
       qk_min(  17)=0.219272632E+00
       w_qmin(  17)=0.190759703E-03
       qk_min(  18)=0.222258868E+00
       w_qmin(  18)=0.224894342E-03
       qk_min(  19)=0.225245104E+00
       w_qmin(  19)=0.262475065E-03
       qk_min(  20)=0.228231340E+00
       w_qmin(  20)=0.303599320E-03
       qk_min(  21)=0.231217576E+00
       w_qmin(  21)=0.348358980E-03
       qk_min(  22)=0.234203812E+00
       w_qmin(  22)=0.396841130E-03
       qk_min(  23)=0.237190048E+00
       w_qmin(  23)=0.449128509E-03
       qk_min(  24)=0.240176284E+00
       w_qmin(  24)=0.505299824E-03
       qk_min(  25)=0.243162520E+00
       w_qmin(  25)=0.565430276E-03
       qk_min(  26)=0.246148756E+00
       w_qmin(  26)=0.629592245E-03
       qk_min(  27)=0.249134992E+00
       w_qmin(  27)=0.697855178E-03
       qk_min(  28)=0.252121228E+00
       w_qmin(  28)=0.770286144E-03
       qk_min(  29)=0.255107464E+00
       w_qmin(  29)=0.846949458E-03
       qk_min(  30)=0.258093700E+00
       w_qmin(  30)=0.927908699E-03
       qk_min(  31)=0.261079936E+00
       w_qmin(  31)=0.101322560E-02
       qk_min(  32)=0.264066172E+00
       w_qmin(  32)=0.110295947E-02
       qk_min(  33)=0.267052408E+00
       w_qmin(  33)=0.119716975E-02
       qk_min(  34)=0.270038643E+00
       w_qmin(  34)=0.129591394E-02
       qk_min(  35)=0.273024879E+00
       w_qmin(  35)=0.139924978E-02
       qk_min(  36)=0.276011115E+00
       w_qmin(  36)=0.150723308E-02
       qk_min(  37)=0.278997351E+00
       w_qmin(  37)=0.161992059E-02
       qk_min(  38)=0.281983587E+00
       w_qmin(  38)=0.173736751E-02
       qk_min(  39)=0.284969823E+00
       w_qmin(  39)=0.185962999E-02
       qk_min(  40)=0.287956059E+00
       w_qmin(  40)=0.198676282E-02
       qk_min(  41)=0.290942295E+00
       w_qmin(  41)=0.211882209E-02
       qk_min(  42)=0.293928531E+00
       w_qmin(  42)=0.225586359E-02
       qk_min(  43)=0.296914767E+00
       w_qmin(  43)=0.239794313E-02
       qk_min(  44)=0.299901003E+00
       w_qmin(  44)=0.254511773E-02
       qk_min(  45)=0.302887239E+00
       w_qmin(  45)=0.269744335E-02
       qk_min(  46)=0.305873475E+00
       w_qmin(  46)=0.285497820E-02
       qk_min(  47)=0.308859711E+00
       w_qmin(  47)=0.301777964E-02
       qk_min(  48)=0.311845947E+00
       w_qmin(  48)=0.318590749E-02
       qk_min(  49)=0.314832183E+00
       w_qmin(  49)=0.335942016E-02
       qk_min(  50)=0.317818419E+00
       w_qmin(  50)=0.353837928E-02
       qk_min(  51)=0.320804655E+00
       w_qmin(  51)=0.372284553E-02
       qk_min(  52)=0.323790891E+00
       w_qmin(  52)=0.391288143E-02
       qk_min(  53)=0.326777127E+00
       w_qmin(  53)=0.410855094E-02
       qk_min(  54)=0.329763363E+00
       w_qmin(  54)=0.430991813E-02
       qk_min(  55)=0.332749599E+00
       w_qmin(  55)=0.451704899E-02
       qk_min(  56)=0.335735835E+00
       w_qmin(  56)=0.473001097E-02
       qk_min(  57)=0.338722071E+00
       w_qmin(  57)=0.494887236E-02
       qk_min(  58)=0.341708307E+00
       w_qmin(  58)=0.517370225E-02
       qk_min(  59)=0.344694543E+00
       w_qmin(  59)=0.540457306E-02
       qk_min(  60)=0.347680779E+00
       w_qmin(  60)=0.564155599E-02
       qk_min(  61)=0.350667015E+00
       w_qmin(  61)=0.588472630E-02
       qk_min(  62)=0.353653251E+00
       w_qmin(  62)=0.613416049E-02
       qk_min(  63)=0.356639487E+00
       w_qmin(  63)=0.638993469E-02
       qk_min(  64)=0.359625722E+00
       w_qmin(  64)=0.665212895E-02
       qk_min(  65)=0.362611958E+00
       w_qmin(  65)=0.692082301E-02
       qk_min(  66)=0.365598194E+00
       w_qmin(  66)=0.719610191E-02
       qk_min(  67)=0.368584430E+00
       w_qmin(  67)=0.747804819E-02
       qk_min(  68)=0.371570666E+00
       w_qmin(  68)=0.776674972E-02
       qk_min(  69)=0.374556902E+00
       w_qmin(  69)=0.806229452E-02
       qk_min(  70)=0.377543138E+00
       w_qmin(  70)=0.836477374E-02
       qk_min(  71)=0.380529374E+00
       w_qmin(  71)=0.867427853E-02
       qk_min(  72)=0.383515610E+00
       w_qmin(  72)=0.899090516E-02
       qk_min(  73)=0.386501846E+00
       w_qmin(  73)=0.931475228E-02
       qk_min(  74)=0.389488082E+00
       w_qmin(  74)=0.964591670E-02
       qk_min(  75)=0.392474318E+00
       w_qmin(  75)=0.998450149E-02
       qk_min(  76)=0.395460554E+00
       w_qmin(  76)=0.103306093E-01
       qk_min(  77)=0.398446790E+00
       w_qmin(  77)=0.106843494E-01
       qk_min(  78)=0.401433026E+00
       w_qmin(  78)=0.110458285E-01
       qk_min(  79)=0.404419262E+00
       w_qmin(  79)=0.114151601E-01
       qk_min(  80)=0.407405498E+00
       w_qmin(  80)=0.117924571E-01
       qk_min(  81)=0.410391734E+00
       w_qmin(  81)=0.121778394E-01
       qk_min(  82)=0.413377970E+00
       w_qmin(  82)=0.125714246E-01
       qk_min(  83)=0.416364206E+00
       w_qmin(  83)=0.129733356E-01
       qk_min(  84)=0.419350442E+00
       w_qmin(  84)=0.133837007E-01
       qk_min(  85)=0.422336678E+00
       w_qmin(  85)=0.138026463E-01
       qk_min(  86)=0.425322914E+00
       w_qmin(  86)=0.142303061E-01
       qk_min(  87)=0.428309150E+00
       w_qmin(  87)=0.146668143E-01
       qk_min(  88)=0.431295386E+00
       w_qmin(  88)=0.151123120E-01
       qk_min(  89)=0.434281622E+00
       w_qmin(  89)=0.155669372E-01
       qk_min(  90)=0.437267858E+00
       w_qmin(  90)=0.160308391E-01
       qk_min(  91)=0.440254094E+00
       w_qmin(  91)=0.165041642E-01
       qk_min(  92)=0.443240330E+00
       w_qmin(  92)=0.169870675E-01
       qk_min(  93)=0.446226566E+00
       w_qmin(  93)=0.174797039E-01
       qk_min(  94)=0.449212802E+00
       w_qmin(  94)=0.179822341E-01
       qk_min(  95)=0.452199037E+00
       w_qmin(  95)=0.184948231E-01
       qk_min(  96)=0.455185273E+00
       w_qmin(  96)=0.190176372E-01
       qk_min(  97)=0.458171509E+00
       w_qmin(  97)=0.195508546E-01
       qk_min(  98)=0.461157745E+00
       w_qmin(  98)=0.200946463E-01
       qk_min(  99)=0.464143981E+00
       w_qmin(  99)=0.206491980E-01
       qk_min( 100)=0.467130217E+00
       w_qmin( 100)=0.212146893E-01
       qk_min( 101)=0.470116453E+00
       w_qmin( 101)=0.217913198E-01
       qk_min( 102)=0.473102689E+00
       w_qmin( 102)=0.223792750E-01
       qk_min( 103)=0.476088925E+00
       w_qmin( 103)=0.229787614E-01
       qk_min( 104)=0.479075161E+00
       w_qmin( 104)=0.235899856E-01
       qk_min( 105)=0.482061397E+00
       w_qmin( 105)=0.242131536E-01
       qk_min( 106)=0.485047633E+00
       w_qmin( 106)=0.248484858E-01
       qk_min( 107)=0.488033869E+00
       w_qmin( 107)=0.254961973E-01
       qk_min( 108)=0.491020105E+00
       w_qmin( 108)=0.261565231E-01
       qk_min( 109)=0.494006341E+00
       w_qmin( 109)=0.268296892E-01
       qk_min( 110)=0.496992577E+00
       w_qmin( 110)=0.275159405E-01
       qk_min( 111)=0.499978813E+00
       w_qmin( 111)=0.282155152E-01
       qk_min( 112)=0.502965049E+00
       w_qmin( 112)=0.289286682E-01
       qk_min( 113)=0.505951285E+00
       w_qmin( 113)=0.296556595E-01
       qk_min( 114)=0.508937521E+00
       w_qmin( 114)=0.303967464E-01
       qk_min( 115)=0.511923757E+00
       w_qmin( 115)=0.311522080E-01
       qk_min( 116)=0.514909993E+00
       w_qmin( 116)=0.319223193E-01
       qk_min( 117)=0.517896229E+00
       w_qmin( 117)=0.327073638E-01
       qk_min( 118)=0.520882465E+00
       w_qmin( 118)=0.335076407E-01
       qk_min( 119)=0.523868701E+00
       w_qmin( 119)=0.343234413E-01
       qk_min( 120)=0.526854937E+00
       w_qmin( 120)=0.351550793E-01
       qk_min( 121)=0.529841173E+00
       w_qmin( 121)=0.360028709E-01
       qk_min( 122)=0.532827409E+00
       w_qmin( 122)=0.368671391E-01
       qk_min( 123)=0.535813645E+00
       w_qmin( 123)=0.377482165E-01
       qk_min( 124)=0.538799881E+00
       w_qmin( 124)=0.386464530E-01
       qk_min( 125)=0.541786116E+00
       w_qmin( 125)=0.395621863E-01
       qk_min( 126)=0.544772352E+00
       w_qmin( 126)=0.404957853E-01
       qk_min( 127)=0.547758588E+00
       w_qmin( 127)=0.414476179E-01
       qk_min( 128)=0.550744824E+00
       w_qmin( 128)=0.424180670E-01
       qk_min( 129)=0.553731060E+00
       w_qmin( 129)=0.434075151E-01
       qk_min( 130)=0.556717296E+00
       w_qmin( 130)=0.444163612E-01
       qk_min( 131)=0.559703532E+00
       w_qmin( 131)=0.454450215E-01
       qk_min( 132)=0.562689768E+00
       w_qmin( 132)=0.464939133E-01
       qk_min( 133)=0.565676004E+00
       w_qmin( 133)=0.475634622E-01
       qk_min( 134)=0.568662240E+00
       w_qmin( 134)=0.486541189E-01
       qk_min( 135)=0.571648476E+00
       w_qmin( 135)=0.497663408E-01
       qk_min( 136)=0.574634712E+00
       w_qmin( 136)=0.509005855E-01
       qk_min( 137)=0.577620948E+00
       w_qmin( 137)=0.520573344E-01
       qk_min( 138)=0.580607184E+00
       w_qmin( 138)=0.532370803E-01
       qk_min( 139)=0.583593420E+00
       w_qmin( 139)=0.544403257E-01
       qk_min( 140)=0.586579656E+00
       w_qmin( 140)=0.556675790E-01
       qk_min( 141)=0.589565892E+00
       w_qmin( 141)=0.569193795E-01
       qk_min( 142)=0.592552128E+00
       w_qmin( 142)=0.581962668E-01
       qk_min( 143)=0.595538364E+00
       w_qmin( 143)=0.594988003E-01
       qk_min( 144)=0.598524600E+00
       w_qmin( 144)=0.608275376E-01
       qk_min( 145)=0.601510836E+00
       w_qmin( 145)=0.621830754E-01
       qk_min( 146)=0.604497072E+00
       w_qmin( 146)=0.635660133E-01
       qk_min( 147)=0.607483308E+00
       w_qmin( 147)=0.649769525E-01
       qk_min( 148)=0.610469544E+00
       w_qmin( 148)=0.664165334E-01
       qk_min( 149)=0.613455780E+00
       w_qmin( 149)=0.678853833E-01
       qk_min( 150)=0.616442016E+00
       w_qmin( 150)=0.693841730E-01
       qk_min( 151)=0.619428252E+00
       w_qmin( 151)=0.709135700E-01
       qk_min( 152)=0.622414488E+00
       w_qmin( 152)=0.724742457E-01
       qk_min( 153)=0.625400724E+00
       w_qmin( 153)=0.740669312E-01
       qk_min( 154)=0.628386960E+00
       w_qmin( 154)=0.756923306E-01
       qk_min( 155)=0.631373196E+00
       w_qmin( 155)=0.773511675E-01
       qk_min( 156)=0.634359431E+00
       w_qmin( 156)=0.790441825E-01
       qk_min( 157)=0.637345667E+00
       w_qmin( 157)=0.807721648E-01
       qk_min( 158)=0.640331903E+00
       w_qmin( 158)=0.825358806E-01
       qk_min( 159)=0.643318139E+00
       w_qmin( 159)=0.843360846E-01
       qk_min( 160)=0.646304375E+00
       w_qmin( 160)=0.861736137E-01
       qk_min( 161)=0.649290611E+00
       w_qmin( 161)=0.880492855E-01
       qk_min( 162)=0.652276847E+00
       w_qmin( 162)=0.899639046E-01
       qk_min( 163)=0.655263083E+00
       w_qmin( 163)=0.919183180E-01
       qk_min( 164)=0.658249319E+00
       w_qmin( 164)=0.939133758E-01
       qk_min( 165)=0.661235555E+00
       w_qmin( 165)=0.959499343E-01
       qk_min( 166)=0.664221791E+00
       w_qmin( 166)=0.980288583E-01
       qk_min( 167)=0.667208027E+00
       w_qmin( 167)=0.100150999E+00
       qk_min( 168)=0.670194263E+00
       w_qmin( 168)=0.102317274E+00
       qk_min( 169)=0.673180499E+00
       w_qmin( 169)=0.104528527E+00
       qk_min( 170)=0.676166735E+00
       w_qmin( 170)=0.106785656E+00
       qk_min( 171)=0.679152971E+00
       w_qmin( 171)=0.109089540E+00
       qk_min( 172)=0.682139207E+00
       w_qmin( 172)=0.111441042E+00
       qk_min( 173)=0.685125443E+00
       w_qmin( 173)=0.113841040E+00
       qk_min( 174)=0.688111679E+00
       w_qmin( 174)=0.116290378E+00
       qk_min( 175)=0.691097915E+00
       w_qmin( 175)=0.118789907E+00
       qk_min( 176)=0.694084151E+00
       w_qmin( 176)=0.121340451E+00
       qk_min( 177)=0.697070387E+00
       w_qmin( 177)=0.123942829E+00
       qk_min( 178)=0.700056623E+00
       w_qmin( 178)=0.126597840E+00
       qk_min( 179)=0.703042859E+00
       w_qmin( 179)=0.129306173E+00
       qk_min( 180)=0.706029095E+00
       w_qmin( 180)=0.132068605E+00
       qk_min( 181)=0.709015331E+00
       w_qmin( 181)=0.134885806E+00
       qk_min( 182)=0.712001567E+00
       w_qmin( 182)=0.137758367E+00
       qk_min( 183)=0.714987803E+00
       w_qmin( 183)=0.140686902E+00
       qk_min( 184)=0.717974039E+00
       w_qmin( 184)=0.143671942E+00
       qk_min( 185)=0.720960275E+00
       w_qmin( 185)=0.146713860E+00
       qk_min( 186)=0.723946510E+00
       w_qmin( 186)=0.149813125E+00
       qk_min( 187)=0.726932746E+00
       w_qmin( 187)=0.152970028E+00
       qk_min( 188)=0.729918982E+00
       w_qmin( 188)=0.156184758E+00
       qk_min( 189)=0.732905218E+00
       w_qmin( 189)=0.159457427E+00
       qk_min( 190)=0.735891454E+00
       w_qmin( 190)=0.162788087E+00
       qk_min( 191)=0.738877690E+00
       w_qmin( 191)=0.166176643E+00
       qk_min( 192)=0.741863926E+00
       w_qmin( 192)=0.169622853E+00
       qk_min( 193)=0.744850162E+00
       w_qmin( 193)=0.173126448E+00
       qk_min( 194)=0.747836398E+00
       w_qmin( 194)=0.176686927E+00
       qk_min( 195)=0.750822634E+00
       w_qmin( 195)=0.180303670E+00
       qk_min( 196)=0.753808870E+00
       w_qmin( 196)=0.183975975E+00
       qk_min( 197)=0.756795106E+00
       w_qmin( 197)=0.187702937E+00
       qk_min( 198)=0.759781342E+00
       w_qmin( 198)=0.191483510E+00
       qk_min( 199)=0.762767578E+00
       w_qmin( 199)=0.195316488E+00
       qk_min( 200)=0.765753814E+00
       w_qmin( 200)=0.199200519E+00
       qk_min( 201)=0.768740050E+00
       w_qmin( 201)=0.203134049E+00
       qk_min( 202)=0.771726286E+00
       w_qmin( 202)=0.207115438E+00
       qk_min( 203)=0.774712522E+00
       w_qmin( 203)=0.211142814E+00
       qk_min( 204)=0.777698758E+00
       w_qmin( 204)=0.215214139E+00
       qk_min( 205)=0.780684994E+00
       w_qmin( 205)=0.219327332E+00
       qk_min( 206)=0.783671230E+00
       w_qmin( 206)=0.223479985E+00
       qk_min( 207)=0.786657466E+00
       w_qmin( 207)=0.227669736E+00
       qk_min( 208)=0.789643702E+00
       w_qmin( 208)=0.231893987E+00
       qk_min( 209)=0.792629938E+00
       w_qmin( 209)=0.236149965E+00
       qk_min( 210)=0.795616174E+00
       w_qmin( 210)=0.240434886E+00
       qk_min( 211)=0.798602410E+00
       w_qmin( 211)=0.244745863E+00
       qk_min( 212)=0.801588646E+00
       w_qmin( 212)=0.249079798E+00
       qk_min( 213)=0.804574882E+00
       w_qmin( 213)=0.253433668E+00
       qk_min( 214)=0.807561118E+00
       w_qmin( 214)=0.257804281E+00
       qk_min( 215)=0.810547354E+00
       w_qmin( 215)=0.262188468E+00
       qk_min( 216)=0.813533590E+00
       w_qmin( 216)=0.266583002E+00
       qk_min( 217)=0.816519825E+00
       w_qmin( 217)=0.270984691E+00
       qk_min( 218)=0.819506061E+00
       w_qmin( 218)=0.275390257E+00
       qk_min( 219)=0.822492297E+00
       w_qmin( 219)=0.279796413E+00
       qk_min( 220)=0.825478533E+00
       w_qmin( 220)=0.284200109E+00
       qk_min( 221)=0.828464769E+00
       w_qmin( 221)=0.288598073E+00
       qk_min( 222)=0.831451005E+00
       w_qmin( 222)=0.292987345E+00
       qk_min( 223)=0.834437241E+00
       w_qmin( 223)=0.297364830E+00
       qk_min( 224)=0.837423477E+00
       w_qmin( 224)=0.301727636E+00
       qk_min( 225)=0.840409713E+00
       w_qmin( 225)=0.306072843E+00
       qk_min( 226)=0.843395949E+00
       w_qmin( 226)=0.310397878E+00
       qk_min( 227)=0.846382185E+00
       w_qmin( 227)=0.314699925E+00
       qk_min( 228)=0.849368421E+00
       w_qmin( 228)=0.318976548E+00
       qk_min( 229)=0.852354657E+00
       w_qmin( 229)=0.323225378E+00
       qk_min( 230)=0.855340893E+00
       w_qmin( 230)=0.327444127E+00
       qk_min( 231)=0.858327129E+00
       w_qmin( 231)=0.331630590E+00
       qk_min( 232)=0.861313365E+00
       w_qmin( 232)=0.335782836E+00
       qk_min( 233)=0.864299601E+00
       w_qmin( 233)=0.339898887E+00
       qk_min( 234)=0.867285837E+00
       w_qmin( 234)=0.343977060E+00
       qk_min( 235)=0.870272073E+00
       w_qmin( 235)=0.348015655E+00
       qk_min( 236)=0.873258309E+00
       w_qmin( 236)=0.352013175E+00
       qk_min( 237)=0.876244545E+00
       w_qmin( 237)=0.355968259E+00
       qk_min( 238)=0.879230781E+00
       w_qmin( 238)=0.359879663E+00
       qk_min( 239)=0.882217017E+00
       w_qmin( 239)=0.363746132E+00
       qk_min( 240)=0.885203253E+00
       w_qmin( 240)=0.367566742E+00
       qk_min( 241)=0.888189489E+00
       w_qmin( 241)=0.371340555E+00
       qk_min( 242)=0.891175725E+00
       w_qmin( 242)=0.375066698E+00
       qk_min( 243)=0.894161961E+00
       w_qmin( 243)=0.378744459E+00
       qk_min( 244)=0.897148197E+00
       w_qmin( 244)=0.382373302E+00
       qk_min( 245)=0.900134433E+00
       w_qmin( 245)=0.385952648E+00
       qk_min( 246)=0.903120669E+00
       w_qmin( 246)=0.389481991E+00
       qk_min( 247)=0.906106904E+00
       w_qmin( 247)=0.392961024E+00
       qk_min( 248)=0.909093140E+00
       w_qmin( 248)=0.396389433E+00
       qk_min( 249)=0.912079376E+00
       w_qmin( 249)=0.399767068E+00
       qk_min( 250)=0.915065612E+00
       w_qmin( 250)=0.403093691E+00
       qk_min( 251)=0.918051848E+00
       w_qmin( 251)=0.406369271E+00
       qk_min( 252)=0.921038084E+00
       w_qmin( 252)=0.409593802E+00
       qk_min( 253)=0.924024320E+00
       w_qmin( 253)=0.412767258E+00
       qk_min( 254)=0.927010556E+00
       w_qmin( 254)=0.415889743E+00
       qk_min( 255)=0.929996792E+00
       w_qmin( 255)=0.418961428E+00
       qk_min( 256)=0.932983028E+00
       w_qmin( 256)=0.421982410E+00
       qk_min( 257)=0.935969264E+00
       w_qmin( 257)=0.424952958E+00
       qk_min( 258)=0.938955500E+00
       w_qmin( 258)=0.427873267E+00
       qk_min( 259)=0.941941736E+00
       w_qmin( 259)=0.430743698E+00
       qk_min( 260)=0.944927972E+00
       w_qmin( 260)=0.433564510E+00
       qk_min( 261)=0.947914208E+00
       w_qmin( 261)=0.436336023E+00
       qk_min( 262)=0.950900444E+00
       w_qmin( 262)=0.439058676E+00
       qk_min( 263)=0.953886680E+00
       w_qmin( 263)=0.441732832E+00
       qk_min( 264)=0.956872916E+00
       w_qmin( 264)=0.444358892E+00
       qk_min( 265)=0.959859152E+00
       w_qmin( 265)=0.446937204E+00
       qk_min( 266)=0.962845388E+00
       w_qmin( 266)=0.449468407E+00
       qk_min( 267)=0.965831624E+00
       w_qmin( 267)=0.451952818E+00
       qk_min( 268)=0.968817860E+00
       w_qmin( 268)=0.454390902E+00
       qk_min( 269)=0.971804096E+00
       w_qmin( 269)=0.456783271E+00
       qk_min( 270)=0.974790332E+00
       w_qmin( 270)=0.459130352E+00
       qk_min( 271)=0.977776568E+00
       w_qmin( 271)=0.461432624E+00
       qk_min( 272)=0.980762804E+00
       w_qmin( 272)=0.463690551E+00
       qk_min( 273)=0.983749040E+00
       w_qmin( 273)=0.465904807E+00
       qk_min( 274)=0.986735276E+00
       w_qmin( 274)=0.468075779E+00
       qk_min( 275)=0.989721512E+00
       w_qmin( 275)=0.470204073E+00
       qk_min( 276)=0.992707748E+00
       w_qmin( 276)=0.472290120E+00
       qk_min( 277)=0.995693984E+00
       w_qmin( 277)=0.474334590E+00
       qk_min( 278)=0.998680219E+00
       w_qmin( 278)=0.476337904E+00
       qk_min( 279)=0.100166646E+01
       w_qmin( 279)=0.478300583E+00
       qk_min( 280)=0.100465269E+01
       w_qmin( 280)=0.480223343E+00
       qk_min( 281)=0.100763893E+01
       w_qmin( 281)=0.482106413E+00
       qk_min( 282)=0.101062516E+01
       w_qmin( 282)=0.483950519E+00
       qk_min( 283)=0.101361140E+01
       w_qmin( 283)=0.485756267E+00
       qk_min( 284)=0.101659764E+01
       w_qmin( 284)=0.487523890E+00
       qk_min( 285)=0.101958387E+01
       w_qmin( 285)=0.489254119E+00
       qk_min( 286)=0.102257011E+01
       w_qmin( 286)=0.490947531E+00
       qk_min( 287)=0.102555634E+01
       w_qmin( 287)=0.492604494E+00
       qk_min( 288)=0.102854258E+01
       w_qmin( 288)=0.494225644E+00
       qk_min( 289)=0.103152882E+01
       w_qmin( 289)=0.495811230E+00
       qk_min( 290)=0.103451505E+01
       w_qmin( 290)=0.497361969E+00
       qk_min( 291)=0.103750129E+01
       w_qmin( 291)=0.498878451E+00
       qk_min( 292)=0.104048752E+01
       w_qmin( 292)=0.500360890E+00
       qk_min( 293)=0.104347376E+01
       w_qmin( 293)=0.501809882E+00
       qk_min( 294)=0.104645999E+01
       w_qmin( 294)=0.503226110E+00
       qk_min( 295)=0.104944623E+01
       w_qmin( 295)=0.504609696E+00
       qk_min( 296)=0.105243247E+01
       w_qmin( 296)=0.505961288E+00
       qk_min( 297)=0.105541870E+01
       w_qmin( 297)=0.507281527E+00
       qk_min( 298)=0.105840494E+01
       w_qmin( 298)=0.508570536E+00
       qk_min( 299)=0.106139117E+01
       w_qmin( 299)=0.509828925E+00
       qk_min( 300)=0.106437741E+01
       w_qmin( 300)=0.511057261E+00
       qk_min( 301)=0.106736365E+01
       w_qmin( 301)=0.512255761E+00
       qk_min( 302)=0.107034988E+01
       w_qmin( 302)=0.513425076E+00
       qk_min( 303)=0.107333612E+01
       w_qmin( 303)=0.514565403E+00
       qk_min( 304)=0.107632235E+01
       w_qmin( 304)=0.515677265E+00
       qk_min( 305)=0.107930859E+01
       w_qmin( 305)=0.516761289E+00
       qk_min( 306)=0.108229483E+01
       w_qmin( 306)=0.517817503E+00
       qk_min( 307)=0.108528106E+01
       w_qmin( 307)=0.518846546E+00
       qk_min( 308)=0.108826730E+01
       w_qmin( 308)=0.519848956E+00
       qk_min( 309)=0.109125353E+01
       w_qmin( 309)=0.520824818E+00
       qk_min( 310)=0.109423977E+01
       w_qmin( 310)=0.521774666E+00
       qk_min( 311)=0.109722601E+01
       w_qmin( 311)=0.522699091E+00
       qk_min( 312)=0.110021224E+01
       w_qmin( 312)=0.523598083E+00
       qk_min( 313)=0.110319848E+01
       w_qmin( 313)=0.524472371E+00
       qk_min( 314)=0.110618471E+01
       w_qmin( 314)=0.525322030E+00
       qk_min( 315)=0.110917095E+01
       w_qmin( 315)=0.526147561E+00
       qk_min( 316)=0.111215719E+01
       w_qmin( 316)=0.526949538E+00
       qk_min( 317)=0.111514342E+01
       w_qmin( 317)=0.527727861E+00
       qk_min( 318)=0.111812966E+01
       w_qmin( 318)=0.528483162E+00
       qk_min( 319)=0.112111589E+01
       w_qmin( 319)=0.529215940E+00
       qk_min( 320)=0.112410213E+01
       w_qmin( 320)=0.529926124E+00
       qk_min( 321)=0.112708837E+01
       w_qmin( 321)=0.530614345E+00
       qk_min( 322)=0.113007460E+01
       w_qmin( 322)=0.531281004E+00
       qk_min( 323)=0.113306084E+01
       w_qmin( 323)=0.531926110E+00
       qk_min( 324)=0.113604707E+01
       w_qmin( 324)=0.532550319E+00
       qk_min( 325)=0.113903331E+01
       w_qmin( 325)=0.533153644E+00
       qk_min( 326)=0.114201955E+01
       w_qmin( 326)=0.533736493E+00
       qk_min( 327)=0.114500578E+01
       w_qmin( 327)=0.534299452E+00
       qk_min( 328)=0.114799202E+01
       w_qmin( 328)=0.534842376E+00
       qk_min( 329)=0.115097825E+01
       w_qmin( 329)=0.535365735E+00
       qk_min( 330)=0.115396449E+01
       w_qmin( 330)=0.535870152E+00
       qk_min( 331)=0.115695073E+01
       w_qmin( 331)=0.536355393E+00
       qk_min( 332)=0.115993696E+01
       w_qmin( 332)=0.536821972E+00
       qk_min( 333)=0.116292320E+01
       w_qmin( 333)=0.537270389E+00
       qk_min( 334)=0.116590943E+01
       w_qmin( 334)=0.537700432E+00
       qk_min( 335)=0.116889567E+01
       w_qmin( 335)=0.538112944E+00
       qk_min( 336)=0.117188191E+01
       w_qmin( 336)=0.538507669E+00
       qk_min( 337)=0.117486814E+01
       w_qmin( 337)=0.538885100E+00
       qk_min( 338)=0.117785438E+01
       w_qmin( 338)=0.539245715E+00
       qk_min( 339)=0.118084061E+01
       w_qmin( 339)=0.539589356E+00
       qk_min( 340)=0.118382685E+01
       w_qmin( 340)=0.539916440E+00
       qk_min( 341)=0.118681309E+01
       w_qmin( 341)=0.540227405E+00
       qk_min( 342)=0.118979932E+01
       w_qmin( 342)=0.540522272E+00
       qk_min( 343)=0.119278556E+01
       w_qmin( 343)=0.540801246E+00
       qk_min( 344)=0.119577179E+01
       w_qmin( 344)=0.541064908E+00
       qk_min( 345)=0.119875803E+01
       w_qmin( 345)=0.541313059E+00
       qk_min( 346)=0.120174427E+01
       w_qmin( 346)=0.541546322E+00
       qk_min( 347)=0.120473050E+01
       w_qmin( 347)=0.541764479E+00
       qk_min( 348)=0.120771674E+01
       w_qmin( 348)=0.541968063E+00
       qk_min( 349)=0.121070297E+01
       w_qmin( 349)=0.542157403E+00
       qk_min( 350)=0.121368921E+01
       w_qmin( 350)=0.542332339E+00
       qk_min( 351)=0.121667545E+01
       w_qmin( 351)=0.542493333E+00
       qk_min( 352)=0.121966168E+01
       w_qmin( 352)=0.542640711E+00
       qk_min( 353)=0.122264792E+01
       w_qmin( 353)=0.542774338E+00
       qk_min( 354)=0.122563415E+01
       w_qmin( 354)=0.542894618E+00
       qk_min( 355)=0.122862039E+01
       w_qmin( 355)=0.543001938E+00
       qk_min( 356)=0.123160662E+01
       w_qmin( 356)=0.543096128E+00
       qk_min( 357)=0.123459286E+01
       w_qmin( 357)=0.543177701E+00
       qk_min( 358)=0.123757910E+01
       w_qmin( 358)=0.543246522E+00
       qk_min( 359)=0.124056533E+01
       w_qmin( 359)=0.543302958E+00
       qk_min( 360)=0.124355157E+01
       w_qmin( 360)=0.543347446E+00
       qk_min( 361)=0.124653780E+01
       w_qmin( 361)=0.543379671E+00
       qk_min( 362)=0.124952404E+01
       w_qmin( 362)=0.543400084E+00
       qk_min( 363)=0.125251028E+01
       w_qmin( 363)=0.543409059E+00
       qk_min( 364)=0.125549651E+01
       w_qmin( 364)=0.543406400E+00
       qk_min( 365)=0.125848275E+01
       w_qmin( 365)=0.543392428E+00
       qk_min( 366)=0.126146898E+01
       w_qmin( 366)=0.543367507E+00
       qk_min( 367)=0.126445522E+01
       w_qmin( 367)=0.543331472E+00
       qk_min( 368)=0.126744146E+01
       w_qmin( 368)=0.543284860E+00
       qk_min( 369)=0.127042769E+01
       w_qmin( 369)=0.543227392E+00
       qk_min( 370)=0.127341393E+01
       w_qmin( 370)=0.543159509E+00
       qk_min( 371)=0.127640016E+01
       w_qmin( 371)=0.543081523E+00
       qk_min( 372)=0.127938640E+01
       w_qmin( 372)=0.542993193E+00
       qk_min( 373)=0.128237264E+01
       w_qmin( 373)=0.542894873E+00
       qk_min( 374)=0.128535887E+01
       w_qmin( 374)=0.542786982E+00
       qk_min( 375)=0.128834511E+01
       w_qmin( 375)=0.542669147E+00
       qk_min( 376)=0.129133134E+01
       w_qmin( 376)=0.542541860E+00
       qk_min( 377)=0.129431758E+01
       w_qmin( 377)=0.542405370E+00
       qk_min( 378)=0.129730382E+01
       w_qmin( 378)=0.542259464E+00
       qk_min( 379)=0.130029005E+01
       w_qmin( 379)=0.542104674E+00
       qk_min( 380)=0.130327629E+01
       w_qmin( 380)=0.541940738E+00
       qk_min( 381)=0.130626252E+01
       w_qmin( 381)=0.541767974E+00
       qk_min( 382)=0.130924876E+01
       w_qmin( 382)=0.541586815E+00
       qk_min( 383)=0.131223500E+01
       w_qmin( 383)=0.541396898E+00
       qk_min( 384)=0.131522123E+01
       w_qmin( 384)=0.541198618E+00
       qk_min( 385)=0.131820747E+01
       w_qmin( 385)=0.540992289E+00
       qk_min( 386)=0.132119370E+01
       w_qmin( 386)=0.540777604E+00
       qk_min( 387)=0.132417994E+01
       w_qmin( 387)=0.540554936E+00
       qk_min( 388)=0.132716618E+01
       w_qmin( 388)=0.540324661E+00
       qk_min( 389)=0.133015241E+01
       w_qmin( 389)=0.540086396E+00
       qk_min( 390)=0.133313865E+01
       w_qmin( 390)=0.539840736E+00
       qk_min( 391)=0.133612488E+01
       w_qmin( 391)=0.539587426E+00
       qk_min( 392)=0.133911112E+01
       w_qmin( 392)=0.539326728E+00
       qk_min( 393)=0.134209736E+01
       w_qmin( 393)=0.539058968E+00
       qk_min( 394)=0.134508359E+01
       w_qmin( 394)=0.538783923E+00
       qk_min( 395)=0.134806983E+01
       w_qmin( 395)=0.538501861E+00
       qk_min( 396)=0.135105606E+01
       w_qmin( 396)=0.538213123E+00
       qk_min( 397)=0.135404230E+01
       w_qmin( 397)=0.537917360E+00
       qk_min( 398)=0.135702854E+01
       w_qmin( 398)=0.537614984E+00
       qk_min( 399)=0.136001477E+01
       w_qmin( 399)=0.537306213E+00
       qk_min( 400)=0.136300101E+01
       w_qmin( 400)=0.536990875E+00
       qk_min( 401)=0.136598724E+01
       w_qmin( 401)=0.536669385E+00
       qk_min( 402)=0.136897348E+01
       w_qmin( 402)=0.536341436E+00
       qk_min( 403)=0.137195972E+01
       w_qmin( 403)=0.536007425E+00
       qk_min( 404)=0.137494595E+01
       w_qmin( 404)=0.535667621E+00
       qk_min( 405)=0.137793219E+01
       w_qmin( 405)=0.535321655E+00
       qk_min( 406)=0.138091842E+01
       w_qmin( 406)=0.534969909E+00
       qk_min( 407)=0.138390466E+01
       w_qmin( 407)=0.534612668E+00
       qk_min( 408)=0.138689090E+01
       w_qmin( 408)=0.534249689E+00
       qk_min( 409)=0.138987713E+01
       w_qmin( 409)=0.533881140E+00
       qk_min( 410)=0.139286337E+01
       w_qmin( 410)=0.533507407E+00
       qk_min( 411)=0.139584960E+01
       w_qmin( 411)=0.533128187E+00
       qk_min( 412)=0.139883584E+01
       w_qmin( 412)=0.532744012E+00
       qk_min( 413)=0.140182208E+01
       w_qmin( 413)=0.532354469E+00
       qk_min( 414)=0.140480831E+01
       w_qmin( 414)=0.531959918E+00
       qk_min( 415)=0.140779455E+01
       w_qmin( 415)=0.531560676E+00
       qk_min( 416)=0.141078078E+01
       w_qmin( 416)=0.531156400E+00
       qk_min( 417)=0.141376702E+01
       w_qmin( 417)=0.530747389E+00
       qk_min( 418)=0.141675325E+01
       w_qmin( 418)=0.530333888E+00
       qk_min( 419)=0.141973949E+01
       w_qmin( 419)=0.529915692E+00
       qk_min( 420)=0.142272573E+01
       w_qmin( 420)=0.529492975E+00
       qk_min( 421)=0.142571196E+01
       w_qmin( 421)=0.529066088E+00
       qk_min( 422)=0.142869820E+01
       w_qmin( 422)=0.528634723E+00
       qk_min( 423)=0.143168443E+01
       w_qmin( 423)=0.528199314E+00
       qk_min( 424)=0.143467067E+01
       w_qmin( 424)=0.527759608E+00
       qk_min( 425)=0.143765691E+01
       w_qmin( 425)=0.527315821E+00
       qk_min( 426)=0.144064314E+01
       w_qmin( 426)=0.526868363E+00
       qk_min( 427)=0.144362938E+01
       w_qmin( 427)=0.526416729E+00
       qk_min( 428)=0.144661561E+01
       w_qmin( 428)=0.525961361E+00
       qk_min( 429)=0.144960185E+01
       w_qmin( 429)=0.525502422E+00
       qk_min( 430)=0.145258809E+01
       w_qmin( 430)=0.525039686E+00
       qk_min( 431)=0.145557432E+01
       w_qmin( 431)=0.524573407E+00
       qk_min( 432)=0.145856056E+01
       w_qmin( 432)=0.524103775E+00
       qk_min( 433)=0.146154679E+01
       w_qmin( 433)=0.523630614E+00
       qk_min( 434)=0.146453303E+01
       w_qmin( 434)=0.523154084E+00
       qk_min( 435)=0.146751927E+01
       w_qmin( 435)=0.522674555E+00
       qk_min( 436)=0.147050550E+01
       w_qmin( 436)=0.522191583E+00
       qk_min( 437)=0.147349174E+01
       w_qmin( 437)=0.521705726E+00
       qk_min( 438)=0.147647797E+01
       w_qmin( 438)=0.521216638E+00
       qk_min( 439)=0.147946421E+01
       w_qmin( 439)=0.520724590E+00
       qk_min( 440)=0.148245045E+01
       w_qmin( 440)=0.520229833E+00
       qk_min( 441)=0.148543668E+01
       w_qmin( 441)=0.519732090E+00
       qk_min( 442)=0.148842292E+01
       w_qmin( 442)=0.519231598E+00
       qk_min( 443)=0.149140915E+01
       w_qmin( 443)=0.518728584E+00
       qk_min( 444)=0.149439539E+01
       w_qmin( 444)=0.518222796E+00
       qk_min( 445)=0.149738163E+01
       w_qmin( 445)=0.517714472E+00
       qk_min( 446)=0.150036786E+01
       w_qmin( 446)=0.517203847E+00
       qk_min( 447)=0.150335410E+01
       w_qmin( 447)=0.516690634E+00
       qk_min( 448)=0.150634033E+01
       w_qmin( 448)=0.516175256E+00
       qk_min( 449)=0.150932657E+01
       w_qmin( 449)=0.515657396E+00
       qk_min( 450)=0.151231281E+01
       w_qmin( 450)=0.515137294E+00
       qk_min( 451)=0.151529904E+01
       w_qmin( 451)=0.514615279E+00
       qk_min( 452)=0.151828528E+01
       w_qmin( 452)=0.514091024E+00
       qk_min( 453)=0.152127151E+01
       w_qmin( 453)=0.513564681E+00
       qk_min( 454)=0.152425775E+01
       w_qmin( 454)=0.513036626E+00
       qk_min( 455)=0.152724399E+01
       w_qmin( 455)=0.512506443E+00
       qk_min( 456)=0.153023022E+01
       w_qmin( 456)=0.511974495E+00
       qk_min( 457)=0.153321646E+01
       w_qmin( 457)=0.511440948E+00
       qk_min( 458)=0.153620269E+01
       w_qmin( 458)=0.510905468E+00
       qk_min( 459)=0.153918893E+01
       w_qmin( 459)=0.510368563E+00
       qk_min( 460)=0.154217517E+01
       w_qmin( 460)=0.509829851E+00
       qk_min( 461)=0.154516140E+01
       w_qmin( 461)=0.509289616E+00
       qk_min( 462)=0.154814764E+01
       w_qmin( 462)=0.508748070E+00
       qk_min( 463)=0.155113387E+01
       w_qmin( 463)=0.508204964E+00
       qk_min( 464)=0.155412011E+01
       w_qmin( 464)=0.507660398E+00
       qk_min( 465)=0.155710635E+01
       w_qmin( 465)=0.507114773E+00
       qk_min( 466)=0.156009258E+01
       w_qmin( 466)=0.506567686E+00
       qk_min( 467)=0.156307882E+01
       w_qmin( 467)=0.506019400E+00
       qk_min( 468)=0.156606505E+01
       w_qmin( 468)=0.505470120E+00
       qk_min( 469)=0.156905129E+01
       w_qmin( 469)=0.504919604E+00
       qk_min( 470)=0.157203753E+01
       w_qmin( 470)=0.504368231E+00
       qk_min( 471)=0.157502376E+01
       w_qmin( 471)=0.503815639E+00
       qk_min( 472)=0.157801000E+01
       w_qmin( 472)=0.503262117E+00
       qk_min( 473)=0.158099623E+01
       w_qmin( 473)=0.502707951E+00
       qk_min( 474)=0.158398247E+01
       w_qmin( 474)=0.502152789E+00
       qk_min( 475)=0.158696871E+01
       w_qmin( 475)=0.501596826E+00
       qk_min( 476)=0.158995494E+01
       w_qmin( 476)=0.501040318E+00
       qk_min( 477)=0.159294118E+01
       w_qmin( 477)=0.500482902E+00
       qk_min( 478)=0.159592741E+01
       w_qmin( 478)=0.499924909E+00
       qk_min( 479)=0.159891365E+01
       w_qmin( 479)=0.499366604E+00
       qk_min( 480)=0.160189988E+01
       w_qmin( 480)=0.498807513E+00
       qk_min( 481)=0.160488612E+01
       w_qmin( 481)=0.498248118E+00
       qk_min( 482)=0.160787236E+01
       w_qmin( 482)=0.497688151E+00
       qk_min( 483)=0.161085859E+01
       w_qmin( 483)=0.497127798E+00
       qk_min( 484)=0.161384483E+01
       w_qmin( 484)=0.496567364E+00
       qk_min( 485)=0.161683106E+01
       w_qmin( 485)=0.496006401E+00
       qk_min( 486)=0.161981730E+01
       w_qmin( 486)=0.495445240E+00
       qk_min( 487)=0.162280354E+01
       w_qmin( 487)=0.494884131E+00
       qk_min( 488)=0.162578977E+01
       w_qmin( 488)=0.494322710E+00
       qk_min( 489)=0.162877601E+01
       w_qmin( 489)=0.493761162E+00
       qk_min( 490)=0.163176224E+01
       w_qmin( 490)=0.493199736E+00
       qk_min( 491)=0.163474848E+01
       w_qmin( 491)=0.492638233E+00
       qk_min( 492)=0.163773472E+01
       w_qmin( 492)=0.492076930E+00
       qk_min( 493)=0.164072095E+01
       w_qmin( 493)=0.491515523E+00
       qk_min( 494)=0.164370719E+01
       w_qmin( 494)=0.490954318E+00
       qk_min( 495)=0.164669342E+01
       w_qmin( 495)=0.490393470E+00
       qk_min( 496)=0.164967966E+01
       w_qmin( 496)=0.489832666E+00
       qk_min( 497)=0.165266590E+01
       w_qmin( 497)=0.489272158E+00
       qk_min( 498)=0.165565213E+01
       w_qmin( 498)=0.488712218E+00
       qk_min( 499)=0.165863837E+01
       w_qmin( 499)=0.488152449E+00
       qk_min( 500)=0.166162460E+01
       w_qmin( 500)=0.487593148E+00
       qk_min( 501)=0.166461084E+01
       w_qmin( 501)=0.487034420E+00
       qk_min( 502)=0.166759708E+01
       w_qmin( 502)=0.486476058E+00
       qk_min( 503)=0.167058331E+01
       w_qmin( 503)=0.485918525E+00
       qk_min( 504)=0.167356955E+01
       w_qmin( 504)=0.485361316E+00
       qk_min( 505)=0.167655578E+01
       w_qmin( 505)=0.484804867E+00
       qk_min( 506)=0.167954202E+01
       w_qmin( 506)=0.484249219E+00
       qk_min( 507)=0.168252826E+01
       w_qmin( 507)=0.483694177E+00
       qk_min( 508)=0.168551449E+01
       w_qmin( 508)=0.483139964E+00
       qk_min( 509)=0.168850073E+01
       w_qmin( 509)=0.482586753E+00
       qk_min( 510)=0.169148696E+01
       w_qmin( 510)=0.482034212E+00
       qk_min( 511)=0.169447320E+01
       w_qmin( 511)=0.481482667E+00
       qk_min( 512)=0.169745944E+01
       w_qmin( 512)=0.480932319E+00
       qk_min( 513)=0.170044567E+01
       w_qmin( 513)=0.480382781E+00
       qk_min( 514)=0.170343191E+01
       w_qmin( 514)=0.479834505E+00
       qk_min( 515)=0.170641814E+01
       w_qmin( 515)=0.479287223E+00
       qk_min( 516)=0.170940438E+01
       w_qmin( 516)=0.478741097E+00
       qk_min( 517)=0.171239062E+01
       w_qmin( 517)=0.478196318E+00
       qk_min( 518)=0.171537685E+01
       w_qmin( 518)=0.477652681E+00
       qk_min( 519)=0.171836309E+01
       w_qmin( 519)=0.477110392E+00
       qk_min( 520)=0.172134932E+01
       w_qmin( 520)=0.476569620E+00
       qk_min( 521)=0.172433556E+01
       w_qmin( 521)=0.476030097E+00
       qk_min( 522)=0.172732180E+01
       w_qmin( 522)=0.475492019E+00
       qk_min( 523)=0.173030803E+01
       w_qmin( 523)=0.474955653E+00
       qk_min( 524)=0.173329427E+01
       w_qmin( 524)=0.474420744E+00
       qk_min( 525)=0.173628050E+01
       w_qmin( 525)=0.473887491E+00
       qk_min( 526)=0.173926674E+01
       w_qmin( 526)=0.473355784E+00
       qk_min( 527)=0.174225298E+01
       w_qmin( 527)=0.472825855E+00
       qk_min( 528)=0.174523921E+01
       w_qmin( 528)=0.472297828E+00
       qk_min( 529)=0.174822545E+01
       w_qmin( 529)=0.471771411E+00
       qk_min( 530)=0.175121168E+01
       w_qmin( 530)=0.471246842E+00
       qk_min( 531)=0.175419792E+01
       w_qmin( 531)=0.470724461E+00
       qk_min( 532)=0.175718416E+01
       w_qmin( 532)=0.470203938E+00
       qk_min( 533)=0.176017039E+01
       w_qmin( 533)=0.469685388E+00
       qk_min( 534)=0.176315663E+01
       w_qmin( 534)=0.469169078E+00
       qk_min( 535)=0.176614286E+01
       w_qmin( 535)=0.468654804E+00
       qk_min( 536)=0.176912910E+01
       w_qmin( 536)=0.468142913E+00
       qk_min( 537)=0.177211533E+01
       w_qmin( 537)=0.467633090E+00
       qk_min( 538)=0.177510157E+01
       w_qmin( 538)=0.467125676E+00
       qk_min( 539)=0.177808781E+01
       w_qmin( 539)=0.466620768E+00
       qk_min( 540)=0.178107404E+01
       w_qmin( 540)=0.466118218E+00
       qk_min( 541)=0.178406028E+01
       w_qmin( 541)=0.465618229E+00
       qk_min( 542)=0.178704651E+01
       w_qmin( 542)=0.465120957E+00
       qk_min( 543)=0.179003275E+01
       w_qmin( 543)=0.464626218E+00
       qk_min( 544)=0.179301899E+01
       w_qmin( 544)=0.464134308E+00
       qk_min( 545)=0.179600522E+01
       w_qmin( 545)=0.463645353E+00
       qk_min( 546)=0.179899146E+01
       w_qmin( 546)=0.463159109E+00
       qk_min( 547)=0.180197769E+01
       w_qmin( 547)=0.462676061E+00
       qk_min( 548)=0.180496393E+01
       w_qmin( 548)=0.462195894E+00
       qk_min( 549)=0.180795017E+01
       w_qmin( 549)=0.461718899E+00
       qk_min( 550)=0.181093640E+01
       w_qmin( 550)=0.461245308E+00
       qk_min( 551)=0.181392264E+01
       w_qmin( 551)=0.460774970E+00
       qk_min( 552)=0.181690887E+01
       w_qmin( 552)=0.460308043E+00
       qk_min( 553)=0.181989511E+01
       w_qmin( 553)=0.459844769E+00
       qk_min( 554)=0.182288135E+01
       w_qmin( 554)=0.459385040E+00
       qk_min( 555)=0.182586758E+01
       w_qmin( 555)=0.458928991E+00
       qk_min( 556)=0.182885382E+01
       w_qmin( 556)=0.458477051E+00
       qk_min( 557)=0.183184005E+01
       w_qmin( 557)=0.458028918E+00
       qk_min( 558)=0.183482629E+01
       w_qmin( 558)=0.457584775E+00
       qk_min( 559)=0.183781253E+01
       w_qmin( 559)=0.457145134E+00
       qk_min( 560)=0.184079876E+01
       w_qmin( 560)=0.456709548E+00
       qk_min( 561)=0.184378500E+01
       w_qmin( 561)=0.456278691E+00
       qk_min( 562)=0.184677123E+01
       w_qmin( 562)=0.455852184E+00
       qk_min( 563)=0.184975747E+01
       w_qmin( 563)=0.455430482E+00
       qk_min( 564)=0.185274371E+01
       w_qmin( 564)=0.455013800E+00
       qk_min( 565)=0.185572994E+01
       w_qmin( 565)=0.454601965E+00
       qk_min( 566)=0.185871618E+01
       w_qmin( 566)=0.454195356E+00
       qk_min( 567)=0.186170241E+01
       w_qmin( 567)=0.453794181E+00
       qk_min( 568)=0.186468865E+01
       w_qmin( 568)=0.453398279E+00
       qk_min( 569)=0.186767489E+01
       w_qmin( 569)=0.453007960E+00
       qk_min( 570)=0.187066112E+01
       w_qmin( 570)=0.452623488E+00
       qk_min( 571)=0.187364736E+01
       w_qmin( 571)=0.452244727E+00
       qk_min( 572)=0.187663359E+01
       w_qmin( 572)=0.451872140E+00
       qk_min( 573)=0.187961983E+01
       w_qmin( 573)=0.451505427E+00
       qk_min( 574)=0.188260607E+01
       w_qmin( 574)=0.451144991E+00
       qk_min( 575)=0.188559230E+01
       w_qmin( 575)=0.450790965E+00
       qk_min( 576)=0.188857854E+01
       w_qmin( 576)=0.450443152E+00
       qk_min( 577)=0.189156477E+01
       w_qmin( 577)=0.450101827E+00
       qk_min( 578)=0.189455101E+01
       w_qmin( 578)=0.449767107E+00
       qk_min( 579)=0.189753725E+01
       w_qmin( 579)=0.449438788E+00
       qk_min( 580)=0.190052348E+01
       w_qmin( 580)=0.449116874E+00
       qk_min( 581)=0.190350972E+01
       w_qmin( 581)=0.448801707E+00
       qk_min( 582)=0.190649595E+01
       w_qmin( 582)=0.448492636E+00
       qk_min( 583)=0.190948219E+01
       w_qmin( 583)=0.448190064E+00
       qk_min( 584)=0.191246843E+01
       w_qmin( 584)=0.447893625E+00
       qk_min( 585)=0.191545466E+01
       w_qmin( 585)=0.447603192E+00
       qk_min( 586)=0.191844090E+01
       w_qmin( 586)=0.447318726E+00
       qk_min( 587)=0.192142713E+01
       w_qmin( 587)=0.447039829E+00
       qk_min( 588)=0.192441337E+01
       w_qmin( 588)=0.446766438E+00
       qk_min( 589)=0.192739961E+01
       w_qmin( 589)=0.446498400E+00
       qk_min( 590)=0.193038584E+01
       w_qmin( 590)=0.446235093E+00
       qk_min( 591)=0.193337208E+01
       w_qmin( 591)=0.445976490E+00
       qk_min( 592)=0.193635831E+01
       w_qmin( 592)=0.445722525E+00
       qk_min( 593)=0.193934455E+01
       w_qmin( 593)=0.445472463E+00
       qk_min( 594)=0.194233079E+01
       w_qmin( 594)=0.445226434E+00
       qk_min( 595)=0.194531702E+01
       w_qmin( 595)=0.444983799E+00
       qk_min( 596)=0.194830326E+01
       w_qmin( 596)=0.444744462E+00
       qk_min( 597)=0.195128949E+01
       w_qmin( 597)=0.444508317E+00
       qk_min( 598)=0.195427573E+01
       w_qmin( 598)=0.444274700E+00
       qk_min( 599)=0.195726196E+01
       w_qmin( 599)=0.444043671E+00
       qk_min( 600)=0.196024820E+01
       w_qmin( 600)=0.443815056E+00
       qk_min( 601)=0.196323444E+01
       w_qmin( 601)=0.443588401E+00
       qk_min( 602)=0.196622067E+01
       w_qmin( 602)=0.443363527E+00
       qk_min( 603)=0.196920691E+01
       w_qmin( 603)=0.443140523E+00
       qk_min( 604)=0.197219314E+01
       w_qmin( 604)=0.442918871E+00
       qk_min( 605)=0.197517938E+01
       w_qmin( 605)=0.442698666E+00
       qk_min( 606)=0.197816562E+01
       w_qmin( 606)=0.442479554E+00
       qk_min( 607)=0.198115185E+01
       w_qmin( 607)=0.442261491E+00
       qk_min( 608)=0.198413809E+01
       w_qmin( 608)=0.442044628E+00
       qk_min( 609)=0.198712432E+01
       w_qmin( 609)=0.441828483E+00
       qk_min( 610)=0.199011056E+01
       w_qmin( 610)=0.441613110E+00
       qk_min( 611)=0.199309680E+01
       w_qmin( 611)=0.441398474E+00
       qk_min( 612)=0.199608303E+01
       w_qmin( 612)=0.441184405E+00
       qk_min( 613)=0.199906927E+01
       w_qmin( 613)=0.440970956E+00
       qk_min( 614)=0.200205550E+01
       w_qmin( 614)=0.440758017E+00
       qk_min( 615)=0.200504174E+01
       w_qmin( 615)=0.440545795E+00
       qk_min( 616)=0.200802798E+01
       w_qmin( 616)=0.440333851E+00
       qk_min( 617)=0.201101421E+01
       w_qmin( 617)=0.440122243E+00
       qk_min( 618)=0.201400045E+01
       w_qmin( 618)=0.439911113E+00
       qk_min( 619)=0.201698668E+01
       w_qmin( 619)=0.439700340E+00
       qk_min( 620)=0.201997292E+01
       w_qmin( 620)=0.439490007E+00
       qk_min( 621)=0.202295916E+01
       w_qmin( 621)=0.439280374E+00
       qk_min( 622)=0.202594539E+01
       w_qmin( 622)=0.439070781E+00
       qk_min( 623)=0.202893163E+01
       w_qmin( 623)=0.438861636E+00
       qk_min( 624)=0.203191786E+01
       w_qmin( 624)=0.438652871E+00
       qk_min( 625)=0.203490410E+01
       w_qmin( 625)=0.438444506E+00
       qk_min( 626)=0.203789034E+01
       w_qmin( 626)=0.438236817E+00
       qk_min( 627)=0.204087657E+01
       w_qmin( 627)=0.438029311E+00
       qk_min( 628)=0.204386281E+01
       w_qmin( 628)=0.437822142E+00
       qk_min( 629)=0.204684904E+01
       w_qmin( 629)=0.437615533E+00
       qk_min( 630)=0.204983528E+01
       w_qmin( 630)=0.437409289E+00
       qk_min( 631)=0.205282152E+01
       w_qmin( 631)=0.437203640E+00
       qk_min( 632)=0.205580775E+01
       w_qmin( 632)=0.436998714E+00
       qk_min( 633)=0.205879399E+01
       w_qmin( 633)=0.436794000E+00
       qk_min( 634)=0.206178022E+01
       w_qmin( 634)=0.436589930E+00
       qk_min( 635)=0.206476646E+01
       w_qmin( 635)=0.436386245E+00
       qk_min( 636)=0.206775270E+01
       w_qmin( 636)=0.436183265E+00
       qk_min( 637)=0.207073893E+01
       w_qmin( 637)=0.435981103E+00
       qk_min( 638)=0.207372517E+01
       w_qmin( 638)=0.435779330E+00
       qk_min( 639)=0.207671140E+01
       w_qmin( 639)=0.435578204E+00
       qk_min( 640)=0.207969764E+01
       w_qmin( 640)=0.435377642E+00
       qk_min( 641)=0.208268388E+01
       w_qmin( 641)=0.435177798E+00
       qk_min( 642)=0.208567011E+01
       w_qmin( 642)=0.434978644E+00
       qk_min( 643)=0.208865635E+01
       w_qmin( 643)=0.434780538E+00
       qk_min( 644)=0.209164258E+01
       w_qmin( 644)=0.434582760E+00
       qk_min( 645)=0.209462882E+01
       w_qmin( 645)=0.434385882E+00
       qk_min( 646)=0.209761506E+01
       w_qmin( 646)=0.434189713E+00
       qk_min( 647)=0.210060129E+01
       w_qmin( 647)=0.433994294E+00
       qk_min( 648)=0.210358753E+01
       w_qmin( 648)=0.433799996E+00
       qk_min( 649)=0.210657376E+01
       w_qmin( 649)=0.433606223E+00
       qk_min( 650)=0.210956000E+01
       w_qmin( 650)=0.433413285E+00
       qk_min( 651)=0.211254624E+01
       w_qmin( 651)=0.433221235E+00
       qk_min( 652)=0.211553247E+01
       w_qmin( 652)=0.433030134E+00
       qk_min( 653)=0.211851871E+01
       w_qmin( 653)=0.432839735E+00
       qk_min( 654)=0.212150494E+01
       w_qmin( 654)=0.432650622E+00
       qk_min( 655)=0.212449118E+01
       w_qmin( 655)=0.432462082E+00
       qk_min( 656)=0.212747742E+01
       w_qmin( 656)=0.432274568E+00
       qk_min( 657)=0.213046365E+01
       w_qmin( 657)=0.432087939E+00
       qk_min( 658)=0.213344989E+01
       w_qmin( 658)=0.431902265E+00
       qk_min( 659)=0.213643612E+01
       w_qmin( 659)=0.431717845E+00
       qk_min( 660)=0.213942236E+01
       w_qmin( 660)=0.431534055E+00
       qk_min( 661)=0.214240859E+01
       w_qmin( 661)=0.431351397E+00
       qk_min( 662)=0.214539483E+01
       w_qmin( 662)=0.431169650E+00
       qk_min( 663)=0.214838107E+01
       w_qmin( 663)=0.430988927E+00
       qk_min( 664)=0.215136730E+01
       w_qmin( 664)=0.430809172E+00
       qk_min( 665)=0.215435354E+01
       w_qmin( 665)=0.430630735E+00
       qk_min( 666)=0.215733977E+01
       w_qmin( 666)=0.430453213E+00
       qk_min( 667)=0.216032601E+01
       w_qmin( 667)=0.430276570E+00
       qk_min( 668)=0.216331225E+01
       w_qmin( 668)=0.430100990E+00
       qk_min( 669)=0.216629848E+01
       w_qmin( 669)=0.429926630E+00
       qk_min( 670)=0.216928472E+01
       w_qmin( 670)=0.429753425E+00
       qk_min( 671)=0.217227095E+01
       w_qmin( 671)=0.429581112E+00
       qk_min( 672)=0.217525719E+01
       w_qmin( 672)=0.429409943E+00
       qk_min( 673)=0.217824343E+01
       w_qmin( 673)=0.429239857E+00
       qk_min( 674)=0.218122966E+01
       w_qmin( 674)=0.429070807E+00
       qk_min( 675)=0.218421590E+01
       w_qmin( 675)=0.428902973E+00
       qk_min( 676)=0.218720213E+01
       w_qmin( 676)=0.428736351E+00
       qk_min( 677)=0.219018837E+01
       w_qmin( 677)=0.428570846E+00
       qk_min( 678)=0.219317461E+01
       w_qmin( 678)=0.428406256E+00
       qk_min( 679)=0.219616084E+01
       w_qmin( 679)=0.428242836E+00
       qk_min( 680)=0.219914708E+01
       w_qmin( 680)=0.428080605E+00
       qk_min( 681)=0.220213331E+01
       w_qmin( 681)=0.427919878E+00
       qk_min( 682)=0.220511955E+01
       w_qmin( 682)=0.427759946E+00
       qk_min( 683)=0.220810579E+01
       w_qmin( 683)=0.427601150E+00
       qk_min( 684)=0.221109202E+01
       w_qmin( 684)=0.427443439E+00
       qk_min( 685)=0.221407826E+01
       w_qmin( 685)=0.427287027E+00
       qk_min( 686)=0.221706449E+01
       w_qmin( 686)=0.427131635E+00
       qk_min( 687)=0.222005073E+01
       w_qmin( 687)=0.426977770E+00
       qk_min( 688)=0.222303697E+01
       w_qmin( 688)=0.426824760E+00
       qk_min( 689)=0.222602320E+01
       w_qmin( 689)=0.426672886E+00
       qk_min( 690)=0.222900944E+01
       w_qmin( 690)=0.426522195E+00
       qk_min( 691)=0.223199567E+01
       w_qmin( 691)=0.426372740E+00
       qk_min( 692)=0.223498191E+01
       w_qmin( 692)=0.426224617E+00
       qk_min( 693)=0.223796815E+01
       w_qmin( 693)=0.426077496E+00
       qk_min( 694)=0.224095438E+01
       w_qmin( 694)=0.425931487E+00
       qk_min( 695)=0.224394062E+01
       w_qmin( 695)=0.425786782E+00
       qk_min( 696)=0.224692685E+01
       w_qmin( 696)=0.425643151E+00
       qk_min( 697)=0.224991309E+01
       w_qmin( 697)=0.425500645E+00
       qk_min( 698)=0.225289933E+01
       w_qmin( 698)=0.425359677E+00
       qk_min( 699)=0.225588556E+01
       w_qmin( 699)=0.425219639E+00
       qk_min( 700)=0.225887180E+01
       w_qmin( 700)=0.425080800E+00
       qk_min( 701)=0.226185803E+01
       w_qmin( 701)=0.424943155E+00
       qk_min( 702)=0.226484427E+01
       w_qmin( 702)=0.424806756E+00
       qk_min( 703)=0.226783051E+01
       w_qmin( 703)=0.424671810E+00
       qk_min( 704)=0.227081674E+01
       w_qmin( 704)=0.424537798E+00
       qk_min( 705)=0.227380298E+01
       w_qmin( 705)=0.424404996E+00
       qk_min( 706)=0.227678921E+01
       w_qmin( 706)=0.424273362E+00
       qk_min( 707)=0.227977545E+01
       w_qmin( 707)=0.424143043E+00
       qk_min( 708)=0.228276169E+01
       w_qmin( 708)=0.424013824E+00
       qk_min( 709)=0.228574792E+01
       w_qmin( 709)=0.423886073E+00
       qk_min( 710)=0.228873416E+01
       w_qmin( 710)=0.423759328E+00
       qk_min( 711)=0.229172039E+01
       w_qmin( 711)=0.423633786E+00
       qk_min( 712)=0.229470663E+01
       w_qmin( 712)=0.423509341E+00
       qk_min( 713)=0.229769287E+01
       w_qmin( 713)=0.423386290E+00
       qk_min( 714)=0.230067910E+01
       w_qmin( 714)=0.423264588E+00
       qk_min( 715)=0.230366534E+01
       w_qmin( 715)=0.423143812E+00
       qk_min( 716)=0.230665157E+01
       w_qmin( 716)=0.423024215E+00
       qk_min( 717)=0.230963781E+01
       w_qmin( 717)=0.422905902E+00
       qk_min( 718)=0.231262405E+01
       w_qmin( 718)=0.422788711E+00
       qk_min( 719)=0.231561028E+01
       w_qmin( 719)=0.422672726E+00
       qk_min( 720)=0.231859652E+01
       w_qmin( 720)=0.422558089E+00
       qk_min( 721)=0.232158275E+01
       w_qmin( 721)=0.422444571E+00
       qk_min( 722)=0.232456899E+01
       w_qmin( 722)=0.422332205E+00
       qk_min( 723)=0.232755522E+01
       w_qmin( 723)=0.422220969E+00
       qk_min( 724)=0.233054146E+01
       w_qmin( 724)=0.422110965E+00
       qk_min( 725)=0.233352770E+01
       w_qmin( 725)=0.422002384E+00
       qk_min( 726)=0.233651393E+01
       w_qmin( 726)=0.421894768E+00
       qk_min( 727)=0.233950017E+01
       w_qmin( 727)=0.421788434E+00
       qk_min( 728)=0.234248640E+01
       w_qmin( 728)=0.421683337E+00
       qk_min( 729)=0.234547264E+01
       w_qmin( 729)=0.421579294E+00
       qk_min( 730)=0.234845888E+01
       w_qmin( 730)=0.421476620E+00
       qk_min( 731)=0.235144511E+01
       w_qmin( 731)=0.421375333E+00
       qk_min( 732)=0.235443135E+01
       w_qmin( 732)=0.421275019E+00
       qk_min( 733)=0.235741758E+01
       w_qmin( 733)=0.421176101E+00
       qk_min( 734)=0.236040382E+01
       w_qmin( 734)=0.421078203E+00
       qk_min( 735)=0.236339006E+01
       w_qmin( 735)=0.420981495E+00
       qk_min( 736)=0.236637629E+01
       w_qmin( 736)=0.420886382E+00
       qk_min( 737)=0.236936253E+01
       w_qmin( 737)=0.420792297E+00
       qk_min( 738)=0.237234876E+01
       w_qmin( 738)=0.420699280E+00
       qk_min( 739)=0.237533500E+01
       w_qmin( 739)=0.420607569E+00
       qk_min( 740)=0.237832124E+01
       w_qmin( 740)=0.420516928E+00
       qk_min( 741)=0.238130747E+01
       w_qmin( 741)=0.420427695E+00
       qk_min( 742)=0.238429371E+01
       w_qmin( 742)=0.420339885E+00
       qk_min( 743)=0.238727994E+01
       w_qmin( 743)=0.420252997E+00
       qk_min( 744)=0.239026618E+01
       w_qmin( 744)=0.420167183E+00
       qk_min( 745)=0.239325242E+01
       w_qmin( 745)=0.420082598E+00
       qk_min( 746)=0.239623865E+01
       w_qmin( 746)=0.419999125E+00
       qk_min( 747)=0.239922489E+01
       w_qmin( 747)=0.419917219E+00
       qk_min( 748)=0.240221112E+01
       w_qmin( 748)=0.419835982E+00
       qk_min( 749)=0.240519736E+01
       w_qmin( 749)=0.419756184E+00
       qk_min( 750)=0.240818360E+01
       w_qmin( 750)=0.419677361E+00
       qk_min( 751)=0.241116983E+01
       w_qmin( 751)=0.419599736E+00
       qk_min( 752)=0.241415607E+01
       w_qmin( 752)=0.419523267E+00
       qk_min( 753)=0.241714230E+01
       w_qmin( 753)=0.419448087E+00
       qk_min( 754)=0.242012854E+01
       w_qmin( 754)=0.419373828E+00
       qk_min( 755)=0.242311478E+01
       w_qmin( 755)=0.419300734E+00
       qk_min( 756)=0.242610101E+01
       w_qmin( 756)=0.419228732E+00
       qk_min( 757)=0.242908725E+01
       w_qmin( 757)=0.419157933E+00
       qk_min( 758)=0.243207348E+01
       w_qmin( 758)=0.419088375E+00
       qk_min( 759)=0.243505972E+01
       w_qmin( 759)=0.419019684E+00
       qk_min( 760)=0.243804596E+01
       w_qmin( 760)=0.418952174E+00
       qk_min( 761)=0.244103219E+01
       w_qmin( 761)=0.418885790E+00
       qk_min( 762)=0.244401843E+01
       w_qmin( 762)=0.418820479E+00
       qk_min( 763)=0.244700466E+01
       w_qmin( 763)=0.418756407E+00
       qk_min( 764)=0.244999090E+01
       w_qmin( 764)=0.418693566E+00
       qk_min( 765)=0.245297714E+01
       w_qmin( 765)=0.418631618E+00
       qk_min( 766)=0.245596337E+01
       w_qmin( 766)=0.418570856E+00
       qk_min( 767)=0.245894961E+01
       w_qmin( 767)=0.418511290E+00
       qk_min( 768)=0.246193584E+01
       w_qmin( 768)=0.418452610E+00
       qk_min( 769)=0.246492208E+01
       w_qmin( 769)=0.418395480E+00
       qk_min( 770)=0.246790832E+01
       w_qmin( 770)=0.418339268E+00
       qk_min( 771)=0.247089455E+01
       w_qmin( 771)=0.418284168E+00
       qk_min( 772)=0.247388079E+01
       w_qmin( 772)=0.418230108E+00
       qk_min( 773)=0.247686702E+01
       w_qmin( 773)=0.418177159E+00
       qk_min( 774)=0.247985326E+01
       w_qmin( 774)=0.418125396E+00
       qk_min( 775)=0.248283950E+01
       w_qmin( 775)=0.418074852E+00
       qk_min( 776)=0.248582573E+01
       w_qmin( 776)=0.418025387E+00
       qk_min( 777)=0.248881197E+01
       w_qmin( 777)=0.417976794E+00
       qk_min( 778)=0.249179820E+01
       w_qmin( 778)=0.417929350E+00
       qk_min( 779)=0.249478444E+01
       w_qmin( 779)=0.417883051E+00
       qk_min( 780)=0.249777068E+01
       w_qmin( 780)=0.417838039E+00
       qk_min( 781)=0.250075691E+01
       w_qmin( 781)=0.417793664E+00
       qk_min( 782)=0.250374315E+01
       w_qmin( 782)=0.417750561E+00
       qk_min( 783)=0.250672938E+01
       w_qmin( 783)=0.417708424E+00
       qk_min( 784)=0.250971562E+01
       w_qmin( 784)=0.417667126E+00
       qk_min( 785)=0.251270185E+01
       w_qmin( 785)=0.417627067E+00
       qk_min( 786)=0.251568809E+01
       w_qmin( 786)=0.417588083E+00
       qk_min( 787)=0.251867433E+01
       w_qmin( 787)=0.417549905E+00
       qk_min( 788)=0.252166056E+01
       w_qmin( 788)=0.417512662E+00
       qk_min( 789)=0.252464680E+01
       w_qmin( 789)=0.417476336E+00
       qk_min( 790)=0.252763303E+01
       w_qmin( 790)=0.417440995E+00
       qk_min( 791)=0.253061927E+01
       w_qmin( 791)=0.417406728E+00
       qk_min( 792)=0.253360551E+01
       w_qmin( 792)=0.417373123E+00
       qk_min( 793)=0.253659174E+01
       w_qmin( 793)=0.417340483E+00
       qk_min( 794)=0.253957798E+01
       w_qmin( 794)=0.417308720E+00
       qk_min( 795)=0.254256421E+01
       w_qmin( 795)=0.417277690E+00
       qk_min( 796)=0.254555045E+01
       w_qmin( 796)=0.417247521E+00
       qk_min( 797)=0.254853669E+01
       w_qmin( 797)=0.417218395E+00
       qk_min( 798)=0.255152292E+01
       w_qmin( 798)=0.417189876E+00
       qk_min( 799)=0.255450916E+01
       w_qmin( 799)=0.417162311E+00
       qk_min( 800)=0.255749539E+01
       w_qmin( 800)=0.417135566E+00
       qk_min( 801)=0.256048163E+01
       w_qmin( 801)=0.417109384E+00
       qk_min( 802)=0.256346787E+01
       w_qmin( 802)=0.417084383E+00
       qk_min( 803)=0.256645410E+01
       w_qmin( 803)=0.417059850E+00
       qk_min( 804)=0.256944034E+01
       w_qmin( 804)=0.417036166E+00
       qk_min( 805)=0.257242657E+01
       w_qmin( 805)=0.417013231E+00
       qk_min( 806)=0.257541281E+01
       w_qmin( 806)=0.416991167E+00
       qk_min( 807)=0.257839905E+01
       w_qmin( 807)=0.416969663E+00
       qk_min( 808)=0.258138528E+01
       w_qmin( 808)=0.416949225E+00
       qk_min( 809)=0.258437152E+01
       w_qmin( 809)=0.416929324E+00
       qk_min( 810)=0.258735775E+01
       w_qmin( 810)=0.416910287E+00
       qk_min( 811)=0.259034399E+01
       w_qmin( 811)=0.416891862E+00
       qk_min( 812)=0.259333023E+01
       w_qmin( 812)=0.416874231E+00
       qk_min( 813)=0.259631646E+01
       w_qmin( 813)=0.416857581E+00
       qk_min( 814)=0.259930270E+01
       w_qmin( 814)=0.416841411E+00
       qk_min( 815)=0.260228893E+01
       w_qmin( 815)=0.416826065E+00
       qk_min( 816)=0.260527517E+01
       w_qmin( 816)=0.416811573E+00
       qk_min( 817)=0.260826141E+01
       w_qmin( 817)=0.416797766E+00
       qk_min( 818)=0.261124764E+01
       w_qmin( 818)=0.416784748E+00
       qk_min( 819)=0.261423388E+01
       w_qmin( 819)=0.416772741E+00
       qk_min( 820)=0.261722011E+01
       w_qmin( 820)=0.416761080E+00
       qk_min( 821)=0.262020635E+01
       w_qmin( 821)=0.416750303E+00
       qk_min( 822)=0.262319259E+01
       w_qmin( 822)=0.416740399E+00
       qk_min( 823)=0.262617882E+01
       w_qmin( 823)=0.416731210E+00
       qk_min( 824)=0.262916506E+01
       w_qmin( 824)=0.416722860E+00
       qk_min( 825)=0.263215129E+01
       w_qmin( 825)=0.416715056E+00
       qk_min( 826)=0.263513753E+01
       w_qmin( 826)=0.416708245E+00
       qk_min( 827)=0.263812377E+01
       w_qmin( 827)=0.416701947E+00
       qk_min( 828)=0.264111000E+01
       w_qmin( 828)=0.416696476E+00
       qk_min( 829)=0.264409624E+01
       w_qmin( 829)=0.416691823E+00
       qk_min( 830)=0.264708247E+01
       w_qmin( 830)=0.416688036E+00
       qk_min( 831)=0.265006871E+01
       w_qmin( 831)=0.416684648E+00
       qk_min( 832)=0.265305495E+01
       w_qmin( 832)=0.416682209E+00
       qk_min( 833)=0.265604118E+01
       w_qmin( 833)=0.416680447E+00
       qk_min( 834)=0.265902742E+01
       w_qmin( 834)=0.416679214E+00
       qk_min( 835)=0.266201365E+01
       w_qmin( 835)=0.416679097E+00
       qk_min( 836)=0.266499989E+01
       w_qmin( 836)=0.416679330E+00
       qk_min( 837)=0.266798613E+01
       w_qmin( 837)=0.416680213E+00
       qk_min( 838)=0.267097236E+01
       w_qmin( 838)=0.416681818E+00
       qk_min( 839)=0.267395860E+01
       w_qmin( 839)=0.416684115E+00
       qk_min( 840)=0.267694483E+01
       w_qmin( 840)=0.416687081E+00
       qk_min( 841)=0.267993107E+01
       w_qmin( 841)=0.416690710E+00
       qk_min( 842)=0.268291730E+01
       w_qmin( 842)=0.416694802E+00
       qk_min( 843)=0.268590354E+01
       w_qmin( 843)=0.416699643E+00
       qk_min( 844)=0.268888978E+01
       w_qmin( 844)=0.416704850E+00
       qk_min( 845)=0.269187601E+01
       w_qmin( 845)=0.416710788E+00
       qk_min( 846)=0.269486225E+01
       w_qmin( 846)=0.416717432E+00
       qk_min( 847)=0.269784848E+01
       w_qmin( 847)=0.416724497E+00
       qk_min( 848)=0.270083472E+01
       w_qmin( 848)=0.416731992E+00
       qk_min( 849)=0.270382096E+01
       w_qmin( 849)=0.416740184E+00
       qk_min( 850)=0.270680719E+01
       w_qmin( 850)=0.416748775E+00
       qk_min( 851)=0.270979343E+01
       w_qmin( 851)=0.416757893E+00
       qk_min( 852)=0.271277966E+01
       w_qmin( 852)=0.416767655E+00
       qk_min( 853)=0.271576590E+01
       w_qmin( 853)=0.416777899E+00
       qk_min( 854)=0.271875214E+01
       w_qmin( 854)=0.416788451E+00
       qk_min( 855)=0.272173837E+01
       w_qmin( 855)=0.416799633E+00
       qk_min( 856)=0.272472461E+01
       w_qmin( 856)=0.416811034E+00
       qk_min( 857)=0.272771084E+01
       w_qmin( 857)=0.416823228E+00
       qk_min( 858)=0.273069708E+01
       w_qmin( 858)=0.416835674E+00
       qk_min( 859)=0.273368332E+01
       w_qmin( 859)=0.416848655E+00
       qk_min( 860)=0.273666955E+01
       w_qmin( 860)=0.416861949E+00
       qk_min( 861)=0.273965579E+01
       w_qmin( 861)=0.416875863E+00
       qk_min( 862)=0.274264202E+01
       w_qmin( 862)=0.416890002E+00
       qk_min( 863)=0.274562826E+01
       w_qmin( 863)=0.416904920E+00
       qk_min( 864)=0.274861450E+01
       w_qmin( 864)=0.416919949E+00
       qk_min( 865)=0.275160073E+01
       w_qmin( 865)=0.416935558E+00
       qk_min( 866)=0.275458697E+01
       w_qmin( 866)=0.416951582E+00
       qk_min( 867)=0.275757320E+01
       w_qmin( 867)=0.416967951E+00
       qk_min( 868)=0.276055944E+01
       w_qmin( 868)=0.416984879E+00
       qk_min( 869)=0.276354568E+01
       w_qmin( 869)=0.417002279E+00
       qk_min( 870)=0.276653191E+01
       w_qmin( 870)=0.417019985E+00
       qk_min( 871)=0.276951815E+01
       w_qmin( 871)=0.417038117E+00
       qk_min( 872)=0.277250438E+01
       w_qmin( 872)=0.417056789E+00
       qk_min( 873)=0.277549062E+01
       w_qmin( 873)=0.417075838E+00
       qk_min( 874)=0.277847686E+01
       w_qmin( 874)=0.417095599E+00
       qk_min( 875)=0.278146309E+01
       w_qmin( 875)=0.417115495E+00
       qk_min( 876)=0.278444933E+01
       w_qmin( 876)=0.417135994E+00
       qk_min( 877)=0.278743556E+01
       w_qmin( 877)=0.417156892E+00
       qk_min( 878)=0.279042180E+01
       w_qmin( 878)=0.417178162E+00
       qk_min( 879)=0.279340804E+01
       w_qmin( 879)=0.417200124E+00
       qk_min( 880)=0.279639427E+01
       w_qmin( 880)=0.417222614E+00
       qk_min( 881)=0.279938051E+01
       w_qmin( 881)=0.417245462E+00
       qk_min( 882)=0.280236674E+01
       w_qmin( 882)=0.417268802E+00
       qk_min( 883)=0.280535298E+01
       w_qmin( 883)=0.417292628E+00
       qk_min( 884)=0.280833922E+01
       w_qmin( 884)=0.417317034E+00
       qk_min( 885)=0.281132545E+01
       w_qmin( 885)=0.417342034E+00
       qk_min( 886)=0.281431169E+01
       w_qmin( 886)=0.417367399E+00
       qk_min( 887)=0.281729792E+01
       w_qmin( 887)=0.417393442E+00
       qk_min( 888)=0.282028416E+01
       w_qmin( 888)=0.417419853E+00
       qk_min( 889)=0.282327040E+01
       w_qmin( 889)=0.417446876E+00
       qk_min( 890)=0.282625663E+01
       w_qmin( 890)=0.417474595E+00
       qk_min( 891)=0.282924287E+01
       w_qmin( 891)=0.417502880E+00
       qk_min( 892)=0.283222910E+01
       w_qmin( 892)=0.417531699E+00
       qk_min( 893)=0.283521534E+01
       w_qmin( 893)=0.417560957E+00
       qk_min( 894)=0.283820158E+01
       w_qmin( 894)=0.417590881E+00
       qk_min( 895)=0.284118781E+01
       w_qmin( 895)=0.417621424E+00
       qk_min( 896)=0.284417405E+01
       w_qmin( 896)=0.417652705E+00
       qk_min( 897)=0.284716028E+01
       w_qmin( 897)=0.417684488E+00
       qk_min( 898)=0.285014652E+01
       w_qmin( 898)=0.417716795E+00
       qk_min( 899)=0.285313276E+01
       w_qmin( 899)=0.417749770E+00
       qk_min( 900)=0.285611899E+01
       w_qmin( 900)=0.417783415E+00
       qk_min( 901)=0.285910523E+01
       w_qmin( 901)=0.417817650E+00
       qk_min( 902)=0.286209146E+01
       w_qmin( 902)=0.417852756E+00
       qk_min( 903)=0.286507770E+01
       w_qmin( 903)=0.417888196E+00
       qk_min( 904)=0.286806393E+01
       w_qmin( 904)=0.417924287E+00
       qk_min( 905)=0.287105017E+01
       w_qmin( 905)=0.417961065E+00
       qk_min( 906)=0.287403641E+01
       w_qmin( 906)=0.417998507E+00
       qk_min( 907)=0.287702264E+01
       w_qmin( 907)=0.418036673E+00
       qk_min( 908)=0.288000888E+01
       w_qmin( 908)=0.418075365E+00
       qk_min( 909)=0.288299511E+01
       w_qmin( 909)=0.418114664E+00
       qk_min( 910)=0.288598135E+01
       w_qmin( 910)=0.418154584E+00
       qk_min( 911)=0.288896759E+01
       w_qmin( 911)=0.418195201E+00
       qk_min( 912)=0.289195382E+01
       w_qmin( 912)=0.418236414E+00
       qk_min( 913)=0.289494006E+01
       w_qmin( 913)=0.418278474E+00
       qk_min( 914)=0.289792629E+01
       w_qmin( 914)=0.418320916E+00
       qk_min( 915)=0.290091253E+01
       w_qmin( 915)=0.418363923E+00
       qk_min( 916)=0.290389877E+01
       w_qmin( 916)=0.418407579E+00
       qk_min( 917)=0.290688500E+01
       w_qmin( 917)=0.418451921E+00
       qk_min( 918)=0.290987124E+01
       w_qmin( 918)=0.418496849E+00
       qk_min( 919)=0.291285747E+01
       w_qmin( 919)=0.418542473E+00
       qk_min( 920)=0.291584371E+01
       w_qmin( 920)=0.418588574E+00
       qk_min( 921)=0.291882995E+01
       w_qmin( 921)=0.418635091E+00
       qk_min( 922)=0.292181618E+01
       w_qmin( 922)=0.418682232E+00
       qk_min( 923)=0.292480242E+01
       w_qmin( 923)=0.418730162E+00
       qk_min( 924)=0.292778865E+01
       w_qmin( 924)=0.418778683E+00
       qk_min( 925)=0.293077489E+01
       w_qmin( 925)=0.418827430E+00
       qk_min( 926)=0.293376113E+01
       w_qmin( 926)=0.418877012E+00
       qk_min( 927)=0.293674736E+01
       w_qmin( 927)=0.418926996E+00
       qk_min( 928)=0.293973360E+01
       w_qmin( 928)=0.418977476E+00
       qk_min( 929)=0.294271983E+01
       w_qmin( 929)=0.419028829E+00
       qk_min( 930)=0.294570607E+01
       w_qmin( 930)=0.419080234E+00
       qk_min( 931)=0.294869231E+01
       w_qmin( 931)=0.419132533E+00
       qk_min( 932)=0.295167854E+01
       w_qmin( 932)=0.419185102E+00
       qk_min( 933)=0.295466478E+01
       w_qmin( 933)=0.419238273E+00
       qk_min( 934)=0.295765101E+01
       w_qmin( 934)=0.419291930E+00
       qk_min( 935)=0.296063725E+01
       w_qmin( 935)=0.419346255E+00
       qk_min( 936)=0.296362349E+01
       w_qmin( 936)=0.419400775E+00
       qk_min( 937)=0.296660972E+01
       w_qmin( 937)=0.419455873E+00
       qk_min( 938)=0.296959596E+01
       w_qmin( 938)=0.419511534E+00
       qk_min( 939)=0.297258219E+01
       w_qmin( 939)=0.419567570E+00
       qk_min( 940)=0.297556843E+01
       w_qmin( 940)=0.419624267E+00
       qk_min( 941)=0.297855467E+01
       w_qmin( 941)=0.419681302E+00
       qk_min( 942)=0.298154090E+01
       w_qmin( 942)=0.419738750E+00
       qk_min( 943)=0.298452714E+01
       w_qmin( 943)=0.419796740E+00
       qk_min( 944)=0.298751337E+01
       w_qmin( 944)=0.419855138E+00
       qk_min( 945)=0.299049961E+01
       w_qmin( 945)=0.419913962E+00
       qk_min( 946)=0.299348585E+01
       w_qmin( 946)=0.419973372E+00
       qk_min( 947)=0.299647208E+01
       w_qmin( 947)=0.420033065E+00
       qk_min( 948)=0.299945832E+01
       w_qmin( 948)=0.420093276E+00
       qk_min( 949)=0.300244455E+01
       w_qmin( 949)=0.420153950E+00
       qk_min( 950)=0.300543079E+01
       w_qmin( 950)=0.420215081E+00
       qk_min( 951)=0.300841703E+01
       w_qmin( 951)=0.420276759E+00
       qk_min( 952)=0.301140326E+01
       w_qmin( 952)=0.420338767E+00
       qk_min( 953)=0.301438950E+01
       w_qmin( 953)=0.420401147E+00
       qk_min( 954)=0.301737573E+01
       w_qmin( 954)=0.420464001E+00
       qk_min( 955)=0.302036197E+01
       w_qmin( 955)=0.420527289E+00
       qk_min( 956)=0.302334821E+01
       w_qmin( 956)=0.420591146E+00
       qk_min( 957)=0.302633444E+01
       w_qmin( 957)=0.420655400E+00
       qk_min( 958)=0.302932068E+01
       w_qmin( 958)=0.420719993E+00
       qk_min( 959)=0.303230691E+01
       w_qmin( 959)=0.420785033E+00
       qk_min( 960)=0.303529315E+01
       w_qmin( 960)=0.420850599E+00
       qk_min( 961)=0.303827939E+01
       w_qmin( 961)=0.420916605E+00
       qk_min( 962)=0.304126562E+01
       w_qmin( 962)=0.420983049E+00
       qk_min( 963)=0.304425186E+01
       w_qmin( 963)=0.421049871E+00
       qk_min( 964)=0.304723809E+01
       w_qmin( 964)=0.421117200E+00
       qk_min( 965)=0.305022433E+01
       w_qmin( 965)=0.421185009E+00
       qk_min( 966)=0.305321056E+01
       w_qmin( 966)=0.421253111E+00
       qk_min( 967)=0.305619680E+01
       w_qmin( 967)=0.421321834E+00
       qk_min( 968)=0.305918304E+01
       w_qmin( 968)=0.421390933E+00
       qk_min( 969)=0.306216927E+01
       w_qmin( 969)=0.421460475E+00
       qk_min( 970)=0.306515551E+01
       w_qmin( 970)=0.421530529E+00
       qk_min( 971)=0.306814174E+01
       w_qmin( 971)=0.421600942E+00
       qk_min( 972)=0.307112798E+01
       w_qmin( 972)=0.421671707E+00
       qk_min( 973)=0.307411422E+01
       w_qmin( 973)=0.421743221E+00
       qk_min( 974)=0.307710045E+01
       w_qmin( 974)=0.421815057E+00
       qk_min( 975)=0.308008669E+01
       w_qmin( 975)=0.421887359E+00
       qk_min( 976)=0.308307292E+01
       w_qmin( 976)=0.421960061E+00
       qk_min( 977)=0.308605916E+01
       w_qmin( 977)=0.422033248E+00
       qk_min( 978)=0.308904540E+01
       w_qmin( 978)=0.422106903E+00
       qk_min( 979)=0.309203163E+01
       w_qmin( 979)=0.422181094E+00
       qk_min( 980)=0.309501787E+01
       w_qmin( 980)=0.422255771E+00
       qk_min( 981)=0.309800410E+01
       w_qmin( 981)=0.422330828E+00
       qk_min( 982)=0.310099034E+01
       w_qmin( 982)=0.422406359E+00
       qk_min( 983)=0.310397658E+01
       w_qmin( 983)=0.422482453E+00
       qk_min( 984)=0.310696281E+01
       w_qmin( 984)=0.422559102E+00
       qk_min( 985)=0.310994905E+01
       w_qmin( 985)=0.422636031E+00
       qk_min( 986)=0.311293528E+01
       w_qmin( 986)=0.422713635E+00
       qk_min( 987)=0.311592152E+01
       w_qmin( 987)=0.422791554E+00
       qk_min( 988)=0.311890776E+01
       w_qmin( 988)=0.422870031E+00
       qk_min( 989)=0.312189399E+01
       w_qmin( 989)=0.422949021E+00
       qk_min( 990)=0.312488023E+01
       w_qmin( 990)=0.423028519E+00
       qk_min( 991)=0.312786646E+01
       w_qmin( 991)=0.423108604E+00
       qk_min( 992)=0.313085270E+01
       w_qmin( 992)=0.423189024E+00
       qk_min( 993)=0.313383894E+01
       w_qmin( 993)=0.423269934E+00
       qk_min( 994)=0.313682517E+01
       w_qmin( 994)=0.423351608E+00
       qk_min( 995)=0.313981141E+01
       w_qmin( 995)=0.423433622E+00
       qk_min( 996)=0.314279764E+01
       w_qmin( 996)=0.423516074E+00
       qk_min( 997)=0.314578388E+01
       w_qmin( 997)=0.423599178E+00
       qk_min( 998)=0.314877012E+01
       w_qmin( 998)=0.423682779E+00
       qk_min( 999)=0.315175635E+01
       w_qmin( 999)=0.423766988E+00
       qk_min(1000)=0.315474259E+01
       w_qmin(1000)=0.423851339E+00
       qk_min(1001)=0.315772882E+01
       w_qmin(1001)=0.423936557E+00
       END
       SUBROUTINE InitA1TABKKpi(qk_min    ,w_qmin    )
       PARAMETER (Nq=1001)
       real*8               qk_min(Nq),w_qmin(Nq)
 ! this routine was automatically created by  a1-tabler
       qk_min(   1)=0.171492857E+00
       w_qmin(   1)=0.000000000E+00
       qk_min(   2)=0.174479093E+00
       w_qmin(   2)=0.000000000E+00
       qk_min(   3)=0.177465328E+00
       w_qmin(   3)=0.000000000E+00
       qk_min(   4)=0.180451564E+00
       w_qmin(   4)=0.000000000E+00
       qk_min(   5)=0.183437800E+00
       w_qmin(   5)=0.000000000E+00
       qk_min(   6)=0.186424036E+00
       w_qmin(   6)=0.000000000E+00
       qk_min(   7)=0.189410272E+00
       w_qmin(   7)=0.000000000E+00
       qk_min(   8)=0.192396508E+00
       w_qmin(   8)=0.000000000E+00
       qk_min(   9)=0.195382744E+00
       w_qmin(   9)=0.000000000E+00
       qk_min(  10)=0.198368980E+00
       w_qmin(  10)=0.000000000E+00
       qk_min(  11)=0.201355216E+00
       w_qmin(  11)=0.000000000E+00
       qk_min(  12)=0.204341452E+00
       w_qmin(  12)=0.000000000E+00
       qk_min(  13)=0.207327688E+00
       w_qmin(  13)=0.000000000E+00
       qk_min(  14)=0.210313924E+00
       w_qmin(  14)=0.000000000E+00
       qk_min(  15)=0.213300160E+00
       w_qmin(  15)=0.000000000E+00
       qk_min(  16)=0.216286396E+00
       w_qmin(  16)=0.000000000E+00
       qk_min(  17)=0.219272632E+00
       w_qmin(  17)=0.000000000E+00
       qk_min(  18)=0.222258868E+00
       w_qmin(  18)=0.000000000E+00
       qk_min(  19)=0.225245104E+00
       w_qmin(  19)=0.000000000E+00
       qk_min(  20)=0.228231340E+00
       w_qmin(  20)=0.000000000E+00
       qk_min(  21)=0.231217576E+00
       w_qmin(  21)=0.000000000E+00
       qk_min(  22)=0.234203812E+00
       w_qmin(  22)=0.000000000E+00
       qk_min(  23)=0.237190048E+00
       w_qmin(  23)=0.000000000E+00
       qk_min(  24)=0.240176284E+00
       w_qmin(  24)=0.000000000E+00
       qk_min(  25)=0.243162520E+00
       w_qmin(  25)=0.000000000E+00
       qk_min(  26)=0.246148756E+00
       w_qmin(  26)=0.000000000E+00
       qk_min(  27)=0.249134992E+00
       w_qmin(  27)=0.000000000E+00
       qk_min(  28)=0.252121228E+00
       w_qmin(  28)=0.000000000E+00
       qk_min(  29)=0.255107464E+00
       w_qmin(  29)=0.000000000E+00
       qk_min(  30)=0.258093700E+00
       w_qmin(  30)=0.000000000E+00
       qk_min(  31)=0.261079936E+00
       w_qmin(  31)=0.000000000E+00
       qk_min(  32)=0.264066172E+00
       w_qmin(  32)=0.000000000E+00
       qk_min(  33)=0.267052408E+00
       w_qmin(  33)=0.000000000E+00
       qk_min(  34)=0.270038643E+00
       w_qmin(  34)=0.000000000E+00
       qk_min(  35)=0.273024879E+00
       w_qmin(  35)=0.000000000E+00
       qk_min(  36)=0.276011115E+00
       w_qmin(  36)=0.000000000E+00
       qk_min(  37)=0.278997351E+00
       w_qmin(  37)=0.000000000E+00
       qk_min(  38)=0.281983587E+00
       w_qmin(  38)=0.000000000E+00
       qk_min(  39)=0.284969823E+00
       w_qmin(  39)=0.000000000E+00
       qk_min(  40)=0.287956059E+00
       w_qmin(  40)=0.000000000E+00
       qk_min(  41)=0.290942295E+00
       w_qmin(  41)=0.000000000E+00
       qk_min(  42)=0.293928531E+00
       w_qmin(  42)=0.000000000E+00
       qk_min(  43)=0.296914767E+00
       w_qmin(  43)=0.000000000E+00
       qk_min(  44)=0.299901003E+00
       w_qmin(  44)=0.000000000E+00
       qk_min(  45)=0.302887239E+00
       w_qmin(  45)=0.000000000E+00
       qk_min(  46)=0.305873475E+00
       w_qmin(  46)=0.000000000E+00
       qk_min(  47)=0.308859711E+00
       w_qmin(  47)=0.000000000E+00
       qk_min(  48)=0.311845947E+00
       w_qmin(  48)=0.000000000E+00
       qk_min(  49)=0.314832183E+00
       w_qmin(  49)=0.000000000E+00
       qk_min(  50)=0.317818419E+00
       w_qmin(  50)=0.000000000E+00
       qk_min(  51)=0.320804655E+00
       w_qmin(  51)=0.000000000E+00
       qk_min(  52)=0.323790891E+00
       w_qmin(  52)=0.000000000E+00
       qk_min(  53)=0.326777127E+00
       w_qmin(  53)=0.000000000E+00
       qk_min(  54)=0.329763363E+00
       w_qmin(  54)=0.000000000E+00
       qk_min(  55)=0.332749599E+00
       w_qmin(  55)=0.000000000E+00
       qk_min(  56)=0.335735835E+00
       w_qmin(  56)=0.000000000E+00
       qk_min(  57)=0.338722071E+00
       w_qmin(  57)=0.000000000E+00
       qk_min(  58)=0.341708307E+00
       w_qmin(  58)=0.000000000E+00
       qk_min(  59)=0.344694543E+00
       w_qmin(  59)=0.000000000E+00
       qk_min(  60)=0.347680779E+00
       w_qmin(  60)=0.000000000E+00
       qk_min(  61)=0.350667015E+00
       w_qmin(  61)=0.000000000E+00
       qk_min(  62)=0.353653251E+00
       w_qmin(  62)=0.000000000E+00
       qk_min(  63)=0.356639487E+00
       w_qmin(  63)=0.000000000E+00
       qk_min(  64)=0.359625722E+00
       w_qmin(  64)=0.000000000E+00
       qk_min(  65)=0.362611958E+00
       w_qmin(  65)=0.000000000E+00
       qk_min(  66)=0.365598194E+00
       w_qmin(  66)=0.000000000E+00
       qk_min(  67)=0.368584430E+00
       w_qmin(  67)=0.000000000E+00
       qk_min(  68)=0.371570666E+00
       w_qmin(  68)=0.000000000E+00
       qk_min(  69)=0.374556902E+00
       w_qmin(  69)=0.000000000E+00
       qk_min(  70)=0.377543138E+00
       w_qmin(  70)=0.000000000E+00
       qk_min(  71)=0.380529374E+00
       w_qmin(  71)=0.000000000E+00
       qk_min(  72)=0.383515610E+00
       w_qmin(  72)=0.000000000E+00
       qk_min(  73)=0.386501846E+00
       w_qmin(  73)=0.000000000E+00
       qk_min(  74)=0.389488082E+00
       w_qmin(  74)=0.000000000E+00
       qk_min(  75)=0.392474318E+00
       w_qmin(  75)=0.000000000E+00
       qk_min(  76)=0.395460554E+00
       w_qmin(  76)=0.000000000E+00
       qk_min(  77)=0.398446790E+00
       w_qmin(  77)=0.000000000E+00
       qk_min(  78)=0.401433026E+00
       w_qmin(  78)=0.000000000E+00
       qk_min(  79)=0.404419262E+00
       w_qmin(  79)=0.000000000E+00
       qk_min(  80)=0.407405498E+00
       w_qmin(  80)=0.000000000E+00
       qk_min(  81)=0.410391734E+00
       w_qmin(  81)=0.000000000E+00
       qk_min(  82)=0.413377970E+00
       w_qmin(  82)=0.000000000E+00
       qk_min(  83)=0.416364206E+00
       w_qmin(  83)=0.000000000E+00
       qk_min(  84)=0.419350442E+00
       w_qmin(  84)=0.000000000E+00
       qk_min(  85)=0.422336678E+00
       w_qmin(  85)=0.000000000E+00
       qk_min(  86)=0.425322914E+00
       w_qmin(  86)=0.000000000E+00
       qk_min(  87)=0.428309150E+00
       w_qmin(  87)=0.000000000E+00
       qk_min(  88)=0.431295386E+00
       w_qmin(  88)=0.000000000E+00
       qk_min(  89)=0.434281622E+00
       w_qmin(  89)=0.000000000E+00
       qk_min(  90)=0.437267858E+00
       w_qmin(  90)=0.000000000E+00
       qk_min(  91)=0.440254094E+00
       w_qmin(  91)=0.000000000E+00
       qk_min(  92)=0.443240330E+00
       w_qmin(  92)=0.000000000E+00
       qk_min(  93)=0.446226566E+00
       w_qmin(  93)=0.000000000E+00
       qk_min(  94)=0.449212802E+00
       w_qmin(  94)=0.000000000E+00
       qk_min(  95)=0.452199037E+00
       w_qmin(  95)=0.000000000E+00
       qk_min(  96)=0.455185273E+00
       w_qmin(  96)=0.000000000E+00
       qk_min(  97)=0.458171509E+00
       w_qmin(  97)=0.000000000E+00
       qk_min(  98)=0.461157745E+00
       w_qmin(  98)=0.000000000E+00
       qk_min(  99)=0.464143981E+00
       w_qmin(  99)=0.000000000E+00
       qk_min( 100)=0.467130217E+00
       w_qmin( 100)=0.000000000E+00
       qk_min( 101)=0.470116453E+00
       w_qmin( 101)=0.000000000E+00
       qk_min( 102)=0.473102689E+00
       w_qmin( 102)=0.000000000E+00
       qk_min( 103)=0.476088925E+00
       w_qmin( 103)=0.000000000E+00
       qk_min( 104)=0.479075161E+00
       w_qmin( 104)=0.000000000E+00
       qk_min( 105)=0.482061397E+00
       w_qmin( 105)=0.000000000E+00
       qk_min( 106)=0.485047633E+00
       w_qmin( 106)=0.000000000E+00
       qk_min( 107)=0.488033869E+00
       w_qmin( 107)=0.000000000E+00
       qk_min( 108)=0.491020105E+00
       w_qmin( 108)=0.000000000E+00
       qk_min( 109)=0.494006341E+00
       w_qmin( 109)=0.000000000E+00
       qk_min( 110)=0.496992577E+00
       w_qmin( 110)=0.000000000E+00
       qk_min( 111)=0.499978813E+00
       w_qmin( 111)=0.000000000E+00
       qk_min( 112)=0.502965049E+00
       w_qmin( 112)=0.000000000E+00
       qk_min( 113)=0.505951285E+00
       w_qmin( 113)=0.000000000E+00
       qk_min( 114)=0.508937521E+00
       w_qmin( 114)=0.000000000E+00
       qk_min( 115)=0.511923757E+00
       w_qmin( 115)=0.000000000E+00
       qk_min( 116)=0.514909993E+00
       w_qmin( 116)=0.000000000E+00
       qk_min( 117)=0.517896229E+00
       w_qmin( 117)=0.000000000E+00
       qk_min( 118)=0.520882465E+00
       w_qmin( 118)=0.000000000E+00
       qk_min( 119)=0.523868701E+00
       w_qmin( 119)=0.000000000E+00
       qk_min( 120)=0.526854937E+00
       w_qmin( 120)=0.000000000E+00
       qk_min( 121)=0.529841173E+00
       w_qmin( 121)=0.000000000E+00
       qk_min( 122)=0.532827409E+00
       w_qmin( 122)=0.000000000E+00
       qk_min( 123)=0.535813645E+00
       w_qmin( 123)=0.000000000E+00
       qk_min( 124)=0.538799881E+00
       w_qmin( 124)=0.000000000E+00
       qk_min( 125)=0.541786116E+00
       w_qmin( 125)=0.000000000E+00
       qk_min( 126)=0.544772352E+00
       w_qmin( 126)=0.000000000E+00
       qk_min( 127)=0.547758588E+00
       w_qmin( 127)=0.000000000E+00
       qk_min( 128)=0.550744824E+00
       w_qmin( 128)=0.000000000E+00
       qk_min( 129)=0.553731060E+00
       w_qmin( 129)=0.000000000E+00
       qk_min( 130)=0.556717296E+00
       w_qmin( 130)=0.000000000E+00
       qk_min( 131)=0.559703532E+00
       w_qmin( 131)=0.000000000E+00
       qk_min( 132)=0.562689768E+00
       w_qmin( 132)=0.000000000E+00
       qk_min( 133)=0.565676004E+00
       w_qmin( 133)=0.000000000E+00
       qk_min( 134)=0.568662240E+00
       w_qmin( 134)=0.000000000E+00
       qk_min( 135)=0.571648476E+00
       w_qmin( 135)=0.000000000E+00
       qk_min( 136)=0.574634712E+00
       w_qmin( 136)=0.000000000E+00
       qk_min( 137)=0.577620948E+00
       w_qmin( 137)=0.000000000E+00
       qk_min( 138)=0.580607184E+00
       w_qmin( 138)=0.000000000E+00
       qk_min( 139)=0.583593420E+00
       w_qmin( 139)=0.000000000E+00
       qk_min( 140)=0.586579656E+00
       w_qmin( 140)=0.000000000E+00
       qk_min( 141)=0.589565892E+00
       w_qmin( 141)=0.000000000E+00
       qk_min( 142)=0.592552128E+00
       w_qmin( 142)=0.000000000E+00
       qk_min( 143)=0.595538364E+00
       w_qmin( 143)=0.000000000E+00
       qk_min( 144)=0.598524600E+00
       w_qmin( 144)=0.000000000E+00
       qk_min( 145)=0.601510836E+00
       w_qmin( 145)=0.000000000E+00
       qk_min( 146)=0.604497072E+00
       w_qmin( 146)=0.000000000E+00
       qk_min( 147)=0.607483308E+00
       w_qmin( 147)=0.000000000E+00
       qk_min( 148)=0.610469544E+00
       w_qmin( 148)=0.000000000E+00
       qk_min( 149)=0.613455780E+00
       w_qmin( 149)=0.000000000E+00
       qk_min( 150)=0.616442016E+00
       w_qmin( 150)=0.000000000E+00
       qk_min( 151)=0.619428252E+00
       w_qmin( 151)=0.000000000E+00
       qk_min( 152)=0.622414488E+00
       w_qmin( 152)=0.000000000E+00
       qk_min( 153)=0.625400724E+00
       w_qmin( 153)=0.000000000E+00
       qk_min( 154)=0.628386960E+00
       w_qmin( 154)=0.000000000E+00
       qk_min( 155)=0.631373196E+00
       w_qmin( 155)=0.000000000E+00
       qk_min( 156)=0.634359431E+00
       w_qmin( 156)=0.000000000E+00
       qk_min( 157)=0.637345667E+00
       w_qmin( 157)=0.000000000E+00
       qk_min( 158)=0.640331903E+00
       w_qmin( 158)=0.000000000E+00
       qk_min( 159)=0.643318139E+00
       w_qmin( 159)=0.000000000E+00
       qk_min( 160)=0.646304375E+00
       w_qmin( 160)=0.000000000E+00
       qk_min( 161)=0.649290611E+00
       w_qmin( 161)=0.000000000E+00
       qk_min( 162)=0.652276847E+00
       w_qmin( 162)=0.000000000E+00
       qk_min( 163)=0.655263083E+00
       w_qmin( 163)=0.000000000E+00
       qk_min( 164)=0.658249319E+00
       w_qmin( 164)=0.000000000E+00
       qk_min( 165)=0.661235555E+00
       w_qmin( 165)=0.000000000E+00
       qk_min( 166)=0.664221791E+00
       w_qmin( 166)=0.000000000E+00
       qk_min( 167)=0.667208027E+00
       w_qmin( 167)=0.000000000E+00
       qk_min( 168)=0.670194263E+00
       w_qmin( 168)=0.000000000E+00
       qk_min( 169)=0.673180499E+00
       w_qmin( 169)=0.000000000E+00
       qk_min( 170)=0.676166735E+00
       w_qmin( 170)=0.000000000E+00
       qk_min( 171)=0.679152971E+00
       w_qmin( 171)=0.000000000E+00
       qk_min( 172)=0.682139207E+00
       w_qmin( 172)=0.000000000E+00
       qk_min( 173)=0.685125443E+00
       w_qmin( 173)=0.000000000E+00
       qk_min( 174)=0.688111679E+00
       w_qmin( 174)=0.000000000E+00
       qk_min( 175)=0.691097915E+00
       w_qmin( 175)=0.000000000E+00
       qk_min( 176)=0.694084151E+00
       w_qmin( 176)=0.000000000E+00
       qk_min( 177)=0.697070387E+00
       w_qmin( 177)=0.000000000E+00
       qk_min( 178)=0.700056623E+00
       w_qmin( 178)=0.000000000E+00
       qk_min( 179)=0.703042859E+00
       w_qmin( 179)=0.000000000E+00
       qk_min( 180)=0.706029095E+00
       w_qmin( 180)=0.000000000E+00
       qk_min( 181)=0.709015331E+00
       w_qmin( 181)=0.000000000E+00
       qk_min( 182)=0.712001567E+00
       w_qmin( 182)=0.000000000E+00
       qk_min( 183)=0.714987803E+00
       w_qmin( 183)=0.000000000E+00
       qk_min( 184)=0.717974039E+00
       w_qmin( 184)=0.000000000E+00
       qk_min( 185)=0.720960275E+00
       w_qmin( 185)=0.000000000E+00
       qk_min( 186)=0.723946510E+00
       w_qmin( 186)=0.000000000E+00
       qk_min( 187)=0.726932746E+00
       w_qmin( 187)=0.000000000E+00
       qk_min( 188)=0.729918982E+00
       w_qmin( 188)=0.000000000E+00
       qk_min( 189)=0.732905218E+00
       w_qmin( 189)=0.000000000E+00
       qk_min( 190)=0.735891454E+00
       w_qmin( 190)=0.000000000E+00
       qk_min( 191)=0.738877690E+00
       w_qmin( 191)=0.000000000E+00
       qk_min( 192)=0.741863926E+00
       w_qmin( 192)=0.000000000E+00
       qk_min( 193)=0.744850162E+00
       w_qmin( 193)=0.000000000E+00
       qk_min( 194)=0.747836398E+00
       w_qmin( 194)=0.000000000E+00
       qk_min( 195)=0.750822634E+00
       w_qmin( 195)=0.000000000E+00
       qk_min( 196)=0.753808870E+00
       w_qmin( 196)=0.000000000E+00
       qk_min( 197)=0.756795106E+00
       w_qmin( 197)=0.000000000E+00
       qk_min( 198)=0.759781342E+00
       w_qmin( 198)=0.000000000E+00
       qk_min( 199)=0.762767578E+00
       w_qmin( 199)=0.000000000E+00
       qk_min( 200)=0.765753814E+00
       w_qmin( 200)=0.000000000E+00
       qk_min( 201)=0.768740050E+00
       w_qmin( 201)=0.000000000E+00
       qk_min( 202)=0.771726286E+00
       w_qmin( 202)=0.000000000E+00
       qk_min( 203)=0.774712522E+00
       w_qmin( 203)=0.000000000E+00
       qk_min( 204)=0.777698758E+00
       w_qmin( 204)=0.000000000E+00
       qk_min( 205)=0.780684994E+00
       w_qmin( 205)=0.000000000E+00
       qk_min( 206)=0.783671230E+00
       w_qmin( 206)=0.000000000E+00
       qk_min( 207)=0.786657466E+00
       w_qmin( 207)=0.000000000E+00
       qk_min( 208)=0.789643702E+00
       w_qmin( 208)=0.000000000E+00
       qk_min( 209)=0.792629938E+00
       w_qmin( 209)=0.000000000E+00
       qk_min( 210)=0.795616174E+00
       w_qmin( 210)=0.000000000E+00
       qk_min( 211)=0.798602410E+00
       w_qmin( 211)=0.000000000E+00
       qk_min( 212)=0.801588646E+00
       w_qmin( 212)=0.000000000E+00
       qk_min( 213)=0.804574882E+00
       w_qmin( 213)=0.000000000E+00
       qk_min( 214)=0.807561118E+00
       w_qmin( 214)=0.000000000E+00
       qk_min( 215)=0.810547354E+00
       w_qmin( 215)=0.000000000E+00
       qk_min( 216)=0.813533590E+00
       w_qmin( 216)=0.000000000E+00
       qk_min( 217)=0.816519825E+00
       w_qmin( 217)=0.000000000E+00
       qk_min( 218)=0.819506061E+00
       w_qmin( 218)=0.000000000E+00
       qk_min( 219)=0.822492297E+00
       w_qmin( 219)=0.000000000E+00
       qk_min( 220)=0.825478533E+00
       w_qmin( 220)=0.000000000E+00
       qk_min( 221)=0.828464769E+00
       w_qmin( 221)=0.000000000E+00
       qk_min( 222)=0.831451005E+00
       w_qmin( 222)=0.000000000E+00
       qk_min( 223)=0.834437241E+00
       w_qmin( 223)=0.000000000E+00
       qk_min( 224)=0.837423477E+00
       w_qmin( 224)=0.000000000E+00
       qk_min( 225)=0.840409713E+00
       w_qmin( 225)=0.000000000E+00
       qk_min( 226)=0.843395949E+00
       w_qmin( 226)=0.000000000E+00
       qk_min( 227)=0.846382185E+00
       w_qmin( 227)=0.000000000E+00
       qk_min( 228)=0.849368421E+00
       w_qmin( 228)=0.000000000E+00
       qk_min( 229)=0.852354657E+00
       w_qmin( 229)=0.000000000E+00
       qk_min( 230)=0.855340893E+00
       w_qmin( 230)=0.000000000E+00
       qk_min( 231)=0.858327129E+00
       w_qmin( 231)=0.000000000E+00
       qk_min( 232)=0.861313365E+00
       w_qmin( 232)=0.000000000E+00
       qk_min( 233)=0.864299601E+00
       w_qmin( 233)=0.000000000E+00
       qk_min( 234)=0.867285837E+00
       w_qmin( 234)=0.000000000E+00
       qk_min( 235)=0.870272073E+00
       w_qmin( 235)=0.000000000E+00
       qk_min( 236)=0.873258309E+00
       w_qmin( 236)=0.000000000E+00
       qk_min( 237)=0.876244545E+00
       w_qmin( 237)=0.000000000E+00
       qk_min( 238)=0.879230781E+00
       w_qmin( 238)=0.000000000E+00
       qk_min( 239)=0.882217017E+00
       w_qmin( 239)=0.000000000E+00
       qk_min( 240)=0.885203253E+00
       w_qmin( 240)=0.000000000E+00
       qk_min( 241)=0.888189489E+00
       w_qmin( 241)=0.000000000E+00
       qk_min( 242)=0.891175725E+00
       w_qmin( 242)=0.000000000E+00
       qk_min( 243)=0.894161961E+00
       w_qmin( 243)=0.000000000E+00
       qk_min( 244)=0.897148197E+00
       w_qmin( 244)=0.000000000E+00
       qk_min( 245)=0.900134433E+00
       w_qmin( 245)=0.000000000E+00
       qk_min( 246)=0.903120669E+00
       w_qmin( 246)=0.000000000E+00
       qk_min( 247)=0.906106904E+00
       w_qmin( 247)=0.000000000E+00
       qk_min( 248)=0.909093140E+00
       w_qmin( 248)=0.000000000E+00
       qk_min( 249)=0.912079376E+00
       w_qmin( 249)=0.000000000E+00
       qk_min( 250)=0.915065612E+00
       w_qmin( 250)=0.000000000E+00
       qk_min( 251)=0.918051848E+00
       w_qmin( 251)=0.000000000E+00
       qk_min( 252)=0.921038084E+00
       w_qmin( 252)=0.000000000E+00
       qk_min( 253)=0.924024320E+00
       w_qmin( 253)=0.000000000E+00
       qk_min( 254)=0.927010556E+00
       w_qmin( 254)=0.000000000E+00
       qk_min( 255)=0.929996792E+00
       w_qmin( 255)=0.000000000E+00
       qk_min( 256)=0.932983028E+00
       w_qmin( 256)=0.000000000E+00
       qk_min( 257)=0.935969264E+00
       w_qmin( 257)=0.000000000E+00
       qk_min( 258)=0.938955500E+00
       w_qmin( 258)=0.000000000E+00
       qk_min( 259)=0.941941736E+00
       w_qmin( 259)=0.000000000E+00
       qk_min( 260)=0.944927972E+00
       w_qmin( 260)=0.000000000E+00
       qk_min( 261)=0.947914208E+00
       w_qmin( 261)=0.000000000E+00
       qk_min( 262)=0.950900444E+00
       w_qmin( 262)=0.000000000E+00
       qk_min( 263)=0.953886680E+00
       w_qmin( 263)=0.000000000E+00
       qk_min( 264)=0.956872916E+00
       w_qmin( 264)=0.000000000E+00
       qk_min( 265)=0.959859152E+00
       w_qmin( 265)=0.000000000E+00
       qk_min( 266)=0.962845388E+00
       w_qmin( 266)=0.000000000E+00
       qk_min( 267)=0.965831624E+00
       w_qmin( 267)=0.000000000E+00
       qk_min( 268)=0.968817860E+00
       w_qmin( 268)=0.000000000E+00
       qk_min( 269)=0.971804096E+00
       w_qmin( 269)=0.000000000E+00
       qk_min( 270)=0.974790332E+00
       w_qmin( 270)=0.000000000E+00
       qk_min( 271)=0.977776568E+00
       w_qmin( 271)=0.000000000E+00
       qk_min( 272)=0.980762804E+00
       w_qmin( 272)=0.000000000E+00
       qk_min( 273)=0.983749040E+00
       w_qmin( 273)=0.000000000E+00
       qk_min( 274)=0.986735276E+00
       w_qmin( 274)=0.000000000E+00
       qk_min( 275)=0.989721512E+00
       w_qmin( 275)=0.000000000E+00
       qk_min( 276)=0.992707748E+00
       w_qmin( 276)=0.000000000E+00
       qk_min( 277)=0.995693984E+00
       w_qmin( 277)=0.000000000E+00
       qk_min( 278)=0.998680219E+00
       w_qmin( 278)=0.000000000E+00
       qk_min( 279)=0.100166646E+01
       w_qmin( 279)=0.000000000E+00
       qk_min( 280)=0.100465269E+01
       w_qmin( 280)=0.000000000E+00
       qk_min( 281)=0.100763893E+01
       w_qmin( 281)=0.000000000E+00
       qk_min( 282)=0.101062516E+01
       w_qmin( 282)=0.000000000E+00
       qk_min( 283)=0.101361140E+01
       w_qmin( 283)=0.000000000E+00
       qk_min( 284)=0.101659764E+01
       w_qmin( 284)=0.000000000E+00
       qk_min( 285)=0.101958387E+01
       w_qmin( 285)=0.000000000E+00
       qk_min( 286)=0.102257011E+01
       w_qmin( 286)=0.000000000E+00
       qk_min( 287)=0.102555634E+01
       w_qmin( 287)=0.000000000E+00
       qk_min( 288)=0.102854258E+01
       w_qmin( 288)=0.000000000E+00
       qk_min( 289)=0.103152882E+01
       w_qmin( 289)=0.000000000E+00
       qk_min( 290)=0.103451505E+01
       w_qmin( 290)=0.000000000E+00
       qk_min( 291)=0.103750129E+01
       w_qmin( 291)=0.000000000E+00
       qk_min( 292)=0.104048752E+01
       w_qmin( 292)=0.000000000E+00
       qk_min( 293)=0.104347376E+01
       w_qmin( 293)=0.000000000E+00
       qk_min( 294)=0.104645999E+01
       w_qmin( 294)=0.000000000E+00
       qk_min( 295)=0.104944623E+01
       w_qmin( 295)=0.000000000E+00
       qk_min( 296)=0.105243247E+01
       w_qmin( 296)=0.000000000E+00
       qk_min( 297)=0.105541870E+01
       w_qmin( 297)=0.000000000E+00
       qk_min( 298)=0.105840494E+01
       w_qmin( 298)=0.000000000E+00
       qk_min( 299)=0.106139117E+01
       w_qmin( 299)=0.000000000E+00
       qk_min( 300)=0.106437741E+01
       w_qmin( 300)=0.000000000E+00
       qk_min( 301)=0.106736365E+01
       w_qmin( 301)=0.000000000E+00
       qk_min( 302)=0.107034988E+01
       w_qmin( 302)=0.000000000E+00
       qk_min( 303)=0.107333612E+01
       w_qmin( 303)=0.000000000E+00
       qk_min( 304)=0.107632235E+01
       w_qmin( 304)=0.000000000E+00
       qk_min( 305)=0.107930859E+01
       w_qmin( 305)=0.000000000E+00
       qk_min( 306)=0.108229483E+01
       w_qmin( 306)=0.000000000E+00
       qk_min( 307)=0.108528106E+01
       w_qmin( 307)=0.000000000E+00
       qk_min( 308)=0.108826730E+01
       w_qmin( 308)=0.000000000E+00
       qk_min( 309)=0.109125353E+01
       w_qmin( 309)=0.000000000E+00
       qk_min( 310)=0.109423977E+01
       w_qmin( 310)=0.000000000E+00
       qk_min( 311)=0.109722601E+01
       w_qmin( 311)=0.000000000E+00
       qk_min( 312)=0.110021224E+01
       w_qmin( 312)=0.000000000E+00
       qk_min( 313)=0.110319848E+01
       w_qmin( 313)=0.000000000E+00
       qk_min( 314)=0.110618471E+01
       w_qmin( 314)=0.000000000E+00
       qk_min( 315)=0.110917095E+01
       w_qmin( 315)=0.000000000E+00
       qk_min( 316)=0.111215719E+01
       w_qmin( 316)=0.000000000E+00
       qk_min( 317)=0.111514342E+01
       w_qmin( 317)=0.000000000E+00
       qk_min( 318)=0.111812966E+01
       w_qmin( 318)=0.000000000E+00
       qk_min( 319)=0.112111589E+01
       w_qmin( 319)=0.000000000E+00
       qk_min( 320)=0.112410213E+01
       w_qmin( 320)=0.000000000E+00
       qk_min( 321)=0.112708837E+01
       w_qmin( 321)=0.000000000E+00
       qk_min( 322)=0.113007460E+01
       w_qmin( 322)=0.000000000E+00
       qk_min( 323)=0.113306084E+01
       w_qmin( 323)=0.000000000E+00
       qk_min( 324)=0.113604707E+01
       w_qmin( 324)=0.000000000E+00
       qk_min( 325)=0.113903331E+01
       w_qmin( 325)=0.000000000E+00
       qk_min( 326)=0.114201955E+01
       w_qmin( 326)=0.000000000E+00
       qk_min( 327)=0.114500578E+01
       w_qmin( 327)=0.000000000E+00
       qk_min( 328)=0.114799202E+01
       w_qmin( 328)=0.000000000E+00
       qk_min( 329)=0.115097825E+01
       w_qmin( 329)=0.000000000E+00
       qk_min( 330)=0.115396449E+01
       w_qmin( 330)=0.000000000E+00
       qk_min( 331)=0.115695073E+01
       w_qmin( 331)=0.000000000E+00
       qk_min( 332)=0.115993696E+01
       w_qmin( 332)=0.000000000E+00
       qk_min( 333)=0.116292320E+01
       w_qmin( 333)=0.000000000E+00
       qk_min( 334)=0.116590943E+01
       w_qmin( 334)=0.000000000E+00
       qk_min( 335)=0.116889567E+01
       w_qmin( 335)=0.000000000E+00
       qk_min( 336)=0.117188191E+01
       w_qmin( 336)=0.000000000E+00
       qk_min( 337)=0.117486814E+01
       w_qmin( 337)=0.000000000E+00
       qk_min( 338)=0.117785438E+01
       w_qmin( 338)=0.000000000E+00
       qk_min( 339)=0.118084061E+01
       w_qmin( 339)=0.000000000E+00
       qk_min( 340)=0.118382685E+01
       w_qmin( 340)=0.000000000E+00
       qk_min( 341)=0.118681309E+01
       w_qmin( 341)=0.000000000E+00
       qk_min( 342)=0.118979932E+01
       w_qmin( 342)=0.000000000E+00
       qk_min( 343)=0.119278556E+01
       w_qmin( 343)=0.000000000E+00
       qk_min( 344)=0.119577179E+01
       w_qmin( 344)=0.000000000E+00
       qk_min( 345)=0.119875803E+01
       w_qmin( 345)=0.000000000E+00
       qk_min( 346)=0.120174427E+01
       w_qmin( 346)=0.000000000E+00
       qk_min( 347)=0.120473050E+01
       w_qmin( 347)=0.000000000E+00
       qk_min( 348)=0.120771674E+01
       w_qmin( 348)=0.000000000E+00
       qk_min( 349)=0.121070297E+01
       w_qmin( 349)=0.000000000E+00
       qk_min( 350)=0.121368921E+01
       w_qmin( 350)=0.000000000E+00
       qk_min( 351)=0.121667545E+01
       w_qmin( 351)=0.000000000E+00
       qk_min( 352)=0.121966168E+01
       w_qmin( 352)=0.000000000E+00
       qk_min( 353)=0.122264792E+01
       w_qmin( 353)=0.000000000E+00
       qk_min( 354)=0.122563415E+01
       w_qmin( 354)=0.000000000E+00
       qk_min( 355)=0.122862039E+01
       w_qmin( 355)=0.000000000E+00
       qk_min( 356)=0.123160662E+01
       w_qmin( 356)=0.000000000E+00
       qk_min( 357)=0.123459286E+01
       w_qmin( 357)=0.000000000E+00
       qk_min( 358)=0.123757910E+01
       w_qmin( 358)=0.000000000E+00
       qk_min( 359)=0.124056533E+01
       w_qmin( 359)=0.000000000E+00
       qk_min( 360)=0.124355157E+01
       w_qmin( 360)=0.000000000E+00
       qk_min( 361)=0.124653780E+01
       w_qmin( 361)=0.000000000E+00
       qk_min( 362)=0.124952404E+01
       w_qmin( 362)=0.000000000E+00
       qk_min( 363)=0.125251028E+01
       w_qmin( 363)=0.000000000E+00
       qk_min( 364)=0.125549651E+01
       w_qmin( 364)=0.000000000E+00
       qk_min( 365)=0.125848275E+01
       w_qmin( 365)=0.000000000E+00
       qk_min( 366)=0.126146898E+01
       w_qmin( 366)=0.000000000E+00
       qk_min( 367)=0.126445522E+01
       w_qmin( 367)=0.000000000E+00
       qk_min( 368)=0.126744146E+01
       w_qmin( 368)=0.000000000E+00
       qk_min( 369)=0.127042769E+01
       w_qmin( 369)=0.000000000E+00
       qk_min( 370)=0.127341393E+01
       w_qmin( 370)=0.000000000E+00
       qk_min( 371)=0.127640016E+01
       w_qmin( 371)=0.467385957E-11
       qk_min( 372)=0.127938640E+01
       w_qmin( 372)=0.336675634E-09
       qk_min( 373)=0.128237264E+01
       w_qmin( 373)=0.180589425E-08
       qk_min( 374)=0.128535887E+01
       w_qmin( 374)=0.521060960E-08
       qk_min( 375)=0.128834511E+01
       w_qmin( 375)=0.112917588E-07
       qk_min( 376)=0.129133134E+01
       w_qmin( 376)=0.207377627E-07
       qk_min( 377)=0.129431758E+01
       w_qmin( 377)=0.341888871E-07
       qk_min( 378)=0.129730382E+01
       w_qmin( 378)=0.522393827E-07
       qk_min( 379)=0.130029005E+01
       w_qmin( 379)=0.754446882E-07
       qk_min( 380)=0.130327629E+01
       w_qmin( 380)=0.104319866E-06
       qk_min( 381)=0.130626252E+01
       w_qmin( 381)=0.139349343E-06
       qk_min( 382)=0.130924876E+01
       w_qmin( 382)=0.180986475E-06
       qk_min( 383)=0.131223500E+01
       w_qmin( 383)=0.229654233E-06
       qk_min( 384)=0.131522123E+01
       w_qmin( 384)=0.285754285E-06
       qk_min( 385)=0.131820747E+01
       w_qmin( 385)=0.349669848E-06
       qk_min( 386)=0.132119370E+01
       w_qmin( 386)=0.421756559E-06
       qk_min( 387)=0.132417994E+01
       w_qmin( 387)=0.502363262E-06
       qk_min( 388)=0.132716618E+01
       w_qmin( 388)=0.591822610E-06
       qk_min( 389)=0.133015241E+01
       w_qmin( 389)=0.690452213E-06
       qk_min( 390)=0.133313865E+01
       w_qmin( 390)=0.798572685E-06
       qk_min( 391)=0.133612488E+01
       w_qmin( 391)=0.916480606E-06
       qk_min( 392)=0.133911112E+01
       w_qmin( 392)=0.104448389E-05
       qk_min( 393)=0.134209736E+01
       w_qmin( 393)=0.118288886E-05
       qk_min( 394)=0.134508359E+01
       w_qmin( 394)=0.133198624E-05
       qk_min( 395)=0.134806983E+01
       w_qmin( 395)=0.149208590E-05
       qk_min( 396)=0.135105606E+01
       w_qmin( 396)=0.166349876E-05
       qk_min( 397)=0.135404230E+01
       w_qmin( 397)=0.184652662E-05
       qk_min( 398)=0.135702854E+01
       w_qmin( 398)=0.204149727E-05
       qk_min( 399)=0.136001477E+01
       w_qmin( 399)=0.224874847E-05
       qk_min( 400)=0.136300101E+01
       w_qmin( 400)=0.246860209E-05
       qk_min( 401)=0.136598724E+01
       w_qmin( 401)=0.270143320E-05
       qk_min( 402)=0.136897348E+01
       w_qmin( 402)=0.294758441E-05
       qk_min( 403)=0.137195972E+01
       w_qmin( 403)=0.320745294E-05
       qk_min( 404)=0.137494595E+01
       w_qmin( 404)=0.348143917E-05
       qk_min( 405)=0.137793219E+01
       w_qmin( 405)=0.376994513E-05
       qk_min( 406)=0.138091842E+01
       w_qmin( 406)=0.407341659E-05
       qk_min( 407)=0.138390466E+01
       w_qmin( 407)=0.439231612E-05
       qk_min( 408)=0.138689090E+01
       w_qmin( 408)=0.472710275E-05
       qk_min( 409)=0.138987713E+01
       w_qmin( 409)=0.507828748E-05
       qk_min( 410)=0.139286337E+01
       w_qmin( 410)=0.544639815E-05
       qk_min( 411)=0.139584960E+01
       w_qmin( 411)=0.583197441E-05
       qk_min( 412)=0.139883584E+01
       w_qmin( 412)=0.623560306E-05
       qk_min( 413)=0.140182208E+01
       w_qmin( 413)=0.665787419E-05
       qk_min( 414)=0.140480831E+01
       w_qmin( 414)=0.709942081E-05
       qk_min( 415)=0.140779455E+01
       w_qmin( 415)=0.756091702E-05
       qk_min( 416)=0.141078078E+01
       w_qmin( 416)=0.804302920E-05
       qk_min( 417)=0.141376702E+01
       w_qmin( 417)=0.854648533E-05
       qk_min( 418)=0.141675325E+01
       w_qmin( 418)=0.907205680E-05
       qk_min( 419)=0.141973949E+01
       w_qmin( 419)=0.962050306E-05
       qk_min( 420)=0.142272573E+01
       w_qmin( 420)=0.101926529E-04
       qk_min( 421)=0.142571196E+01
       w_qmin( 421)=0.107893739E-04
       qk_min( 422)=0.142869820E+01
       w_qmin( 422)=0.114115273E-04
       qk_min( 423)=0.143168443E+01
       w_qmin( 423)=0.120600685E-04
       qk_min( 424)=0.143467067E+01
       w_qmin( 424)=0.127359231E-04
       qk_min( 425)=0.143765691E+01
       w_qmin( 425)=0.134401138E-04
       qk_min( 426)=0.144064314E+01
       w_qmin( 426)=0.141736967E-04
       qk_min( 427)=0.144362938E+01
       w_qmin( 427)=0.149377099E-04
       qk_min( 428)=0.144661561E+01
       w_qmin( 428)=0.157332917E-04
       qk_min( 429)=0.144960185E+01
       w_qmin( 429)=0.165616119E-04
       qk_min( 430)=0.145258809E+01
       w_qmin( 430)=0.174238532E-04
       qk_min( 431)=0.145557432E+01
       w_qmin( 431)=0.183212644E-04
       qk_min( 432)=0.145856056E+01
       w_qmin( 432)=0.192551614E-04
       qk_min( 433)=0.146154679E+01
       w_qmin( 433)=0.202268169E-04
       qk_min( 434)=0.146453303E+01
       w_qmin( 434)=0.212376507E-04
       qk_min( 435)=0.146751927E+01
       w_qmin( 435)=0.222890955E-04
       qk_min( 436)=0.147050550E+01
       w_qmin( 436)=0.233825795E-04
       qk_min( 437)=0.147349174E+01
       w_qmin( 437)=0.245196597E-04
       qk_min( 438)=0.147647797E+01
       w_qmin( 438)=0.257018674E-04
       qk_min( 439)=0.147946421E+01
       w_qmin( 439)=0.269308165E-04
       qk_min( 440)=0.148245045E+01
       w_qmin( 440)=0.282082199E-04
       qk_min( 441)=0.148543668E+01
       w_qmin( 441)=0.295357207E-04
       qk_min( 442)=0.148842292E+01
       w_qmin( 442)=0.309151355E-04
       qk_min( 443)=0.149140915E+01
       w_qmin( 443)=0.323482962E-04
       qk_min( 444)=0.149439539E+01
       w_qmin( 444)=0.338370198E-04
       qk_min( 445)=0.149738163E+01
       w_qmin( 445)=0.353832483E-04
       qk_min( 446)=0.150036786E+01
       w_qmin( 446)=0.369890317E-04
       qk_min( 447)=0.150335410E+01
       w_qmin( 447)=0.386563318E-04
       qk_min( 448)=0.150634033E+01
       w_qmin( 448)=0.403872831E-04
       qk_min( 449)=0.150932657E+01
       w_qmin( 449)=0.421839819E-04
       qk_min( 450)=0.151231281E+01
       w_qmin( 450)=0.440487213E-04
       qk_min( 451)=0.151529904E+01
       w_qmin( 451)=0.459837761E-04
       qk_min( 452)=0.151828528E+01
       w_qmin( 452)=0.479914037E-04
       qk_min( 453)=0.152127151E+01
       w_qmin( 453)=0.500740716E-04
       qk_min( 454)=0.152425775E+01
       w_qmin( 454)=0.522342430E-04
       qk_min( 455)=0.152724399E+01
       w_qmin( 455)=0.544743776E-04
       qk_min( 456)=0.153023022E+01
       w_qmin( 456)=0.567971198E-04
       qk_min( 457)=0.153321646E+01
       w_qmin( 457)=0.592052189E-04
       qk_min( 458)=0.153620269E+01
       w_qmin( 458)=0.617012423E-04
       qk_min( 459)=0.153918893E+01
       w_qmin( 459)=0.642881487E-04
       qk_min( 460)=0.154217517E+01
       w_qmin( 460)=0.669687132E-04
       qk_min( 461)=0.154516140E+01
       w_qmin( 461)=0.697458622E-04
       qk_min( 462)=0.154814764E+01
       w_qmin( 462)=0.726227948E-04
       qk_min( 463)=0.155113387E+01
       w_qmin( 463)=0.756024033E-04
       qk_min( 464)=0.155412011E+01
       w_qmin( 464)=0.786879653E-04
       qk_min( 465)=0.155710635E+01
       w_qmin( 465)=0.818828300E-04
       qk_min( 466)=0.156009258E+01
       w_qmin( 466)=0.851901657E-04
       qk_min( 467)=0.156307882E+01
       w_qmin( 467)=0.886134930E-04
       qk_min( 468)=0.156606505E+01
       w_qmin( 468)=0.921564139E-04
       qk_min( 469)=0.156905129E+01
       w_qmin( 469)=0.958223661E-04
       qk_min( 470)=0.157203753E+01
       w_qmin( 470)=0.996152064E-04
       qk_min( 471)=0.157502376E+01
       w_qmin( 471)=0.103538575E-03
       qk_min( 472)=0.157801000E+01
       w_qmin( 472)=0.107596429E-03
       qk_min( 473)=0.158099623E+01
       w_qmin( 473)=0.111792857E-03
       qk_min( 474)=0.158398247E+01
       w_qmin( 474)=0.116131714E-03
       qk_min( 475)=0.158696871E+01
       w_qmin( 475)=0.120617363E-03
       qk_min( 476)=0.158995494E+01
       w_qmin( 476)=0.125254114E-03
       qk_min( 477)=0.159294118E+01
       w_qmin( 477)=0.130046282E-03
       qk_min( 478)=0.159592741E+01
       w_qmin( 478)=0.134998333E-03
       qk_min( 479)=0.159891365E+01
       w_qmin( 479)=0.140115085E-03
       qk_min( 480)=0.160189988E+01
       w_qmin( 480)=0.145401211E-03
       qk_min( 481)=0.160488612E+01
       w_qmin( 481)=0.150861547E-03
       qk_min( 482)=0.160787236E+01
       w_qmin( 482)=0.156501094E-03
       qk_min( 483)=0.161085859E+01
       w_qmin( 483)=0.162324989E-03
       qk_min( 484)=0.161384483E+01
       w_qmin( 484)=0.168338617E-03
       qk_min( 485)=0.161683106E+01
       w_qmin( 485)=0.174547287E-03
       qk_min( 486)=0.161981730E+01
       w_qmin( 486)=0.180956528E-03
       qk_min( 487)=0.162280354E+01
       w_qmin( 487)=0.187572222E-03
       qk_min( 488)=0.162578977E+01
       w_qmin( 488)=0.194400044E-03
       qk_min( 489)=0.162877601E+01
       w_qmin( 489)=0.201446006E-03
       qk_min( 490)=0.163176224E+01
       w_qmin( 490)=0.208716558E-03
       qk_min( 491)=0.163474848E+01
       w_qmin( 491)=0.216217754E-03
       qk_min( 492)=0.163773472E+01
       w_qmin( 492)=0.223956400E-03
       qk_min( 493)=0.164072095E+01
       w_qmin( 493)=0.231938975E-03
       qk_min( 494)=0.164370719E+01
       w_qmin( 494)=0.240172480E-03
       qk_min( 495)=0.164669342E+01
       w_qmin( 495)=0.248664180E-03
       qk_min( 496)=0.164967966E+01
       w_qmin( 496)=0.257421259E-03
       qk_min( 497)=0.165266590E+01
       w_qmin( 497)=0.266451259E-03
       qk_min( 498)=0.165565213E+01
       w_qmin( 498)=0.275762188E-03
       qk_min( 499)=0.165863837E+01
       w_qmin( 499)=0.285361592E-03
       qk_min( 500)=0.166162460E+01
       w_qmin( 500)=0.295258054E-03
       qk_min( 501)=0.166461084E+01
       w_qmin( 501)=0.305460085E-03
       qk_min( 502)=0.166759708E+01
       w_qmin( 502)=0.315976296E-03
       qk_min( 503)=0.167058331E+01
       w_qmin( 503)=0.326815911E-03
       qk_min( 504)=0.167356955E+01
       w_qmin( 504)=0.337987974E-03
       qk_min( 505)=0.167655578E+01
       w_qmin( 505)=0.349502332E-03
       qk_min( 506)=0.167954202E+01
       w_qmin( 506)=0.361369081E-03
       qk_min( 507)=0.168252826E+01
       w_qmin( 507)=0.373598144E-03
       qk_min( 508)=0.168551449E+01
       w_qmin( 508)=0.386200398E-03
       qk_min( 509)=0.168850073E+01
       w_qmin( 509)=0.399186814E-03
       qk_min( 510)=0.169148696E+01
       w_qmin( 510)=0.412568599E-03
       qk_min( 511)=0.169447320E+01
       w_qmin( 511)=0.426357630E-03
       qk_min( 512)=0.169745944E+01
       w_qmin( 512)=0.440566135E-03
       qk_min( 513)=0.170044567E+01
       w_qmin( 513)=0.455206668E-03
       qk_min( 514)=0.170343191E+01
       w_qmin( 514)=0.470292258E-03
       qk_min( 515)=0.170641814E+01
       w_qmin( 515)=0.485836517E-03
       qk_min( 516)=0.170940438E+01
       w_qmin( 516)=0.501853310E-03
       qk_min( 517)=0.171239062E+01
       w_qmin( 517)=0.518357580E-03
       qk_min( 518)=0.171537685E+01
       w_qmin( 518)=0.535364141E-03
       qk_min( 519)=0.171836309E+01
       w_qmin( 519)=0.552888807E-03
       qk_min( 520)=0.172134932E+01
       w_qmin( 520)=0.570948256E-03
       qk_min( 521)=0.172433556E+01
       w_qmin( 521)=0.589558827E-03
       qk_min( 522)=0.172732180E+01
       w_qmin( 522)=0.608738722E-03
       qk_min( 523)=0.173030803E+01
       w_qmin( 523)=0.628506301E-03
       qk_min( 524)=0.173329427E+01
       w_qmin( 524)=0.648880463E-03
       qk_min( 525)=0.173628050E+01
       w_qmin( 525)=0.669881735E-03
       qk_min( 526)=0.173926674E+01
       w_qmin( 526)=0.691529804E-03
       qk_min( 527)=0.174225298E+01
       w_qmin( 527)=0.713847476E-03
       qk_min( 528)=0.174523921E+01
       w_qmin( 528)=0.736857210E-03
       qk_min( 529)=0.174822545E+01
       w_qmin( 529)=0.760582053E-03
       qk_min( 530)=0.175121168E+01
       w_qmin( 530)=0.785047285E-03
       qk_min( 531)=0.175419792E+01
       w_qmin( 531)=0.810279012E-03
       qk_min( 532)=0.175718416E+01
       w_qmin( 532)=0.836303084E-03
       qk_min( 533)=0.176017039E+01
       w_qmin( 533)=0.863148241E-03
       qk_min( 534)=0.176315663E+01
       w_qmin( 534)=0.890844941E-03
       qk_min( 535)=0.176614286E+01
       w_qmin( 535)=0.919422666E-03
       qk_min( 536)=0.176912910E+01
       w_qmin( 536)=0.948914732E-03
       qk_min( 537)=0.177211533E+01
       w_qmin( 537)=0.979354556E-03
       qk_min( 538)=0.177510157E+01
       w_qmin( 538)=0.101077768E-02
       qk_min( 539)=0.177808781E+01
       w_qmin( 539)=0.104322263E-02
       qk_min( 540)=0.178107404E+01
       w_qmin( 540)=0.107672683E-02
       qk_min( 541)=0.178406028E+01
       w_qmin( 541)=0.111133165E-02
       qk_min( 542)=0.178704651E+01
       w_qmin( 542)=0.114708143E-02
       qk_min( 543)=0.179003275E+01
       w_qmin( 543)=0.118402019E-02
       qk_min( 544)=0.179301899E+01
       w_qmin( 544)=0.122219555E-02
       qk_min( 545)=0.179600522E+01
       w_qmin( 545)=0.126165782E-02
       qk_min( 546)=0.179899146E+01
       w_qmin( 546)=0.130245871E-02
       qk_min( 547)=0.180197769E+01
       w_qmin( 547)=0.134465404E-02
       qk_min( 548)=0.180496393E+01
       w_qmin( 548)=0.138830074E-02
       qk_min( 549)=0.180795017E+01
       w_qmin( 549)=0.143345950E-02
       qk_min( 550)=0.181093640E+01
       w_qmin( 550)=0.148019468E-02
       qk_min( 551)=0.181392264E+01
       w_qmin( 551)=0.152857157E-02
       qk_min( 552)=0.181690887E+01
       w_qmin( 552)=0.157866061E-02
       qk_min( 553)=0.181989511E+01
       w_qmin( 553)=0.163053610E-02
       qk_min( 554)=0.182288135E+01
       w_qmin( 554)=0.168427316E-02
       qk_min( 555)=0.182586758E+01
       w_qmin( 555)=0.173995289E-02
       qk_min( 556)=0.182885382E+01
       w_qmin( 556)=0.179765952E-02
       qk_min( 557)=0.183184005E+01
       w_qmin( 557)=0.185747904E-02
       qk_min( 558)=0.183482629E+01
       w_qmin( 558)=0.191950230E-02
       qk_min( 559)=0.183781253E+01
       w_qmin( 559)=0.198382641E-02
       qk_min( 560)=0.184079876E+01
       w_qmin( 560)=0.205054436E-02
       qk_min( 561)=0.184378500E+01
       w_qmin( 561)=0.211976203E-02
       qk_min( 562)=0.184677123E+01
       w_qmin( 562)=0.219158092E-02
       qk_min( 563)=0.184975747E+01
       w_qmin( 563)=0.226610813E-02
       qk_min( 564)=0.185274371E+01
       w_qmin( 564)=0.234345373E-02
       qk_min( 565)=0.185572994E+01
       w_qmin( 565)=0.242372713E-02
       qk_min( 566)=0.185871618E+01
       w_qmin( 566)=0.250704165E-02
       qk_min( 567)=0.186170241E+01
       w_qmin( 567)=0.259350806E-02
       qk_min( 568)=0.186468865E+01
       w_qmin( 568)=0.268323895E-02
       qk_min( 569)=0.186767489E+01
       w_qmin( 569)=0.277634333E-02
       qk_min( 570)=0.187066112E+01
       w_qmin( 570)=0.287292872E-02
       qk_min( 571)=0.187364736E+01
       w_qmin( 571)=0.297309477E-02
       qk_min( 572)=0.187663359E+01
       w_qmin( 572)=0.307693952E-02
       qk_min( 573)=0.187961983E+01
       w_qmin( 573)=0.318454804E-02
       qk_min( 574)=0.188260607E+01
       w_qmin( 574)=0.329599762E-02
       qk_min( 575)=0.188559230E+01
       w_qmin( 575)=0.341135647E-02
       qk_min( 576)=0.188857854E+01
       w_qmin( 576)=0.353067138E-02
       qk_min( 577)=0.189156477E+01
       w_qmin( 577)=0.365398018E-02
       qk_min( 578)=0.189455101E+01
       w_qmin( 578)=0.378129690E-02
       qk_min( 579)=0.189753725E+01
       w_qmin( 579)=0.391261592E-02
       qk_min( 580)=0.190052348E+01
       w_qmin( 580)=0.404791404E-02
       qk_min( 581)=0.190350972E+01
       w_qmin( 581)=0.418713972E-02
       qk_min( 582)=0.190649595E+01
       w_qmin( 582)=0.433021438E-02
       qk_min( 583)=0.190948219E+01
       w_qmin( 583)=0.447704605E-02
       qk_min( 584)=0.191246843E+01
       w_qmin( 584)=0.462750363E-02
       qk_min( 585)=0.191545466E+01
       w_qmin( 585)=0.478143724E-02
       qk_min( 586)=0.191844090E+01
       w_qmin( 586)=0.493868030E-02
       qk_min( 587)=0.192142713E+01
       w_qmin( 587)=0.509903314E-02
       qk_min( 588)=0.192441337E+01
       w_qmin( 588)=0.526228436E-02
       qk_min( 589)=0.192739961E+01
       w_qmin( 589)=0.542820514E-02
       qk_min( 590)=0.193038584E+01
       w_qmin( 590)=0.559654903E-02
       qk_min( 591)=0.193337208E+01
       w_qmin( 591)=0.576706235E-02
       qk_min( 592)=0.193635831E+01
       w_qmin( 592)=0.593948743E-02
       qk_min( 593)=0.193934455E+01
       w_qmin( 593)=0.611355365E-02
       qk_min( 594)=0.194233079E+01
       w_qmin( 594)=0.628900004E-02
       qk_min( 595)=0.194531702E+01
       w_qmin( 595)=0.646555578E-02
       qk_min( 596)=0.194830326E+01
       w_qmin( 596)=0.664296362E-02
       qk_min( 597)=0.195128949E+01
       w_qmin( 597)=0.682098036E-02
       qk_min( 598)=0.195427573E+01
       w_qmin( 598)=0.699935510E-02
       qk_min( 599)=0.195726196E+01
       w_qmin( 599)=0.717785742E-02
       qk_min( 600)=0.196024820E+01
       w_qmin( 600)=0.735627661E-02
       qk_min( 601)=0.196323444E+01
       w_qmin( 601)=0.753439863E-02
       qk_min( 602)=0.196622067E+01
       w_qmin( 602)=0.771203465E-02
       qk_min( 603)=0.196920691E+01
       w_qmin( 603)=0.788901168E-02
       qk_min( 604)=0.197219314E+01
       w_qmin( 604)=0.806515833E-02
       qk_min( 605)=0.197517938E+01
       w_qmin( 605)=0.824033396E-02
       qk_min( 606)=0.197816562E+01
       w_qmin( 606)=0.841439593E-02
       qk_min( 607)=0.198115185E+01
       w_qmin( 607)=0.858722323E-02
       qk_min( 608)=0.198413809E+01
       w_qmin( 608)=0.875871271E-02
       qk_min( 609)=0.198712432E+01
       w_qmin( 609)=0.892875728E-02
       qk_min( 610)=0.199011056E+01
       w_qmin( 610)=0.909727135E-02
       qk_min( 611)=0.199309680E+01
       w_qmin( 611)=0.926418915E-02
       qk_min( 612)=0.199608303E+01
       w_qmin( 612)=0.942943341E-02
       qk_min( 613)=0.199906927E+01
       w_qmin( 613)=0.959294908E-02
       qk_min( 614)=0.200205550E+01
       w_qmin( 614)=0.975468924E-02
       qk_min( 615)=0.200504174E+01
       w_qmin( 615)=0.991462229E-02
       qk_min( 616)=0.200802798E+01
       w_qmin( 616)=0.100726954E-01
       qk_min( 617)=0.201101421E+01
       w_qmin( 617)=0.102288893E-01
       qk_min( 618)=0.201400045E+01
       w_qmin( 618)=0.103831838E-01
       qk_min( 619)=0.201698668E+01
       w_qmin( 619)=0.105355603E-01
       qk_min( 620)=0.201997292E+01
       w_qmin( 620)=0.106860098E-01
       qk_min( 621)=0.202295916E+01
       w_qmin( 621)=0.108345319E-01
       qk_min( 622)=0.202594539E+01
       w_qmin( 622)=0.109811024E-01
       qk_min( 623)=0.202893163E+01
       w_qmin( 623)=0.111257341E-01
       qk_min( 624)=0.203191786E+01
       w_qmin( 624)=0.112684298E-01
       qk_min( 625)=0.203490410E+01
       w_qmin( 625)=0.114091923E-01
       qk_min( 626)=0.203789034E+01
       w_qmin( 626)=0.115480445E-01
       qk_min( 627)=0.204087657E+01
       w_qmin( 627)=0.116849705E-01
       qk_min( 628)=0.204386281E+01
       w_qmin( 628)=0.118199965E-01
       qk_min( 629)=0.204684904E+01
       w_qmin( 629)=0.119531346E-01
       qk_min( 630)=0.204983528E+01
       w_qmin( 630)=0.120844006E-01
       qk_min( 631)=0.205282152E+01
       w_qmin( 631)=0.122138135E-01
       qk_min( 632)=0.205580775E+01
       w_qmin( 632)=0.123414023E-01
       qk_min( 633)=0.205879399E+01
       w_qmin( 633)=0.124671592E-01
       qk_min( 634)=0.206178022E+01
       w_qmin( 634)=0.125911196E-01
       qk_min( 635)=0.206476646E+01
       w_qmin( 635)=0.127132995E-01
       qk_min( 636)=0.206775270E+01
       w_qmin( 636)=0.128337247E-01
       qk_min( 637)=0.207073893E+01
       w_qmin( 637)=0.129524278E-01
       qk_min( 638)=0.207372517E+01
       w_qmin( 638)=0.130694050E-01
       qk_min( 639)=0.207671140E+01
       w_qmin( 639)=0.131846907E-01
       qk_min( 640)=0.207969764E+01
       w_qmin( 640)=0.132983123E-01
       qk_min( 641)=0.208268388E+01
       w_qmin( 641)=0.134102859E-01
       qk_min( 642)=0.208567011E+01
       w_qmin( 642)=0.135206395E-01
       qk_min( 643)=0.208865635E+01
       w_qmin( 643)=0.136294069E-01
       qk_min( 644)=0.209164258E+01
       w_qmin( 644)=0.137365862E-01
       qk_min( 645)=0.209462882E+01
       w_qmin( 645)=0.138422124E-01
       qk_min( 646)=0.209761506E+01
       w_qmin( 646)=0.139463108E-01
       qk_min( 647)=0.210060129E+01
       w_qmin( 647)=0.140489012E-01
       qk_min( 648)=0.210358753E+01
       w_qmin( 648)=0.141500199E-01
       qk_min( 649)=0.210657376E+01
       w_qmin( 649)=0.142496662E-01
       qk_min( 650)=0.210956000E+01
       w_qmin( 650)=0.143478715E-01
       qk_min( 651)=0.211254624E+01
       w_qmin( 651)=0.144446632E-01
       qk_min( 652)=0.211553247E+01
       w_qmin( 652)=0.145400583E-01
       qk_min( 653)=0.211851871E+01
       w_qmin( 653)=0.146340828E-01
       qk_min( 654)=0.212150494E+01
       w_qmin( 654)=0.147267691E-01
       qk_min( 655)=0.212449118E+01
       w_qmin( 655)=0.148181103E-01
       qk_min( 656)=0.212747742E+01
       w_qmin( 656)=0.149081453E-01
       qk_min( 657)=0.213046365E+01
       w_qmin( 657)=0.149968931E-01
       qk_min( 658)=0.213344989E+01
       w_qmin( 658)=0.150843719E-01
       qk_min( 659)=0.213643612E+01
       w_qmin( 659)=0.151706188E-01
       qk_min( 660)=0.213942236E+01
       w_qmin( 660)=0.152556261E-01
       qk_min( 661)=0.214240859E+01
       w_qmin( 661)=0.153394238E-01
       qk_min( 662)=0.214539483E+01
       w_qmin( 662)=0.154220366E-01
       qk_min( 663)=0.214838107E+01
       w_qmin( 663)=0.155034801E-01
       qk_min( 664)=0.215136730E+01
       w_qmin( 664)=0.155837748E-01
       qk_min( 665)=0.215435354E+01
       w_qmin( 665)=0.156629538E-01
       qk_min( 666)=0.215733977E+01
       w_qmin( 666)=0.157410074E-01
       qk_min( 667)=0.216032601E+01
       w_qmin( 667)=0.158179681E-01
       qk_min( 668)=0.216331225E+01
       w_qmin( 668)=0.158938541E-01
       qk_min( 669)=0.216629848E+01
       w_qmin( 669)=0.159686826E-01
       qk_min( 670)=0.216928472E+01
       w_qmin( 670)=0.160424840E-01
       qk_min( 671)=0.217227095E+01
       w_qmin( 671)=0.161152479E-01
       qk_min( 672)=0.217525719E+01
       w_qmin( 672)=0.161870103E-01
       qk_min( 673)=0.217824343E+01
       w_qmin( 673)=0.162577788E-01
       qk_min( 674)=0.218122966E+01
       w_qmin( 674)=0.163275787E-01
       qk_min( 675)=0.218421590E+01
       w_qmin( 675)=0.163964232E-01
       qk_min( 676)=0.218720213E+01
       w_qmin( 676)=0.164643391E-01
       qk_min( 677)=0.219018837E+01
       w_qmin( 677)=0.165313205E-01
       qk_min( 678)=0.219317461E+01
       w_qmin( 678)=0.165973881E-01
       qk_min( 679)=0.219616084E+01
       w_qmin( 679)=0.166625641E-01
       qk_min( 680)=0.219914708E+01
       w_qmin( 680)=0.167268572E-01
       qk_min( 681)=0.220213331E+01
       w_qmin( 681)=0.167903009E-01
       qk_min( 682)=0.220511955E+01
       w_qmin( 682)=0.168528775E-01
       qk_min( 683)=0.220810579E+01
       w_qmin( 683)=0.169146152E-01
       qk_min( 684)=0.221109202E+01
       w_qmin( 684)=0.169755246E-01
       qk_min( 685)=0.221407826E+01
       w_qmin( 685)=0.170356213E-01
       qk_min( 686)=0.221706449E+01
       w_qmin( 686)=0.170949171E-01
       qk_min( 687)=0.222005073E+01
       w_qmin( 687)=0.171534415E-01
       qk_min( 688)=0.222303697E+01
       w_qmin( 688)=0.172111728E-01
       qk_min( 689)=0.222602320E+01
       w_qmin( 689)=0.172681398E-01
       qk_min( 690)=0.222900944E+01
       w_qmin( 690)=0.173243582E-01
       qk_min( 691)=0.223199567E+01
       w_qmin( 691)=0.173798301E-01
       qk_min( 692)=0.223498191E+01
       w_qmin( 692)=0.174345971E-01
       qk_min( 693)=0.223796815E+01
       w_qmin( 693)=0.174886299E-01
       qk_min( 694)=0.224095438E+01
       w_qmin( 694)=0.175419632E-01
       qk_min( 695)=0.224394062E+01
       w_qmin( 695)=0.175946025E-01
       qk_min( 696)=0.224692685E+01
       w_qmin( 696)=0.176465629E-01
       qk_min( 697)=0.224991309E+01
       w_qmin( 697)=0.176978608E-01
       qk_min( 698)=0.225289933E+01
       w_qmin( 698)=0.177485202E-01
       qk_min( 699)=0.225588556E+01
       w_qmin( 699)=0.177985286E-01
       qk_min( 700)=0.225887180E+01
       w_qmin( 700)=0.178478959E-01
       qk_min( 701)=0.226185803E+01
       w_qmin( 701)=0.178966537E-01
       qk_min( 702)=0.226484427E+01
       w_qmin( 702)=0.179448012E-01
       qk_min( 703)=0.226783051E+01
       w_qmin( 703)=0.179923681E-01
       qk_min( 704)=0.227081674E+01
       w_qmin( 704)=0.180393329E-01
       qk_min( 705)=0.227380298E+01
       w_qmin( 705)=0.180857222E-01
       qk_min( 706)=0.227678921E+01
       w_qmin( 706)=0.181315403E-01
       qk_min( 707)=0.227977545E+01
       w_qmin( 707)=0.181767980E-01
       qk_min( 708)=0.228276169E+01
       w_qmin( 708)=0.182215072E-01
       qk_min( 709)=0.228574792E+01
       w_qmin( 709)=0.182656858E-01
       qk_min( 710)=0.228873416E+01
       w_qmin( 710)=0.183093106E-01
       qk_min( 711)=0.229172039E+01
       w_qmin( 711)=0.183524032E-01
       qk_min( 712)=0.229470663E+01
       w_qmin( 712)=0.183949727E-01
       qk_min( 713)=0.229769287E+01
       w_qmin( 713)=0.184370255E-01
       qk_min( 714)=0.230067910E+01
       w_qmin( 714)=0.184785801E-01
       qk_min( 715)=0.230366534E+01
       w_qmin( 715)=0.185196119E-01
       qk_min( 716)=0.230665157E+01
       w_qmin( 716)=0.185601483E-01
       qk_min( 717)=0.230963781E+01
       w_qmin( 717)=0.186001946E-01
       qk_min( 718)=0.231262405E+01
       w_qmin( 718)=0.186397515E-01
       qk_min( 719)=0.231561028E+01
       w_qmin( 719)=0.186788294E-01
       qk_min( 720)=0.231859652E+01
       w_qmin( 720)=0.187174615E-01
       qk_min( 721)=0.232158275E+01
       w_qmin( 721)=0.187556096E-01
       qk_min( 722)=0.232456899E+01
       w_qmin( 722)=0.187933082E-01
       qk_min( 723)=0.232755522E+01
       w_qmin( 723)=0.188305624E-01
       qk_min( 724)=0.233054146E+01
       w_qmin( 724)=0.188673824E-01
       qk_min( 725)=0.233352770E+01
       w_qmin( 725)=0.189037949E-01
       qk_min( 726)=0.233651393E+01
       w_qmin( 726)=0.189397778E-01
       qk_min( 727)=0.233950017E+01
       w_qmin( 727)=0.189753492E-01
       qk_min( 728)=0.234248640E+01
       w_qmin( 728)=0.190105254E-01
       qk_min( 729)=0.234547264E+01
       w_qmin( 729)=0.190453153E-01
       qk_min( 730)=0.234845888E+01
       w_qmin( 730)=0.190797269E-01
       qk_min( 731)=0.235144511E+01
       w_qmin( 731)=0.191137845E-01
       qk_min( 732)=0.235443135E+01
       w_qmin( 732)=0.191474604E-01
       qk_min( 733)=0.235741758E+01
       w_qmin( 733)=0.191807822E-01
       qk_min( 734)=0.236040382E+01
       w_qmin( 734)=0.192137587E-01
       qk_min( 735)=0.236339006E+01
       w_qmin( 735)=0.192463875E-01
       qk_min( 736)=0.236637629E+01
       w_qmin( 736)=0.192786988E-01
       qk_min( 737)=0.236936253E+01
       w_qmin( 737)=0.193106570E-01
       qk_min( 738)=0.237234876E+01
       w_qmin( 738)=0.193422911E-01
       qk_min( 739)=0.237533500E+01
       w_qmin( 739)=0.193736007E-01
       qk_min( 740)=0.237832124E+01
       w_qmin( 740)=0.194045867E-01
       qk_min( 741)=0.238130747E+01
       w_qmin( 741)=0.194352575E-01
       qk_min( 742)=0.238429371E+01
       w_qmin( 742)=0.194656293E-01
       qk_min( 743)=0.238727994E+01
       w_qmin( 743)=0.194956639E-01
       qk_min( 744)=0.239026618E+01
       w_qmin( 744)=0.195253908E-01
       qk_min( 745)=0.239325242E+01
       w_qmin( 745)=0.195548023E-01
       qk_min( 746)=0.239623865E+01
       w_qmin( 746)=0.195839100E-01
       qk_min( 747)=0.239922489E+01
       w_qmin( 747)=0.196127197E-01
       qk_min( 748)=0.240221112E+01
       w_qmin( 748)=0.196412099E-01
       qk_min( 749)=0.240519736E+01
       w_qmin( 749)=0.196693885E-01
       qk_min( 750)=0.240818360E+01
       w_qmin( 750)=0.196972667E-01
       qk_min( 751)=0.241116983E+01
       w_qmin( 751)=0.197248463E-01
       qk_min( 752)=0.241415607E+01
       w_qmin( 752)=0.197521245E-01
       qk_min( 753)=0.241714230E+01
       w_qmin( 753)=0.197791298E-01
       qk_min( 754)=0.242012854E+01
       w_qmin( 754)=0.198058300E-01
       qk_min( 755)=0.242311478E+01
       w_qmin( 755)=0.198322465E-01
       qk_min( 756)=0.242610101E+01
       w_qmin( 756)=0.198583839E-01
       qk_min( 757)=0.242908725E+01
       w_qmin( 757)=0.198842496E-01
       qk_min( 758)=0.243207348E+01
       w_qmin( 758)=0.199098720E-01
       qk_min( 759)=0.243505972E+01
       w_qmin( 759)=0.199352178E-01
       qk_min( 760)=0.243804596E+01
       w_qmin( 760)=0.199603180E-01
       qk_min( 761)=0.244103219E+01
       w_qmin( 761)=0.199851764E-01
       qk_min( 762)=0.244401843E+01
       w_qmin( 762)=0.200098100E-01
       qk_min( 763)=0.244700466E+01
       w_qmin( 763)=0.200342154E-01
       qk_min( 764)=0.244999090E+01
       w_qmin( 764)=0.200584302E-01
       qk_min( 765)=0.245297714E+01
       w_qmin( 765)=0.200824225E-01
       qk_min( 766)=0.245596337E+01
       w_qmin( 766)=0.201062163E-01
       qk_min( 767)=0.245894961E+01
       w_qmin( 767)=0.201298320E-01
       qk_min( 768)=0.246193584E+01
       w_qmin( 768)=0.201532643E-01
       qk_min( 769)=0.246492208E+01
       w_qmin( 769)=0.201765465E-01
       qk_min( 770)=0.246790832E+01
       w_qmin( 770)=0.201996489E-01
       qk_min( 771)=0.247089455E+01
       w_qmin( 771)=0.202225967E-01
       qk_min( 772)=0.247388079E+01
       w_qmin( 772)=0.202453918E-01
       qk_min( 773)=0.247686702E+01
       w_qmin( 773)=0.202680411E-01
       qk_min( 774)=0.247985326E+01
       w_qmin( 774)=0.202905521E-01
       qk_min( 775)=0.248283950E+01
       w_qmin( 775)=0.203129378E-01
       qk_min( 776)=0.248582573E+01
       w_qmin( 776)=0.203351772E-01
       qk_min( 777)=0.248881197E+01
       w_qmin( 777)=0.203572820E-01
       qk_min( 778)=0.249179820E+01
       w_qmin( 778)=0.203792539E-01
       qk_min( 779)=0.249478444E+01
       w_qmin( 779)=0.204011005E-01
       qk_min( 780)=0.249777068E+01
       w_qmin( 780)=0.204228238E-01
       qk_min( 781)=0.250075691E+01
       w_qmin( 781)=0.204444094E-01
       qk_min( 782)=0.250374315E+01
       w_qmin( 782)=0.204658525E-01
       qk_min( 783)=0.250672938E+01
       w_qmin( 783)=0.204871663E-01
       qk_min( 784)=0.250971562E+01
       w_qmin( 784)=0.205083419E-01
       qk_min( 785)=0.251270185E+01
       w_qmin( 785)=0.205293734E-01
       qk_min( 786)=0.251568809E+01
       w_qmin( 786)=0.205502826E-01
       qk_min( 787)=0.251867433E+01
       w_qmin( 787)=0.205710219E-01
       qk_min( 788)=0.252166056E+01
       w_qmin( 788)=0.205916076E-01
       qk_min( 789)=0.252464680E+01
       w_qmin( 789)=0.206120492E-01
       qk_min( 790)=0.252763303E+01
       w_qmin( 790)=0.206323276E-01
       qk_min( 791)=0.253061927E+01
       w_qmin( 791)=0.206524614E-01
       qk_min( 792)=0.253360551E+01
       w_qmin( 792)=0.206724200E-01
       qk_min( 793)=0.253659174E+01
       w_qmin( 793)=0.206922093E-01
       qk_min( 794)=0.253957798E+01
       w_qmin( 794)=0.207118259E-01
       qk_min( 795)=0.254256421E+01
       w_qmin( 795)=0.207312799E-01
       qk_min( 796)=0.254555045E+01
       w_qmin( 796)=0.207505596E-01
       qk_min( 797)=0.254853669E+01
       w_qmin( 797)=0.207696884E-01
       qk_min( 798)=0.255152292E+01
       w_qmin( 798)=0.207886275E-01
       qk_min( 799)=0.255450916E+01
       w_qmin( 799)=0.208073986E-01
       qk_min( 800)=0.255749539E+01
       w_qmin( 800)=0.208260043E-01
       qk_min( 801)=0.256048163E+01
       w_qmin( 801)=0.208444481E-01
       qk_min( 802)=0.256346787E+01
       w_qmin( 802)=0.208627364E-01
       qk_min( 803)=0.256645410E+01
       w_qmin( 803)=0.208808634E-01
       qk_min( 804)=0.256944034E+01
       w_qmin( 804)=0.208988348E-01
       qk_min( 805)=0.257242657E+01
       w_qmin( 805)=0.209166573E-01
       qk_min( 806)=0.257541281E+01
       w_qmin( 806)=0.209343399E-01
       qk_min( 807)=0.257839905E+01
       w_qmin( 807)=0.209518904E-01
       qk_min( 808)=0.258138528E+01
       w_qmin( 808)=0.209693308E-01
       qk_min( 809)=0.258437152E+01
       w_qmin( 809)=0.209866286E-01
       qk_min( 810)=0.258735775E+01
       w_qmin( 810)=0.210038226E-01
       qk_min( 811)=0.259034399E+01
       w_qmin( 811)=0.210209134E-01
       qk_min( 812)=0.259333023E+01
       w_qmin( 812)=0.210379105E-01
       qk_min( 813)=0.259631646E+01
       w_qmin( 813)=0.210548386E-01
       qk_min( 814)=0.259930270E+01
       w_qmin( 814)=0.210716804E-01
       qk_min( 815)=0.260228893E+01
       w_qmin( 815)=0.210884458E-01
       qk_min( 816)=0.260527517E+01
       w_qmin( 816)=0.211051674E-01
       qk_min( 817)=0.260826141E+01
       w_qmin( 817)=0.211218344E-01
       qk_min( 818)=0.261124764E+01
       w_qmin( 818)=0.211384647E-01
       qk_min( 819)=0.261423388E+01
       w_qmin( 819)=0.211550795E-01
       qk_min( 820)=0.261722011E+01
       w_qmin( 820)=0.211716577E-01
       qk_min( 821)=0.262020635E+01
       w_qmin( 821)=0.211882200E-01
       qk_min( 822)=0.262319259E+01
       w_qmin( 822)=0.212047806E-01
       qk_min( 823)=0.262617882E+01
       w_qmin( 823)=0.212213377E-01
       qk_min( 824)=0.262916506E+01
       w_qmin( 824)=0.212379084E-01
       qk_min( 825)=0.263215129E+01
       w_qmin( 825)=0.212544766E-01
       qk_min( 826)=0.263513753E+01
       w_qmin( 826)=0.212710628E-01
       qk_min( 827)=0.263812377E+01
       w_qmin( 827)=0.212876624E-01
       qk_min( 828)=0.264111000E+01
       w_qmin( 828)=0.213042839E-01
       qk_min( 829)=0.264409624E+01
       w_qmin( 829)=0.213209307E-01
       qk_min( 830)=0.264708247E+01
       w_qmin( 830)=0.213376112E-01
       qk_min( 831)=0.265006871E+01
       w_qmin( 831)=0.213542997E-01
       qk_min( 832)=0.265305495E+01
       w_qmin( 832)=0.213710178E-01
       qk_min( 833)=0.265604118E+01
       w_qmin( 833)=0.213877536E-01
       qk_min( 834)=0.265902742E+01
       w_qmin( 834)=0.214045146E-01
       qk_min( 835)=0.266201365E+01
       w_qmin( 835)=0.214212992E-01
       qk_min( 836)=0.266499989E+01
       w_qmin( 836)=0.214380876E-01
       qk_min( 837)=0.266798613E+01
       w_qmin( 837)=0.214548824E-01
       qk_min( 838)=0.267097236E+01
       w_qmin( 838)=0.214716909E-01
       qk_min( 839)=0.267395860E+01
       w_qmin( 839)=0.214884922E-01
       qk_min( 840)=0.267694483E+01
       w_qmin( 840)=0.215052997E-01
       qk_min( 841)=0.267993107E+01
       w_qmin( 841)=0.215221154E-01
       qk_min( 842)=0.268291730E+01
       w_qmin( 842)=0.215388959E-01
       qk_min( 843)=0.268590354E+01
       w_qmin( 843)=0.215556562E-01
       qk_min( 844)=0.268888978E+01
       w_qmin( 844)=0.215723984E-01
       qk_min( 845)=0.269187601E+01
       w_qmin( 845)=0.215891019E-01
       qk_min( 846)=0.269486225E+01
       w_qmin( 846)=0.216057878E-01
       qk_min( 847)=0.269784848E+01
       w_qmin( 847)=0.216224100E-01
       qk_min( 848)=0.270083472E+01
       w_qmin( 848)=0.216389850E-01
       qk_min( 849)=0.270382096E+01
       w_qmin( 849)=0.216555078E-01
       qk_min( 850)=0.270680719E+01
       w_qmin( 850)=0.216719646E-01
       qk_min( 851)=0.270979343E+01
       w_qmin( 851)=0.216883554E-01
       qk_min( 852)=0.271277966E+01
       w_qmin( 852)=0.217046927E-01
       qk_min( 853)=0.271576590E+01
       w_qmin( 853)=0.217209473E-01
       qk_min( 854)=0.271875214E+01
       w_qmin( 854)=0.217371074E-01
       qk_min( 855)=0.272173837E+01
       w_qmin( 855)=0.217531950E-01
       qk_min( 856)=0.272472461E+01
       w_qmin( 856)=0.217692009E-01
       qk_min( 857)=0.272771084E+01
       w_qmin( 857)=0.217851223E-01
       qk_min( 858)=0.273069708E+01
       w_qmin( 858)=0.218009415E-01
       qk_min( 859)=0.273368332E+01
       w_qmin( 859)=0.218166714E-01
       qk_min( 860)=0.273666955E+01
       w_qmin( 860)=0.218323080E-01
       qk_min( 861)=0.273965579E+01
       w_qmin( 861)=0.218478544E-01
       qk_min( 862)=0.274264202E+01
       w_qmin( 862)=0.218633094E-01
       qk_min( 863)=0.274562826E+01
       w_qmin( 863)=0.218786805E-01
       qk_min( 864)=0.274861450E+01
       w_qmin( 864)=0.218939532E-01
       qk_min( 865)=0.275160073E+01
       w_qmin( 865)=0.219091338E-01
       qk_min( 866)=0.275458697E+01
       w_qmin( 866)=0.219242331E-01
       qk_min( 867)=0.275757320E+01
       w_qmin( 867)=0.219392516E-01
       qk_min( 868)=0.276055944E+01
       w_qmin( 868)=0.219542071E-01
       qk_min( 869)=0.276354568E+01
       w_qmin( 869)=0.219690692E-01
       qk_min( 870)=0.276653191E+01
       w_qmin( 870)=0.219838618E-01
       qk_min( 871)=0.276951815E+01
       w_qmin( 871)=0.219985889E-01
       qk_min( 872)=0.277250438E+01
       w_qmin( 872)=0.220132414E-01
       qk_min( 873)=0.277549062E+01
       w_qmin( 873)=0.220278473E-01
       qk_min( 874)=0.277847686E+01
       w_qmin( 874)=0.220424160E-01
       qk_min( 875)=0.278146309E+01
       w_qmin( 875)=0.220569198E-01
       qk_min( 876)=0.278444933E+01
       w_qmin( 876)=0.220713889E-01
       qk_min( 877)=0.278743556E+01
       w_qmin( 877)=0.220858345E-01
       qk_min( 878)=0.279042180E+01
       w_qmin( 878)=0.221002366E-01
       qk_min( 879)=0.279340804E+01
       w_qmin( 879)=0.221146318E-01
       qk_min( 880)=0.279639427E+01
       w_qmin( 880)=0.221290153E-01
       qk_min( 881)=0.279938051E+01
       w_qmin( 881)=0.221433895E-01
       qk_min( 882)=0.280236674E+01
       w_qmin( 882)=0.221577555E-01
       qk_min( 883)=0.280535298E+01
       w_qmin( 883)=0.221721332E-01
       qk_min( 884)=0.280833922E+01
       w_qmin( 884)=0.221865214E-01
       qk_min( 885)=0.281132545E+01
       w_qmin( 885)=0.222009499E-01
       qk_min( 886)=0.281431169E+01
       w_qmin( 886)=0.222154001E-01
       qk_min( 887)=0.281729792E+01
       w_qmin( 887)=0.222298815E-01
       qk_min( 888)=0.282028416E+01
       w_qmin( 888)=0.222444115E-01
       qk_min( 889)=0.282327040E+01
       w_qmin( 889)=0.222589981E-01
       qk_min( 890)=0.282625663E+01
       w_qmin( 890)=0.222736348E-01
       qk_min( 891)=0.282924287E+01
       w_qmin( 891)=0.222883526E-01
       qk_min( 892)=0.283222910E+01
       w_qmin( 892)=0.223031301E-01
       qk_min( 893)=0.283521534E+01
       w_qmin( 893)=0.223179987E-01
       qk_min( 894)=0.283820158E+01
       w_qmin( 894)=0.223329387E-01
       qk_min( 895)=0.284118781E+01
       w_qmin( 895)=0.223479838E-01
       qk_min( 896)=0.284417405E+01
       w_qmin( 896)=0.223631308E-01
       qk_min( 897)=0.284716028E+01
       w_qmin( 897)=0.223783701E-01
       qk_min( 898)=0.285014652E+01
       w_qmin( 898)=0.223937147E-01
       qk_min( 899)=0.285313276E+01
       w_qmin( 899)=0.224091676E-01
       qk_min( 900)=0.285611899E+01
       w_qmin( 900)=0.224247395E-01
       qk_min( 901)=0.285910523E+01
       w_qmin( 901)=0.224404371E-01
       qk_min( 902)=0.286209146E+01
       w_qmin( 902)=0.224562550E-01
       qk_min( 903)=0.286507770E+01
       w_qmin( 903)=0.224721864E-01
       qk_min( 904)=0.286806393E+01
       w_qmin( 904)=0.224882382E-01
       qk_min( 905)=0.287105017E+01
       w_qmin( 905)=0.225044056E-01
       qk_min( 906)=0.287403641E+01
       w_qmin( 906)=0.225207191E-01
       qk_min( 907)=0.287702264E+01
       w_qmin( 907)=0.225371498E-01
       qk_min( 908)=0.288000888E+01
       w_qmin( 908)=0.225536978E-01
       qk_min( 909)=0.288299511E+01
       w_qmin( 909)=0.225703591E-01
       qk_min( 910)=0.288598135E+01
       w_qmin( 910)=0.225871495E-01
       qk_min( 911)=0.288896759E+01
       w_qmin( 911)=0.226040409E-01
       qk_min( 912)=0.289195382E+01
       w_qmin( 912)=0.226210497E-01
       qk_min( 913)=0.289494006E+01
       w_qmin( 913)=0.226381881E-01
       qk_min( 914)=0.289792629E+01
       w_qmin( 914)=0.226554163E-01
       qk_min( 915)=0.290091253E+01
       w_qmin( 915)=0.226727412E-01
       qk_min( 916)=0.290389877E+01
       w_qmin( 916)=0.226901695E-01
       qk_min( 917)=0.290688500E+01
       w_qmin( 917)=0.227076836E-01
       qk_min( 918)=0.290987124E+01
       w_qmin( 918)=0.227252900E-01
       qk_min( 919)=0.291285747E+01
       w_qmin( 919)=0.227429806E-01
       qk_min( 920)=0.291584371E+01
       w_qmin( 920)=0.227607376E-01
       qk_min( 921)=0.291882995E+01
       w_qmin( 921)=0.227785663E-01
       qk_min( 922)=0.292181618E+01
       w_qmin( 922)=0.227964612E-01
       qk_min( 923)=0.292480242E+01
       w_qmin( 923)=0.228144070E-01
       qk_min( 924)=0.292778865E+01
       w_qmin( 924)=0.228324315E-01
       qk_min( 925)=0.293077489E+01
       w_qmin( 925)=0.228504829E-01
       qk_min( 926)=0.293376113E+01
       w_qmin( 926)=0.228685835E-01
       qk_min( 927)=0.293674736E+01
       w_qmin( 927)=0.228867053E-01
       qk_min( 928)=0.293973360E+01
       w_qmin( 928)=0.229048652E-01
       qk_min( 929)=0.294271983E+01
       w_qmin( 929)=0.229230603E-01
       qk_min( 930)=0.294570607E+01
       w_qmin( 930)=0.229412501E-01
       qk_min( 931)=0.294869231E+01
       w_qmin( 931)=0.229594679E-01
       qk_min( 932)=0.295167854E+01
       w_qmin( 932)=0.229776894E-01
       qk_min( 933)=0.295466478E+01
       w_qmin( 933)=0.229959141E-01
       qk_min( 934)=0.295765101E+01
       w_qmin( 934)=0.230141333E-01
       qk_min( 935)=0.296063725E+01
       w_qmin( 935)=0.230323669E-01
       qk_min( 936)=0.296362349E+01
       w_qmin( 936)=0.230505769E-01
       qk_min( 937)=0.296660972E+01
       w_qmin( 937)=0.230687806E-01
       qk_min( 938)=0.296959596E+01
       w_qmin( 938)=0.230869536E-01
       qk_min( 939)=0.297258219E+01
       w_qmin( 939)=0.231051289E-01
       qk_min( 940)=0.297556843E+01
       w_qmin( 940)=0.231232707E-01
       qk_min( 941)=0.297855467E+01
       w_qmin( 941)=0.231413792E-01
       qk_min( 942)=0.298154090E+01
       w_qmin( 942)=0.231594641E-01
       qk_min( 943)=0.298452714E+01
       w_qmin( 943)=0.231775183E-01
       qk_min( 944)=0.298751337E+01
       w_qmin( 944)=0.231955454E-01
       qk_min( 945)=0.299049961E+01
       w_qmin( 945)=0.232135373E-01
       qk_min( 946)=0.299348585E+01
       w_qmin( 946)=0.232315001E-01
       qk_min( 947)=0.299647208E+01
       w_qmin( 947)=0.232494192E-01
       qk_min( 948)=0.299945832E+01
       w_qmin( 948)=0.232673042E-01
       qk_min( 949)=0.300244455E+01
       w_qmin( 949)=0.232851580E-01
       qk_min( 950)=0.300543079E+01
       w_qmin( 950)=0.233029694E-01
       qk_min( 951)=0.300841703E+01
       w_qmin( 951)=0.233207514E-01
       qk_min( 952)=0.301140326E+01
       w_qmin( 952)=0.233384939E-01
       qk_min( 953)=0.301438950E+01
       w_qmin( 953)=0.233561936E-01
       qk_min( 954)=0.301737573E+01
       w_qmin( 954)=0.233738681E-01
       qk_min( 955)=0.302036197E+01
       w_qmin( 955)=0.233914984E-01
       qk_min( 956)=0.302334821E+01
       w_qmin( 956)=0.234091148E-01
       qk_min( 957)=0.302633444E+01
       w_qmin( 957)=0.234266923E-01
       qk_min( 958)=0.302932068E+01
       w_qmin( 958)=0.234442482E-01
       qk_min( 959)=0.303230691E+01
       w_qmin( 959)=0.234617565E-01
       qk_min( 960)=0.303529315E+01
       w_qmin( 960)=0.234792575E-01
       qk_min( 961)=0.303827939E+01
       w_qmin( 961)=0.234967327E-01
       qk_min( 962)=0.304126562E+01
       w_qmin( 962)=0.235141912E-01
       qk_min( 963)=0.304425186E+01
       w_qmin( 963)=0.235316186E-01
       qk_min( 964)=0.304723809E+01
       w_qmin( 964)=0.235490383E-01
       qk_min( 965)=0.305022433E+01
       w_qmin( 965)=0.235664448E-01
       qk_min( 966)=0.305321056E+01
       w_qmin( 966)=0.235838473E-01
       qk_min( 967)=0.305619680E+01
       w_qmin( 967)=0.236012407E-01
       qk_min( 968)=0.305918304E+01
       w_qmin( 968)=0.236186374E-01
       qk_min( 969)=0.306216927E+01
       w_qmin( 969)=0.236360168E-01
       qk_min( 970)=0.306515551E+01
       w_qmin( 970)=0.236534103E-01
       qk_min( 971)=0.306814174E+01
       w_qmin( 971)=0.236708031E-01
       qk_min( 972)=0.307112798E+01
       w_qmin( 972)=0.236882071E-01
       qk_min( 973)=0.307411422E+01
       w_qmin( 973)=0.237056372E-01
       qk_min( 974)=0.307710045E+01
       w_qmin( 974)=0.237230797E-01
       qk_min( 975)=0.308008669E+01
       w_qmin( 975)=0.237405326E-01
       qk_min( 976)=0.308307292E+01
       w_qmin( 976)=0.237580036E-01
       qk_min( 977)=0.308605916E+01
       w_qmin( 977)=0.237755219E-01
       qk_min( 978)=0.308904540E+01
       w_qmin( 978)=0.237930688E-01
       qk_min( 979)=0.309203163E+01
       w_qmin( 979)=0.238106525E-01
       qk_min( 980)=0.309501787E+01
       w_qmin( 980)=0.238282891E-01
       qk_min( 981)=0.309800410E+01
       w_qmin( 981)=0.238459565E-01
       qk_min( 982)=0.310099034E+01
       w_qmin( 982)=0.238636744E-01
       qk_min( 983)=0.310397658E+01
       w_qmin( 983)=0.238814546E-01
       qk_min( 984)=0.310696281E+01
       w_qmin( 984)=0.238993028E-01
       qk_min( 985)=0.310994905E+01
       w_qmin( 985)=0.239171986E-01
       qk_min( 986)=0.311293528E+01
       w_qmin( 986)=0.239351669E-01
       qk_min( 987)=0.311592152E+01
       w_qmin( 987)=0.239532177E-01
       qk_min( 988)=0.311890776E+01
       w_qmin( 988)=0.239713333E-01
       qk_min( 989)=0.312189399E+01
       w_qmin( 989)=0.239895409E-01
       qk_min( 990)=0.312488023E+01
       w_qmin( 990)=0.240078380E-01
       qk_min( 991)=0.312786646E+01
       w_qmin( 991)=0.240262264E-01
       qk_min( 992)=0.313085270E+01
       w_qmin( 992)=0.240447070E-01
       qk_min( 993)=0.313383894E+01
       w_qmin( 993)=0.240632838E-01
       qk_min( 994)=0.313682517E+01
       w_qmin( 994)=0.240819687E-01
       qk_min( 995)=0.313981141E+01
       w_qmin( 995)=0.241007713E-01
       qk_min( 996)=0.314279764E+01
       w_qmin( 996)=0.241196804E-01
       qk_min( 997)=0.314578388E+01
       w_qmin( 997)=0.241387052E-01
       qk_min( 998)=0.314877012E+01
       w_qmin( 998)=0.241578507E-01
       qk_min( 999)=0.315175635E+01
       w_qmin( 999)=0.241771224E-01
       qk_min(1000)=0.315474259E+01
       w_qmin(1000)=0.241965235E-01
       qk_min(1001)=0.315772882E+01
       w_qmin(1001)=0.242160684E-01
       END
      SUBROUTINE CH3PISET(JJ)
C information on 3 pion sub-channel under construction obtained
C J3PI=1 means 3 prong
C J3PI=2 means 1 prong
C to be initialized in routine DPHSAA of tauola.f
      COMMON /CHANOPT/ J3PI
      INTEGER          J3PI
!
      BIND(C,NAME='tauola_chanopt') /CHANOPT/
!
      J3PI=JJ
      end
      SUBROUTINE CH3PIGET(JJ)
C information on 3 pion sub-channel under construction obtained
C J3PI=1 means 3 prong
C J3PI=2 means 1 prong
C to be initialized in routine DPHSAA of tauola.f
      COMMON /CHANOPT/ J3PI
      INTEGER          J3PI
!
      BIND(C,NAME='tauola_chanopt') /CHANOPT/
!
      IF (J3PI.EQ.1.OR.J3PI.EQ.2) THEN
        JJ=J3PI
      ELSE
       CALL TAUOLA_PRINT(
     $   'FROM value_parameter.f CH3PIGET, wrong J3PI='//CHAR(0))
!      write(*,*) 'FROM value_parameter.f CH3PIGET, wrong J3PI=',J3PI
       CALL TAUOLA_STOP()
      ENDIF
      end
      SUBROUTINE GETFF2PIRHO(JJ)
      IMPLICIT NONE
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      INTEGER JJ
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      JJ = FF2PIRHO
      END
      FUNCTION COEFrr(I,J)
C clebsh gordan (or so ...)  coefs for 3 scalar final states
      implicit none
C TAUOLA RChL COEF(I,J) =  COEFr(I,J)
      REAL COEFr(1:5,0:7)
      REAL COEFrr
      DATA PI /3.141592653589793238462643/
      REAL PI
      DATA ICONT /0/
      INTEGER ICONT
      INTEGER I,J
      REAL FPIr
C initialization of FPI matrix defined in ...
C FPIc is to be used with cleo initialization
C FPIr is to be used with RChL initialization
C actual choice is made in ???
      DATA  FPIr /92.4E-3/
C initialization of COEF matrix defined in ...
C COEFc is to be used with cleo initialization
C COEFr is to be used with RChL initialization
      IF (ICONT.EQ.0) THEN
       ICONT=1
C
C********* COEFr(I,J) *******
       COEFr(1,0)= 1.
       COEFr(2,0)= -1.
       COEFr(3,0)= 0.
       COEFr(4,0)= 1.
       COEFr(5,0)= 0.
       COEFr(1,1)= 1.
       COEFr(2,1)= -1.
       COEFr(3,1)= 0.
       COEFr(4,1)= 1.
       COEFr(5,1)= 1.
C
       COEFr(1,2)=1.
       COEFr(2,2)= -1.
       COEFr(3,2)= 0.0
       COEFr(4,2)= 1.
       COEFr(5,2)=1.
C
       COEFr(1,3)= 0.
       COEFr(2,3)= 1.
       COEFr(3,3)= -1.
       COEFr(4,3)= 1.
       COEFr(5,3)= - 1.
C
       COEFr(1,4)= 1.0/SQRT(2.)/3.0
       COEFr(2,4)=-1.0/SQRT(2.)/3.0
       COEFr(3,4)= 0.0
       COEFr(4,4)= 0.0
       COEFr(5,4)= 0.0
C
       COEFr(1,5)=-SQRT(2.)/3.0
       COEFr(2,5)= SQRT(2.)/3.0
       COEFr(3,5)= 0.0
       COEFr(4,5)= 0.0
       COEFr(5,5)=-SQRT(2.)
C
       COEFr(1,6)= 1./3.
       COEFr(2,6)=-2./3.
       COEFr(3,6)= 2./3.
       COEFr(4,6)= 0.0
       COEFr(5,6)=-2.0
C
       COEFr(1,7)= 0.0
       COEFr(2,7)= 0.0
       COEFr(3,7)= 0.0
       COEFr(4,7)= 0.0
       COEFr(5,7)=-SQRT(2.0/3.0)
      ENDIF
      COEFrr=COEFr(I,J)
      END
      subroutine rchl_parameters(KAK)
      implicit none
C==============================================================================
C  Initialization, of '../parameter.inc' common block group
C
C  KAK may be equal to JAK of TAUOLA namespace, but it is not always the case
C  Hard-coded  fit parameters:
C  rho, rhoprime, f2(1275), f0(1186), sigma(made up!)
C  The value of both the mass and width of resonances are taken
C  from  fit to ALEPH data  (ref [1], Set 1)
C  References: [1] arXiv: 0911.4436 [hep-ph] D.  Gomez Dumm et al
C              [2] arXiv: 0911.2640 [hep-ph] D.  Gomez Dumm et al.
C              [3] P Roig, talk PhiPsi2011, Novosibirsk
C              [4] arXiv:0807.4883 [hep-ph] Diogo R. Boito et al.
C              [5] arXiv:0803.1786 [hep-ph] M. Jamin et al.
C  WARNING:    some of parameters require RERUN of da1wid_tot_rho1_gauss.f
C              pretabulating Q dependent a1 width,
C              directory RChL-currents/tabler/a1
C==============================================================================
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      INTEGER KAK
      DATA IWARM/0/
      INTEGER IWARM
      INTEGER          J3PI
      COMMON /CHANOPT/ J3PI
!
      BIND(C,NAME='tauola_chanopt') /CHANOPT/
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
         IF(KAK.EQ.4) THEN
C  /MASS_RES/; resonances parameters initialization:
C                           ! at present only for two pion mode non-default
c                           ! values are used:
      mro = 0.77554d0
      mrho1 = 1.453d0
      grho1 = 0.50155D0
c   /PAR_RHOPRIME/; parameters of rho' and rho''
C                   used for 2 pion form factor, reference [3]
      COEF_GA =  0.14199D0
      COEF_DE = -0.12623D0
      phi_1   = -0.17377D0
      phi_2   =  0.27632D0
      grho2   =  0.41786D0
      mrho2   =  1.8105d0
         ELSE IF(KAK.EQ.5) THEN
      MRO   = 0.771849d0      !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MRHO1 = 1.35d0     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GRHO1 = 0.448379d0 !0.473287d0       !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      ELSE
      MRO   = 0.775     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MRHO1 = 1.465     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GRHO1 = 0.4       !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
c   /PAR_RHOPRIME/; parameters of rho' and rho''
C                   used for 2 kaon form factor, reference [3]
c                   FOR THE MOMENT THEIR NUMERICAL VALUES COINCIDE WITH
c                   ONES FOR THE TWO PION MODE !!!!
      COEF_GA =  0.14199D0
      COEF_DE = -0.12623D0
      phi_1   = -0.17377D0
      phi_2   =  0.27632D0
      grho2   =  0.41786D0
      mrho2   =  1.8105d0
         ENDIF
        IF(KAK.EQ.70)  THEN    ! non default values to be used
                               ! for KPI MODE NO FSR INTERACTION
c   /PAR_KPI/; parameters for Kpi mode, reference [4], table 4, row2
      MKST        = 0.943d0
      MKSTPR      = 1.374D0
      GAMMA_KST   = 0.06672d0
      GAMMA_KSTPR = 0.240d0
      GAMMA_RCHT  =-0.039d0
        ELSE IF(KAK.EQ.71)  THEN     ! non default values to be used
                                     ! for KPI MODE WITH FSR INTERACTION
c   parameters for Kpi mode, reference [5]
      MKST = 0.8953d0
      GAMMA_KST = 0.0475d0
      MKSTPR = 1.307d0
      GAMMA_KSTPR = 0.206d0
      GAMMA_RCHT = -0.043d0
        ELSE
C   /MASS_SCAL/; stable particles - final scalars
      Mksp  = 0.89166d0   !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      Mks0  = 0.89610d0   !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MKST = (Mksp +Mks0)/2.
      MKSTPR = 1.374d0    !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GAMMA_KST   = 0.06672
      GAMMA_KSTPR = 0.240
      GAMMA_RCHT  = -0.043 !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
        ENDIF
C  /RCHT_3PI/; model parameters; their value are from fit,
c             reference [1], set 1
C              CHANGE OF THEIR VALUES REQUIRES
C              RERUN /tabler/a1/da1wid_tot_rho1_gauss.f
      IF(KAK.EQ.5) THEN
      FPI_RPT  = 0.091337d0
      FV_RPT   = 0.168652d0
      FA_RPT   = 0.131425d0
      BETA_RHO = -0.318551d0
      ELSE
      FPI_RPT  = 0.0924
      FV_RPT   = 0.18
      FA_RPT   = 0.149
      BETA_RHO =  -0.25
      ENDIF
      FK_RPT = FPI_RPT*1.198d0
      GV_RPT   = FPI_RPT*FPI_RPT/FV_RPT
c$$$c   It has to be used for a new parametrization of rho1 for 3pions,
C$$$c   that is not checked yet
c$$$c        IF(KAK.EQ.5)  THEN  ! high energy behaviour imposes these relations
c$$$c      GV_RPT   = 0.066
c$$$c      FV1_RPT = 0.18D0
c$$$c      GV1_RPT = (FPI_RPT*FPI_RPT- FV_RPT*GV_RPT)/FV1_RPT
c$$$c        ELSE
c$$$c      GV_RPT   = FPI_RPT*FPI_RPT/FV_RPT
c$$$c      ENDIF
c   /SCAL_3PI/; parameters of sigma meson for 3 pion mode
C* Parameteres for the sigma contribution, using BW for sigma
      IF(KAK.EQ.5)  THEN
          IF (J3PI.EQ.1) THEN
      alpsig = -8.795938d0
      betasig = 9.763701d0
      gamsig =  1.264263d0
      delsig =  0.656762d0
      rsigma =  1.866913d0
           ELSE IF (J3PI.EQ.2) THEN
      alpsig = 1.139486d0*0.63d0
      betasig = 1.139486d0*0.63d0
      gamsig =  0.889769d0*0.63d0
      delsig =  0.889769d0*0.63d0
      rsigma =  0.000013d0
           ENDIF
       ENDIF
C  /MASS_RES/
      IF(KAK.EQ.5)  THEN
      MMA1  = 1.091865d0      !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
          IF (J3PI.EQ.1) THEN
      MSIG  = 0.487512d0
      GSIG  = 0.70d0
          ELSE IF(J3PI.EQ.2) THEN
      MSIG  = 0.55d0
      GSIG  = 0.7d0
          ENDIF
      ELSE
      MMA1 = 1.12
      MSIG = 0.475
      GSIG = 0.550
      ENDIF
        call rchl_REparam(0,IWARM,KAK)
         IF (IWARM.EQ.1) RETURN  ! parameters below do not need
         IWARM=1                 ! re-initialization
C  /MASS_RES/
      GRO   = 0.149d0     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MF2   = 1.275d0
      GF2   = 0.185d0
      MF0   = 1.186d0
      GF0   = 0.350d0
      MSG   = 0.860d0
      GSG   = 0.880d0
      MPHI  = 1.019d0     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GPHI  = 0.0042d0    !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MOM   = 0.781940d0  !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GOM   = 0.00843d0   !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
C    /RES_MIXING_RCHT/; a parameter defines w-phi angle mixing
      THETA = 35.*PI/180.
C    /FF0SCKPI/ a parameter normalized FFSC_KPI
      F00 = 0.972
C   /MASS_SCAL/; stable particles - final scalars
C              CHANGE OF THEIR VALUES (useful for some tests) REQUIRES,
C              RERUN /tabler/a1/da1wid_tot_rho1_gauss.f
      MPIZ    = 0.1349766d0           !PKORB(1,7)    ! NEUTRAL PION MASS
      MPIC    = 0.13957018d0          !PKORB(1,8)    ! CHARGED PION MASS
      MMPI_AV = (MPIZ+2.*MPIC)/3.d0
      MKZ     = 0.497648d0            !PKORB(1,12)    ! NEUTRAL KAON MASS
      MKC     = 0.493677d0            !PKORB(1,11)    ! CHARGED KAON MASS
      MMK     = (MKC+MKZ)/2.d0
      MTAU    = 1.777
      MNUTA   = 0.001
      META    = 0.547d0
c   /PAR_KKPI/; parameters to describe KKpi modes, reference [2]
C              CHANGE OF THEIR VALUES REQUIRES
C              RERUN /tabler/a1/da1wid_tot_rho1_gauss.f
      G2    =  mro/(192.*pi*pi*sqrt(2.)*FV_RPT)*3.
      G13   = -2.*g2
      G4    = -0.72
      G5    = -0.6-2.*g4
      C125  =  0.
      C1256 = -3/96./pi**2*FV_RPT*MRO/SQRT(2.)/FPI_RPT**2
      C1235 =  0.
      C4    = -0.07
      D123  =  0.05
      D3    = -MRO**2/(64.*PI*PI*FPI_RPT**2)
c   /PAR_KPI/; parameters to describe Kpi mode, reference [4]
      Ht0 = -1.2400398216503017D-2
C     Ht0 = !!!!! TO ADD A FORMULAE FOR Ht0 (Jamin's email) !!!!!
      lap_KPI = 24.66e-3
      lapp_KPI = 11.99e-4
      c1_KPI = lap_KPI/mpic**2
      c2_KPI = (lapp_KPI - lap_kpi**2)/2.d0/mpic**4
c /KPISC_EM/; parameters for Kpi scalar FF from
c  http://arxiv.org/pdf/1103.4855.pdf
      lnC = 0.20193d0
      lambda0 = 0.013139d0
c   /SCAL_3PI/; parameters of sigma meson for 3 pion mode
      a00_3piscal = 0.220
      b00_3piscal = 0.268/mmpi_av**2
      c00_3piscal = -0.0139/mmpi_av**4
      d00_3piscal = -0.00139/mmpi_av**6
      x00_3piscal = 36.77*mmpi_av**2
      a02_3piscal = -0.0444
      b02_3piscal = -0.0857/mmpi_av**2
      c02_3piscal = -0.00221/mmpi_av**4
      d02_3piscal = -0.000129/mmpi_av**6
      x02_3piscal = -21.62*mmpi_av**2
      MMF0 = 0.441
c /SCAL_3PI/; parameters for the scalar part 3 pion modes
c   Pablo private
      ALPHA0_3PI = 1.
      ALPHA1_3PI = 1.
      GAMMA0_3PI = 1.
      GAMMA1_3PI = 1.
C FFVEC: dipswitch for Final State interaction in two scalar modes
C     with FSI (default FFVEC =1) and
C     without FSI (FFVEC =0)
      FFVEC = 1
C  FFKPIVEC : parameter to choose the parametrization for
C             vector Kpi form factor with FSI effects
C    FFKPIVEC = 0 parametrization Eqs.(17),(18) of [4]
C    FFKPIVEC = 1 parametrization Eq.(5) of [5]
C    FFKPIVEC = 2 parmetrization [4], total result
      FFKPIVEC = 2
C  FFKPISCAL : parameter to choose the parametrization for
C             scalar Kpi form factor with FSI effects
C  FFKPISCAL = 0 no scalar contribution
C  FFKPISCAL = 1 parametrization of Mathias Jamin,adopted his private code
C  FFKPISCAL = 2 parametrization of Emilie Passerman,
C               adopted her private code []
      FFKPISCAL = 1
C FFKKVEC: dipswitch for K0K- mode
C     with rho' and rho'' (FFKKVEC =1) and
C     without rho' and rho''  (default FFKKVEC =0)
      FFKKVEC = 0
C FF3PISCAL: dipswitch for the scalar contribution for 3 pion modes
C     with the scalar contribution ( default FF3PISCAL = 2)
c     FF3PISCAL = 2 BW parametrization for sigma meson
c     FF3PISCAL = 1 simplified RCHT results
C     FF3PISCAL =0  no sigma contribution
      FF3PISCAL = 2
C  Implemetation of another parametrization rho1, not checked yet by tests
C FF3PIRHOPR: dipswitch for the parametrization for rho' contribution
C     For 3 pion modes
C     general parametrization  ( default FF3RHOPR =1) and
C     simplified  (FF3PIRHOPR =0)
      FF3PIRHOPR = 0
C FF2PIRHO: dipswitch for the two pion form factor (default FF2PIRHO = 1)
C FF2PIRHO =1   RChL parametrization
C FF2PIRHO = 2  Belle parametrization,
C               all parameters par (1...11) of fit are free
C FF2PIRHO = 3  Belle parametrization,
C               parameters of fit are free
C               except for fixed par(1)=F_pi(0)=1
      FF2PIRHO =2
C FCOUL: dipswitch for the Coulomb interaction
C FCOUL = 1  with
C FCOUL = 0  without
      FCOUL = 0
      call rchl_REparam(1,IWARM,KAK)
      return
      end
      subroutine rchl_REparam(IMODE,IWARM,KAK)
C==============================================================================
C  common blocks and constants
C==============================================================================
      COMMON/RCHT_3PI/        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      DOUBLE PRECISION        FPI_RPT,FV_RPT,GV_RPT,FA_RPT,BETA_RHO,FK_RPT
     &                       ,FV1_RPT,GV1_RPT
      COMMON/MASS_RES/        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      DOUBLE PRECISION        MRO,GRO,MRHO1,GRHO1,MMA1,MKSP,MKS0
     &                       ,MF2,GF2,MF0,GF0,MSG,GSG,MOM,GOM,MPHI,GPHI
     &                       ,MKST,MSIG,GSIG
      COMMON/MASS_STABLE/     MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        MPIZ,MPIC,MMPI_AV,MKZ,MKC,MMK,MTAU,MNUTA
     &                       ,META
      DOUBLE PRECISION        PI
      DATA                    PI /3.141592653589793238462643D0/
      COMPLEX                 i
      DATA                    i/(0.,1.)/
      COMMON/PAR_KKPI/        G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      REAL                    G2,G13,G4,G5,C125,C1256,C1235,C4,D123,D3
      COMMON/RES_MIXING_RCHT/ THETA
      REAL                    THETA
      COMMON/PAR_RHOPRIME/    COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      DOUBLE PRECISION        COEF_GA,COEF_DE,phi_1,phi_2,grho2,mrho2
      COMMON/PAR_KPI/         MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      DOUBLE PRECISION        MKSTPR,GAMMA_KST,GAMMA_KSTPR,GAMMA_RCHT,Ht0
     &                       ,lap_KPI,lapp_KPI,c1_KPI,c2_KPI
      COMMON/VERSION_VECFF/   FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      INTEGER                 FFVEC,FFKPIVEC,FFKKVEC,FFKPISCAL,FF3PISCAL
     &                       ,FF3PIRHOPR,FFKKPIRHOPR,FFKKPIKPR,FF2PIRHO
     &                       ,FCOUL
      COMMON/FF0SCKPI/        F00
      REAL                    F00
      COMMON/KPISC_EM/        lnC,lambda0
      DOUBLE PRECISION        lnC,lambda0
      COMMON/SCAL_3PI/        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      DOUBLE PRECISION        a00_3piscal,b00_3piscal,c00_3piscal
     &                       ,d00_3piscal,x00_3piscal,a02_3piscal
     &                       ,b02_3piscal,c02_3piscal,d02_3piscal
     &                       ,x02_3piscal,ALPHA0_3PI,ALPHA1_3PI
     &                       ,GAMMA0_3PI,GAMMA1_3PI,MMF0
     &                       ,alpsig,betasig,gamsig,delsig,rsigma
      common / PARAMS / P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,IUSE
      INTEGER                                                                  IUSE
      DOUBLE PRECISION  P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16
      DATA IUSE /0/
!
      BIND(C,NAME='tauola_ff0sckpi') /FF0SCKPI/
      BIND(C,NAME='tauola_kpisc_em') /KPISC_EM/
      BIND(C,NAME='tauola_mass_res') /MASS_RES/
      BIND(C,NAME='tauola_mass_stable') /MASS_STABLE/
      BIND(C,NAME='tauola_params') /PARAMS/
      BIND(C,NAME='tauola_par_kkpi') /PAR_KKPI/
      BIND(C,NAME='tauola_par_kpi') /PAR_KPI/
      BIND(C,NAME='tauola_par_rhoprime') /PAR_RHOPRIME/
      BIND(C,NAME='tauola_rcht_3pi') /RCHT_3PI/
      BIND(C,NAME='tauola_res_mixing_rcht') /RES_MIXING_RCHT/
      BIND(C,NAME='tauola_scal_3pi') /SCAL_3PI/
      BIND(C,NAME='tauola_version_vecff') /VERSION_VECFF/
!
      IF(IUSE.EQ.0) RETURN
      IF (IMODE.EQ.-1) THEN
        IWARM=IWARM
      ELSE
C FF3PISCAL: dipswitch for the scalar contribution for 3 pion modes
C     with the scalar contribution ( default FF3PISCAL = 2)
c     FF3PISCAL = 2 BW parametrization for sigma meson
c     FF3PISCAL = 1 simplified RCHT results
C     FF3PISCAL =0  no sigma contribution
c      FF3PISCAL = 2
C     CANDIDATES FOR PARAMETERS TO FIT with default values
C* Parameteres for the sigma contribution, using BW for sigma
      alpsig  = P1
      betasig = P2
      gamsig  = P3
      delsig  = P4
      rsigma  = P5
      MRO   = P6     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MRHO1 = P7     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      GRHO1 = P8       !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
C  /MASS_RES/
      GRO   = P9     !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MMA1  = P10    !CHANGE REQUIRES RERUN da1wid_tot_rho1_gauss.f
      MSIG  = P11
      GSIG  = P12
C  /RCHT_3PI/; model parameters; their value are from fit,
c             reference [1], set 1
C              CHANGE OF THEIR VALUES REQUIRES
C              RERUN /tabler/a1/da1wid_tot_rho1_gauss.f
      FPI_RPT  = P13
      FV_RPT   = P14
      FA_RPT   = P15
      BETA_RHO = P16
      FK_RPT = FPI_RPT*1.198d0
      GV_RPT   = FPI_RPT*FPI_RPT/FV_RPT
      ENDIF
      return
      end
      SUBROUTINE IFGFACT(IKEY,ITYPE,IINIT)
C stores and manages options for   wid_a1_fit
C ITYPE=0,1,2 means zero value, from full table, from GFACT plus table for KKpi
C IINIT is for GFACT, whether it  has to reinitialize itself or just be used
C IKEY 0,1,2 means  initialization  use by GFUNCT (reini) other use
C IKEY -1 mieans initialization with recalculation og GFUNCT constats
      DATA ITYPE0,IINIT0 /1,0/
      SAVE ITYPE0,IINIT0
      IF    (IKEY.EQ.-1) THEN
        ITYPE0=ITYPE
        IINIT0=IINIT
      ELSEIF(IKEY.EQ.0) THEN
        ITYPE0=ITYPE
      ELSEIF(IKEY.EQ.1) THEN
        IINIT=IINIT0
        IINIT0=1
      ELSEIF(IKEY.EQ.2) THEN
        ITYPE= ITYPE0
      ENDIF
      END
      real function wid_a1_fit(qq)
      implicit none
      real                     qq
C..............................................................
C.    Output:    a1 width  as function of qq (formula (29) of ref [1])
C.               but to speed up execution linear interpolation
C.               from numerical table is used. From read in table
C.               linear interpolation (extrapolation) is used.
C.    Input:     qq  [GeV**2]
C.               table wida1_qq_tot_table.txt prepared in
C.               ./a1-tabler directory
C.    References:[1] arXiv:0911.4436 (hep-ph) D.Gomez Dumm et al.(tau ->3pi nu)
C.    Called    : files f3pi_rcht.f  fkkpi.f fkk0pi0.f
C                       ./demo-standalone/tautestroman.f
C..............................................................
      integer Nq
      PARAMETER (Nq=1001)
      integer ik
      common /a1_width/ qmax,qmin,qk_min,    w_qmin,    del_qq
      real*8            qmax,qmin,qk_min(Nq),w_qmin(Nq),del_qq
      integer kq
      character path*120
      real*8 aq,bq,qk_max,w_qmax
      real*8 xqq
      INTEGER IFIRST,ITYPE
      DATA    IFIRST,ITYPE /0,1/
      INTEGER  IMODE,IINIT
!      DATA IMODE /0/
C. READING TABLE AND INITIALIZATION
C. ================================
!
      BIND(C,NAME='tauola_a1_width') /a1_width/
!
      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
C Initialization from automatically created routine
        call InitA1TAB(qk_min    ,w_qmin    )
        qmax = qk_min(Nq)
        qmin = qk_min(1)
        del_qq = (qmax-qmin)/float(Nq-1)
      ENDIF
      CALL IFGFACT(2,ITYPE,IINIT)
C ALTERNATIVE CALCULATIONS FOR SPECIAL PURPOSES
C. ============================================
      IF (ITYPE.EQ.0) THEN
       wid_a1_fit=0
       RETURN
      ELSEIF (ITYPE.EQ.2) THEN
       xqq=qq
       wid_a1_fit=wid_a1_fitKKpi(qq)+GFACT(xqq)
!       write(*,*) 'uwaga', wid_a1_fit, wid_a1_fitKKpi(qq)
       RETURN
      ENDIF
C. INTEPOLATION, for extrapolation values at ends are used
C. =======================================================
      if(qq.gt.qmin.and.qq.le.qmax) then
        kq = (qq-qmin)/del_qq
        kq = kq +1
        qk_max = qk_min(kq+1)
        w_qmax = w_qmin(kq+1)
        aq = (w_qmax-w_qmin(kq))/(qk_max-qk_min(kq))
        bq = (w_qmax*qk_min(kq) -w_qmin(kq)*qk_max)
     $      /(qk_min(kq)-qk_max)
        wid_a1_fit = aq*qq+bq
      elseif(qq.ge.qmax) then               ! above maximun
        wid_a1_fit = w_qmin(Nq)
      else                                  ! below minimum
        wid_a1_fit = w_qmin(1)
      endif
      return
      end
      real function wid_a1_fitKKpi(qq)
      implicit none
      real                     qq
C..............................................................
C.    Output:    KKpi part ofa1 width  as function of qq
C.               but to speed up execution linear interpolation
C.               from numerical table is used. From read in table
C.               linear interpolation (extrapolation) is used.
C.    Input:     qq  [GeV**2]
C..............................................................
      integer Nq
      PARAMETER (Nq=1001)
      integer ik
      real*8 qmax,qmin,qk_min(Nq),w_qmin(Nq),del_qq
      integer kq
      character path*120
      real*8 aq,bq,qk_max,w_qmax
      INTEGER IFIRST
      DATA IFIRST/0/
      INTEGER  IMODE
!      DATA IMODE /0/
      COMMON   /IMODE/  IMODE      ! 0 is for calculation of G , 1 is for spectra, 2 is for use of G
!
      BIND(C,NAME='tauola_imode') /IMODE/
!
      IF (IMODE.EQ.0) THEN
       wid_a1_fitKKpi=0
      ENDIF
C. READING TABLE AND INITIALIZATION
C. ================================
      IF (IFIRST.EQ.0) THEN
        IFIRST = 1
C Initialization from automatically created routine
        call InitA1TABKKpi(qk_min    ,w_qmin    )
        qmax = qk_min(Nq)
        qmin = qk_min(1)
        del_qq = (qmax-qmin)/float(Nq-1)
      ENDIF
C. INTEPOLATION, for extrapolation values at ends are used
C. =======================================================
      if(qq.gt.qmin.and.qq.le.qmax) then
        kq = (qq-qmin)/del_qq
        kq = kq +1
        qk_max = qk_min(kq+1)
        w_qmax = w_qmin(kq+1)
        aq = (w_qmax-w_qmin(kq))/(qk_max-qk_min(kq))
        bq = (w_qmax*qk_min(kq) -w_qmin(kq)*qk_max)
     $      /(qk_min(kq)-qk_max)
        wid_a1_fitKKpi = aq*qq+bq
      elseif(qq.ge.qmax) then               ! above maximun
        wid_a1_fitKKpi = w_qmin(Nq)
      else                                  ! below minimum
        wid_a1_fitKKpi = w_qmin(1)
      endif
      return
      end
      SUBROUTINE CHOICE(MNUM,RR,ICHAN,PROB1,PROB2,PROB3,
     $            AMRX,GAMRX,AMRA,GAMRA,AMRB,GAMRB)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      AMROP=1.1
      GAMROP=0.36
      AMOM=.782
      GAMOM=0.0084
C     XXXXA CORRESPOND TO S2 CHANNEL !
      IF(MNUM.EQ.0) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =AMA1
       GAMRX=GAMA1
       AMRA =AMRO
       GAMRA=GAMRO
       AMRB =AMRO
       GAMRB=GAMRO
      ELSEIF(MNUM.EQ.1) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =1.57
       GAMRX=0.9
       AMRB =AMKST
       GAMRB=GAMKST
       AMRA =AMRO
       GAMRA=GAMRO
      ELSEIF(MNUM.EQ.2) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =1.57
       GAMRX=0.9
       AMRB =AMKST
       GAMRB=GAMKST
       AMRA =AMRO
       GAMRA=GAMRO
      ELSEIF(MNUM.EQ.3) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =1.27
       GAMRX=0.3
       AMRA =AMKST
       GAMRA=GAMKST
       AMRB =AMKST
       GAMRB=GAMKST
      ELSEIF(MNUM.EQ.4) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =1.27
       GAMRX=0.3
       AMRA =AMKST
       GAMRA=GAMKST
       AMRB =AMKST
       GAMRB=GAMKST
      ELSEIF(MNUM.EQ.5) THEN
       PROB1=0.5
       PROB2=0.5
       AMRX =1.27
       GAMRX=0.3
       AMRA =AMKST
       GAMRA=GAMKST
       AMRB =AMRO
       GAMRB=GAMRO
      ELSEIF(MNUM.EQ.6) THEN
       PROB1=0.4
       PROB2=0.4
       AMRX =1.27
       GAMRX=0.3
       AMRA =AMRO
       GAMRA=GAMRO
       AMRB =AMKST
       GAMRB=GAMKST
      ELSEIF(MNUM.EQ.7) THEN
       PROB1=0.0
       PROB2=1.0
       AMRX =1.27
       GAMRX=0.9
       AMRA =AMRO
       GAMRA=GAMRO
       AMRB =AMRO
       GAMRB=GAMRO
      ELSEIF(MNUM.EQ.8) THEN
       PROB1=0.0
       PROB2=1.0
       AMRX =AMROP
       GAMRX=GAMROP
       AMRB =AMOM
       GAMRB=GAMOM
       AMRA =AMRO
       GAMRA=GAMRO
      ELSEIF(MNUM.EQ.101) THEN
       PROB1=.35
       PROB2=.35
       AMRX =1.2
       GAMRX=.46
       AMRB =AMOM
       GAMRB=GAMOM
       AMRA =AMOM
       GAMRA=GAMOM
      ELSEIF(MNUM.EQ.102) THEN
       PROB1=0.0
       PROB2=0.0
       AMRX =1.4
       GAMRX=.6
       AMRB =AMOM
       GAMRB=GAMOM
       AMRA =AMOM
       GAMRA=GAMOM
      ELSE
       PROB1=0.0
       PROB2=0.0
       AMRX =AMA1
       GAMRX=GAMA1
       AMRA =AMRO
       GAMRA=GAMRO
       AMRB =AMRO
       GAMRB=GAMRO
      ENDIF
C
      IF    (RR.LE.PROB1) THEN
       ICHAN=1
      ELSEIF(RR.LE.(PROB1+PROB2)) THEN
       ICHAN=2
        AX   =AMRA
        GX   =GAMRA
        AMRA =AMRB
        GAMRA=GAMRB
        AMRB =AX
        GAMRB=GX
        PX   =PROB1
        PROB1=PROB2
        PROB2=PX
      ELSE
       ICHAN=3
      ENDIF
C
      PROB3=1.0-PROB1-PROB2
      END
      SUBROUTINE INITDK
* ----------------------------------------------------------------------
*     INITIALISATION OF TAU DECAY PARAMETERS  and routines
*
*     called by : KORALZ
* ----------------------------------------------------------------------
      COMMON / DECPAR / GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      REAL*4            GFERMI,GV,GA,CCABIB,SCABIB,GAMEL
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
      COMMON / TAUBRA / GAMPRT(30),JLIST(30),NCHAN
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      PARAMETER (NMODE=15,NM1=0,NM2=1,NM3=8,NM4=2,NM5=1,NM6=3)
      COMMON / TAUDCD /IDFFIN(9,NMODE),MULPIK(NMODE)
     &                ,NAMES
      CHARACTER NAMES(NMODE)*31
      CHARACTER OLDNAMES(7)*31
      CHARACTER*80 bxINIT
      PARAMETER (
     $  bxINIT ='(1x,1h*,g17.8,            16x, a31,a4,a4, 1x,1h*)'
     $ )
      REAL*4 PI
*
*
* LIST OF BRANCHING RATIOS
CAM normalised to e nu nutau channel
CAM                  enu   munu   pinu  rhonu   A1nu   Knu    K*nu   pi
CAM   DATA JLIST  /    1,     2,     3,     4,     5,     6,     7,
*AM   DATA GAMPRT /1.000,0.9730,0.6054,1.2432,0.8432,0.0432,O.O811,0.616
*AM
*AM  multipion decays
*
*    conventions of particles names
*                 K-,P-,K+,  K0,P-,KB,  K-,P0,K0
*                  3, 1,-3  , 4, 1,-4  , 3, 2, 4  ,
*                 P0,P0,K-,  K-,P-,P+,  P-,KB,P0
*                  2, 2, 3  , 3, 1,-1  , 1,-4, 2  ,
*                 ET,P-,P0   P-,P0,GM
*                  9, 1, 2  , 1, 2, 8
*
C
!
      BIND(C,NAME='tauola_decpar') /DECPAR/
      BIND(C,NAME='tauola_parmas') /PARMAS/
      BIND(C,NAME='tauola_taubra') /TAUBRA/
      BIND(C,NAME='tauola_taudcd') /TAUDCD/
      BIND(C,NAME='tauola_taukle') /TAUKLE/
!
      DIMENSION NOPIK(6,NMODE),NPIK(NMODE)
*AM   outgoing multiplicity and flavors of multi-pion /multi-K modes
      DATA   NPIK  /                4,                    4,
     1                              5,                    5,
     2                              6,                    6,
     3                              3,                    3,
     4                              3,                    3,
     5                              3,                    3,
     6                              3,                    3,
     7                              2                         /
      DATA  NOPIK / -1,-1, 1, 2, 0, 0,     2, 2, 2,-1, 0, 0,
     1              -1,-1, 1, 2, 2, 0,    -1,-1,-1, 1, 1, 0,
     2              -1,-1,-1, 1, 1, 2,    -1,-1, 1, 2, 2, 2,
     3              -3,-1, 3, 0, 0, 0,    -4,-1, 4, 0, 0, 0,
     4              -3, 2,-4, 0, 0, 0,     2, 2,-3, 0, 0, 0,
     5              -3,-1, 1, 0, 0, 0,    -1, 4, 2, 0, 0, 0,
     6               9,-1, 2, 0, 0, 0,    -1, 2, 8, 0, 0, 0,
C AJWMOD fix sign bug, 2/22/99
     7              -3,-4, 0, 0, 0, 0                         /
* LIST OF BRANCHING RATIOS
      NCHAN = NMODE + 7
      DO 1 I = 1,30
      IF (I.LE.NCHAN) THEN
        JLIST(I) = I
        IF(I.EQ. 1) GAMPRT(I) =0.1800
        IF(I.EQ. 2) GAMPRT(I) =0.1751
        IF(I.EQ. 3) GAMPRT(I) =0.1110
        IF(I.EQ. 4) GAMPRT(I) =0.2515
        IF(I.EQ. 5) GAMPRT(I) =0.1790
        IF(I.EQ. 6) GAMPRT(I) =0.0071
        IF(I.EQ. 7) GAMPRT(I) =0.0134
        IF(I.EQ. 8) GAMPRT(I) =0.0450
        IF(I.EQ. 9) GAMPRT(I) =0.0100
        IF(I.EQ.10) GAMPRT(I) =0.0009
        IF(I.EQ.11) GAMPRT(I) =0.0004
        IF(I.EQ.12) GAMPRT(I) =0.0003
        IF(I.EQ.13) GAMPRT(I) =0.0005
        IF(I.EQ.14) GAMPRT(I) =0.0015
        IF(I.EQ.15) GAMPRT(I) =0.0015
        IF(I.EQ.16) GAMPRT(I) =0.0015
        IF(I.EQ.17) GAMPRT(I) =0.0005
        IF(I.EQ.18) GAMPRT(I) =0.0050
        IF(I.EQ.19) GAMPRT(I) =0.0055
        IF(I.EQ.20) GAMPRT(I) =0.0017
        IF(I.EQ.21) GAMPRT(I) =0.0013
        IF(I.EQ.22) GAMPRT(I) =0.0010
        IF(I.EQ. 1) OLDNAMES(I)='  TAU-  -->   E-               '
        IF(I.EQ. 2) OLDNAMES(I)='  TAU-  -->  MU-               '
        IF(I.EQ. 3) OLDNAMES(I)='  TAU-  -->  PI-               '
        IF(I.EQ. 4) OLDNAMES(I)='  TAU-  -->  PI-, PI0          '
        IF(I.EQ. 5) OLDNAMES(I)='  TAU-  -->  A1- (two subch)   '
        IF(I.EQ. 6) OLDNAMES(I)='  TAU-  -->   K-               '
        IF(I.EQ. 7) OLDNAMES(I)='  TAU-  -->  K*- (two subch)   '
        IF(I.EQ. 8) NAMES(I-7)='  TAU-  --> 2PI-,  PI0,  PI+   '
        IF(I.EQ. 9) NAMES(I-7)='  TAU-  --> 3PI0,        PI-   '
        IF(I.EQ.10) NAMES(I-7)='  TAU-  --> 2PI-,  PI+, 2PI0   '
        IF(I.EQ.11) NAMES(I-7)='  TAU-  --> 3PI-, 2PI+,        '
        IF(I.EQ.12) NAMES(I-7)='  TAU-  --> 3PI-, 2PI+,  PI0   '
        IF(I.EQ.13) NAMES(I-7)='  TAU-  --> 2PI-,  PI+, 3PI0   '
        IF(I.EQ.14) NAMES(I-7)='  TAU-  -->  K-, PI-,  K+      '
        IF(I.EQ.15) NAMES(I-7)='  TAU-  -->  K0, PI-, K0B      '
        IF(I.EQ.16) NAMES(I-7)='  TAU-  -->  K-,  K0, PI0      '
        IF(I.EQ.17) NAMES(I-7)='  TAU-  --> PI0  PI0   K-      '
        IF(I.EQ.18) NAMES(I-7)='  TAU-  -->  K-  PI-  PI+      '
        IF(I.EQ.19) NAMES(I-7)='  TAU-  --> PI-  K0B  PI0      '
        IF(I.EQ.20) NAMES(I-7)='  TAU-  --> ETA  PI-  PI0      '
        IF(I.EQ.21) NAMES(I-7)='  TAU-  --> PI-  PI0  GAM      '
        IF(I.EQ.22) NAMES(I-7)='  TAU-  -->  K-  K0            '
      ELSE
        JLIST(I) = 0
        GAMPRT(I) = 0.
      ENDIF
   1  CONTINUE
      DO I=1,NMODE
        MULPIK(I)=NPIK(I)
        DO J=1,MULPIK(I)
         IDFFIN(J,I)=NOPIK(J,I)
        ENDDO
      ENDDO
*
*
* --- COEFFICIENTS TO FIX RATIO OF:
* --- A1 3CHARGED/ A1 1CHARGED 2 NEUTRALS MATRIX ELEMENTS (MASLESS LIM.)
* --- PROBABILITY OF K0 TO BE KS
* --- PROBABILITY OF K0B TO BE KS
* --- RATIO OF COEFFICIENTS FOR K*--> K0 PI-
* --- ALL COEFFICENTS SHOULD BE IN THE RANGE (0.0,1.0)
* --- THEY MEANING IS PROBABILITY OF THE FIRST CHOICE ONLY IF ONE
* --- NEGLECTS MASS-PHASE SPACE EFFECTS
      BRA1=0.5
      BRK0=0.5
      BRK0B=0.5
      BRKS=0.6667
*
      GFERMI = 1.16637E-5
      CCABIB = 0.975
      GV     = 1.0
      GA     =-1.0
* ZW 13.04.89 HERE WAS AN ERROR
      SCABIB = SQRT(1.-CCABIB**2)
      PI =4.*ATAN(1.)
      GAMEL  = GFERMI**2*AMTAU**5/(192*PI**3)
*
*      CALL DEXAY(-1,pol1)
*
      RETURN
      END
      FUNCTION DCDMAS(IDENT)
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      IF      (IDENT.EQ. 1) THEN
        APKMAS=AMPI
      ELSEIF  (IDENT.EQ.-1) THEN
        APKMAS=AMPI
      ELSEIF  (IDENT.EQ. 2) THEN
        APKMAS=AMPIZ
      ELSEIF  (IDENT.EQ.-2) THEN
        APKMAS=AMPIZ
      ELSEIF  (IDENT.EQ. 3) THEN
        APKMAS=AMK
      ELSEIF  (IDENT.EQ.-3) THEN
        APKMAS=AMK
      ELSEIF  (IDENT.EQ. 4) THEN
        APKMAS=AMKZ
      ELSEIF  (IDENT.EQ.-4) THEN
        APKMAS=AMKZ
      ELSEIF  (IDENT.EQ. 8) THEN
        APKMAS=0.0001
      ELSEIF  (IDENT.EQ.-8) THEN
        APKMAS=0.0001
      ELSEIF  (IDENT.EQ. 9) THEN
        APKMAS=0.5488
      ELSEIF  (IDENT.EQ.-9) THEN
        APKMAS=0.5488
      ELSE
        CALL TAUOLA_PRINT('STOP IN APKMAS, WRONG IDENT='//CHAR(0))
!       PRINT *, 'STOP IN APKMAS, WRONG IDENT=',IDENT
        CALL TAUOLA_STOP()
      ENDIF
      DCDMAS=APKMAS
      END
      FUNCTION LUNPIK(ID,ISGN)
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      REAL*4 XIO(1)
!
      BIND(C,NAME='tauola_taukle') /TAUKLE/
!
      IDENT=ID*ISGN
      IF      (IDENT.EQ. 1) THEN
        IPKDEF=-211
      ELSEIF  (IDENT.EQ.-1) THEN
        IPKDEF= 211
      ELSEIF  (IDENT.EQ. 2) THEN
        IPKDEF=111
      ELSEIF  (IDENT.EQ.-2) THEN
        IPKDEF=111
      ELSEIF  (IDENT.EQ. 3) THEN
        IPKDEF=-321
      ELSEIF  (IDENT.EQ.-3) THEN
        IPKDEF= 321
      ELSEIF  (IDENT.EQ. 4) THEN
*
* K0 --> K0_LONG (IS 130) / K0_SHORT (IS 310) = 1/1
        CALL RANMAR(XIO,1)
        IF (XIO(1).GT.BRK0) THEN
          IPKDEF= 130
        ELSE
          IPKDEF= 310
        ENDIF
      ELSEIF  (IDENT.EQ.-4) THEN
*
* K0B--> K0_LONG (IS 130) / K0_SHORT (IS 310) = 1/1
        CALL RANMAR(XIO,1)
        IF (XIO(1).GT.BRK0B) THEN
          IPKDEF= 130
        ELSE
          IPKDEF= 310
        ENDIF
      ELSEIF  (IDENT.EQ. 8) THEN
        IPKDEF= 22
      ELSEIF  (IDENT.EQ.-8) THEN
        IPKDEF= 22
      ELSEIF  (IDENT.EQ. 9) THEN
        IPKDEF= 221
      ELSEIF  (IDENT.EQ.-9) THEN
        IPKDEF= 221
      ELSE
        CALL TAUOLA_PRINT('STOP IN IPKDEF, WRONG IDENT='//CHAR(0))
!       PRINT *, 'STOP IN IPKDEF, WRONG IDENT=',IDENT
        CALL TAUOLA_STOP()
      ENDIF
      LUNPIK=IPKDEF
      END
      SUBROUTINE TAURDF(KTO)
C THIS ROUTINE CAN BE CALLED BEFORE ANY TAU+ OR TAU- EVENT IS GENERATED
C IT CAN BE USED TO GENERATE TAU+ AND TAU- SAMPLES OF DIFFERENT
C CONTENTS
      COMMON / TAUKLE / BRA1,BRK0,BRK0B,BRKS
      REAL*4            BRA1,BRK0,BRK0B,BRKS
      COMMON / TAUBRA / GAMPRT(30),JLIST(30),NCHAN
C Subroutine TAURDF is disabled
!
      BIND(C,NAME='tauola_taubra') /TAUBRA/
      BIND(C,NAME='tauola_taukle') /TAUKLE/
!
      RETURN
      IF (KTO.EQ.1) THEN
C     ==================
C AJWMOD: Set the BRs for (A1+ -> rho+ pi0) and (K*+ -> K0 pi+)
      BRA1 = PKORB(4,1)
      BRKS = PKORB(4,3)
      BRK0  = PKORB(4,5)
      BRK0B  = PKORB(4,6)
      ELSE
C     ====
C AJWMOD: Set the BRs for (A1+ -> rho+ pi0) and (K*+ -> K0 pi+)
      BRA1 = PKORB(4,2)
      BRKS = PKORB(4,4)
      BRK0  = PKORB(4,5)
      BRK0B  = PKORB(4,6)
      ENDIF
C     =====
      END
      SUBROUTINE INIPHY(XK00)
* ----------------------------------------------------------------------
*     INITIALISATION OF PARAMETERS
*     USED IN QED and/or GSW ROUTINES
* ----------------------------------------------------------------------
      COMMON / QEDPRM /ALFINV,ALFPI,XK0
      REAL*8           ALFINV,ALFPI,XK0
      REAL*8 PI8,XK00
*
!
      BIND(C,NAME='tauola_qedprm') /QEDPRM/
!
      PI8    = 4.D0*DATAN(1.D0)
      ALFINV = 137.03604D0
      ALFPI  = 1D0/(ALFINV*PI8)
      XK0=XK00
      END
      SUBROUTINE INIMAS
C ----------------------------------------------------------------------
C     INITIALISATION OF MASSES
C
C     called by : KORALZ
C ----------------------------------------------------------------------
      COMMON / PARMAS / AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
*
      REAL*4            AMTAU,AMNUTA,AMEL,AMNUE,AMMU,AMNUMU
     *                 ,AMPIZ,AMPI,AMRO,GAMRO,AMA1,GAMA1
     *                 ,AMK,AMKZ,AMKST,GAMKST
C
C IN-COMING / OUT-GOING  FERMION MASSES
!
      BIND(C,NAME='tauola_parmas') /PARMAS/
!
      AMTAU  = 1.7842
C --- let us update tau mass ...
      AMTAU  = 1.777
      AMNUTA = 0.010
      AMEL   = 0.0005111
      AMNUE  = 0.0
      AMMU   = 0.105659
      AMNUMU = 0.0
*
* MASSES USED IN TAU DECAYS
      AMPIZ  = 0.134964
      AMPI   = 0.139568
      AMRO   = 0.773
      GAMRO  = 0.145
*C    GAMRO  = 0.666
      AMA1   = 1.251
      GAMA1  = 0.599
      AMK    = 0.493667
      AMKZ   = 0.49772
      AMKST  = 0.8921
      GAMKST = 0.0513
C
C
C IN-COMING / OUT-GOING  FERMION MASSES
!!      AMNUTA = PKORB(1,2)
!!      AMNUE  = PKORB(1,4)
!!      AMNUMU = PKORB(1,6)
C
C MASSES USED IN TAU DECAYS  Cleo settings
!!      AMPIZ  = PKORB(1,7)
!!      AMPI   = PKORB(1,8)
!!      AMRO   = PKORB(1,9)
!!      GAMRO  = PKORB(2,9)
      AMA1   = 1.275   !! PKORB(1,10)
      GAMA1  = 0.615   !! PKORB(2,10)
!!      AMK    = PKORB(1,11)
!!      AMKZ   = PKORB(1,12)
!!      AMKST  = PKORB(1,13)
!!      GAMKST = PKORB(2,13)
C
      RETURN
      END
      SUBROUTINE TRALO4(KTO,P,Q,AMS)
      INTEGER KTO
      REAL*4 P(4),Q(4),AMS
      END SUBROUTINE TRALO4
      END SUBROUTINE TAUOLA_DECAY