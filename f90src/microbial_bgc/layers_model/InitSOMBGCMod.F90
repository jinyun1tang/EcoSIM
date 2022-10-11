module InitSOMBGCMOD
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use MicrobialDataType
  use SOMDataType
  use GridConsts
  use minimathmod, only : AZMAX1
  use SoilPhysDataType
  use FlagDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use EcosimConst
  use SurfLitterDataType
  use SoilPropertyDataType
  use GridDataType
  use MicBGCPars, only : micpar

  implicit none

  private


  real(r8), allocatable :: CORGCX(:)  !C concentations from OM complexes
  real(r8), allocatable :: CORGNX(:)  !N concentations from OM complexes
  real(r8), allocatable :: CORGPX(:)  !P concentations from OM complexes

  public :: InitSOMVars
  public :: InitSOMProfile
  public :: InitSOMConsts
  public :: InitSOMBGC
  public :: DestructSOMBGC
  contains
!------------------------------------------------------------------------------------------

  subroutine InitSOMBGC(nmicbguilds)
  implicit none
  integer, intent(in) :: nmicbguilds

  call MicPar%Init(nmicbguilds)

  JGnio  => micpar%JGnio
  JGnfo  => micpar%JGnfo
  JGniA  => micpar%JGniA
  JGnfA  => micpar%JGnfA
  NMICBSA= micpar%NMICBSA
  NMICBSO= micpar%NMICBSO

  JG=micpar%jguilds
  allocate(CORGCX(0:jcplx1))
  allocate(CORGNX(0:jcplx1))
  allocate(CORGPX(0:jcplx1))
  end subroutine InitSOMBGC

!------------------------------------------------------------------------------------------

  subroutine InitSOMConsts
  use MicBGCPars, only : micpar
  implicit none

  call MicPar%SetPars

  end subroutine InitSOMConsts
!------------------------------------------------------------------------------------------

  subroutine InitSOMVars(L,NY,NX,FCX)

  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: FCX

  real(r8), parameter :: DCKR=0.25_r8,DCKM=2.5E+04_r8

  integer :: K,M,KK,N,NGL,NN
  real(r8) :: OC,ON,OP,X
  real(r8) :: OMC1,OMN1,OMP1
  real(r8) :: FOSCI,FOSNI,FOSPI
  real(r8) :: RNT,RPT
  real(r8) :: FRNT,FRPT
  real(r8) :: TOSNI,TOSPI,TOSCI
  real(r8) :: OSCI(0:jcplx1),OSNI(0:jcplx1),OSPI(0:jcplx1)
  real(r8) :: OSCM(0:jcplx1)
  real(r8) :: TOSCK(0:jcplx1),TOSNK(0:jcplx1),TOSPK(0:jcplx1)
  real(r8) :: RC
  real(r8) :: CNOSCT(0:jcplx1), CPOSCT(0:jcplx1)
  real(r8) :: OSCX(0:jcplx1)
  real(r8) :: OSNX(0:jcplx1)
  real(r8) :: OSPX(0:jcplx1)

! begin_execution

  associate(                  &
    CNOMC   => micpar%CNOMC  ,&
    CPOMC   => micpar%CPOMC  ,&
    OHCK    => micpar%OHCK   ,&
    OMCK    => micpar%OMCK   ,&
    OQCK    => micpar%OQCK   ,&
    ORCK    => micpar%ORCK   ,&
    ORCI    => micpar%ORCI   ,&
    OMCI    => micpar%OMCI   ,&
    CNRH    => micpar%CNRH   ,&
    CPRH    => micpar%CPRH   ,&
    CNOFC   => micpar%CNOFC  ,&
    CPOFC   => micpar%CPOFC  ,&
    OMCF    => micpar%OMCF   ,&
    OMCA    => micpar%OMCA    &
  )
  DO 975 K=0,2
    CNOSCT(K)=0.0_r8
    CPOSCT(K)=0.0_r8
    IF(RSC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      RNT=0.0_r8
      RPT=0.0_r8
      DO 970 M=1,jsken
        RNT=RNT+RSC(K,L,NY,NX)*CFOSC(M,K,L,NY,NX)*CNOFC(M,K)
        RPT=RPT+RSC(K,L,NY,NX)*CFOSC(M,K,L,NY,NX)*CPOFC(M,K)
970   CONTINUE
      FRNT=RSN(K,L,NY,NX)/RNT
      FRPT=RSP(K,L,NY,NX)/RPT
      DO 960 M=1,jsken
        CNOSC(M,K,L,NY,NX)=CNOFC(M,K)*FRNT
        CPOSC(M,K,L,NY,NX)=CPOFC(M,K)*FRPT
        CNOSCT(K)=CNOSCT(K)+CFOSC(M,K,L,NY,NX)*CNOSC(M,K,L,NY,NX)
        CPOSCT(K)=CPOSCT(K)+CFOSC(M,K,L,NY,NX)*CPOSC(M,K,L,NY,NX)
960   CONTINUE
    ELSE
      DO 965 M=1,jsken
        CNOSC(M,K,L,NY,NX)=CNRH(K)
        CPOSC(M,K,L,NY,NX)=CPRH(K)
965   CONTINUE
      CNOSCT(K)=CNRH(K)
      CPOSCT(K)=CPRH(K)
    ENDIF
975 CONTINUE
  DO 990 K=3,jcplx1
    CNOSCT(K)=0.0_r8
    CPOSCT(K)=0.0_r8
    IF(CORGCX(K).GT.ZERO)THEN
      DO 985 M=1,jsken
        CNOSC(M,K,L,NY,NX)=CORGNX(K)/CORGCX(K)
        CPOSC(M,K,L,NY,NX)=CORGPX(K)/CORGCX(K)
        CNOSCT(K)=CNOSCT(K)+CFOSC(M,K,L,NY,NX)*CNOSC(M,K,L,NY,NX)
        CPOSCT(K)=CPOSCT(K)+CFOSC(M,K,L,NY,NX)*CPOSC(M,K,L,NY,NX)
985   CONTINUE
    ELSE
      DO 980 M=1,jsken
        CNOSC(M,K,L,NY,NX)=CNRH(K)
        CPOSC(M,K,L,NY,NX)=CPRH(K)
980   CONTINUE
      CNOSCT(K)=CNRH(K)
      CPOSCT(K)=CPRH(K)
    ENDIF
990 CONTINUE
!
!     MICROBIAL BIOMASS,RESIDUE, DOC, ADSORBED
!
!     OSCI,OSNI,OSPI=initial SOC,SON,SOP mass in each complex (g)
!     OMCK,ORCK,OQCK,OHCK=fractions of SOC in biomass,litter,DOC adsorbed C
!     OSCM=total biomass in each complex (g)
!     DCKR,DCKM=parameters to initialize microbial biomass from SOC
!
  TOSCI=0.0_r8
  TOSNI=0.0_r8
  TOSPI=0.0_r8
  DO 995 K=0,jcplx1
    IF(L.EQ.0)THEN
      KK=K
    ELSE
      KK=4
    ENDIF
    IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
      OSCI(K)=CORGCX(K)*BKVL(L,NY,NX)
      OSNI(K)=CORGNX(K)*BKVL(L,NY,NX)
      OSPI(K)=CORGPX(K)*BKVL(L,NY,NX)
    ELSE
      OSCI(K)=CORGCX(K)*VOLT(L,NY,NX)
      OSNI(K)=CORGNX(K)*VOLT(L,NY,NX)
      OSPI(K)=CORGPX(K)*VOLT(L,NY,NX)
    ENDIF
    TOSCK(K)=OMCK(K)+ORCK(K)+OQCK(K)+OHCK(K)
    TOSNK(K)=ORCK(K)*CNRH(K)+OQCK(K)*CNOSCT(KK)+OHCK(K)*CNOSCT(KK)
    TOSPK(K)=ORCK(K)*CPRH(K)+OQCK(K)*CPOSCT(KK)+OHCK(K)*CPOSCT(KK)
    do NGL=1,JG
      TOSNK(K)=TOSNK(K)+OMCI(1,K)*CNOMC(1,NGL,1,K)+OMCI(2,K)*CNOMC(2,NGL,1,K)
      TOSPK(K)=TOSPK(K)+OMCI(1,K)*CPOMC(1,NGL,1,K)+OMCI(2,K)*CPOMC(2,NGL,1,K)
    enddo
    TOSCI=TOSCI+OSCI(K)*TOSCK(K)
    TOSNI=TOSNI+OSCI(K)*TOSNK(K)
    TOSPI=TOSPI+OSCI(K)*TOSPK(K)
    OSCX(K)=0.0_r8
    OSNX(K)=0.0_r8
    OSPX(K)=0.0_r8
995 CONTINUE

  DO 8995 K=0,jcplx1
    IF(L.EQ.0)THEN
      OSCM(K)=DCKR*CORGCX(K)*BKVL(L,NY,NX)
      X=0.0_r8
      KK=K
      FOSCI=1.0
      FOSNI=1.0
      FOSPI=1.0
    ELSE
      IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
        IF(K.LE.2)THEN
          OSCM(K)=DCKR*CORGCX(K)*BKVL(L,NY,NX)
        ELSE
          OSCM(K)=FCX*CORGCX(K)*BKVL(L,NY,NX)*DCKM/(CORGCX(4)+DCKM)
        ENDIF
      ELSE
        IF(K.LE.2)THEN
          OSCM(K)=DCKR*CORGCX(K)*VOLT(L,NY,NX)
        ELSE
          OSCM(K)=FCX*CORGCX(K)*VOLT(L,NY,NX)*DCKM/(CORGCX(4)+DCKM)
        ENDIF
      ENDIF
      !     IF(L.EQ.NU(NY,NX))THEN
      !     WRITE(*,2424)'OSCM',NX,NY,L,K,OSCM(K),OSCI(K),CORGCX(K)
      !    2,BKVL(L,NY,NX),CORGCX(K)*BKVL(L,NY,NX),FCX
!2424  FORMAT(A8,4I4,12E12.4)
      !     ENDIF
      X=1.0
      KK=4
      IF(TOSCI.GT.ZEROS(NY,NX))THEN
        FOSCI=AMIN1(1.0,OSCI(KK)/TOSCI)
        FOSNI=AMIN1(1.0,OSCI(KK)*CNOSCT(KK)/TOSNI)
        FOSPI=AMIN1(1.0,OSCI(KK)*CPOSCT(KK)/TOSPI)
      ELSE
        FOSCI=0.0_r8
        FOSNI=0.0_r8
        FOSPI=0.0_r8
      ENDIF
    ENDIF
!
!     MICROBIAL C, N AND P
!
!     OMC,OMN,OMP=microbial C,N,P
!     OMCI=microbial biomass content in litter
!     OMCF,OMCA=hetero,autotrophic biomass composition in litter
!     CNOMC,CPOMC=maximum N:C and P:C ratios in microbial biomass
!     OSCX,OSNX,OSPX=remaining unallocated SOC,SON,SOP
!  The reason that initialization of complex 5 microbes is repated for each
! complex is because complex 5 is shared by the other complexes
    DO 7990 N=1,NFGs
      DO NGL=1,JG
        DO 7985 M=1,3
          OMCff(M,NGL,N,L,NY,NX)=0.0_r8
          OMNff(M,NGL,N,L,NY,NX)=0.0_r8
          OMPff(M,NGL,N,L,NY,NX)=0.0_r8
7985    CONTINUE
      enddo
7990  CONTINUE

    DO 8990 N=1,NFGs
      do NGL=1,JG
        DO 8991 M=1,3
          OMC1=AZMAX1(OSCM(K)*OMCI(M,K)*OMCF(N)*FOSCI)
          OMN1=AZMAX1(OMC1*CNOMC(M,NGL,N,K)*FOSNI)
          OMP1=AZMAX1(OMC1*CPOMC(M,NGL,N,K)*FOSPI)
          OMC(M,NGL,N,K,L,NY,NX)=OMC1
          OMN(M,NGL,N,K,L,NY,NX)=OMN1
          OMP(M,NGL,N,K,L,NY,NX)=OMP1
          OSCX(KK)=OSCX(KK)+OMC1
          OSNX(KK)=OSNX(KK)+OMN1
          OSPX(KK)=OSPX(KK)+OMP1
          DO 8992 NN=1,NFGs
            OMCff(M,NGL,NN,L,NY,NX)=OMCff(M,NGL,NN,L,NY,NX)+OMC1*OMCA(NN)
            OMNff(M,NGL,NN,L,NY,NX)=OMNff(M,NGL,NN,L,NY,NX)+OMN1*OMCA(NN)
            OMPff(M,NGL,NN,L,NY,NX)=OMPff(M,NGL,NN,L,NY,NX)+OMP1*OMCA(NN)
            OSCX(KK)=OSCX(KK)+OMC1*OMCA(NN)
            OSNX(KK)=OSNX(KK)+OMN1*OMCA(NN)
            OSPX(KK)=OSPX(KK)+OMP1*OMCA(NN)
8992      CONTINUE
8991    CONTINUE
      enddo
8990  CONTINUE
!
!     MICROBIAL RESIDUE C, N AND P
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     ORCI=allocation of microbial residue to kinetic components
!  X is an indicator of surface residual layer
    DO 8985 M=1,2
      ORC(M,K,L,NY,NX)=X*AZMAX1(OSCM(K)*ORCI(M,K)*FOSCI)
      ORN(M,K,L,NY,NX)=AZMAX1(ORC(M,K,L,NY,NX)*CNOMC(M,1,1,K)*FOSNI)
      ORP(M,K,L,NY,NX)=AZMAX1(ORC(M,K,L,NY,NX)*CPOMC(M,1,1,K)*FOSPI)
      OSCX(KK)=OSCX(KK)+ORC(M,K,L,NY,NX)
      OSNX(KK)=OSNX(KK)+ORN(M,K,L,NY,NX)
      OSPX(KK)=OSPX(KK)+ORP(M,K,L,NY,NX)
8985  CONTINUE
!
!     DOC, DON AND DOP
!
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores (g)
!     OQCH,OQNH,OQPH,OQAH=DOC,DON,DOP,acetate in macropores (g)
!
    OQC(K,L,NY,NX)=X*AZMAX1(OSCM(K)*OQCK(K)*FOSCI)
    OQN(K,L,NY,NX)=AZMAX1(OQC(K,L,NY,NX)*CNOSCT(KK)*FOSNI)
    OQP(K,L,NY,NX)=AZMAX1(OQC(K,L,NY,NX)*CPOSCT(KK)*FOSPI)
    OQA(K,L,NY,NX)=0.0_r8
    OQCH(K,L,NY,NX)=0.0_r8
    OQNH(K,L,NY,NX)=0.0_r8
    OQPH(K,L,NY,NX)=0.0_r8
    OQAH(K,L,NY,NX)=0.0_r8
    OSCX(KK)=OSCX(KK)+OQC(K,L,NY,NX)
    OSNX(KK)=OSNX(KK)+OQN(K,L,NY,NX)
    OSPX(KK)=OSPX(KK)+OQP(K,L,NY,NX)
!
!     ADSORBED C, N AND P
!
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!
    OHC(K,L,NY,NX)=X*AZMAX1(OSCM(K)*OHCK(K)*FOSCI)
    OHN(K,L,NY,NX)=AZMAX1(OHC(K,L,NY,NX)*CNOSCT(KK)*FOSNI)
    OHP(K,L,NY,NX)=AZMAX1(OHC(K,L,NY,NX)*CPOSCT(KK)*FOSPI)
    OHA(K,L,NY,NX)=0.0_r8
    OSCX(KK)=OSCX(KK)+OHC(K,L,NY,NX)+OHA(K,L,NY,NX)
    OSNX(KK)=OSNX(KK)+OHN(K,L,NY,NX)
    OSPX(KK)=OSPX(KK)+OHP(K,L,NY,NX)
!
!     HUMUS C, N AND P
!
!     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP

    DO 8980 M=1,jsken
      OSC(M,K,L,NY,NX)=AZMAX1(CFOSC(M,K,L,NY,NX)*(OSCI(K)-OSCX(K)))
      IF(CNOSCT(K).GT.ZERO)THEN
        OSN(M,K,L,NY,NX)=AZMAX1(CFOSC(M,K,L,NY,NX)*CNOSC(M,K,L,NY,NX) &
          /CNOSCT(K)*(OSNI(K)-OSNX(K)))
      ELSE
        OSN(M,K,L,NY,NX)=0.0_r8
      ENDIF
      IF(CPOSCT(K).GT.ZERO)THEN
        OSP(M,K,L,NY,NX)=AZMAX1(CFOSC(M,K,L,NY,NX)*CPOSC(M,K,L,NY,NX) &
          /CPOSCT(K)*(OSPI(K)-OSPX(K)))
      ELSE
        OSP(M,K,L,NY,NX)=0.0_r8
      ENDIF
      IF(K.EQ.0)THEN
        do NGL=1,JG
          OSA(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)*OMCI(1+(NGL-1)*3,K)
        ENDDO
      ELSE
        OSA(M,K,L,NY,NX)=OSC(M,K,L,NY,NX)
      ENDIF
8980  CONTINUE
8995  CONTINUE
!
!     ADD ALL LITTER,POC,HUMUS COMPONENTS TO GET TOTAL SOC
!
  OC=0.0_r8
  ON=0.0_r8
  OP=0.0_r8
  RC=0.0_r8
  IF(L.EQ.0)THEN
    DO 6975 K=0,jcplx1
      RC0(K,NY,NX)=0.0_r8
6975  CONTINUE
    RC0ff(NY,NX)=0._r8
  ENDIF

  DO 6990 K=0,jcplx1
    DO  N=1,NFGs
      do NGL=1,JG
        ROXYS(NGL,N,K,L,NY,NX)=0.0_r8
        RVMX4(NGL,N,K,L,NY,NX)=0.0_r8
        RVMX3(NGL,N,K,L,NY,NX)=0.0_r8
        RVMX2(NGL,N,K,L,NY,NX)=0.0_r8
        RVMX1(NGL,N,K,L,NY,NX)=0.0_r8
        RINHO(NGL,N,K,L,NY,NX)=0.0_r8
        RINHB(NGL,N,K,L,NY,NX)=0.0_r8
        RINOO(NGL,N,K,L,NY,NX)=0.0_r8
        RIPOO(NGL,N,K,L,NY,NX)=0.0_r8
        RINOB(NGL,N,K,L,NY,NX)=0.0_r8
        RIPBO(NGL,N,K,L,NY,NX)=0.0_r8
        RIPO1(NGL,N,K,L,NY,NX)=0.0_r8
        RIPB1(NGL,N,K,L,NY,NX)=0.0_r8
        IF(L.EQ.0)THEN
          RINHOR(NGL,N,K,NY,NX)=0.0_r8
          RINOOR(NGL,N,K,NY,NX)=0.0_r8
          RIPOOR(NGL,N,K,NY,NX)=0.0_r8
        ENDIF
        DO  M=1,3
          OC=OC+OMC(M,NGL,N,K,L,NY,NX)
          ON=ON+OMN(M,NGL,N,K,L,NY,NX)
          OP=OP+OMP(M,NGL,N,K,L,NY,NX)
          IF(K.LE.2)THEN
            RC=RC+OMC(M,NGL,N,K,L,NY,NX)
          ENDIF
          RC0(K,NY,NX)=RC0(K,NY,NX)+OMC(M,NGL,N,K,L,NY,NX)
        ENDDO
      ENDDO
    enddo
6990  CONTINUE

    DO  N=1,NFGs
      do NGL=1,JG
        ROXYSff(NGL,N,L,NY,NX)=0.0_r8
        RVMX4ff(NGL,N,L,NY,NX)=0.0_r8
        RVMX3ff(NGL,N,L,NY,NX)=0.0_r8
        RVMB3ff(NGL,N,L,NY,NX)=0.0_r8
        RVMX2ff(NGL,N,L,NY,NX)=0.0_r8
        RVMX1ff(NGL,N,L,NY,NX)=0.0_r8
        RINHOff(NGL,N,L,NY,NX)=0.0_r8
        RINHBff(NGL,N,L,NY,NX)=0.0_r8
        RINOOff(NGL,N,L,NY,NX)=0.0_r8
        RINOBff(NGL,N,L,NY,NX)=0.0_r8
        RIPOOff(NGL,N,L,NY,NX)=0.0_r8
        RIPBOff(NGL,N,L,NY,NX)=0.0_r8
        RIPO1ff(NGL,N,L,NY,NX)=0.0_r8
        RIPB1ff(NGL,N,L,NY,NX)=0.0_r8
        IF(L.EQ.0)THEN
          RINHORff(NGL,N,NY,NX)=0.0_r8
          RINOORff(NGL,N,NY,NX)=0.0_r8
          RIPOORff(NGL,N,NY,NX)=0.0_r8
        ENDIF
        DO  M=1,3
          OC=OC+OMCff(M,NGL,N,L,NY,NX)
          ON=ON+OMNff(M,NGL,N,L,NY,NX)
          OP=OP+OMPff(M,NGL,N,L,NY,NX)
          RC0ff(NY,NX)=RC0(K,NY,NX)+OMCff(M,NGL,N,L,NY,NX)
        ENDDO
      ENDDO
    enddo

  DO 6995 K=0,jcplx1
    DO 6985 M=1,2
      OC=OC+ORC(M,K,L,NY,NX)
      ON=ON+ORN(M,K,L,NY,NX)
      OP=OP+ORP(M,K,L,NY,NX)
      IF(K.LE.2)THEN
        RC=RC+ORC(M,K,L,NY,NX)
      ENDIF
      IF(L.EQ.0)THEN
        RC0(K,NY,NX)=RC0(K,NY,NX)+ORC(M,K,L,NY,NX)
      ENDIF
6985  CONTINUE
    OC=OC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
      +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
    ON=ON+OQN(K,L,NY,NX)+OQNH(K,L,NY,NX)+OHN(K,L,NY,NX)
    OP=OP+OQP(K,L,NY,NX)+OQPH(K,L,NY,NX)+OHP(K,L,NY,NX)
    OC=OC+OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)
    IF(K.LE.2)THEN
      RC=RC+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX)+OHC(K,L,NY,NX) &
        +OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
      RC=RC+OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)
    ENDIF
    IF(L.EQ.0)THEN
      RC0(K,NY,NX)=RC0(K,NY,NX)+OQC(K,L,NY,NX)+OQCH(K,L,NY,NX) &
        +OHC(K,L,NY,NX)+OQA(K,L,NY,NX)+OQAH(K,L,NY,NX)+OHA(K,L,NY,NX)
    ENDIF
    DO 6980 M=1,jsken
      OC=OC+OSC(M,K,L,NY,NX)
      ON=ON+OSN(M,K,L,NY,NX)
      OP=OP+OSP(M,K,L,NY,NX)
      IF(K.LE.2)THEN
        RC=RC+OSC(M,K,L,NY,NX)
      ENDIF
      IF(L.EQ.0)THEN
        RC0(K,NY,NX)=RC0(K,NY,NX)+OSC(M,K,L,NY,NX)
      ENDIF
6980  CONTINUE
6995  CONTINUE
  ORGC(L,NY,NX)=OC
  ORGCX(L,NY,NX)=ORGC(L,NY,NX)
  ORGR(L,NY,NX)=RC
  ORGN(L,NY,NX)=ON
  end associate
  end subroutine InitSOMVars

!------------------------------------------------------------------------------------------
  subroutine InitSurfResiduKinetiComponent(L,NY,NX)
  implicit none
  integer, intent(in) :: L, NY,NX
!     begin_execution
  IF(L.EQ.0)THEN
    !
    !     CFOSC=fraction of litter in protein(1),nonstructural(2)
    !     cellulose(3) and lignin(4)
    !
    !     PREVIOUS COARSE WOODY RESIDUE
    !
    CFOSC(1:4,0,L,NY,NX)=real((/0.00,0.045,0.660,0.295/),r8)

    !
    !     MAIZE
    !
    IF(IXTYP(1,NY,NX).EQ.1)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.080,0.245,0.613,0.062/),r8)
      !
      !     WHEAT
      !
    ELSEIF(IXTYP(1,NY,NX).EQ.2)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.125,0.171,0.560,0.144/),r8)
!
      !     SOYBEAN
      !
    ELSEIF(IXTYP(1,NY,NX).EQ.3)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.138,0.426,0.316,0.120/),r8)

      !     NEW STRAW
!
    ELSEIF(IXTYP(1,NY,NX).EQ.4)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.036,0.044,0.767,0.153/),r8)

      !     OLD STRAW
!
    ELSEIF(IXTYP(1,NY,NX).EQ.5)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.075,0.125,0.550,0.250/),r8)

!
      !     COMPOST
!
    ELSEIF(IXTYP(1,NY,NX).EQ.6)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.143,0.015,0.640,0.202/),r8)
!
      !     GREEN MANURE
!
    ELSEIF(IXTYP(1,NY,NX).EQ.7)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.202,0.013,0.560,0.225/),r8)

      !     NEW DECIDUOUS FOREST
!
    ELSEIF(IXTYP(1,NY,NX).EQ.8)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.070,0.41,0.36,0.16/),r8)

!
      !     NEW CONIFEROUS FOREST
!
    ELSEIF(IXTYP(1,NY,NX).EQ.9)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.07,0.25,0.38,0.30/),r8)

      !     OLD DECIDUOUS FOREST
!
    ELSEIF(IXTYP(1,NY,NX).EQ.10)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.02,0.06,0.34,0.58/),r8)
!
      !     OLD CONIFEROUS FOREST
!
    ELSEIF(IXTYP(1,NY,NX).EQ.11)THEN
      CFOSC(1:4,1,L,NY,NX)=real((/0.02,0.06,0.34,0.58/),r8)
      !     DEFAULT
!
    ELSE
      CFOSC(1:4,1,L,NY,NX)=real((/0.075,0.125,0.550,0.250/),r8)

    ENDIF
!
    !     PREVIOUS COARSE (K=0) AND FINE (K=1) ROOTS
!
  ELSE
    CFOSC(1:4,0,L,NY,NX)=real((/0.00,0.00,0.20,0.80/),r8)
    CFOSC(1:4,1,L,NY,NX)=real((/0.02,0.06,0.34,0.58/),r8)
  ENDIF
  end subroutine InitSurfResiduKinetiComponent
!------------------------------------------------------------------------------------------
  subroutine InitManureKinetiComponent(L,NY,NX)

  implicit none
  integer, intent(in) :: L, NY, NX
  !     begin_execution
  !
  !     RUMINANT
!
  IF(IXTYP(2,NY,NX).EQ.1)THEN
    CFOSC(1:4,2,L,NY,NX)=real((/0.036,0.044,0.630,0.290/),r8)

!
    !     NON-RUMINANT
!
  ELSEIF(IXTYP(2,NY,NX).EQ.2)THEN
    CFOSC(1:4,2,L,NY,NX)=real((/0.138,0.401,0.316,0.145/),r8)
!
!     OTHER
!
  ELSE
    CFOSC(1:4,2,L,NY,NX)=real((/0.138,0.401,0.316,0.145/),r8)
  ENDIF
  end subroutine InitManureKinetiComponent
!------------------------------------------------------------------------------------------
  subroutine InitPOMKinetiComponent(L,NY,NX,HCX,TORGL,CDPTHG,FCX,CORGCM)
  implicit none
  integer, intent(in)  :: L, NY, NX
  real(r8), intent(in) :: HCX
  real(r8), intent(in) :: TORGL
  real(r8), intent(in) :: CDPTHG
  real(r8), intent(out):: FCX
  real(r8), intent(out):: CORGCM
  real(r8) :: FCO,FCY,FC0,FC1

! begin_execution

!
! CFOSC=siNGLe kinetic fraction in POM
!
  IF(L.NE.0)THEN
    CFOSC(1:4,3,L,NY,NX)=real((/1.0,0.0,0.0,0.0/),r8)

!
!  HUMUS PARTITIONED TO DIFFERENT FRACTIONS
!  BASED ON SOC ACCUMULATION ABOVE EACH LAYER
!
!  NATURAL SOILS
!
    IF(ISOILR(NY,NX).EQ.0)THEN
!
      !     DRYLAND SOIL
      !
      !     CORGC,FORGC=SOC,minimum SOC for organic soil(g Mg-1)
      !     DPTH,DTBLZ=depth to layer midpoint,external water table(m)
      !     FC0=partitioning to less resistant component at DPTH=0
      !     FCX=reduction in FC0 at DPTH
      !     CORGCX,CORGNX,CORGPX=C,N,P concentations in humus
!
      IF(CORGC(L,NY,NX).LE.FORGC.OR.DPTH(L,NY,NX).LE.DTBLZ(NY,NX) &
        +CDPTH(NU(NY,NX),NY,NX)-CDPTHG)THEN
        FCY=0.60
        IF(CORGCX(4).GT.1.0E-32)THEN
          FC0=FCY*EXP(-5.0*(AMIN1(CORGNX(4),10.0*CORGPX(4))/CORGCX(4)))
        ELSE
          FCO=FCY
        ENDIF
        FCX=EXP(HCX*TORGL)
        !     WETLAND
!
      ELSE
        FCY=0.60
        IF(CORGCX(4).GT.1.0E-32)THEN
          FC0=FCY*EXP(-5.0*(AMIN1(CORGNX(4),10.0*CORGPX(4))/CORGCX(4)))
        ELSE
          FCO=FCY
        ENDIF
!     FCX=(EXP(HCX*TORGL))**0.5
        FCX=EXP(HCX*TORGL)
      ENDIF
!
      !     RECONSTRUCTED SOILS
!
    ELSE
      FCY=0.30
      IF(CORGCX(4).GT.1.0E-32)THEN
        FC0=FCY*EXP(-5.0*(AMIN1(CORGNX(4),10.0*CORGPX(4))/CORGCX(4)))
      ELSE
        FCO=FCY
      ENDIF
      FCX=1.0
    ENDIF
!
!   PARTITION HUMUS
!
!   CFOSC=fraction of humus in less(1),more(2) resistant component
!
    FC1=FC0*FCX
    CFOSC(1,4,L,NY,NX)=FC1
    CFOSC(2,4,L,NY,NX)=1.0-FC1
    CFOSC(3,4,L,NY,NX)=0.00_r8
    CFOSC(4,4,L,NY,NX)=0.00_r8
!
!   MICROBIAL DETRITUS ALLOCATED TO HUMUS MAINTAINS
!   HUMUS PARTITIONING TO COMPONENTS
!
!   CFOMC=fraction of microbial litter allocated to humus components
!
    CFOMC(1,L,NY,NX)=3.0*FC1/(2.0*FC1+1.0)
    CFOMC(2,L,NY,NX)=1.0-CFOMC(1,L,NY,NX)
  ENDIF

  IF(L.GT.0)THEN
    IF(BKDS(L,NY,NX).GT.ZERO)THEN
      CORGCM=AMIN1(0.55E+06_r8,(CORGCX(1)+CORGCX(2)+CORGCX(3)+CORGCX(4)))/0.55_r8
    else
      CORGCM=0._r8
    endif
  endif
  end subroutine InitPOMKinetiComponent

!------------------------------------------------------------------------------------------

  subroutine InitSOMProfile(L,NY,NX,HCX,TORGL,CDPTHG,CORGCM,FCX)

  implicit none
  integer,  intent(in)  :: L,NY,NX
  real(r8), intent(in)  :: HCX
  real(r8), intent(in)  :: TORGL
  real(r8), intent(in)  :: CDPTHG
  real(r8), intent(out) :: CORGCM
  real(r8), intent(out) :: FCX

  call InitLitterProfile(L,NY,NX)

  !     SURFACE RESIDUE KINETIC COMPONENTS
  call InitSurfResiduKinetiComponent(L,NY,NX)
  !
  !     ANIMAL MANURE
  call InitManureKinetiComponent(L,NY,NX)
  !
  !     POM
  call InitPOMKinetiComponent(L,NY,NX,HCX,TORGL,CDPTHG,FCX,CORGCM)

  end subroutine InitSOMProfile

!------------------------------------------------------------------------------------------

  subroutine InitLitterProfile(L,NY,NX)

  implicit none
  integer, intent(in) :: L,NY,NX

  real(r8) :: CORGCZ,CORGNZ,CORGPZ,CORGRZ
  real(r8) :: scal
! begin_execution
  associate(                  &
    CNRH    => micpar%CNRH   ,&
    CPRH    => micpar%CPRH    &
  )
  IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
    scal=AREA(3,L,NY,NX)/BKVL(L,NY,NX)
    CORGCX(0:2)=RSC(0:2,L,NY,NX)*scal
    CORGNX(0:2)=RSN(0:2,L,NY,NX)*scal
    CORGPX(0:2)=RSP(0:2,L,NY,NX)*scal
  ELSE
    scal=AREA(3,L,NY,NX)/VOLT(L,NY,NX)
    CORGCX(0:2)=RSC(0:2,L,NY,NX)*scal
    CORGNX(0:2)=RSN(0:2,L,NY,NX)*scal
    CORGPX(0:2)=RSP(0:2,L,NY,NX)*scal
  ENDIF
    !
    !     ALLOCATE SOC TO POC(3) AND HUMUS(4)
    !
    IF(L.GT.0)THEN
      CORGCZ=CORGC(L,NY,NX)
      CORGRZ=CORGR(L,NY,NX)
      CORGNZ=CORGN(L,NY,NX)
      CORGPZ=CORGP(L,NY,NX)
      IF(CORGCZ.GT.ZERO)THEN
        CORGCX(3)=CORGRZ
        CORGCX(4)=AZMAX1(CORGCZ-CORGCX(3))
        CORGNX(3)=AMIN1(CNRH(3)*CORGCX(3),CORGNZ)
        CORGNX(4)=AZMAX1(CORGNZ-CORGNX(3))
        CORGPX(3)=AMIN1(CPRH(3)*CORGCX(3),CORGPZ)
        CORGPX(4)=AZMAX1(CORGPZ-CORGPX(3))
      ELSE
        CORGCX(3)=0.0_r8
        CORGCX(4)=0.0_r8
        CORGNX(3)=0.0_r8
        CORGNX(4)=0.0_r8
        CORGPX(3)=0.0_r8
        CORGPX(4)=0.0_r8
      ENDIF
    ELSE
      CORGCX(3)=0.0_r8
      CORGCX(4)=0.0_r8
      CORGNX(3)=0.0_r8
      CORGNX(4)=0.0_r8
      CORGPX(3)=0.0_r8
      CORGPX(4)=0.0_r8
    ENDIF
  end associate
  end subroutine InitLitterProfile


!------------------------------------------------------------------------------------------

  subroutine DestructSOMBGC
  use abortutils, only : destroy
  implicit none

  call destroy(CORGCX)
  call destroy(CORGNX)
  call destroy(CORGPX)

  end subroutine DestructSOMBGC
end module InitSOMBGCMOD
