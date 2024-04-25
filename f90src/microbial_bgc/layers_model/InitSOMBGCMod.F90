module InitSOMBGCMOD
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use MicrobialDataType
  use SOMDataType
  use GridConsts
  use minimathmod, only : AZMAX1
  use SoilPhysDataType
  use FlagDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use EcosimConst
  use EcoSIMConfig, only : nlbiomcp => NumLiveMicrbCompts, ndbiomcp=> NumDeadMicrbCompts
  use SurfLitterDataType
  use SoilPropertyDataType
  use GridDataType
  use EcoSiMParDataMod, only : micpar

  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__
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
  NumMicrobAutrophCmplx= micpar%NumMicrobAutrophCmplx
  NumMicrbHetetrophCmplx= micpar%NumMicrbHetetrophCmplx
  NumLiveHeterBioms =micpar%NumLiveHeterBioms
  NumLiveAutoBioms  =micpar%NumLiveAutoBioms
  allocate(CORGCX(1:jcplx))
  allocate(CORGNX(1:jcplx))
  allocate(CORGPX(1:jcplx))
  end subroutine InitSOMBGC

!------------------------------------------------------------------------------------------

  subroutine InitSOMConsts
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
  real(r8) :: OSCI(1:jcplx),OSNI(1:jcplx),OSPI(1:jcplx)
  real(r8) :: OSCM(1:jcplx)
  real(r8) :: TOSCK(1:jcplx),TOSNK(1:jcplx),TOSPK(1:jcplx)
  real(r8) :: RC
  real(r8) :: CNOSCT(1:jcplx), CPOSCT(1:jcplx)
  real(r8) :: OSCX(1:jcplx)
  real(r8) :: OSNX(1:jcplx)
  real(r8) :: OSPX(1:jcplx)
  real(r8) :: tglds
  integer :: MID
  ! begin_execution

  associate(                  &
    rNCOMCa  => micpar%rNCOMCa ,&
    rPCOMCa               => micpar%rPCOMCa               , &
    nlbiomcp              => micpar%nlbiomcp              , &
    NumMicrobAutrophCmplx => micpar%NumMicrobAutrophCmplx , &
    k_humus => micpar%k_humus, &
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

  D975: DO K=1,micpar%NumOfLitrCmplxs
    CNOSCT(K)=0.0_r8
    CPOSCT(K)=0.0_r8
    IF(RSC(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      RNT=0.0_r8
      RPT=0.0_r8
      D970: DO M=1,jsken
        RNT=RNT+RSC(K,L,NY,NX)*CFOSC(M,K,L,NY,NX)*CNOFC(M,K)
        RPT=RPT+RSC(K,L,NY,NX)*CFOSC(M,K,L,NY,NX)*CPOFC(M,K)
      ENDDO D970
      FRNT=RSN(K,L,NY,NX)/RNT
      FRPT=RSP(K,L,NY,NX)/RPT
      D960: DO M=1,jsken
        CNOSC(M,K,L,NY,NX)=CNOFC(M,K)*FRNT
        CPOSC(M,K,L,NY,NX)=CPOFC(M,K)*FRPT
        CNOSCT(K)=CNOSCT(K)+CFOSC(M,K,L,NY,NX)*CNOSC(M,K,L,NY,NX)
        CPOSCT(K)=CPOSCT(K)+CFOSC(M,K,L,NY,NX)*CPOSC(M,K,L,NY,NX)
      ENDDO D960
    ELSE
      D965: DO M=1,jsken
        CNOSC(M,K,L,NY,NX)=CNRH(K)
        CPOSC(M,K,L,NY,NX)=CPRH(K)
      ENDDO D965
      CNOSCT(K)=CNRH(K)
      CPOSCT(K)=CPRH(K)
    ENDIF
  ENDDO D975

  D990: DO K=micpar%NumOfLitrCmplxs+1,jcplx1
    CNOSCT(K)=0.0_r8
    CPOSCT(K)=0.0_r8
    IF(CORGCX(K).GT.ZERO)THEN
      D985: DO M=1,jsken
        CNOSC(M,K,L,NY,NX)=CORGNX(K)/CORGCX(K)
        CPOSC(M,K,L,NY,NX)=CORGPX(K)/CORGCX(K)
        CNOSCT(K)=CNOSCT(K)+CFOSC(M,K,L,NY,NX)*CNOSC(M,K,L,NY,NX)
        CPOSCT(K)=CPOSCT(K)+CFOSC(M,K,L,NY,NX)*CPOSC(M,K,L,NY,NX)
      ENDDO D985
    ELSE
      D980: DO M=1,jsken
        CNOSC(M,K,L,NY,NX)=CNRH(K)
        CPOSC(M,K,L,NY,NX)=CPRH(K)
      ENDDO D980
      CNOSCT(K)=CNRH(K)
      CPOSCT(K)=CPRH(K)
    ENDIF
  ENDDO D990
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
  D995: DO K=1,jcplx
    IF(L.EQ.0)THEN
      KK=K
    ELSE
      !humus complex
      KK=micpar%k_humus
    ENDIF
    IF(SoilMicPMassLayer(L,NY,NX).GT.ZEROS(NY,NX))THEN
      OSCI(K)=CORGCX(K)*SoilMicPMassLayer(L,NY,NX)
      OSNI(K)=CORGNX(K)*SoilMicPMassLayer(L,NY,NX)
      OSPI(K)=CORGPX(K)*SoilMicPMassLayer(L,NY,NX)
    ELSE
      OSCI(K)=CORGCX(K)*VGeomLayer(L,NY,NX)
      OSNI(K)=CORGNX(K)*VGeomLayer(L,NY,NX)
      OSPI(K)=CORGPX(K)*VGeomLayer(L,NY,NX)
    ENDIF
    TOSCK(K)=OMCK(K)+ORCK(K)+OQCK(K)+OHCK(K)
    TOSNK(K)=ORCK(K)*CNRH(K)+OQCK(K)*CNOSCT(KK)+OHCK(K)*CNOSCT(KK)
    TOSPK(K)=ORCK(K)*CPRH(K)+OQCK(K)*CPOSCT(KK)+OHCK(K)*CPOSCT(KK)
!   based on aerobic heterotrophs

    TOSNK(K)=TOSNK(K)+OMCI(1,K)*rNCOMCa(1,1,K)+OMCI(2,K)*rNCOMCa(2,1,K)
    TOSPK(K)=TOSPK(K)+OMCI(1,K)*rPCOMCa(1,1,K)+OMCI(2,K)*rPCOMCa(2,1,K)

    TOSCI=TOSCI+OSCI(K)*TOSCK(K)
    TOSNI=TOSNI+OSCI(K)*TOSNK(K)
    TOSPI=TOSPI+OSCI(K)*TOSPK(K)
    OSCX(K)=0.0_r8
    OSNX(K)=0.0_r8
    OSPX(K)=0.0_r8
  ENDDO D995

  D8995: DO K=1,jcplx
    IF(L.EQ.0)THEN
      OSCM(K)=DCKR*CORGCX(K)*SoilMicPMassLayer(L,NY,NX)
      X=0.0_r8
      KK=K
      FOSCI=1.0
      FOSNI=1.0
      FOSPI=1.0
    ELSE
      IF(SoilMicPMassLayer(L,NY,NX).GT.ZEROS(NY,NX))THEN
        IF(K.LE.micpar%NumOfLitrCmplxs)THEN
          OSCM(K)=DCKR*CORGCX(K)*SoilMicPMassLayer(L,NY,NX)
        ELSE
          OSCM(K)=FCX*CORGCX(K)*SoilMicPMassLayer(L,NY,NX)*DCKM/(CORGCX(k_humus)+DCKM)
        ENDIF
      ELSE
        IF(K.LE.micpar%NumOfLitrCmplxs)THEN
          OSCM(K)=DCKR*CORGCX(K)*VGeomLayer(L,NY,NX)
        ELSE
          OSCM(K)=FCX*CORGCX(K)*VGeomLayer(L,NY,NX)*DCKM/(CORGCX(k_humus)+DCKM)
        ENDIF
      ENDIF
      X=1.0_r8
      KK=micpar%k_humus
      IF(TOSCI.GT.ZEROS(NY,NX))THEN
        FOSCI=AMIN1(1.0_r8,OSCI(KK)/TOSCI)
        FOSNI=AMIN1(1.0_r8,OSCI(KK)*CNOSCT(KK)/TOSNI)
        FOSPI=AMIN1(1.0_r8,OSCI(KK)*CPOSCT(KK)/TOSPI)
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
!     rNCOMC,rPCOMC=maximum N:C and P:C ratios in microbial biomass
!     OSCX,OSNX,OSPX=remaining unallocated SOC,SON,SOP
!  The reason that initialization of complex 5 microbes is repated for each
! complex is because complex 5 is shared by the other complexes
    OMEauto(1:NumPlantChemElms,1:NumLiveAutoBioms,L,NY,NX)=0._r8

    D8990: DO N=1,NumMicbFunGroups
      tglds=JGnfo(N)-JGnio(N)+1._r8
      D8991: DO M=1,nlbiomcp
        OMC1=AZMAX1(OSCM(K)*OMCI(M,K)*OMCF(N)*FOSCI)
        OMN1=AZMAX1(OMC1*rNCOMCa(M,N,K)*FOSNI)
        OMP1=AZMAX1(OMC1*rPCOMCa(M,N,K)*FOSPI)
        do NGL=JGnio(N),JGnfo(N)
          MID=micpar%get_micb_id(M,NGL)
          OMEheter(ielmc,MID,K,L,NY,NX)=OMC1/tglds
          OMEheter(ielmn,MID,K,L,NY,NX)=OMN1/tglds
          OMEheter(ielmp,MID,K,L,NY,NX)=OMP1/tglds
        ENDDO
        OSCX(KK)=OSCX(KK)+OMC1
        OSNX(KK)=OSNX(KK)+OMN1
        OSPX(KK)=OSPX(KK)+OMP1
        D8992: DO NN=1,NumMicbFunGroups
          tglds=JGnfA(N)-JGniA(N)+1._r8
          do NGL=JGniA(N),JGnfA(N)
            MID=micpar%get_micb_id(M,NGL)
            OMEauto(ielmc,MID,L,NY,NX)=OMEauto(ielmc,MID,L,NY,NX)+OMC1*OMCA(NN)/tglds
            OMEauto(ielmn,MID,L,NY,NX)=OMEauto(ielmn,MID,L,NY,NX)+OMN1*OMCA(NN)/tglds
            OMEauto(ielmp,MID,L,NY,NX)=OMEauto(ielmp,MID,L,NY,NX)+OMP1*OMCA(NN)/tglds
          ENDDO
          OSCX(KK)=OSCX(KK)+OMC1*OMCA(NN)
          OSNX(KK)=OSNX(KK)+OMN1*OMCA(NN)
          OSPX(KK)=OSPX(KK)+OMP1*OMCA(NN)
        ENDDO D8992
      ENDDO D8991
    ENDDO D8990
!
!     MICROBIAL RESIDUE C, N AND P
!
!     ORC,ORN,ORP=microbial residue C,N,P
!     ORCI=allocation of microbial residue to kinetic components
!  X is an indicator of surface residual layer
    D8985: DO M=1,ndbiomcp
      ORM(ielmc,M,K,L,NY,NX)=X*AZMAX1(OSCM(K)*ORCI(M,K)*FOSCI)
      ORM(ielmn,M,K,L,NY,NX)=AZMAX1(ORM(ielmc,M,K,L,NY,NX)*rNCOMCa(M,1,K)*FOSNI)
      ORM(ielmp,M,K,L,NY,NX)=AZMAX1(ORM(ielmc,M,K,L,NY,NX)*rPCOMCa(M,1,K)*FOSPI)
      OSCX(KK)=OSCX(KK)+ORM(ielmc,M,K,L,NY,NX)
      OSNX(KK)=OSNX(KK)+ORM(ielmn,M,K,L,NY,NX)
      OSPX(KK)=OSPX(KK)+ORM(ielmp,M,K,L,NY,NX)
    ENDDO D8985
!
!     DOC, DON AND DOP
!
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores (g)
!     OQCH,OQNH,OQPH,OQAH=DOC,DON,DOP,acetate in macropores (g)
!
    DOM(idom_doc,K,L,NY,NX)=X*AZMAX1(OSCM(K)*OQCK(K)*FOSCI)
    DOM(idom_don,K,L,NY,NX)=AZMAX1(DOM(idom_doc,K,L,NY,NX)*CNOSCT(KK)*FOSNI)
    DOM(idom_dop,K,L,NY,NX)=AZMAX1(DOM(idom_doc,K,L,NY,NX)*CPOSCT(KK)*FOSPI)
    DOM(idom_acetate,K,L,NY,NX)=0.0_r8
    DOM_Macp(idom_doc,K,L,NY,NX)=0.0_r8
    DOM_Macp(idom_don,K,L,NY,NX)=0.0_r8
    DOM_Macp(idom_dop,K,L,NY,NX)=0.0_r8
    DOM_Macp(idom_acetate,K,L,NY,NX)=0.0_r8
    OSCX(KK)=OSCX(KK)+DOM(idom_doc,K,L,NY,NX)
    OSNX(KK)=OSNX(KK)+DOM(idom_don,K,L,NY,NX)
    OSPX(KK)=OSPX(KK)+DOM(idom_dop,K,L,NY,NX)
!
!     ADSORBED C, N AND P
!
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!
    OHM(ielmc,K,L,NY,NX)=X*AZMAX1(OSCM(K)*OHCK(K)*FOSCI)
    OHM(ielmn,K,L,NY,NX)=AZMAX1(OHM(ielmc,K,L,NY,NX)*CNOSCT(KK)*FOSNI)
    OHM(ielmp,K,L,NY,NX)=AZMAX1(OHM(ielmc,K,L,NY,NX)*CPOSCT(KK)*FOSPI)
    OHM(idom_acetate,K,L,NY,NX)=0.0_r8
    OSCX(KK)=OSCX(KK)+OHM(ielmc,K,L,NY,NX)+OHM(idom_acetate,K,L,NY,NX)
    OSNX(KK)=OSNX(KK)+OHM(ielmn,K,L,NY,NX)
    OSPX(KK)=OSPX(KK)+OHM(ielmp,K,L,NY,NX)
!
!     HUMUS C, N AND P
!
!     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP

    D8980: DO M=1,jsken
      OSM(ielmc,M,K,L,NY,NX)=AZMAX1(CFOSC(M,K,L,NY,NX)*(OSCI(K)-OSCX(K)))
      IF(CNOSCT(K).GT.ZERO)THEN
        OSM(ielmn,M,K,L,NY,NX)=AZMAX1(CFOSC(M,K,L,NY,NX)*CNOSC(M,K,L,NY,NX) &
          /CNOSCT(K)*(OSNI(K)-OSNX(K)))
      ELSE
        OSM(ielmn,M,K,L,NY,NX)=0.0_r8
      ENDIF
      IF(CPOSCT(K).GT.ZERO)THEN
        OSM(ielmp,M,K,L,NY,NX)=AZMAX1(CFOSC(M,K,L,NY,NX)*CPOSC(M,K,L,NY,NX) &
          /CPOSCT(K)*(OSPI(K)-OSPX(K)))
      ELSE
        OSM(ielmp,M,K,L,NY,NX)=0.0_r8
      ENDIF
      IF(K.EQ.micpar%k_woody_litr)THEN
        OSA(M,K,L,NY,NX)=OSM(ielmc,M,K,L,NY,NX)*OMCI(1,K)
      ELSE
        OSA(M,K,L,NY,NX)=OSM(ielmc,M,K,L,NY,NX)
      ENDIF
    ENDDO D8980
  ENDDO D8995
!
!     ADD ALL LITTER,POC,HUMUS COMPONENTS TO GET TOTAL SOC
!
  OC=0.0_r8
  ON=0.0_r8
  OP=0.0_r8
  RC=0.0_r8
  IF(L.EQ.0)THEN
    RC0(:,NY,NX)=0.0_r8
    RC0ff(NY,NX)=0._r8
  ENDIF

  ROXYS(:,:,L,NY,NX)=0.0_r8
  RVMX4(:,:,L,NY,NX)=0.0_r8
  RVMX3(:,:,L,NY,NX)=0.0_r8
  RVMX2(:,:,L,NY,NX)=0.0_r8
  RVMX1(:,:,L,NY,NX)=0.0_r8
  RINHO(:,:,L,NY,NX)=0.0_r8
  RINHB(:,:,L,NY,NX)=0.0_r8
  RINOO(:,:,L,NY,NX)=0.0_r8
  RIPOO(:,:,L,NY,NX)=0.0_r8
  RINOB(:,:,L,NY,NX)=0.0_r8
  RIPBO(:,:,L,NY,NX)=0.0_r8
  RIPO1(:,:,L,NY,NX)=0.0_r8
  RIPB1(:,:,L,NY,NX)=0.0_r8
  IF(L.EQ.0)THEN
    RINHOR(:,:,NY,NX)=0.0_r8
    RINOOR(:,:,NY,NX)=0.0_r8
    RIPOOR(:,:,NY,NX)=0.0_r8
  ENDIF

  D6990: DO K=1,jcplx
    DO  N=1,NumMicbFunGroups
      do NGL=JGnio(n),JGnfo(n)
        DO  M=1,nlbiomcp
          OC=OC+OMEheter(ielmc,MID,K,L,NY,NX)
          ON=ON+OMEheter(ielmn,MID,K,L,NY,NX)
          OP=OP+OMEheter(ielmp,MID,K,L,NY,NX)
          IF(K.LE.micpar%NumOfLitrCmplxs)THEN
            RC=RC+OMEheter(ielmc,MID,K,L,NY,NX)
          ENDIF
          RC0(K,NY,NX)=RC0(K,NY,NX)+OMEheter(ielmc,MID,K,L,NY,NX)
        ENDDO
      ENDDO
    enddo
  ENDDO D6990

  ROXYSff(:,L,NY,NX)=0.0_r8
  RNH3OxidAutor(:,L,NY,NX)=0.0_r8
  RVMX3ff(:,L,NY,NX)=0.0_r8
  RVMB3ff(:,L,NY,NX)=0.0_r8
  RNO2OxidAutor(:,L,NY,NX)=0.0_r8
  RVMX1ff(:,L,NY,NX)=0.0_r8
  RINHOff(:,L,NY,NX)=0.0_r8
  RINHBff(:,L,NY,NX)=0.0_r8
  RINOOff(:,L,NY,NX)=0.0_r8
  RINOBff(:,L,NY,NX)=0.0_r8
  RIPOOff(:,L,NY,NX)=0.0_r8
  RIPBOff(:,L,NY,NX)=0.0_r8
  RIPO1ff(:,L,NY,NX)=0.0_r8
  RIPB1ff(:,L,NY,NX)=0.0_r8
  IF(L.EQ.0)THEN
    RINHORff(:,NY,NX)=0.0_r8
    RINOORff(:,NY,NX)=0.0_r8
    RIPOORff(:,NY,NX)=0.0_r8
  ENDIF
    DO  N=1,NumMicbFunGroups
      do NGL=JGniA(n),JGnfA(n)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          OC=OC+OMEauto(ielmc,MID,L,NY,NX)
          ON=ON+OMEauto(ielmn,MID,L,NY,NX)
          OP=OP+OMEauto(ielmp,MID,L,NY,NX)
          RC0ff(NY,NX)=RC0ff(NY,NX)+OMEauto(ielmc,MID,L,NY,NX)
        ENDDO
      ENDDO
    enddo

  D6995: DO K=1,jcplx
    D6985: DO M=1,ndbiomcp
      OC=OC+ORM(ielmc,M,K,L,NY,NX)
      ON=ON+ORM(ielmn,M,K,L,NY,NX)
      OP=OP+ORM(ielmp,M,K,L,NY,NX)
      IF(K.LE.micpar%NumOfLitrCmplxs)THEN
        RC=RC+ORM(ielmc,M,K,L,NY,NX)
      ENDIF
      IF(L.EQ.0)THEN
        RC0(K,NY,NX)=RC0(K,NY,NX)+ORM(ielmc,M,K,L,NY,NX)
      ENDIF
    ENDDO D6985
    OC=OC+DOM(idom_doc,K,L,NY,NX)+DOM_Macp(idom_doc,K,L,NY,NX)+OHM(ielmc,K,L,NY,NX) &
      +DOM(idom_acetate,K,L,NY,NX)+DOM_Macp(idom_acetate,K,L,NY,NX)+OHM(idom_acetate,K,L,NY,NX)
    ON=ON+DOM(idom_don,K,L,NY,NX)+DOM_Macp(idom_don,K,L,NY,NX)+OHM(ielmn,K,L,NY,NX)
    OP=OP+DOM(idom_dop,K,L,NY,NX)+DOM_Macp(idom_dop,K,L,NY,NX)+OHM(ielmp,K,L,NY,NX)
    OC=OC+DOM(idom_acetate,K,L,NY,NX)+DOM_Macp(idom_acetate,K,L,NY,NX)
    IF(K.LE.micpar%NumOfLitrCmplxs)THEN
      RC=RC+DOM(idom_doc,K,L,NY,NX)+DOM_Macp(idom_doc,K,L,NY,NX)+OHM(ielmc,K,L,NY,NX) &
        +DOM(idom_acetate,K,L,NY,NX)+DOM_Macp(idom_acetate,K,L,NY,NX)+OHM(idom_acetate,K,L,NY,NX)
      RC=RC+DOM(idom_acetate,K,L,NY,NX)+DOM_Macp(idom_acetate,K,L,NY,NX)
    ENDIF
    IF(L.EQ.0)THEN
      RC0(K,NY,NX)=RC0(K,NY,NX)+DOM(idom_doc,K,L,NY,NX)+DOM_Macp(idom_doc,K,L,NY,NX) &
        +OHM(ielmc,K,L,NY,NX)+DOM(idom_acetate,K,L,NY,NX)+DOM_Macp(idom_acetate,K,L,NY,NX)+OHM(idom_acetate,K,L,NY,NX)
    ENDIF
    D6980: DO M=1,jsken
      OC=OC+OSM(ielmc,M,K,L,NY,NX)
      ON=ON+OSM(ielmn,M,K,L,NY,NX)
      OP=OP+OSM(ielmp,M,K,L,NY,NX)
      IF(K.LE.micpar%NumOfLitrCmplxs)THEN
        RC=RC+OSM(ielmc,M,K,L,NY,NX)
      ENDIF
      IF(L.EQ.0)THEN
        RC0(K,NY,NX)=RC0(K,NY,NX)+OSM(ielmc,M,K,L,NY,NX)
      ENDIF
    ENDDO D6980
  ENDDO D6995
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
  associate(                               &
   k_woody_litr  => micpar%k_woody_litr  , &
   k_fine_litr   => micpar%k_fine_litr   , &
   k_manure      => micpar%k_manure        &
  )
  IF(L.EQ.0)THEN
    !
    !     CFOSC=fraction of litter in protein(1),nonstructural(2)
    !     cellulose(3) and lignin(4)
    !
    !     PREVIOUS COARSE WOODY RESIDUE
    !
    CFOSC(1:jsken,k_woody_litr,L,NY,NX)=real((/0.00,0.045,0.660,0.295/),r8)

    !
    !     MAIZE
    !
    IF(IXTYP(1,NY,NX).EQ.1)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.080,0.245,0.613,0.062/),r8)
      !
      !     WHEAT
      !
    ELSEIF(IXTYP(1,NY,NX).EQ.2)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.125,0.171,0.560,0.144/),r8)
!
      !     SOYBEAN
      !
    ELSEIF(IXTYP(1,NY,NX).EQ.3)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.138,0.426,0.316,0.120/),r8)

      !     NEW STRAW
!
    ELSEIF(IXTYP(1,NY,NX).EQ.4)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.036,0.044,0.767,0.153/),r8)

      !     OLD STRAW
!
    ELSEIF(IXTYP(1,NY,NX).EQ.5)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.075,0.125,0.550,0.250/),r8)

!
      !     COMPOST
!
    ELSEIF(IXTYP(1,NY,NX).EQ.6)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.143,0.015,0.640,0.202/),r8)
!
      !     GREEN MANURE
!
    ELSEIF(IXTYP(1,NY,NX).EQ.7)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.202,0.013,0.560,0.225/),r8)

      !     NEW DECIDUOUS FOREST
!
    ELSEIF(IXTYP(1,NY,NX).EQ.8)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.070,0.41,0.36,0.16/),r8)

!
      !     NEW CONIFEROUS FOREST
!
    ELSEIF(IXTYP(1,NY,NX).EQ.9)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.07,0.25,0.38,0.30/),r8)

      !     OLD DECIDUOUS FOREST
!
    ELSEIF(IXTYP(1,NY,NX).EQ.10)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.02,0.06,0.34,0.58/),r8)
!
      !     OLD CONIFEROUS FOREST
!
    ELSEIF(IXTYP(1,NY,NX).EQ.11)THEN
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.02,0.06,0.34,0.58/),r8)
      !     DEFAULT
!
    ELSE
      CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.075,0.125,0.550,0.250/),r8)

    ENDIF
!
    !     PREVIOUS COARSE (K=0) AND FINE (K=1) ROOTS
!
  ELSE
    CFOSC(1:jsken,k_woody_litr,L,NY,NX)=real((/0.00,0.00,0.20,0.80/),r8)
    CFOSC(1:jsken,k_fine_litr,L,NY,NX)=real((/0.02,0.06,0.34,0.58/),r8)
  ENDIF
  end associate
  end subroutine InitSurfResiduKinetiComponent
!------------------------------------------------------------------------------------------
  subroutine InitManureKinetiComponent(L,NY,NX)

  implicit none
  integer, intent(in) :: L, NY, NX
  !     begin_execution
  associate(                      &
    k_manure  => micpar%k_manure, &
    jsken     => micpar%jsken     &
  )
  !
  !     RUMINANT
!
  IF(IXTYP(2,NY,NX).EQ.1)THEN
    CFOSC(1:jsken,k_manure,L,NY,NX)=real((/0.036,0.044,0.630,0.290/),r8)

!
    !     NON-RUMINANT
!
  ELSEIF(IXTYP(2,NY,NX).EQ.2)THEN
    CFOSC(1:jsken,k_manure,L,NY,NX)=real((/0.138,0.401,0.316,0.145/),r8)
!
!     OTHER
!
  ELSE
    CFOSC(1:jsken,k_manure,L,NY,NX)=real((/0.138,0.401,0.316,0.145/),r8)
  ENDIF
  end associate
  end subroutine InitManureKinetiComponent
!------------------------------------------------------------------------------------------
  subroutine InitPOMKinetiComponent(L,NY,NX,HCX,TORGL,LandScape1stSoiLayDepth,FCX,CORGCM)
  implicit none
  integer, intent(in)  :: L, NY, NX
  real(r8), intent(in) :: HCX
  real(r8), intent(in) :: TORGL
  real(r8), intent(in) :: LandScape1stSoiLayDepth
  real(r8), intent(out):: FCX
  real(r8), intent(out):: CORGCM
  real(r8) :: FCY,FC0,FC1

! begin_execution
  associate(                      &
    k_POM => micpar%k_POM      ,  &
    k_humus  => micpar%k_humus ,  &
    iprotein => micpar%iprotein,  &
    icarbhyro=> micpar%icarbhyro, &
    icellulos=> micpar%icellulos, &
    ilignin  => micpar%ilignin  , &
    jsken => micpar%jsken         &
  )
!
! CFOSC=siNGLe kinetic fraction in POM
!
  IF(L.NE.0)THEN
    CFOSC(1:jsken,k_POM,L,NY,NX)=real((/1.0,0.0,0.0,0.0/),r8)

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
      !     DPTH,ExtWaterTablet0=depth to layer midpoint,external water table(m)
      !     FC0=partitioning to less resistant component at DPTH=0
      !     FCX=reduction in FC0 at DPTH
      !     CORGCX,CORGNX,CORGPX=C,N,P concentations in humus
!
      IF(CORGC(L,NY,NX).LE.FORGC.OR.SoiDepthMidLay(L,NY,NX).LE.ExtWaterTablet0(NY,NX) &
        +CumDepth2LayerBottom(NU(NY,NX),NY,NX)-LandScape1stSoiLayDepth)THEN
        FCY=0.60_r8
        IF(CORGCX(k_humus).GT.1.0E-32_r8)THEN
          FC0=FCY*EXP(-5.0_r8*(AMIN1(CORGNX(k_humus), &
            10.0_r8*CORGPX(k_humus))/CORGCX(k_humus)))
        ELSE
          FC0=FCY
        ENDIF
        FCX=EXP(HCX*TORGL)
        !     WETLAND
!
      ELSE
        FCY=0.60_r8
        IF(CORGCX(k_humus).GT.1.0E-32_r8)THEN
          FC0=FCY*EXP(-5.0_r8*(AMIN1(CORGNX(k_humus), &
            10.0_r8*CORGPX(k_humus))/CORGCX(k_humus)))
        ELSE
          FC0=FCY
        ENDIF
!     FCX=(EXP(HCX*TORGL))**0.5_r8
        FCX=EXP(HCX*TORGL)
      ENDIF
!
      !     RECONSTRUCTED SOILS
!
    ELSE
      FCY=0.30_r8
      IF(CORGCX(k_humus).GT.1.0E-32_r8)THEN
        FC0=FCY*EXP(-5.0_r8*(AMIN1(CORGNX(k_humus), &
          10.0_r8*CORGPX(k_humus))/CORGCX(k_humus)))
      ELSE
        FC0=FCY
      ENDIF
      FCX=1.0_r8
    ENDIF
!
!   PARTITION HUMUS
!
!   CFOSC=fraction of humus in less(1),more(2) resistant component
!
    FC1=FC0*FCX
    CFOSC(iprotein,k_humus,L,NY,NX)=FC1
    CFOSC(icarbhyro,k_humus,L,NY,NX)=1.0_r8-FC1
    CFOSC(icellulos,k_humus,L,NY,NX)=0.00_r8
    CFOSC(ilignin,k_humus,L,NY,NX)=0.00_r8
!
!   MICROBIAL DETRITUS ALLOCATED TO HUMUS MAINTAINS
!   HUMUS PARTITIONING TO COMPONENTS
!
!   CFOMC=fraction of microbial litter allocated to humus components
!
    CFOMC(1,L,NY,NX)=3.0_r8*FC1/(2.0_r8*FC1+1.0_r8)
    CFOMC(2,L,NY,NX)=1.0_r8-CFOMC(1,L,NY,NX)
  ENDIF

  IF(L.GT.0)THEN
    IF(SoiBulkDensity(L,NY,NX).GT.ZERO)THEN
      CORGCM=AMIN1(orgcden,(CORGCX(1)+CORGCX(2)+CORGCX(3)+CORGCX(k_humus)))/0.55_r8
    else
      CORGCM=0._r8
    endif
  endif
  end associate
  end subroutine InitPOMKinetiComponent

!------------------------------------------------------------------------------------------

  subroutine InitSOMProfile(L,NY,NX,HCX,TORGL,LandScape1stSoiLayDepth,CORGCM,FCX)

  implicit none
  integer,  intent(in)  :: L,NY,NX
  real(r8), intent(in)  :: HCX
  real(r8), intent(in)  :: TORGL
  real(r8), intent(in)  :: LandScape1stSoiLayDepth
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
  call InitPOMKinetiComponent(L,NY,NX,HCX,TORGL,LandScape1stSoiLayDepth,FCX,CORGCM)

  end subroutine InitSOMProfile

!------------------------------------------------------------------------------------------

  subroutine InitLitterProfile(L,NY,NX)

  implicit none
  integer, intent(in) :: L,NY,NX

  real(r8) :: CORGCZ,CORGNZ,CORGPZ,CORGRZ
  real(r8) :: scal
! begin_execution
  associate(                     &
    NumOfLitrCmplxs=> micpar%NumOfLitrCmplxs , &
    k_POM   => micpar%k_POM    , &
    k_humus => micpar%k_humus  , &
    CNRH    => micpar%CNRH     , &
    CPRH    => micpar%CPRH       &
  )
  IF(SoilMicPMassLayer(L,NY,NX).GT.ZEROS(NY,NX))THEN
    scal=AREA(3,L,NY,NX)/SoilMicPMassLayer(L,NY,NX)
    CORGCX(1:NumOfLitrCmplxs)=RSC(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGNX(1:NumOfLitrCmplxs)=RSN(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGPX(1:NumOfLitrCmplxs)=RSP(1:NumOfLitrCmplxs,L,NY,NX)*scal
  ELSE
    scal=AREA(3,L,NY,NX)/VGeomLayer(L,NY,NX)
    CORGCX(1:NumOfLitrCmplxs)=RSC(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGNX(1:NumOfLitrCmplxs)=RSN(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGPX(1:NumOfLitrCmplxs)=RSP(1:NumOfLitrCmplxs,L,NY,NX)*scal
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
        CORGCX(k_POM)=CORGRZ
        CORGCX(k_humus)=AZMAX1(CORGCZ-CORGCX(3))
        CORGNX(k_POM)=AMIN1(CNRH(3)*CORGCX(3),CORGNZ)
        CORGNX(k_humus)=AZMAX1(CORGNZ-CORGNX(3))
        CORGPX(k_POM)=AMIN1(CPRH(3)*CORGCX(3),CORGPZ)
        CORGPX(k_humus)=AZMAX1(CORGPZ-CORGPX(3))
      ELSE
        CORGCX(k_POM)=0.0_r8
        CORGCX(k_humus)=0.0_r8
        CORGNX(k_POM)=0.0_r8
        CORGNX(k_humus)=0.0_r8
        CORGPX(k_POM)=0.0_r8
        CORGPX(k_humus)=0.0_r8
      ENDIF
    ELSE
      CORGCX(k_POM)=0.0_r8
      CORGCX(k_humus)=0.0_r8
      CORGNX(k_POM)=0.0_r8
      CORGNX(k_humus)=0.0_r8
      CORGPX(k_POM)=0.0_r8
      CORGPX(k_humus)=0.0_r8
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
