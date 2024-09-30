module InitSOMBGCMOD
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use EcoSIMConfig,     only: nlbiomcp => NumLiveMicrbCompts, ndbiomcp=> NumDeadMicrbCompts
  use SoilBGCNLayMod,   only: sumorgmlayl,sumLitrOMLayL, sumMicBiomLayL
  use minimathmod,      only: AZMAX1
  use EcoSiMParDataMod, only: micpar
  use MicrobialDataType
  use SOMDataType
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use EcosimConst
  use SurfLitterDataType
  use SoilPropertyDataType
  use GridDataType
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
  NumHetetrMicCmplx= micpar%NumHetetrMicCmplx
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

  real(r8), parameter :: DCKR=0.25_r8
  real(r8), parameter :: DCKM=2.5E+04_r8

  integer :: K,M,KK,N,NGL,NN
  real(r8) :: OC,ON,OP,X
  real(r8) :: OME1(1:NumPlantChemElms)
  real(r8) :: FOSCI,FOSNI,FOSPI
  real(r8) :: RNT,RPT
  real(r8) :: FRNT,FRPT
  real(r8) :: TOSNI,TOSPI,TOSCI
  real(r8) :: OSCI(1:jcplx),OSNI(1:jcplx),OSPI(1:jcplx)
  real(r8) :: OSCM(1:jcplx)
  real(r8) :: TOSCK(1:jcplx),TOSNK(1:jcplx),TOSPK(1:jcplx)
  real(r8) :: CNOSCT(1:jcplx), CPOSCT(1:jcplx)
  real(r8) :: OSCX(1:jcplx)
  real(r8) :: OSNX(1:jcplx)
  real(r8) :: OSPX(1:jcplx)
  real(r8) :: litrOM(NumPlantChemElms)
  real(r8) :: ORGM(NumPlantChemElms)
  real(r8) :: tglds
  integer  :: MID,NE
  ! begin_execution

  associate(                                               &
    rNCOMCa               => micpar%rNCOMCa,               &
    rPCOMCa               => micpar%rPCOMCa,               &
    nlbiomcp              => micpar%nlbiomcp,              &
    NumMicrobAutrophCmplx => micpar%NumMicrobAutrophCmplx, &
    k_humus               => micpar%k_humus,               &
    OHCK                  => micpar%OHCK,                  &
    OMCK                  => micpar%OMCK,                  &
    OQCK                  => micpar%OQCK,                  &
    ORCK                  => micpar%ORCK,                  &
    ORCI                  => micpar%ORCI,                  &
    OMCI                  => micpar%OMCI,                  &
    CNRH                  => micpar%CNRH,                  &
    CPRH                  => micpar%CPRH,                  &
    CNOFC                 => micpar%CNOFC,                 &
    CPOFC                 => micpar%CPOFC,                 &
    OMCF                  => micpar%OMCF,                  &
    OMCA                  => micpar%OMCA                   &
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
        CNOSC(M,K,L,NY,NX) = CNOFC(M,K)*FRNT
        CPOSC(M,K,L,NY,NX) = CPOFC(M,K)*FRPT
        CNOSCT(K)          = CNOSCT(K)+CFOSC(M,K,L,NY,NX)*CNOSC(M,K,L,NY,NX)
        CPOSCT(K)          = CPOSCT(K)+CFOSC(M,K,L,NY,NX)*CPOSC(M,K,L,NY,NX)
      ENDDO D960
    ELSE
      D965: DO M=1,jsken
        CNOSC(M,K,L,NY,NX)=CNRH(K)
        CPOSC(M,K,L,NY,NX)=CPRH(K)
      ENDDO D965
      CNOSCT(K)=CNRH(K)
      CPOSCT(K)=CPRH(K)
    ENDIF
!    write(*,*)K,L,CNOSCT(K),CPOSCT(K)    
  ENDDO D975

  D990: DO K=micpar%NumOfLitrCmplxs+1,jcplx
    CNOSCT(K)=0.0_r8
    CPOSCT(K)=0.0_r8
    IF(CORGCX(K).GT.ZERO)THEN
      D985: DO M=1,jsken
        CNOSC(M,K,L,NY,NX) = CORGNX(K)/CORGCX(K)
        CPOSC(M,K,L,NY,NX) = CORGPX(K)/CORGCX(K)
        CNOSCT(K)          = CNOSCT(K)+CFOSC(M,K,L,NY,NX)*CNOSC(M,K,L,NY,NX)
        CPOSCT(K)          = CPOSCT(K)+CFOSC(M,K,L,NY,NX)*CPOSC(M,K,L,NY,NX)
      ENDDO D985
    ELSE
      D980: DO M=1,jsken
        CNOSC(M,K,L,NY,NX)=CNRH(K)
        CPOSC(M,K,L,NY,NX)=CPRH(K)
      ENDDO D980
      CNOSCT(K)=CNRH(K)
      CPOSCT(K)=CPRH(K)
    ENDIF
!    write(*,*)K,L,CNOSCT(K),CPOSCT(K),jcplx
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
    IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
      OSCI(K)=CORGCX(K)*VLSoilMicPMass_vr(L,NY,NX)
      OSNI(K)=CORGNX(K)*VLSoilMicPMass_vr(L,NY,NX)
      OSPI(K)=CORGPX(K)*VLSoilMicPMass_vr(L,NY,NX)
    ELSE
      OSCI(K)=CORGCX(K)*VGeomLayer_vr(L,NY,NX)
      OSNI(K)=CORGNX(K)*VGeomLayer_vr(L,NY,NX)
      OSPI(K)=CORGPX(K)*VGeomLayer_vr(L,NY,NX)
    ENDIF
!    write(*,*)L,K,CORGCX(K),VLSoilMicPMass_vr(L,NY,NX),VGeomLayer_vr(L,NY,NX)
    TOSCK(K)=OMCK(K)+ORCK(K)+OQCK(K)+OHCK(K)
    TOSNK(K)=ORCK(K)*CNRH(K)+OQCK(K)*CNOSCT(KK)+OHCK(K)*CNOSCT(KK)
    TOSPK(K)=ORCK(K)*CPRH(K)+OQCK(K)*CPOSCT(KK)+OHCK(K)*CPOSCT(KK)
!   based on aerobic heterotrophs

    TOSNK(K)=TOSNK(K)+OMCI(1,K)*rNCOMCa(1,1,K)+OMCI(2,K)*rNCOMCa(2,1,K)
    TOSPK(K)=TOSPK(K)+OMCI(1,K)*rPCOMCa(1,1,K)+OMCI(2,K)*rPCOMCa(2,1,K)

    TOSCI   = TOSCI+OSCI(K)*TOSCK(K)
    TOSNI   = TOSNI+OSCI(K)*TOSNK(K)
    TOSPI   = TOSPI+OSCI(K)*TOSPK(K)
    OSCX(K) = 0.0_r8
    OSNX(K) = 0.0_r8
    OSPX(K) = 0.0_r8
  ENDDO D995

  D8995: DO K=1,jcplx
    IF(L.EQ.0)THEN
      OSCM(K) = AMIN1(DCKR,1._r8)*CORGCX(K)*VLSoilMicPMass_vr(L,NY,NX)
      X       = 0.0_r8
      KK      = K
      FOSCI   = 1.0_r8
      FOSNI   = 1.0_r8
      FOSPI   = 1.0_r8
    ELSE
      IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
        IF(K.LE.micpar%NumOfLitrCmplxs)THEN
          OSCM(K)=AMIN1(DCKR,1._r8)*CORGCX(K)*VLSoilMicPMass_vr(L,NY,NX)
        ELSE
          OSCM(K)=AMIN1(FCX,1._r8)*CORGCX(K)*VLSoilMicPMass_vr(L,NY,NX)*DCKM/(CORGCX(k_humus)+DCKM)
        ENDIF
      ELSE
        IF(K.LE.micpar%NumOfLitrCmplxs)THEN
          OSCM(K)=AMIN1(DCKR,1._r8)*CORGCX(K)*VGeomLayer_vr(L,NY,NX)
        ELSE
          OSCM(K)=AMIN1(FCX,1._r8)*CORGCX(K)*VGeomLayer_vr(L,NY,NX)*DCKM/(CORGCX(k_humus)+DCKM)
        ENDIF
      ENDIF
      X=1.0_r8
      KK=micpar%k_humus
      IF(TOSCI.GT.ZEROS(NY,NX))THEN
        FOSCI=AMIN1(1.0_r8,OSCI(KK)/TOSCI)
        FOSNI=AMIN1(1.0_r8,OSCI(KK)*CNOSCT(KK)/TOSNI)
        FOSPI=AMIN1(1.0_r8,OSCI(KK)*CPOSCT(KK)/TOSPI)
!        write(*,*)K,L,OSCI(KK),CNOSCT(KK),CPOSCT(KK),TOSCI,TOSNI,TOSPI
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
    mBiomeAutor_vr(1:NumPlantChemElms,1:NumLiveAutoBioms,L,NY,NX)=0._r8

    D8990: DO N=1,NumMicbFunGrupsPerCmplx
      tglds=JGnfo(N)-JGnio(N)+1._r8
      D8991: DO M=1,nlbiomcp
        OME1(ielmc)=AZMAX1(OSCM(K)*OMCI(M,K)*OMCF(N)*FOSCI)
        OME1(ielmn)=AZMAX1(OME1(ielmc)*rNCOMCa(M,N,K)*FOSNI)
        OME1(ielmp)=AZMAX1(OME1(ielmc)*rPCOMCa(M,N,K)*FOSPI)
!        write(*,*)L,M,OME1(ielmc),rNCOMCa(M,N,K),rPCOMCa(M,N,K)
        do NGL=JGnio(N),JGnfo(N)
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            mBiomeHeter_vr(NE,MID,K,L,NY,NX)=OME1(NE)/tglds
          ENDDO
        ENDDO
        OSCX(KK)=OSCX(KK)+OME1(ielmc)
        OSNX(KK)=OSNX(KK)+OME1(ielmn)
        OSPX(KK)=OSPX(KK)+OME1(ielmp)
        D8992: DO NN=1,NumMicbFunGrupsPerCmplx
          tglds=JGnfA(N)-JGniA(N)+1._r8
          do NGL=JGniA(N),JGnfA(N)
            MID=micpar%get_micb_id(M,NGL)
            DO NE=1,NumPlantChemElms
              mBiomeAutor_vr(NE,MID,L,NY,NX)=mBiomeAutor_vr(NE,MID,L,NY,NX)+OME1(NE)*OMCA(NN)/tglds
            ENDDO
          ENDDO
          OSCX(KK)=OSCX(KK)+OME1(ielmc)*OMCA(NN)
          OSNX(KK)=OSNX(KK)+OME1(ielmn)*OMCA(NN)
          OSPX(KK)=OSPX(KK)+OME1(ielmp)*OMCA(NN)
!          write(*,*)L,K,OME1(ielmc)*OMCA(NN),OME1(ielmn)*OMCA(NN),OME1(ielmp)*OMCA(NN)
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
      OMBioResdu_vr(ielmc,M,K,L,NY,NX)=X*AZMAX1(OSCM(K)*ORCI(M,K)*FOSCI)
      OMBioResdu_vr(ielmn,M,K,L,NY,NX)=AZMAX1(OMBioResdu_vr(ielmc,M,K,L,NY,NX)*rNCOMCa(M,1,K)*FOSNI)
      OMBioResdu_vr(ielmp,M,K,L,NY,NX)=AZMAX1(OMBioResdu_vr(ielmc,M,K,L,NY,NX)*rPCOMCa(M,1,K)*FOSPI)
      OSCX(KK)=OSCX(KK)+OMBioResdu_vr(ielmc,M,K,L,NY,NX)
      OSNX(KK)=OSNX(KK)+OMBioResdu_vr(ielmn,M,K,L,NY,NX)
      OSPX(KK)=OSPX(KK)+OMBioResdu_vr(ielmp,M,K,L,NY,NX)
    ENDDO D8985
!
!     DOC, DON AND DOP
!
!     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores (g)
!     OQCH,OQNH,OQPH,OQAH=DOC,DON,DOP,acetate in macropores (g)
!
    DOM_vr(idom_doc,K,L,NY,NX)               = X*AZMAX1(OSCM(K)*OQCK(K)*FOSCI)
    DOM_vr(idom_don,K,L,NY,NX)               = AZMAX1(DOM_vr(idom_doc,K,L,NY,NX)*CNOSCT(KK)*FOSNI)
    DOM_vr(idom_dop,K,L,NY,NX)               = AZMAX1(DOM_vr(idom_doc,K,L,NY,NX)*CPOSCT(KK)*FOSPI)
    DOM_vr(idom_acetate,K,L,NY,NX)           = 0.0_r8
    DOM_MacP_vr(idom_beg:idom_end,K,L,NY,NX) = 0.0_r8
    OSCX(KK)                                 = OSCX(KK)+DOM_vr(idom_doc,K,L,NY,NX)
    OSNX(KK)                                 = OSNX(KK)+DOM_vr(idom_don,K,L,NY,NX)
    OSPX(KK)                                 = OSPX(KK)+DOM_vr(idom_dop,K,L,NY,NX)
!
!     ADSORBED C, N AND P
!
!     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!
    SorbedOM_vr(ielmc,K,L,NY,NX)=X*AZMAX1(OSCM(K)*OHCK(K)*FOSCI)
    SorbedOM_vr(ielmn,K,L,NY,NX)=AZMAX1(SorbedOM_vr(ielmc,K,L,NY,NX)*CNOSCT(KK)*FOSNI)
    SorbedOM_vr(ielmp,K,L,NY,NX)=AZMAX1(SorbedOM_vr(ielmc,K,L,NY,NX)*CPOSCT(KK)*FOSPI)
    SorbedOM_vr(idom_acetate,K,L,NY,NX)=0.0_r8
    OSCX(KK)=OSCX(KK)+SorbedOM_vr(ielmc,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)
    OSNX(KK)=OSNX(KK)+SorbedOM_vr(ielmn,K,L,NY,NX)
    OSPX(KK)=OSPX(KK)+SorbedOM_vr(ielmp,K,L,NY,NX)
!
!     HUMUS C, N AND P
!
!     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP

    D8980: DO M=1,jsken
      SolidOM_vr(ielmc,M,K,L,NY,NX)=AZMAX1(CFOSC(M,K,L,NY,NX)*(OSCI(K)-OSCX(K)))
      IF(CNOSCT(K).GT.ZERO)THEN
        SolidOM_vr(ielmn,M,K,L,NY,NX)=AZMAX1(CFOSC(M,K,L,NY,NX)*CNOSC(M,K,L,NY,NX) &
          /CNOSCT(K)*(OSNI(K)-OSNX(K)))
      ELSE
        SolidOM_vr(ielmn,M,K,L,NY,NX)=0.0_r8
      ENDIF
      IF(CPOSCT(K).GT.ZERO)THEN
        SolidOM_vr(ielmp,M,K,L,NY,NX)=AZMAX1(CFOSC(M,K,L,NY,NX)*CPOSC(M,K,L,NY,NX) &
          /CPOSCT(K)*(OSPI(K)-OSPX(K)))
      ELSE
        SolidOM_vr(ielmp,M,K,L,NY,NX)=0.0_r8
      ENDIF

      IF(K.EQ.micpar%k_woody_litr)THEN
        SolidOMAct_vr(M,K,L,NY,NX)=SolidOM_vr(ielmc,M,K,L,NY,NX)*OMCI(1,K)
      ELSE
        SolidOMAct_vr(M,K,L,NY,NX)=SolidOM_vr(ielmc,M,K,L,NY,NX)
      ENDIF
    ENDDO D8980
  ENDDO D8995
!
!     ADD ALL LITTER,POC,HUMUS COMPONENTS TO GET TOTAL SOC
!

  RO2DmndHetert(:,:,L,NY,NX)             = 0.0_r8
  RNO3ReduxDmndSoilHeter_vr(:,:,L,NY,NX) = 0.0_r8
  RNO2DmndReduxSoilHeter_vr(:,:,L,NY,NX) = 0.0_r8
  RN2ODmndReduxHeter_vr(:,:,L,NY,NX)     = 0.0_r8
  RNH4DmndSoilHeter_vr(:,:,L,NY,NX)      = 0.0_r8
  RNH4DmndBandHeter_vr(:,:,L,NY,NX)      = 0.0_r8
  RNO3DmndSoilHeter_vr(:,:,L,NY,NX)      = 0.0_r8
  RH2PO4DmndSoilHeter_vr(:,:,L,NY,NX)    = 0.0_r8
  RNO3DmndBandHeter_vr(:,:,L,NY,NX)      = 0.0_r8
  RH2PO4DmndBandHeter_vr(:,:,L,NY,NX)    = 0.0_r8
  RH1PO4DmndSoilHeter_vr(:,:,L,NY,NX)    = 0.0_r8
  RH1PO4DmndBandHeter_vr(:,:,L,NY,NX)    = 0.0_r8
  IF(L.EQ.0)THEN
    RNH4DmndLitrHeter_col(:,:,NY,NX)   = 0.0_r8
    RNO3DmndLitrHeter_col(:,:,NY,NX)   = 0.0_r8
    RH2PO4DmndLitrHeter_col(:,:,NY,NX) = 0.0_r8
  ENDIF

  RO2DmndAutort_vr(:,L,NY,NX)       = 0.0_r8
  RNH3OxidAutor(:,L,NY,NX)          = 0.0_r8
  RNO2OxidAutor(:,L,NY,NX)          = 0.0_r8
  RN2ODmndReduxAutor_vr(:,L,NY,NX)  = 0.0_r8
  RNH4UptkSoilAutor_vr(:,L,NY,NX)   = 0.0_r8
  RNH4UptkBandAutor_vr(:,L,NY,NX)   = 0.0_r8
  RNO3UptkSoilAutor_vr(:,L,NY,NX)   = 0.0_r8
  RNO3UptkBandAutor_vr(:,L,NY,NX)   = 0.0_r8
  RH2PO4UptkSoilAutor_vr(:,L,NY,NX) = 0.0_r8
  RH2PO4UptkBandAutor_vr(:,L,NY,NX) = 0.0_r8
  RH1PO4UptkSoilAutor_vr(:,L,NY,NX) = 0.0_r8
  RH1PO4UptkBandAutor_vr(:,L,NY,NX) = 0.0_r8

  IF(L.EQ.0)THEN
    RNH4UptkLitrAutor_col(:,NY,NX)   = 0.0_r8
    RNO3UptkLitrAutor_col(:,NY,NX)   = 0.0_r8
    RH2PO4UptkLitrAutor_col(:,NY,NX) = 0.0_r8
  ENDIF
  
  call sumORGMLayL(L,NY,NX,ORGM)

  SoilOrgM_vr(1:NumPlantChemElms,L,NY,NX)=ORGM(1:NumPlantChemElms)
  ORGCX_vr(L,NY,NX)=SoilOrgM_vr(ielmc,L,NY,NX)
    
  call sumLitrOMLayL(L,NY,NX,litrOM)

  OMLitrC_vr(L,NY,NX)=litrOM(ielmc)

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
    CFOSC(1:jsken,k_woody_litr,L,NY,NX) = real((/0.00,0.00,0.20,0.80/),r8)
    CFOSC(1:jsken,k_fine_litr,L,NY,NX)  = real((/0.02,0.06,0.34,0.58/),r8)
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
  associate(                           &
    k_POM       => micpar%k_POM,       &
    k_humus     => micpar%k_humus,     &
    k_fine_litr => micpar%k_fine_litr, &
    k_manure    => micpar%k_manure,    &
    iprotein    => micpar%iprotein,    &
    icarbhyro   => micpar%icarbhyro,   &
    icellulos   => micpar%icellulos,   &
    ilignin     => micpar%ilignin,     &
    jsken       => micpar%jsken        &
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
      IF(CSoilOrgM_vr(ielmc,L,NY,NX).LE.FORGC .OR. SoiDepthMidLay_vr(L,NY,NX).LE.ExtWaterTablet0(NY,NX) &
        +CumDepz2LayerBot_vr(NU(NY,NX),NY,NX)-LandScape1stSoiLayDepth)THEN
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
        FC0=FCY*EXP(-5.0_r8*(AMIN1(CORGNX(k_humus), 10.0_r8*CORGPX(k_humus))/CORGCX(k_humus)))
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
    CFOSC(iprotein,k_humus,L,NY,NX)  = FC1
    CFOSC(icarbhyro,k_humus,L,NY,NX) = 1.0_r8-FC1
    CFOSC(icellulos,k_humus,L,NY,NX) = 0.00_r8
    CFOSC(ilignin,k_humus,L,NY,NX)   = 0.00_r8
!
!   MICROBIAL DETRITUS ALLOCATED TO HUMUS MAINTAINS
!   HUMUS PARTITIONING TO COMPONENTS
!
!   ElmAllocmatMicrblitr2POM_vr=fraction of microbial litter allocated to humus components
!
    ElmAllocmatMicrblitr2POM_vr(1,L,NY,NX)=3.0_r8*FC1/(2.0_r8*FC1+1.0_r8)
    ElmAllocmatMicrblitr2POM_vr(2,L,NY,NX)=1.0_r8-ElmAllocmatMicrblitr2POM_vr(1,L,NY,NX)
  ENDIF

  IF(L.GT.0)THEN
    IF(SoiBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
      CORGCM=AMIN1(orgcden,(CORGCX(k_fine_litr)+CORGCX(k_manure)+CORGCX(k_POM)+CORGCX(k_humus)))/0.55_r8
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
  associate(                                   &
    NumOfLitrCmplxs => micpar%NumOfLitrCmplxs, &
    k_POM           => micpar%k_POM,           &
    k_humus         => micpar%k_humus,         &
    CNRH            => micpar%CNRH,            &
    CPRH            => micpar%CPRH             &
  )
  IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
    scal                      = AREA(3,L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX)
    CORGCX(1:NumOfLitrCmplxs) = RSC(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGNX(1:NumOfLitrCmplxs) = RSN(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGPX(1:NumOfLitrCmplxs) = RSP(1:NumOfLitrCmplxs,L,NY,NX)*scal
  ELSE
    scal                      = AREA(3,L,NY,NX)/VGeomLayer_vr(L,NY,NX)
    CORGCX(1:NumOfLitrCmplxs) = RSC(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGNX(1:NumOfLitrCmplxs) = RSN(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGPX(1:NumOfLitrCmplxs) = RSP(1:NumOfLitrCmplxs,L,NY,NX)*scal
  ENDIF
    !
    !     ALLOCATE SOC TO POC(3) AND HUMUS(4)
    !
  IF(L.GT.0)THEN
    CORGCZ=CSoilOrgM_vr(ielmc,L,NY,NX)
    CORGRZ=COMLitrC_vr(L,NY,NX)
    CORGNZ=CSoilOrgM_vr(ielmn,L,NY,NX)
    CORGPZ=CSoilOrgM_vr(ielmp,L,NY,NX)
    IF(CORGCZ.GT.ZERO)THEN
      CORGCX(k_POM)   = CORGRZ
      CORGCX(k_humus) = AZMAX1(CORGCZ-CORGCX(k_POM))
      CORGNX(k_POM)   = AMIN1(CNRH(k_POM)*CORGCX(k_POM),CORGNZ)
      CORGNX(k_humus) = AZMAX1(CORGNZ-CORGNX(k_POM))
      CORGPX(k_POM)   = AMIN1(CPRH(k_POM)*CORGCX(k_POM),CORGPZ)
      CORGPX(k_humus) = AZMAX1(CORGPZ-CORGPX(k_POM))

    ELSE
      CORGCX(k_POM)   = 0.0_r8
      CORGCX(k_humus) = 0.0_r8
      CORGNX(k_POM)   = 0.0_r8
      CORGNX(k_humus) = 0.0_r8
      CORGPX(k_POM)   = 0.0_r8
      CORGPX(k_humus) = 0.0_r8
    ENDIF
  ELSE
    CORGCX(k_POM)   = 0.0_r8
    CORGCX(k_humus) = 0.0_r8
    CORGNX(k_POM)   = 0.0_r8
    CORGNX(k_humus) = 0.0_r8
    CORGPX(k_POM)   = 0.0_r8
    CORGPX(k_humus) = 0.0_r8
  ENDIF
!  write(*,*)L,CORGCX
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
