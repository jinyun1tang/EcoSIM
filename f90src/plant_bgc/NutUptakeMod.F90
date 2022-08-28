module NutUptakeMod

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use StomatesMod   , only : stomates
  use minimathmod  , only : safe_adb,vapsat,test_aneb
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use PlantAPIData
  use RootGasMod
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = __FILE__
  public :: NutO2Uptake
  contains

!------------------------------------------------------------------------

  subroutine NutO2Uptake(NZ,FDMP,OSTRN,OSTRD,PATH,RRADL,FPQ,FPP,FRTDPX,RTARR)
  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in):: FDMP
  real(r8), intent(in) :: PATH(2,JZ1),RRADL(2,JZ1),FPQ(2,JZ1,JP1),FPP(2,JZ1,JP1)
  real(r8), intent(in) :: FRTDPX(JZ1,JP1),RTARR(2,JZ1)
  real(r8), intent(inout) :: OSTRN,OSTRD

  call CanopyNH3Flux(NZ,FDMP)
!
!     ROOT(N=1) AD MYCORRHIZAL(N=2) O2 AND NUTRIENT UPTAKE
!
  call RootMycoO2NutrientUptake(NZ,OSTRN,OSTRD,PATH,RRADL,&
    FPQ,FPP,FRTDPX,RTARR)
  end subroutine NutO2Uptake
!------------------------------------------------------------------------

  subroutine CanopyNH3Flux(NZ,FDMP)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in):: FDMP
  real(r8) :: FNH3P
  real(r8) :: CNH3P
  real(r8) :: SNH3P
  real(r8) :: ZPOOLB
  integer :: NB

  associate(                            &
    TCCs1     =>  plt_ew%TCCs1        , &
    NUs1      =>  plt_site%NUs1       , &
    AREA3s1   =>  plt_site%AREA3s1    , &
    CNH3Es1   =>  plt_site%CNH3Es1    , &
    RNH3Bs1   =>  plt_rbgc%RNH3Bs1    , &
    ZEROPs1   =>  plt_biom%ZEROPs1    , &
    WTLSBs1   =>  plt_biom%WTLSBs1    , &
    ZPOOLs1   =>  plt_biom%ZPOOLs1    , &
    CCPOLBs1  =>  plt_biom%CCPOLBs1   , &
    CZPOLBs1  =>  plt_biom%CZPOLBs1   , &
    CPPOLBs1  =>  plt_biom%CPPOLBs1   , &
    RCs1      =>  plt_photo%RCs1      , &
    RAs1      =>  plt_photo%RAs1      , &
    ARLFBs1   =>  plt_morph%ARLFBs1   , &
    NBRs1     =>  plt_morph%NBRs1     , &
    FRADPs1   =>  plt_rad%FRADPs1     , &
    ARLFPs1   =>  plt_morph%ARLFPs1     &
  )
  !
  !     NH3 EXCHANGE BETWEEN CANOPY AND ATMOSPHERE FROM NH3
  !     CONCENTRATION DIFFERENCES 'CNH3E' (ATMOSPHERE FROM 'READS') AND
  !     'CNH3P' (CANOPY), AND FROM STOMATAL + BOUNDARY LAYER RESISTANCE
  !
  !     SNH3P,SNH3X=NH3 solubility at TCC, 25 oC
  !     TCC=canopy temperature (oC)
  !     FDMP,FNH3P=canopy dry matter content,NH3 concentration
  !     ARLFB,ARLFP=branch,canopy leaf area
  !     CNH3P,CNH3E=gaseous NH3 concentration in branch,atmosphere
  !     CZPOLB,ZPOOLB=nonstplt_rbgc%RUCtural N concentration,content in branch
  !     RNH3B=NH3 flux between atmosphere and branch
  !     RA,RC=canopy boundary layer,stomatal resistance
  !     FRADP=fraction of radiation received by each PFT canopy
  !
  SNH3P=SNH3X*EXP(0.513-0.0171*TCCs1(NZ))
  FNH3P=1.0E-04*FDMP
  DO 105 NB=1,NBRs1(NZ)
    IF(WTLSBs1(NB,NZ).GT.ZEROPs1(NZ) &
      .AND.ARLFBs1(NB,NZ).GT.ZEROPs1(NZ) &
      .AND.ARLFPs1(NZ).GT.ZEROPs1(NZ))THEN
      CNH3P=AMAX1(0.0,FNH3P*CZPOLBs1(NB,NZ)/SNH3P)
      ZPOOLB=AMAX1(0.0,ZPOOLs1(NB,NZ))
      RNH3Bs1(NB,NZ)=AMIN1(0.1*ZPOOLB,AMAX1((CNH3Es1-CNH3P)/(RAs1(NZ)+RCs1(NZ)) &
        *FRADPs1(NZ)*AREA3s1(NUs1)*ARLFBs1(NB,NZ)/ARLFPs1(NZ),-0.1*ZPOOLB))
      ELSE
      RNH3Bs1(NB,NZ)=0.0_r8
      ENDIF

105   CONTINUE
  end associate
  end subroutine CanopyNH3Flux

!------------------------------------------------------------------------------------------

  subroutine RootMycoO2NutrientUptake(NZ,OSTRN,OSTRD,PATH,RRADL,&
    FPQ,FPP,FRTDPX,RTARR)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: PATH(2,JZ1),RRADL(2,JZ1),FPQ(2,JZ1,JP1),FPP(2,JZ1,JP1)
  real(r8), intent(in) :: FRTDPX(JZ1,JP1),RTARR(2,JZ1)
  real(r8), intent(inout) :: OSTRN,OSTRD
  real(r8) :: TFOXYX
  real(r8) :: FCUP,FZUP,FPUP,FWSRT,UPWTRP,UPWTRH,FOXYX,RUPOXT
  integer :: N,L
!     begin_execution
  associate(                             &
    THETWs1 =>  plt_soilchem%THETWs1   , &
    VOLXs1  =>  plt_soilchem%VOLXs1    , &
    ZEROS2s1=>  plt_site%ZEROS2s1      , &
    NUs1    =>  plt_site%NUs1          , &
    ZEROs1  =>  plt_site%ZEROs1        , &
    ZEROPs1 =>  plt_biom%ZEROPs1       , &
    ROXYPs1 =>  plt_rbgc%ROXYPs1       , &
    WFRs1   =>  plt_rbgc%WFRs1         , &
    RCO2Ms1 =>  plt_rbgc%RCO2Ms1       , &
    RTLGPs1 =>  plt_morph%RTLGPs1      , &
    MYs1    =>  plt_morph%MYs1         , &
    RTDNPs1 =>  plt_morph%RTDNPs1      , &
    RTVLWs1 =>  plt_morph%RTVLWs1      , &
    NIs1    =>  plt_morph%NIs1           &
  )


  call ZeroUptake(NZ)

  DO 955 N=1,MYs1(NZ)
    DO 950 L=NUs1,NIs1(NZ)
      IF(VOLXs1(L).GT.ZEROS2s1 &
        .AND.RTDNPs1(N,L,NZ).GT.ZEROs1 &
        .AND.RTVLWs1(N,L,NZ).GT.ZEROPs1(NZ) &
        .AND.THETWs1(L).GT.ZEROs1)THEN
        TFOXYX=0.0_r8
        call GetUptakeCapcity(N,L,NZ,FPQ,FPP,FCUP,FZUP,FPUP,&
          FWSRT,UPWTRP,UPWTRH,FOXYX)

        TFOXYX=TFOXYX+FOXYX
!
!     ROOT O2 DEMAND CALCULATED FROM O2 NON-LIMITED RESPIRATION RATE
!
!     ROXYP=O2 demand
!     RCO2M=respiration unlimited by O2
!     RTVLW=root or myco aqueous volume
!     FOXYX=fraction of total O2 demand from previous hour
!
        ROXYPs1(N,L,NZ)=2.667*RCO2Ms1(N,L,NZ)

        call RootSoilGasExchange(N,L,NZ,RRADL,FPQ,FRTDPX,RTARR,UPWTRH,&
          FOXYX,RUPOXT)

        OSTRD=OSTRD+ROXYPs1(N,L,NZ)
        OSTRN=OSTRN+RUPOXT

        call RootExudates(N,L,NZ)
!
!     NUTRIENT UPTAKE
!
!     WFR=constraint by O2 consumption on all biological processes
!     FCUP=limitation to active uptake respiration from CPOOLR
!     FWSRT=protein concentration relative to 5%
!     RTLGP=root,myco length per plant
!
        IF(WFRs1(N,L,NZ).GT.ZEROs1.AND.FCUP.GT.ZEROs1.AND.FWSRT.GT.ZEROs1 &
          .AND.RTLGPs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
!
!     FZUP=limitn to active uptake respiration from CZPOLR
!
          call UptakeMineralNitrogen(N,L,NZ,PATH,RRADL,FPQ,FPP,&
            RTARR,FCUP,FZUP,FWSRT,UPWTRP)
!
!     FPUP=limitn to active uptake respiration from CPPOLR
!
          call UptakeMineralPhosporhus(N,L,NZ,PATH,RRADL,FPQ,FPP,&
            RTARR,FCUP,FPUP,FWSRT,UPWTRP)
        ENDIF
      ELSE

      ENDIF

      call SumupNutrientUptake(N,L,NZ)

950 CONTINUE
955 CONTINUE
  end associate
  end subroutine RootMycoO2NutrientUptake
!------------------------------------------------------------------------
  subroutine ZeroUptake(NZ)

  implicit none
  integer, intent(in) :: NZ

  integer :: K, L1,L2,NN
  !     begin_execution

  L1=plt_site%NUs1;L2=plt_morph%NIs1(NZ);NN=plt_morph%MYs1(NZ)

  plt_rbgc%RCOFLAs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%ROXFLAs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RCHFLAs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RN2FLAs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RNHFLAs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RCODFAs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%ROXDFAs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RCHDFAs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RN2DFAs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RNHDFAs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RCO2Ss1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPOXSs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPCHSs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPN2Ss1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPN3Ss1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RCO2Ps1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPOXPs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RDFOMCs1(1:NN,0:jcplx11,L1:L2,NZ)=0.0_r8
  plt_rbgc%RDFOMNs1(1:NN,0:jcplx11,L1:L2,NZ)=0.0_r8
  plt_rbgc%RDFOMPs1(1:NN,0:jcplx11,L1:L2,NZ)=0.0_r8
  plt_rbgc%WFRs1(1:NN,L1:L2,NZ)=1.0
  plt_rbgc%RUNNHPs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPNH4s1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUONH4s1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCNH4s1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUNNBPs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPNHBs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUONHBs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCNHBs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUNNOPs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPNO3s1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUONO3s1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCNO3s1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUNNXPs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPNOBs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUONOBs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCNOBs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPP2Ps1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPH2Ps1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUOH2Ps1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCH2Ps1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPP2Bs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPH2Bs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUOH2Bs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCH2Bs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPP1Ps1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPH1Ps1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUOH1Ps1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCH1Ps1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPP1Bs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPH1Bs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUOH1Bs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCH1Bs1(1:NN,L1:L2,NZ)=0.0_r8
  plt_bgcr%RUPNFs1(L1:L2,NZ)=0.0_r8
  end subroutine ZeroUptake

!------------------------------------------------------------------------0

  subroutine UptakeMineralPhosporhus(N,L,NZ,PATH,RRADL,FPQ,FPP,&
    RTARR,FCUP,FPUP,FWSRT,UPWTRP)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  real(r8), intent(in):: PATH(2,JZ1),RRADL(2,JZ1),FPQ(2,JZ1,JP1),FPP(2,JZ1,JP1)
  real(r8), intent(in):: RTARR(2,JZ1),FCUP,FPUP,FWSRT,UPWTRP
  real(r8) :: TFPO4X,TFP14X,TFPOBX,TFP1BX
  real(r8) :: DIFFL
  real(r8) :: FP14X,FP1BX
  real(r8) :: FPO4X,FPOBX
  real(r8) :: POSGX
  real(r8) :: PATHL
  associate(                               &
    TORTs1    =>   plt_site%TORTs1       , &
    ZEROSs1   =>   plt_site%ZEROSs1      , &
    ZERO2s1   =>   plt_site%ZERO2s1      , &
    RUPP1Bs1  =>   plt_rbgc%RUPP1Bs1     , &
    RUPP1Ps1  =>   plt_rbgc%RUPP1Ps1     , &
    RUPP2Ps1  =>   plt_rbgc%RUPP2Ps1     , &
    RUPP2Bs1  =>   plt_rbgc%RUPP2Bs1     , &
    RP1BYs1   =>   plt_bgcr%RP1BYs1      , &
    RPOBYs1   =>   plt_bgcr%RPOBYs1      , &
    RP14Ys1   =>   plt_bgcr%RP14Ys1      , &
    RPO4Ys1   =>   plt_bgcr%RPO4Ys1      , &
    POSGLs1   =>   plt_soilchem%POSGLs1    &
  )
  TFPO4X=0.0_r8
  TFPOBX=0.0_r8
  TFP14X=0.0_r8
  TFP1BX=0.0_r8

!     begin_execution
  IF(RPO4Ys1(L).GT.ZEROSs1)THEN
    FPO4X=AMAX1(FPP(N,L,NZ),RUPP2Ps1(N,L,NZ)/RPO4Ys1(L))
  ELSE
    FPO4X=FPQ(N,L,NZ)
  ENDIF
  IF(RPOBYs1(L).GT.ZEROSs1)THEN
    FPOBX=AMAX1(FPP(N,L,NZ),RUPP2Bs1(N,L,NZ)/RPOBYs1(L))
  ELSE
    FPOBX=FPQ(N,L,NZ)
  ENDIF
  IF(RP14Ys1(L).GT.ZEROSs1)THEN
    FP14X=AMAX1(FPP(N,L,NZ),RUPP1Ps1(N,L,NZ)/RP14Ys1(L))
  ELSE
    FP14X=FPQ(N,L,NZ)
  ENDIF
  IF(RP1BYs1(L).GT.ZEROSs1)THEN
    FP1BX=AMAX1(FPP(N,L,NZ),RUPP1Bs1(N,L,NZ)/RP1BYs1(L))
  ELSE
    FP1BX=FPQ(N,L,NZ)
  ENDIF

  TFPO4X=TFPO4X+FPO4X
  TFP14X=TFP14X+FP14X

  TFPOBX=TFPOBX+FPOBX
  TFP1BX=TFP1BX+FP1BX

  IF(FPUP.GT.ZERO2s1)THEN
    !
    !     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF H2PO4,HPO4
    !     FROM SOIL TO ROOT
    !
    !     POSGL=PO4 diffusivity
    !     TORT=soil tortuosity
    !     PATH=path length of water and nutrient uptake
    !     RRADL=root radius
    !     DIFFL=PO4 diffusion per plant
    !
    POSGX=POSGLs1(L)*TORTs1(NPH,L)
    PATHL=AMIN1(PATH(N,L),RRADL(N,L)+SQRT(2.0*POSGX))
    DIFFL=POSGX*safe_adb(RTARR(N,L),LOG(PATHL/RRADL(N,L)))

    call UptakeH2PO4(N,L,NZ,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FWSRT,UPWTRP)

    call UptakeHPO4(N,L,NZ,DIFFL,FP14X,FP1BX,FCUP,FPUP,FWSRT,UPWTRP)

  ENDIF
  end associate
  end subroutine UptakeMineralPhosporhus
!------------------------------------------------------------------------

  subroutine UptakeNO3(N,L,NZ,FNO3X,FNOBX,PATH,RRADL,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  real(r8), intent(in):: FNO3X,FNOBX,PATH(2,JZ1),RRADL(2,JZ1),RTARR(2,JZ1)
  real(r8), intent(in):: FCUP,FZUP,FWSRT,UPWTRP
  real(r8) :: B,C
  real(r8) :: BP,CP
  real(r8) :: DIFFL
  real(r8) :: DIFNO3,DIFNOB
  real(r8) :: PATHL
  real(r8) :: RMFNO3,RTKNO3,RTKNOP,RMFNOB,RTKNOB,RTKNPB
  real(r8) :: UPMX,UPMXP,X,Y
  real(r8) :: ZOSGX,ZNO3M,ZNO3X,ZNOBM,ZNOBX
! begin_execution
  associate(                              &
    PPs1      =>  plt_site%PPs1         , &
    TORTs1    =>  plt_site%TORTs1       , &
    ZEROs1    =>  plt_site%ZEROs1       , &
    RTARPs1   =>  plt_morph%RTARPs1     , &
    TFN4s1    =>  plt_pheno%TFN4s1      , &
    WFRs1     =>  plt_rbgc%WFRs1        , &
    UPMXZOs1  =>  plt_rbgc%UPMXZOs1     , &
    UPMNZOs1  =>  plt_rbgc%UPMNZOs1     , &
    UPKMZOs1  =>  plt_rbgc%UPKMZOs1     , &
    RUCNOBs1  =>  plt_rbgc%RUCNOBs1     , &
    RUNNOPs1  =>  plt_rbgc%RUNNOPs1     , &
    RUCNO3s1  =>  plt_rbgc%RUCNO3s1     , &
    RUPNO3s1  =>  plt_rbgc%RUPNO3s1     , &
    RUNNXPs1  =>  plt_rbgc%RUNNXPs1     , &
    RUONOBs1  =>  plt_rbgc%RUONOBs1     , &
    RUPNOBs1  =>  plt_rbgc%RUPNOBs1     , &
    RUONO3s1  =>  plt_rbgc%RUONO3s1     , &
    VLNOBs1   =>  plt_soilchem%VLNOBs1  , &
    CNO3Bs1   =>  plt_soilchem%CNO3Bs1  , &
    ZNO3Ss1   =>  plt_soilchem%ZNO3Ss1  , &
    ZOSGLs1   =>  plt_soilchem%ZOSGLs1  , &
    CNO3Ss1   =>  plt_soilchem%CNO3Ss1  , &
    VOLWs1    =>  plt_soilchem%VOLWs1   , &
    VLNO3s1   =>  plt_soilchem%VLNO3s1  , &
    ZNO3Bs1   =>  plt_soilchem%ZNO3Bs1    &
  )
!
! PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NO3
! FROM SOIL TO ROOT
!
! ZOSGL=NO3 diffusivity
! TORT=soil tortuosity
! RRADL=root radius
! PATH=path length of water and nutrient uptake
! DIFFL=NO3 diffusion per plant
!
  ZOSGX=ZOSGLs1(L)*TORTs1(NPH,L)
  PATHL=AMIN1(PATH(N,L),RRADL(N,L)+SQRT(2.0*ZOSGX))
  DIFFL=ZOSGX*safe_adb(RTARR(N,L),LOG(PATHL/RRADL(N,L)))
  !
  ! NO3 UPTAKE IN NON-BAND SOIL ZONE
  !
  !  VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
  !     CNO3S=NO3 concentration in non-band
  !     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFNO3=soil-root convective NO3 flux per plant in non-band
  !     DIFNO3=soil-root NO3 diffusion per plant in non-band
  !
  IF(VLNO3s1(L).GT.ZEROs1.AND.CNO3Ss1(L) &
    .GT.UPMNZOs1(N,NZ))THEN
    RMFNO3=UPWTRP*CNO3Ss1(L)*VLNO3s1(L)
    DIFNO3=DIFFL*VLNO3s1(L)
    !
    !     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=max NO3 uptake in non-band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     TFN4=temperature function for root growth
    !     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=UPMXZOs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLNO3s1(L)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFNO3=soil-root convective N03 flux per plant in non-band
    !     DIFNO3=soil-root N03 diffusion per plant in non-band
    !     CNO3S=NO3 concentration in non-band
    !     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
    !     RTKNO3,RTKNOP=NO3 uptake per plant in non-band lmtd,unlmtd by O2
    !     ZNO3M,ZNO3X=minimum,maximum NO3 available for uptake in non-band
    !     FNO3X=fraction of total NH4 uptake in non-band by root,myco populn
    !     RUNNOP,RUPNO3=NO3 uptake in non-band unlimited,limited by NO3
    !     RUONO3=NO3 uptake in non-band unlimited by O2
    !     RUCNO3=NO3 uptake in non-band unlimited by nonstructural C
!
    X=(DIFNO3+RMFNO3)*CNO3Ss1(L)
    Y=DIFNO3*UPMNZOs1(N,NZ)
    B=-UPMX-DIFNO3*UPKMZOs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKNO3=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNO3*UPKMZOs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKNOP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNO3M=UPMNZOs1(N,NZ)*VOLWs1(L)*VLNO3s1(L)
    ZNO3X=AMAX1(0.0,FNO3X*(ZNO3Ss1(L)-ZNO3M))
    RUNNOPs1(N,L,NZ)=AMAX1(0.0,RTKNO3*PPs1(NZ))
    RUPNO3s1(N,L,NZ)=AMIN1(ZNO3X,RUNNOPs1(N,L,NZ))
    RUONO3s1(N,L,NZ)=AMIN1(ZNO3X,AMAX1(0.0 &
      ,RTKNOP*PPs1(NZ)))
    RUCNO3s1(N,L,NZ)=RUPNO3s1(N,L,NZ)/FCUP
  ENDIF
  !
  !     NO3 UPTAKE IN BAND SOIL ZONE
  !
  !     VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
  !     CNO3B=NO3 concentration in band
  !     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFNOB=soil-root convective NO3 flux per plant in band
  !     DIFNOB=soil-root NO3 diffusion per plant in band
  !

  IF(VLNOBs1(L).GT.ZEROs1.AND.CNO3Bs1(L) &
    .GT.UPMNZOs1(N,NZ))THEN
    RMFNOB=UPWTRP*CNO3Bs1(L)*VLNOBs1(L)
    DIFNOB=DIFFL*VLNOBs1(L)
    !
    !     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum NO3 uptake in band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     TFN4=temperature function for root growth
    !     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=UPMXZOs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLNOBs1(L)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFNOB=soil-root convective NO3 flux per plant in band
    !     DIFNOB=soil-root NO3 diffusion per plant in band
    !     CNO3B=NH4 concentration in band
    !     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
    !     RTKNOB,RTKNPB=NO3 uptake per plant in band lmtd,unlmtd by O2
    !     ZNOBM,ZNOBX=minimum,maximum NO3 available for uptake in band
    !     FNOBX=fraction of total NO3 uptake in band by root,myco populn
    !     RUNNXP,RUPNOB=NO3 uptake in band unlimited,limited by NH4
    !     RUONOB=NO3 uptake in band unlimited by O2
    !     RUCNOB=NO3 uptake in band unlimited by nonstructural C
    !
    X=(DIFNOB+RMFNOB)*CNO3Bs1(L)
    Y=DIFNOB*UPMNZOs1(N,NZ)
    B=-UPMX-DIFNOB*UPKMZOs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKNOB=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNOB*UPKMZOs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKNPB=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNOBM=UPMNZOs1(N,NZ)*VOLWs1(L)*VLNOBs1(L)
    ZNOBX=AMAX1(0.0,FNOBX*(ZNO3Bs1(L)-ZNOBM))
    RUNNXPs1(N,L,NZ)=AMAX1(0.0,RTKNOB*PPs1(NZ))
    RUPNOBs1(N,L,NZ)=AMIN1(ZNOBX,RUNNXPs1(N,L,NZ))
    RUONOBs1(N,L,NZ)=AMIN1(ZNOBX,AMAX1(0.0,RTKNPB*PPs1(NZ)))
    RUCNOBs1(N,L,NZ)=RUPNOBs1(N,L,NZ)/FCUP
  ENDIF
  end associate
  end subroutine UptakeNO3
!------------------------------------------------------------------------
  subroutine UptakeNH4(N,L,NZ,FNH4X,FNHBX,PATH,RRADL,RTARR,&
    FCUP,FZUP,FWSRT,UPWTRP)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: FNH4X,FNHBX,PATH(2,JZ1),RRADL(2,JZ1)
  real(r8), intent(in) :: RTARR(2,JZ1)
  real(r8), intent(in) :: FCUP,FZUP,FWSRT,UPWTRP
  real(r8) :: B,C
  real(r8) :: BP,CP
  real(r8) :: DIFFL
  real(r8) :: DIFNH4,DIFNHB
  real(r8) :: PATHL
  real(r8) :: RMFNH4,RTKNH4,RTKNHP,RMFNHB,RTKNHB,RTKNBP
  real(r8) :: UPMX,UPMXP,X,Y
  real(r8) :: ZNHBX,ZNSGX,ZNH4M,ZNH4X,ZNHBM
! begin_execution
  associate(                              &
    PPs1      =>  plt_site%PPs1         , &
    ZEROs1    =>  plt_site%ZEROs1       , &
    TORTs1    =>  plt_site%TORTs1       , &
    TFN4s1    =>  plt_pheno%TFN4s1      , &
    WFRs1     =>  plt_rbgc%WFRs1        , &
    UPMNZHs1  =>  plt_rbgc%UPMNZHs1     , &
    UPMXZHs1  =>  plt_rbgc%UPMXZHs1     , &
    UPKMZHs1  =>  plt_rbgc%UPKMZHs1     , &
    RUPNHBs1  =>  plt_rbgc%RUPNHBs1     , &
    RUCNH4s1  =>  plt_rbgc%RUCNH4s1     , &
    RUONHBs1  =>  plt_rbgc%RUONHBs1     , &
    RUONH4s1  =>  plt_rbgc%RUONH4s1     , &
    RUNNHPs1  =>  plt_rbgc%RUNNHPs1     , &
    RUCNHBs1  =>  plt_rbgc%RUCNHBs1     , &
    RUNNBPs1  =>  plt_rbgc%RUNNBPs1     , &
    RUPNH4s1  =>  plt_rbgc%RUPNH4s1     , &
    RTARPs1   =>  plt_morph%RTARPs1     , &
    ZNSGLs1   =>  plt_soilchem%ZNSGLs1  , &
    VOLWs1    =>  plt_soilchem%VOLWs1   , &
    VLNHBs1   =>  plt_soilchem%VLNHBs1  , &
    ZNH4Bs1   =>  plt_soilchem%ZNH4Bs1  , &
    CNH4Bs1   =>  plt_soilchem%CNH4Bs1  , &
    VLNH4s1   =>  plt_soilchem%VLNH4s1  , &
    ZNH4Ss1   =>  plt_soilchem%ZNH4Ss1  , &
    CNH4Ss1   =>  plt_soilchem%CNH4Ss1    &
  )
! ZNSGL=NH4 diffusivity
! TORT=soil tortuosity
! PATH=path length of water and nutrient uptake
! RRADL=root radius
! DIFFL=NH4 diffusion per plant
!
  ZNSGX=ZNSGLs1(L)*TORTs1(NPH,L)
  PATHL=AMIN1(PATH(N,L),RRADL(N,L)+SQRT(2.0*ZNSGX))
  DIFFL=ZNSGX*safe_adb(RTARR(N,L),LOG(PATHL/RRADL(N,L)))
!
! NH4 UPTAKE IN NON-BAND SOIL ZONE
!
! VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
! CNH4S=NH4 concentration in non-band
! UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
! UPWTRP=root water uptake per plant
! RMFNH4=soil-root convective NH4 flux per plant in non-band
! DIFNH4=soil-root NH4 diffusion per plant in non-band
!
  IF(VLNH4s1(L).GT.ZEROs1.AND.CNH4Ss1(L).GT.UPMNZHs1(N,NZ))THEN
    RMFNH4=UPWTRP*CNH4Ss1(L)*VLNH4s1(L)
    DIFNH4=DIFFL*VLNH4s1(L)
!
!   NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
!   AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
!
!   UPMXP,UPMX=max NH4 uptake in non-band unlimited,limited by O2
!   RTARP=root surface area per plant from grosub.f
!   FWSRT=protein concentration relative to 5%
!   TFN4=temperature function for root growth
!   FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
!   WFR=constraint by O2 consumption on all biological processes
!
    UPMXP=UPMXZHs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLNH4s1(L)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
!
!   SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
!   SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
!   WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!   RMFNH4=soil-root convective NH4 flux per plant in non-band
!   DIFNH4=soil-root NH4 diffusion per plant in non-band
!   CNH4S=NH4 concentration in non-band
!   UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
!   RTKNH4,RTKNHP=NH4 uptake per plant in non-band lmtd,unlmtd by O2
!   ZNH4M,ZNH4X=minimum,maximum NH4 available for uptake in non-band
!   FNH4X=fraction of total NH4 uptake in non-band by root,myco populn
!   RUNNHP,RUPNH4=NH4 uptake in non-band unlimited,limited by NH4
!   RUONH4=NH4 uptake in non-band unlimited by O2
!   RUCNH4=NH4 uptake in non-band unlimited by nonstructural C
!
    X=(DIFNH4+RMFNH4)*CNH4Ss1(L)
    Y=DIFNH4*UPMNZHs1(N,NZ)
    B=-UPMX-DIFNH4*UPKMZHs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKNH4=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNH4*UPKMZHs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKNHP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNH4M=UPMNZHs1(N,NZ)*VOLWs1(L)*VLNH4s1(L)
    ZNH4X=AMAX1(0.0,FNH4X*(ZNH4Ss1(L)-ZNH4M))
    RUNNHPs1(N,L,NZ)=AMAX1(0.0,RTKNH4*PPs1(NZ))
    RUPNH4s1(N,L,NZ)=AMIN1(ZNH4X,RUNNHPs1(N,L,NZ))
    RUONH4s1(N,L,NZ)=AMIN1(ZNH4X,AMAX1(0.0 &
      ,RTKNHP*PPs1(NZ)))
    RUCNH4s1(N,L,NZ)=RUPNH4s1(N,L,NZ)/FCUP
  ENDIF
!
! NH4 UPTAKE IN BAND SOIL ZONE
!
! VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
! CNH4B=NH4 concentration in band
! UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
! UPWTRP=root water uptake per plant
! RMFNHB=soil-root convective NH4 flux per plant in band
! DIFNHB=soil-root NH4 diffusion per plant in band
!

  IF(VLNHBs1(L).GT.ZEROs1.AND.CNH4Bs1(L) &
    .GT.UPMNZHs1(N,NZ))THEN
    RMFNHB=UPWTRP*CNH4Bs1(L)*VLNHBs1(L)
    DIFNHB=DIFFL*VLNHBs1(L)
!
!   NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
!   AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
!
!   UPMXP,UPMX=maximum NH4 uptake in band unlimited,limited by O2
!   RTARP=root surface area per plant from grosub.f
!   FWSRT=protein concentration relative to 5%
!   TFN4=temperature function for root growth
!   FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
!   WFR=constraint by O2 consumption on all biological processes
!
    UPMXP=UPMXZHs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLNHBs1(L)*AMIN1(FCUP,FZUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
!
!   SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
!   SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
!   WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!   RMFNHB=soil-root convective NH4 flux per plant in band
!   DIFNHB=soil-root NH4 diffusion per plant in band
!   CNH4B=NH4 concentration in band
!   UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
!   RTKNHB,RTKNBP=NH4 uptake per plant in band lmtd,unlmtd by O2
!   ZNHBM,ZNHBX=minimum,maximum NH4 available for uptake in band
!   FNHBX=fraction of total NH4 uptake in band by root,myco populn
!   RUNNBP,RUPNHB=NH4 uptake in band unlimited,limited by NH4
!   RUONHB=NH4 uptake in band unlimited by O2
!   RUCNHB=NH4 uptake in band unlimited by nonstructural C
!
    X=(DIFNHB+RMFNHB)*CNH4Bs1(L)
    Y=DIFNHB*UPMNZHs1(N,NZ)
    B=-UPMX-DIFNHB*UPKMZHs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKNHB=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFNHB*UPKMZHs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKNBP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    ZNHBM=UPMNZHs1(N,NZ)*VOLWs1(L)*VLNHBs1(L)
    ZNHBX=AMAX1(0.0,FNHBX*(ZNH4Bs1(L)-ZNHBM))
    RUNNBPs1(N,L,NZ)=AMAX1(0.0,RTKNHB*PPs1(NZ))
    RUPNHBs1(N,L,NZ)=AMIN1(ZNHBX,RUNNBPs1(N,L,NZ))
    RUONHBs1(N,L,NZ)=AMIN1(ZNHBX,AMAX1(0.0 &
      ,RTKNBP*PPs1(NZ)))
    RUCNHBs1(N,L,NZ)=RUPNHBs1(N,L,NZ)/FCUP
  ENDIF
  end associate
  end subroutine UptakeNH4

!------------------------------------------------------------------------
  subroutine UptakeHPO4(N,L,NZ,DIFFL,FP14X,FP1BX,FCUP,FPUP,FWSRT,UPWTRP)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: DIFFL
  real(r8), intent(in) :: FP14X,FP1BX
  real(r8), intent(in) :: FCUP,FPUP,FWSRT,UPWTRP
  real(r8) :: B,C
  real(r8) :: BP,CP
  real(r8) :: DIFH1P
  real(r8) :: DIFH1B
  real(r8) :: H1POM,H1POX,H1PXM,H1PXB
  real(r8) :: RMFH1P,RTKH1P,RTKHP1,RTKH1B,RTKHB1
  real(r8) :: RMFH2B
  real(r8) :: UPMX,UPMXP
  real(r8) :: X,Y
  !     begin_execution
  associate(                             &
    PPs1      => plt_site%PPs1         , &
    ZEROs1    => plt_site%ZEROs1       , &
    TFN4s1    => plt_pheno%TFN4s1      , &
    WFRs1     => plt_rbgc%WFRs1        , &
    UPMNPOs1  => plt_rbgc%UPMNPOs1     , &
    UPKMPOs1  => plt_rbgc%UPKMPOs1     , &
    UPMXPOs1  => plt_rbgc%UPMXPOs1     , &
    RUPH1Bs1  => plt_rbgc%RUPH1Bs1     , &
    RUCH1Ps1  => plt_rbgc%RUCH1Ps1     , &
    RUOH1Ps1  => plt_rbgc%RUOH1Ps1     , &
    RUPP1Ps1  => plt_rbgc%RUPP1Ps1     , &
    RUPH1Ps1  => plt_rbgc%RUPH1Ps1     , &
    RUCH1Bs1  => plt_rbgc%RUCH1Bs1     , &
    RUOH1Bs1  => plt_rbgc%RUOH1Bs1     , &
    RUPP1Bs1  => plt_rbgc%RUPP1Bs1     , &
    RTARPs1   => plt_morph%RTARPs1     , &
    VOLWs1    => plt_soilchem%VOLWs1   , &
    VLPO4s1   => plt_soilchem%VLPO4s1  , &
    CH1P4s1   => plt_soilchem%CH1P4s1  , &
    H1PO4s1   => plt_soilchem%H1PO4s1  , &
    VLPOBs1   => plt_soilchem%VLPOBs1  , &
    CH1P4Bs1  => plt_soilchem%CH1P4Bs1 , &
    H1POBs1   => plt_soilchem%H1POBs1    &
  )
  !
  !     HPO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH1P4=HPO4 concentration in non-band
  !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFH1P=soil-root convective HPO4 flux per plant in non-band
  !     DIFH1P=soil-root HPO4 diffusion per plant in non-band
!
    IF(VLPO4s1(L).GT.ZEROs1.AND.CH1P4s1(L) &
      .GT.UPMNPOs1(N,NZ))THEN
      RMFH1P=UPWTRP*CH1P4s1(L)*VLPO4s1(L)
      DIFH1P=DIFFL*VLPO4s1(L)
    !
    !     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=max HPO4 uptake in non-band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     TFN4=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=0.1*UPMXPOs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLPO4s1(L)*AMIN1(FCUP,FPUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH1P=soil-root convective HPO4 flux per plant in non-band
    !     DIFH1P=soil-root HPO4 diffusion per plant in non-band
    !     CH1P4=HPO4 concentration in non-band
    !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH1P,RTKHP1=HPO4 uptake per plant in non-band lmtd,unlmtd by O2
    !     H1POM,H1POX=minimum,maximum HPO4 available for uptake in non-band
    !     FP14X=fraction of total HPO4 uptake in non-band by root,myco populn
    !     RUPP1P,RUPH1P=HPO4 uptake in non-band unlimited,limited by HPO4
    !     RUOH1P=HPO4 uptake in non-band unlimited by O2
    !     RUCH1P=HPO4 uptake in non-band unlimited by nonstructural C
!
    X=(DIFH1P+RMFH1P)*CH1P4s1(L)
    Y=DIFH1P*UPMNPOs1(N,NZ)
    B=-UPMX-DIFH1P*UPKMPOs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKH1P=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFH1P*UPKMPOs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKHP1=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    H1POM=UPMNPOs1(N,NZ)*VOLWs1(L)*VLPO4s1(L)
    H1POX=AMAX1(0.0,FP14X*(H1PO4s1(L)-H1POM))
    RUPP1Ps1(N,L,NZ)=AMAX1(0.0,RTKH1P*PPs1(NZ))
    RUPH1Ps1(N,L,NZ)=AMIN1(H1POX,RUPP1Ps1(N,L,NZ))
    RUOH1Ps1(N,L,NZ)=AMIN1(H1POX,AMAX1(0.0,RTKHP1*PPs1(NZ)))
    RUCH1Ps1(N,L,NZ)=RUPH1Ps1(N,L,NZ)/FCUP

  ENDIF
  !
  !     HPO4 UPTAKE IN BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH1P4B=HPO4 concentration in band
  !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFH1B=soil-root convective HPO4 flux per plant in band
  !     DIFH1B=soil-root HPO4 diffusion per plant in band
  !
  IF(VLPOBs1(L).GT.ZEROs1.AND.CH1P4Bs1(L) &
    .GT.UPMNPOs1(N,NZ))THEN
    RMFH2B=UPWTRP*CH1P4Bs1(L)*VLPOBs1(L)
    DIFH1B=DIFFL*VLPOBs1(L)
    !
    !     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum HPO4 uptake in band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     TFN4=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=0.1*UPMXPOs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLPOBs1(L)*AMIN1(FCUP,FPUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH1B=soil-root convective HPO4 flux per plant in band
    !     DIFH1B=soil-root HPO4 diffusion per plant in band
    !     CH1P4B=HPO4 concentration in band
    !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH1B,RTKHB1=HPO4 uptake per plant in band lmtd,unlmtd by O2
    !     H1PXM,H1PXB=minimum,maximum HPO4 available for uptake in band
    !     FP1BX=fraction of total HPO4 uptake in band by root,myco populn
    !     RUPP1B,RUPH1B=HPO4 uptake in band unlimited,limited by H2PO4
    !     RUOH1B=HPO4 uptake in band unlimited by O2
    !     RUCH1B=HPO4 uptake in band unlimited by nonstructural C
!
    X=(DIFH1B+RMFH2B)*CH1P4Bs1(L)
    Y=DIFH1B*UPMNPOs1(N,NZ)
    B=-UPMX-DIFH1B*UPKMPOs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKH1B=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFH1B*UPKMPOs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKHB1=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    H1PXM=UPMNPOs1(N,NZ)*VOLWs1(L)*VLPOBs1(L)
    H1PXB=AMAX1(0.0,FP1BX*(H1POBs1(L)-H1PXM))
    RUPP1Bs1(N,L,NZ)=AMAX1(0.0,RTKH1B*PPs1(NZ))
    RUPH1Bs1(N,L,NZ)=AMIN1(H1PXB,RUPP1Bs1(N,L,NZ))
    RUOH1Bs1(N,L,NZ)=AMIN1(H1PXB,AMAX1(0.0,RTKHB1*PPs1(NZ)))
    RUCH1Bs1(N,L,NZ)=RUPH1Bs1(N,L,NZ)/FCUP
  ENDIF
  end associate
  end subroutine UptakeHPO4
!------------------------------------------------------------------------

  subroutine UptakeH2PO4(N,L,NZ,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FWSRT,UPWTRP)

  implicit none
  integer,  intent(in) :: N,L
  integer,  intent(in) :: NZ
  real(r8), intent(in) :: DIFFL
  real(r8), intent(in) :: FPO4X,FPOBX
  real(r8), intent(in) :: FCUP,FPUP,FWSRT,UPWTRP
  real(r8) :: B,C
  real(r8) :: BP,CP
  real(r8) :: DIFH2P,DIFH2B
  real(r8) :: H2POM,H2POX,H2PXM,H2PXB
  real(r8) :: RTKHPB,RMFH2B,RMFH2P,RTKH2P,RTKHPP,RTKH2B
  real(r8) :: UPMX,UPMXP,X,Y
  !
  associate(                             &
    PPs1      => plt_site%PPs1         , &
    ZEROs1    => plt_site%ZEROs1       , &
    WFRs1     => plt_rbgc%WFRs1        , &
    UPKMPOs1  => plt_rbgc%UPKMPOs1     , &
    UPMXPOs1  => plt_rbgc%UPMXPOs1     , &
    UPMNPOs1  => plt_rbgc%UPMNPOs1     , &
    RUCH2Bs1  => plt_rbgc%RUCH2Bs1     , &
    RUOH2Bs1  => plt_rbgc%RUOH2Bs1     , &
    RUPH2Bs1  => plt_rbgc%RUPH2Bs1     , &
    RUOH2Ps1  => plt_rbgc%RUOH2Ps1     , &
    RUPP2Bs1  => plt_rbgc%RUPP2Bs1     , &
    RUCH2Ps1  => plt_rbgc%RUCH2Ps1     , &
    RUPH2Ps1  => plt_rbgc%RUPH2Ps1     , &
    RUPP2Ps1  => plt_rbgc%RUPP2Ps1     , &
    TFN4s1    => plt_pheno%TFN4s1      , &
    RTARPs1   => plt_morph%RTARPs1     , &
    H2POBs1   => plt_soilchem%H2POBs1  , &
    VOLWs1    => plt_soilchem%VOLWs1   , &
    H2PO4s1   => plt_soilchem%H2PO4s1  , &
    CH2P4s1   => plt_soilchem%CH2P4s1  , &
    VLPO4s1   => plt_soilchem%VLPO4s1  , &
    CH2P4Bs1  => plt_soilchem%CH2P4Bs1 , &
    VLPOBs1   => plt_soilchem%VLPOBs1    &
  )
  !     H2PO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH2P4=H2PO4 concentration in non-band
  !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
  !     UPWTRP=root water uptake per plant
  !     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
  !     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
!
    IF(VLPO4s1(L).GT.ZEROs1.AND.CH2P4s1(L) &
      .GT.UPMNPOs1(N,NZ))THEN
      RMFH2P=UPWTRP*CH2P4s1(L)*VLPO4s1(L)
      DIFH2P=DIFFL*VLPO4s1(L)
      !
      !     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
      !     AND FROM ROOT SURFACE AREA, C AND P CONSTRAINTS CALCULATED ABOVE
      !
      !     UPMXP,UPMX=max H2PO4 uptake in non-band unlimited,limited by O2
      !     RTARP=root surface area per plant from grosub.f
      !     FWSRT=protein concentration relative to 5%
      !     TFN4=temperature function for root growth
      !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
      !     WFR=constraint by O2 consumption on all biological processes
!
      UPMXP=UPMXPOs1(N,NZ)*RTARPs1(N,L,NZ) &
        *FWSRT*TFN4s1(L,NZ)*VLPO4s1(L)*AMIN1(FCUP,FPUP)
      UPMX=UPMXP*WFRs1(N,L,NZ)
      !
      !     SOLUTION FOR MASS FLOW + DIFFUSION OF H2PO4 IN AQUEOUS PHASE OF
      !     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY
      !     COMPETITION WITH OTHER ROOT AND MICROBIAL POPULATIONS
      !
      !     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
      !     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
      !     CH2P4=H2PO4 concentration in non-band
      !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
      !     RTKH2P,RTKHPP=H2PO4 uptake per plant in non-band lmtd,unlmtd by O2
      !     H2POM,H2POX=minimum,maximum H2PO4 available for uptake in non-band
      !     FPO4X=fraction of total H2PO4 uptake in non-band by root,myco populn
      !     RUPP2P,RUPH2P=H2PO4 uptake in non-band unlimited,limited by H2PO4
      !     RUOH2P=H2PO4 uptake in non-band unlimited by O2
      !     RUCH2P=H2PO4 uptake in non-band unlimited by nonstructural C
!
      X=(DIFH2P+RMFH2P)*CH2P4s1(L)
      Y=DIFH2P*UPMNPOs1(N,NZ)
      B=-UPMX-DIFH2P*UPKMPOs1(N,NZ)-X+Y
      C=(X-Y)*UPMX
      RTKH2P=(-B-SQRT(B*B-4.0*C))/2.0
      BP=-UPMXP-DIFH2P*UPKMPOs1(N,NZ)-X+Y
      CP=(X-Y)*UPMXP
      RTKHPP=(-BP-SQRT(BP*BP-4.0*CP))/2.0
      H2POM=UPMNPOs1(N,NZ)*VOLWs1(L)*VLPO4s1(L)
      H2POX=AMAX1(0.0,FPO4X*(H2PO4s1(L)-H2POM))
      RUPP2Ps1(N,L,NZ)=AMAX1(0.0,RTKH2P*PPs1(NZ))
      RUPH2Ps1(N,L,NZ)=AMIN1(H2POX,RUPP2Ps1(N,L,NZ))
      RUOH2Ps1(N,L,NZ)=AMIN1(H2POX,AMAX1(0.0_r8,RTKHPP*PPs1(NZ)))
      RUCH2Ps1(N,L,NZ)=RUPH2Ps1(N,L,NZ)/FCUP
    ENDIF
    !
    !     H2PO4 UPTAKE IN BAND SOIL ZONE
    !
    !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
    !     CH2P4B=H2PO4 concentration in band
    !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
    !     UPWTRP=root water uptake per plant
    !     RMFH2B=soil-root convective H2PO4 flux per plant in band
    !     DIFH2B=soil-root H2PO4 diffusion per plant in band
    !

  IF(VLPOBs1(L).GT.ZEROs1.AND.CH2P4Bs1(L) &
    .GT.UPMNPOs1(N,NZ))THEN
    RMFH2B=UPWTRP*CH2P4Bs1(L)*VLPOBs1(L)
    DIFH2B=DIFFL*VLPOBs1(L)
    !
    !     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum H2PO4 uptake in band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     TFN4=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=UPMXPOs1(N,NZ)*RTARPs1(N,L,NZ) &
      *FWSRT*TFN4s1(L,NZ)*VLPOBs1(L)*AMIN1(FCUP,FPUP)
    UPMX=UPMXP*WFRs1(N,L,NZ)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF PO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH2B=soil-root convective H2PO4 flux per plant in band
    !     DIFH2B=soil-root H2PO4 diffusion per plant in band
    !     CH2P4B=H2PO4 concentration in band
    !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH2B,RTKHPB=H2PO4 uptake per plant in band lmtd,unlmtd by O2
    !     H2PXM,H2PXB=minimum,maximum H2PO4 available for uptake in band
    !     FPOBX=fraction of total H2PO4 uptake in band by root,myco populn
    !     RUPP2B,RUPH2B=H2PO4 uptake in band unlimited,limited by H2PO4
    !     RUOH2B=H2PO4 uptake in band unlimited by O2
    !     RUCH2B=H2PO4 uptake in band unlimited by nonstructural C
!
    X=(DIFH2B+RMFH2B)*CH2P4Bs1(L)
    Y=DIFH2B*UPMNPOs1(N,NZ)
    B=-UPMX-DIFH2B*UPKMPOs1(N,NZ)-X+Y
    C=(X-Y)*UPMX
    RTKH2B=(-B-SQRT(B*B-4.0*C))/2.0
    BP=-UPMXP-DIFH2B*UPKMPOs1(N,NZ)-X+Y
    CP=(X-Y)*UPMXP
    RTKHPB=(-BP-SQRT(BP*BP-4.0*CP))/2.0
    H2PXM=UPMNPOs1(N,NZ)*VOLWs1(L)*VLPOBs1(L)
    H2PXB=AMAX1(0.0,FPOBX*(H2POBs1(L)-H2PXM))
    RUPP2Bs1(N,L,NZ)=AMAX1(0.0,RTKH2B*PPs1(NZ))
    RUPH2Bs1(N,L,NZ)=AMIN1(H2PXB,RUPP2Bs1(N,L,NZ))
    RUOH2Bs1(N,L,NZ)=AMIN1(H2PXB,AMAX1(0.0,RTKHPB*PPs1(NZ)))
    RUCH2Bs1(N,L,NZ)=RUPH2Bs1(N,L,NZ)/FCUP
  ENDIF
  end associate
  end subroutine UptakeH2PO4

!------------------------------------------------------------------------

  subroutine UptakeMineralNitrogen(N,L,NZ,PATH,RRADL,FPQ,&
    FPP,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: PATH(2,JZ1),RRADL(2,JZ1),FPQ(2,JZ1,JP1)
  real(r8), intent(in) :: FPP(2,JZ1,JP1),RTARR(2,JZ1)
  real(r8), intent(in) :: FCUP,FZUP,FWSRT,UPWTRP
  real(r8) :: FNO3X,FNOBX,FNH4X,FNHBX
  real(r8) :: TFNH4X,TFNO3X,TFNHBX,TFNOBX

!     begin_execution
  associate(                          &
    ZERO2s1   =>  plt_site%ZERO2s1  , &
    ZEROSs1   =>  plt_site%ZEROSs1  , &
    RUNNXPs1  =>  plt_rbgc%RUNNXPs1 , &
    RUNNHPs1  =>  plt_rbgc%RUNNHPs1 , &
    RUNNBPs1  =>  plt_rbgc%RUNNBPs1 , &
    RUNNOPs1  =>  plt_rbgc%RUNNOPs1 , &
    RNO3Ys1   =>  plt_bgcr%RNO3Ys1  , &
    RNHBYs1   =>  plt_bgcr%RNHBYs1  , &
    RNH4Ys1   =>  plt_bgcr%RNH4Ys1  , &
    RN3BYs1   =>  plt_bgcr%RN3BYs1    &
  )
  TFNH4X=0.0_r8
  TFNHBX=0.0_r8
  TFNO3X=0.0_r8
  TFNOBX=0.0_r8


  IF(RNH4Ys1(L).GT.ZEROSs1)THEN
    FNH4X=AMAX1(FPP(N,L,NZ),RUNNHPs1(N,L,NZ)/RNH4Ys1(L))
  ELSE
    FNH4X=FPQ(N,L,NZ)
  ENDIF
  IF(RNHBYs1(L).GT.ZEROSs1)THEN
    FNHBX=AMAX1(FPP(N,L,NZ),RUNNBPs1(N,L,NZ)/RNHBYs1(L))
  ELSE
    FNHBX=FPQ(N,L,NZ)
  ENDIF

  IF(RNO3Ys1(L).GT.ZEROSs1)THEN
    FNO3X=AMAX1(FPP(N,L,NZ),RUNNOPs1(N,L,NZ)/RNO3Ys1(L))
  ELSE
    FNO3X=FPQ(N,L,NZ)
  ENDIF
  IF(RN3BYs1(L).GT.ZEROSs1)THEN
    FNOBX=AMAX1(FPP(N,L,NZ),RUNNXPs1(N,L,NZ)/RN3BYs1(L))
  ELSE
    FNOBX=FPQ(N,L,NZ)
  ENDIF

  TFNH4X=TFNH4X+FNH4X
  TFNO3X=TFNO3X+FNO3X
  TFNHBX=TFNHBX+FNHBX
  TFNOBX=TFNOBX+FNOBX

  IF(FZUP.GT.ZERO2s1)THEN
!
    !     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NH4,NO3
    !     FROM SOIL TO ROOT
    !
    call UptakeNH4(N,L,NZ,FNH4X,FNHBX,PATH,RRADL,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

    call UptakeNO3(N,L,NZ,FNO3X,FNOBX,PATH,RRADL,RTARR,FCUP,FZUP,FWSRT,UPWTRP)

  ENDIF
  end associate
  end subroutine UptakeMineralNitrogen

!------------------------------------------------------------------------

  subroutine GetUptakeCapcity(N,L,NZ,FPQ,FPP,FCUP,FZUP,FPUP,&
    FWSRT,UPWTRP,UPWTRH,FOXYX)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  REAL(R8), INTENT(IN):: FPQ(2,JZ1,JP1),FPP(2,JZ1,JP1)
  real(r8), intent(out):: FCUP,FZUP,FPUP,FWSRT,UPWTRP,UPWTRH,FOXYX

  associate(                          &
    ROXYYs1   => plt_bgcr%ROXYYs1   , &
    RCO2Ns1   => plt_rbgc%RCO2Ns1   , &
    ROXYPs1   => plt_rbgc%ROXYPs1   , &
    PPs1      => plt_site%PPs1      , &
    ZEROSs1   => plt_site%ZEROSs1   , &
    ZEROs1    => plt_site%ZEROs1    , &
    UPWTRs1   => plt_ew%UPWTRs1     , &
    CWSRTs1   => plt_allom%CWSRTs1  , &
    ZEROPs1   => plt_biom%ZEROPs1   , &
    CWSRTLs1  => plt_biom%CWSRTLs1  , &
    WSRTLs1   => plt_biom%WSRTLs1   , &
    CPPOLRs1  => plt_biom%CPPOLRs1  , &
    CZPOLRs1  => plt_biom%CZPOLRs1  , &
    CPOOLRs1  => plt_biom%CPOOLRs1  , &
    CCPOLRs1  => plt_biom%CCPOLRs1  , &
    WTRTLs1   => plt_biom%WTRTLs1     &
  )
  !
  !     UPTAKE CAPACITY 'FWSRT' DEPENDS ON ROOT,MYCORRHIZAL
  !     PROTEIN CONTENT RELATIVE TO 5% FOR WHICH ACTIVE UPTAKE
  !     PARAMETERS ARE DEFINED
  !
  !     CWSRTL,CWSRT=current,maximum protein concentration
  !     WSRTL,WTRTL=protein content,mass
  !     FWSRT=protein concentration relative to 5%
  !
  IF(WTRTLs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
    CWSRTLs1(N,L,NZ)=AMIN1(CWSRTs1(NZ),WSRTLs1(N,L,NZ)/WTRTLs1(N,L,NZ))
    FWSRT=CWSRTLs1(N,L,NZ)/0.05
  ELSE
    CWSRTLs1(N,L,NZ)=CWSRTs1(NZ)
    FWSRT=1.0
  ENDIF
  !
  !     RESPIRATION CONSTRAINT ON UPTAKE FROM NON-STRUCTURAL C
  !
  !     RCO2N=total respiration from CPOOLR
  !     FCUP=limitation to active uptake respiration from CPOOLR
  !     CPOOLR=nonstructural C content
  !
  IF(RCO2Ns1(N,L,NZ).GT.ZEROPs1(NZ))THEN
    FCUP=AMAX1(0.0,AMIN1(1.0,0.25*safe_adb(CPOOLRs1(N,L,NZ) &
      ,RCO2Ns1(N,L,NZ))))
  ELSE
    FCUP=0.0_r8
  ENDIF
  !
  !     FEEDBACK CONSTRAINT ON N UPTAKE FROM NON-STRUCTURAL N AND P
  !
  !     FZUP,FPUP=limitn to active uptake respiration from CZPOLR,CPPOLR
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition effect on N,P uptake
  !     UPWTRH=water uptake at time step for gas flux calculations
  !
  IF(CCPOLRs1(N,L,NZ).GT.ZEROs1)THEN
    FZUP=AMIN1(safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ)+CZPOLRs1(N,L,NZ)/ZCKI) &
      ,safe_adb(CPPOLRs1(N,L,NZ),CPPOLRs1(N,L,NZ)+CZPOLRs1(N,L,NZ)/ZPKI))
    FPUP=AMIN1(safe_adb(CCPOLRs1(N,L,NZ),CCPOLRs1(N,L,NZ)+CPPOLRs1(N,L,NZ)/PCKI) &
      ,safe_adb(CZPOLRs1(N,L,NZ),CZPOLRs1(N,L,NZ)+CPPOLRs1(N,L,NZ)/PZKI))
  ELSE
    FZUP=0.0_r8
    FPUP=0.0_r8
  ENDIF
  !NN=0
  UPWTRP=AMAX1(0.0,-UPWTRs1(N,L,NZ)/PPs1(NZ))
  UPWTRH=UPWTRP*XNPG
  !
  !     FACTORS CONSTRAINING O2 AND NUTRIENT UPTAKE AMONG
  !     COMPETING ROOT,MYCORRHIZAL AND MICROBIAL POPULATIONS
  !     IN BAND AND NON-BAND SOIL ZONES FROM DEMAND CALCULATED
  !     IN PREVIOUS HOUR
  !
  !     ROXYY=O2 demand by all microbial,root,myco populations
  !     ROXYP=O2 demand by each root,myco population
  !     FOXYX=fraction of ROXYY by each root,myco population
  !     RNH4Y=NH4 demand in non-band by all microbial,root,myco populations
  !     RUNNHP=NH4 demand in non-band by each root,myco population
  !     FNH4X=fraction of RNH4Y by each root,myco populn
  !     RNHBY=NH4 demand in band by all microbial,root,myco populations
  !     RUNNBP=NH4 demand in band by each root,myco population
  !     FNHBX=fraction of RNHBY by each root,myco populn
  !     RNO3Y=NO3 demand in non-band by all microbial,root,myco populations
  !     RUNNOP=NO3 demand in non-band by each root,myco population
  !     FNO3X=fraction of RNO3Y by each root,myco populn
  !     RN3BY=NO3 demand in band by all microbial,root,myco populations
  !     RUNNXB=NO3 demand in band by each root,myco population
  !     FNOBX=fraction of RN3BY by each root,myco populn
  !     RPO4Y=H2PO4 demand in non-band by all microbial,root,myco populations
  !     RUPP2P=H2PO4 demand in non-band by each root,myco population
  !     FPO4X=fraction of RPO4Y by each root,myco populn
  !     RPOBY=H2PO4 demand in band by all microbial,root,myco populations
  !     RUPP2B=H2PO4 demand in band by each root,myco population
  !     FPOBX=fraction of RPOBY by each root,myco populn
  !     RP14Y=HPO4 demand in non-band by all microbial,root,myco populations
  !     RUPP1P=HPO4 demand in non-band by each root,myco population
  !     FP14X=fraction of RP14Y by each root,myco populn
  !     RP1BY=HPO4 demand in band by all microbial,root,myco populations
  !     RUPP1B=HPO4 demand in band by each root,myco population
  !     FP1BX=fraction of RP1BY by each root,myco populn
  !     FPP=minimum uptake fraction
  !     FPQ=PFT fraction of biome root mass
  !
  IF(ROXYYs1(L).GT.ZEROSs1)THEN
    FOXYX=AMAX1(FPP(N,L,NZ),ROXYPs1(N,L,NZ)/ROXYYs1(L))
  ELSE
    FOXYX=FPQ(N,L,NZ)
  ENDIF
  end associate
  end subroutine GetUptakeCapcity
!------------------------------------------------------------------------

  subroutine RootExudates(N,L,NZ)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ

  real(r8) :: CPOOLX,CPOOLT
  real(r8) :: PPOOLX,ZPOOLX
  real(r8) :: VOLWK,VOLWT
  real(r8) :: XFRC,XFRN,XFRP
  integer :: K
  !     begin_execution
  associate(                           &
    CPOOLRs1=>  plt_biom%CPOOLRs1    , &
    ZPOOLRs1=>  plt_biom%ZPOOLRs1    , &
    PPOOLRs1=>  plt_biom%PPOOLRs1    , &
    ZEROPs1 =>  plt_biom%ZEROPs1     , &
    ZEROSs1 =>  plt_site%ZEROSs1     , &
    ZEROS2s1=>  plt_site%ZEROS2s1    , &
    VOLWMs1 =>  plt_site%VOLWMs1     , &
    RTVLWs1 =>  plt_morph%RTVLWs1    , &
    RDFOMCs1=>  plt_rbgc%RDFOMCs1    , &
    RDFOMNs1=>  plt_rbgc%RDFOMNs1    , &
    RDFOMPs1=>  plt_rbgc%RDFOMPs1    , &
    FOSRHs1 =>  plt_soilchem%FOSRHs1 , &
    OQNs1   =>  plt_soilchem%OQNs1   , &
    OQPs1   =>  plt_soilchem%OQPs1   , &
    OQCs1   =>  plt_soilchem%OQCs1     &
  )
  !
  !     ROOT EXUDATION OF C, N AND P DEPENDS ON CONCN DIFFERENCES
  !     BETWEEN ROOT NON-STRUCTURAL POOLS AND SOIL DISSOLVED POOLS
  !
  !     VOLWMM=soil micropore water volume
  !     FOSRH=fraction of total SOC in each substrate K from nitro.f
  !     RTVLW=root aqueous volume
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P in root,myco
  !     XFRC,XFRN,XFRP=nonstructural C,N,P exchg at root-soil DOC equilibrium
  !     OQC=soil DOC
  !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
  !     FEXUC,FEXUN,FEXUP=rate constant for root C,N,P exudation
  !     TLEC,TSHC=total fluxes x blr for calculating canopy air temperature,
  !     vapor pressure in watsub.f
  !      EFLXC,SFLXC=canopylatent,sensible heat fluxes
  !      RA=canopy boundary layer resistance
  !     OSTR=O2 stress indicator
  !
  DO 195 K=0,jcplx11
    VOLWK=VOLWMs1(NPH,L)*FOSRHs1(K,L)
    IF(VOLWK.GT.ZEROS2s1 &
      .AND.RTVLWs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
      VOLWT=VOLWK+RTVLWs1(N,L,NZ)
      CPOOLX=AMIN1(1.25E+03*RTVLWs1(N,L,NZ),CPOOLRs1(N,L,NZ))
      XFRC=(OQCs1(K,L)*RTVLWs1(N,L,NZ)-CPOOLX*VOLWK)/VOLWT
      RDFOMCs1(N,K,L,NZ)=FEXUC*XFRC
      IF(OQCs1(K,L).GT.ZEROSs1 &
        .AND.CPOOLRs1(N,L,NZ).GT.ZEROPs1(NZ))THEN
        CPOOLT=OQCs1(K,L)+CPOOLRs1(N,L,NZ)
        ZPOOLX=0.1*ZPOOLRs1(N,L,NZ)
        PPOOLX=0.1*PPOOLRs1(N,L,NZ)
        XFRN=(OQNs1(K,L)*CPOOLRs1(N,L,NZ)-ZPOOLX*OQCs1(K,L))/CPOOLT
        XFRP=(OQPs1(K,L)*CPOOLRs1(N,L,NZ)-PPOOLX*OQCs1(K,L))/CPOOLT
        RDFOMNs1(N,K,L,NZ)=FEXUN*XFRN
        RDFOMPs1(N,K,L,NZ)=FEXUP*XFRP
      ELSE
        RDFOMNs1(N,K,L,NZ)=0.0_r8
        RDFOMPs1(N,K,L,NZ)=0.0_r8
      ENDIF
    ELSE
      RDFOMCs1(N,K,L,NZ)=0.0_r8
      RDFOMNs1(N,K,L,NZ)=0.0_r8
      RDFOMPs1(N,K,L,NZ)=0.0_r8
    ENDIF

195   CONTINUE
  end associate
  end subroutine RootExudates
!------------------------------------------------------------------------

  subroutine SumupNutrientUptake(N,L,NZ)

  implicit none
  integer, intent(in) :: N, L
  integer, intent(in) :: NZ

  integer :: K
  !     begin_execution
  associate(                              &
    XOQCSs1    =>  plt_bgcr%XOQCSs1     , &
    XOQNSs1    =>  plt_bgcr%XOQNSs1     , &
    XOQPSs1    =>  plt_bgcr%XOQPSs1     , &
    RDFOMCs1   =>  plt_rbgc%RDFOMCs1    , &
    RDFOMNs1   =>  plt_rbgc%RDFOMNs1    , &
    RDFOMPs1   =>  plt_rbgc%RDFOMPs1    , &
    UPOMCs1    =>  plt_rbgc%UPOMCs1     , &
    UPOMNs1    =>  plt_rbgc%UPOMNs1     , &
    UPOMPs1    =>  plt_rbgc%UPOMPs1     , &
    UPH1Ps1    =>  plt_rbgc%UPH1Ps1     , &
    UPH2Ps1    =>  plt_rbgc%UPH2Ps1     , &
    UPNH4s1    =>  plt_rbgc%UPNH4s1     , &
    UPNO3s1    =>  plt_rbgc%UPNO3s1     , &
    RUPH2Bs1   =>  plt_rbgc%RUPH2Bs1    , &
    RUPNHBs1   =>  plt_rbgc%RUPNHBs1    , &
    RUPNOBs1   =>  plt_rbgc%RUPNOBs1    , &
    RUPNH4s1   =>  plt_rbgc%RUPNH4s1    , &
    RUPNO3s1   =>  plt_rbgc%RUPNO3s1    , &
    RUPH2Ps1   =>  plt_rbgc%RUPH2Ps1    , &
    RUPH1Ps1   =>  plt_rbgc%RUPH1Ps1    , &
    RUPH1Bs1   =>  plt_rbgc%RUPH1Bs1      &
  )
  !
  !     TOTAL C,N,P EXCHANGE BETWEEN ROOTS AND SOIL
  !
  !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
  !     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
  !     XOQCS,XOQNZ,XOQPS=accumulated change in DOC,DON,DOP from nitro.f
  !     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
  !     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
  !     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
  !
  DO 295 K=0,jcplx11
    UPOMCs1(NZ)=UPOMCs1(NZ)+RDFOMCs1(N,K,L,NZ)
    UPOMNs1(NZ)=UPOMNs1(NZ)+RDFOMNs1(N,K,L,NZ)
    UPOMPs1(NZ)=UPOMPs1(NZ)+RDFOMPs1(N,K,L,NZ)
    XOQCSs1(K,L)=XOQCSs1(K,L)-RDFOMCs1(N,K,L,NZ)
    XOQNSs1(K,L)=XOQNSs1(K,L)-RDFOMNs1(N,K,L,NZ)
    XOQPSs1(K,L)=XOQPSs1(K,L)-RDFOMPs1(N,K,L,NZ)
295   CONTINUE
  UPNH4s1(NZ)=UPNH4s1(NZ)+RUPNH4s1(N,L,NZ)+RUPNHBs1(N,L,NZ)
  UPNO3s1(NZ)=UPNO3s1(NZ)+RUPNO3s1(N,L,NZ)+RUPNOBs1(N,L,NZ)
  UPH2Ps1(NZ)=UPH2Ps1(NZ)+RUPH2Ps1(N,L,NZ)+RUPH2Bs1(N,L,NZ)
  UPH1Ps1(NZ)=UPH1Ps1(NZ)+RUPH1Ps1(N,L,NZ)+RUPH1Bs1(N,L,NZ)
  end associate
  end subroutine SumupNutrientUptake

end module NutUptakeMod
