module NutUptakeMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use StomatesMod   , only : stomates
  use minimathmod  , only : safe_adb,vapsat,AZMAX1
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use PlantAPIData
  use PlantMathFuncMod
  use RootGasMod
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  public :: PopPlantNutientO2Uptake
  contains

!------------------------------------------------------------------------

  subroutine PopPlantNutientO2Uptake(NZ,FDMP,PopPlantO2Uptake,PopPlantO2Demand,&
    PATH,FineRootRadius,FracPRoot4Uptake,MinFracPRoot4Uptake,FracSoiLayByPrimRoot,RootAreaDivRadius)
  !
  !DESCRIPTION
  !doing plant population level nutrient, and O2 uptake
  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in):: FDMP
  real(r8), intent(in) :: PATH(2,JZ1),FineRootRadius(2,JZ1),FracPRoot4Uptake(2,JZ1,JP1)
  real(r8), intent(in) :: MinFracPRoot4Uptake(2,JZ1,JP1)
  real(r8), intent(in) :: FracSoiLayByPrimRoot(JZ1,JP1),RootAreaDivRadius(2,JZ1)
  real(r8), intent(inout) :: PopPlantO2Uptake,PopPlantO2Demand

  call CanopyNH3Flux(NZ,FDMP)
!
!     ROOT(N=1) AD MYCORRHIZAL(N=2) O2 AND NUTRIENT UPTAKE
!
  call RootMycoO2NutrientUptake(NZ,PopPlantO2Uptake,PopPlantO2Demand,PATH,FineRootRadius,&
    FracPRoot4Uptake,MinFracPRoot4Uptake,FracSoiLayByPrimRoot,RootAreaDivRadius)

  end subroutine PopPlantNutientO2Uptake
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

  associate(                        &
    TCelciusCanopy    =>  plt_ew%TCelciusCanopy       , &
    NU      =>  plt_site%NU       , &
    AREA3   =>  plt_site%AREA3    , &
    CNH3E   =>  plt_site%CNH3E    , &
    RNH3B   =>  plt_rbgc%RNH3B    , &
    ZEROP   =>  plt_biom%ZEROP    , &
    LeafPetioleBiomassC_brch   =>  plt_biom%LeafPetioleBiomassC_brch    , &
    EPOOL   =>  plt_biom%EPOOL    , &
    CEPOLB  =>  plt_biom%CEPOLB   , &
    CanPStomaResistH2O      =>  plt_photo%CanPStomaResistH2O      , &
    CanPbndlResist     =>  plt_photo%CanPbndlResist     , &
    CanopyBranchLeafA_pft   =>  plt_morph%CanopyBranchLeafA_pft   , &
    NumOfBranches_pft     =>  plt_morph%NumOfBranches_pft     , &
    FracPARByCanP   =>  plt_rad%FracPARByCanP     , &
    CanopyLeafA_pft   =>  plt_morph%CanopyLeafA_pft     &
  )
  !
  !     NH3 EXCHANGE BETWEEN CANOPY AND ATMOSPHERE FROM NH3
  !     CONCENTRATION DIFFERENCES 'CNH3E' (ATMOSPHERE FROM 'READS') AND
  !     'CNH3P' (CANOPY), AND FROM STOMATAL + BOUNDARY LAYER RESISTANCE
  !
  !     SNH3P,SNH3X=NH3 solubility at TCelciusCanopy, 25 oC
  !     TCelciusCanopy=canopy temperature (oC)
  !     FDMP,FNH3P=canopy dry matter content,NH3 concentration
  !     CanopyBranchLeafA_pft,CanopyLeafA_pft=branch,canopy leaf area
  !     CNH3P,CNH3E=gaseous NH3 concentration in branch,atmosphere
  !     CZPOLB,ZPOOLB=nonstplt_rbgc%RUCtural N concentration,content in branch
  !     RNH3B=NH3 flux between atmosphere and branch
  !     RA,CanPStomaResistH2O=canopy boundary layer,stomatal resistance
  !     FracPARByCanP=fraction of radiation received by each PFT canopy
  !
  SNH3P=SNH3X*EXP(0.513_r8-0.0171_r8*TCelciusCanopy(NZ))
  FNH3P=1.0E-04_r8*FDMP
  D105: DO NB=1,NumOfBranches_pft(NZ)
    IF(LeafPetioleBiomassC_brch(NB,NZ).GT.ZEROP(NZ).AND.CanopyBranchLeafA_pft(NB,NZ).GT.ZEROP(NZ) &
      .AND.CanopyLeafA_pft(NZ).GT.ZEROP(NZ))THEN
      CNH3P=AZMAX1(FNH3P*CEPOLB(ielmn,NB,NZ)/SNH3P)
      ZPOOLB=AZMAX1(EPOOL(ielmn,NB,NZ))
      RNH3B(NB,NZ)=AMIN1(0.1_r8*ZPOOLB,AMAX1((CNH3E-CNH3P)/(CanPbndlResist(NZ)+CanPStomaResistH2O(NZ)) &
        *FracPARByCanP(NZ)*AREA3(NU)*CanopyBranchLeafA_pft(NB,NZ)/CanopyLeafA_pft(NZ),-0.1_r8*ZPOOLB))
    ELSE
      RNH3B(NB,NZ)=0.0_r8
    ENDIF

  ENDDO D105
  end associate
  end subroutine CanopyNH3Flux

!------------------------------------------------------------------------------------------

  subroutine RootMycoO2NutrientUptake(NZ,PopPlantO2Uptake,PopPlantO2Demand,PATH,FineRootRadius,&
    FracPRoot4Uptake,MinFracPRoot4Uptake,FracSoiLayByPrimRoot,RootAreaDivRadius)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in) :: PATH(2,JZ1),FineRootRadius(2,JZ1),FracPRoot4Uptake(2,JZ1,JP1)
  real(r8),  intent(in) :: MinFracPRoot4Uptake(2,JZ1,JP1)
  real(r8), intent(in) :: FracSoiLayByPrimRoot(JZ1,JP1),RootAreaDivRadius(2,JZ1)
  real(r8), intent(inout) :: PopPlantO2Uptake,PopPlantO2Demand
  real(r8) :: TFOXYX
  real(r8) :: FCUP,FZUP,FPUP,FWSRT,PerPlantRootH2OUptake,CumPerPlantRootH2OUptake,FOXYX,PopPlantO2Uptake_vr
  integer :: N,L
!     begin_execution
  associate(                             &
    THETW =>  plt_soilchem%THETW   , &
    VLSoilPoreMicP  =>  plt_soilchem%VLSoilPoreMicP    , &
    ZEROS2=>  plt_site%ZEROS2      , &
    NU    =>  plt_site%NU          , &
    ZERO  =>  plt_site%ZERO        , &
    ZEROP =>  plt_biom%ZEROP       , &
    ROXYP =>  plt_rbgc%ROXYP       , &
    WFR   =>  plt_rbgc%WFR         , &
    RCO2M =>  plt_rbgc%RCO2M       , &
    RootLenPerP =>  plt_morph%RootLenPerP      , &
    MY    =>  plt_morph%MY         , &
    RootLenDensNLP =>  plt_morph%RootLenDensNLP      , &
    RTVLW =>  plt_morph%RTVLW      , &
    NI    =>  plt_morph%NI           &
  )


  call ZeroUptake(NZ)

  D955: DO N=1,MY(NZ)
    D950: DO L=NU,NI(NZ)
      IF(VLSoilPoreMicP(L).GT.ZEROS2.AND.RootLenDensNLP(N,L,NZ).GT.ZERO &
        .AND.RTVLW(N,L,NZ).GT.ZEROP(NZ).AND.THETW(L).GT.ZERO)THEN
        TFOXYX=0.0_r8
        call GetUptakeCapcity(N,L,NZ,FracPRoot4Uptake,MinFracPRoot4Uptake,FCUP,FZUP,FPUP,&
          FWSRT,PerPlantRootH2OUptake,CumPerPlantRootH2OUptake,FOXYX)

        TFOXYX=TFOXYX+FOXYX
!
!     ROOT O2 DEMAND CALCULATED FROM O2 NON-LIMITED RESPIRATION RATE
!
!     ROXYP=O2 demand, g O2
!     RCO2M=respiration unlimited by O2
!     RTVLW=root or myco aqueous volume
!     FOXYX=fraction of total O2 demand from previous hour
!
        ROXYP(N,L,NZ)=2.667_r8*RCO2M(N,L,NZ)

        call RootSoilGasExchange(N,L,NZ,FineRootRadius,FracPRoot4Uptake,FracSoiLayByPrimRoot,&
          RootAreaDivRadius,CumPerPlantRootH2OUptake,FOXYX,PopPlantO2Uptake_vr)

        PopPlantO2Demand=PopPlantO2Demand+ROXYP(N,L,NZ)
        PopPlantO2Uptake=PopPlantO2Uptake+PopPlantO2Uptake_vr

        call RootExudates(N,L,NZ)
!
!     NUTRIENT UPTAKE
!
!     WFR=constraint by O2 consumption on all biological processes
!     FCUP=limitation to active uptake respiration from CPOOLR
!     FWSRT=protein concentration relative to 5%
!     RootLenPerP=root,myco length per plant
!
        IF(WFR(N,L,NZ).GT.ZERO.AND.FCUP.GT.ZERO.AND.FWSRT.GT.ZERO &
          .AND.RootLenPerP(N,L,NZ).GT.ZEROP(NZ))THEN
!
!     FZUP=limitn to active uptake respiration from CZPOLR
!
          call UptakeMineralNitrogen(N,L,NZ,PATH,FineRootRadius,FracPRoot4Uptake,MinFracPRoot4Uptake,&
            RootAreaDivRadius,FCUP,FZUP,FWSRT,PerPlantRootH2OUptake)
!
!     FPUP=limitn to active uptake respiration from CPPOLR
!
          call UptakeMineralPhosporhus(N,L,NZ,PATH,FineRootRadius,FracPRoot4Uptake,MinFracPRoot4Uptake,&
            RootAreaDivRadius,FCUP,FPUP,FWSRT,PerPlantRootH2OUptake)
        ENDIF
      ELSE

      ENDIF

      call SumupNutrientUptake(N,L,NZ)

    ENDDO D950
  ENDDO D955
  end associate
  end subroutine RootMycoO2NutrientUptake
!------------------------------------------------------------------------
  subroutine ZeroUptake(NZ)

  implicit none
  integer, intent(in) :: NZ

  integer :: K, L1,L2,NN
  !     begin_execution

  L1=plt_site%NU;L2=plt_morph%NI(NZ);NN=plt_morph%MY(NZ)

  plt_rbgc%trcg_air2root_flx_pft_vr(idg_beg:idg_end-1,1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%trcg_Root_DisEvap_flx_vr(idg_beg:idg_end-1,1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RCO2S(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPOXS(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPCHS(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPN2S(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPN3S(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RCO2P(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPOXP(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RDFOME(1:NumOfPlantChemElements,1:NN,1:jcplx,L1:L2,NZ)=0.0_r8
  plt_rbgc%WFR(1:NN,L1:L2,NZ)=1.0
  plt_rbgc%RUNNHP(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPNH4(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUONH4(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCNH4(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUNNBP(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPNHB(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUONHB(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCNHB(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUNNOP(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPNO3(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUONO3(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCNO3(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUNNXP(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPNOB(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUONOB(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCNOB(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPP2P(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPH2P(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUOH2P(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCH2P(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPP2B(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPH2B(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUOH2B(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCH2B(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPP1P(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPH1P(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUOH1P(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCH1P(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPP1B(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUPH1B(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUOH1B(1:NN,L1:L2,NZ)=0.0_r8
  plt_rbgc%RUCH1B(1:NN,L1:L2,NZ)=0.0_r8
  plt_bgcr%RUPNF(L1:L2,NZ)=0.0_r8
  end subroutine ZeroUptake

!------------------------------------------------------------------------0

  subroutine UptakeMineralPhosporhus(N,L,NZ,PATH,FineRootRadius,FracPRoot4Uptake,MinFracPRoot4Uptake,&
    RootAreaDivRadius,FCUP,FPUP,FWSRT,PerPlantRootH2OUptake)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  real(r8), intent(in):: PATH(2,JZ1),FineRootRadius(2,JZ1),FracPRoot4Uptake(2,JZ1,JP1)
  real(r8), intent(in) :: MinFracPRoot4Uptake(2,JZ1,JP1)
  real(r8), intent(in):: RootAreaDivRadius(2,JZ1),FCUP,FPUP,FWSRT,PerPlantRootH2OUptake
  real(r8) :: TFPO4X,TFP14X,TFPOBX,TFP1BX
  real(r8) :: DIFFL
  real(r8) :: FP14X,FP1BX
  real(r8) :: FPO4X,FPOBX
  real(r8) :: POSGX
  real(r8) :: PATHL
  associate(                               &
    TortMicPM    =>   plt_site%TortMicPM       , &
    ZEROS   =>   plt_site%ZEROS      , &
    ZERO2   =>   plt_site%ZERO2      , &
    RUPP1B  =>   plt_rbgc%RUPP1B     , &
    RUPP1P  =>   plt_rbgc%RUPP1P     , &
    RUPP2P  =>   plt_rbgc%RUPP2P     , &
    RUPP2B  =>   plt_rbgc%RUPP2B     , &
    SolDifc =>   plt_soilchem%SolDifc, &
    RP1BY   =>   plt_bgcr%RP1BY      , &
    RPOBY   =>   plt_bgcr%RPOBY      , &
    RP14Y   =>   plt_bgcr%RP14Y      , &
    RPO4Y   =>   plt_bgcr%RPO4Y        &
  )
  TFPO4X=0.0_r8
  TFPOBX=0.0_r8
  TFP14X=0.0_r8
  TFP1BX=0.0_r8

!     begin_execution
  IF(RPO4Y(L).GT.ZEROS)THEN
    FPO4X=AMAX1(MinFracPRoot4Uptake(N,L,NZ),RUPP2P(N,L,NZ)/RPO4Y(L))
  ELSE
    FPO4X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RPOBY(L).GT.ZEROS)THEN
    FPOBX=AMAX1(MinFracPRoot4Uptake(N,L,NZ),RUPP2B(N,L,NZ)/RPOBY(L))
  ELSE
    FPOBX=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RP14Y(L).GT.ZEROS)THEN
    FP14X=AMAX1(MinFracPRoot4Uptake(N,L,NZ),RUPP1P(N,L,NZ)/RP14Y(L))
  ELSE
    FP14X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RP1BY(L).GT.ZEROS)THEN
    FP1BX=AMAX1(MinFracPRoot4Uptake(N,L,NZ),RUPP1B(N,L,NZ)/RP1BY(L))
  ELSE
    FP1BX=FracPRoot4Uptake(N,L,NZ)
  ENDIF

  TFPO4X=TFPO4X+FPO4X
  TFP14X=TFP14X+FP14X

  TFPOBX=TFPOBX+FPOBX
  TFP1BX=TFP1BX+FP1BX

  IF(FPUP.GT.ZERO2)THEN
    !
    !     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF H2PO4,HPO4
    !     FROM SOIL TO ROOT
    !
    !     POSGL=PO4 diffusivity
    !     TortMicPM=soil tortuosity
    !     PATH=path length of water and nutrient uptake
    !     FineRootRadius=root radius
    !     DIFFL=PO4 diffusion per plant
    !
    POSGX=SolDifc(ids_H1PO4,L)*TortMicPM(NPH,L)
    PATHL=AMIN1(PATH(N,L),FineRootRadius(N,L)+SQRT(2.0*POSGX))
    DIFFL=POSGX*safe_adb(RootAreaDivRadius(N,L),LOG(PATHL/FineRootRadius(N,L)))

    call UptakeH2PO4(N,L,NZ,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FWSRT,PerPlantRootH2OUptake)

    call UptakeHPO4(N,L,NZ,DIFFL,FP14X,FP1BX,FCUP,FPUP,FWSRT,PerPlantRootH2OUptake)

  ENDIF
  end associate
  end subroutine UptakeMineralPhosporhus
!------------------------------------------------------------------------

  subroutine UptakeNO3(N,L,NZ,FNO3X,FNOBX,PATH,FineRootRadius,RootAreaDivRadius,FCUP,FZUP,&
    FWSRT,PerPlantRootH2OUptake)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  real(r8), intent(in):: FNO3X,FNOBX,PATH(2,JZ1),FineRootRadius(2,JZ1),RootAreaDivRadius(2,JZ1)
  real(r8), intent(in):: FCUP,FZUP,FWSRT,PerPlantRootH2OUptake

  real(r8) :: DIFFL
  real(r8) :: DIFNO3,DIFNOB
  real(r8) :: PATHL
  real(r8) :: RMFNO3,RTKNO3,RTKNOP,RMFNOB,RTKNOB,RTKNPB
  real(r8) :: UPMX,UPMXP
  real(r8) :: ZOSGX,ZNO3M,ZNO3X,ZNOBM,ZNOBX
  type(PlantSoluteUptakeConfig_type) :: PlantSoluteUptakeConfig

! begin_execution
  associate(                              &
    pftPlantPopulation      =>  plt_site%pftPlantPopulation         , &
    TortMicPM    =>  plt_site%TortMicPM       , &
    ZERO    =>  plt_site%ZERO       , &
    RTARP   =>  plt_morph%RTARP     , &
    fTgrowRootP    =>  plt_pheno%fTgrowRootP      , &
    WFR     =>  plt_rbgc%WFR        , &
    UPMXZO  =>  plt_rbgc%UPMXZO     , &
    UPMNZO  =>  plt_rbgc%UPMNZO     , &
    UPKMZO  =>  plt_rbgc%UPKMZO     , &
    RUCNOB  =>  plt_rbgc%RUCNOB     , &
    RUNNOP  =>  plt_rbgc%RUNNOP     , &
    RUCNO3  =>  plt_rbgc%RUCNO3     , &
    RUPNO3  =>  plt_rbgc%RUPNO3     , &
    RUNNXP  =>  plt_rbgc%RUNNXP     , &
    RUONOB  =>  plt_rbgc%RUONOB     , &
    RUPNOB  =>  plt_rbgc%RUPNOB     , &
    RUONO3  =>  plt_rbgc%RUONO3     , &
    SolDifc =>  plt_soilchem%SolDifc, &
    trc_solml   => plt_soilchem%trc_solml  , &
    trcs_VLN   =>  plt_soilchem%trcs_VLN  , &
    trc_solcl  =>  plt_soilchem%trc_solcl , &
    VLWatMicP    =>  plt_soilchem%VLWatMicP     &
  )
!
! PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NO3
! FROM SOIL TO ROOT
!
! ZOSGL=NO3 diffusivity
! TortMicPM=soil tortuosity
! FineRootRadius=root radius
! PATH=path length of water and nutrient uptake
! DIFFL=NO3 diffusion per plant
!
  ZOSGX=SolDifc(ids_NO3,L)*TortMicPM(NPH,L)
  PATHL=AMIN1(PATH(N,L),FineRootRadius(N,L)+SQRT(2.0_r8*ZOSGX))
  DIFFL=ZOSGX*safe_adb(RootAreaDivRadius(N,L),LOG(PATHL/FineRootRadius(N,L)))
  !
  ! NO3 UPTAKE IN NON-BAND SOIL ZONE
  !
  !  VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
  !     CNO3S=NO3 concentration in non-band
  !     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFNO3=soil-root convective NO3 flux per plant in non-band
  !     DIFNO3=soil-root NO3 diffusion per plant in non-band
  !
  IF(trcs_VLN(ids_NO3,L).GT.ZERO.AND.trc_solcl(ids_NO3,L).GT.UPMNZO(N,NZ))THEN
    RMFNO3=PerPlantRootH2OUptake*trcs_VLN(ids_NO3,L)
    DIFNO3=DIFFL*trcs_VLN(ids_NO3,L)
    !
    !     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=max NO3 uptake in non-band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     fTgrowRootP=temperature function for root growth
    !     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=UPMXZO(N,NZ)*RTARP(N,L,NZ) &
      *FWSRT*fTgrowRootP(L,NZ)*trcs_VLN(ids_NO3,L)*AMIN1(FCUP,FZUP)

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
    ZNO3M=UPMNZO(N,NZ)*VLWatMicP(L)*trcs_VLN(ids_NO3,L)
    ZNO3X=AZMAX1(FNO3X*(trc_solml(ids_NO3,L)-ZNO3M))

    PlantSoluteUptakeConfig%SoluteConcMin=UPMNZO(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx = RMFNO3
    PlantSoluteUptakeConfig%SolDifusFlx = DIFNO3
    PlantSoluteUptakeConfig%UptakeRateMax = UPMXP
    PlantSoluteUptakeConfig%O2Stress  = WFR(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc =trc_solml(ids_NO3,L)
    PlantSoluteUptakeConfig%SoluteMassMax=ZNO3X
    PlantSoluteUptakeConfig%CAvailStress=FCUP
    PlantSoluteUptakeConfig%PlantPopulation= pftPlantPopulation(NZ)
    PlantSoluteUptakeConfig%SoluteKM=UPKMZO(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RUNNOP(N,L,NZ),RUPNO3(N,L,NZ),&
      RUPNOB(N,L,NZ),RUCNO3(N,L,NZ))

  ENDIF
  !
  !     NO3 UPTAKE IN BAND SOIL ZONE
  !
  !     VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
  !     CNO3B=NO3 concentration in band
  !     UPMXZO,UPKMZO,UPMNZO=NO3 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFNOB=soil-root convective NO3 flux per plant in band
  !     DIFNOB=soil-root NO3 diffusion per plant in band
  !

  IF(trcs_VLN(ids_NO3B,L).GT.ZERO.AND.trc_solcl(ids_NO3B,L).GT.UPMNZO(N,NZ))THEN
    RMFNOB=PerPlantRootH2OUptake*trcs_VLN(ids_NO3B,L)
    DIFNOB=DIFFL*trcs_VLN(ids_NO3B,L)
    !
    !     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum NO3 uptake in band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     fTgrowRootP=temperature function for root growth
    !     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=UPMXZO(N,NZ)*RTARP(N,L,NZ) &
      *FWSRT*fTgrowRootP(L,NZ)*trcs_VLN(ids_NO3B,L)*AMIN1(FCUP,FZUP)
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
    ZNOBM=UPMNZO(N,NZ)*VLWatMicP(L)*trcs_VLN(ids_NO3B,L)
    ZNOBX=AZMAX1(FNOBX*(trc_solml(ids_NO3B,L)-ZNOBM))

    PlantSoluteUptakeConfig%SoluteConcMin=UPMNZO(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx = RMFNOB
    PlantSoluteUptakeConfig%SolDifusFlx = DIFNOB
    PlantSoluteUptakeConfig%UptakeRateMax = UPMXP
    PlantSoluteUptakeConfig%O2Stress  = WFR(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc =trc_solcl(ids_NO3B,L)
    PlantSoluteUptakeConfig%SoluteMassMax=ZNOBX
    PlantSoluteUptakeConfig%CAvailStress=FCUP
    PlantSoluteUptakeConfig%PlantPopulation= pftPlantPopulation(NZ)
    PlantSoluteUptakeConfig%SoluteKM=UPKMZO(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RUNNXP(N,L,NZ),RUONOB(N,L,NZ),&
      RUPNOB(N,L,NZ),RUCNOB(N,L,NZ))

  ENDIF
  end associate
  end subroutine UptakeNO3
!------------------------------------------------------------------------
  subroutine UptakeNH4(N,L,NZ,FNH4X,FNHBX,PATH,FineRootRadius,RootAreaDivRadius,&
    FCUP,FZUP,FWSRT,PerPlantRootH2OUptake)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: FNH4X,FNHBX,PATH(2,JZ1),FineRootRadius(2,JZ1)
  real(r8), intent(in) :: RootAreaDivRadius(2,JZ1)
  real(r8), intent(in) :: FCUP,FZUP,FWSRT,PerPlantRootH2OUptake
  real(r8) :: DIFFL
  real(r8) :: DIFNH4,DIFNHB
  real(r8) :: PATHL
  real(r8) :: RMFNH4,RTKNH4,RTKNHP,RMFNHB,RTKNHB,RTKNBP
  real(r8) :: UPMX,UPMXP
  real(r8) :: ZNHBX,ZNSGX,ZNH4M,ZNH4X,ZNHBM
  type(PlantSoluteUptakeConfig_type) :: PlantSoluteUptakeConfig
! begin_execution
  associate(                              &
    pftPlantPopulation      =>  plt_site%pftPlantPopulation         , &
    ZERO    =>  plt_site%ZERO       , &
    TortMicPM    =>  plt_site%TortMicPM       , &
    fTgrowRootP    =>  plt_pheno%fTgrowRootP      , &
    WFR     =>  plt_rbgc%WFR        , &
    UPMNZH  =>  plt_rbgc%UPMNZH     , &
    UPMXZH  =>  plt_rbgc%UPMXZH     , &
    UPKMZH  =>  plt_rbgc%UPKMZH     , &
    RUPNHB  =>  plt_rbgc%RUPNHB     , &
    RUCNH4  =>  plt_rbgc%RUCNH4     , &
    RUONHB  =>  plt_rbgc%RUONHB     , &
    RUONH4  =>  plt_rbgc%RUONH4     , &
    RUNNHP  =>  plt_rbgc%RUNNHP     , &
    RUCNHB  =>  plt_rbgc%RUCNHB     , &
    RUNNBP  =>  plt_rbgc%RUNNBP     , &
    RUPNH4  =>  plt_rbgc%RUPNH4     , &
    RTARP   =>  plt_morph%RTARP     , &
    SolDifc =>  plt_soilchem%SolDifc, &
    VLWatMicP    =>  plt_soilchem%VLWatMicP   , &
    trcs_VLN   =>  plt_soilchem%trcs_VLN  , &
    trc_solml  =>  plt_soilchem%trc_solml , &
    trc_solcl   =>  plt_soilchem%trc_solcl    &
  )
! ZNSGL=NH4 diffusivity
! TortMicPM=soil tortuosity
! PATH=path length of water and nutrient uptake
! FineRootRadius=root radius
! DIFFL=NH4 diffusion per plant
!
  ZNSGX=SolDifc(idg_NH3,L)*TortMicPM(NPH,L)
  PATHL=AMIN1(PATH(N,L),FineRootRadius(N,L)+SQRT(2.0*ZNSGX))
  DIFFL=ZNSGX*safe_adb(RootAreaDivRadius(N,L),LOG(PATHL/FineRootRadius(N,L)))
!
! NH4 UPTAKE IN NON-BAND SOIL ZONE
!
! VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
! CNH4S=NH4 concentration in non-band
! UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
! PerPlantRootH2OUptake=root water uptake per plant
! RMFNH4=soil-root convective NH4 flux per plant in non-band
! DIFNH4=soil-root NH4 diffusion per plant in non-band
!
  IF(trcs_VLN(ids_NH4,L).GT.ZERO.AND.trc_solcl(ids_NH4,L).GT.UPMNZH(N,NZ))THEN
    RMFNH4=PerPlantRootH2OUptake*trcs_VLN(ids_NH4,L)
    DIFNH4=DIFFL*trcs_VLN(ids_NH4,L)
!
!   NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
!   AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
!
!   UPMXP,UPMX=max NH4 uptake in non-band unlimited,limited by O2
!   RTARP=root surface area per plant from grosub.f
!   FWSRT=protein concentration relative to 5%
!   fTgrowRootP=temperature function for root growth
!   FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
!   WFR=constraint by O2 consumption on all biological processes
!
    UPMXP=UPMXZH(N,NZ)*RTARP(N,L,NZ) &
      *FWSRT*fTgrowRootP(L,NZ)*trcs_VLN(ids_NH4,L)*AMIN1(FCUP,FZUP)
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
    ZNH4M=UPMNZH(N,NZ)*VLWatMicP(L)*trcs_VLN(ids_NH4,L)
    ZNH4X=AZMAX1(FNH4X*(trc_solml(ids_NH4,L)-ZNH4M))

    PlantSoluteUptakeConfig%SoluteConcMin=UPMNZH(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx = RMFNH4
    PlantSoluteUptakeConfig%SolDifusFlx = DIFNH4
    PlantSoluteUptakeConfig%UptakeRateMax = UPMXP
    PlantSoluteUptakeConfig%O2Stress  = WFR(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc =trc_solcl(ids_NH4,L)
    PlantSoluteUptakeConfig%SoluteMassMax=ZNH4X
    PlantSoluteUptakeConfig%CAvailStress=FCUP
    PlantSoluteUptakeConfig%PlantPopulation= pftPlantPopulation(NZ)
    PlantSoluteUptakeConfig%SoluteKM=UPKMZH(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RUNNHP(N,L,NZ),RUONH4(N,L,NZ),&
      RUPNH4(N,L,NZ),RUCNH4(N,L,NZ))

  ENDIF
!
! NH4 UPTAKE IN BAND SOIL ZONE
!
! VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
! CNH4B=NH4 concentration in band
! UPMXZH,UPKMZH,UPMNZH=NH4 max uptake,Km,min concn from PFT file
! PerPlantRootH2OUptake=root water uptake per plant
! RMFNHB=soil-root convective NH4 flux per plant in band
! DIFNHB=soil-root NH4 diffusion per plant in band
!

  IF(trcs_VLN(ids_NH4B,L).GT.ZERO.AND.trc_solcl(ids_NH4B,L).GT.UPMNZH(N,NZ))THEN
    RMFNHB=PerPlantRootH2OUptake*trcs_VLN(ids_NH4B,L)
    DIFNHB=DIFFL*trcs_VLN(ids_NH4B,L)
!
!   NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
!   AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
!
!   UPMXP,UPMX=maximum NH4 uptake in band unlimited,limited by O2
!   RTARP=root surface area per plant from grosub.f
!   FWSRT=protein concentration relative to 5%
!   fTgrowRootP=temperature function for root growth
!   FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
!   WFR=constraint by O2 consumption on all biological processes
!
    UPMXP=UPMXZH(N,NZ)*RTARP(N,L,NZ) &
      *FWSRT*fTgrowRootP(L,NZ)*trcs_VLN(ids_NH4B,L)*AMIN1(FCUP,FZUP)
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

    ZNHBM=UPMNZH(N,NZ)*VLWatMicP(L)*trcs_VLN(ids_NH4B,L)
    ZNHBX=AZMAX1(FNHBX*(trc_solml(ids_NH4B,L)-ZNHBM))

    PlantSoluteUptakeConfig%SoluteConcMin=UPMNZH(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx = RMFNHB
    PlantSoluteUptakeConfig%SolDifusFlx = DIFNHB
    PlantSoluteUptakeConfig%UptakeRateMax = UPMXP
    PlantSoluteUptakeConfig%O2Stress  = WFR(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc =trc_solcl(ids_NH4B,L)
    PlantSoluteUptakeConfig%SoluteMassMax=ZNHBX
    PlantSoluteUptakeConfig%CAvailStress=FCUP
    PlantSoluteUptakeConfig%PlantPopulation= pftPlantPopulation(NZ)
    PlantSoluteUptakeConfig%SoluteKM=UPKMZH(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RUNNBP(N,L,NZ),RUONHB(N,L,NZ),&
      RUPNHB(N,L,NZ),RUCNHB(N,L,NZ))

  ENDIF
  end associate
  end subroutine UptakeNH4

!------------------------------------------------------------------------
  subroutine UptakeHPO4(N,L,NZ,DIFFL,FP14X,FP1BX,FCUP,FPUP,FWSRT,PerPlantRootH2OUptake)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: DIFFL
  real(r8), intent(in) :: FP14X,FP1BX
  real(r8), intent(in) :: FCUP,FPUP,FWSRT,PerPlantRootH2OUptake
  real(r8) :: B,C
  real(r8) :: BP,CP
  real(r8) :: DIFH1P
  real(r8) :: DIFH1B
  real(r8) :: H1POM,H1POX,H1PXM,H1PXB
  real(r8) :: RMFH1P,RTKH1P,RTKHP1,RTKH1B,RTKHB1
  real(r8) :: RMFH2B
  real(r8) :: UPMX,UPMXP
  real(r8) :: X,Y
  type(PlantSoluteUptakeConfig_type) :: PlantSoluteUptakeConfig
  !     begin_execution
  associate(                             &
    pftPlantPopulation      => plt_site%pftPlantPopulation         , &
    ZERO    => plt_site%ZERO       , &
    fTgrowRootP    => plt_pheno%fTgrowRootP      , &
    WFR     => plt_rbgc%WFR        , &
    UPMNPO  => plt_rbgc%UPMNPO     , &
    UPKMPO  => plt_rbgc%UPKMPO     , &
    UPMXPO  => plt_rbgc%UPMXPO     , &
    RUPH1B  => plt_rbgc%RUPH1B     , &
    RUCH1P  => plt_rbgc%RUCH1P     , &
    RUOH1P  => plt_rbgc%RUOH1P     , &
    RUPP1P  => plt_rbgc%RUPP1P     , &
    RUPH1P  => plt_rbgc%RUPH1P     , &
    RUCH1B  => plt_rbgc%RUCH1B     , &
    RUOH1B  => plt_rbgc%RUOH1B     , &
    RUPP1B  => plt_rbgc%RUPP1B     , &
    RTARP   => plt_morph%RTARP     , &
    VLWatMicP    => plt_soilchem%VLWatMicP   , &
    trcs_VLN   => plt_soilchem%trcs_VLN  , &
    trc_solcl   => plt_soilchem%trc_solcl  , &
    trc_solml   => plt_soilchem%trc_solml    &

  )
  !
  !     HPO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH1P4=HPO4 concentration in non-band
  !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant > 0._r8
  !     RMFH1P=soil-root convective HPO4 flux per plant in non-band
  !     DIFH1P=soil-root HPO4 diffusion per plant in non-band > 0._r8
!
    IF(trcs_VLN(ids_H1PO4,L).GT.ZERO.AND.trc_solcl(ids_H1PO4,L).GT.UPMNPO(N,NZ))THEN
      RMFH1P=PerPlantRootH2OUptake*trcs_VLN(ids_H1PO4,L)
      DIFH1P=DIFFL*trcs_VLN(ids_H1PO4,L)
    !
    !     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=max HPO4 uptake in non-band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     fTgrowRootP=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=0.1_r8*UPMXPO(N,NZ)*RTARP(N,L,NZ) &
      *FWSRT*fTgrowRootP(L,NZ)*trcs_VLN(ids_H1PO4,L)*AMIN1(FCUP,FPUP)
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

    H1POM=UPMNPO(N,NZ)*VLWatMicP(L)*trcs_VLN(ids_H1PO4,L)
    H1POX=AZMAX1(FP14X*(trc_solml(ids_H1PO4,L)-H1POM))

    PlantSoluteUptakeConfig%SoluteConcMin=UPMNPO(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx = RMFH1P
    PlantSoluteUptakeConfig%SolDifusFlx = DIFH1P
    PlantSoluteUptakeConfig%UptakeRateMax = UPMXP
    PlantSoluteUptakeConfig%O2Stress  = WFR(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc =trc_solcl(ids_H1PO4,L)
    PlantSoluteUptakeConfig%SoluteMassMax=H1POX
    PlantSoluteUptakeConfig%CAvailStress=FCUP
    PlantSoluteUptakeConfig%PlantPopulation= pftPlantPopulation(NZ)
    PlantSoluteUptakeConfig%SoluteKM=UPKMPO(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RUPP1P(N,L,NZ),RUOH1P(N,L,NZ),&
      RUPH1P(N,L,NZ),RUCH1P(N,L,NZ))

  ENDIF
  !
  !     HPO4 UPTAKE IN BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH1P4B=HPO4 concentration in band
  !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFH1B=soil-root convective HPO4 flux per plant in band
  !     DIFH1B=soil-root HPO4 diffusion per plant in band
  !
  IF(trcs_VLN(ids_H1PO4B,L).GT.ZERO.AND.trc_solcl(ids_H1PO4B,L).GT.UPMNPO(N,NZ))THEN
    RMFH2B=PerPlantRootH2OUptake*trcs_VLN(ids_H1PO4B,L)
    DIFH1B=DIFFL*trcs_VLN(ids_H1PO4B,L)
    !
    !     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum HPO4 uptake in band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     fTgrowRootP=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=0.1_r8*UPMXPO(N,NZ)*RTARP(N,L,NZ) &
      *FWSRT*fTgrowRootP(L,NZ)*trcs_VLN(ids_H1PO4B,L)*AMIN1(FCUP,FPUP)
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

    H1PXM=UPMNPO(N,NZ)*VLWatMicP(L)*trcs_VLN(ids_H1PO4B,L)
    H1PXB=AZMAX1(FP1BX*(trc_solml(ids_H1PO4B,L)-H1PXM))

    PlantSoluteUptakeConfig%SoluteConcMin=UPMNPO(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx = RMFH2B
    PlantSoluteUptakeConfig%SolDifusFlx = DIFH1B
    PlantSoluteUptakeConfig%UptakeRateMax = UPMXP
    PlantSoluteUptakeConfig%O2Stress  = WFR(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc =trc_solcl(ids_H1PO4B,L)
    PlantSoluteUptakeConfig%SoluteMassMax=H1PXB
    PlantSoluteUptakeConfig%CAvailStress=FCUP
    PlantSoluteUptakeConfig%PlantPopulation= pftPlantPopulation(NZ)
    PlantSoluteUptakeConfig%SoluteKM=UPKMPO(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RUPP1B(N,L,NZ),RUOH1B(N,L,NZ),&
      RUPH1B(N,L,NZ),RUCH1B(N,L,NZ))

  ENDIF
  end associate
  end subroutine UptakeHPO4
!------------------------------------------------------------------------

  subroutine UptakeH2PO4(N,L,NZ,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FWSRT,PerPlantRootH2OUptake)

  implicit none
  integer,  intent(in) :: N,L
  integer,  intent(in) :: NZ
  real(r8), intent(in) :: DIFFL
  real(r8), intent(in) :: FPO4X,FPOBX
  real(r8), intent(in) :: FCUP,FPUP,FWSRT,PerPlantRootH2OUptake

  real(r8) :: DIFH2P,DIFH2B
  real(r8) :: H2POM,H2POX,H2PXM,H2PXB
  real(r8) :: RTKHPB,RMFH2B,RMFH2P,RTKH2P,RTKHPP,RTKH2B
  real(r8) :: UPMX,UPMXP
  type(PlantSoluteUptakeConfig_type) :: PlantSoluteUptakeConfig
  !
  associate(                             &
    pftPlantPopulation      => plt_site%pftPlantPopulation         , &
    ZERO    => plt_site%ZERO       , &
    WFR     => plt_rbgc%WFR        , &
    UPKMPO  => plt_rbgc%UPKMPO     , &
    UPMXPO  => plt_rbgc%UPMXPO     , &
    UPMNPO  => plt_rbgc%UPMNPO     , &
    RUCH2B  => plt_rbgc%RUCH2B     , &
    RUOH2B  => plt_rbgc%RUOH2B     , &
    RUPH2B  => plt_rbgc%RUPH2B     , &
    RUOH2P  => plt_rbgc%RUOH2P     , &
    RUPP2B  => plt_rbgc%RUPP2B     , &
    RUCH2P  => plt_rbgc%RUCH2P     , &
    RUPH2P  => plt_rbgc%RUPH2P     , &
    RUPP2P  => plt_rbgc%RUPP2P     , &
    fTgrowRootP    => plt_pheno%fTgrowRootP      , &
    RTARP   => plt_morph%RTARP     , &
    trc_solml   => plt_soilchem%trc_solml  , &
    VLWatMicP    => plt_soilchem%VLWatMicP   , &
    trc_solcl   => plt_soilchem%trc_solcl  , &
    trcs_VLN   => plt_soilchem%trcs_VLN  &
  )
  !     H2PO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH2P4=H2PO4 concentration in non-band
  !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
  !     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
!
    IF(trcs_VLN(ids_H1PO4,L).GT.ZERO.AND.trc_solcl(ids_H2PO4,L).GT.UPMNPO(N,NZ))THEN
      RMFH2P=PerPlantRootH2OUptake*trcs_VLN(ids_H1PO4,L)
      DIFH2P=DIFFL*trcs_VLN(ids_H1PO4,L)
      !
      !     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
      !     AND FROM ROOT SURFACE AREA, C AND P CONSTRAINTS CALCULATED ABOVE
      !
      !     UPMXP,UPMX=max H2PO4 uptake in non-band unlimited,limited by O2
      !     RTARP=root surface area per plant from grosub.f
      !     FWSRT=protein concentration relative to 5%
      !     fTgrowRootP=temperature function for root growth
      !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
      !     WFR=constraint by O2 consumption on all biological processes
!
      UPMXP=UPMXPO(N,NZ)*RTARP(N,L,NZ) &
        *FWSRT*fTgrowRootP(L,NZ)*trcs_VLN(ids_H1PO4,L)*AMIN1(FCUP,FPUP)
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

      H2POM=UPMNPO(N,NZ)*VLWatMicP(L)*trcs_VLN(ids_H1PO4,L)
      H2POX=AZMAX1(FPO4X*(trc_solml(ids_H2PO4,L)-H2POM))

      PlantSoluteUptakeConfig%SoluteConcMin=UPMNPO(N,NZ)
      PlantSoluteUptakeConfig%SolAdvFlx = RMFH2P
      PlantSoluteUptakeConfig%SolDifusFlx = DIFH2P
      PlantSoluteUptakeConfig%UptakeRateMax = UPMXP
      PlantSoluteUptakeConfig%O2Stress  = WFR(N,L,NZ)
      PlantSoluteUptakeConfig%SoluteConc =trc_solcl(ids_H2PO4,L)
      PlantSoluteUptakeConfig%SoluteMassMax=H2POX
      PlantSoluteUptakeConfig%CAvailStress=FCUP
      PlantSoluteUptakeConfig%PlantPopulation= pftPlantPopulation(NZ)
      PlantSoluteUptakeConfig%SoluteKM=UPKMPO(N,NZ)

      call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RUPP2P(N,L,NZ),RUOH2P(N,L,NZ),&
        RUPH2P(N,L,NZ),RUCH2P(N,L,NZ))

    ENDIF
    !
    !     H2PO4 UPTAKE IN BAND SOIL ZONE
    !
    !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
    !     CH2P4B=H2PO4 concentration in band
    !     UPMXPO,UPKMPO,UPMNPO=H2PO4 max uptake,Km,min concn from PFT file
    !     PerPlantRootH2OUptake=root water uptake per plant
    !     RMFH2B=soil-root convective H2PO4 flux per plant in band
    !     DIFH2B=soil-root H2PO4 diffusion per plant in band
    !

  IF(trcs_VLN(ids_H1PO4B,L).GT.ZERO.AND.trc_solcl(ids_H2PO4B,L).GT.UPMNPO(N,NZ))THEN
    RMFH2B=PerPlantRootH2OUptake*trcs_VLN(ids_H1PO4B,L)
    DIFH2B=DIFFL*trcs_VLN(ids_H1PO4B,L)
    !
    !     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum H2PO4 uptake in band unlimited,limited by O2
    !     RTARP=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     fTgrowRootP=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     WFR=constraint by O2 consumption on all biological processes
    !
    UPMXP=UPMXPO(N,NZ)*RTARP(N,L,NZ) &
      *FWSRT*fTgrowRootP(L,NZ)*trcs_VLN(ids_H1PO4B,L)*AMIN1(FCUP,FPUP)

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

    H2PXM=UPMNPO(N,NZ)*VLWatMicP(L)*trcs_VLN(ids_H1PO4B,L)
    H2PXB=AZMAX1(FPOBX*(trc_solml(ids_H2PO4B,L)-H2PXM))

    PlantSoluteUptakeConfig%SoluteConcMin=UPMNPO(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx = RMFH2B
    PlantSoluteUptakeConfig%SolDifusFlx = DIFH2B
    PlantSoluteUptakeConfig%UptakeRateMax = UPMXP
    PlantSoluteUptakeConfig%O2Stress  = WFR(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc =trc_solcl(ids_H2PO4B,L)
    PlantSoluteUptakeConfig%SoluteMassMax=H2PXB
    PlantSoluteUptakeConfig%CAvailStress=FCUP
    PlantSoluteUptakeConfig%PlantPopulation= pftPlantPopulation(NZ)
    PlantSoluteUptakeConfig%SoluteKM=UPKMPO(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RUPP2B(N,L,NZ),RUOH2B(N,L,NZ),&
      RUPH2B(N,L,NZ),RUCH2B(N,L,NZ))

  ENDIF
  end associate
  end subroutine UptakeH2PO4

!------------------------------------------------------------------------

  subroutine UptakeMineralNitrogen(N,L,NZ,PATH,FineRootRadius,FracPRoot4Uptake,&
    MinFracPRoot4Uptake,RootAreaDivRadius,FCUP,FZUP,FWSRT,PerPlantRootH2OUptake)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: PATH(2,JZ1),FineRootRadius(2,JZ1),FracPRoot4Uptake(2,JZ1,JP1)
  real(r8), intent(in) :: MinFracPRoot4Uptake(2,JZ1,JP1),RootAreaDivRadius(2,JZ1)
  real(r8), intent(in) :: FCUP,FZUP,FWSRT,PerPlantRootH2OUptake
  real(r8) :: FNO3X,FNOBX,FNH4X,FNHBX
  real(r8) :: TFNH4X,TFNO3X,TFNHBX,TFNOBX

!     begin_execution
  associate(                          &
    ZERO2   =>  plt_site%ZERO2  , &
    ZEROS   =>  plt_site%ZEROS  , &
    RUNNXP  =>  plt_rbgc%RUNNXP , &
    RUNNHP  =>  plt_rbgc%RUNNHP , &
    RUNNBP  =>  plt_rbgc%RUNNBP , &
    RUNNOP  =>  plt_rbgc%RUNNOP , &
    RNO3Y   =>  plt_bgcr%RNO3Y  , &
    RNHBY   =>  plt_bgcr%RNHBY  , &
    RNH4Y   =>  plt_bgcr%RNH4Y  , &
    RN3BY   =>  plt_bgcr%RN3BY    &
  )
  TFNH4X=0.0_r8
  TFNHBX=0.0_r8
  TFNO3X=0.0_r8
  TFNOBX=0.0_r8


  IF(RNH4Y(L).GT.ZEROS)THEN
    FNH4X=AMAX1(MinFracPRoot4Uptake(N,L,NZ),RUNNHP(N,L,NZ)/RNH4Y(L))
  ELSE
    FNH4X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RNHBY(L).GT.ZEROS)THEN
    FNHBX=AMAX1(MinFracPRoot4Uptake(N,L,NZ),RUNNBP(N,L,NZ)/RNHBY(L))
  ELSE
    FNHBX=FracPRoot4Uptake(N,L,NZ)
  ENDIF

  IF(RNO3Y(L).GT.ZEROS)THEN
    FNO3X=AMAX1(MinFracPRoot4Uptake(N,L,NZ),RUNNOP(N,L,NZ)/RNO3Y(L))
  ELSE
    FNO3X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RN3BY(L).GT.ZEROS)THEN
    FNOBX=AMAX1(MinFracPRoot4Uptake(N,L,NZ),RUNNXP(N,L,NZ)/RN3BY(L))
  ELSE
    FNOBX=FracPRoot4Uptake(N,L,NZ)
  ENDIF

  TFNH4X=TFNH4X+FNH4X
  TFNO3X=TFNO3X+FNO3X
  TFNHBX=TFNHBX+FNHBX
  TFNOBX=TFNOBX+FNOBX

  IF(FZUP.GT.ZERO2)THEN
!
    !     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NH4,NO3
    !     FROM SOIL TO ROOT
    !
    call UptakeNH4(N,L,NZ,FNH4X,FNHBX,PATH,FineRootRadius,RootAreaDivRadius,FCUP,FZUP,FWSRT,&
      PerPlantRootH2OUptake)

    call UptakeNO3(N,L,NZ,FNO3X,FNOBX,PATH,FineRootRadius,RootAreaDivRadius,FCUP,FZUP,FWSRT,&
      PerPlantRootH2OUptake)

  ENDIF
  end associate
  end subroutine UptakeMineralNitrogen

!------------------------------------------------------------------------

  subroutine GetUptakeCapcity(N,L,NZ,FracPRoot4Uptake,MinFracPRoot4Uptake,FCUP,FZUP,FPUP,&
    FWSRT,PerPlantRootH2OUptake,CumPerPlantRootH2OUptake,FOXYX)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  REAL(R8), INTENT(IN):: FracPRoot4Uptake(2,JZ1,JP1),MinFracPRoot4Uptake(2,JZ1,JP1)
  real(r8), intent(out):: FCUP,FZUP,FPUP,FWSRT,PerPlantRootH2OUptake,CumPerPlantRootH2OUptake,FOXYX

  associate(                          &
    ROXYY   => plt_bgcr%ROXYY   , &
    RCO2N   => plt_rbgc%RCO2N   , &
    ROXYP   => plt_rbgc%ROXYP   , &
    pftPlantPopulation      => plt_site%pftPlantPopulation      , &
    ZEROS   => plt_site%ZEROS   , &
    ZERO    => plt_site%ZERO    , &
    AllPlantRootH2OUptake_vr   => plt_ew%AllPlantRootH2OUptake_vr     , &
    RootFracRemobilizableBiom   => plt_allom%RootFracRemobilizableBiom  , &
    ZEROP   => plt_biom%ZEROP   , &
    RootProteinConc_pftvr  => plt_biom%RootProteinConc_pftvr  , &
    WSRTL   => plt_biom%WSRTL   , &
    EPOOLR  => plt_biom%EPOOLR  , &
    RootNonstructElementConcpft_vr  => plt_biom%RootNonstructElementConcpft_vr  , &
    WTRTL   => plt_biom%WTRTL     &
  )
  !
  !     UPTAKE CAPACITY 'FWSRT' DEPENDS ON ROOT,MYCORRHIZAL
  !     PROTEIN CONTENT RELATIVE TO 5% FOR WHICH ACTIVE UPTAKE
  !     PARAMETERS ARE DEFINED
  !
  !     RootProteinConc_pftvr,RootFracRemobilizableBiom=current,maximum protein concentration
  !     WSRTL,WTRTL=protein content,mass
  !     FWSRT=protein concentration relative to 5%
  !
  IF(WTRTL(N,L,NZ).GT.ZEROP(NZ))THEN
    RootProteinConc_pftvr(N,L,NZ)=AMIN1(RootFracRemobilizableBiom(NZ),WSRTL(N,L,NZ)/WTRTL(N,L,NZ))
    FWSRT=RootProteinConc_pftvr(N,L,NZ)/0.05_r8
  ELSE
    RootProteinConc_pftvr(N,L,NZ)=RootFracRemobilizableBiom(NZ)
    FWSRT=1.0_r8
  ENDIF
  !
  !     RESPIRATION CONSTRAINT ON UPTAKE FROM NON-STRUCTURAL C
  !
  !     RCO2N=total respiration from CPOOLR
  !     FCUP=limitation to active uptake respiration from CPOOLR
  !     CPOOLR=nonstructural C content
  !
  IF(RCO2N(N,L,NZ).GT.ZEROP(NZ))THEN
    FCUP=AZMAX1(AMIN1(1.0_r8,0.25_r8*safe_adb(EPOOLR(ielmc,N,L,NZ),RCO2N(N,L,NZ))))
  ELSE
    FCUP=0.0_r8
  ENDIF
  !
  !     FEEDBACK CONSTRAINT ON N UPTAKE FROM NON-STRUCTURAL N AND P
  !
  !     FZUP,FPUP=limitn to active uptake respiration from CZPOLR,CPPOLR
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition effect on N,P uptake
  !     CumPerPlantRootH2OUptake=water uptake at time step for gas flux calculations
  !
  IF(RootNonstructElementConcpft_vr(ielmc,N,L,NZ).GT.ZERO)THEN
    FZUP=AMIN1(safe_adb(RootNonstructElementConcpft_vr(ielmc,N,L,NZ),RootNonstructElementConcpft_vr(ielmc,N,L,NZ) &
      +RootNonstructElementConcpft_vr(ielmn,N,L,NZ)/ZCKI) &
      ,safe_adb(RootNonstructElementConcpft_vr(ielmp,N,L,NZ),RootNonstructElementConcpft_vr(ielmp,N,L,NZ) &
       +RootNonstructElementConcpft_vr(ielmn,N,L,NZ)/ZPKI))
    FPUP=AMIN1(safe_adb(RootNonstructElementConcpft_vr(ielmc,N,L,NZ),RootNonstructElementConcpft_vr(ielmc,N,L,NZ) &
      +RootNonstructElementConcpft_vr(ielmp,N,L,NZ)/PCKI) &
      ,safe_adb(RootNonstructElementConcpft_vr(ielmn,N,L,NZ),RootNonstructElementConcpft_vr(ielmn,N,L,NZ) &
      +RootNonstructElementConcpft_vr(ielmp,N,L,NZ)/PZKI))
  ELSE
    FZUP=0.0_r8
    FPUP=0.0_r8
  ENDIF
  !NN=0
  PerPlantRootH2OUptake=AZMAX1(-AllPlantRootH2OUptake_vr(N,L,NZ)/pftPlantPopulation(NZ))
  CumPerPlantRootH2OUptake=PerPlantRootH2OUptake*dts_gas
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
  !     MinFracPRoot4Uptake=minimum uptake fraction
  !     FracPRoot4Uptake=PFT fraction of biome root mass
  !
  IF(ROXYY(L).GT.ZEROS)THEN
    FOXYX=AMAX1(MinFracPRoot4Uptake(N,L,NZ),ROXYP(N,L,NZ)/ROXYY(L))
  ELSE
    FOXYX=FracPRoot4Uptake(N,L,NZ)
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
  real(r8) :: VLWatMicPK,VLWatMicPT
  real(r8) :: XFRC,XFRN,XFRP
  integer :: K
  !     begin_execution
  associate(                       &
    EPOOLR=>  plt_biom%EPOOLR    , &
    ZEROP =>  plt_biom%ZEROP     , &
    ZEROS =>  plt_site%ZEROS     , &
    ZEROS2=>  plt_site%ZEROS2    , &
    VLWatMicPM =>  plt_site%VLWatMicPM     , &
    RTVLW =>  plt_morph%RTVLW    , &
    RDFOME=>  plt_rbgc%RDFOME    , &
    FOSRH =>  plt_soilchem%FOSRH , &
    DOM   =>  plt_soilchem%DOM     &
  )
  !
  !     ROOT EXUDATION OF C, N AND P DEPENDS ON CONCN DIFFERENCES
  !     BETWEEN ROOT NON-STRUCTURAL POOLS AND SOIL DISSOLVED POOLS
  !
  !     VLWatMicPMM=soil micropore water volume
  !     FOSRH=fraction of total SOC in each substrate K from nitro.f
  !     RTVLW=root aqueous volume
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P in root,myco
  !     XFRC,XFRN,XFRP=nonstructural C,N,P exchg at root-soil DOC equilibrium
  !     OQC=soil DOC
  !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
  !     FEXUC,FEXUN,FEXUP=rate constant for root C,N,P exudation
  !     Canopy_Heat_Latent_col,Canopy_Heat_Sens_col=total fluxes x blr for calculating canopy air temperature,
  !     vapor pressure in watsub.f
  !      EvapTransHeatP,SFLXC=canopylatent,sensible heat fluxes
  !      RA=canopy boundary layer resistance
  !     OSTR=O2 stress indicator
  !
  D195: DO K=1,jcplx
    VLWatMicPK=VLWatMicPM(NPH,L)*FOSRH(K,L)
    IF(VLWatMicPK.GT.ZEROS2.AND.RTVLW(N,L,NZ).GT.ZEROP(NZ))THEN
      VLWatMicPT=VLWatMicPK+RTVLW(N,L,NZ)
      CPOOLX=AMIN1(1.25E+03_r8*RTVLW(N,L,NZ),EPOOLR(ielmc,N,L,NZ))
      XFRC=(DOM(idom_doc,K,L)*RTVLW(N,L,NZ)-CPOOLX*VLWatMicPK)/VLWatMicPT
      RDFOME(ielmc,N,K,L,NZ)=FEXUC*XFRC
      IF(DOM(idom_doc,K,L).GT.ZEROS.AND.EPOOLR(ielmc,N,L,NZ).GT.ZEROP(NZ))THEN
        CPOOLT=DOM(idom_doc,K,L)+EPOOLR(ielmc,N,L,NZ)
        ZPOOLX=0.1_r8*EPOOLR(ielmn,N,L,NZ)
        PPOOLX=0.1_r8*EPOOLR(ielmp,N,L,NZ)
        XFRN=(DOM(idom_don,K,L)*EPOOLR(ielmc,N,L,NZ)-ZPOOLX*DOM(idom_doc,K,L))/CPOOLT
        XFRP=(DOM(idom_dop,K,L)*EPOOLR(ielmc,N,L,NZ)-PPOOLX*DOM(idom_doc,K,L))/CPOOLT
        RDFOME(ielmn,N,K,L,NZ)=FEXUN*XFRN
        RDFOME(ielmp,N,K,L,NZ)=FEXUP*XFRP
      ELSE
        RDFOME(ielmn,N,K,L,NZ)=0.0_r8
        RDFOME(ielmp,N,K,L,NZ)=0.0_r8
      ENDIF
    ELSE
      RDFOME(1:NumOfPlantChemElements,N,K,L,NZ)=0.0_r8
    ENDIF

  ENDDO D195
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
    RDOM_micb_flx    =>  plt_bgcr%RDOM_micb_flx     , &
    RDFOME   =>  plt_rbgc%RDFOME    , &
    UPOME    =>  plt_rbgc%UPOME     , &
    UPH1P    =>  plt_rbgc%UPH1P     , &
    UPH2P    =>  plt_rbgc%UPH2P     , &
    UPNH4    =>  plt_rbgc%UPNH4     , &
    UPNO3    =>  plt_rbgc%UPNO3     , &
    RUPH2B   =>  plt_rbgc%RUPH2B    , &
    RUPNHB   =>  plt_rbgc%RUPNHB    , &
    RUPNOB   =>  plt_rbgc%RUPNOB    , &
    RUPNH4   =>  plt_rbgc%RUPNH4    , &
    RUPNO3   =>  plt_rbgc%RUPNO3    , &
    RUPH2P   =>  plt_rbgc%RUPH2P    , &
    RUPH1P   =>  plt_rbgc%RUPH1P    , &
    RUPH1B   =>  plt_rbgc%RUPH1B      &
  )
  !
  !     TOTAL C,N,P EXCHANGE BETWEEN ROOTS AND SOIL
  !
  !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
  !     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
  !     RUPNH4,RUPNHB,RUPN03,RUPNOB=uptake from non-band,band of NH4,NO3
  !     RUPH2P,RUPH2B,RUPH1P,RUPH1B=uptake from non-band,band of H2PO4,HPO4
  !     UPNH4,UPNO3,UPH2P,UPH1P=PFT uptake of NH4,NO3,H2PO4,HPO4
  !
  D295: DO K=1,jcplx
    UPOME(ielmc,NZ)=UPOME(ielmc,NZ)+RDFOME(ielmc,N,K,L,NZ)
    UPOME(ielmn,NZ)=UPOME(ielmn,NZ)+RDFOME(ielmn,N,K,L,NZ)
    UPOME(ielmp,NZ)=UPOME(ielmp,NZ)+RDFOME(ielmp,N,K,L,NZ)
    RDOM_micb_flx(idom_doc,K,L)=RDOM_micb_flx(idom_doc,K,L)-RDFOME(ielmc,N,K,L,NZ)
    RDOM_micb_flx(idom_don,K,L)=RDOM_micb_flx(idom_don,K,L)-RDFOME(ielmn,N,K,L,NZ)
    RDOM_micb_flx(idom_dop,K,L)=RDOM_micb_flx(idom_dop,K,L)-RDFOME(ielmp,N,K,L,NZ)
  ENDDO D295
  UPNH4(NZ)=UPNH4(NZ)+RUPNH4(N,L,NZ)+RUPNHB(N,L,NZ)
  UPNO3(NZ)=UPNO3(NZ)+RUPNO3(N,L,NZ)+RUPNOB(N,L,NZ)
  UPH2P(NZ)=UPH2P(NZ)+RUPH2P(N,L,NZ)+RUPH2B(N,L,NZ)
  UPH1P(NZ)=UPH1P(NZ)+RUPH1P(N,L,NZ)+RUPH1B(N,L,NZ)
  end associate
  end subroutine SumupNutrientUptake

end module NutUptakeMod
