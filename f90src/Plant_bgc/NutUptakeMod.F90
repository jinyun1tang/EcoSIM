module NutUptakeMod

  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: safe_adb, vapsat, AZMAX1,dssign
  use TracerPropMod, only: gas_solubility
  use DebugToolMod, only : PrintInfo
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
  public :: PlantNutientO2Uptake
  public :: ZeroNutrientUptake
  contains

!------------------------------------------------------------------------

  subroutine PlantNutientO2Uptake(I,J,NZ,FDMP, PathLen_pvr,FineRootRadius,FracPRoot4Uptake,&
    MinFracPRoot4Uptake_pvr,FracSoiLayByPrimRoot,RootAreaDivRadius_vr)
  !
  !DESCRIPTION
  !doing plant population level nutrient, and O2 uptake
  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in):: FDMP
  real(r8), intent(in) :: PathLen_pvr(jroots,JZ1),FineRootRadius(jroots,JZ1),FracPRoot4Uptake(jroots,JZ1,JP1)
  real(r8), intent(in) :: MinFracPRoot4Uptake_pvr(jroots,JZ1,JP1)
  real(r8), intent(in) :: FracSoiLayByPrimRoot(JZ1,JP1),RootAreaDivRadius_vr(jroots,JZ1)
  
  real(r8)  :: PopPlantO2Uptake,PopPlantO2Demand

  associate(                                                  &
    PlantO2Stress_pft      => plt_pheno%PlantO2Stress_pft,    &
    ZERO4Groth_pft         => plt_biom%ZERO4Groth_pft         &    
  )

  call FoliarNutrientInterception(NZ,FDMP)
!
!     ROOT(N=1) AD MYCORRHIZAL(N=2) O2 AND NUTRIENT UPTAKE
!
  call RootMycoO2NutrientUptake(I,J,NZ,PathLen_pvr,FineRootRadius,&
    FracPRoot4Uptake,MinFracPRoot4Uptake_pvr,FracSoiLayByPrimRoot,RootAreaDivRadius_vr &
    ,PopPlantO2Uptake,PopPlantO2Demand)

  call SumNutrientUptake(NZ)

  IF(PopPlantO2Demand.GT.ZERO4Groth_pft(NZ))THEN
    PlantO2Stress_pft(NZ)=AZMAX1(PopPlantO2Uptake)/PopPlantO2Demand
  ELSE
    PlantO2Stress_pft(NZ)=1.0_r8
  ENDIF    
  end associate
  end subroutine PlantNutientO2Uptake
!------------------------------------------------------------------------

  subroutine FoliarNutrientInterception(NZ,FDMP)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in):: FDMP

  character(len=*), parameter :: subname='FoliarNutrientInterception'
  real(r8) :: FNH3P
  real(r8) :: CNH3P
  real(r8) :: SNH3P
  real(r8) :: ZPOOLB
  integer :: NB

  associate(                                                         &
    TdegCCanopy_pft           => plt_ew%TdegCCanopy_pft,             &
    NU                        => plt_site%NU,                        &
    AREA3                     => plt_site%AREA3,                     &
    NH3Dep2Can_brch           => plt_rbgc%NH3Dep2Can_brch,           &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,            &
    AtmGasc                   => plt_site%AtmGasc,                   &
    LeafPetolBiomassC_brch    => plt_biom%LeafPetolBiomassC_brch,    &
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch,      &
    LeafPetoNonstElmConc_brch => plt_biom%LeafPetoNonstElmConc_brch, &
    CanPStomaResistH2O_pft    => plt_photo%CanPStomaResistH2O_pft,   &
    CanopyBndlResist_pft      => plt_photo%CanopyBndlResist_pft,     &
    LeafAreaLive_brch         => plt_morph%LeafAreaLive_brch,        &
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft,        &
    FracPARads2Canopy_pft     => plt_rad%FracPARads2Canopy_pft,      &
    CanopyLeafArea_pft        => plt_morph%CanopyLeafArea_pft        &
  )

  call PrintInfo('beg '//subname)
  !
  !     NH3 EXCHANGE BETWEEN CANOPY AND ATMOSPHERE FROM NH3
  !     CONCENTRATION DIFFERENCES 'AtmGasc(idg_NH3)' (ATMOSPHERE FROM 'READS') AND
  !     'CNH3P' (CANOPY), AND FROM STOMATAL + BOUNDARY LAYER RESISTANCE
  !
  !     SNH3P =NH3 solubility at TdegCCanopy_pft
  !     TdegCCanopy_pft=canopy temperature (oC)
  !     FDMP,FNH3P=canopy dry matter content,NH3 concentration
  !     LeafAreaLive_brch,CanopyLeafArea_pft=branch,canopy leaf area
  !     CNH3P,AtmGasc(idg_NH3)=gaseous NH3 concentration in branch,atmosphere
  !     CZPOLB,ZPOOLB=nonstplt_rbgc%RUCtural N concentration,content in branch
  !     NH3Dep2Can_brch=NH3 flux between atmosphere and branch
  !     RA,CanPStomaResistH2O_pft=canopy boundary layer,stomatal resistance
  !     FracPARads2Canopy_pft=fraction of radiation received by each PFT canopy
  !
  SNH3P = gas_solubility(idg_NH3,TdegCCanopy_pft(NZ))
  FNH3P = 1.0E-04_r8*FDMP
  if(FracPARads2Canopy_pft(NZ).GT.ZERO4Groth_pft(NZ))then
    D105: DO NB=1,NumOfBranches_pft(NZ)
      IF(LeafPetolBiomassC_brch(NB,NZ).GT.ZERO4Groth_pft(NZ).AND.LeafAreaLive_brch(NB,NZ).GT.ZERO4Groth_pft(NZ) &
        .AND.CanopyLeafArea_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
        CNH3P                  = AZMAX1(FNH3P*LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/SNH3P)
        ZPOOLB                 = AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))
        NH3Dep2Can_brch(NB,NZ) = AMIN1(0.1_r8*ZPOOLB &
          ,AMAX1((AtmGasc(idg_NH3)-CNH3P)/(CanopyBndlResist_pft(NZ)+CanPStomaResistH2O_pft(NZ)) &
          *FracPARads2Canopy_pft(NZ)*AREA3(NU)*LeafAreaLive_brch(NB,NZ)/CanopyLeafArea_pft(NZ),-0.1_r8*ZPOOLB))
      ELSE
        NH3Dep2Can_brch(NB,NZ)=0.0_r8
      ENDIF
    ENDDO D105
  ELSE
    NH3Dep2Can_brch(1:NumOfBranches_pft(NZ),NZ)=0.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine FoliarNutrientInterception

!------------------------------------------------------------------------------------------

  subroutine RootMycoO2NutrientUptake(I,J,NZ,PathLen_pvr,FineRootRadius,FracPRoot4Uptake,&
    MinFracPRoot4Uptake_pvr,FracSoiLayByPrimRoot,RootAreaDivRadius_vr,PopPlantO2Uptake,PopPlantO2Demand)

  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: PathLen_pvr(jroots,JZ1)
  real(r8), intent(in) :: FineRootRadius(jroots,JZ1)
  real(r8), intent(in) :: FracPRoot4Uptake(jroots,JZ1,JP1)
  real(r8), intent(in) :: MinFracPRoot4Uptake_pvr(jroots,JZ1,JP1)   !minimum active fraction of roots for substance uptake from soil [none]
  real(r8), intent(in) :: FracSoiLayByPrimRoot(JZ1,JP1)
  real(r8), intent(in) :: RootAreaDivRadius_vr(jroots,JZ1)
  real(r8), intent(out) :: PopPlantO2Uptake
  real(r8), intent(out) :: PopPlantO2Demand

  character(len=*), parameter :: subname='RootMycoO2NutrientUptake'
  real(r8) :: FCUP,FZUP,FPUP,FWSRT,PerPlantRootH2OUptake
  real(r8) :: dtPerPlantRootH2OUptake,FOXYX,PopPlantO2Uptake_vr
  real(r8) :: RootCO2Ar,RootCO2ArB
  integer :: N,L,idg
!     begin_execution
  associate(                                                      &
    THETW_vr                => plt_soilchem%THETW_vr,             &
    VLSoilPoreMicP_vr       => plt_soilchem%VLSoilPoreMicP_vr,    &
    ZEROS2                  => plt_site%ZEROS2,                   &
    NU                      => plt_site%NU,                       &
    NK                      => plt_site%NK,                       &    
    ZERO                    => plt_site%ZERO,                     &
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft,           &
    RootO2Dmnd4Resp_pvr     => plt_rbgc%RootO2Dmnd4Resp_pvr,      &
    RAutoRootO2Limter_rpvr  => plt_rbgc%RAutoRootO2Limter_rpvr,   &
    RootRespPotent_pvr      => plt_rbgc%RootRespPotent_pvr,       &
    RootLenPerPlant_pvr     => plt_morph%RootLenPerPlant_pvr,     &
    MY                      => plt_morph%MY,                      &
    RootLenDensPerPlant_pvr => plt_morph%RootLenDensPerPlant_pvr, &
    trcs_deadroot2soil_pvr  => plt_rbgc%trcs_deadroot2soil_pvr,   &
    trcg_rootml_pvr         => plt_rbgc%trcg_rootml_pvr,          &
    trcs_rootml_pvr         => plt_rbgc%trcs_rootml_pvr,          &
    RootVH2O_pvr            => plt_morph%RootVH2O_pvr,            &
    RootCO2Ar2Soil_pvr      => plt_rbgc%RootCO2Ar2Soil_pvr,       &
    RootCO2Ar2Root_pvr      => plt_rbgc%RootCO2Ar2Root_pvr,       &
    MaxSoiL4Root_pft        => plt_morph%MaxSoiL4Root_pft         &
  )

  call PrintInfo('beg '//subname)

  PopPlantO2Uptake = 0._r8
  PopPlantO2Demand = 0._r8
  RootCO2Ar        = 0._r8
  RootCO2ArB       = 0._r8
  trcs_deadroot2soil_pvr(:,:,NZ) = 0._r8  

  D950: DO L=NU,NK
    IF(VLSoilPoreMicP_vr(L).GT.ZEROS2 .AND. THETW_vr(L).GT.ZERO) then

      D955: DO N  = 1, MY(NZ)
        if(RootLenDensPerPlant_pvr(N,L,NZ).GT.ZERO .AND. RootVH2O_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN

          call GetUptakeCapcity(N,L,NZ,FracPRoot4Uptake,MinFracPRoot4Uptake_pvr,FCUP,FZUP,FPUP,&
            FWSRT,PerPlantRootH2OUptake,dtPerPlantRootH2OUptake,FOXYX)

!
!     ROOT O2 DEMAND CALCULATED FROM O2 NON-LIMITED RESPIRATION RATE
!
!     RootO2Dmnd4Resp_pvr=O2 demand, g O2
!     RootRespPotent_pvr=respiration unlimited by O2
!     RootVH2O_pvr=root or myco aqueous volume
!     FOXYX=fraction of total O2 demand from previous hour
!
          RootO2Dmnd4Resp_pvr(N,L,NZ)=2.667_r8*RootRespPotent_pvr(N,L,NZ)

          call RootSoilGasExchange(I,J,N,L,NZ,FineRootRadius,FracPRoot4Uptake,FracSoiLayByPrimRoot,&
            RootAreaDivRadius_vr,dtPerPlantRootH2OUptake,FOXYX,PopPlantO2Uptake_vr)
        
          PopPlantO2Demand = PopPlantO2Demand+RootO2Dmnd4Resp_pvr(N,L,NZ)
          PopPlantO2Uptake = PopPlantO2Uptake+PopPlantO2Uptake_vr
       
          call RootExudates(I,J,N,L,NZ)
!
!     NUTRIENT UPTAKE
!
!     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
!     FCUP=limitation to active uptake respiration from CPOOLR
!     FWSRT=protein concentration relative to 5%
!     RootLenPerPlant_pvr=root,myco length per plant
!
          IF(RAutoRootO2Limter_rpvr(N,L,NZ).GT.ZERO .AND. FCUP.GT.ZERO .AND. FWSRT.GT.ZERO &
            .AND. RootLenPerPlant_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
!
!     FZUP=limitn to active uptake respiration from CZPOLR
!         
!          if(I>176)print*,'uptakemin'
            call UptakeMineralNitrogen(I,J,N,L,NZ,PathLen_pvr,FineRootRadius,FracPRoot4Uptake,MinFracPRoot4Uptake_pvr,&
              RootAreaDivRadius_vr,FCUP,FZUP,FWSRT,PerPlantRootH2OUptake)
!
!     FPUP=limitn to active uptake respiration from CPPOLR
!         
!          if(I>176)print*,'uptakeppp'
            call UptakeMineralPhosporhus(N,L,NZ,PathLen_pvr,FineRootRadius,FracPRoot4Uptake,MinFracPRoot4Uptake_pvr,&
              RootAreaDivRadius_vr,FCUP,FPUP,FWSRT,PerPlantRootH2OUptake)
          ENDIF
          RootCO2Ar=RootCO2Ar-plt_rbgc%RootCO2AutorX_pvr(N,L,NZ)
        ENDIF
!        if(I==140 .and. J<=2)write(116,*)'rootgas',(I*1000+J)*100+N,L,plt_rbgc%trcg_air2root_flx_pvr(idg_CH4,N,L,NZ)
      ENDDO D955      
    ELSE
      D956: DO N  = 1, MY(NZ)    
        RootCO2ArB=RootCO2ArB-plt_rbgc%RootCO2AutorX_pvr(N,L,NZ)
        RootCO2Ar2Soil_pvr(L,NZ)=RootCO2Ar2Soil_pvr(L,NZ)-plt_rbgc%RootCO2AutorX_pvr(N,L,NZ)
        DO idg=idg_beg,idg_NH3
          trcs_deadroot2soil_pvr(idg,L,NZ) = trcs_deadroot2soil_pvr(idg,L,NZ) + trcg_rootml_pvr(idg,N,L,NZ)
          trcs_deadroot2soil_pvr(idg,L,NZ) = trcs_deadroot2soil_pvr(idg,L,NZ) + trcs_rootml_pvr(idg,N,L,NZ)
          trcg_rootml_pvr(idg,N,L,NZ)      = 0._r8
          trcs_rootml_pvr(idg,N,L,NZ)      = 0._r8
        ENDDO
      ENDDO D956
    ENDIF
  ENDDO D950

  call PrintInfo('end '//subname)
  end associate
  end subroutine RootMycoO2NutrientUptake
!------------------------------------------------------------------------
  subroutine ZeroNutrientUptake(NZ)

  implicit none
  integer, intent(in) :: NZ

  integer :: K, L1,L2,NN
  !     begin_execution

  L1=plt_site%NU;L2=plt_site%NK;NN=plt_morph%MY(NZ)
  plt_rbgc%RootCO2Ar2Soil_pvr(:,NZ)       = 0._r8
  plt_rbgc%RootCO2Ar2Root_pvr(:,NZ)       = 0._r8
  plt_rbgc%trcg_air2root_flx_pvr(idg_beg:idg_NH3,1:NN,L1:L2,NZ)        = 0.0_r8
  plt_rbgc%trcg_Root_gas2aqu_flx_vr(idg_beg:idg_NH3,1:NN,L1:L2,NZ)     = 0.0_r8
  plt_rbgc%RootUptkSoiSol_pvr(idg_beg:idg_end,1:NN,L1:L2,NZ)              = 0.0_r8
  plt_rbgc%RootCO2Emis_pvr(1:NN,L1:L2,NZ)                                = 0.0_r8
  plt_rbgc%RootMycoExudEUptk_pvr(1:NumPlantChemElms,1:NN,1:jcplx,L1:L2,NZ) = 0.0_r8
  plt_rbgc%RAutoRootO2Limter_rpvr(1:NN,L1:L2,NZ)                          = 1.0_r8
  plt_rbgc%RootNH4DmndSoil_pvr(1:NN,L1:L2,NZ)                            = 0.0_r8
  plt_rbgc%RootNutUptake_pvr(ids_NH4,1:NN,L1:L2,NZ)                      = 0.0_r8
  plt_rbgc%RootOUlmNutUptake_pvr(ids_NH4,1:NN,L1:L2,NZ)                  = 0.0_r8
  plt_rbgc%RootCUlmNutUptake_pvr(ids_NH4,1:NN,L1:L2,NZ)                  = 0.0_r8
  plt_rbgc%RootNH4DmndBand_pvr(1:NN,L1:L2,NZ)                            = 0.0_r8
  plt_rbgc%RootNutUptake_pvr(ids_NH4B,1:NN,L1:L2,NZ)                     = 0.0_r8
  plt_rbgc%RootOUlmNutUptake_pvr(ids_NH4B,1:NN,L1:L2,NZ)                 = 0.0_r8
  plt_rbgc%RootCUlmNutUptake_pvr(ids_NH4B,1:NN,L1:L2,NZ)                 = 0.0_r8
  plt_rbgc%RootNO3DmndSoil_pvr(1:NN,L1:L2,NZ)                            = 0.0_r8
  plt_rbgc%RootNutUptake_pvr(ids_NO3,1:NN,L1:L2,NZ)                      = 0.0_r8
  plt_rbgc%RootOUlmNutUptake_pvr(ids_NO3,1:NN,L1:L2,NZ)                  = 0.0_r8
  plt_rbgc%RootCUlmNutUptake_pvr(ids_NO3,1:NN,L1:L2,NZ)                  = 0.0_r8
  plt_rbgc%RootNO3DmndBand_pvr(1:NN,L1:L2,NZ)                            = 0.0_r8
  plt_rbgc%RootNutUptake_pvr(ids_NO3B,1:NN,L1:L2,NZ)                     = 0.0_r8
  plt_rbgc%RootOUlmNutUptake_pvr(ids_NO3B,1:NN,L1:L2,NZ)                 = 0.0_r8
  plt_rbgc%RootCUlmNutUptake_pvr(ids_NO3B,1:NN,L1:L2,NZ)                 = 0.0_r8
  plt_rbgc%RootH2PO4DmndSoil_pvr(1:NN,L1:L2,NZ)                          = 0.0_r8
  plt_rbgc%RootNutUptake_pvr(ids_H2PO4,1:NN,L1:L2,NZ)                    = 0.0_r8
  plt_rbgc%RootOUlmNutUptake_pvr(ids_H2PO4,1:NN,L1:L2,NZ)                = 0.0_r8
  plt_rbgc%RootCUlmNutUptake_pvr(ids_H2PO4,1:NN,L1:L2,NZ)                = 0.0_r8
  plt_rbgc%RootH2PO4DmndBand_pvr(1:NN,L1:L2,NZ)                          = 0.0_r8
  plt_rbgc%RootNutUptake_pvr(ids_H2PO4B,1:NN,L1:L2,NZ)                   = 0.0_r8
  plt_rbgc%RootOUlmNutUptake_pvr(ids_H2PO4B,1:NN,L1:L2,NZ)               = 0.0_r8
  plt_rbgc%RootCUlmNutUptake_pvr(ids_H2PO4B,1:NN,L1:L2,NZ)               = 0.0_r8
  plt_rbgc%RootH1PO4DmndSoil_pvr(1:NN,L1:L2,NZ)                          = 0.0_r8
  plt_rbgc%RootNutUptake_pvr(ids_H1PO4,1:NN,L1:L2,NZ)                    = 0.0_r8
  plt_rbgc%RootOUlmNutUptake_pvr(ids_H1PO4,1:NN,L1:L2,NZ)                = 0.0_r8
  plt_rbgc%RootCUlmNutUptake_pvr(ids_H1PO4,1:NN,L1:L2,NZ)                = 0.0_r8
  plt_rbgc%RootH1PO4DmndBand_pvr(1:NN,L1:L2,NZ)                          = 0.0_r8
  plt_rbgc%RootNutUptake_pvr(ids_H1PO4B,1:NN,L1:L2,NZ)                   = 0.0_r8
  plt_rbgc%RootOUlmNutUptake_pvr(ids_H1PO4B,1:NN,L1:L2,NZ)               = 0.0_r8
  plt_rbgc%RootCUlmNutUptake_pvr(ids_H1PO4B,1:NN,L1:L2,NZ)               = 0.0_r8
  plt_bgcr%RootN2Fix_pvr(L1:L2,NZ)                                       = 0.0_r8
  end subroutine ZeroNutrientUptake

!------------------------------------------------------------------------

  subroutine UptakeMineralPhosporhus(N,L,NZ,PathLen_pvr,FineRootRadius,FracPRoot4Uptake,MinFracPRoot4Uptake_pvr,&
    RootAreaDivRadius_vr,FCUP,FPUP,FWSRT,PerPlantRootH2OUptake)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  real(r8), intent(in):: PathLen_pvr(jroots,JZ1),FineRootRadius(jroots,JZ1),FracPRoot4Uptake(jroots,JZ1,JP1)
  real(r8), intent(in) :: MinFracPRoot4Uptake_pvr(jroots,JZ1,JP1)  !the minimum root fraction involved in substrate uptake
  real(r8), intent(in):: RootAreaDivRadius_vr(jroots,JZ1),FCUP,FPUP,FWSRT,PerPlantRootH2OUptake
  real(r8) :: TFPO4X,TFP14X,TFPOBX,TFP1BX
  real(r8) :: DIFFL
  real(r8) :: FP14X,FP1BX
  real(r8) :: FPO4X,FPOBX
  real(r8) :: POSGX
  real(r8) :: PATHL
  associate(                                                       &
    TortMicPM_vr             => plt_site%TortMicPM_vr,             &
    ZEROS                    => plt_site%ZEROS,                    &
    ZERO2                    => plt_site%ZERO2,                    &
    RootH1PO4DmndBand_pvr    => plt_rbgc%RootH1PO4DmndBand_pvr,    &
    RootH1PO4DmndSoil_pvr    => plt_rbgc%RootH1PO4DmndSoil_pvr,    &
    RootH2PO4DmndSoil_pvr    => plt_rbgc%RootH2PO4DmndSoil_pvr,    &
    RootH2PO4DmndBand_pvr    => plt_rbgc%RootH2PO4DmndBand_pvr,    &
    SoluteDifusvty_vr        => plt_soilchem%SoluteDifusvty_vr,    &
    RH1PO4EcoDmndBandPrev_vr => plt_bgcr%RH1PO4EcoDmndBandPrev_vr, &
    RH2PO4EcoDmndBandPrev_vr => plt_bgcr%RH2PO4EcoDmndBandPrev_vr, &
    RH1PO4EcoDmndSoilPrev_vr => plt_bgcr%RH1PO4EcoDmndSoilPrev_vr, &
    RH2PO4EcoDmndSoilPrev_vr => plt_bgcr%RH2PO4EcoDmndSoilPrev_vr  &
  )
  TFPO4X=0.0_r8
  TFPOBX=0.0_r8
  TFP14X=0.0_r8
  TFP1BX=0.0_r8

!     begin_execution
  IF(RH2PO4EcoDmndSoilPrev_vr(L).GT.ZEROS)THEN
    FPO4X=AMAX1(MinFracPRoot4Uptake_pvr(N,L,NZ),RootH2PO4DmndSoil_pvr(N,L,NZ)/RH2PO4EcoDmndSoilPrev_vr(L))
  ELSE
    FPO4X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RH2PO4EcoDmndBandPrev_vr(L).GT.ZEROS)THEN
    FPOBX=AMAX1(MinFracPRoot4Uptake_pvr(N,L,NZ),RootH2PO4DmndBand_pvr(N,L,NZ)/RH2PO4EcoDmndBandPrev_vr(L))
  ELSE
    FPOBX=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RH1PO4EcoDmndSoilPrev_vr(L).GT.ZEROS)THEN
    FP14X=AMAX1(MinFracPRoot4Uptake_pvr(N,L,NZ),RootH1PO4DmndSoil_pvr(N,L,NZ)/RH1PO4EcoDmndSoilPrev_vr(L))
  ELSE
    FP14X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RH1PO4EcoDmndBandPrev_vr(L).GT.ZEROS)THEN
    FP1BX=AMAX1(MinFracPRoot4Uptake_pvr(N,L,NZ),RootH1PO4DmndBand_pvr(N,L,NZ)/RH1PO4EcoDmndBandPrev_vr(L))
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
    !     TortMicPM_vr=soil tortuosity
    !     PathLen_pvr=path length of water and nutrient uptake
    !     FineRootRadius=root radius
    !     DIFFL=PO4 diffusion per plant
    !
    POSGX=SoluteDifusvty_vr(ids_H1PO4,L)*TortMicPM_vr(NPH,L)
    PATHL=AMIN1(PathLen_pvr(N,L),FineRootRadius(N,L)+SQRT(2.0*POSGX))
    DIFFL=POSGX*safe_adb(RootAreaDivRadius_vr(N,L),LOG(PATHL/FineRootRadius(N,L)))

    call UptakeH2PO4(N,L,NZ,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FWSRT,PerPlantRootH2OUptake)

    call UptakeHPO4(N,L,NZ,DIFFL,FP14X,FP1BX,FCUP,FPUP,FWSRT,PerPlantRootH2OUptake)

  ENDIF
  end associate
  end subroutine UptakeMineralPhosporhus
!------------------------------------------------------------------------

  subroutine UptakeNO3(I,J,N,L,NZ,FNO3X,FNOBX,PathLen_pvr,FineRootRadius,RootAreaDivRadius_vr,FCUP,FZUP,&
    FWSRT,PerPlantRootH2OUptake)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  real(r8), intent(in):: FNO3X,FNOBX
  real(r8), intent(in) :: PathLen_pvr(jroots,JZ1)
  real(r8), intent(in) :: FineRootRadius(jroots,JZ1)
  real(r8), intent(in) :: RootAreaDivRadius_vr(jroots,JZ1)
  real(r8), intent(in):: FCUP,FZUP,FWSRT,PerPlantRootH2OUptake

  real(r8) :: DIFFL
  real(r8) :: DIFNO3,DIFNOB
  real(r8) :: PATHL
  real(r8) :: RMFNO3,RTKNO3,RTKNOP,RMFNOB,RTKNOB,RTKNPB
  real(r8) :: UPMX,UPMXP
  real(r8) :: ZOSGX,ZNO3M,ZNO3X,ZNOBM,ZNOBX
  type(PlantSoluteUptakeConfig_type) :: PlantSoluteUptakeConfig
  logical :: ldebug

! begin_execution
  associate(                                                   &
    PlantPopulation_pft    => plt_site%PlantPopulation_pft,    &
    TortMicPM_vr           => plt_site%TortMicPM_vr,           &
    ZERO                   => plt_site%ZERO,                   &
    RootAreaPerPlant_pvr   => plt_morph%RootAreaPerPlant_pvr,  &
    fTgrowRootP_vr         => plt_pheno%fTgrowRootP_vr,        &
    RAutoRootO2Limter_rpvr => plt_rbgc%RAutoRootO2Limter_rpvr, &
    VmaxNO3Root_pft        => plt_rbgc%VmaxNO3Root_pft,        &
    CminNO3Root_pft        => plt_rbgc%CminNO3Root_pft,        &
    KmNO3Root_pft          => plt_rbgc%KmNO3Root_pft,          &
    RootCUlmNutUptake_pvr  => plt_rbgc%RootCUlmNutUptake_pvr,  &
    RootNO3DmndSoil_pvr    => plt_rbgc%RootNO3DmndSoil_pvr,    &
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr,      &
    RootNO3DmndBand_pvr    => plt_rbgc%RootNO3DmndBand_pvr,    &
    RootOUlmNutUptake_pvr  => plt_rbgc%RootOUlmNutUptake_pvr,  &
    SoluteDifusvty_vr      => plt_soilchem%SoluteDifusvty_vr,  &
    trcs_solml_vr          => plt_soilchem%trcs_solml_vr,      &
    trcs_VLN_vr            => plt_soilchem%trcs_VLN_vr,        &
    trc_solcl_vr           => plt_soilchem%trc_solcl_vr,       &
    VLWatMicP_vr           => plt_soilchem%VLWatMicP_vr        &
  )
  ldebug=.false.
!
! PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NO3
! FROM SOIL TO ROOT
!
! ZOSGL=NO3 diffusivity
! TortMicPM_vr=soil tortuosity
! FineRootRadius=root radius
! PATH=path length of water and nutrient uptake
! DIFFL=NO3 diffusion per plant
!
  ZOSGX=SoluteDifusvty_vr(ids_NO3,L)*TortMicPM_vr(NPH,L)
  PATHL=AMIN1(PathLen_pvr(N,L),FineRootRadius(N,L)+SQRT(2.0_r8*ZOSGX))
  DIFFL=ZOSGX*safe_adb(RootAreaDivRadius_vr(N,L),LOG(PATHL/FineRootRadius(N,L)))
  !
  ! NO3 UPTAKE IN NON-BAND SOIL ZONE
  !
  !  VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
  !     CNO3S=NO3 concentration in non-band
  !     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFNO3=soil-root convective NO3 flux per plant in non-band
  !     DIFNO3=soil-root NO3 diffusion per plant in non-band
  !
  IF(trcs_VLN_vr(ids_NO3,L).GT.ZERO.AND.trc_solcl_vr(ids_NO3,L).GT.CminNO3Root_pft(N,NZ))THEN
    RMFNO3 = PerPlantRootH2OUptake*trcs_VLN_vr(ids_NO3,L)
    DIFNO3 = DIFFL*trcs_VLN_vr(ids_NO3,L)
    !
    !     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=max NO3 uptake in non-band unlimited,limited by O2
    !     RootAreaPerPlant_pvr=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     fTgrowRootP_vr=temperature function for root growth
    !     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
    !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
    !
    UPMXP=VmaxNO3Root_pft(N,NZ)*RootAreaPerPlant_pvr(N,L,NZ) &
      *FWSRT*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_NO3,L)*AMIN1(FCUP,FZUP)

    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFNO3=soil-root convective N03 flux per plant in non-band
    !     DIFNO3=soil-root N03 diffusion per plant in non-band
    !     CNO3S=NO3 concentration in non-band
    !     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake,Km,min concn from PFT file
    !     RTKNO3,RTKNOP=NO3 uptake per plant in non-band lmtd,unlmtd by O2
    !     ZNO3M,ZNO3X=minimum,maximum NO3 available for uptake in non-band
    !     FNO3X=fraction of total NH4 uptake in non-band by root,myco populn
    !     RootNO3DmndSoil_pvr,RootNO3Uptake_pvr=NO3 uptake in non-band unlimited,limited by NO3
    !     RootOUlmNutUptake_pvr=NO3 uptake in non-band unlimited by O2
    !     RootCUlmNutUptake_pvr=NO3 uptake in non-band unlimited by nonstructural C
!
    ZNO3M=CminNO3Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_NO3,L)
    ZNO3X=AZMAX1(FNO3X*(trcs_solml_vr(ids_NO3,L)-ZNO3M))

    PlantSoluteUptakeConfig%SoluteConcMin   = CminNO3Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFNO3
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFNO3
    PlantSoluteUptakeConfig%UptakeRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_NO3,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = ZNO3X
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopulation_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmNO3Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootNO3DmndSoil_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_NO3,N,L,NZ),&
      RootCUlmNutUptake_pvr(ids_NO3,N,L,NZ),RootNutUptake_pvr(ids_NO3,N,L,NZ),ldebug)

  ENDIF
  !
  !     NO3 UPTAKE IN BAND SOIL ZONE
  !
  !     VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
  !     CNO3B=NO3 concentration in band
  !     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFNOB=soil-root convective NO3 flux per plant in band
  !     DIFNOB=soil-root NO3 diffusion per plant in band
  !

  IF(trcs_VLN_vr(ids_NO3B,L).GT.ZERO.AND.trc_solcl_vr(ids_NO3B,L).GT.CminNO3Root_pft(N,NZ))THEN
    RMFNOB = PerPlantRootH2OUptake*trcs_VLN_vr(ids_NO3B,L)
    DIFNOB = DIFFL*trcs_VLN_vr(ids_NO3B,L)
    !
    !     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum NO3 uptake in band unlimited,limited by O2
    !     RootAreaPerPlant_pvr=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     fTgrowRootP_vr=temperature function for root growth
    !     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
    !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
    !
    UPMXP=VmaxNO3Root_pft(N,NZ)*RootAreaPerPlant_pvr(N,L,NZ) &
      *FWSRT*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_NO3B,L)*AMIN1(FCUP,FZUP)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFNOB=soil-root convective NO3 flux per plant in band
    !     DIFNOB=soil-root NO3 diffusion per plant in band
    !     CNO3B=NH4 concentration in band
    !     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake,Km,min concn from PFT file
    !     RTKNOB,RTKNPB=NO3 uptake per plant in band lmtd,unlmtd by O2
    !     ZNOBM,ZNOBX=minimum,maximum NO3 available for uptake in band
    !     FNOBX=fraction of total NO3 uptake in band by root,myco populn
    !     RootNO3DmndBand_pvr,RootNO3BUptake_pvr=NO3 uptake in band unlimited,limited by NH4
    !     RootOUlmNutUptake_pvr=NO3 uptake in band unlimited by O2
    !     RootCUlmNutUptake_pvr=NO3 uptake in band unlimited by nonstructural C
    !
    ZNOBM = CminNO3Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_NO3B,L)
    ZNOBX = AZMAX1(FNOBX*(trcs_solml_vr(ids_NO3B,L)-ZNOBM))

    PlantSoluteUptakeConfig%SoluteConcMin   = CminNO3Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFNOB
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFNOB
    PlantSoluteUptakeConfig%UptakeRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_NO3B,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = ZNOBX
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopulation_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmNO3Root_pft(N,NZ)
  
    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootNO3DmndBand_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_NO3B,N,L,NZ),&
      RootCUlmNutUptake_pvr(ids_NO3B,N,L,NZ),RootNutUptake_pvr(ids_NO3B,N,L,NZ))

  ENDIF
  end associate
  end subroutine UptakeNO3
!------------------------------------------------------------------------
  subroutine UptakeNH4(N,L,NZ,FNH4X,FNHBX,PathLen_pvr,FineRootRadius,RootAreaDivRadius_vr,&
    FCUP,FZUP,FWSRT,PerPlantRootH2OUptake)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: FNH4X,FNHBX,PathLen_pvr(jroots,JZ1),FineRootRadius(jroots,JZ1)
  real(r8), intent(in) :: RootAreaDivRadius_vr(jroots,JZ1)
  real(r8), intent(in) :: FCUP,FZUP,FWSRT,PerPlantRootH2OUptake
  real(r8) :: DIFFL
  real(r8) :: DIFNH4,DIFNHB
  real(r8) :: PATHL
  real(r8) :: RMFNH4,RTKNH4,RTKNHP,RMFNHB,RTKNHB,RTKNBP
  real(r8) :: UPMX,UPMXP
  real(r8) :: ZNHBX,ZNSGX,ZNH4M,ZNH4X,ZNHBM
  type(PlantSoluteUptakeConfig_type) :: PlantSoluteUptakeConfig
! begin_execution
  associate(                                                   &
    PlantPopulation_pft    => plt_site%PlantPopulation_pft,    &
    ZERO                   => plt_site%ZERO,                   &
    TortMicPM_vr           => plt_site%TortMicPM_vr,           &
    fTgrowRootP_vr         => plt_pheno%fTgrowRootP_vr,        &
    RAutoRootO2Limter_rpvr => plt_rbgc%RAutoRootO2Limter_rpvr, &
    CMinNH4Root_pft        => plt_rbgc%CMinNH4Root_pft,        &
    VmaxNH4Root_pft        => plt_rbgc%VmaxNH4Root_pft,        &
    KmNH4Root_pft          => plt_rbgc%KmNH4Root_pft,          &
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr,      &
    RootCUlmNutUptake_pvr  => plt_rbgc%RootCUlmNutUptake_pvr,  &
    RootOUlmNutUptake_pvr  => plt_rbgc%RootOUlmNutUptake_pvr,  &
    RootNH4DmndSoil_pvr    => plt_rbgc%RootNH4DmndSoil_pvr,    &
    RootNH4DmndBand_pvr    => plt_rbgc%RootNH4DmndBand_pvr,    &
    RootAreaPerPlant_pvr   => plt_morph%RootAreaPerPlant_pvr,  &
    SoluteDifusvty_vr      => plt_soilchem%SoluteDifusvty_vr,  &
    VLWatMicP_vr           => plt_soilchem%VLWatMicP_vr,       &
    trcs_VLN_vr            => plt_soilchem%trcs_VLN_vr,        &
    trcs_solml_vr           => plt_soilchem%trcs_solml_vr,       &
    trc_solcl_vr           => plt_soilchem%trc_solcl_vr        &
  )
! ZNSGL=NH4 diffusivity
! TortMicPM_vr=soil tortuosity
! PATH=path length of water and nutrient uptake
! FineRootRadius=root radius
! DIFFL=NH4 diffusion per plant
!
  ZNSGX = SoluteDifusvty_vr(idg_NH3,L)*TortMicPM_vr(NPH,L)
  PATHL = AMIN1(PathLen_pvr(N,L),FineRootRadius(N,L)+SQRT(2.0*ZNSGX))
  DIFFL = ZNSGX*safe_adb(RootAreaDivRadius_vr(N,L),LOG(PATHL/FineRootRadius(N,L)))
!
! NH4 UPTAKE IN NON-BAND SOIL ZONE
!
! VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
! CNH4S=NH4 concentration in non-band
! VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake,Km,min concn from PFT file
! PerPlantRootH2OUptake=root water uptake per plant
! RMFNH4=soil-root convective NH4 flux per plant in non-band
! DIFNH4=soil-root NH4 diffusion per plant in non-band
!
  IF(trcs_VLN_vr(ids_NH4,L).GT.ZERO.AND.trc_solcl_vr(ids_NH4,L).GT.CMinNH4Root_pft(N,NZ))THEN
    RMFNH4=PerPlantRootH2OUptake*trcs_VLN_vr(ids_NH4,L)
    DIFNH4=DIFFL*trcs_VLN_vr(ids_NH4,L)
!
!   NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
!   AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
!
!   UPMXP,UPMX=max NH4 uptake in non-band unlimited,limited by O2
!   RootAreaPerPlant_pvr=root surface area per plant from grosub.f
!   FWSRT=protein concentration relative to 5%
!   fTgrowRootP_vr=temperature function for root growth
!   FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
!   RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
!
    UPMXP=VmaxNH4Root_pft(N,NZ)*RootAreaPerPlant_pvr(N,L,NZ) &
      *FWSRT*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_NH4,L)*AMIN1(FCUP,FZUP)
!
!   SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
!   SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
!   WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!   RMFNH4=soil-root convective NH4 flux per plant in non-band
!   DIFNH4=soil-root NH4 diffusion per plant in non-band
!   CNH4S=NH4 concentration in non-band
!   VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake,Km,min concn from PFT file
!   RTKNH4,RTKNHP=NH4 uptake per plant in non-band lmtd,unlmtd by O2
!   ZNH4M,ZNH4X=minimum,maximum NH4 available for uptake in non-band
!   FNH4X=fraction of total NH4 uptake in non-band by root,myco populn
!   RootNH4DmndSoil_pvr,RootNH4Uptake_pvr=NH4 uptake in non-band unlimited,limited by NH4
!   RootOUlmNutUptake_pvr=NH4 uptake in non-band unlimited by O2
!   RootCUlmNutUptake_pvr=NH4 uptake in non-band unlimited by nonstructural C
!
    ZNH4M = CMinNH4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_NH4,L)
    ZNH4X = AZMAX1(FNH4X*(trcs_solml_vr(ids_NH4,L)-ZNH4M))

    PlantSoluteUptakeConfig%SoluteConcMin   = CMinNH4Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFNH4
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFNH4
    PlantSoluteUptakeConfig%UptakeRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_NH4,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = ZNH4X
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopulation_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmNH4Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootNH4DmndSoil_pvr(N,L,NZ),&
      RootOUlmNutUptake_pvr(ids_NH4,N,L,NZ),RootCUlmNutUptake_pvr(ids_NH4,N,L,NZ),&
      RootNutUptake_pvr(ids_NH4,N,L,NZ))

  ENDIF
!
! NH4 UPTAKE IN BAND SOIL ZONE
!
! VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
! CNH4B=NH4 concentration in band
! VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake,Km,min concn from PFT file
! PerPlantRootH2OUptake=root water uptake per plant
! RMFNHB=soil-root convective NH4 flux per plant in band
! DIFNHB=soil-root NH4 diffusion per plant in band
!

  IF(trcs_VLN_vr(ids_NH4B,L).GT.ZERO.AND.trc_solcl_vr(ids_NH4B,L).GT.CMinNH4Root_pft(N,NZ))THEN
    RMFNHB=PerPlantRootH2OUptake*trcs_VLN_vr(ids_NH4B,L)
    DIFNHB=DIFFL*trcs_VLN_vr(ids_NH4B,L)
!
!   NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
!   AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
!
!   UPMXP,UPMX=maximum NH4 uptake in band unlimited,limited by O2
!   RootAreaPerPlant_pvr=root surface area per plant from grosub.f
!   FWSRT=protein concentration relative to 5%
!   fTgrowRootP_vr=temperature function for root growth
!   FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
!   RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
!
    UPMXP=VmaxNH4Root_pft(N,NZ)*RootAreaPerPlant_pvr(N,L,NZ) &
      *FWSRT*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_NH4B,L)*AMIN1(FCUP,FZUP)
!
!   SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
!   SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
!   WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!   RMFNHB=soil-root convective NH4 flux per plant in band
!   DIFNHB=soil-root NH4 diffusion per plant in band
!   CNH4B=NH4 concentration in band
!   VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake,Km,min concn from PFT file
!   RTKNHB,RTKNBP=NH4 uptake per plant in band lmtd,unlmtd by O2
!   ZNHBM,ZNHBX=minimum,maximum NH4 available for uptake in band
!   FNHBX=fraction of total NH4 uptake in band by root,myco populn
!   RootNH4DmndBand_pvr,RootNH4BUptake_pvr=NH4 uptake in band unlimited,limited by NH4
!   RootOUlmNutUptake_pvr=NH4 uptake in band unlimited by O2
!   RootCUlmNutUptake_pvr=NH4 uptake in band unlimited by nonstructural C
!

    ZNHBM=CMinNH4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_NH4B,L)
    ZNHBX=AZMAX1(FNHBX*(trcs_solml_vr(ids_NH4B,L)-ZNHBM))

    PlantSoluteUptakeConfig%SoluteConcMin   = CMinNH4Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFNHB
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFNHB
    PlantSoluteUptakeConfig%UptakeRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_NH4B,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = ZNHBX
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopulation_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmNH4Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootNH4DmndBand_pvr(N,L,NZ),&
      RootOUlmNutUptake_pvr(ids_NH4B,N,L,NZ),RootCUlmNutUptake_pvr(ids_NH4B,N,L,NZ),&
      RootNutUptake_pvr(ids_NH4B,N,L,NZ))

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
  associate(                                                   &
    PlantPopulation_pft    => plt_site%PlantPopulation_pft,    &
    ZERO                   => plt_site%ZERO,                   &
    fTgrowRootP_vr         => plt_pheno%fTgrowRootP_vr,        &
    RAutoRootO2Limter_rpvr => plt_rbgc%RAutoRootO2Limter_rpvr, &
    CMinPO4Root_pft        => plt_rbgc%CMinPO4Root_pft,        &
    KmPO4Root_pft          => plt_rbgc%KmPO4Root_pft,          &
    VmaxPO4Root_pft        => plt_rbgc%VmaxPO4Root_pft,        &
    RootCUlmNutUptake_pvr  => plt_rbgc%RootCUlmNutUptake_pvr,  &
    RootH1PO4DmndSoil_pvr  => plt_rbgc%RootH1PO4DmndSoil_pvr,  &
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr,      &
    RootOUlmNutUptake_pvr  => plt_rbgc%RootOUlmNutUptake_pvr,  &
    RootH1PO4DmndBand_pvr  => plt_rbgc%RootH1PO4DmndBand_pvr,  &
    RootAreaPerPlant_pvr   => plt_morph%RootAreaPerPlant_pvr,  &
    VLWatMicP_vr           => plt_soilchem%VLWatMicP_vr,       &
    trcs_VLN_vr            => plt_soilchem%trcs_VLN_vr,        &
    trc_solcl_vr           => plt_soilchem%trc_solcl_vr,       &
    trcs_solml_vr          => plt_soilchem%trcs_solml_vr        &
  )
  !
  !     HPO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH1P4=HPO4 concentration in non-band
  !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant > 0._r8
  !     RMFH1P=soil-root convective HPO4 flux per plant in non-band
  !     DIFH1P=soil-root HPO4 diffusion per plant in non-band > 0._r8
!
    IF(trcs_VLN_vr(ids_H1PO4,L).GT.ZERO.AND.trc_solcl_vr(ids_H1PO4,L).GT.CMinPO4Root_pft(N,NZ))THEN
      RMFH1P = PerPlantRootH2OUptake*trcs_VLN_vr(ids_H1PO4,L)
      DIFH1P = DIFFL*trcs_VLN_vr(ids_H1PO4,L)
    !
    !     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=max HPO4 uptake in non-band unlimited,limited by O2
    !     RootAreaPerPlant_pvr=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     fTgrowRootP_vr=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
    !
    UPMXP=0.1_r8*VmaxPO4Root_pft(N,NZ)*RootAreaPerPlant_pvr(N,L,NZ) &
      *FWSRT*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_H1PO4,L)*AMIN1(FCUP,FPUP)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH1P=soil-root convective HPO4 flux per plant in non-band
    !     DIFH1P=soil-root HPO4 diffusion per plant in non-band
    !     CH1P4=HPO4 concentration in non-band
    !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH1P,RTKHP1=HPO4 uptake per plant in non-band lmtd,unlmtd by O2
    !     H1POM,H1POX=minimum,maximum HPO4 available for uptake in non-band
    !     FP14X=fraction of total HPO4 uptake in non-band by root,myco populn
    !     RootH1PO4DmndSoil_pvr,RootNutUptake_pvr=HPO4 uptake in non-band unlimited,limited by HPO4
    !     RootOUlmNutUptake_pvr(ids_H1PO4,=HPO4 uptake in non-band unlimited by O2
    !     RootCUlmNutUptake_pvr=HPO4 uptake in non-band unlimited by nonstructural C

    H1POM=CMinPO4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_H1PO4,L)
    H1POX=AZMAX1(FP14X*(trcs_solml_vr(ids_H1PO4,L)-H1POM))

    PlantSoluteUptakeConfig%SoluteConcMin   = CMinPO4Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFH1P
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFH1P
    PlantSoluteUptakeConfig%UptakeRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_H1PO4,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = H1POX
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopulation_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmPO4Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootH1PO4DmndSoil_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_H1PO4,N,L,NZ),&
      RootCUlmNutUptake_pvr(ids_H1PO4,N,L,NZ),RootNutUptake_pvr(ids_H1PO4,N,L,NZ))

  ENDIF
  !
  !     HPO4 UPTAKE IN BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH1P4B=HPO4 concentration in band
  !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFH1B=soil-root convective HPO4 flux per plant in band
  !     DIFH1B=soil-root HPO4 diffusion per plant in band
  !
  IF(trcs_VLN_vr(ids_H1PO4B,L).GT.ZERO.AND.trc_solcl_vr(ids_H1PO4B,L).GT.CMinPO4Root_pft(N,NZ))THEN
    RMFH2B = PerPlantRootH2OUptake*trcs_VLN_vr(ids_H1PO4B,L)
    DIFH1B = DIFFL*trcs_VLN_vr(ids_H1PO4B,L)
    !
    !     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum HPO4 uptake in band unlimited,limited by O2
    !     RootAreaPerPlant_pvr=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     fTgrowRootP_vr=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
    !
    UPMXP=0.1_r8*VmaxPO4Root_pft(N,NZ)*RootAreaPerPlant_pvr(N,L,NZ) &
      *FWSRT*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_H1PO4B,L)*AMIN1(FCUP,FPUP)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH1B=soil-root convective HPO4 flux per plant in band
    !     DIFH1B=soil-root HPO4 diffusion per plant in band
    !     CH1P4B=HPO4 concentration in band
    !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH1B,RTKHB1=HPO4 uptake per plant in band lmtd,unlmtd by O2
    !     H1PXM,H1PXB=minimum,maximum HPO4 available for uptake in band
    !     FP1BX=fraction of total HPO4 uptake in band by root,myco populn
    !     RootH1PO4DmndBand_pvr,RootH1PO4BUptake_pvr=HPO4 uptake in band unlimited,limited by H2PO4
    !     RootOUlmNutUptake_pvr=HPO4 uptake in band unlimited by O2
    !     RootCUlmNutUptake_pvr=HPO4 uptake in band unlimited by nonstructural C

    H1PXM = CMinPO4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_H1PO4B,L)
    H1PXB = AZMAX1(FP1BX*(trcs_solml_vr(ids_H1PO4B,L)-H1PXM))

    PlantSoluteUptakeConfig%SoluteConcMin   = CMinPO4Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFH2B
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFH1B
    PlantSoluteUptakeConfig%UptakeRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_H1PO4B,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = H1PXB
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopulation_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmPO4Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootH1PO4DmndBand_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ),&
      RootCUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ),RootNutUptake_pvr(ids_H1PO4B,N,L,NZ))

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
  associate(                                                   &
    PlantPopulation_pft    => plt_site%PlantPopulation_pft,    &
    ZERO                   => plt_site%ZERO,                   &
    RAutoRootO2Limter_rpvr => plt_rbgc%RAutoRootO2Limter_rpvr, &
    KmPO4Root_pft          => plt_rbgc%KmPO4Root_pft,          &
    VmaxPO4Root_pft        => plt_rbgc%VmaxPO4Root_pft,        &
    CMinPO4Root_pft        => plt_rbgc%CMinPO4Root_pft,        &
    RootCUlmNutUptake_pvr  => plt_rbgc%RootCUlmNutUptake_pvr,  &
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr,      &
    RootOUlmNutUptake_pvr  => plt_rbgc%RootOUlmNutUptake_pvr,  &
    RootH2PO4DmndBand_pvr  => plt_rbgc%RootH2PO4DmndBand_pvr,  &
    RootH2PO4DmndSoil_pvr  => plt_rbgc%RootH2PO4DmndSoil_pvr,  &
    fTgrowRootP_vr         => plt_pheno%fTgrowRootP_vr,        &
    RootAreaPerPlant_pvr   => plt_morph%RootAreaPerPlant_pvr,  &
    trcs_solml_vr           => plt_soilchem%trcs_solml_vr,       &
    VLWatMicP_vr           => plt_soilchem%VLWatMicP_vr,       &
    trc_solcl_vr           => plt_soilchem%trc_solcl_vr,       &
    trcs_VLN_vr            => plt_soilchem%trcs_VLN_vr         &
  )
  !     H2PO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH2P4=H2PO4 concentration in non-band
  !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
  !     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
!
    IF(trcs_VLN_vr(ids_H1PO4,L).GT.ZERO.AND.trc_solcl_vr(ids_H2PO4,L).GT.CMinPO4Root_pft(N,NZ))THEN
      RMFH2P=PerPlantRootH2OUptake*trcs_VLN_vr(ids_H1PO4,L)
      DIFH2P=DIFFL*trcs_VLN_vr(ids_H1PO4,L)
      !
      !     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
      !     AND FROM ROOT SURFACE AREA, C AND P CONSTRAINTS CALCULATED ABOVE
      !
      !     UPMXP,UPMX=max H2PO4 uptake in non-band unlimited,limited by O2
      !     RootAreaPerPlant_pvr=root surface area per plant from grosub.f
      !     FWSRT=protein concentration relative to 5%
      !     fTgrowRootP_vr=temperature function for root growth
      !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
      !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
!
      UPMXP=VmaxPO4Root_pft(N,NZ)*RootAreaPerPlant_pvr(N,L,NZ) &
        *FWSRT*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_H1PO4,L)*AMIN1(FCUP,FPUP)
      !
      !     SOLUTION FOR MASS FLOW + DIFFUSION OF H2PO4 IN AQUEOUS PHASE OF
      !     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY
      !     COMPETITION WITH OTHER ROOT AND MICROBIAL POPULATIONS
      !
      !     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
      !     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
      !     CH2P4=H2PO4 concentration in non-band
      !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
      !     RTKH2P,RTKHPP=H2PO4 uptake per plant in non-band lmtd,unlmtd by O2
      !     H2POM,H2POX=minimum,maximum H2PO4 available for uptake in non-band
      !     FPO4X=fraction of total H2PO4 uptake in non-band by root,myco populn
      !     RootH2PO4DmndSoil_pvr,RootH2PO4Uptake_pvr=H2PO4 uptake in non-band unlimited,limited by H2PO4
      !     RootOUlmNutUptake_pvr=H2PO4 uptake in non-band unlimited by O2
      !     RootCUlmNutUptake_pvr=H2PO4 uptake in non-band unlimited by nonstructural C

      H2POM = CMinPO4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_H1PO4,L)
      H2POX = AZMAX1(FPO4X*(trcs_solml_vr(ids_H2PO4,L)-H2POM))

      PlantSoluteUptakeConfig%SoluteConcMin   = CMinPO4Root_pft(N,NZ)
      PlantSoluteUptakeConfig%SolAdvFlx       = RMFH2P
      PlantSoluteUptakeConfig%SolDifusFlx     = DIFH2P
      PlantSoluteUptakeConfig%UptakeRateMax   = UPMXP
      PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
      PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_H2PO4,L)
      PlantSoluteUptakeConfig%SoluteMassMax   = H2POX
      PlantSoluteUptakeConfig%CAvailStress    = FCUP
      PlantSoluteUptakeConfig%PlantPopulation = PlantPopulation_pft(NZ)
      PlantSoluteUptakeConfig%SoluteKM        = KmPO4Root_pft(N,NZ)

      call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootH2PO4DmndSoil_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_H2PO4,N,L,NZ),&
        RootCUlmNutUptake_pvr(ids_H2PO4,N,L,NZ),RootNutUptake_pvr(ids_H2PO4,N,L,NZ))

    ENDIF
    !
    !     H2PO4 UPTAKE IN BAND SOIL ZONE
    !
    !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
    !     CH2P4B=H2PO4 concentration in band
    !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
    !     PerPlantRootH2OUptake=root water uptake per plant
    !     RMFH2B=soil-root convective H2PO4 flux per plant in band
    !     DIFH2B=soil-root H2PO4 diffusion per plant in band
    !

  IF(trcs_VLN_vr(ids_H1PO4B,L).GT.ZERO.AND.trc_solcl_vr(ids_H2PO4B,L).GT.CMinPO4Root_pft(N,NZ))THEN
    RMFH2B=PerPlantRootH2OUptake*trcs_VLN_vr(ids_H1PO4B,L)
    DIFH2B=DIFFL*trcs_VLN_vr(ids_H1PO4B,L)
    !
    !     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum H2PO4 uptake in band unlimited,limited by O2
    !     RootAreaPerPlant_pvr=root surface area per plant from grosub.f
    !     FWSRT=protein concentration relative to 5%
    !     fTgrowRootP_vr=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
    !
    UPMXP=VmaxPO4Root_pft(N,NZ)*RootAreaPerPlant_pvr(N,L,NZ) &
      *FWSRT*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_H1PO4B,L)*AMIN1(FCUP,FPUP)

    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF PO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH2B=soil-root convective H2PO4 flux per plant in band
    !     DIFH2B=soil-root H2PO4 diffusion per plant in band
    !     CH2P4B=H2PO4 concentration in band
    !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH2B,RTKHPB=H2PO4 uptake per plant in band lmtd,unlmtd by O2
    !     H2PXM,H2PXB=minimum,maximum H2PO4 available for uptake in band
    !     FPOBX=fraction of total H2PO4 uptake in band by root,myco populn
    !     RootH2PO4DmndBand_pvr,RootNutUptake_pvr=H2PO4 uptake in band unlimited,limited by H2PO4
    !     RootOUlmNutUptake_pvr=H2PO4 uptake in band unlimited by O2
    !     RootCUlmNutUptake_pvr=H2PO4 uptake in band unlimited by nonstructural C

    H2PXM=CMinPO4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_H1PO4B,L)
    H2PXB=AZMAX1(FPOBX*(trcs_solml_vr(ids_H2PO4B,L)-H2PXM))

    PlantSoluteUptakeConfig%SoluteConcMin   = CMinPO4Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFH2B
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFH2B
    PlantSoluteUptakeConfig%UptakeRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_H2PO4B,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = H2PXB
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopulation_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmPO4Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootH2PO4DmndBand_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ),&
      RootCUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ),RootNutUptake_pvr(ids_H2PO4B,N,L,NZ))

  ENDIF
  end associate
  end subroutine UptakeH2PO4

!------------------------------------------------------------------------

  subroutine UptakeMineralNitrogen(I,J,N,L,NZ,PathLen_pvr,FineRootRadius,FracPRoot4Uptake,&
    MinFracPRoot4Uptake_pvr,RootAreaDivRadius_vr,FCUP,FZUP,FWSRT,PerPlantRootH2OUptake)

  implicit none
  integer , intent(in) :: I,J
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: PathLen_pvr(jroots,JZ1),FineRootRadius(jroots,JZ1),FracPRoot4Uptake(jroots,JZ1,JP1)
  real(r8), intent(in) :: MinFracPRoot4Uptake_pvr(jroots,JZ1,JP1)
  real(r8), intent(in) :: RootAreaDivRadius_vr(jroots,JZ1)
  real(r8), intent(in) :: FCUP,FZUP,FWSRT,PerPlantRootH2OUptake
  real(r8) :: FNO3X,FNOBX,FNH4X,FNHBX
  real(r8) :: TFNH4X,TFNO3X,TFNHBX,TFNOBX

!     begin_execution
  associate(                                                   &
    ZERO2                  => plt_site%ZERO2,                  &
    ZEROS                  => plt_site%ZEROS,                  &
    RootNO3DmndBand_pvr    => plt_rbgc%RootNO3DmndBand_pvr,    &
    RootNH4DmndSoil_pvr    => plt_rbgc%RootNH4DmndSoil_pvr,    &
    RootNH4DmndBand_pvr    => plt_rbgc%RootNH4DmndBand_pvr,    &
    RootNO3DmndSoil_pvr    => plt_rbgc%RootNO3DmndSoil_pvr,    &
    RNO3EcoDmndSoilPrev_vr => plt_bgcr%RNO3EcoDmndSoilPrev_vr, &
    RNH4EcoDmndBandPrev_vr => plt_bgcr%RNH4EcoDmndBandPrev_vr, &
    RNH4EcoDmndSoilPrev_vr => plt_bgcr%RNH4EcoDmndSoilPrev_vr, &
    RNO3EcoDmndBandPrev_vr => plt_bgcr%RNO3EcoDmndBandPrev_vr  &
  )
  TFNH4X=0.0_r8
  TFNHBX=0.0_r8
  TFNO3X=0.0_r8
  TFNOBX=0.0_r8

  IF(RNH4EcoDmndSoilPrev_vr(L).GT.ZEROS)THEN
    FNH4X=AMAX1(MinFracPRoot4Uptake_pvr(N,L,NZ),RootNH4DmndSoil_pvr(N,L,NZ)/RNH4EcoDmndSoilPrev_vr(L))
  ELSE
    FNH4X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RNH4EcoDmndBandPrev_vr(L).GT.ZEROS)THEN
    FNHBX=AMAX1(MinFracPRoot4Uptake_pvr(N,L,NZ),RootNH4DmndBand_pvr(N,L,NZ)/RNH4EcoDmndBandPrev_vr(L))
  ELSE
    FNHBX=FracPRoot4Uptake(N,L,NZ)
  ENDIF

  IF(RNO3EcoDmndSoilPrev_vr(L).GT.ZEROS)THEN
    FNO3X=AMAX1(MinFracPRoot4Uptake_pvr(N,L,NZ),RootNO3DmndSoil_pvr(N,L,NZ)/RNO3EcoDmndSoilPrev_vr(L))
  ELSE
    FNO3X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RNO3EcoDmndBandPrev_vr(L).GT.ZEROS)THEN
    FNOBX=AMAX1(MinFracPRoot4Uptake_pvr(N,L,NZ),RootNO3DmndBand_pvr(N,L,NZ)/RNO3EcoDmndBandPrev_vr(L))
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

    call UptakeNH4(N,L,NZ,FNH4X,FNHBX,PathLen_pvr,FineRootRadius,RootAreaDivRadius_vr,FCUP,FZUP,FWSRT,&
      PerPlantRootH2OUptake)

    call UptakeNO3(I,J,N,L,NZ,FNO3X,FNOBX,PathLen_pvr,FineRootRadius,RootAreaDivRadius_vr,FCUP,FZUP,FWSRT,&
      PerPlantRootH2OUptake)

  ENDIF
  end associate
  end subroutine UptakeMineralNitrogen

!------------------------------------------------------------------------

  subroutine GetUptakeCapcity(N,L,NZ,FracPRoot4Uptake,MinFracPRoot4Uptake_pvr,FCUP,FZUP,FPUP,&
    FWSRT,PerPlantRootH2OUptake,dtPerPlantRootH2OUptake,FOXYX)

  implicit none
  integer, intent(in)   :: N,L
  integer, intent(in)   :: NZ
  REAL(R8), INTENT(IN)  :: FracPRoot4Uptake(jroots,JZ1,JP1)
  real(r8), intent(in)  :: MinFracPRoot4Uptake_pvr(jroots,JZ1,JP1)
  real(r8), intent(out) :: FCUP,FZUP,FPUP,FWSRT
  real(r8), intent(out) :: PerPlantRootH2OUptake,dtPerPlantRootH2OUptake
  real(r8), intent(out) :: FOXYX

  character(len=*), parameter :: subname='GetUptakeCapcity'
  real(r8) :: dss
  associate(                                                          &
    RO2EcoDmndPrev_vr         => plt_bgcr%RO2EcoDmndPrev_vr,          &
    RootCO2EmisPot_pvr        => plt_rbgc%RootCO2EmisPot_pvr,         &
    RootO2Dmnd4Resp_pvr       => plt_rbgc%RootO2Dmnd4Resp_pvr,        &
    PlantPopulation_pft       => plt_site%PlantPopulation_pft,        &
    ZEROS                     => plt_site%ZEROS,                      &
    ZERO                      => plt_site%ZERO,                       &
    AllPlantRootH2OLoss_vr    => plt_ew%AllPlantRootH2OLoss_vr,       &
    RootFracRemobilizableBiom => plt_allom%RootFracRemobilizableBiom, &
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft,             &
    RootProteinConc_rpvr      => plt_biom%RootProteinConc_rpvr,       &
    RootProteinC_pvr          => plt_biom%RootProteinC_pvr,           &
    RootMycoNonstElms_rpvr    => plt_biom%RootMycoNonstElms_rpvr,     &
    RootNonstructElmConc_rpvr => plt_biom%RootNonstructElmConc_rpvr,  &
    RootMycoActiveBiomC_pvr   => plt_biom%RootMycoActiveBiomC_pvr     &
  )
  !
  !     UPTAKE CAPACITY 'FWSRT' DEPENDS ON ROOT,MYCORRHIZAL
  !     PROTEIN CONTENT RELATIVE TO 5% FOR WHICH ACTIVE UPTAKE
  !     PARAMETERS ARE DEFINED
  !
  !     RootProteinConc_rpvr,RootFracRemobilizableBiom=current,maximum protein concentration
  !     RootProteinC_pvr,WTRTL=protein content,mass
  !     FWSRT=protein concentration relative to 5%
  !
  call PrintInfo('beg '//subname)
  IF(RootMycoActiveBiomC_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
    RootProteinConc_rpvr(N,L,NZ)=AMIN1(RootFracRemobilizableBiom(NZ),RootProteinC_pvr(N,L,NZ)/RootMycoActiveBiomC_pvr(N,L,NZ))
    FWSRT=RootProteinConc_rpvr(N,L,NZ)/0.05_r8
  ELSE
    RootProteinConc_rpvr(N,L,NZ)=RootFracRemobilizableBiom(NZ)
    FWSRT=1.0_r8
  ENDIF
  !
  !     RESPIRATION CONSTRAINT ON UPTAKE FROM NON-STRUCTURAL C
  !
  !     RootCO2EmisPot_pvr=total respiration from CPOOLR
  !     FCUP=limitation to active uptake respiration from CPOOLR
  !     CPOOLR=nonstructural C content
  !
  IF(RootCO2EmisPot_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
    FCUP=AZMAX1(AMIN1(1.0_r8,0.25_r8*safe_adb(RootMycoNonstElms_rpvr(ielmc,N,L,NZ),RootCO2EmisPot_pvr(N,L,NZ))))
  ELSE
    FCUP=0.0_r8
  ENDIF
  !
  !     FEEDBACK CONSTRAINT ON N UPTAKE FROM NON-STRUCTURAL N AND P
  !
  !     FZUP,FPUP=limitn to active uptake respiration from CZPOLR,CPPOLR
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition effect on N,P uptake
  !     dtPerPlantRootH2OUptake=water uptake at time step for gas flux calculations
  !
  IF(RootNonstructElmConc_rpvr(ielmc,N,L,NZ).GT.ZERO)THEN
    FZUP=AMIN1(safe_adb(RootNonstructElmConc_rpvr(ielmc,N,L,NZ),RootNonstructElmConc_rpvr(ielmc,N,L,NZ) &
      +RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/ZCKI) &
      ,safe_adb(RootNonstructElmConc_rpvr(ielmp,N,L,NZ),RootNonstructElmConc_rpvr(ielmp,N,L,NZ) &
       +RootNonstructElmConc_rpvr(ielmn,N,L,NZ)/ZPKI))
    FPUP=AMIN1(safe_adb(RootNonstructElmConc_rpvr(ielmc,N,L,NZ),RootNonstructElmConc_rpvr(ielmc,N,L,NZ) &
      +RootNonstructElmConc_rpvr(ielmp,N,L,NZ)/PCKI) &
      ,safe_adb(RootNonstructElmConc_rpvr(ielmn,N,L,NZ),RootNonstructElmConc_rpvr(ielmn,N,L,NZ) &
      +RootNonstructElmConc_rpvr(ielmp,N,L,NZ)/PZKI))
  ELSE
    FZUP=0.0_r8
    FPUP=0.0_r8
  ENDIF
  !NN=0
  PerPlantRootH2OUptake=AZMAX1(-AllPlantRootH2OLoss_vr(N,L,NZ)/PlantPopulation_pft(NZ))
  dtPerPlantRootH2OUptake=PerPlantRootH2OUptake*dts_gas
  !
  !     FACTORS CONSTRAINING O2 AND NUTRIENT UPTAKE AMONG
  !     COMPETING ROOT,MYCORRHIZAL AND MICROBIAL POPULATIONS
  !     IN BAND AND NON-BAND SOIL ZONES FROM DEMAND CALCULATED
  !     IN PREVIOUS HOUR
  !
  !
  IF(RO2EcoDmndPrev_vr(L).GT.ZEROS)THEN
    FOXYX=AMAX1(MinFracPRoot4Uptake_pvr(N,L,NZ),RootO2Dmnd4Resp_pvr(N,L,NZ)/RO2EcoDmndPrev_vr(L))
  ELSE
    FOXYX=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine GetUptakeCapcity
!------------------------------------------------------------------------

  subroutine RootExudates(I,J,N,L,NZ)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ

  real(r8) :: CPOOLX,CPOOLT
  real(r8) :: PPOOLX,ZPOOLX
  real(r8) :: VLWatMicPK,VLWatMicPT
  real(r8) :: XFRE(NumPlantChemElms)
  real(r8) :: RootExudE,scal
  integer :: K,NE
  !     begin_execution
  associate(                                                   &
    RootMycoNonstElms_rpvr => plt_biom%RootMycoNonstElms_rpvr, &
    ZERO4Groth_pft         => plt_biom%ZERO4Groth_pft,         &
    ZEROS                  => plt_site%ZEROS,                  &
    ZEROS2                 => plt_site%ZEROS2,                 &
    VLWatMicPM_vr          => plt_site%VLWatMicPM_vr,          &
    RootVH2O_pvr           => plt_morph%RootVH2O_pvr,          &
    RootMycoExudEUptk_pvr  => plt_rbgc%RootMycoExudEUptk_pvr,  &
    FracBulkSOMC_vr        => plt_soilchem%FracBulkSOMC_vr,    &
    DOM_vr                 => plt_soilchem%DOM_vr              &
  )
  !
  !     ROOT EXUDATION OF C, N AND P DEPENDS ON CONCN DIFFERENCES
  !     BETWEEN ROOT NON-STRUCTURAL POOLS AND SOIL DISSOLVED POOLS
  !
  !     VLWatMicPMM=soil micropore water volume
  !     FracBulkSOMC_vr=fraction of total SOC in each substrate K from nitro.f
  !     RootVH2O_pvr=root aqueous volume
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P in root,myco
  !     XFRC,XFRN,XFRP=nonstructural C,N,P exchg at root-soil DOC equilibrium
  !     OQC=soil DOC
  !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
  !     FEXUC,FEXUN,FEXUP=rate constant for root C,N,P exudation
  !     Air_Heat_Latent_store_col,Air_Heat_Sens_store_col=total fluxes x blr for calculating canopy air temperature,
  !     vapor pressure in watsub.f
  !      EvapTransLHeat_pft,SFLXC=canopylatent,sensible heat fluxes
  !      RA=canopy boundary layer resistance
  !     OSTR=O2 stress indicator
  !
  D195: DO K=1,jcplx
    VLWatMicPK=VLWatMicPM_vr(NPH,L)*FracBulkSOMC_vr(K,L)
    IF(VLWatMicPK.GT.ZEROS2.AND.RootVH2O_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
      VLWatMicPT=VLWatMicPK+RootVH2O_pvr(N,L,NZ)
      CPOOLX      = AMIN1(1.25E+03_r8*RootVH2O_pvr(N,L,NZ),RootMycoNonstElms_rpvr(ielmc,N,L,NZ))
      XFRE(ielmc) = (DOM_vr(idom_doc,K,L)*RootVH2O_pvr(N,L,NZ)-CPOOLX*VLWatMicPK)/VLWatMicPT
      !XFRC, >0 positive into plants, 
      RootMycoExudEUptk_pvr(ielmc,N,K,L,NZ)=FEXUDE(ielmc)*XFRE(ielmc)
      IF(DOM_vr(idom_doc,K,L).GT.ZEROS .AND. RootMycoNonstElms_rpvr(ielmc,N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
        CPOOLT                                = DOM_vr(idom_doc,K,L)+RootMycoNonstElms_rpvr(ielmc,N,L,NZ)
        ZPOOLX                                = 0.1_r8*RootMycoNonstElms_rpvr(ielmn,N,L,NZ)
        PPOOLX                                = 0.1_r8*RootMycoNonstElms_rpvr(ielmp,N,L,NZ)
        XFRE(ielmn)                           = (DOM_vr(idom_don,K,L)*RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-ZPOOLX*DOM_vr(idom_doc,K,L))/CPOOLT
        XFRE(ielmp)                           = (DOM_vr(idom_dop,K,L)*RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-PPOOLX*DOM_vr(idom_doc,K,L))/CPOOLT
        RootMycoExudEUptk_pvr(ielmn,N,K,L,NZ) = FEXUDE(ielmn)*XFRE(ielmn)
        RootMycoExudEUptk_pvr(ielmp,N,K,L,NZ) = FEXUDE(ielmp)*XFRE(ielmp)
      ELSE
        RootMycoExudEUptk_pvr(ielmn,N,K,L,NZ)=0.0_r8
        RootMycoExudEUptk_pvr(ielmp,N,K,L,NZ)=0.0_r8
      ENDIF
    ELSE
      RootMycoExudEUptk_pvr(1:NumPlantChemElms,N,K,L,NZ)=0.0_r8
    ENDIF

  ENDDO D195
!  if(I>176)print*,'rootexudd195',L
  !avoid excessive exudation 
  DO NE=1,NumPlantChemElms
    RootExudE=0._r8
    DO K=1,jcplx
      RootExudE=RootExudE+RootMycoExudEUptk_pvr(NE,N,K,L,NZ)
    ENDDO

    if(RootExudE<0._r8)then
      if(RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootExudE<0._r8)then
        scal=RootMycoNonstElms_rpvr(NE,N,L,NZ)/(-RootExudE)*0.999999_r8
        DO K=1,jcplx
          RootMycoExudEUptk_pvr(NE,N,K,L,NZ)=RootMycoExudEUptk_pvr(NE,N,K,L,NZ)*scal
        ENDDO
      endif
    endif
    
  ENDDO      
!  if(I>176)print*,'exexud'
  end associate
  end subroutine RootExudates
!------------------------------------------------------------------------

  subroutine SumNutrientUptake(NZ)

  implicit none
  integer, intent(in) :: NZ
  real(r8) :: DOM_uptk(1:NumPlantChemElms)
  real(r8) :: scal
  integer :: L  
  integer :: K,NE,N
  !     begin_execution
  associate(                                                   &
    NU                     => plt_site%NU,                     &
    ZERO                   => plt_site%ZERO,                   &
    ZEROS2                 => plt_site%ZEROS2,                 &
    THETW_vr               => plt_soilchem%THETW_vr,           &
    VLSoilPoreMicP_vr      => plt_soilchem%VLSoilPoreMicP_vr,  &
    REcoDOMProd_vr         => plt_bgcr%REcoDOMProd_vr,         &
    RootMycoNonstElms_rpvr => plt_biom%RootMycoNonstElms_rpvr, &
    RootMycoExudEUptk_pvr  => plt_rbgc%RootMycoExudEUptk_pvr,  &
    RootMycoExudElms_pft   => plt_rbgc%RootMycoExudElms_pft,   &
    RootHPO4Uptake_pft     => plt_rbgc%RootHPO4Uptake_pft,     &
    RootH2PO4Uptake_pft    => plt_rbgc%RootH2PO4Uptake_pft,    &
    RootNH4Uptake_pft      => plt_rbgc%RootNH4Uptake_pft,      &
    RootNO3Uptake_pft      => plt_rbgc%RootNO3Uptake_pft,      &
    DOM_vr                 => plt_soilchem%DOM_vr,             &
    MY                     => plt_morph%MY,                    &
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr,      &
    MaxSoiL4Root_pft       => plt_morph%MaxSoiL4Root_pft       &
  )
  !
  !     TOTAL C,N,P EXCHANGE BETWEEN ROOTS AND SOIL
  !
  !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
  !     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
  !     RootNH4Uptake_pvr,RootNH4BUptake_pvr,RUPN03,RootNO3BUptake_pvr=uptake from non-band,band of NH4,NO3
  !     RootH2PO4Uptake_pvr,RootNutUptake_pvr,RootNutUptake_pvr,RootH1PO4BUptake_pvr=uptake from non-band,band of H2PO4,HPO4
  !     RootNH4Uptake_pft,RootNO3Uptake_pft,RootH2PO4Uptake_pft,RootHPO4Uptake_pft=PFT uptake of NH4,NO3,H2PO4,HPO4
  !
  D9501: DO L=NU,MaxSoiL4Root_pft(NZ)
    IF(VLSoilPoreMicP_vr(L).GT.ZEROS2 .AND. THETW_vr(L).GT.ZERO)THEN  
      D295: DO K=1,jcplx
        DOM_uptk=0._r8
        DO N=1,MY(NZ)
          !positve means uptake from soil
          DO NE=1,NumPlantChemElms
            DOM_uptk(NE)=DOM_uptk(NE)+RootMycoExudEUptk_pvr(NE,N,K,L,NZ)
          ENDDO
        ENDDO
        if(any(DOM_uptk/=0._r8))then
          DO NE=1,NumPlantChemElms
            if(DOM_uptk(NE)>DOM_vr(NE,K,L))then
              scal=DOM_uptk(NE)/DOM_vr(NE,K,L)
              DO N=1,MY(NZ)
                RootMycoExudEUptk_pvr(NE,N,K,L,NZ)=RootMycoExudEUptk_pvr(NE,N,K,L,NZ)*scal
              ENDDO
              DOM_uptk(NE)=DOM_uptk(NE)*scal
            endif
            DOM_vr(NE,K,L)=DOM_vr(NE,K,L)-DOM_uptk(NE)
          ENDDO

          DO NE=1,NumPlantChemElms
            REcoDOMProd_vr(NE,K,L)=REcoDOMProd_vr(NE,K,L)-DOM_uptk(NE)
          ENDDO

          DO N=1,MY(NZ)
            DO NE=1,NumPlantChemElms
              RootMycoExudElms_pft(NE,NZ)       = RootMycoExudElms_pft(NE,NZ)+RootMycoExudEUptk_pvr(NE,N,K,L,NZ)
              RootMycoNonstElms_rpvr(NE,N,L,NZ) = RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootMycoExudEUptk_pvr(NE,N,K,L,NZ)
            ENDDO
          ENDDO
        endif
      ENDDO D295
      
      DO N=1,MY(NZ)
        RootNH4Uptake_pft(NZ)   = RootNH4Uptake_pft(NZ)  +(RootNutUptake_pvr(ids_NH4,N,L,NZ)+RootNutUptake_pvr(ids_NH4B,N,L,NZ))
        RootNO3Uptake_pft(NZ)   = RootNO3Uptake_pft(NZ)  +(RootNutUptake_pvr(ids_NO3,N,L,NZ)+RootNutUptake_pvr(ids_NO3B,N,L,NZ))
        RootH2PO4Uptake_pft(NZ) = RootH2PO4Uptake_pft(NZ) +(RootNutUptake_pvr(ids_H2PO4,N,L,NZ)+RootNutUptake_pvr(ids_H2PO4B,N,L,NZ))
        RootHPO4Uptake_pft(NZ)  = RootHPO4Uptake_pft(NZ)  +(RootNutUptake_pvr(ids_H1PO4,N,L,NZ)+RootNutUptake_pvr(ids_H1PO4B,N,L,NZ))
      ENDDO        
    endif
  ENDDO  D9501

  end associate
  end subroutine SumNutrientUptake

end module NutUptakeMod
