module ExtractsMod
!!
!Description:
!     THIS SUBROUTINE AGGREGATES ALL SOIL-PLANT C,N,P EXCHANGES
!     FROM 'UPTAKE' AMD 'GROSUB' AND SENDS RESULTS TO 'REDIST'
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcosimConst
  use GrosubPars
  use PlantAPIData
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: extracts
  contains

  SUBROUTINE extracts(I,J)
!     execution begins here
  implicit none

  integer, intent(in) :: I, J

  integer :: NZ

  call TotalLitrFall()

  DO NZ=1,plt_site%NP
    IF(plt_pheno%IsPlantActive_pft(NZ).EQ.iPlantIsActive)THEN

      call TotalLeafArea(NZ)

      call TotalGasandSoluteUptake(NZ)

      call CanopyFluxesandFixation(NZ)

    ENDIF
  ENDDO
  RETURN
  END subroutine extracts
!------------------------------------------------------------------------------------------

  subroutine TotalLitrFall()

  implicit none
  integer :: NZ,L,K,M
  integer :: NE
  associate(                             &
   NP0                      => plt_site%NP0        , &
   WGLFT                    => plt_biom%WGLFT      , &
   StandingDeadChemElm_col  => plt_biom%StandingDeadChemElm_col   , &
   StandingDeadChemElms_pft => plt_biom%StandingDeadChemElms_pft     , &
   LitrFallChemElm_pft      => plt_bgcr%LitrFallChemElm_pft      , &
   LitrFallChemElm_col      => plt_bgcr%LitrFallChemElm_col      , &
   LitrfalChemElemnts_vr    => plt_bgcr%LitrfalChemElemnts_vr       , &
   LitrFallChemElm_pvr      => plt_bgcr%LitrFallChemElm_pvr       , &
   MaxSoiL4Root             => plt_morph%MaxSoiL4Root       , &
   CanopyStemA_lyr          => plt_morph%CanopyStemA_lyr     , &
   CanopyLAgrid_lyr         =>  plt_morph%CanopyLAgrid_lyr    , &
   StemArea_grd             =>  plt_morph%StemArea_grd    , &
   CanopyLeafArea_grd       =>  plt_morph%CanopyLeafArea_grd      &
  )
  DO NZ=1,NP0
!
!   TOTAL LitrFall OF ALL PLANT SPECIES
!
!   ZCSNC,ZZSNC,ZPSNC=total C,N,P LitrFall
!   HCSNC,HZSNC,HPSNC=hourly PFT C,N,P LitrFall from grosub.f
!   StandingDeadChemElm_col=total standing dead C,N,P mass
!   WTSTG=PFT standing dead C,N,P mass
!   LitrFallChemElm_pvr,=cumulative PFT C,N,P LitrFall from grosub.f
!   LitrfalChemElemnts_vr,=cumulative total C,N,P LitrFall
!
    DO NE=1,NumPlantChemElms
      LitrFallChemElm_col(NE)=LitrFallChemElm_col(NE)+LitrFallChemElm_pft(NE,NZ)
      StandingDeadChemElm_col(NE)=StandingDeadChemElm_col(NE)+StandingDeadChemElms_pft(NE,NZ)
    ENDDO

    DO  L=0,MaxSoiL4Root(NZ)
      DO K=1,pltpar%NumOfPlantLitrCmplxs
        DO NE=1,NumPlantChemElms
          DO  M=1,pltpar%jsken
            LitrfalChemElemnts_vr(NE,M,K,L)=LitrfalChemElemnts_vr(NE,M,K,L)+LitrFallChemElm_pvr(NE,M,K,L,NZ)
          enddo
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CanopyLeafArea_grd=0._r8
  StemArea_grd=0._r8
  DO  L=1,NumOfCanopyLayers1
    CanopyLAgrid_lyr(L)=0._r8
    WGLFT(L)=0._r8
    CanopyStemA_lyr(L)=0._r8
  ENDDO
  end associate
  end subroutine TotalLitrFall
!------------------------------------------------------------------------------------------

  subroutine TotalLeafArea(NZ)
!
!     TOTAL LEAF AREA OF ALL PLANT SPECIES
!
!     CanopyLAgrid_lyr,CanopyStemA_lyr=total leaf,stalk area of combined canopy layer
!     CanopyLeafApft_lyr,CanopyStemApft_lyr=PFT leaf,stalk area in canopy layer
!     WGLFT=total leaf C of combined canopy layer
!     CanopyLeafCpft_lyr=PFT leaf C in canopy layer
!
  implicit none
  integer, intent(in) :: NZ
  integer :: L
  associate(                              &
    CanopyLeafCpft_lyr    => plt_biom%CanopyLeafCpft_lyr      , &
    WGLFT                 => plt_biom%WGLFT      , &
    CanopyLAgrid_lyr      =>  plt_morph%CanopyLAgrid_lyr    , &
    CanopyStemApft_lyr    =>  plt_morph%CanopyStemApft_lyr    , &
    CanopyStemA_lyr       => plt_morph%CanopyStemA_lyr     , &
    CanopyLeafApft_lyr    => plt_morph%CanopyLeafApft_lyr       &
  )
  DO L=1,NumOfCanopyLayers1
    CanopyLAgrid_lyr(L)=CanopyLAgrid_lyr(L)+CanopyLeafApft_lyr(L,NZ)
    WGLFT(L)=WGLFT(L)+CanopyLeafCpft_lyr(L,NZ)
    CanopyStemA_lyr(L)=CanopyStemA_lyr(L)+CanopyStemApft_lyr(L,NZ)
  ENDDO
  end associate
  end subroutine TotalLeafArea
!------------------------------------------------------------------------------------------

  subroutine TotalGasandSoluteUptake(NZ)
!
!     TOTAL GAS AND SOLUTE UPTAKE BY ALL PLANT SPECIES
!

  implicit none
  integer, intent(in) :: NZ

  integer :: N,L,K,NTG,NE

  associate(                       &
    NU                       => plt_site%NU     , &
    AREA3                    => plt_site%AREA3  , &
    PlantPopulation_pft      => plt_site%PlantPopulation_pft     , &
    RUPP1B                   => plt_rbgc%RUPP1B , &
    RUPP2B                   => plt_rbgc%RUPP2B , &
    RUNNXP                   => plt_rbgc%RUNNXP , &
    RUPOXP                   => plt_rbgc%RUPOXP , &
    trcg_Root_DisEvap_flx_vr => plt_rbgc%trcg_Root_DisEvap_flx_vr , &
    trcg_air2root_flx_pft_vr => plt_rbgc%trcg_air2root_flx_pft_vr , &
    RCO2P                    => plt_rbgc%RCO2P  , &
    RUPGasSol_vr             => plt_rbgc%RUPGasSol_vr , &
    RootNutUptake_pvr        => plt_rbgc%RootNutUptake_pvr , &
    trcg_air2root_flx_vr     => plt_rbgc%trcg_air2root_flx_vr , &
    trcg_TLP                 => plt_rbgc%trcg_TLP , &
    ROXYP                    => plt_rbgc%ROXYP  , &
    RDFOME                   => plt_rbgc%RDFOME , &
    RUNNHP                   => plt_rbgc%RUNNHP , &
    RUNNOP                   => plt_rbgc%RUNNOP , &
    RUPP2P                   => plt_rbgc%RUPP2P , &
    RUNNBP                   => plt_rbgc%RUNNBP , &
    RUPP1P                   => plt_rbgc%RUPP1P , &
    trcs_plant_uptake_vr     => plt_rbgc%trcs_plant_uptake_vr , &
    RNO3X                    => plt_bgcr%RNO3X  , &
    RNH4X                    => plt_bgcr%RNH4X  , &
    RPO4X                    => plt_bgcr%RPO4X  , &
    RN3BX                    => plt_bgcr%RN3BX  , &
    RP14X                    => plt_bgcr%RP14X  , &
    RNHBX                    => plt_bgcr%RNHBX  , &
    ROXYX                    => plt_bgcr%ROXYX  , &
    TDFOME                   => plt_bgcr%TDFOME , &
    TUPOXP                   => plt_bgcr%TUPOXP , &
    TCO2P                    => plt_bgcr%TCO2P  , &
    RPOBX                    => plt_bgcr%RPOBX  , &
    RP1BX                    => plt_bgcr%RP1BX  , &
    TKS                      => plt_ew%TKS      , &
    THeatRootUptake          => plt_ew%THeatRootUptake    , &
    GridPlantRootH2OUptake_vr=> plt_ew%GridPlantRootH2OUptake_vr   , &
    AllPlantRootH2OUptake_vr => plt_ew%AllPlantRootH2OUptake_vr    , &
    trcg_rootml_vr           => plt_rbgc%trcg_rootml_vr,&
    trcs_rootml_vr           => plt_rbgc%trcs_rootml_vr, &
    RootLenDensPerPlant_pvr  => plt_morph%RootLenDensPerPlant_pvr , &
    RTDNT                    => plt_morph%RTDNT , &
    MY                       => plt_morph%MY    , &
    MaxSoiL4Root             => plt_morph%MaxSoiL4Root     &
  )

  DO N=1,MY(NZ)
    DO L=NU,MaxSoiL4Root(NZ)
!
!     TOTAL ROOT DENSITY
!
!     RTDNT=total root length density
!     RootLenDensPerPlant_pvr=PFT root length density per plant
!     AllPlantRootH2OUptake_vr=total water uptake
!     AllPlantRootH2OUptake_vr=PFT root water uptake
!     THeatRootUptake=total convective heat in root water uptake
!     TKS=soil temperature
!     PP=PFT population, this is dynamic, and can goes to zero
!
      IF(N.EQ.ipltroot)THEN
        RTDNT(L)=RTDNT(L)+RootLenDensPerPlant_pvr(N,L,NZ)*PlantPopulation_pft(NZ)/AREA3(L)
      ENDIF
!
!     TOTAL WATER UPTAKE
!
      GridPlantRootH2OUptake_vr(L)=GridPlantRootH2OUptake_vr(L)+AllPlantRootH2OUptake_vr(N,L,NZ)
      THeatRootUptake(L)=THeatRootUptake(L)+AllPlantRootH2OUptake_vr(N,L,NZ)*cpw*TKS(L)
!
!     ROOT GAS CONTENTS FROM FLUXES IN 'UPTAKE'
!
!     *A,*P=PFT root gaseous, aqueous gas content
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     R*FLA=root gaseous-atmosphere CO2 exchange
!     R*DFA=root aqueous-gaseous CO2 exchange
!
      DO NTG=idg_beg,idg_NH3
        trcg_rootml_vr(NTG,N,L,NZ)=trcg_rootml_vr(NTG,N,L,NZ) &
          +trcg_air2root_flx_pft_vr(NTG,N,L,NZ)-trcg_Root_DisEvap_flx_vr(NTG,N,L,NZ)
      ENDDO

      trcs_rootml_vr(idg_CO2,N,L,NZ)=trcs_rootml_vr(idg_CO2,N,L,NZ)+trcg_Root_DisEvap_flx_vr(idg_CO2,N,L,NZ)+RCO2P(N,L,NZ)
      trcs_rootml_vr(idg_O2,N,L,NZ)=trcs_rootml_vr(idg_O2,N,L,NZ)+trcg_Root_DisEvap_flx_vr(idg_O2,N,L,NZ)-RUPOXP(N,L,NZ)
      trcs_rootml_vr(idg_CH4,N,L,NZ)=trcs_rootml_vr(idg_CH4,N,L,NZ) &
        +trcg_Root_DisEvap_flx_vr(idg_CH4,N,L,NZ)+RUPGasSol_vr(idg_CH4,N,L,NZ)
      trcs_rootml_vr(idg_N2O,N,L,NZ)=trcs_rootml_vr(idg_N2O,N,L,NZ) &
        +trcg_Root_DisEvap_flx_vr(idg_N2O,N,L,NZ)+RUPGasSol_vr(idg_N2O,N,L,NZ)
      trcs_rootml_vr(idg_NH3,N,L,NZ)=trcs_rootml_vr(idg_NH3,N,L,NZ) &
        +trcg_Root_DisEvap_flx_vr(idg_NH3,N,L,NZ)+RUPGasSol_vr(idg_NH3,N,L,NZ)+RUPGasSol_vr(idg_NH3B,N,L,NZ)
      trcs_rootml_vr(idg_H2,N,L,NZ)=trcs_rootml_vr(idg_H2,N,L,NZ) &
        +trcg_Root_DisEvap_flx_vr(idg_H2,N,L,NZ)+RUPGasSol_vr(idg_H2,N,L,NZ)
!
!     TOTAL ROOT GAS CONTENTS
!
!     TL*P=total root gas content, dissolved + gaseous phase
!     *A,*P=PFT root gaseous, aqueous gas content
!
      DO NTG=idg_beg,idg_end-1
        trcg_TLP(NTG,L)=trcg_TLP(NTG,L)+trcs_rootml_vr(NTG,N,L,NZ)+trcg_rootml_vr(NTG,N,L,NZ)
      ENDDO
!
!     TOTAL ROOT BOUNDARY GAS FLUXES
!
!     T*FLA=total root gaseous-atmosphere CO2 exchange
!     R*FLA=PFT root gaseous-atmosphere CO2 exchange
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     TUP*S,TUP*B=total root-soil gas, solute exchange in non-band,band
!     RUP*S,RUP*B*=PFT root-soil gas, solute exchange in non-band,band
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     solute code:NH4=NH4,NO3=NO3,H2P=H2PO4,H1P=H1PO4 in non-band
!                :NHB=NH4,NOB=NO3,H2B=H2PO4,H1B=H1PO4 in band
!
      DO NTG=idg_beg,idg_end-1
        trcg_air2root_flx_vr(NTG,L)=trcg_air2root_flx_vr(NTG,L)+trcg_air2root_flx_pft_vr(NTG,N,L,NZ)
      ENDDO

      TCO2P(L)=TCO2P(L)-RCO2P(N,L,NZ)
      TUPOXP(L)=TUPOXP(L)+RUPOXP(N,L,NZ)
      trcs_plant_uptake_vr(idg_CO2,L)=trcs_plant_uptake_vr(idg_CO2,L)+RUPGasSol_vr(idg_CO2,N,L,NZ)
      trcs_plant_uptake_vr(idg_O2,L)=trcs_plant_uptake_vr(idg_O2,L)+RUPGasSol_vr(idg_O2,N,L,NZ)
      trcs_plant_uptake_vr(idg_CH4,L)=trcs_plant_uptake_vr(idg_CH4,L)+RUPGasSol_vr(idg_CH4,N,L,NZ)
      trcs_plant_uptake_vr(idg_N2O,L)=trcs_plant_uptake_vr(idg_N2O,L)+RUPGasSol_vr(idg_N2O,N,L,NZ)
      trcs_plant_uptake_vr(idg_NH3,L)=trcs_plant_uptake_vr(idg_NH3,L)+RUPGasSol_vr(idg_NH3,N,L,NZ)
      trcs_plant_uptake_vr(idg_NH3B,L)=trcs_plant_uptake_vr(idg_NH3B,L)+RUPGasSol_vr(idg_NH3B,N,L,NZ)
      trcs_plant_uptake_vr(idg_H2,L)=trcs_plant_uptake_vr(idg_H2,L)+RUPGasSol_vr(idg_H2,N,L,NZ)
      trcs_plant_uptake_vr(idg_N2,L)=trcs_plant_uptake_vr(idg_N2,L)+RUPGasSol_vr(idg_N2,N,L,NZ)

      trcs_plant_uptake_vr(ids_NH4,L)=trcs_plant_uptake_vr(ids_NH4,L)+RootNutUptake_pvr(ids_NH4,N,L,NZ)
      trcs_plant_uptake_vr(ids_NO3,L)=trcs_plant_uptake_vr(ids_NO3,L)+RootNutUptake_pvr(ids_NO3,N,L,NZ)
      trcs_plant_uptake_vr(ids_H2PO4,L)=trcs_plant_uptake_vr(ids_H2PO4,L)+RootNutUptake_pvr(ids_H2PO4,N,L,NZ)
      trcs_plant_uptake_vr(ids_H1PO4,L)=trcs_plant_uptake_vr(ids_H1PO4,L)+RootNutUptake_pvr(ids_H1PO4,N,L,NZ)
      trcs_plant_uptake_vr(ids_NH4B,L)=trcs_plant_uptake_vr(ids_NH4B,L)+RootNutUptake_pvr(ids_NH4B,N,L,NZ)
      trcs_plant_uptake_vr(ids_NO3B,L)=trcs_plant_uptake_vr(ids_NO3B,L)+RootNutUptake_pvr(ids_NO3B,N,L,NZ)
      trcs_plant_uptake_vr(ids_H2PO4B,L)=trcs_plant_uptake_vr(ids_H2PO4B,L)+RootNutUptake_pvr(ids_H2PO4B,N,L,NZ)
      trcs_plant_uptake_vr(ids_H1PO4B,L)=trcs_plant_uptake_vr(ids_H1PO4B,L)+RootNutUptake_pvr(ids_H1PO4B,N,L,NZ)
!
!     TOTAL ROOT C,N,P EXUDATION
!
!     TDFOME=total nonstructl C,N,P exchange
!     RDFOMC,RDFOMN,RDFOMP=PFT nonstructl C,N,P exchange
!
      DO K=1,jcplx
        DO NE=1,NumPlantChemElms
          TDFOME(NE,K,L)=TDFOME(NE,K,L)-RDFOME(ielmc,N,K,L,NZ)
        ENDDO
      ENDDO
!
!     TOTAL ROOT O2, NH4, NO3, PO4 UPTAKE CONTRIBUTES TO
!     TOTAL ROOT + MICROBIAL UPTAKE USED TO CALCULATE
!     COMPETITION CONSTRAINTS
!
!     ROXYX=O2 demand by all microbial,root,myco populations
!     RNH4X=NH4 demand in non-band by all microbial,root,myco populations
!     RNO3X=NO3 demand in non-band by all microbial,root,myco populations
!     RPO4X=H2PO4 demand in non-band by all microbial,root,myco populations
!     RP14X=HPO4 demand in non-band by all microbial,root,myco populations
!     RNHBX=NH4 demand in band by all microbial,root,myco populations
!     RN3BX=NO3 demand in band by all microbial,root,myco populations
!     RPOBX=H2PO4 demand in band by all microbial,root,myco populations
!     RP1BX=HPO4 demand in band by all microbial,root,myco populations
!     ROXYP=O2 demand by each root,myco population
!     RUNNHP=NH4 demand in non-band by each root population
!     RUNNOP=NO3 demand in non-band by each root population
!     RUPP2P=H2PO4 demand in non-band by each root population
!     RUPP1P=HPO4 demand in non-band by each root population
!     RUNNBP=NH4 demand in band by each root population
!     RUNNXB=NO3 demand in band by each root population
!     RUPP2B=H2PO4 demand in band by each root population
!     RUPP1B=HPO4 demand in band by each root population
!
      ROXYX(L)=ROXYX(L)+ROXYP(N,L,NZ)
      RNH4X(L)=RNH4X(L)+RUNNHP(N,L,NZ)
      RNO3X(L)=RNO3X(L)+RUNNOP(N,L,NZ)
      RPO4X(L)=RPO4X(L)+RUPP2P(N,L,NZ)
      RP14X(L)=RP14X(L)+RUPP1P(N,L,NZ)
      RNHBX(L)=RNHBX(L)+RUNNBP(N,L,NZ)
      RN3BX(L)=RN3BX(L)+RUNNXP(N,L,NZ)
      RPOBX(L)=RPOBX(L)+RUPP2B(N,L,NZ)
      RP1BX(L)=RP1BX(L)+RUPP1B(N,L,NZ)
    ENDDO
  ENDDO
  end associate
  end subroutine TotalGasandSoluteUptake
!------------------------------------------------------------------------------------------

  subroutine CanopyFluxesandFixation(NZ)
!
!     TOTAL ROOT N2 FIXATION BY ALL PLANT SPECIES
!
!     TRootN2Fix_pft=total root N2 fixation
!     RootN2Fix_pvr=PFT root N2 fixation
!

  implicit none
  integer, intent(in) :: NZ
  integer :: L, NE,NB,NTG
  real(r8) :: ENGYC

  associate(                       &
    PlantElemntStoreLandscape  => plt_site%PlantElemntStoreLandscape  , &
    ElmntBalanceCum_pft        => plt_site%ElmntBalanceCum_pft   , &
    NH3EmiCum_pft              => plt_bgcr%NH3EmiCum_pft  , &
    NH3Dep2Can_pft             => plt_bgcr%NH3Dep2Can_pft  , &
    Canopy_NEE_col             => plt_bgcr%Canopy_NEE_col  , &
    LitrFallChemElm_col        => plt_bgcr%LitrFallChemElm_col  , &
    RootGasLossDisturb_pft     => plt_bgcr%RootGasLossDisturb_pft, &
    RootN2Fix_pvr              => plt_bgcr%RootN2Fix_pvr  , &
    CO2NetFix_pft              => plt_bgcr%CO2NetFix_pft   , &
    ETCanopy_pft               => plt_ew%ETCanopy_pft    , &
    TH2GZ                      => plt_bgcr%TH2GZ  , &
    trcs_plant_uptake_vr       => plt_rbgc%trcs_plant_uptake_vr, &    
    NH3Dep2_brch               => plt_rbgc%NH3Dep2_brch  , &
    PlantRootSoilChemNetX_pft  => plt_rbgc%PlantRootSoilChemNetX_pft , &
    TRootGasLossDisturb_pft    => plt_rbgc%TRootGasLossDisturb_pft  , &
    Transpiration_pft          => plt_ew%Transpiration_pft      , &
    PrecIntcptByCanopy_pft     => plt_ew%PrecIntcptByCanopy_pft     , &
    VapXAir2Canopy_pft         => plt_ew%VapXAir2Canopy_pft   , &
    WatByPCanopy               => plt_ew%WatByPCanopy    , &
    CanopyWater_pft            => plt_ew%CanopyWater_pft    , &
    Eco_Heat_Grnd_col          => plt_ew%Eco_Heat_Grnd_col      , &
    HeatXAir2PCan              => plt_ew%HeatXAir2PCan    , &
    EvapTransHeat_pft          => plt_ew%EvapTransHeat_pft    , &
    CanWatg                    => plt_ew%CanWatg   , &
    TKC                        => plt_ew%TKC      , &
    TKS                        => plt_ew%TKS      , &
    ENGYX                      => plt_ew%ENGYX    , &
    Eco_Heat_Sens_col          => plt_ew%Eco_Heat_Sens_col      , &
    VapXAir2CanG               => plt_ew%VapXAir2CanG   , &
    TENGYC                     => plt_ew%TENGYC   , &
    TEVAPP                     => plt_ew%TEVAPP   , &
    THFLXC                     => plt_ew%THFLXC   , &
    LWRadCanG                  => plt_ew%LWRadCanG    , &
    TairK                      => plt_ew%TairK      , &
    HeatStorCanP               => plt_ew%HeatStorCanP    , &
    Eco_Heat_Latent_col        => plt_ew%Eco_Heat_Latent_col      , &
    CanH2OHeldVg               => plt_ew%CanH2OHeldVg   , &
    NU                         => plt_site%NU     , &
    StemArea_grd               => plt_morph%StemArea_grd , &
    CanopyLeafArea_grd         => plt_morph%CanopyLeafArea_grd , &
    MaxSoiL4Root               => plt_morph%MaxSoiL4Root   , &
    NumOfBranches_pft          => plt_morph%NumOfBranches_pft   , &
    CanopyStemA_pft            => plt_morph%CanopyStemA_pft , &
    CanopyLeafArea_pft         => plt_morph%CanopyLeafArea_pft , &
    RadNet2CanP                => plt_rad%RadNet2CanP    , &
    LWRadCanP                  => plt_rad%LWRadCanP   , &
    Eco_NetRad_col             => plt_rad%Eco_NetRad_col       &
  )
  DO L=NU,MaxSoiL4Root(NZ)
    trcs_plant_uptake_vr(idg_N2,L)=trcs_plant_uptake_vr(idg_N2,L)+RootN2Fix_pvr(L,NZ)
  ENDDO
!
!     TOTAL ENERGY, WATER, CO2 FLUXES
!
!     Eco_NetRad_col=total net SW+LW absorbed by canopy
!     RadNet2CanP=PFT net SW+LW absorbed by canopy
!     Eco_Heat_Latent_col=total canopy latent heat flux
!     EvapTransHeat_pft=PFT canopy latent heat flux
!     Eco_Heat_Sens_col=total canopy sensible heat flux
!     HeatXAir2PCan=PFT canopy sensible heat flux
!     Eco_Heat_Grnd_col=total canopy storage heat flux
!     HeatStorCanP=PFT canopy storage heat flux
!     Canopy_NEE_col=total net CO2 fixation
!     CO2NetFix_pft=PFT net CO2 fixation
!     CanWatg,CanH2OHeldVg=total water volume in canopy,on canopy surfaces
!     CanopyWater_pft,WatByPCanopy=PFT water volume in canopy,on canopy surfaces
!     TEVAPP,VapXAir2CanG=total water flux to,from canopy,canopy surfaces
!     VapXAir2PCan,Transpiration_pft=water flux to,from canopy surfaces, inside canopy
!     TENGYC=total canopy water heat content
!     ENGYC=PFT canopy water heat content
!     CanopyLeafArea_grd,StemArea_grd=total leaf,stalk area
!     CanopyLeafArea_pft,CanopyStemA_pft=PFT leaf,stalk area
!     ZCSNC,ZZSNC,ZPSNC=total net root-soil C,N,P exchange
!     HCUPTK,HZUPTK,HPUPTK=PFT net root-soil C,N,P exchange
!     TBALC,TBALN,TBALP=total C,N,P balance
!     BALC,BALN,BALP=PFT C,N,P balance
!     TRootGasLossDisturb_pft=total loss of root CO2, O2, CH4, N2O, NH3, H2
!     RootGasLossDisturb_pft=PFT loss of root CO2, O2, CH4, N2O, NH3, H2
!
  Eco_NetRad_col=Eco_NetRad_col+RadNet2CanP(NZ)
  Eco_Heat_Latent_col=Eco_Heat_Latent_col+EvapTransHeat_pft(NZ)
  Eco_Heat_Sens_col=Eco_Heat_Sens_col+HeatXAir2PCan(NZ)
  Eco_Heat_Grnd_col=Eco_Heat_Grnd_col+HeatStorCanP(NZ)
  Canopy_NEE_col=Canopy_NEE_col+CO2NetFix_pft(NZ)
  ETCanopy_pft(NZ)=ETCanopy_pft(NZ)+Transpiration_pft(NZ)+VapXAir2Canopy_pft(NZ)
  CanWatg=CanWatg+CanopyWater_pft(NZ)
  CanH2OHeldVg=CanH2OHeldVg+WatByPCanopy(NZ)
  TEVAPP=TEVAPP+Transpiration_pft(NZ)+VapXAir2Canopy_pft(NZ)
  VapXAir2CanG=VapXAir2CanG+VapXAir2Canopy_pft(NZ)
  ENGYC=cpw*(WatByPCanopy(NZ)+PrecIntcptByCanopy_pft(NZ)+VapXAir2Canopy_pft(NZ))*TKC(NZ)
  TENGYC=TENGYC+ENGYC
  THFLXC=THFLXC+ENGYC-ENGYX(NZ)-(PrecIntcptByCanopy_pft(NZ)*cpw*TairK)
  ENGYX(NZ)=ENGYC
  LWRadCanG=LWRadCanG+LWRadCanP(NZ)
  CanopyLeafArea_grd=CanopyLeafArea_grd+CanopyLeafArea_pft(NZ)
  StemArea_grd=StemArea_grd+CanopyStemA_pft(NZ)
  DO NE=1,NumPlantChemElms
    LitrFallChemElm_col(NE)=LitrFallChemElm_col(NE)-PlantRootSoilChemNetX_pft(NE,NZ)
    PlantElemntStoreLandscape(NE)=PlantElemntStoreLandscape(NE)+ElmntBalanceCum_pft(NE,NZ)
  ENDDO

  DO NTG=idg_beg,idg_end-1
    TRootGasLossDisturb_pft(NTG)=TRootGasLossDisturb_pft(NTG)+RootGasLossDisturb_pft(NTG,NZ)
  ENDDO
!
!     TOTAL CANOPY NH3 EXCHANGE AND EXUDATION
!
!     NH3Dep2_brch,NH3Dep2Can_pft=PFT NH3 flux between atmosphere and branch,canopy
!     NH3EmiCum_pft=total NH3 flux between atmosphere and canopy
!
  NH3Dep2Can_pft(NZ)=0._r8
  DO NB=1,NumOfBranches_pft(NZ)
    NH3Dep2Can_pft(NZ)=NH3Dep2Can_pft(NZ)+NH3Dep2_brch(NB,NZ)
    NH3EmiCum_pft(NZ)=NH3EmiCum_pft(NZ)+NH3Dep2_brch(NB,NZ)
  ENDDO
  end associate
  end subroutine CanopyFluxesandFixation

  end module ExtractsMod
