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

  call TotalLitterfall()

  DO NZ=1,plt_site%NP
    IF(plt_pheno%IsPlantActive(NZ).EQ.iPlantIsActive)THEN

      call TotalLeafArea(NZ)

      call TotalGasandSoluteUptake(NZ)

      call CanopyFluxesandFixation(NZ)

    ENDIF
  ENDDO
  RETURN
  END subroutine extracts
!------------------------------------------------------------------------------------------

  subroutine TotalLitterfall()

  implicit none
  integer :: NZ,L,K,M
  integer :: NE
  associate(                             &
   NP0      => plt_site%NP0        , &
   WGLFT    => plt_biom%WGLFT      , &
   StandingDeadChemElmnt_col  => plt_biom%StandingDeadChemElmnt_col   , &
   StandingDeadChemElmnts_pft   => plt_biom%StandingDeadChemElmnts_pft     , &
   LitterFallChemElmnt_pft    => plt_bgcr%LitterFallChemElmnt_pft      , &
   LitterFallChemElmnt_col    => plt_bgcr%LitterFallChemElmnt_col      , &
   LitrfalChemElemnts_vr     => plt_bgcr%LitrfalChemElemnts_vr       , &
   LitterFallChemElmnt_pftvr     => plt_bgcr%LitterFallChemElmnt_pftvr       , &
   NI       => plt_morph%NI        , &
   CanopyStemA_lyr    => plt_morph%CanopyStemA_lyr     , &
   CanopyLAgrid_lyr    =>  plt_morph%CanopyLAgrid_lyr    , &
   StemAreag    =>  plt_morph%StemAreag    , &
   CanopyLA_grd    =>  plt_morph%CanopyLA_grd      &
  )
  DO NZ=1,NP0
!
!   TOTAL LITTERFALL OF ALL PLANT SPECIES
!
!   ZCSNC,ZZSNC,ZPSNC=total C,N,P litterfall
!   HCSNC,HZSNC,HPSNC=hourly PFT C,N,P litterfall from grosub.f
!   StandingDeadChemElmnt_col=total standing dead C,N,P mass
!   WTSTG=PFT standing dead C,N,P mass
!   LitterFallChemElmnt_pftvr,=cumulative PFT C,N,P litterfall from grosub.f
!   LitrfalChemElemnts_vr,=cumulative total C,N,P litterfall
!
    DO NE=1,NumOfPlantChemElmnts
      LitterFallChemElmnt_col(NE)=LitterFallChemElmnt_col(NE)+LitterFallChemElmnt_pft(NE,NZ)
      StandingDeadChemElmnt_col(NE)=StandingDeadChemElmnt_col(NE)+StandingDeadChemElmnts_pft(NE,NZ)
    ENDDO

    DO  L=0,NI(NZ)
      DO K=1,pltpar%NumOfPlantLitrCmplxs
        DO NE=1,NumOfPlantChemElmnts
          DO  M=1,pltpar%jsken
            LitrfalChemElemnts_vr(NE,M,K,L)=LitrfalChemElemnts_vr(NE,M,K,L)+LitterFallChemElmnt_pftvr(NE,M,K,L,NZ)
          enddo
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CanopyLA_grd=0._r8
  StemAreag=0._r8
  DO  L=1,NumOfCanopyLayers1
    CanopyLAgrid_lyr(L)=0._r8
    WGLFT(L)=0._r8
    CanopyStemA_lyr(L)=0._r8
  ENDDO
  end associate
  end subroutine TotalLitterfall
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
    WGLFT    => plt_biom%WGLFT      , &
    CanopyLAgrid_lyr    =>  plt_morph%CanopyLAgrid_lyr    , &
    CanopyStemApft_lyr    =>  plt_morph%CanopyStemApft_lyr    , &
    CanopyStemA_lyr    => plt_morph%CanopyStemA_lyr     , &
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
    NU    => plt_site%NU     , &
    AREA3 => plt_site%AREA3  , &
    pftPlantPopulation    => plt_site%pftPlantPopulation     , &
    RUPP1B=> plt_rbgc%RUPP1B , &
    RUPP2B=> plt_rbgc%RUPP2B , &
    RUNNXP=> plt_rbgc%RUNNXP , &
    trcg_Root_DisEvap_flx_vr=> plt_rbgc%trcg_Root_DisEvap_flx_vr , &
    trcg_air2root_flx_pft_vr=> plt_rbgc%trcg_air2root_flx_pft_vr , &
    RCO2P => plt_rbgc%RCO2P  , &
    RUPCHS=> plt_rbgc%RUPCHS , &
    RUPOXP=> plt_rbgc%RUPOXP , &
    RUPN2S=> plt_rbgc%RUPN2S , &
    RUPN3B=> plt_rbgc%RUPN3B , &
    RUPHGS=> plt_rbgc%RUPHGS , &
    RUPN3S=> plt_rbgc%RUPN3S , &
    RUPNOB=> plt_rbgc%RUPNOB , &
    RCO2S => plt_rbgc%RCO2S  , &
    RUPOXS=> plt_rbgc%RUPOXS , &
    RootNH4Uptake_pvr=> plt_rbgc%RootNH4Uptake_pvr , &
    RootNO3Uptake_pvr=> plt_rbgc%RootNO3Uptake_pvr , &
    RootHPO4Uptake_pvr=> plt_rbgc%RootHPO4Uptake_pvr , &
    RootH2PO4Uptake_pvr=> plt_rbgc%RootH2PO4Uptake_pvr , &
    RUPNHB=> plt_rbgc%RUPNHB , &
    RUPH2B=> plt_rbgc%RUPH2B , &
    RUPH1B=> plt_rbgc%RUPH1B , &
    trcg_air2root_flx_vr=> plt_rbgc%trcg_air2root_flx_vr , &
    trcg_TLP=> plt_rbgc%trcg_TLP , &
    ROXYP => plt_rbgc%ROXYP  , &
    RDFOME=> plt_rbgc%RDFOME , &
    RUNNHP=> plt_rbgc%RUNNHP , &
    RUNNOP=> plt_rbgc%RUNNOP , &
    RUPP2P=> plt_rbgc%RUPP2P , &
    RUNNBP=> plt_rbgc%RUNNBP , &
    RUPP1P=> plt_rbgc%RUPP1P , &
    trcs_plant_uptake_vr=> plt_rbgc%trcs_plant_uptake_vr , &
    RNO3X => plt_bgcr%RNO3X  , &
    RNH4X => plt_bgcr%RNH4X  , &
    RPO4X => plt_bgcr%RPO4X  , &
    RN3BX => plt_bgcr%RN3BX  , &
    RP14X => plt_bgcr%RP14X  , &
    RNHBX => plt_bgcr%RNHBX  , &
    ROXYX => plt_bgcr%ROXYX  , &
    TDFOME=> plt_bgcr%TDFOME , &
    TUPOXP=> plt_bgcr%TUPOXP , &
    TCO2P => plt_bgcr%TCO2P  , &
    RPOBX => plt_bgcr%RPOBX  , &
    RP1BX => plt_bgcr%RP1BX  , &
    TKS   => plt_ew%TKS      , &
    THeatRootUptake => plt_ew%THeatRootUptake    , &
    GridPlantRootH2OUptake_vr=> plt_ew%GridPlantRootH2OUptake_vr   , &
    AllPlantRootH2OUptake_vr => plt_ew%AllPlantRootH2OUptake_vr    , &
    trcg_rootml  => plt_rbgc%trcg_rootml,&
    trcs_rootml => plt_rbgc%trcs_rootml, &
    RootLenthDensPerPopu_pvr => plt_morph%RootLenthDensPerPopu_pvr , &
    RTDNT => plt_morph%RTDNT , &
    MY    => plt_morph%MY    , &
    NI    => plt_morph%NI      &
  )

  DO N=1,MY(NZ)
    DO L=NU,NI(NZ)
!
!     TOTAL ROOT DENSITY
!
!     RTDNT=total root length density
!     RootLenthDensPerPopu_pvr=PFT root length density per plant
!     AllPlantRootH2OUptake_vr=total water uptake
!     AllPlantRootH2OUptake_vr=PFT root water uptake
!     THeatRootUptake=total convective heat in root water uptake
!     TKS=soil temperature
!     PP=PFT population, this is dynamic, and can goes to zero
!
      IF(N.EQ.ipltroot)THEN
        RTDNT(L)=RTDNT(L)+RootLenthDensPerPopu_pvr(N,L,NZ)*pftPlantPopulation(NZ)/AREA3(L)
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
        trcg_rootml(NTG,N,L,NZ)=trcg_rootml(NTG,N,L,NZ)+trcg_air2root_flx_pft_vr(NTG,N,L,NZ)-trcg_Root_DisEvap_flx_vr(NTG,N,L,NZ)
      ENDDO

      trcs_rootml(idg_CO2,N,L,NZ)=trcs_rootml(idg_CO2,N,L,NZ)+trcg_Root_DisEvap_flx_vr(idg_CO2,N,L,NZ)+RCO2P(N,L,NZ)
      trcs_rootml(idg_O2,N,L,NZ)=trcs_rootml(idg_O2,N,L,NZ)+trcg_Root_DisEvap_flx_vr(idg_O2,N,L,NZ)-RUPOXP(N,L,NZ)
      trcs_rootml(idg_CH4,N,L,NZ)=trcs_rootml(idg_CH4,N,L,NZ)+trcg_Root_DisEvap_flx_vr(idg_CH4,N,L,NZ)+RUPCHS(N,L,NZ)
      trcs_rootml(idg_N2O,N,L,NZ)=trcs_rootml(idg_N2O,N,L,NZ)+trcg_Root_DisEvap_flx_vr(idg_N2O,N,L,NZ)+RUPN2S(N,L,NZ)
      trcs_rootml(idg_NH3,N,L,NZ)=trcs_rootml(idg_NH3,N,L,NZ)+trcg_Root_DisEvap_flx_vr(idg_NH3,N,L,NZ)+RUPN3S(N,L,NZ)+RUPN3B(N,L,NZ)
      trcs_rootml(idg_H2,N,L,NZ)=trcs_rootml(idg_H2,N,L,NZ)+trcg_Root_DisEvap_flx_vr(idg_H2,N,L,NZ)+RUPHGS(N,L,NZ)
!
!     TOTAL ROOT GAS CONTENTS
!
!     TL*P=total root gas content, dissolved + gaseous phase
!     *A,*P=PFT root gaseous, aqueous gas content
!
      DO NTG=idg_beg,idg_end-1
        trcg_TLP(NTG,L)=trcg_TLP(NTG,L)+trcs_rootml(NTG,N,L,NZ)+trcg_rootml(NTG,N,L,NZ)
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
      trcs_plant_uptake_vr(idg_CO2,L)=trcs_plant_uptake_vr(idg_CO2,L)+RCO2S(N,L,NZ)
      trcs_plant_uptake_vr(idg_O2,L)=trcs_plant_uptake_vr(idg_O2,L)+RUPOXS(N,L,NZ)
      trcs_plant_uptake_vr(idg_CH4,L)=trcs_plant_uptake_vr(idg_CH4,L)+RUPCHS(N,L,NZ)
      trcs_plant_uptake_vr(idg_N2O,L)=trcs_plant_uptake_vr(idg_N2O,L)+RUPN2S(N,L,NZ)
      trcs_plant_uptake_vr(idg_NH3,L)=trcs_plant_uptake_vr(idg_NH3,L)+RUPN3S(N,L,NZ)
      trcs_plant_uptake_vr(idg_NH3B,L)=trcs_plant_uptake_vr(idg_NH3B,L)+RUPN3B(N,L,NZ)
      trcs_plant_uptake_vr(idg_H2,L)=trcs_plant_uptake_vr(idg_H2,L)+RUPHGS(N,L,NZ)

      trcs_plant_uptake_vr(ids_NH4,L)=trcs_plant_uptake_vr(ids_NH4,L)+RootNH4Uptake_pvr(N,L,NZ)
      trcs_plant_uptake_vr(ids_NO3,L)=trcs_plant_uptake_vr(ids_NO3,L)+RootNO3Uptake_pvr(N,L,NZ)
      trcs_plant_uptake_vr(ids_H2PO4,L)=trcs_plant_uptake_vr(ids_H2PO4,L)+RootH2PO4Uptake_pvr(N,L,NZ)
      trcs_plant_uptake_vr(ids_H1PO4,L)=trcs_plant_uptake_vr(ids_H1PO4,L)+RootHPO4Uptake_pvr(N,L,NZ)
      trcs_plant_uptake_vr(ids_NH4B,L)=trcs_plant_uptake_vr(ids_NH4B,L)+RUPNHB(N,L,NZ)
      trcs_plant_uptake_vr(ids_NO3B,L)=trcs_plant_uptake_vr(ids_NO3B,L)+RUPNOB(N,L,NZ)
      trcs_plant_uptake_vr(ids_H2PO4B,L)=trcs_plant_uptake_vr(ids_H2PO4B,L)+RUPH2B(N,L,NZ)
      trcs_plant_uptake_vr(ids_H1PO4B,L)=trcs_plant_uptake_vr(ids_H1PO4B,L)+RUPH1B(N,L,NZ)
!
!     TOTAL ROOT C,N,P EXUDATION
!
!     TDFOME=total nonstructl C,N,P exchange
!     RDFOMC,RDFOMN,RDFOMP=PFT nonstructl C,N,P exchange
!
      DO K=1,jcplx
        DO NE=1,NumOfPlantChemElmnts
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
    PlantElemntStoreLandscape => plt_site%PlantElemntStoreLandscape  , &
    ElmntBalanceCum_pft  => plt_site%ElmntBalanceCum_pft   , &
    NH3EmiCum_pft => plt_bgcr%NH3EmiCum_pft  , &
    RNH3C => plt_bgcr%RNH3C  , &
    Canopy_NEE_col => plt_bgcr%Canopy_NEE_col  , &
    LitterFallChemElmnt_col => plt_bgcr%LitterFallChemElmnt_col  , &
    RootGasLoss_disturb => plt_bgcr%RootGasLoss_disturb, &
    RootN2Fix_pvr => plt_bgcr%RootN2Fix_pvr  , &
    CO2NetFix_pft  => plt_bgcr%CO2NetFix_pft   , &
    ETCanP => plt_ew%ETCanP    , &
    TH2GZ => plt_bgcr%TH2GZ  , &
    trcs_plant_uptake_vr => plt_rbgc%trcs_plant_uptake_vr, &    
    RNH3B => plt_rbgc%RNH3B  , &
    PlantRootSoilChemNetX_pft=> plt_rbgc%PlantRootSoilChemNetX_pft , &
    TRootGasLoss_disturb => plt_rbgc%TRootGasLoss_disturb  , &
    PTrans    => plt_ew%PTrans       , &
    PrecIntcptByCanP  => plt_ew%PrecIntcptByCanP     , &
    VapXAir2PCan => plt_ew%VapXAir2PCan    , &
    WatByPCanopy => plt_ew%WatByPCanopy    , &
    CanWatP => plt_ew%CanWatP    , &
    Eco_Heat_Grnd_col   => plt_ew%Eco_Heat_Grnd_col      , &
    HeatXAir2PCan => plt_ew%HeatXAir2PCan    , &
    EvapTransHeatP => plt_ew%EvapTransHeatP    , &
    CanWatg=> plt_ew%CanWatg   , &
    TKC   => plt_ew%TKC      , &
    TKS   => plt_ew%TKS      , &
    ENGYX => plt_ew%ENGYX    , &
    Eco_Heat_Sens_col   => plt_ew%Eco_Heat_Sens_col      , &
    VapXAir2CanG=> plt_ew%VapXAir2CanG   , &
    TENGYC=> plt_ew%TENGYC   , &
    TEVAPP=> plt_ew%TEVAPP   , &
    THFLXC=> plt_ew%THFLXC   , &
    LWRadCanG => plt_ew%LWRadCanG    , &
    TairK   => plt_ew%TairK      , &
    HeatStorCanP => plt_ew%HeatStorCanP    , &
    Eco_Heat_Latent_col   => plt_ew%Eco_Heat_Latent_col      , &
    CanH2OHeldVg=> plt_ew%CanH2OHeldVg   , &
    NU    => plt_site%NU     , &
    StemAreag => plt_morph%StemAreag , &
    CanopyLA_grd => plt_morph%CanopyLA_grd , &
    NI    => plt_morph%NI    , &
    NumOfBranches_pft   => plt_morph%NumOfBranches_pft   , &
    CanopyStemA_pft => plt_morph%CanopyStemA_pft , &
    CanopyLeafA_pft => plt_morph%CanopyLeafA_pft , &
    RadNet2CanP  => plt_rad%RadNet2CanP    , &
    LWRadCanP => plt_rad%LWRadCanP   , &
    Eco_NetRad_col   => plt_rad%Eco_NetRad_col       &
  )
  DO L=NU,NI(NZ)
    trcs_plant_uptake_vr(idg_N2,L)=trcs_plant_uptake_vr(idg_N2,L)+RootN2Fix_pvr(L,NZ)
  ENDDO
!
!     TOTAL ENERGY, WATER, CO2 FLUXES
!
!     Eco_NetRad_col=total net SW+LW absorbed by canopy
!     RadNet2CanP=PFT net SW+LW absorbed by canopy
!     Eco_Heat_Latent_col=total canopy latent heat flux
!     EvapTransHeatP=PFT canopy latent heat flux
!     Eco_Heat_Sens_col=total canopy sensible heat flux
!     HeatXAir2PCan=PFT canopy sensible heat flux
!     Eco_Heat_Grnd_col=total canopy storage heat flux
!     HeatStorCanP=PFT canopy storage heat flux
!     Canopy_NEE_col=total net CO2 fixation
!     CO2NetFix_pft=PFT net CO2 fixation
!     CanWatg,CanH2OHeldVg=total water volume in canopy,on canopy surfaces
!     CanWatP,WatByPCanopy=PFT water volume in canopy,on canopy surfaces
!     TEVAPP,VapXAir2CanG=total water flux to,from canopy,canopy surfaces
!     VapXAir2PCan,PTrans=water flux to,from canopy surfaces, inside canopy
!     TENGYC=total canopy water heat content
!     ENGYC=PFT canopy water heat content
!     CanopyLA_grd,StemAreag=total leaf,stalk area
!     CanopyLeafA_pft,CanopyStemA_pft=PFT leaf,stalk area
!     ZCSNC,ZZSNC,ZPSNC=total net root-soil C,N,P exchange
!     HCUPTK,HZUPTK,HPUPTK=PFT net root-soil C,N,P exchange
!     TBALC,TBALN,TBALP=total C,N,P balance
!     BALC,BALN,BALP=PFT C,N,P balance
!     TRootGasLoss_disturb=total loss of root CO2, O2, CH4, N2O, NH3, H2
!     RootGasLoss_disturb=PFT loss of root CO2, O2, CH4, N2O, NH3, H2
!
  Eco_NetRad_col=Eco_NetRad_col+RadNet2CanP(NZ)
  Eco_Heat_Latent_col=Eco_Heat_Latent_col+EvapTransHeatP(NZ)
  Eco_Heat_Sens_col=Eco_Heat_Sens_col+HeatXAir2PCan(NZ)
  Eco_Heat_Grnd_col=Eco_Heat_Grnd_col+HeatStorCanP(NZ)
  Canopy_NEE_col=Canopy_NEE_col+CO2NetFix_pft(NZ)
  ETCanP(NZ)=ETCanP(NZ)+PTrans(NZ)+VapXAir2PCan(NZ)
  CanWatg=CanWatg+CanWatP(NZ)
  CanH2OHeldVg=CanH2OHeldVg+WatByPCanopy(NZ)
  TEVAPP=TEVAPP+PTrans(NZ)+VapXAir2PCan(NZ)
  VapXAir2CanG=VapXAir2CanG+VapXAir2PCan(NZ)
  ENGYC=cpw*(WatByPCanopy(NZ)+PrecIntcptByCanP(NZ)+VapXAir2PCan(NZ))*TKC(NZ)
  TENGYC=TENGYC+ENGYC
  THFLXC=THFLXC+ENGYC-ENGYX(NZ)-(PrecIntcptByCanP(NZ)*cpw*TairK)
  ENGYX(NZ)=ENGYC
  LWRadCanG=LWRadCanG+LWRadCanP(NZ)
  CanopyLA_grd=CanopyLA_grd+CanopyLeafA_pft(NZ)
  StemAreag=StemAreag+CanopyStemA_pft(NZ)
  DO NE=1,NumOfPlantChemElmnts
    LitterFallChemElmnt_col(NE)=LitterFallChemElmnt_col(NE)-PlantRootSoilChemNetX_pft(NE,NZ)
    PlantElemntStoreLandscape(NE)=PlantElemntStoreLandscape(NE)+ElmntBalanceCum_pft(NE,NZ)
  ENDDO

  DO NTG=idg_beg,idg_end-1
    TRootGasLoss_disturb(NTG)=TRootGasLoss_disturb(NTG)+RootGasLoss_disturb(NTG,NZ)
  ENDDO
!
!     TOTAL CANOPY NH3 EXCHANGE AND EXUDATION
!
!     RNH3B,RNH3C=PFT NH3 flux between atmosphere and branch,canopy
!     NH3EmiCum_pft=total NH3 flux between atmosphere and canopy
!
  RNH3C(NZ)=0._r8
  DO NB=1,NumOfBranches_pft(NZ)
    RNH3C(NZ)=RNH3C(NZ)+RNH3B(NB,NZ)
    NH3EmiCum_pft(NZ)=NH3EmiCum_pft(NZ)+RNH3B(NB,NZ)
  ENDDO
  end associate
  end subroutine CanopyFluxesandFixation

  end module ExtractsMod
