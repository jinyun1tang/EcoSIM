module ExtractsMod
!!
!Description:
!     THIS SUBROUTINE AGGREGATES ALL SOIL-PLANT C,N,P EXCHANGES
!     FROM 'UPTAKE' AMD 'GROSUB' AND SENDS RESULTS TO 'REDIST'
!
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: AZMAX1
  use PlantBalMod, only : SumPlantRootGas
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

  call SumPlantRootGas(I,J)

  DO NZ=1,plt_site%NP
    IF(plt_pheno%IsPlantActive_pft(NZ).EQ.iActive)THEN

      call TotalLeafArea(NZ)

      call TotalGasandSoluteUptake(I,J,NZ)

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
  associate(                                                        &
   NP0                       => plt_site%NP0,                       &
   tCanLeafC_cl              => plt_biom%tCanLeafC_cl,              &
   StandingDeadStrutElms_col => plt_biom%StandingDeadStrutElms_col, &
   StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft,    &
   LitrfalStrutElms_pft      => plt_bgcr%LitrfalStrutElms_pft,      &
   LitrFallStrutElms_col     => plt_bgcr%LitrFallStrutElms_col,     &
   LitrfalStrutElms_vr       => plt_bgcr%LitrfalStrutElms_vr,       &
   LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr,      &
   MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,         &
   CanopyStemAareZ_col       => plt_morph%CanopyStemAareZ_col,      &
   CanopyLeafAareZ_col       => plt_morph%CanopyLeafAareZ_col,      &
   StemArea_col              => plt_morph%StemArea_col,             &
   CanopyLeafArea_col        => plt_morph%CanopyLeafArea_col        &
  )
  DO NZ=1,NP0
!
!   TOTAL LitrFall OF ALL PLANT SPECIES
!
!   ZCSNC,ZZSNC,ZPSNC=total C,N,P LitrFall
!   HCSNC,HZSNC,HPSNC=hourly PFT C,N,P LitrFall from grosub.f
!   StandingDeadStrutElms_col=total standing dead C,N,P mass
!   WTSTG=PFT standing dead C,N,P mass
!   LitrfalStrutElms_pvr,=cumulative PFT C,N,P LitrFall from grosub.f
!   LitrfalStrutElms_vr,=cumulative total C,N,P LitrFall
!
    DO NE=1,NumPlantChemElms
      LitrFallStrutElms_col(NE)=LitrFallStrutElms_col(NE)+LitrfalStrutElms_pft(NE,NZ)
      StandingDeadStrutElms_col(NE)=StandingDeadStrutElms_col(NE)+StandDeadStrutElms_pft(NE,NZ)
    ENDDO

    DO  L=0,MaxSoiL4Root_pft(NZ)
      DO K=1,pltpar%NumOfPlantLitrCmplxs
        DO NE=1,NumPlantChemElms
          DO  M=1,pltpar%jsken
            LitrfalStrutElms_vr(NE,M,K,L)=LitrfalStrutElms_vr(NE,M,K,L)+AZMAX1(LitrfalStrutElms_pvr(NE,M,K,L,NZ))
            if(LitrfalStrutElms_vr(NE,M,K,L)<0._r8)then
            write(*,*)'extract',K,L,LitrfalStrutElms_vr(NE,M,K,L),LitrfalStrutElms_pvr(NE,M,K,L,NZ)
            endif
          enddo
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CanopyLeafArea_col=0._r8
  StemArea_col=0._r8
  DO  L=1,NumOfCanopyLayers1
    CanopyLeafAareZ_col(L)=0._r8
    tCanLeafC_cl(L)=0._r8
    CanopyStemAareZ_col(L)=0._r8
  ENDDO
  end associate
  end subroutine TotalLitrFall
!------------------------------------------------------------------------------------------

  subroutine TotalLeafArea(NZ)
!
!     TOTAL LEAF AREA OF ALL PLANT SPECIES
!
!     CanopyLeafAareZ_col,CanopyStemAareZ_col=total leaf,stalk area of combined canopy layer
!     CanopyLeafAreaZ_pft,CanopyStemAreaZ_pft=PFT leaf,stalk area in canopy layer
!     tCanLeafC_cl=total leaf C of combined canopy layer
!     CanopyLeafCLyr_pft=PFT leaf C in canopy layer
!
  implicit none
  integer, intent(in) :: NZ
  integer :: L
  associate(                                              &
    CanopyLeafCLyr_pft  => plt_biom%CanopyLeafCLyr_pft,   &
    tCanLeafC_cl        => plt_biom%tCanLeafC_cl,         &
    CanopyLeafAareZ_col => plt_morph%CanopyLeafAareZ_col, &
    CanopyStemAreaZ_pft => plt_morph%CanopyStemAreaZ_pft, &
    CanopyStemAareZ_col => plt_morph%CanopyStemAareZ_col, &
    CanopyLeafAreaZ_pft => plt_morph%CanopyLeafAreaZ_pft  &
  )
  DO L=1,NumOfCanopyLayers1
    CanopyLeafAareZ_col(L)=CanopyLeafAareZ_col(L)+CanopyLeafAreaZ_pft(L,NZ)
    tCanLeafC_cl(L)=tCanLeafC_cl(L)+CanopyLeafCLyr_pft(L,NZ)
    CanopyStemAareZ_col(L)=CanopyStemAareZ_col(L)+CanopyStemAreaZ_pft(L,NZ)
  ENDDO
  end associate
  end subroutine TotalLeafArea
!------------------------------------------------------------------------------------------

  subroutine TotalGasandSoluteUptake(I,J,NZ)
!
!     TOTAL GAS AND SOLUTE UPTAKE BY ALL PLANT SPECIES
!

  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: N,L,K,idg,NE,ids
 
  associate(                                                       &
    NU                       => plt_site%NU,                       &
    NK                       => plt_site%NK,                       &    
    MaxNumRootLays           => plt_site%MaxNumRootLays         ,  &
    AREA3                    => plt_site%AREA3,                    &
    PlantPopulation_pft      => plt_site%PlantPopulation_pft,      &
    RootH1PO4DmndBand_pvr    => plt_rbgc%RootH1PO4DmndBand_pvr,    &
    RootH2PO4DmndBand_pvr    => plt_rbgc%RootH2PO4DmndBand_pvr,    &
    RootNO3DmndBand_pvr      => plt_rbgc%RootNO3DmndBand_pvr,      &
    RootO2Uptk_pvr           => plt_rbgc%RootO2Uptk_pvr,           &
    trcg_Root_gas2aqu_flx_vr => plt_rbgc%trcg_Root_gas2aqu_flx_vr, &
    trcg_air2root_flx_pvr    => plt_rbgc%trcg_air2root_flx_pvr,    &
    RootCO2Emis_pvr          => plt_rbgc%RootCO2Emis_pvr,          &
    RootUptkSoiSol_pvr       => plt_rbgc%RootUptkSoiSol_pvr,       &
    RootNutUptake_pvr        => plt_rbgc%RootNutUptake_pvr,        &
    trcg_air2root_flx_vr     => plt_rbgc%trcg_air2root_flx_vr,     &
    RootO2Dmnd4Resp_pvr      => plt_rbgc%RootO2Dmnd4Resp_pvr,      &
    RootMycoExudEUptk_pvr    => plt_rbgc%RootMycoExudEUptk_pvr,    &
    RootNH4DmndSoil_pvr      => plt_rbgc%RootNH4DmndSoil_pvr,      &
    RootNO3DmndSoil_pvr      => plt_rbgc%RootNO3DmndSoil_pvr,      &
    RootH2PO4DmndSoil_pvr    => plt_rbgc%RootH2PO4DmndSoil_pvr,    &
    RootNH4DmndBand_pvr      => plt_rbgc%RootNH4DmndBand_pvr,      &
    RootH1PO4DmndSoil_pvr    => plt_rbgc%RootH1PO4DmndSoil_pvr,    &
    trcs_plant_uptake_vr     => plt_rbgc%trcs_plant_uptake_vr,     &
    RootN2Fix_pvr            => plt_bgcr%RootN2Fix_pvr,            &
    REcoNO3DmndSoil_vr       => plt_bgcr%REcoNO3DmndSoil_vr,       &
    REcoNH4DmndSoil_vr       => plt_bgcr%REcoNH4DmndSoil_vr,       &
    REcoH2PO4DmndSoil_vr     => plt_bgcr%REcoH2PO4DmndSoil_vr,     &
    REcoNO3DmndBand_vr       => plt_bgcr%REcoNO3DmndBand_vr,       &
    REcoH1PO4DmndSoil_vr     => plt_bgcr%REcoH1PO4DmndSoil_vr,     &
    REcoNH4DmndBand_vr       => plt_bgcr%REcoNH4DmndBand_vr,       &
    REcoO2DmndResp_vr        => plt_bgcr%REcoO2DmndResp_vr,        &
    tRootMycoExud2Soil_vr    => plt_bgcr%tRootMycoExud2Soil_vr,    &
    RUptkRootO2_vr           => plt_bgcr%RUptkRootO2_vr,           &
    RootCO2Emis2Root_vr     => plt_bgcr%RootCO2Emis2Root_vr,     &
    REcoH2PO4DmndBand_vr     => plt_bgcr%REcoH2PO4DmndBand_vr,     &
    REcoH1PO4DmndBand_vr     => plt_bgcr%REcoH1PO4DmndBand_vr,     &
    TKCanopy_pft             => plt_ew%TKCanopy_pft,               &
    TKS_vr                   => plt_ew%TKS_vr,                     &
    THeatLossRoot2Soil_vr    => plt_ew%THeatLossRoot2Soil_vr,      &
    TPlantRootH2OLoss_vr     => plt_ew%TPlantRootH2OLoss_vr,       &
    AllPlantRootH2OLoss_vr   => plt_ew%AllPlantRootH2OLoss_vr,     &
    RootLenDensPerPlant_pvr  => plt_morph%RootLenDensPerPlant_pvr, &
    totRootLenDens_vr        => plt_morph%totRootLenDens_vr,       &
    MY                       => plt_morph%MY,                      &
    MaxSoiL4Root_pft         => plt_morph%MaxSoiL4Root_pft         &
  )
  
  DO L=NU,MaxNumRootLays
    DO N=1,MY(NZ)  
!
!     TOTAL ROOT DENSITY
!
!     totRootLenDens_vr=total root length density
!     RootLenDensPerPlant_pvr=PFT root length density per plant
!     AllPlantRootH2OLoss_vr=total water uptake
!     AllPlantRootH2OLoss_vr=PFT root water uptake
!     THeatLossRoot2Soil_vr=total convective heat in root water uptake
!     TKS=soil temperature
!     PP=PFT population, this is dynamic, and can goes to zero
!
      IF(N.EQ.ipltroot)THEN
        totRootLenDens_vr(L)=totRootLenDens_vr(L)+RootLenDensPerPlant_pvr(N,L,NZ)*PlantPopulation_pft(NZ)/AREA3(L)
      ENDIF
!
!     TOTAL WATER UPTAKE
!
      TPlantRootH2OLoss_vr(L) = TPlantRootH2OLoss_vr(L)+AllPlantRootH2OLoss_vr(N,L,NZ)

      !water lose from canopy to soil
      if(AllPlantRootH2OLoss_vr(N,L,NZ)>0._r8)then
        THeatLossRoot2Soil_vr(L)     = THeatLossRoot2Soil_vr(L)+AllPlantRootH2OLoss_vr(N,L,NZ)*cpw*TKCanopy_pft(NZ)
      !water lose from soil to canopy  
      else
        THeatLossRoot2Soil_vr(L)     = THeatLossRoot2Soil_vr(L)+AllPlantRootH2OLoss_vr(N,L,NZ)*cpw*TKS_vr(L)
      endif
!
!     TOTAL ROOT BOUNDARY GAS FLUXES
!
      DO idg=idg_beg,idg_NH3
        trcg_air2root_flx_vr(idg,L)=trcg_air2root_flx_vr(idg,L)+trcg_air2root_flx_pvr(idg,N,L,NZ)
      ENDDO

      RootCO2Emis2Root_vr(L) = RootCO2Emis2Root_vr(L)+RootCO2Emis_pvr(N,L,NZ)
      RUptkRootO2_vr(L)       = RUptkRootO2_vr(L)+RootO2Uptk_pvr(N,L,NZ)

      DO idg=idg_beg,idg_end
        trcs_plant_uptake_vr(idg,L)=trcs_plant_uptake_vr(idg,L)+RootUptkSoiSol_pvr(idg,N,L,NZ)
      ENDDO

      !Nutrients are sorted according to band|soil
      do ids=ids_NH4B,ids_nuts_end
        trcs_plant_uptake_vr(ids,L)=trcs_plant_uptake_vr(ids,L)+RootNutUptake_pvr(ids,N,L,NZ)
      ENDDO

!
!     TOTAL ROOT C,N,P EXUDATION
!
!     tRootMycoExud2Soil_vr=total nonstructl C,N,P exchange
!     RDFOMC,RDFOMN,RDFOMP=PFT nonstructl C,N,P exchange
!
      DO K=1,jcplx
        DO NE=1,NumPlantChemElms
          tRootMycoExud2Soil_vr(NE,K,L)=tRootMycoExud2Soil_vr(NE,K,L)-RootMycoExudEUptk_pvr(NE,N,K,L,NZ)
        ENDDO
      ENDDO
!
!     TOTAL ROOT O2, NH4, NO3, PO4 UPTAKE CONTRIBUTES TO
!     TOTAL ROOT + MICROBIAL UPTAKE USED TO CALCULATE
!     COMPETITION CONSTRAINTS
!
!
      REcoO2DmndResp_vr(L)    = REcoO2DmndResp_vr(L)+RootO2Dmnd4Resp_pvr(N,L,NZ)
      REcoNH4DmndSoil_vr(L)   = REcoNH4DmndSoil_vr(L)+RootNH4DmndSoil_pvr(N,L,NZ)
      REcoNO3DmndSoil_vr(L)   = REcoNO3DmndSoil_vr(L)+RootNO3DmndSoil_pvr(N,L,NZ)
      REcoH2PO4DmndSoil_vr(L) = REcoH2PO4DmndSoil_vr(L)+RootH2PO4DmndSoil_pvr(N,L,NZ)
      REcoH1PO4DmndSoil_vr(L) = REcoH1PO4DmndSoil_vr(L)+RootH1PO4DmndSoil_pvr(N,L,NZ)
      REcoNH4DmndBand_vr(L)   = REcoNH4DmndBand_vr(L)+RootNH4DmndBand_pvr(N,L,NZ)
      REcoNO3DmndBand_vr(L)   = REcoNO3DmndBand_vr(L)+RootNO3DmndBand_pvr(N,L,NZ)
      REcoH2PO4DmndBand_vr(L) = REcoH2PO4DmndBand_vr(L)+RootH2PO4DmndBand_pvr(N,L,NZ)
      REcoH1PO4DmndBand_vr(L) = REcoH1PO4DmndBand_vr(L)+RootH1PO4DmndBand_pvr(N,L,NZ)
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
  integer :: L, NE,NB,idg
  real(r8) :: ENGYC

  associate(                                                         &
    PlantElemntStoreLandscape => plt_site%PlantElemntStoreLandscape, &
    ElmBalanceCum_pft         => plt_site%ElmBalanceCum_pft,         &
    NH3Emis_CumYr_pft         => plt_bgcr%NH3Emis_CumYr_pft,         &
    Canopy_NEE_col            => plt_bgcr%Canopy_NEE_col,            &
    LitrFallStrutElms_col     => plt_bgcr%LitrFallStrutElms_col,     &
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft,    &
    RootN2Fix_pvr             => plt_bgcr%RootN2Fix_pvr,             &
    CO2NetFix_pft             => plt_bgcr%CO2NetFix_pft,             &
    ETCanopy_CumYr_pft        => plt_ew%ETCanopy_CumYr_pft,          &
    TRootH2Flx_col            => plt_bgcr%TRootH2Flx_col,            &
    PlantRootSoilElmNetX_pft  => plt_rbgc%PlantRootSoilElmNetX_pft,  &
    TRootGasLossDisturb_col   => plt_rbgc%TRootGasLossDisturb_col,   &
    Transpiration_pft         => plt_ew%Transpiration_pft,           &
    PrecIntcptByCanopy_pft    => plt_ew%PrecIntcptByCanopy_pft,      &
    VapXAir2Canopy_pft        => plt_ew%VapXAir2Canopy_pft,          &
    WatHeldOnCanopy_pft       => plt_ew%WatHeldOnCanopy_pft,         &
    VHeatCapCanopy_pft        => plt_ew%VHeatCapCanopy_pft,          &
    CanopyBiomWater_pft       => plt_ew%CanopyBiomWater_pft,         &
    Eco_Heat_GrndSurf_col     => plt_ew%Eco_Heat_GrndSurf_col,       &
    HeatXAir2PCan_pft         => plt_ew%HeatXAir2PCan_pft,           &
    EvapTransLHeat_pft        => plt_ew%EvapTransLHeat_pft,          &
    CanopyWat_col             => plt_ew%CanopyWat_col,               &
    TKC_pft                   => plt_ew%TKC_pft,                     &
    TKS_vr                    => plt_ew%TKS_vr,                      &
    ENGYX_pft                 => plt_ew%ENGYX_pft,                   &
    Eco_Heat_Sens_col         => plt_ew%Eco_Heat_Sens_col,           &
    VapXAir2Canopy_col        => plt_ew%VapXAir2Canopy_col,          &
    CanopyHeatStor_col        => plt_ew%CanopyHeatStor_col,          &
    QVegET_col                => plt_ew%QVegET_col,                  &
    HeatFlx2Canopy_col        => plt_ew%HeatFlx2Canopy_col,          &
    LWRadCanG                 => plt_ew%LWRadCanG,                   &
    TairK                     => plt_ew%TairK,                       &
    HeatStorCanopy_pft        => plt_ew%HeatStorCanopy_pft,          &
    Eco_Heat_Latent_col       => plt_ew%Eco_Heat_Latent_col,         &
    WatHeldOnCanopy_col       => plt_ew%WatHeldOnCanopy_col,         &
    NU                        => plt_site%NU,                        &
    NH3Dep2Can_pft            => plt_bgcr%NH3Dep2Can_pft,            &
    StemArea_col              => plt_morph%StemArea_col,             &
    CanopyLeafArea_col        => plt_morph%CanopyLeafArea_col,       &
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,         &
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft,        &
    CanopyStemArea_pft        => plt_morph%CanopyStemArea_pft,       &
    CanopyLeafArea_pft        => plt_morph%CanopyLeafArea_pft,       &
    RadNet2Canopy_pft         => plt_rad%RadNet2Canopy_pft,          &
    LWRadCanopy_pft           => plt_rad%LWRadCanopy_pft,            &
    Eco_NetRad_col            => plt_rad%Eco_NetRad_col              &
  )
!
!     TOTAL ENERGY, WATER, CO2 FLUXES
!
!     Eco_NetRad_col=total net SW+LW absorbed by canopy
!     RadNet2Canopy_pft=PFT net SW+LW absorbed by canopy
!     Eco_Heat_Latent_col=total canopy latent heat flux
!     EvapTransLHeat_pft=PFT canopy latent heat flux
!     Eco_Heat_Sens_col=total canopy sensible heat flux
!     HeatXAir2PCan_pft=PFT canopy sensible heat flux
!     Eco_Heat_GrndSurf_col=total canopy storage heat flux
!     HeatStorCanopy_pft=PFT canopy storage heat flux
!     Canopy_NEE_col=total net CO2 fixation
!     CO2NetFix_pft=PFT net CO2 fixation
!     CanopyWat_col,WatHeldOnCanopy=total water volume in canopy,on canopy surfaces
!     CanopyBiomWater_pft,WatHeldOnCanopy_pft=PFT water volume in canopy,on canopy surfaces
!     QVegET_col,VapXAir2Canopy_col=total water flux to,from canopy,canopy surfaces
!     VapXAir2PCan,Transpiration_pft=water flux to,from canopy surfaces, inside canopy
!     CanopyHeatStor_col=total canopy water heat content
!     ENGYC=PFT canopy water heat content
!     CanopyLeafArea_col,StemArea_col=total leaf,stalk area
!     CanopyLeafArea_pft,CanopyStemArea_pft=PFT leaf,stalk area
!     ZCSNC,ZZSNC,ZPSNC=total net root-soil C,N,P exchange
!     HCUPTK,HZUPTK,HPUPTK=PFT net root-soil C,N,P exchange
!     TBALC,TBALN,TBALP=total C,N,P balance
!     BALC,BALN,BALP=PFT C,N,P balance
!     TRootGasLossDisturb_col=total loss of root CO2, O2, CH4, N2O, NH3, H2
!     RootGasLossDisturb_pft=PFT loss of root CO2, O2, CH4, N2O, NH3, H2
!
  Eco_NetRad_col         = Eco_NetRad_col+RadNet2Canopy_pft(NZ)
  Eco_Heat_Latent_col    = Eco_Heat_Latent_col+EvapTransLHeat_pft(NZ)
  Eco_Heat_Sens_col      = Eco_Heat_Sens_col+HeatXAir2PCan_pft(NZ)
  Eco_Heat_GrndSurf_col  = Eco_Heat_GrndSurf_col+HeatStorCanopy_pft(NZ)
  Canopy_NEE_col         = Canopy_NEE_col+CO2NetFix_pft(NZ)
  ETCanopy_CumYr_pft(NZ) = ETCanopy_CumYr_pft(NZ)+Transpiration_pft(NZ)+VapXAir2Canopy_pft(NZ)
  CanopyWat_col          = CanopyWat_col+CanopyBiomWater_pft(NZ)
  WatHeldOnCanopy_col    = WatHeldOnCanopy_col+WatHeldOnCanopy_pft(NZ)
  QVegET_col             = QVegET_col+Transpiration_pft(NZ)+VapXAir2Canopy_pft(NZ)
  VapXAir2Canopy_col     = VapXAir2Canopy_col+VapXAir2Canopy_pft(NZ)
  ENGYC                  = VHeatCapCanopy_pft(NZ)*TKC_pft(NZ)
  CanopyHeatStor_col     = CanopyHeatStor_col+ENGYC
  HeatFlx2Canopy_col     = HeatFlx2Canopy_col+ENGYC-ENGYX_pft(NZ)
  ENGYX_pft(NZ)          = ENGYC
  LWRadCanG              = LWRadCanG+LWRadCanopy_pft(NZ)
  CanopyLeafArea_col     = CanopyLeafArea_col+CanopyLeafArea_pft(NZ)
  StemArea_col           = StemArea_col+CanopyStemArea_pft(NZ)

  DO NE=1,NumPlantChemElms
    LitrFallStrutElms_col(NE)     = LitrFallStrutElms_col(NE)-PlantRootSoilElmNetX_pft(NE,NZ)
    PlantElemntStoreLandscape(NE) = PlantElemntStoreLandscape(NE)+ElmBalanceCum_pft(NE,NZ)
  ENDDO

  DO idg=idg_beg,idg_NH3
    TRootGasLossDisturb_col(idg)=TRootGasLossDisturb_col(idg)+RootGasLossDisturb_pft(idg,NZ)
  ENDDO
!
!     TOTAL CANOPY NH3 EXCHANGE AND EXUDATION
!
!     NH3Emis_CumYr_pft=total NH3 flux between atmosphere and canopy
!
    
  NH3Emis_CumYr_pft(NZ)=NH3Emis_CumYr_pft(NZ)+NH3Dep2Can_pft(NZ)

  end associate
  end subroutine CanopyFluxesandFixation

  end module ExtractsMod
