module ExtractsMod
!!
!Description:
!     THIS SUBROUTINE AGGREGATES ALL SOIL-PLANT C,N,P EXCHANGES
!     FROM 'UPTAKE' AMD 'GROSUB' AND SENDS RESULTS TO 'REDIST'
!
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: AZMAX1
  use PlantBalMod,   only: SumPlantRootGas
  use PlantDebugMod, only: PrintRootTracer
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

    call PrintRootTracer(I,J,NZ,'extract')
    IF(plt_pheno%IsPlantActive_pft(NZ).EQ.iActive)THEN

      call CalcTotalLeafArea(NZ)

      call TotalGasandSoluteUptake(I,J,NZ)

      call ExtractCanopyFluxes(I,J,NZ)

    ENDIF
  ENDDO

  RETURN
  END subroutine extracts
!------------------------------------------------------------------------------------------

  subroutine TotalLitrFall()

  implicit none
  integer :: NZ,L,K,M
  integer :: NE
  associate(                                                         &
    LitrfalStrutElms_pft      => plt_bgcr%LitrfalStrutElms_pft,      &  !input
    LitrfalStrutElms_pvr      => plt_bgcr%LitrfalStrutElms_pvr,      &  !input
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft,         &  !input
    NP0                       => plt_site%NP0,                       &  !input
    StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft,    &  !input
    LitrFallStrutElms_col     => plt_bgcr%LitrFallStrutElms_col,     &  !inoput
    LitrfalStrutElms_vr       => plt_bgcr%LitrfalStrutElms_vr,       &  !inoput
    StandingDeadStrutElms_col => plt_biom%StandingDeadStrutElms_col, &  !inoput
    CanopyLeafAareZ_col       => plt_morph%CanopyLeafAareZ_col,      &  !output
    CanopyLeafArea_col        => plt_morph%CanopyLeafArea_col ,      &  !output
    CanopyStemAareZ_col       => plt_morph%CanopyStemAareZ_col,      &  !output
    StemArea_col              => plt_morph%StemArea_col,             &  !output
    tCanLeafC_clyr            => plt_biom%tCanLeafC_clyr             &  !output
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
     CanopyLeafArea_col = 0._r8
     StemArea_col       = 0._r8
  DO L                  = 1, NumOfCanopyLayers1
    CanopyLeafAareZ_col(L) = 0._r8
    tCanLeafC_clyr(L)        = 0._r8
    CanopyStemAareZ_col(L) = 0._r8
  ENDDO
  end associate
  end subroutine TotalLitrFall
!------------------------------------------------------------------------------------------

  subroutine CalcTotalLeafArea(NZ)
!
!     TOTAL LEAF AREA OF ALL PLANT SPECIES
!
!     CanopyLeafAareZ_col,CanopyStemAareZ_col=total leaf,stalk area of combined canopy layer
!     CanopyLeafAreaZ_pft,CanopyStemAreaZ_pft=PFT leaf,stalk area in canopy layer
!     tCanLeafC_clyr=total leaf C of combined canopy layer
!     CanopyLeafCLyr_pft=PFT leaf C in canopy layer
!
  implicit none
  integer, intent(in) :: NZ
  integer :: L
  associate(                                              &
    CanopyLeafAreaZ_pft => plt_morph%CanopyLeafAreaZ_pft, &  !input
    CanopyLeafCLyr_pft  => plt_biom%CanopyLeafCLyr_pft,   &  !input
    CanopyStemAreaZ_pft => plt_morph%CanopyStemAreaZ_pft, &  !input
    CanopyLeafAareZ_col => plt_morph%CanopyLeafAareZ_col, &  !inoput
    CanopyStemAareZ_col => plt_morph%CanopyStemAareZ_col, &  !inoput
    tCanLeafC_clyr      => plt_biom%tCanLeafC_clyr        &  !inoput
  )
  DO L=1,NumOfCanopyLayers1
    CanopyLeafAareZ_col(L)=CanopyLeafAareZ_col(L)+CanopyLeafAreaZ_pft(L,NZ)
    tCanLeafC_clyr(L)=tCanLeafC_clyr(L)+CanopyLeafCLyr_pft(L,NZ)
    CanopyStemAareZ_col(L)=CanopyStemAareZ_col(L)+CanopyStemAreaZ_pft(L,NZ)
  ENDDO
  end associate
  end subroutine CalcTotalLeafArea
!------------------------------------------------------------------------------------------

  subroutine TotalGasandSoluteUptake(I,J,NZ)
!
!     TOTAL GAS AND SOLUTE UPTAKE BY ALL PLANT SPECIES
!

  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: N,L,K,idg,NE,ids
 
  associate(                                                        &
    AREA3                     => plt_site%AREA3,                    &  !input
    AllPlantRootH2OLoss_pvr   => plt_ew%AllPlantRootH2OLoss_pvr,    &  !input
    Myco_pft                    => plt_morph%Myco_pft,                  &  !input
    MaxNumRootLays            => plt_site%MaxNumRootLays,           &  !input
    NU                        => plt_site%NU,                       &  !input
    PlantPopulation_pft       => plt_site%PlantPopulation_pft,      &  !input
    RCO2Emis2Root_pvr         => plt_rbgc%RCO2Emis2Root_pvr,        &  !input
    RootH1PO4DmndBand_pvr     => plt_rbgc%RootH1PO4DmndBand_pvr,    &  !input
    RootH1PO4DmndSoil_pvr     => plt_rbgc%RootH1PO4DmndSoil_pvr,    &  !input
    RootH2PO4DmndBand_pvr     => plt_rbgc%RootH2PO4DmndBand_pvr,    &  !input
    RootH2PO4DmndSoil_pvr     => plt_rbgc%RootH2PO4DmndSoil_pvr,    &  !input
    RootLenDensPerPlant_pvr   => plt_morph%RootLenDensPerPlant_pvr, &  !input
    RootMycoExudEUptk_pvr     => plt_rbgc%RootMycoExudEUptk_pvr,    &  !input
    RootNH4DmndBand_pvr       => plt_rbgc%RootNH4DmndBand_pvr,      &  !input
    RootNH4DmndSoil_pvr       => plt_rbgc%RootNH4DmndSoil_pvr,      &  !input
    RootNO3DmndBand_pvr       => plt_rbgc%RootNO3DmndBand_pvr,      &  !input
    RootNO3DmndSoil_pvr       => plt_rbgc%RootNO3DmndSoil_pvr,      &  !input
    RootNutUptake_pvr         => plt_rbgc%RootNutUptake_pvr,        &  !input
    RootO2Dmnd4Resp_pvr       => plt_rbgc%RootO2Dmnd4Resp_pvr,      &  !input
    RootO2Uptk_pvr            => plt_rbgc%RootO2Uptk_pvr,           &  !input
    RootO2_Xink_pvr           => plt_bgcr%RootO2_Xink_pvr,          &  !input
    RootUptkSoiSol_pvr        => plt_rbgc%RootUptkSoiSol_pvr,       &  !input
    TKCanopy_pft              => plt_ew%TKCanopy_pft,               &  !input
    TKS_vr                    => plt_ew%TKS_vr,                     &  !input
    trcg_air2root_flx_pvr     => plt_rbgc%trcg_air2root_flx_pvr,    &  !input
    REcoH1PO4DmndBand_vr      => plt_bgcr%REcoH1PO4DmndBand_vr,     &  !inoput
    REcoH1PO4DmndSoil_vr      => plt_bgcr%REcoH1PO4DmndSoil_vr,     &  !inoput
    REcoH2PO4DmndBand_vr      => plt_bgcr%REcoH2PO4DmndBand_vr,     &  !inoput
    REcoH2PO4DmndSoil_vr      => plt_bgcr%REcoH2PO4DmndSoil_vr,     &  !inoput
    REcoNH4DmndBand_vr        => plt_bgcr%REcoNH4DmndBand_vr,       &  !inoput
    REcoNH4DmndSoil_vr        => plt_bgcr%REcoNH4DmndSoil_vr,       &  !inoput
    REcoNO3DmndBand_vr        => plt_bgcr%REcoNO3DmndBand_vr,       &  !inoput
    REcoNO3DmndSoil_vr        => plt_bgcr%REcoNO3DmndSoil_vr,       &  !inoput
    REcoO2DmndResp_vr         => plt_bgcr%REcoO2DmndResp_vr,        &  !inoput
    RUptkRootO2_vr            => plt_bgcr%RUptkRootO2_vr,           &  !inoput
    RootCO2Emis2Root_vr       => plt_bgcr%RootCO2Emis2Root_vr,      &  !inoput
    RootO2_Xink_vr            => plt_bgcr%RootO2_Xink_vr,           &  !inoput
    THeatLossRoot2Soil_vr     => plt_ew%THeatLossRoot2Soil_vr,      &  !inoput
    TWaterPlantRoot2Soil_vr   => plt_ew%TWaterPlantRoot2Soil_vr,    &  !inoput
    tRootMycoExud2Soil_vr     => plt_bgcr%tRootMycoExud2Soil_vr,    &  !inoput
    totRootLenDens_vr         => plt_morph%totRootLenDens_vr,       &  !inoput
    trcg_air2root_flx_vr      => plt_rbgc%trcg_air2root_flx_vr,     &  !inoput
    trcs_Soil2plant_uptake_vr => plt_rbgc%trcs_Soil2plant_uptake_vr &  !inoput
  )
  
  DO L=NU,MaxNumRootLays
    DO N=1,Myco_pft(NZ)  
!
!     TOTAL ROOT DENSITY
!
!     totRootLenDens_vr=total root length density
!     RootLenDensPerPlant_pvr=PFT root length density per plant
!     AllPlantRootH2OLoss_pvr=total water uptake
!     AllPlantRootH2OLoss_pvr=PFT root water uptake
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
      TWaterPlantRoot2Soil_vr(L) = TWaterPlantRoot2Soil_vr(L)+AllPlantRootH2OLoss_pvr(N,L,NZ)

      !water lose from canopy to soil
      if(AllPlantRootH2OLoss_pvr(N,L,NZ)>0._r8)then
        THeatLossRoot2Soil_vr(L)     = THeatLossRoot2Soil_vr(L)+AllPlantRootH2OLoss_pvr(N,L,NZ)*cpw*TKCanopy_pft(NZ)
      !water lose from soil to canopy  
      else
        THeatLossRoot2Soil_vr(L)     = THeatLossRoot2Soil_vr(L)+AllPlantRootH2OLoss_pvr(N,L,NZ)*cpw*TKS_vr(L)
      endif
!
!     TOTAL ROOT BOUNDARY GAS FLUXES
!
      DO idg=idg_beg,idg_NH3
        trcg_air2root_flx_vr(idg,L)=trcg_air2root_flx_vr(idg,L)+trcg_air2root_flx_pvr(idg,N,L,NZ)
      ENDDO
      RootO2_Xink_vr(L)      = RootO2_Xink_vr(L) + RootO2_Xink_pvr(N,L,NZ)
      RootCO2Emis2Root_vr(L) = RootCO2Emis2Root_vr(L)+RCO2Emis2Root_pvr(N,L,NZ)
      RUptkRootO2_vr(L)      = RUptkRootO2_vr(L)+RootO2Uptk_pvr(N,L,NZ)

      !(>0) uptake from soil into roots
      DO idg=idg_beg,idg_end
        trcs_Soil2plant_uptake_vr(idg,L)=trcs_Soil2plant_uptake_vr(idg,L)+RootUptkSoiSol_pvr(idg,N,L,NZ)
      ENDDO

      !Nutrients are sorted according to band|soil
      do ids=ids_NH4B,ids_nuts_end
        trcs_Soil2plant_uptake_vr(ids,L)=trcs_Soil2plant_uptake_vr(ids,L)+RootNutUptake_pvr(ids,N,L,NZ)
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

  subroutine ExtractCanopyFluxes(I,J,NZ)
!
!     TOTAL ROOT N2 FIXATION BY ALL PLANT SPECIES
!
!     TRootN2Fix_pft=total root N2 fixation
!     RootN2Fix_pvr=PFT root N2 fixation
!

  implicit none
  integer, intent(in) :: I,J,NZ
  integer :: L, NE,NB,idg
  real(r8) :: ENGYC

  associate(                                                         &
    CO2NetFix_pft             => plt_bgcr%CO2NetFix_pft,             &  !input
    CanopyBiomWater_pft       => plt_ew%CanopyBiomWater_pft,         &  !input
    CanopyLeafArea_pft        => plt_morph%CanopyLeafArea_pft,       &  !input
    CanopyStemArea_pft        => plt_morph%CanopyStemArea_pft,       &  !input
    ElmBalanceCum_pft         => plt_site%ElmBalanceCum_pft,         &  !input
    EvapTransLHeat_pft        => plt_ew%EvapTransLHeat_pft,          &  !input
    HeatStorCanopy_pft        => plt_ew%HeatStorCanopy_pft,          &  !input
    HeatXAir2PCan_pft         => plt_ew%HeatXAir2PCan_pft,           &  !input
    LWRadCanopy_pft           => plt_rad%LWRadCanopy_pft,            &  !input
    NH3Dep2Can_pft            => plt_bgcr%NH3Dep2Can_pft,            &  !input
    PlantRootSoilElmNetX_pft  => plt_rbgc%PlantRootSoilElmNetX_pft,  &  !input
    RadNet2Canopy_pft         => plt_rad%RadNet2Canopy_pft,          &  !input
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft,    &  !input
    TKC_pft                   => plt_ew%TKC_pft,                     &  !input
    Transpiration_pft         => plt_ew%Transpiration_pft,           &  !input
    VHeatCapCanopy_pft        => plt_ew%VHeatCapCanopy_pft,          &  !input
    VapXAir2Canopy_pft        => plt_ew%VapXAir2Canopy_pft,          &  !input
    WatHeldOnCanopy_pft       => plt_ew%WatHeldOnCanopy_pft,         &  !input
    CanopyHeatStor_col        => plt_ew%CanopyHeatStor_col,          &  !inoput
    CanopyLeafArea_col        => plt_morph%CanopyLeafArea_col,       &  !inoput
    CanopyWat_col             => plt_ew%CanopyWat_col,               &  !inoput
    Canopy_NEE_col            => plt_bgcr%Canopy_NEE_col,            &  !inoput
    ENGYX_pft                 => plt_ew%ENGYX_pft,                   &  !inoput
    ETCanopy_CumYr_pft        => plt_ew%ETCanopy_CumYr_pft,          &  !inoput
    Eco_Heat_GrndSurf_col     => plt_ew%Eco_Heat_GrndSurf_col,       &  !inoput
    Eco_Heat_Latent_col       => plt_ew%Eco_Heat_Latent_col,         &  !inoput
    Eco_Heat_Sens_col         => plt_ew%Eco_Heat_Sens_col,           &  !inoput
    Eco_NetRad_col            => plt_rad%Eco_NetRad_col,             &  !inoput
    HeatFlx2Canopy_col        => plt_ew%HeatFlx2Canopy_col,          &  !inoput
    LWRadCanG                 => plt_ew%LWRadCanG,                   &  !inoput
    LitrFallStrutElms_col     => plt_bgcr%LitrFallStrutElms_col,     &  !inoput
    NH3Emis_CumYr_pft         => plt_bgcr%NH3Emis_CumYr_pft,         &  !inoput
    PlantElemntStoreLandscape => plt_site%PlantElemntStoreLandscape, &  !inoput
    QVegET_col                => plt_ew%QVegET_col,                  &  !inoput
    StemArea_col              => plt_morph%StemArea_col,             &  !inoput
    TRootGasLossDisturb_col   => plt_rbgc%TRootGasLossDisturb_col,   &  !inoput
    VapXAir2Canopy_col        => plt_ew%VapXAir2Canopy_col,          &  !inoput
    WatHeldOnCanopy_col       => plt_ew%WatHeldOnCanopy_col          &  !inoput
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
!     Canopy_NEE_col=total net CO2 fixation in the column
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
  end subroutine ExtractCanopyFluxes

  end module ExtractsMod
