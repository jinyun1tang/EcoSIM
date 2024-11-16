module ExtractsMod
!!
!Description:
!     THIS SUBROUTINE AGGREGATES ALL SOIL-PLANT C,N,P EXCHANGES
!     FROM 'UPTAKE' AMD 'GROSUB' AND SENDS RESULTS TO 'REDIST'
!
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: AZMAX1
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
    IF(plt_pheno%IsPlantActive_pft(NZ).EQ.iActive)THEN

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

  subroutine TotalGasandSoluteUptake(NZ)
!
!     TOTAL GAS AND SOLUTE UPTAKE BY ALL PLANT SPECIES
!

  implicit none
  integer, intent(in) :: NZ

  integer :: N,L,K,NTG,NE,idg,ids
 
  associate(                                                        &
    NU                       => plt_site%NU,                       &
    AREA3                    => plt_site%AREA3,                    &
    PlantPopulation_pft      => plt_site%PlantPopulation_pft,      &
    RootH1PO4DmndBand_pvr    => plt_rbgc%RootH1PO4DmndBand_pvr,    &
    RootH2PO4DmndBand_pvr    => plt_rbgc%RootH2PO4DmndBand_pvr,    &
    RootNO3DmndBand_pvr      => plt_rbgc%RootNO3DmndBand_pvr,      &
    RootO2Uptk_pvr           => plt_rbgc%RootO2Uptk_pvr,           &
    trcg_Root_gas2aqu_flx_vr => plt_rbgc%trcg_Root_gas2aqu_flx_vr, &
    trcg_air2root_flx_pvr    => plt_rbgc%trcg_air2root_flx_pvr,    &
    RootCO2Emis_pvr          => plt_rbgc%RootCO2Emis_pvr,          &
    RootUptkSoiSol_vr        => plt_rbgc%RootUptkSoiSol_vr,        &
    RootNutUptake_pvr        => plt_rbgc%RootNutUptake_pvr,        &
    trcg_air2root_flx_vr     => plt_rbgc%trcg_air2root_flx_vr,     &
    trcg_root_vr             => plt_rbgc%trcg_root_vr,             &
    RootO2Dmnd4Resp_pvr      => plt_rbgc%RootO2Dmnd4Resp_pvr,      &
    RootMycoExudElm_pvr      => plt_rbgc%RootMycoExudElm_pvr,      &
    RootNH4DmndSoil_pvr      => plt_rbgc%RootNH4DmndSoil_pvr,      &
    RootNO3DmndSoil_pvr      => plt_rbgc%RootNO3DmndSoil_pvr,      &
    RootH2PO4DmndSoil_pvr    => plt_rbgc%RootH2PO4DmndSoil_pvr,    &
    RootNH4DmndBand_pvr      => plt_rbgc%RootNH4DmndBand_pvr,      &
    RootH1PO4DmndSoil_pvr    => plt_rbgc%RootH1PO4DmndSoil_pvr,    &
    trcs_plant_uptake_vr     => plt_rbgc%trcs_plant_uptake_vr,     &
    REcoNO3DmndSoil_vr       => plt_bgcr%REcoNO3DmndSoil_vr,       &
    REcoNH4DmndSoil_vr       => plt_bgcr%REcoNH4DmndSoil_vr,       &
    REcoH2PO4DmndSoil_vr     => plt_bgcr%REcoH2PO4DmndSoil_vr,     &
    REcoNO3DmndBand_vr       => plt_bgcr%REcoNO3DmndBand_vr,       &
    REcoH1PO4DmndSoil_vr     => plt_bgcr%REcoH1PO4DmndSoil_vr,     &
    REcoNH4DmndBand_vr       => plt_bgcr%REcoNH4DmndBand_vr,       &
    REcoO2DmndResp_vr        => plt_bgcr%REcoO2DmndResp_vr,        &
    tRootMycoExud2Soil_vr    => plt_bgcr%tRootMycoExud2Soil_vr,    &
    tRO2MicrbUptk_vr         => plt_bgcr%tRO2MicrbUptk_vr,         &
    tRootCO2Emis_vr          => plt_bgcr%tRootCO2Emis_vr,          &
    REcoH2PO4DmndBand_vr     => plt_bgcr%REcoH2PO4DmndBand_vr,     &
    REcoH1PO4DmndBand_vr     => plt_bgcr%REcoH1PO4DmndBand_vr,     &
    TKS_vr                   => plt_ew%TKS_vr,                     &
    THeatRootUptake_vr       => plt_ew%THeatRootUptake_vr,         &
    TPlantRootH2OUptake_vr   => plt_ew%TPlantRootH2OUptake_vr,     &
    AllPlantRootH2OUptake_vr => plt_ew%AllPlantRootH2OUptake_vr,   &
    trcg_rootml_pvr          => plt_rbgc%trcg_rootml_pvr,          &
    trcs_rootml_pvr          => plt_rbgc%trcs_rootml_pvr,          &
    RootLenDensPerPlant_pvr  => plt_morph%RootLenDensPerPlant_pvr, &
    totRootLenDens_vr        => plt_morph%totRootLenDens_vr,       &
    MY                       => plt_morph%MY,                      &
    MaxSoiL4Root_pft         => plt_morph%MaxSoiL4Root_pft         &
  )

  trcs_plant_uptake_vr=0._r8
  DO N=1,MY(NZ)
    DO L=NU,MaxSoiL4Root_pft(NZ)
!
!     TOTAL ROOT DENSITY
!
!     totRootLenDens_vr=total root length density
!     RootLenDensPerPlant_pvr=PFT root length density per plant
!     AllPlantRootH2OUptake_vr=total water uptake
!     AllPlantRootH2OUptake_vr=PFT root water uptake
!     THeatRootUptake_vr=total convective heat in root water uptake
!     TKS=soil temperature
!     PP=PFT population, this is dynamic, and can goes to zero
!
      IF(N.EQ.ipltroot)THEN
        totRootLenDens_vr(L)=totRootLenDens_vr(L)+RootLenDensPerPlant_pvr(N,L,NZ)*PlantPopulation_pft(NZ)/AREA3(L)
      ENDIF
!
!     TOTAL WATER UPTAKE
!
      TPlantRootH2OUptake_vr(L)=TPlantRootH2OUptake_vr(L)+AllPlantRootH2OUptake_vr(N,L,NZ)
      THeatRootUptake_vr(L)=THeatRootUptake_vr(L)+AllPlantRootH2OUptake_vr(N,L,NZ)*cpw*TKS_vr(L)
!
!     ROOT GAS CONTENTS FROM FLUXES IN 'UPTAKE'
!
!     *A,*P=PFT root gaseous, aqueous gas content
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     R*FLA=root gaseous-atmosphere CO2 exchange
!     R*DFA=root aqueous-gaseous CO2 exchange
!
      DO NTG=idg_beg,idg_NH3
        trcg_rootml_pvr(NTG,N,L,NZ) = trcg_rootml_pvr(NTG,N,L,NZ)-trcg_Root_gas2aqu_flx_vr(NTG,N,L,NZ)+trcg_air2root_flx_pvr(NTG,N,L,NZ)
        trcs_rootml_pvr(NTG,N,L,NZ) = trcs_rootml_pvr(NTG,N,L,NZ)+trcg_Root_gas2aqu_flx_vr(NTG,N,L,NZ)
      ENDDO

      trcs_rootml_pvr(idg_CO2,N,L,NZ) = trcs_rootml_pvr(idg_CO2,N,L,NZ)+RootCO2Emis_pvr(N,L,NZ)
      trcs_rootml_pvr(idg_O2,N,L,NZ)  = trcs_rootml_pvr(idg_O2,N,L,NZ) -RootO2Uptk_pvr(N,L,NZ)
      trcs_rootml_pvr(idg_CH4,N,L,NZ) = trcs_rootml_pvr(idg_CH4,N,L,NZ)+RootUptkSoiSol_vr(idg_CH4,N,L,NZ)
      trcs_rootml_pvr(idg_N2O,N,L,NZ) = trcs_rootml_pvr(idg_N2O,N,L,NZ)+RootUptkSoiSol_vr(idg_N2O,N,L,NZ)
      trcs_rootml_pvr(idg_H2,N,L,NZ)  = trcs_rootml_pvr(idg_H2,N,L,NZ) +RootUptkSoiSol_vr(idg_H2,N,L,NZ)
      trcs_rootml_pvr(idg_NH3,N,L,NZ) = trcs_rootml_pvr(idg_NH3,N,L,NZ)+RootUptkSoiSol_vr(idg_NH3,N,L,NZ)+RootUptkSoiSol_vr(idg_NH3B,N,L,NZ)
!
!     TOTAL ROOT GAS CONTENTS
!
!     TL*P=total root gas content, dissolved + gaseous phase
!     *A,*P=PFT root gaseous, aqueous gas content
!
      DO NTG=idg_beg,idg_end-1
        trcg_root_vr(NTG,L)=trcg_root_vr(NTG,L)+trcs_rootml_pvr(NTG,N,L,NZ)+trcg_rootml_pvr(NTG,N,L,NZ)
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
        trcg_air2root_flx_vr(NTG,L)=trcg_air2root_flx_vr(NTG,L)+trcg_air2root_flx_pvr(NTG,N,L,NZ)
      ENDDO

      tRootCO2Emis_vr(L) =tRootCO2Emis_vr(L)-RootCO2Emis_pvr(N,L,NZ)
      tRO2MicrbUptk_vr(L)=tRO2MicrbUptk_vr(L)+RootO2Uptk_pvr(N,L,NZ)
      DO idg=idg_beg,idg_end
        trcs_plant_uptake_vr(idg,L)=trcs_plant_uptake_vr(idg,L)+RootUptkSoiSol_vr(idg,N,L,NZ)
      ENDDO

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
          tRootMycoExud2Soil_vr(NE,K,L)=tRootMycoExud2Soil_vr(NE,K,L)-RootMycoExudElm_pvr(NE,N,K,L,NZ)
        ENDDO
      ENDDO
!
!     TOTAL ROOT O2, NH4, NO3, PO4 UPTAKE CONTRIBUTES TO
!     TOTAL ROOT + MICROBIAL UPTAKE USED TO CALCULATE
!     COMPETITION CONSTRAINTS
!
!     REcoO2DmndResp_vr=O2 demand by all microbial,root,myco populations
!     REcoNH4DmndSoil_vr=NH4 demand in non-band by all microbial,root,myco populations
!     REcoNO3DmndSoil_vr=NO3 demand in non-band by all microbial,root,myco populations
!     REcoH2PO4DmndSoil_vr=H2PO4 demand in non-band by all microbial,root,myco populations
!     REcoH1PO4DmndSoil_vr=HPO4 demand in non-band by all microbial,root,myco populations
!     REcoNH4DmndBand_vr=NH4 demand in band by all microbial,root,myco populations
!     REcoNO3DmndBand_vr=NO3 demand in band by all microbial,root,myco populations
!     REcoH2PO4DmndBand_vr=H2PO4 demand in band by all microbial,root,myco populations
!     REcoH1PO4DmndBand_vr=HPO4 demand in band by all microbial,root,myco populations
!     RootO2Dmnd4Resp_pvr=O2 demand by each root,myco population
!     RootNH4DmndSoil_pvr=NH4 demand in non-band by each root population
!     RootNO3DmndSoil_pvr=NO3 demand in non-band by each root population
!     RootH2PO4DmndSoil_pvr=H2PO4 demand in non-band by each root population
!     RootH1PO4DmndSoil_pvr=HPO4 demand in non-band by each root population
!     RootNH4DmndBand_pvr=NH4 demand in band by each root population
!     RUNNXB=NO3 demand in band by each root population
!     RootH2PO4DmndBand_pvr=H2PO4 demand in band by each root population
!     RootH1PO4DmndBand_pvr=HPO4 demand in band by each root population
!
      REcoO2DmndResp_vr(L)=REcoO2DmndResp_vr(L)+RootO2Dmnd4Resp_pvr(N,L,NZ)
      REcoNH4DmndSoil_vr(L)=REcoNH4DmndSoil_vr(L)+RootNH4DmndSoil_pvr(N,L,NZ)
      REcoNO3DmndSoil_vr(L)=REcoNO3DmndSoil_vr(L)+RootNO3DmndSoil_pvr(N,L,NZ)
      REcoH2PO4DmndSoil_vr(L)=REcoH2PO4DmndSoil_vr(L)+RootH2PO4DmndSoil_pvr(N,L,NZ)
      REcoH1PO4DmndSoil_vr(L)=REcoH1PO4DmndSoil_vr(L)+RootH1PO4DmndSoil_pvr(N,L,NZ)
      REcoNH4DmndBand_vr(L)=REcoNH4DmndBand_vr(L)+RootNH4DmndBand_pvr(N,L,NZ)
      REcoNO3DmndBand_vr(L)=REcoNO3DmndBand_vr(L)+RootNO3DmndBand_pvr(N,L,NZ)
      REcoH2PO4DmndBand_vr(L)=REcoH2PO4DmndBand_vr(L)+RootH2PO4DmndBand_pvr(N,L,NZ)
      REcoH1PO4DmndBand_vr(L)=REcoH1PO4DmndBand_vr(L)+RootH1PO4DmndBand_pvr(N,L,NZ)
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
    trcs_plant_uptake_vr      => plt_rbgc%trcs_plant_uptake_vr,      &
    PlantRootSoilElmNetX_pft  => plt_rbgc%PlantRootSoilElmNetX_pft,  &
    TRootGasLossDisturb_pft   => plt_rbgc%TRootGasLossDisturb_pft,   &
    Transpiration_pft         => plt_ew%Transpiration_pft,           &
    PrecIntcptByCanopy_pft    => plt_ew%PrecIntcptByCanopy_pft,      &
    VapXAir2Canopy_pft        => plt_ew%VapXAir2Canopy_pft,          &
    WatByPCanopy_pft          => plt_ew%WatByPCanopy_pft,            &
    CanopyWater_pft           => plt_ew%CanopyWater_pft,             &
    Eco_Heat_GrndSurf_col     => plt_ew%Eco_Heat_GrndSurf_col,       &
    HeatXAir2PCan_pft         => plt_ew%HeatXAir2PCan_pft,           &
    EvapTransHeat_pft         => plt_ew%EvapTransHeat_pft,           &
    CanWat_col                => plt_ew%CanWat_col,                  &
    TKC                       => plt_ew%TKC,                         &
    TKS_vr                    => plt_ew%TKS_vr,                      &
    ENGYX_pft                 => plt_ew%ENGYX_pft,                   &
    Eco_Heat_Sens_col         => plt_ew%Eco_Heat_Sens_col,           &
    VapXAir2Canopy_col        => plt_ew%VapXAir2Canopy_col,          &
    CanopyHeatStor_col        => plt_ew%CanopyHeatStor_col,          &
    QvET_col                  => plt_ew%QvET_col,                    &
    HeatFlx2Canopy_col        => plt_ew%HeatFlx2Canopy_col,          &
    LWRadCanG                 => plt_ew%LWRadCanG,                   &
    TairK                     => plt_ew%TairK,                       &
    HeatStorCanopy_pft        => plt_ew%HeatStorCanopy_pft,          &
    Eco_Heat_Latent_col       => plt_ew%Eco_Heat_Latent_col,         &
    CanH2OHeldVg              => plt_ew%CanH2OHeldVg,                &
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
  DO L=NU,MaxSoiL4Root_pft(NZ)
    trcs_plant_uptake_vr(idg_N2,L)=trcs_plant_uptake_vr(idg_N2,L)+RootN2Fix_pvr(L,NZ)
  ENDDO
!
!     TOTAL ENERGY, WATER, CO2 FLUXES
!
!     Eco_NetRad_col=total net SW+LW absorbed by canopy
!     RadNet2Canopy_pft=PFT net SW+LW absorbed by canopy
!     Eco_Heat_Latent_col=total canopy latent heat flux
!     EvapTransHeat_pft=PFT canopy latent heat flux
!     Eco_Heat_Sens_col=total canopy sensible heat flux
!     HeatXAir2PCan_pft=PFT canopy sensible heat flux
!     Eco_Heat_GrndSurf_col=total canopy storage heat flux
!     HeatStorCanopy_pft=PFT canopy storage heat flux
!     Canopy_NEE_col=total net CO2 fixation
!     CO2NetFix_pft=PFT net CO2 fixation
!     CanWat_col,CanH2OHeldVg=total water volume in canopy,on canopy surfaces
!     CanopyWater_pft,WatByPCanopy_pft=PFT water volume in canopy,on canopy surfaces
!     QvET_col,VapXAir2Canopy_col=total water flux to,from canopy,canopy surfaces
!     VapXAir2PCan,Transpiration_pft=water flux to,from canopy surfaces, inside canopy
!     CanopyHeatStor_col=total canopy water heat content
!     ENGYC=PFT canopy water heat content
!     CanopyLeafArea_col,StemArea_col=total leaf,stalk area
!     CanopyLeafArea_pft,CanopyStemArea_pft=PFT leaf,stalk area
!     ZCSNC,ZZSNC,ZPSNC=total net root-soil C,N,P exchange
!     HCUPTK,HZUPTK,HPUPTK=PFT net root-soil C,N,P exchange
!     TBALC,TBALN,TBALP=total C,N,P balance
!     BALC,BALN,BALP=PFT C,N,P balance
!     TRootGasLossDisturb_pft=total loss of root CO2, O2, CH4, N2O, NH3, H2
!     RootGasLossDisturb_pft=PFT loss of root CO2, O2, CH4, N2O, NH3, H2
!
  Eco_NetRad_col         = Eco_NetRad_col+RadNet2Canopy_pft(NZ)
  Eco_Heat_Latent_col    = Eco_Heat_Latent_col+EvapTransHeat_pft(NZ)
  Eco_Heat_Sens_col      = Eco_Heat_Sens_col+HeatXAir2PCan_pft(NZ)
  Eco_Heat_GrndSurf_col  = Eco_Heat_GrndSurf_col+HeatStorCanopy_pft(NZ)
  Canopy_NEE_col         = Canopy_NEE_col+CO2NetFix_pft(NZ)
  ETCanopy_CumYr_pft(NZ) = ETCanopy_CumYr_pft(NZ)+Transpiration_pft(NZ)+VapXAir2Canopy_pft(NZ)
  CanWat_col             = CanWat_col+CanopyWater_pft(NZ)
  CanH2OHeldVg           = CanH2OHeldVg+WatByPCanopy_pft(NZ)
  QvET_col               = QvET_col+Transpiration_pft(NZ)+VapXAir2Canopy_pft(NZ)
  VapXAir2Canopy_col     = VapXAir2Canopy_col+VapXAir2Canopy_pft(NZ)
  ENGYC                  = cpw*(WatByPCanopy_pft(NZ)+PrecIntcptByCanopy_pft(NZ)+VapXAir2Canopy_pft(NZ))*TKC(NZ)
  CanopyHeatStor_col     = CanopyHeatStor_col+ENGYC
  HeatFlx2Canopy_col     = HeatFlx2Canopy_col+ENGYC-ENGYX_pft(NZ)-(PrecIntcptByCanopy_pft(NZ)*cpw*TairK)
  ENGYX_pft(NZ)          = ENGYC
  LWRadCanG              = LWRadCanG+LWRadCanopy_pft(NZ)
  CanopyLeafArea_col     = CanopyLeafArea_col+CanopyLeafArea_pft(NZ)
  StemArea_col           = StemArea_col+CanopyStemArea_pft(NZ)

  DO NE=1,NumPlantChemElms
    LitrFallStrutElms_col(NE)     = LitrFallStrutElms_col(NE)-PlantRootSoilElmNetX_pft(NE,NZ)
    PlantElemntStoreLandscape(NE) = PlantElemntStoreLandscape(NE)+ElmBalanceCum_pft(NE,NZ)
  ENDDO

  DO NTG=idg_beg,idg_end-1
    TRootGasLossDisturb_pft(NTG)=TRootGasLossDisturb_pft(NTG)+RootGasLossDisturb_pft(NTG,NZ)
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
