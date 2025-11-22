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
  use EcoSIMCtrlMod, only: lverb, ldo_sp_mode
  use DebugToolMod,  only: PrintInfo
  use EcosimConst
  use PlantBGCPars
  use PlantAPIData
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: extracts
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  SUBROUTINE extracts(I,J)
!     execution begins here
  implicit none

  integer, intent(in) :: I, J
  character(len=*), parameter :: subname='extracts'
  integer :: NZ,L

  call PrintInfo('beg '//subname)
  call TotalLitrFall()

  call SumPlantRootGas(I,J)

  DO NZ=1,plt_site%NP

    call PrintRootTracer(I,J,NZ,'extract')
    IF(plt_pheno%IsPlantActive_pft(NZ).EQ.iActive)THEN
      if(.not.ldo_sp_mode)then
        call CalcTotalLeafArea(NZ)

        call TotalGasandSoluteUptake(I,J,NZ)
      endif
      call ExtractCanopyFluxes(I,J,NZ)

    ENDIF
  ENDDO

  call PrintInfo('end '//subname)  
  RETURN
  END subroutine extracts


!----------------------------------------------------------------------------------------------------
  subroutine TotalLitrFall()

  implicit none
  integer :: NZ,L,K,M
  integer :: NE
  associate(                                                          &
    LitrfallElms_pft      => plt_bgcr%LitrfallElms_pft       ,& !input  :plant element LitrFall, [g d-2 h-1]
    LitrfallElms_pvr      => plt_bgcr%LitrfallElms_pvr       ,& !input  :plant LitrFall element, [g d-2 h-1]
    MaxSoiL4Root_pft          => plt_morph%MaxSoiL4Root_pft          ,& !input  :maximum soil layer number for all root axes,[-]
    NP0                       => plt_site%NP0                        ,& !input  :intitial number of plant species,[-]
    StandDeadStrutElms_pft    => plt_biom%StandDeadStrutElms_pft     ,& !input  :standing dead element, [g d-2]
    LitrFallStrutElms_col     => plt_bgcr%LitrFallStrutElms_col      ,& !inoput :total LitrFall structural element mass, [g d-2 h-1]
    LitrfalStrutElms_vr       => plt_bgcr%LitrfalStrutElms_vr        ,& !inoput :total LitrFall element, [g d-2 h-1]
    StandingDeadStrutElms_col => plt_biom%StandingDeadStrutElms_col  ,& !inoput :total standing dead biomass chemical element, [g d-2]
    CanopyLeafAareZ_col       => plt_morph%CanopyLeafAareZ_col       ,& !output :total leaf area, [m2 d-2]
    CanopyLeafArea_col        => plt_morph%CanopyLeafArea_col        ,& !output :grid canopy leaf area, [m2 d-2]
    CanopyStemAareZ_col       => plt_morph%CanopyStemAareZ_col       ,& !output :total stem area, [m2 d-2]
    StemArea_col              => plt_morph%StemArea_col              ,& !output :grid canopy stem area, [m2 d-2]
    tCanLeafC_clyr            => plt_biom%tCanLeafC_clyr              & !output :total leaf carbon mass in canopy layers, [gC d-2]
  )
  DO NZ=1,NP0
!
!   TOTAL LitrFall OF ALL PLANT SPECIES
!
!   ZCSNC,ZZSNC,ZPSNC=total C,N,P LitrFall
!   HCSNC,HZSNC,HPSNC=hourly PFT C,N,P LitrFall from grosub.f
!   StandingDeadStrutElms_col=total standing dead C,N,P mass
!   WTSTG=PFT standing dead C,N,P mass
!   LitrfallElms_pvr,=cumulative PFT C,N,P LitrFall from grosub.f
!   LitrfalStrutElms_vr,=cumulative total C,N,P LitrFall
!
    DO NE=1,NumPlantChemElms
      LitrFallStrutElms_col(NE)=LitrFallStrutElms_col(NE)+LitrfallElms_pft(NE,NZ)
      StandingDeadStrutElms_col(NE)=StandingDeadStrutElms_col(NE)+StandDeadStrutElms_pft(NE,NZ)
    ENDDO

    DO  L=0,MaxSoiL4Root_pft(NZ)
      DO K=1,pltpar%NumOfPlantLitrCmplxs
        DO NE=1,NumPlantChemElms
          DO  M=1,pltpar%jsken
            LitrfalStrutElms_vr(NE,M,K,L)=LitrfalStrutElms_vr(NE,M,K,L)+AZMAX1(LitrfallElms_pvr(NE,M,K,L,NZ))
            if(LitrfalStrutElms_vr(NE,M,K,L)<0._r8)then
              write(*,*)'extract',K,L,LitrfalStrutElms_vr(NE,M,K,L),LitrfallElms_pvr(NE,M,K,L,NZ)
            endif
          enddo
        ENDDO
      ENDDO
    ENDDO
  ENDDO
     CanopyLeafArea_col = 0._r8
     StemArea_col       = 0._r8
  DO L                  = 1, NumCanopyLayers1
    CanopyLeafAareZ_col(L) = 0._r8
    tCanLeafC_clyr(L)        = 0._r8
    CanopyStemAareZ_col(L) = 0._r8
  ENDDO
  end associate
  end subroutine TotalLitrFall

!----------------------------------------------------------------------------------------------------
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
  associate(                                               &
    CanopyLeafAreaZ_pft => plt_morph%CanopyLeafAreaZ_pft  ,& !input  :canopy layer leaf area, [m2 d-2]
    CanopyLeafCLyr_pft  => plt_biom%CanopyLeafCLyr_pft    ,& !input  :canopy layer leaf C, [g d-2]
    CanopyStemAreaZ_pft => plt_morph%CanopyStemAreaZ_pft  ,& !input  :plant canopy layer stem area, [m2 d-2]
    CanopyLeafAareZ_col => plt_morph%CanopyLeafAareZ_col  ,& !inoput :total leaf area, [m2 d-2]
    CanopyStemAareZ_col => plt_morph%CanopyStemAareZ_col  ,& !inoput :total stem area, [m2 d-2]
    tCanLeafC_clyr      => plt_biom%tCanLeafC_clyr         & !inoput :total leaf carbon mass in canopy layers, [gC d-2]
  )
  DO L=1,NumCanopyLayers1
    CanopyLeafAareZ_col(L)=CanopyLeafAareZ_col(L)+CanopyLeafAreaZ_pft(L,NZ)
    tCanLeafC_clyr(L)=tCanLeafC_clyr(L)+CanopyLeafCLyr_pft(L,NZ)
    CanopyStemAareZ_col(L)=CanopyStemAareZ_col(L)+CanopyStemAreaZ_pft(L,NZ)
  ENDDO
  end associate
  end subroutine CalcTotalLeafArea

!----------------------------------------------------------------------------------------------------
  subroutine TotalGasandSoluteUptake(I,J,NZ)
!
!     TOTAL GAS AND SOLUTE UPTAKE BY ALL PLANT SPECIES
!

  implicit none
  integer, intent(in) :: I,J,NZ

  integer :: N,L,K,idg,NE,ids
 
  associate(                                                          &
    AREA3                     => plt_site%AREA3                      ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    RPlantRootH2OUptk_pvr     => plt_ew%RPlantRootH2OUptk_pvr        ,& !input  :root water uptake, [m2 d-2 h-1]
    Myco_pft                  => plt_morph%Myco_pft                  ,& !input  :mycorrhizal type (no or yes),[-]
    MaxNumRootLays            => plt_site%MaxNumRootLays             ,& !input  :maximum root layer number,[-]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    PlantPopulation_pft       => plt_site%PlantPopulation_pft        ,& !input  :plant population, [d-2]
    RCO2Emis2Root_rpvr         => plt_rbgc%RCO2Emis2Root_rpvr          ,& !input  :aqueous CO2 flux from roots to root water, [g d-2 h-1]
    RootH1PO4DmndBand_pvr     => plt_rbgc%RootH1PO4DmndBand_pvr      ,& !input  :HPO4 demand in band by each root population, [g d-2 h-1]
    RootH1PO4DmndSoil_pvr     => plt_rbgc%RootH1PO4DmndSoil_pvr      ,& !input  :HPO4 demand in non-band by each root population, [g d-2 h-1]
    RootH2PO4DmndBand_pvr     => plt_rbgc%RootH2PO4DmndBand_pvr      ,& !input  :root uptake of H2PO4 band, [g d-2 h-1]
    RootH2PO4DmndSoil_pvr     => plt_rbgc%RootH2PO4DmndSoil_pvr      ,& !input  :root uptake of H2PO4 non-band, [g d-2 h-1]
    RootLenDensPerPlant_pvr   => plt_morph%RootLenDensPerPlant_pvr   ,& !input  :root layer length density, [m m-3]
    RootMycoExudEUptk_pvr     => plt_rbgc%RootMycoExudEUptk_pvr      ,& !input  :root uptake (+ve) - exudation (-ve) of DOE, [g d-2 h-1]
    RootNH4DmndBand_pvr       => plt_rbgc%RootNH4DmndBand_pvr        ,& !input  :root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
    RootNH4DmndSoil_pvr       => plt_rbgc%RootNH4DmndSoil_pvr        ,& !input  :root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
    RootNO3DmndBand_pvr       => plt_rbgc%RootNO3DmndBand_pvr        ,& !input  :root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
    RootNO3DmndSoil_pvr       => plt_rbgc%RootNO3DmndSoil_pvr        ,& !input  :root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
    RootNutUptake_pvr         => plt_rbgc%RootNutUptake_pvr          ,& !input  :root uptake of Nutrient band, [g d-2 h-1]
    RootO2Dmnd4Resp_pvr       => plt_rbgc%RootO2Dmnd4Resp_pvr        ,& !input  :root O2 demand from respiration, [g d-2 h-1]
    RootO2Uptk_pvr            => plt_rbgc%RootO2Uptk_pvr             ,& !input  :aqueous O2 flux from roots to root water, [g d-2 h-1]
    RootO2_TotSink_pvr        => plt_bgcr%RootO2_TotSink_pvr         ,& !input  :root O2 sink for autotrophic respiraiton, [gC d-2 h-1]
    RootUptkSoiSol_pvr        => plt_rbgc%RootUptkSoiSol_pvr         ,& !input  :aqueous CO2 flux from roots to soil water, [g d-2 h-1]
    TKCanopy_pft              => plt_ew%TKCanopy_pft                 ,& !input  :canopy temperature, [K]
    TKS_vr                    => plt_ew%TKS_vr                       ,& !input  :mean annual soil temperature, [K]
    trcg_air2root_flx_pvr     => plt_rbgc%trcg_air2root_flx_pvr      ,& !input  :gaseous tracer flux through roots, [g d-2 h-1]
    RootCO2Ar2RootX_rpvr      => plt_rbgc%RootCO2Ar2RootX_rpvr       ,& !input  :root/myco respiration released to root/myco, [gC d-2 h-1]    
    RootCO2Ar2RootX_pvr      => plt_rbgc%RootCO2Ar2RootX_pvr         ,& !input  :root respiration released to root, [gC d-2 h-1]        
    REcoH1PO4DmndBand_vr      => plt_bgcr%REcoH1PO4DmndBand_vr       ,& !inoput :HPO4 demand in band by all microbial, root, myco populations, [gP d-2 h-1]
    REcoH1PO4DmndSoil_vr      => plt_bgcr%REcoH1PO4DmndSoil_vr       ,& !inoput :HPO4 demand in non-band by all microbial, root, myco populations, [gP d-2 h-1]
    REcoH2PO4DmndBand_vr      => plt_bgcr%REcoH2PO4DmndBand_vr       ,& !inoput :total root + microbial PO4 uptake band, [gP d-2 h-1]
    REcoH2PO4DmndSoil_vr      => plt_bgcr%REcoH2PO4DmndSoil_vr       ,& !inoput :total root + microbial PO4 uptake non-band, [gP d-2 h-1]
    REcoNH4DmndBand_vr        => plt_bgcr%REcoNH4DmndBand_vr         ,& !inoput :total root + microbial NH4 uptake band, [gN d-2 h-1]
    REcoNH4DmndSoil_vr        => plt_bgcr%REcoNH4DmndSoil_vr         ,& !inoput :total root + microbial NH4 uptake non-band, [gN d-2 h-1]
    REcoNO3DmndBand_vr        => plt_bgcr%REcoNO3DmndBand_vr         ,& !inoput :total root + microbial NO3 uptake band, [gN d-2 h-1]
    REcoNO3DmndSoil_vr        => plt_bgcr%REcoNO3DmndSoil_vr         ,& !inoput :total root + microbial NO3 uptake non-band, [gN d-2 h-1]
    REcoO2DmndResp_vr         => plt_bgcr%REcoO2DmndResp_vr          ,& !inoput :total root + microbial O2 uptake, [g d-2 h-1]
    RUptkRootO2_vr            => plt_bgcr%RUptkRootO2_vr             ,& !inoput :total root internal O2 flux, [g d-2 h-1]
    RootCO2Emis2Root_vr       => plt_bgcr%RootCO2Emis2Root_vr        ,& !inoput :total root CO2 flux, [gC d-2 h-1]
    RootCO2Emis2Root_pvr      => plt_bgcr%RootCO2Emis2Root_pvr       ,& !inoput :root CO2 flux, [gC d-2 h-1]    
    RootO2_TotSink_vr         => plt_bgcr%RootO2_TotSink_vr          ,& !inoput :all root O2 sink for autotrophic respiraiton, [gC d-2 h-1]
    THeatLossRoot2Soil_vr     => plt_ew%THeatLossRoot2Soil_vr        ,& !inoput :total root heat uptake, [MJ d-2 h-1]
    TWaterPlantRoot2Soil_vr   => plt_ew%TWaterPlantRoot2Soil_vr      ,& !inoput :total root water uptake, [m3 d-2 h-1]
    tRootMycoExud2Soil_vr     => plt_bgcr%tRootMycoExud2Soil_vr      ,& !inoput :total root element exchange, [g d-2 h-1]
    totRootLenDens_vr         => plt_morph%totRootLenDens_vr         ,& !inoput :total root length density, [m m-3]
    trcg_air2root_flx_vr      => plt_rbgc%trcg_air2root_flx_vr       ,& !inoput :total internal root gas flux, [gC d-2 h-1]    
    trcs_Soil2plant_uptake_pvr=> plt_rbgc%trcs_Soil2plant_uptake_pvr ,& !inoput :plant root-soil solute flux non-band, [g d-2 h-1]    
    trcs_Soil2plant_uptake_vr => plt_rbgc%trcs_Soil2plant_uptake_vr   & !inoput :total root-soil solute flux non-band, [g d-2 h-1]
  )
  
  DO L=NU,MaxNumRootLays
    DO N=1,Myco_pft(NZ)  
!
!     TOTAL ROOT DENSITY
!
!     totRootLenDens_vr=total root length density
!     RootLenDensPerPlant_pvr=PFT root length density per plant
!     RPlantRootH2OUptk_pvr=total water uptake
!     RPlantRootH2OUptk_pvr=PFT root water uptake
!     THeatLossRoot2Soil_vr=total convective heat in root water uptake, 
!     TKS=soil temperature
!     PP=PFT population, this is dynamic, and can goes to zero
!
      IF(N.EQ.ipltroot)THEN
        totRootLenDens_vr(L)=totRootLenDens_vr(L)+RootLenDensPerPlant_pvr(N,L,NZ)*PlantPopulation_pft(NZ)/AREA3(L)
      ENDIF
!
!     TOTAL WATER UPTAKE
!
      TWaterPlantRoot2Soil_vr(L) = TWaterPlantRoot2Soil_vr(L)+RPlantRootH2OUptk_pvr(N,L,NZ)

      
      if(RPlantRootH2OUptk_pvr(N,L,NZ)>0._r8)then
        !water lose from canopy to soil
        THeatLossRoot2Soil_vr(L)     = THeatLossRoot2Soil_vr(L)+RPlantRootH2OUptk_pvr(N,L,NZ)*cpw*TKCanopy_pft(NZ)        
      else
        !water lose from soil to canopy  
        THeatLossRoot2Soil_vr(L)     = THeatLossRoot2Soil_vr(L)+RPlantRootH2OUptk_pvr(N,L,NZ)*cpw*TKS_vr(L)
      endif
!
!     TOTAL ROOT BOUNDARY GAS FLUXES
!
      DO idg=idg_beg,idg_NH3
        trcg_air2root_flx_vr(idg,L)=trcg_air2root_flx_vr(idg,L)+trcg_air2root_flx_pvr(idg,N,L,NZ)
      ENDDO
      RootO2_TotSink_vr(L)      = RootO2_TotSink_vr(L) + RootO2_TotSink_pvr(N,L,NZ)  !total O2 uptake by roots to support root autotrophic respiraiton
      RootCO2Emis2Root_vr(L)    = RootCO2Emis2Root_vr(L)+RCO2Emis2Root_rpvr(N,L,NZ)
      RootCO2Emis2Root_pvr(L,NZ)=RootCO2Emis2Root_pvr(L,NZ)+RCO2Emis2Root_rpvr(N,L,NZ)
      RUptkRootO2_vr(L)      = RUptkRootO2_vr(L)+RootO2Uptk_pvr(N,L,NZ)

      !(>0) uptake from soil into roots
      DO idg=idg_beg,idg_end
        trcs_Soil2plant_uptake_vr(idg,L)     = trcs_Soil2plant_uptake_vr(idg,L)+RootUptkSoiSol_pvr(idg,N,L,NZ)
        trcs_Soil2plant_uptake_pvr(idg,L,NZ) = trcs_Soil2plant_uptake_pvr(idg,L,NZ)+RootUptkSoiSol_pvr(idg,N,L,NZ)
      ENDDO
      RootCO2Ar2RootX_pvr(L,NZ) = RootCO2Ar2RootX_pvr(L,NZ)+RootCO2Ar2RootX_rpvr(N,L,NZ)
 
      !Nutrients are sorted according to band|soil
      do ids=ids_NH4B,ids_nuts_end
        trcs_Soil2plant_uptake_vr(ids,L)=trcs_Soil2plant_uptake_vr(ids,L)+RootNutUptake_pvr(ids,N,L,NZ)
        trcs_Soil2plant_uptake_pvr(ids,L,NZ)=trcs_Soil2plant_uptake_pvr(ids,L,NZ)+RootNutUptake_pvr(ids,N,L,NZ)
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

!----------------------------------------------------------------------------------------------------
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

  associate(                                                          &
    CO2NetFix_pft             => plt_bgcr%CO2NetFix_pft              ,& !input  :canopy net CO2 exchange, [gC d-2 h-1]
    CanopyBiomWater_pft       => plt_ew%CanopyBiomWater_pft          ,& !input  :canopy water content, [m3 d-2]
    CanopyLeafArea_pft        => plt_morph%CanopyLeafArea_pft        ,& !input  :plant canopy leaf area, [m2 d-2]
    CanopyStemArea_pft        => plt_morph%CanopyStemArea_pft        ,& !input  :plant stem area, [m2 d-2]
    PlantElmBalCum_pft         => plt_site%PlantElmBalCum_pft          ,& !input  :cumulative plant element balance, [g d-2]
    EvapTransLHeat_pft        => plt_ew%EvapTransLHeat_pft           ,& !input  :canopy latent heat flux, [MJ d-2 h-1]
    HeatStorCanopy_pft        => plt_ew%HeatStorCanopy_pft           ,& !input  :canopy storage heat flux, [MJ d-2 h-1]
    HeatXAir2PCan_pft         => plt_ew%HeatXAir2PCan_pft            ,& !input  :canopy sensible heat flux, [MJ d-2 h-1]
    LWRadCanopy_pft           => plt_rad%LWRadCanopy_pft             ,& !input  :canopy longwave radiation, [MJ d-2 h-1]
    NH3Dep2Can_pft            => plt_bgcr%NH3Dep2Can_pft             ,& !input  :canopy NH3 flux, [gN d-2 h-1]
    PlantRootSoilElmNetX_pft  => plt_rbgc%PlantRootSoilElmNetX_pft   ,& !input  :net root element uptake (+ve) - exudation (-ve), [gC d-2 h-1]
    RadNet2Canopy_pft         => plt_rad%RadNet2Canopy_pft           ,& !input  :canopy net radiation, [MJ d-2 h-1]
    RootGasLossDisturb_pft    => plt_bgcr%RootGasLossDisturb_pft     ,& !input  :gaseous flux fron root disturbance, [g d-2 h-1]
    TKC_pft                   => plt_ew%TKC_pft                      ,& !input  :canopy temperature, [K]
    Transpiration_pft         => plt_ew%Transpiration_pft            ,& !input  :canopy transpiration, [m2 d-2 h-1]
    VHeatCapCanopy_pft        => plt_ew%VHeatCapCanopy_pft           ,& !input  :canopy heat capacity, [MJ d-2 K-1]
    VapXAir2Canopy_pft        => plt_ew%VapXAir2Canopy_pft           ,& !input  :canopy evaporation, [m2 d-2 h-1]
    WatHeldOnCanopy_pft       => plt_ew%WatHeldOnCanopy_pft          ,& !input  :canopy surface water content, [m3 d-2]
    CanopyHeatStor_col        => plt_ew%CanopyHeatStor_col           ,& !inoput :total canopy heat content, [MJ d-2]
    CanopyLeafArea_col        => plt_morph%CanopyLeafArea_col        ,& !inoput :grid canopy leaf area, [m2 d-2]
    CanopyWat_col             => plt_ew%CanopyWat_col                ,& !inoput :total canopy water content stored with dry matter, [m3 d-2]
    Canopy_NEE_col            => plt_bgcr%Canopy_NEE_col             ,& !inoput :total net CO2 fixation, [gC d-2]
    ENGYX_pft                 => plt_ew%ENGYX_pft                    ,& !inoput :canopy heat storage from previous time step, [MJ d-2]
    ETCanopy_CumYr_pft        => plt_ew%ETCanopy_CumYr_pft           ,& !inoput :total transpiration, [m H2O d-2]
    Eco_Heat_GrndSurf_col     => plt_ew%Eco_Heat_GrndSurf_col        ,& !inoput :ecosystem storage heat flux, [MJ d-2 h-1]
    Eco_Heat_Latent_col       => plt_ew%Eco_Heat_Latent_col          ,& !inoput :ecosystem latent heat flux, [MJ d-2 h-1]
    Eco_Heat_Sens_col         => plt_ew%Eco_Heat_Sens_col            ,& !inoput :ecosystem sensible heat flux, [MJ d-2 h-1]
    Eco_NetRad_col            => plt_rad%Eco_NetRad_col              ,& !inoput :ecosystem net radiation, [MJ d-2 h-1]
    HeatFlx2Canopy_col        => plt_ew%HeatFlx2Canopy_col           ,& !inoput :total canopy heat flux, [MJ d-2]
    LWRadCanG                 => plt_ew%LWRadCanG                    ,& !inoput :grid total canopy LW emission, [MJ d-2 h-1]
    LitrFallStrutElms_col     => plt_bgcr%LitrFallStrutElms_col      ,& !inoput :total LitrFall structural element mass, [g d-2 h-1]
    NH3Emis_CumYr_pft         => plt_bgcr%NH3Emis_CumYr_pft          ,& !inoput :total canopy NH3 flux, [gN d-2 ]
    PlantElemntStoreLandscape => plt_site%PlantElemntStoreLandscape  ,& !inoput :total plant element balance, [g d-2]
    QVegET_col                => plt_ew%QVegET_col                   ,& !inoput :total canopy evaporation + transpiration, [m3 d-2]
    StemArea_col              => plt_morph%StemArea_col              ,& !inoput :grid canopy stem area, [m2 d-2]
    TRootGasLossDisturb_col   => plt_rbgc%TRootGasLossDisturb_col    ,& !inoput :total root gas content, [g d-2]
    VapXAir2Canopy_col        => plt_ew%VapXAir2Canopy_col           ,& !inoput :grid canopy evaporation, [m3 d-2]
    WatHeldOnCanopy_col       => plt_ew%WatHeldOnCanopy_col           & !inoput :canopy surface water content, [m3 d-2]
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

  if(ldo_sp_mode)return
  Canopy_NEE_col         = Canopy_NEE_col+CO2NetFix_pft(NZ)
  CanopyLeafArea_col     = CanopyLeafArea_col+CanopyLeafArea_pft(NZ)
  StemArea_col           = StemArea_col+CanopyStemArea_pft(NZ)

  DO NE=1,NumPlantChemElms
    LitrFallStrutElms_col(NE)     = LitrFallStrutElms_col(NE)-PlantRootSoilElmNetX_pft(NE,NZ)
    PlantElemntStoreLandscape(NE) = PlantElemntStoreLandscape(NE)+PlantElmBalCum_pft(NE,NZ)
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
  ![tail]
  end module ExtractsMod
