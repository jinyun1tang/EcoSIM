module RootGasMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use StomatesMod   , only : stomates
  use minimathmod  , only : safe_adb,vapsat,AZMAX1,AZMIN1
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use PlantAPIData
  use TracerIDMod
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  public :: RootSoilGasExchange
  contains
!------------------------------------------------------------------------

  subroutine RootSoilGasExchange(I,J,N,L,NZ,FineRootRadius,FracPRoot4Uptake,FracSoiLayByPrimRoot,&
    RootAreaDivRadius_vr,dtPerPlantRootH2OUptake,FOXYX,PopPlantO2Uptake_vr)

  implicit none
  integer , intent(in) :: I,J,N,L,NZ
  real(r8), intent(in) :: FineRootRadius(jroots,JZ1),FracPRoot4Uptake(jroots,JZ1,JP1)
  real(r8), intent(in) :: FracSoiLayByPrimRoot(JZ1,JP1)
  real(r8), intent(in) :: RootAreaDivRadius_vr(jroots,JZ1)
  real(r8), intent(in) :: dtPerPlantRootH2OUptake,FOXYX
  real(r8), intent(out):: PopPlantO2Uptake_vr
  !local variables
  integer :: M,MX
  real(r8) :: B,C
  real(r8) :: trcg_gmas_loc(idg_beg:idg_end-1)
  real(r8) :: trcaqu_conc_soi_loc(idg_beg:idg_end)
  real(r8) :: trcrootconc_loc(idg_beg:idg_end-1)
  real(r8) :: GasDifc_loc(idg_beg:idg_end-1)
  real(r8) :: SolDifc_loc(idg_beg:idg_end-1)
  real(r8) :: trcg_gcon_loc(idg_beg:idg_end-1)
  real(r8) :: COXYR,GasPX(idg_beg:idg_end-1)
  real(r8) :: DIFOP,DIFLGas(idg_beg:idg_end),DIFOX,DiffusivitySolutEffP
  real(r8) :: DFAGas(idg_beg:idg_end-1),DFGP,RUPOSX
  real(r8) :: trcs_rootml_loc(idg_beg:idg_end-1),trc_solml_loc(idg_beg:idg_end)
  real(r8) :: trc_gasml_loc(idg_beg:idg_end-1)
  real(r8) :: O2AquaDiffusvityP,RTVLWA
  real(r8) :: RTVLWB,ROXYFX,RCO2FX,ROXYLX
  real(r8) :: RDFQSolute(idg_beg:idg_end)
  real(r8) :: RTCR1,RTCR2,RTCRA,RTARRX,RCO2PX,RRADS
  real(r8) :: RMFGas(idg_beg:idg_end),RootOxyUptakePerPlant,RDFOXS
  real(r8) :: RDFOXP,RUPSolute(idg_beg:idg_end)
  real(r8) :: RDXSolute(idg_beg:idg_end),RDFSolute(idg_beg:idg_end)
  real(r8) :: RUPOST,RUPNTX
  real(r8) :: trcg_Root2Soil_flx(idg_beg:idg_end-1)
  real(r8) :: trcg_air2root_flx1(idg_beg:idg_end-1)
  real(r8) :: THETW1,THETM
  real(r8) :: RootOxyDemandPerPlant
  real(r8) :: DisolvedGasVolume(idg_beg:idg_end-1),VLWatMicPMO,VLWatMicPMM,VLsoiAirPMM
  real(r8) :: VOLWSP,VLWatMicPMA,VLWatMicPMB,VOLWSA,VOLWSB,VOLWSolute(idg_beg:idg_end)
  real(r8) :: VOLPNH,VOLPNB
  real(r8) :: X
  real(r8) :: ZH3PA,ZH3PB,ZH3GA,ZH3GB
  integer  :: NTG
!     begin_execution
  associate(                                                       &
    RootStrutElms_pft        => plt_biom%RootStrutElms_pft,        &
    ZERO4Groth_pft           => plt_biom%ZERO4Groth_pft,           &
    PlantPopulation_pft      => plt_site%PlantPopulation_pft,      &
    DPTHZ_vr                 => plt_site%DPTHZ_vr,                 &
    AtmGasc                  => plt_site%AtmGasc,                  &
    ZEROS                    => plt_site%ZEROS,                    &
    ZERO                     => plt_site%ZERO,                     &
    VLWatMicPM_vr            => plt_site%VLWatMicPM_vr,            &
    VLsoiAirPM               => plt_site%VLsoiAirPM,               &
    TortMicPM_vr             => plt_site%TortMicPM_vr,             &
    FILM                     => plt_site%FILM,                     &
    RO2GasXchangePrev_vr     => plt_bgcr%RO2GasXchangePrev_vr,     &
    RCO2GasFlxPrev_vr        => plt_bgcr%RCO2GasFlxPrev_vr,        &
    RO2AquaXchangePrev_vr    => plt_bgcr%RO2AquaXchangePrev_vr,    &
    RootO2Uptk_pvr           => plt_rbgc%RootO2Uptk_pvr,           &
    RAutoRootO2Limter_pvr    => plt_rbgc%RAutoRootO2Limter_pvr,    &
    ZERO4Uptk_pft            => plt_rbgc%ZERO4Uptk_pft,            &
    RootRespPotent_pvr       => plt_rbgc%RootRespPotent_pvr,       &
    RootO2Dmnd4Resp_pvr      => plt_rbgc%RootO2Dmnd4Resp_pvr,      &
    RO2UptkSoilM_vr          => plt_rbgc%RO2UptkSoilM_vr,          &
    RootCO2Emis_pvr          => plt_rbgc%RootCO2Emis_pvr,          &
    trcg_air2root_flx__pvr   => plt_rbgc%trcg_air2root_flx__pvr,   &
    trcg_Root_DisEvap_flx_vr => plt_rbgc%trcg_Root_DisEvap_flx_vr, &
    RUPGasSol_vr             => plt_rbgc%RUPGasSol_vr,             &
    RootCO2Autor_pvr         => plt_rbgc%RootCO2Autor_pvr,         &
    trcg_rootml_pvr          => plt_rbgc%trcg_rootml_pvr,          &
    trcs_rootml_pvr          => plt_rbgc%trcs_rootml_pvr,          &
    TScal4Difsvity_vr        => plt_soilchem%TScal4Difsvity_vr,    &
    trcs_VLN_vr              => plt_soilchem%trcs_VLN_vr,          &
    trc_solml_vr             => plt_soilchem%trc_solml_vr,         &
    trc_gascl_vr             => plt_soilchem%trc_gascl_vr,         &
    GasDifc_vr               => plt_soilchem%GasDifc_vr,           &
    SoluteDifusvty_vr        => plt_soilchem%SoluteDifusvty_vr,    &
    GasSolbility_vr          => plt_soilchem%GasSolbility_vr,      &
    trc_solcl_vr             => plt_soilchem%trc_solcl_vr,         &
    THETY_vr                 => plt_soilchem%THETY_vr,             &
    VLSoilMicP_vr            => plt_soilchem%VLSoilMicP_vr,        &
    THETPM                   => plt_soilchem%THETPM,               &
    trc_gasml_vr             => plt_soilchem%trc_gasml_vr,         &
    DiffusivitySolutEff      => plt_soilchem%DiffusivitySolutEff,  &
    iPlantCalendar_brch      => plt_pheno%iPlantCalendar_brch,     &
    RootPoreTortu4Gas        => plt_morph%RootPoreTortu4Gas,       &
    Root1stRadius_pvr        => plt_morph%Root1stRadius_pvr,       &
    Root2ndAveLen_pvr        => plt_morph%Root2ndAveLen_pvr,       &
    RootPoreVol_pvr          => plt_morph%RootPoreVol_pvr,         &
    RootLenPerPlant_pvr      => plt_morph%RootLenPerPlant_pvr,     &
    Root2ndXNum_pvr          => plt_morph%Root2ndXNum_pvr,         &
    Root2ndRadius_pvr        => plt_morph%Root2ndRadius_pvr,       &
    RootRaidus_rpft          => plt_morph%RootRaidus_rpft,         &
    RootVH2O_pvr             => plt_morph%RootVH2O_pvr,            &
    RootPorosity_pft         => plt_morph%RootPorosity_pft,        &
    Root1stXNumL_pvr         => plt_morph%Root1stXNumL_pvr,        &
    NGTopRootLayer_pft       => plt_morph%NGTopRootLayer_pft,      &
    MainBranchNum_pft        => plt_morph%MainBranchNum_pft        &
  )
  
  IF(RootRespPotent_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ).AND.RootVH2O_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ) &
    .AND.FOXYX.GT.ZERO4Uptk_pft(NZ))THEN
!
!     INITIALIZE VARIABLES USED IN ROOT GAS EXCHANGE
!     (CO2, O2, CH4, N2, N2O, NH3, H2)
!
!     trcg_gmas_loc(idg_CO2),trcs_rootml_loc(idg_CO2),trc_gasml_loc(idg_CO2),trc_solml_loc(idg_CO2)=gaseous,aqueous CO2 in root,soil
!     trcaqu_conc_soi_loc(idg_CH4),trcrootconc_loc(idg_CH4)=aqueous gas concentration in soil,root
!     RTVLWA,RTVLWB=root aqueous volume in non-band,band
!     dts_gas=time step of flux calculation
!     RootOxyDemandPerPlant=O2 demand per plant at time step of flux calculation
!     ROXYFX=net O2 gas flux at time step of flux calculation
!     RCO2FX=net CO2 gas flux at time step of flux calculation
!     ROXYLX=net O2 aqueous flux at time step of flux calculation
!

    DO NTG=idg_beg,idg_end-1
      trcg_gmas_loc(NTG)=AMAX1(ZERO4Groth_pft(NZ),trcg_rootml_pvr(NTG,N,L,NZ))
      trcs_rootml_loc(NTG)=AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_pvr(NTG,N,L,NZ))
    ENDDO

    DO NTG=idg_beg,idg_end
      if(NTG/=idg_O2)trc_solml_loc(NTG)=AMAX1(ZERO4Groth_pft(NZ),trc_solml_vr(NTG,L)*FracPRoot4Uptake(N,L,NZ))
    enddo
    trc_solml_loc(idg_O2)=trc_solml_vr(idg_O2,L)*FOXYX

!  the two lines below may be redundant
!    trc_gasml_loc(idg_CO2)=AMAX1(ZERO4Groth_pft(NZ),trc_gasml_vr(idg_CO2,L)*FracPRoot4Uptake(N,L,NZ))
!    trc_gasml_loc(idg_O2)=AMAX1(ZERO4Groth_pft(NZ),trc_gasml_vr(idg_O2,L)*FOXYX)

    RTVLWA=RootVH2O_pvr(N,L,NZ)*trcs_VLN_vr(ids_NH4,L)
    RTVLWB=RootVH2O_pvr(N,L,NZ)*trcs_VLN_vr(ids_NH4B,L)
    RootOxyDemandPerPlant=RootO2Dmnd4Resp_pvr(N,L,NZ)*dts_gas/PlantPopulation_pft(NZ)
    ROXYFX=RO2GasXchangePrev_vr(L)*FOXYX*dts_gas
    RCO2FX=RCO2GasFlxPrev_vr(L)*FOXYX*dts_gas
    ROXYLX=RO2AquaXchangePrev_vr(L)*FOXYX*dts_gas

!
!     GASEOUS AND AQUEOUS DIFFUSIVITIES IN ROOT AND SOIL
!
!     *SGL1=diffusivity
!     RootPoreTortu4Gas=tortuosity effect of root porosity on diffusivity
!     CG=CO2g,OG=O2g,CH=CH4g,Z2=N2Og,ZH=NH3g,HG=H2g
!     CL=CO2s,OL=O2s,CQ=CH4s,ZV=N2Os,ZN=NH3s,HL=H2s
!
    do NTG=idg_beg,idg_end-1
      GasDifc_loc(NTG)=GasDifc_vr(NTG,L)*dts_gas*RootPoreTortu4Gas(N,NZ)
      SolDifc_loc(NTG)=SoluteDifusvty_vr(NTG,L)*dts_gas*FOXYX
    enddo
    O2AquaDiffusvityP=SoluteDifusvty_vr(idg_O2,L)*dts_gas

    RDFQSolute(idg_beg:idg_end)=0.0_r8
    trcg_Root2Soil_flx(idg_beg:idg_end-1)=0.0_r8
!
!     ROOT CONDUCTANCE TO GAS TRANSFER
!
!     WTRTS=total root,myco mass
!     FracSoiLayByPrimRoot=fraction of each soil layer with primary root
!     RTCR1,RTCR2,RTCRA=cross-sectional area/length of
!     primary,secondary,total root,myco system
!     Root1stXNumL_pvr,Root2ndXNum_pvr=number of root,myco primary,secondary axes
!     Root1stRadius_pvr,Root2ndRadius_pvr=primary,secondary root radius
!     DPTHZ=depth of primary root from surface
!     Root2ndAveLen_pvr=average secondary root length
!
    IF(RootStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ).AND.FracSoiLayByPrimRoot(L,NZ).GT.ZERO)THEN
      RTCR1=AMAX1(PlantPopulation_pft(NZ),Root1stXNumL_pvr(N,L,NZ)) &
        *PICON*Root1stRadius_pvr(N,L,NZ)**2/DPTHZ_vr(L)
      RTCR2=(Root2ndXNum_pvr(N,L,NZ)*PICON*Root2ndRadius_pvr(N,L,NZ)**2 &
        /Root2ndAveLen_pvr(N,L,NZ))/FracSoiLayByPrimRoot(L,NZ)
      IF(RTCR2.GT.RTCR1)THEN
        RTCRA=RTCR1*RTCR2/(RTCR1+RTCR2)
      ELSE
        RTCRA=RTCR1
      ENDIF
    ELSE
      RTCRA=0.0_r8
    ENDIF
!
!     VARIABLES USED TO CALCULATE ROOT GAS TRANSFER
!     BETWEEN AQUEOUS AND GASEOUS PHASES
!
!     RootLenPerPlant_pvr=root,myco length per plant
!     iPlantCalendar_brch(ipltcal_Emerge,=emergence date
!     RTARR=root surface area/radius for uptake
!     RootRaidus_rpft=path length for radial diffusion within root
!     DIFOP=aqueous diffusivity of O2 within root
!     RootVH2O_pvr=root,myco aqueous volume
!     S*L=solubility of gas in water from hour1.f:
!     CO2=CO2,OXY=O2,CH4=CH4,N2O=N2O,NH3=NH3,H2G=H2
!     DF*A=root-atmosphere gas conductance
!     DFGP=rate const for equilibrn of gas concn in gaseous-aqueous phases
!     RCO2PX=root CO2 gas flux at time step for gas flux calculations
!     RootCO2Autor_pvr=root CO2 flux from grosub.f
!

    IF(N.EQ.ipltroot .AND. iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).GT.0 &
      .AND. RootLenPerPlant_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
      RTARRX=RootAreaDivRadius_vr(N,L)/RootRaidus_rpft(N,NZ)
      DIFOP=O2AquaDiffusvityP*RTARRX
      DO NTG=idg_beg,idg_end-1
        DisolvedGasVolume(NTG)=RootVH2O_pvr(N,L,NZ)*GasSolbility_vr(NTG,L)
        DFAGas(NTG)=GasDifc_loc(NTG)*RTCRA
      ENDDO
    ELSE
      RTARRX=0.0_r8
      DIFOP=0.0_r8
      DisolvedGasVolume(idg_beg:idg_end-1)=0.0_r8
      DFAGas(idg_beg:idg_end-1)=0.0_r8
    ENDIF

    DFGP=AMIN1(1.0,XNPD*SQRT(RootPorosity_pft(N,NZ))*TScal4Difsvity_vr(L))
    RCO2PX=-RootCO2Autor_pvr(N,L,NZ)*dts_gas
!
!     SOLVE FOR GAS EXCHANGE IN SOIL AND ROOTS DURING ROOT UPTAKE
!     AT SMALLER TIME STEP NPH
!
    D99: DO M=1,NPH
!
!     AQUEOUS GAS DIFFUSIVITY THROUGH SOIL WATER TO ROOT
!
!     gas code:CO2=CO2,OXY=O2,CH4=CH4,Z2O=N2O,NH3=NH3 non-band,
!     NHB=NH3 band,H2G=H2
!     VLWatMicPMM,VLsoiAirPMM=soil micropore water,air volume
!     FOXYX=root fraction of total O2 demand from previous hour
!     FracPRoot4Uptake=PFT fraction of biome root mass
!     VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
!     VOLX=soil volume excluding rock,macropores
!     THETW1=soil water concentration
!     TortMicPM_vr=soil tortuosity
!     FILM=soil water film thickness
!     FineRootRadius=root radius
!     RRADS=path length for radial diffusion from soil to root
!     DIF*=aqueous diffusivity from soil to root:OL=O2,CL=CH4
!     ZL=N2O,NL=NH3 non-band,NB=NH4 band,HL=H2
!     C*G=soil gaseous concentration
!     VOLW*,VOLP*=VLWatMicPMM,VLsoiAirPMM*gas solubility
!
      VLWatMicPMO=VLWatMicPM_vr(M,L)*FOXYX
      VLWatMicPMM=VLWatMicPM_vr(M,L)*FracPRoot4Uptake(N,L,NZ)
      VLsoiAirPMM=VLsoiAirPM(M,L)*FracPRoot4Uptake(N,L,NZ)
      VOLWSP=RootVH2O_pvr(N,L,NZ)+VLWatMicPMM
      VLWatMicPMA=VLWatMicPMM*trcs_VLN_vr(ids_NH4,L)
      VLWatMicPMB=VLWatMicPMM*trcs_VLN_vr(ids_NH4B,L)
      VOLWSA=RTVLWA+VLWatMicPMA
      VOLWSB=RTVLWB+VLWatMicPMB
      THETW1=AZMAX1(VLWatMicPM_vr(M,L)/VLSoilMicP_vr(L))

      IF(THETW1.GT.THETY_vr(L) .AND. FracPRoot4Uptake(N,L,NZ).GT.ZERO4Uptk_pft(NZ))THEN
        THETM=TortMicPM_vr(M,L)*THETW1
        RRADS=LOG((FILM(M,L)+FineRootRadius(N,L))/FineRootRadius(N,L))
        RTARRX=RootAreaDivRadius_vr(N,L)/RRADS
        do NTG=idg_beg,idg_end-1
          DIFLGas(NTG)=THETM*SolDifc_loc(NTG)*RTARRX
          if(NTG/=idg_O2 .and. NTG/=idg_CO2)then
            trc_gasml_loc(NTG)=trc_gascl_vr(NTG,L)*VLsoiAirPMM
          endif          
        enddo
        DIFLGas(idg_NH3)=DIFLGas(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
        DIFLGas(idg_NH3B)=DIFLGas(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)

        DO NTG=idg_beg,idg_end
          VOLWSolute(NTG)=VLWatMicPMM*GasSolbility_vr(NTG,L)
        ENDDO
        VOLWSolute(idg_NH3)=VOLWSolute(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
        VOLWSolute(idg_NH3B)=VOLWSolute(idg_NH3B)*trcs_VLN_vr(ids_NH4B,L)

        VOLPNH=VLsoiAirPMM*trcs_VLN_vr(ids_NH4,L)
        VOLPNB=VLsoiAirPMM*trcs_VLN_vr(ids_NH4B,L)
!
!     MASS FLOW OF GAS FROM SOIL TO ROOT AT SHORTER TIME STEP NPT
!
!     C*S1=soil aqueous concentration non-band
!     C*B1=soil aqueous concentration band
!     C*A1=root gaseous concentration
!     C*P1=root aqueous concentration
!     ROXYLX=soil net O2 aqueous flux
!     VLWatMicPMM=micropore water volume
!     RootVH2O_pvr,RootPoreVol_pvr=root aqueous,gaseous volume
!     RMF*=soil convective solute flux:COS=CO2,OXS=O2,CHS=CH4,
!     N2S=N2O,NHS=NH3 non-band,NHB=NH3 band,HGS=H2
!     dtPerPlantRootH2OUptake=water uptake
!
        D90: DO MX=1,NPT
          trc_solml_loc(idg_O2)=trc_solml_loc(idg_O2)+ROXYLX
          do NTG=idg_beg,idg_end
            if(NTG/=idg_O2)trcaqu_conc_soi_loc(NTG)=AZMAX1(trc_solml_loc(NTG)/VLWatMicPMM)
          enddo
          trcaqu_conc_soi_loc(idg_O2)=AMIN1(AtmGasc(idg_O2)*GasSolbility_vr(idg_O2,L),AZMAX1(trc_solml_loc(idg_O2)/VLWatMicPMO))
          IF(RootPoreVol_pvr(N,L,NZ).GT.ZERO)THEN
            DO NTG=idg_beg,idg_end-1
              trcg_gcon_loc(NTG)=AZMAX1(trcg_gmas_loc(NTG)/RootPoreVol_pvr(N,L,NZ))
            ENDDO
          ELSE
            trcg_gcon_loc(idg_beg:idg_end-1)=0.0_r8
          ENDIF

          do NTG=idg_beg,idg_end-1
            trcrootconc_loc(NTG)=AZMAX1(trcs_rootml_loc(NTG)/RootVH2O_pvr(N,L,NZ))
          enddo
          trcrootconc_loc(idg_O2)=AMIN1(AtmGasc(idg_O2)*GasSolbility_vr(idg_O2,L),trcrootconc_loc(idg_O2))

          DIFOX=DIFLGas(idg_O2)+DIFOP
          do NTG=idg_beg,idg_end
            RMFGas(NTG)=dtPerPlantRootH2OUptake*trcaqu_conc_soi_loc(NTG)
          enddo
          RMFGas(idg_NH3)=RMFGas(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
          RMFGas(idg_NH3B)=RMFGas(idg_NH3B)*trcs_VLN_vr(ids_NH4B,L)
!
!     SOLUTION FOR MASS FLOW + DIFFUSION OF O2 IN AQUEOUS PHASES OF
!     SOIL AND ROOT = ACTIVE UPTAKE OF O2 BY ROOT
!
!     DIFLGas(idg_O2)=O2 aqueous diffusivity from soil to root
!     dtPerPlantRootH2OUptake=water uptake
!     DIFOP=aqueous diffusivity of O2 within root
!     trcaqu_conc_soi_loc(idg_O2),trcrootconc_loc(idg_O2)=soil,root aqueous O2 concentration
!     RootOxyDemandPerPlant=O2 demand per plant
!     RootOxyUptakePerPlant=root O2 uptake per plant
!     COXYR=aqueous O2 concentration at root surface
!     RDFOXS,RDFOXP=aqueous O2 diffusion per plant:soil-root,within root
!
          X=(DIFLGas(idg_O2)+dtPerPlantRootH2OUptake)*trcaqu_conc_soi_loc(idg_O2)+&
            DIFOP*trcrootconc_loc(idg_O2)

          IF(X.GT.ZERO.AND.trc_solml_loc(idg_O2).GT.ZERO4Groth_pft(NZ))THEN
            !root take up O2
            B=-RootOxyDemandPerPlant-DIFOX*OXKM-X
            C=X*RootOxyDemandPerPlant            
            RootOxyUptakePerPlant=(-B-SQRT(B*B-4.0_r8*C))/2.0_r8

            COXYR=(X-RootOxyUptakePerPlant)/DIFOX
            RDFOXS=RMFGas(idg_O2)+DIFLGas(idg_O2)*(trcaqu_conc_soi_loc(idg_O2)-COXYR)
            RDFOXP=DIFOP*(trcrootconc_loc(idg_O2)-COXYR)
          ELSE
            !root release O2
            X=DIFOP*trcrootconc_loc(idg_O2)
            IF(X.GT.ZERO.AND.trcs_rootml_loc(idg_O2).GT.ZERO4Groth_pft(NZ))THEN
              B=-RootOxyDemandPerPlant-DIFOP*OXKM-X
              C=X*RootOxyDemandPerPlant
              RootOxyUptakePerPlant=(-B-SQRT(B*B-4.0_r8*C))/2.0_r8

              COXYR=(X-RootOxyUptakePerPlant)/DIFOP
              RDFOXS=0.0_r8
              RDFOXP=DIFOP*(trcrootconc_loc(idg_O2)-COXYR)
            ELSE
              RootOxyUptakePerPlant=0.0_r8
              COXYR=0.0_r8
              RDFOXS=0.0_r8
              RDFOXP=0.0_r8
            ENDIF
          ENDIF
!
!     MASS FLOW + DIFFUSIVE EXCHANGE OF OTHER GASES
!     BETWEEN ROOT AND SOIL, CONSTRAINED BY COMPETITION
!     WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!     RUPOSX,-RUPSolute(idg_O2)=aqueous O2 uptake from soil,root
!     PP=PFT population
!     RDFSolute(:),RUPSolute(:)=aqueous gas soil-root diffusion,root uptake, > 0 into soil
!     RMF*=soil convective solute flux
!     DIF*=aqueous diffusivity from soil to root
!     C*S1=soil aqueous concentration non-band
!     C*B1=soil aqueous concentration band
!     C*P1=root aqueous concentration
!
          RUPOSX=RDFOXS*PlantPopulation_pft(NZ)

          RUPSolute(idg_O2)=-RDFOXP*PlantPopulation_pft(NZ)

          RDFSolute(idg_CO2)=RMFGas(idg_CO2)+DIFLGas(idg_CO2)*&
            (trcaqu_conc_soi_loc(idg_CO2)-trcrootconc_loc(idg_CO2))
          RDXSolute(idg_CO2)=(RootVH2O_pvr(N,L,NZ)*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_CO2)) &
            -VLWatMicPMM*AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_loc(idg_CO2)))/VOLWSP

          IF(RDFSolute(idg_CO2).GT.0.0_r8)THEN
            RUPSolute(idg_CO2)=AMIN1(AZMAX1(RDXSolute(idg_CO2)),RDFSolute(idg_CO2) &
              *PlantPopulation_pft(NZ))
          ELSE
            RUPSolute(idg_CO2)=AMAX1(AZMIN1(RDXSolute(idg_CO2)),RDFSolute(idg_CO2) &
              *PlantPopulation_pft(NZ))
          ENDIF

          IF(N.EQ.ipltroot)THEN
            DO NTG=idg_beg,idg_end-1
              if(NTG/=idg_CO2.and.NTG/=idg_NH3.and.NTG/=idg_O2)then
                RDFSolute(NTG)=RMFGas(NTG)+DIFLGas(NTG)*(trcaqu_conc_soi_loc(NTG)-trcrootconc_loc(NTG))
                RDXSolute(NTG)=(RootVH2O_pvr(N,L,NZ)*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(NTG)) &
                  -VLWatMicPMM*AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_loc(NTG)))/VOLWSP
                IF(RDFSolute(NTG).GT.0.0_r8)THEN
                  RUPSolute(NTG)=AMIN1(AZMAX1(RDXSolute(NTG)),RDFSolute(NTG)*PlantPopulation_pft(NZ))
                ELSE
                  RUPSolute(NTG)=AMAX1(AZMIN1(RDXSolute(NTG)),RDFSolute(NTG)*PlantPopulation_pft(NZ))
                ENDIF
              ENDIF
            ENDDO

            RDFSolute(idg_NH3)=RMFGas(idg_NH3)+DIFLGas(idg_NH3)*&
              (trcaqu_conc_soi_loc(idg_NH3)-trcrootconc_loc(idg_NH3))
            IF(VOLWSA.GT.ZERO4Groth_pft(NZ))THEN
              ZH3PA=trcs_rootml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
              RDXSolute(idg_NH3)=(RTVLWA*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_NH3)) &
                -VLWatMicPMA*AMAX1(ZERO4Groth_pft(NZ),ZH3PA))/VOLWSA
            ELSE
              RDXSolute(idg_NH3)=0.0_r8
            ENDIF
            IF(RDFSolute(idg_NH3).GT.0.0_r8)THEN
              RUPSolute(idg_NH3)=AMIN1(AZMAX1(RDXSolute(idg_NH3)),RDFSolute(idg_NH3) &
                *PlantPopulation_pft(NZ))
            ELSE
              RUPSolute(idg_NH3)=AMAX1(AZMIN1(RDXSolute(idg_NH3)),RDFSolute(idg_NH3) &
                *PlantPopulation_pft(NZ))
            ENDIF

            RDFSolute(idg_NH3B)=RMFGas(idg_NH3B)+DIFLGas(idg_NH3B)* &
              (trcaqu_conc_soi_loc(idg_NH3B)-trcrootconc_loc(idg_NH3))
            IF(VOLWSB.GT.ZERO4Groth_pft(NZ))THEN
              ZH3PB=trcs_rootml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)
              RDXSolute(idg_NH3B)=(RTVLWB*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_NH3B)) &
                -VLWatMicPMB*AMAX1(ZERO4Groth_pft(NZ),ZH3PB))/VOLWSB
            ELSE
              RDXSolute(idg_NH3B)=0.0_r8
            ENDIF
            IF(RDFSolute(idg_NH3B).GT.0.0_r8)THEN
              RUPSolute(idg_NH3B)=AMIN1(AZMAX1(RDXSolute(idg_NH3B)),RDFSolute(idg_NH3B) &
                *PlantPopulation_pft(NZ))
            ELSE
              RUPSolute(idg_NH3B)=AMAX1(AZMIN1(RDXSolute(idg_NH3B)),RDFSolute(idg_NH3B) &
                *PlantPopulation_pft(NZ))
            ENDIF

          ELSE
            RUPSolute(idg_beg:idg_end)=0.0_r8
          ENDIF
!
!     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN SOIL
!     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
!     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENTS
!     FROM 'WATSUB'
!
!     THETPM,THETX=air-filled porosity,minimum THETPM
!     R*DFQ=soil gas exchange between gaseous-aqueous phases
!     DiffusivitySolutEff=rate constant for soil gas exchange from watsub.f
!     trc_gasml_loc(:),trc_solml_loc(:)=gaseous,aqueous gas in soil
!     RUPSolute=root aqueous gas uptake
!     ROXYLX=soil net O2 aqueous flux

!     VLWatMicPMM,VLsoiAirPMM=soil micropore water,air volume
!     VOLW*=VLWatMicPMM*gas solubility
!         
          IF(THETPM(M,L).GT.THETX)THEN
            DiffusivitySolutEffP=FracPRoot4Uptake(N,L,NZ)*DiffusivitySolutEff(M,L)
            RDFQSolute(idg_CO2)=DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(idg_CO2))&
              *VOLWSolute(idg_CO2)-(AMAX1(ZEROS,trc_solml_loc(idg_CO2))-RUPSolute(idg_CO2))*VLsoiAirPMM) &
              /(VOLWSolute(idg_CO2)+VLsoiAirPMM)
            RUPOST=RUPOSX-ROXYLX
            RDFQSolute(idg_O2)=DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(idg_O2))*VOLWSolute(idg_O2) &
              -(AMAX1(ZEROS,trc_solml_loc(idg_O2))-RUPOST)*VLsoiAirPMM)/(VOLWSolute(idg_O2)+VLsoiAirPMM)

            IF(N.EQ.ipltroot)THEN
              DO NTG=idg_beg,idg_NH3-1
                if(NTG/=idg_CO2.and.NTG/=idg_O2)then
                  RDFQSolute(NTG)=DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(NTG))*VOLWSolute(NTG) &
                    -(AMAX1(ZEROS,trc_solml_loc(NTG))-RUPSolute(NTG))*VLsoiAirPMM)/(VOLWSolute(NTG)+VLsoiAirPMM)
                endif
              ENDDO

              IF(VOLWSolute(idg_NH3)+VOLPNH.GT.ZERO4Groth_pft(NZ))THEN
                ZH3GA=trc_gasml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
                RDFQSolute(idg_NH3)=AMIN1(RUPSolute(idg_NH3),AMAX1(-RUPSolute(idg_NH3) &
                  ,DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),ZH3GA)*VOLWSolute(idg_NH3) &
                  -(AMAX1(ZEROS,trc_solml_loc(idg_NH3))-RUPSolute(idg_NH3))*VOLPNH) &
                  /(VOLWSolute(idg_NH3)+VOLPNH)))
              ELSE
                RDFQSolute(idg_NH3)=0.0_r8
              ENDIF

              IF(VOLWSolute(idg_NH3B)+VOLPNB.GT.ZERO4Groth_pft(NZ))THEN
                ZH3GB=trc_gasml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)
                RDFQSolute(idg_NH3B)=AMIN1(RUPSolute(idg_NH3),AMAX1(-RUPSolute(idg_NH3) &
                  ,DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),ZH3GB)*VOLWSolute(idg_NH3B)   &
                  -(AMAX1(ZEROS,trc_solml_loc(idg_NH3B))-RUPSolute(idg_NH3B))*VOLPNB)  &
                  /(VOLWSolute(idg_NH3B)+VOLPNB)))
              ELSE
                RDFQSolute(idg_NH3B)=0.0_r8
              ENDIF
            ELSE
              RDFQSolute(idg_beg:idg_end)=0.0_r8
            ENDIF
          ELSE
            RDFQSolute(idg_beg:idg_end)=0.0_r8
          ENDIF
!
!     UPDATE GASEOUS, AQUEOUS GAS CONTENTS AND CONCENTRATIONS
!     FROM GASEOUS-AQUEOUS EXCHANGE, SOIL GAS TRANSFERS
!
          trc_gasml_loc(idg_O2)=trc_gasml_loc(idg_O2)-RDFQSolute(idg_O2)+ROXYFX
          trc_gasml_loc(idg_CO2)=trc_gasml_loc(idg_CO2)-RDFQSolute(idg_CO2)+RCO2FX

          DO NTG=idg_beg,idg_end
            trc_solml_loc(NTG)=trc_solml_loc(NTG)+RDFQSolute(NTG)-RUPSolute(NTG)
          enddo
!
!     GAS TRANSFER THROUGH ROOTS
!
          IF(N.EQ.ipltroot.AND.RootPoreVol_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
            RUPNTX=RUPSolute(idg_NH3)+RUPSolute(idg_NH3B)
!
!     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN ROOTS
!     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
!     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENT
!
!     R*DF1=root gas exchange between gaseous-aqueous phases
!     R*FL1=root gas exchange with atmosphere
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     trcg_gmas_loc(:),trcs_rootml_loc(:)=gaseous,aqueous gas in root
!     RootVH2O_pvr,RootPoreVol_pvr=root aqueous,gaseous volume
!     VOLW*=RTVLW*gas solubility
!     C*E,C*A1=atmosphere,root gas concentration
!     DF*A=root-atmosphere gas conductance
!
            GasPX(idg_CO2)=trcs_rootml_loc(idg_CO2)+RCO2PX
            GasPX(idg_NH3)=trcs_rootml_loc(idg_NH3)+RUPNTX

            DO NTG=idg_beg,idg_end-1
              if(NTG/=idg_CO2.and.NTG/=idg_NH3)GasPX(NTG)=trcs_rootml_loc(NTG)+RUPSolute(NTG)
              trcg_Root2Soil_flx(NTG)=AMAX1(-GasPX(NTG),DFGP*(AMAX1(ZERO4Groth_pft(NZ),trcg_gmas_loc(NTG)) &
                *DisolvedGasVolume(NTG)-GasPX(NTG)*RootPoreVol_pvr(N,L,NZ)) &
                /(DisolvedGasVolume(NTG)+RootPoreVol_pvr(N,L,NZ)))
              !positive into root, negative into atmosphere, assuming specific rate 1/hr
              trcg_air2root_flx1(NTG)=AMIN1(DFAGas(NTG),RootPoreVol_pvr(N,L,NZ)) &
                *(AtmGasc(NTG)-trcg_gcon_loc(NTG))
            enddo
          ELSE
            trcg_Root2Soil_flx(idg_beg:idg_end-1)=0.0_r8
            trcg_air2root_flx1(idg_beg:idg_end-1)=0.0_r8
          ENDIF
!
!     UPDATE ROOT AQUEOUS, GASEOUS GAS CONTENTS AND CONCENTRATIONS
!     FOR ROOT AQUEOUS-GASEOUS, GASEOUS-ATMOSPHERE EXCHANGES
!
          DO NTG=idg_beg,idg_end-1
            trcg_gmas_loc(NTG)=trcg_gmas_loc(NTG)-trcg_Root2Soil_flx(NTG)+trcg_air2root_flx1(NTG)
            trcs_rootml_loc(NTG)=trcs_rootml_loc(NTG)+trcg_Root2Soil_flx(NTG)+RUPSolute(NTG)
          ENDDO
          trcs_rootml_loc(idg_CO2)=trcs_rootml_loc(idg_CO2)+RCO2PX
          trcs_rootml_loc(idg_NH3)=trcs_rootml_loc(idg_NH3)+RUPSolute(idg_NH3B)
!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
!
!     RCO2S=soil-root gas exchange
!
          DO NTG=idg_beg,idg_end
            if(NTG==idg_O2)then
              RUPGasSol_vr(NTG,N,L,NZ)=RUPGasSol_vr(NTG,N,L,NZ)+RUPOSX
            else
              RUPGasSol_vr(NTG,N,L,NZ)=RUPGasSol_vr(NTG,N,L,NZ)+RUPSolute(NTG)
            endif  
          enddo
!
!     ACCUMULATE ROOT-ATMOSPHERE GAS EXCHANGE TO HOURLY TIME SCALE
!
!     R*DFA=root aqueous-gaseous CO2 exchange
!     R*FLA=root gaseous-atmosphere CO2 exchange
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!
          DO NTG=idg_beg,idg_end-1
            trcg_Root_DisEvap_flx_vr(NTG,N,L,NZ)=trcg_Root_DisEvap_flx_vr(NTG,N,L,NZ) &
              +trcg_Root2Soil_flx(NTG)
            trcg_air2root_flx__pvr(NTG,N,L,NZ)=trcg_air2root_flx__pvr(NTG,N,L,NZ) &
              +trcg_air2root_flx1(NTG)
          ENDDO
!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
!
!     RCO2P=root CO2 emission into root
!     RootO2Uptk_pvr=root O2 uptake from root
!     RO2UptkSoilM_vr=total O2 uptake from soil by all microbial,root popns
!
          RootCO2Emis_pvr(N,L,NZ)=RootCO2Emis_pvr(N,L,NZ)+RCO2PX+RUPSolute(idg_CO2)
          RootO2Uptk_pvr(N,L,NZ)=RootO2Uptk_pvr(N,L,NZ)-RUPSolute(idg_O2)
          RO2UptkSoilM_vr(M,L)=RO2UptkSoilM_vr(M,L)+RUPOSX

        ENDDO D90
      ENDIF
    ENDDO D99
!
!     O2 CONSTRAINTS TO ROOT RESPIRATION DEPENDS UPON RATIO
!     OF ROOT O2 UPTAKE 'RO2UptkHeterT' TO ROOT O2 DEMAND 'RootO2Dmnd4Resp_pvr'
!
!     RO2UptkHeterT=O2 uptake from soil+root by each root,myco population
!     RootO2Dmnd4Resp_pvr=O2 demand by each root,myco population
!     RAutoRootO2Limter_pvr=constraint by O2 consumption on all root processes
!     imposed by O2 uptake
!
    PopPlantO2Uptake_vr=RootO2Uptk_pvr(N,L,NZ)+RUPGasSol_vr(idg_O2,N,L,NZ)
    RAutoRootO2Limter_pvr(N,L,NZ)=AMIN1(1.0_r8,AZMAX1(PopPlantO2Uptake_vr/RootO2Dmnd4Resp_pvr(N,L,NZ)))

  ELSE
    PopPlantO2Uptake_vr=0.0_r8
    IF(L.GT.NGTopRootLayer_pft(NZ))THEN
      RAutoRootO2Limter_pvr(N,L,NZ)=RAutoRootO2Limter_pvr(N,L-1,NZ)
    ELSE
      RAutoRootO2Limter_pvr(N,L,NZ)=1.0
    ENDIF
  ENDIF
  end associate
  end subroutine RootSoilGasExchange

end module RootGasMod
