module RootGasMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
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
  real(r8), intent(in) :: dtPerPlantRootH2OUptake   !root water uptake 
  real(r8), intent(in) :: FOXYX
  real(r8), intent(out):: PopPlantO2Uptake_vr
  !local variables
  integer :: M,MX
  real(r8) :: B,C
  real(r8) :: trcg_rootml_loc(idg_beg:idg_end-1)     !vol gaseous mass volatile in roots
  real(r8) :: trcs_rootml_loc(idg_beg:idg_end-1)     !vol aqueous mass volatile in roots
  real(r8) :: trcaqu_conc_soi_loc(idg_beg:idg_end)   !aqueous volatile tracer concentration in soil
  real(r8) :: trc_conc_root_loc(idg_beg:idg_end-1)   !aqueous tracer concentration in roots  [g m-3]
  real(r8) :: GasDifc_loc(idg_beg:idg_end-1)
  real(r8) :: SolDifc_loc(idg_beg:idg_end-1)
  real(r8) :: trcg_gcon_loc(idg_beg:idg_end-1)       !gaseous volatile concentration in roots 
  real(r8) :: COXYR                             !aqueous O2 concentration at root surface used for uptake
  real(r8) :: trcs_maxRootml_loc(idg_beg:idg_end-1)
  real(r8) :: DIFOP                             !aqueous diffusivity of O2 within root
  real(r8) :: DifAqueVolatile(idg_beg:idg_end)  !aqueous diffusivity from soil to root
  real(r8) :: DIFOX,DiffusivitySolutEffP
  real(r8) :: DFAGas(idg_beg:idg_end-1)         !maximum volatile flux from atmosphere to roots
  real(r8) :: DFGP,ROxySoil2Uptk
  real(r8) :: trc_solml_loc(idg_beg:idg_end)
  real(r8) :: trc_gasml_loc(idg_beg:idg_end-1)
  real(r8) :: O2AquaDiffusvityP,RTVLWA
  real(r8) :: RTVLWB,ROXYFX,RCO2FX,ROXYLX
  real(r8) :: RDFQSolute(idg_beg:idg_end)   !gas dissolution into aqueous concentration
  real(r8) :: RTCR1,RTCR2
  real(r8) :: RTCRA,RTARRX
  real(r8) :: RCO2PX                       !root CO2 gas efflux due to root respiration at time step for gas flux calculations
  real(r8) :: RRADS
  real(r8) :: RSolUptkTransp(idg_beg:idg_end)   !volatile tracer uptake by transpiration
  real(r8) :: RootOxyUptakePerPlant
  real(r8) :: ROxySoil2UptkPerPlant            !aqueous O2 uptake flux from soil due to diffusion and transpiration-aided-advection
  real(r8) :: ROxyRoot2UptkPerPlant            !aqueous O2 diffusion uptake flux from inside roots
  real(r8) :: ROxyRoot2Uptk  
  real(r8) :: RootUptkSoiSolute(idg_beg:idg_end)  !root uptake of aqueous volatile from soil to inside roots
  real(r8) :: RDXSolute(idg_beg:idg_end)    !maximum fluxes can be taken by roots
  real(r8) :: RDFSolute(idg_beg:idg_end)
  real(r8) :: RUPOST
  real(r8) :: RUPNTX  !total uptake of NH4 from soil into roots
  real(r8) :: Root_gas2sol_flx(idg_beg:idg_end-1)    !gas dissolution into aqueous phase of the volatile tracers
  real(r8) :: trcg_air2root_flx(idg_beg:idg_end-1)   !diffusion flux of gas from atmosphere to inside roots
  real(r8) :: THETW1,THETM
  real(r8) :: RootOxyDemandPerPlant
  real(r8) :: DisolvedGasVolume(idg_beg:idg_end-1),VLWatMicPMO,VLWatMicPMM,VLsoiAirPMM
  real(r8) :: VOLWSP,VLWatMicPMA,VLWatMicPMB,VOLWSA,VOLWSB,VOLWSolute(idg_beg:idg_end)
  real(r8) :: VOLPNH3,VOLPNH3B
  real(r8) :: X
  real(r8) :: ZH3PA,ZH3PB,ZH3GA,ZH3GB
  integer  :: NTG
!     begin_execution
  associate(                                                       &
    RootStrutElms_pft        => plt_biom%RootStrutElms_pft,        &
    ZERO4Groth_pft           => plt_biom%ZERO4Groth_pft,           &
    PlantPopulation_pft      => plt_site%PlantPopulation_pft,      &
    DPTHZ_vr                 => plt_site%DPTHZ_vr,                 &
    AtmGasc                  => plt_site%AtmGasc,                  &  !in: atmospheric gaseous concentration 
    ZEROS                    => plt_site%ZEROS,                    &
    ZERO                     => plt_site%ZERO,                     &
    VLWatMicPM_vr            => plt_site%VLWatMicPM_vr,            &
    VLsoiAirPM               => plt_site%VLsoiAirPM,               &
    TortMicPM_vr             => plt_site%TortMicPM_vr,             &
    FILM                     => plt_site%FILM,                     &
    RO2GasXchangePrev_vr     => plt_bgcr%RO2GasXchangePrev_vr,     &
    RCO2GasFlxPrev_vr        => plt_bgcr%RCO2GasFlxPrev_vr,        &
    RO2AquaSourcePrev_vr    => plt_bgcr%RO2AquaSourcePrev_vr,    &
    RootO2Uptk_pvr           => plt_rbgc%RootO2Uptk_pvr,           & !out: O2 uptake from O2 inside root
    RAutoRootO2Limter_rpvr    => plt_rbgc%RAutoRootO2Limter_rpvr,    &
    ZERO4Uptk_pft            => plt_rbgc%ZERO4Uptk_pft,            &
    RootRespPotent_pvr       => plt_rbgc%RootRespPotent_pvr,       &
    RootO2Dmnd4Resp_pvr      => plt_rbgc%RootO2Dmnd4Resp_pvr,      &
    RO2UptkSoilM_vr          => plt_rbgc%RO2UptkSoilM_vr,          &
    RootCO2Emis_pvr          => plt_rbgc%RootCO2Emis_pvr,          & !out: total CO2 emitted inside roots
    trcg_air2root_flx_pvr    => plt_rbgc%trcg_air2root_flx_pvr,    &
    trcg_Root_gas2aqu_flx_vr => plt_rbgc%trcg_Root_gas2aqu_flx_vr, &
    RootUptkSoiSol_vr        => plt_rbgc%RootUptkSoiSol_vr,        &  !out: aqueous tracer uptake from soil into roots
    RootCO2Autor_pvr         => plt_rbgc%RootCO2Autor_pvr,         &
    trcg_rootml_pvr          => plt_rbgc%trcg_rootml_pvr,          &
    trcs_rootml_pvr          => plt_rbgc%trcs_rootml_pvr,          &
    TScal4Difsvity_vr        => plt_soilchem%TScal4Difsvity_vr,    &
    trcs_VLN_vr              => plt_soilchem%trcs_VLN_vr,          &
    trc_solml_vr             => plt_soilchem%trc_solml_vr,         &
    trc_gascl_vr             => plt_soilchem%trc_gascl_vr,         &
    GasDifc_vr               => plt_soilchem%GasDifc_vr,           &  !in: gaseous diffusivity of volatile tracers
    SoluteDifusvty_vr        => plt_soilchem%SoluteDifusvty_vr,    &  !in: aqueous diffusivity of volatile tracers
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
!     trcg_rootml_loc(idg_CO2),trcs_rootml_loc(idg_CO2),trc_gasml_loc(idg_CO2),trc_solml_loc(idg_CO2)=gaseous,aqueous CO2 in root,soil
!     trcaqu_conc_soi_loc(idg_CH4),trc_conc_root_loc(idg_CH4)=aqueous gas concentration in soil,root
!     RTVLWA,RTVLWB=root aqueous volume in non-band,band
!     dts_gas=time step of flux calculation
!     RootOxyDemandPerPlant=O2 demand per plant at time step of flux calculation
!     ROXYFX=net O2 gas flux at time step of flux calculation
!     RCO2FX=net CO2 gas flux at time step of flux calculation
!     ROXYLX=net O2 aqueous flux at time step of flux calculation
!

    DO NTG=idg_beg,idg_end-1
      trcg_rootml_loc(NTG) = AMAX1(ZERO4Groth_pft(NZ),trcg_rootml_pvr(NTG,N,L,NZ))
      trcs_rootml_loc(NTG) = AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_pvr(NTG,N,L,NZ))
    ENDDO

    DO NTG=idg_beg,idg_end
      if(NTG==idg_CO2)then
        trc_solml_loc(NTG)=AMAX1(ZERO4Groth_pft(NZ),trc_solml_vr(NTG,L)*FracPRoot4Uptake(N,L,NZ))
      elseif(NTG/=idg_O2)then
        trc_solml_loc(NTG)=trc_solml_vr(NTG,L)*FracPRoot4Uptake(N,L,NZ)
      endif
    enddo
    trc_solml_loc(idg_O2)=trc_solml_vr(idg_O2,L)*FOXYX

!  the two lines below may be redundant
    trc_gasml_loc(idg_CO2) = AMAX1(ZERO4Groth_pft(NZ),trc_gasml_vr(idg_CO2,L)*FracPRoot4Uptake(N,L,NZ))
    trc_gasml_loc(idg_O2)  = AMAX1(ZERO4Groth_pft(NZ),trc_gasml_vr(idg_O2,L)*FOXYX)

    RTVLWA                = RootVH2O_pvr(N,L,NZ)*trcs_VLN_vr(ids_NH4,L)
    RTVLWB                = RootVH2O_pvr(N,L,NZ)*trcs_VLN_vr(ids_NH4B,L)
    RootOxyDemandPerPlant = RootO2Dmnd4Resp_pvr(N,L,NZ)*dts_gas/PlantPopulation_pft(NZ)
    ROXYFX                = RO2GasXchangePrev_vr(L)*FOXYX*dts_gas
    RCO2FX                = RCO2GasFlxPrev_vr(L)*FOXYX*dts_gas
    ROXYLX                = RO2AquaSourcePrev_vr(L)*FOXYX*dts_gas

!
!     GASEOUS AND AQUEOUS DIFFUSIVITIES IN ROOT AND SOIL
!
!     *SGL1=diffusivity
!     RootPoreTortu4Gas=tortuosity effect of root porosity on diffusivity
!     CG=CO2g,OG=O2g,CH=CH4g,Z2=N2Og,ZH=NH3g,HG=H2g
!     CL=CO2s,OL=O2s,CQ=CH4s,ZV=N2Os,ZN=NH3s,HL=H2s
!
    do NTG=idg_beg,idg_end-1
      GasDifc_loc(NTG) = GasDifc_vr(NTG,L)*dts_gas*RootPoreTortu4Gas(N,NZ)
      SolDifc_loc(NTG) = SoluteDifusvty_vr(NTG,L)*dts_gas*FOXYX
    enddo
    O2AquaDiffusvityP = SoluteDifusvty_vr(idg_O2,L)*dts_gas

    RDFQSolute(idg_beg:idg_end)         = 0.0_r8
    Root_gas2sol_flx(idg_beg:idg_end-1) = 0.0_r8
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
    IF(RootStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ) .AND. FracSoiLayByPrimRoot(L,NZ).GT.ZERO)THEN
      RTCR1 = AMAX1(PlantPopulation_pft(NZ),Root1stXNumL_pvr(N,L,NZ))*PICON*Root1stRadius_pvr(N,L,NZ)**2/DPTHZ_vr(L)
      RTCR2 = (Root2ndXNum_pvr(N,L,NZ)*PICON*Root2ndRadius_pvr(N,L,NZ)**2/Root2ndAveLen_pvr(N,L,NZ))/FracSoiLayByPrimRoot(L,NZ)
      IF(RTCR2.GT.RTCR1)THEN
        RTCRA = RTCR1*RTCR2/(RTCR1+RTCR2)
      ELSE
        RTCRA = RTCR1
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
!     DIFOP=
!     RootVH2O_pvr=root,myco aqueous volume
!     S*L=solubility of gas in water from hour1.f:
!     CO2=CO2,OXY=O2,CH4=CH4,N2O=N2O,NH3=NH3,H2G=H2
!     DF*A=root-atmosphere gas conductance
!     DFGP=rate const for equilibrn of gas concn in gaseous-aqueous phases
!     RootCO2Autor_pvr=root CO2 flux from grosub.f
!

    IF(N.EQ.ipltroot .AND. iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).GT.0 &
      .AND. RootLenPerPlant_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
      RTARRX = RootAreaDivRadius_vr(N,L)/RootRaidus_rpft(N,NZ)
      DIFOP  = O2AquaDiffusvityP*RTARRX
      DO NTG    = idg_beg, idg_end-1
        DisolvedGasVolume(NTG) = RootVH2O_pvr(N,L,NZ)*GasSolbility_vr(NTG,L)
        DFAGas(NTG)            = GasDifc_loc(NTG)*RTCRA
      ENDDO
    ELSE
      RTARRX                               = 0.0_r8
      DIFOP                                = 0.0_r8
      DisolvedGasVolume(idg_beg:idg_end-1) = 0.0_r8
      DFAGas(idg_beg:idg_end-1)            = 0.0_r8
      DFAGas(idg_O2)                       = AMIN1(1.e-12_r8,RootPoreVol_pvr(N,L,NZ))
    ENDIF

    DFGP   = AMIN1(1.0,XNPD*SQRT(RootPorosity_pft(N,NZ))*TScal4Difsvity_vr(L))
    RCO2PX = -RootCO2Autor_pvr(N,L,NZ)*dts_gas
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
      VLWatMicPMO = VLWatMicPM_vr(M,L)*FOXYX
      VLWatMicPMM = VLWatMicPM_vr(M,L)*FracPRoot4Uptake(N,L,NZ)
      VLsoiAirPMM = VLsoiAirPM(M,L)*FracPRoot4Uptake(N,L,NZ)
      VOLWSP      = RootVH2O_pvr(N,L,NZ)+VLWatMicPMM
      VLWatMicPMA = VLWatMicPMM*trcs_VLN_vr(ids_NH4,L)
      VLWatMicPMB = VLWatMicPMM*trcs_VLN_vr(ids_NH4B,L)
      VOLWSA      = RTVLWA+VLWatMicPMA
      VOLWSB      = RTVLWB+VLWatMicPMB
      THETW1      = AZMAX1(VLWatMicPM_vr(M,L)/VLSoilMicP_vr(L))

      IF(THETW1.GT.THETY_vr(L) .AND. FracPRoot4Uptake(N,L,NZ).GT.ZERO4Uptk_pft(NZ))THEN
        THETM  = TortMicPM_vr(M,L)*THETW1
        RRADS  = LOG((FILM(M,L)+FineRootRadius(N,L))/FineRootRadius(N,L))
        RTARRX = RootAreaDivRadius_vr(N,L)/RRADS
        do NTG    = idg_beg, idg_end-1
          DifAqueVolatile(NTG)=THETM*SolDifc_loc(NTG)*RTARRX
          if(NTG/=idg_O2 .and. NTG/=idg_CO2)then
            trc_gasml_loc(NTG)=trc_gascl_vr(NTG,L)*VLsoiAirPMM
          endif          
        enddo
        DifAqueVolatile(idg_NH3)  = DifAqueVolatile(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
        DifAqueVolatile(idg_NH3B) = DifAqueVolatile(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)

        DO NTG=idg_beg,idg_end
          VOLWSolute(NTG)=VLWatMicPMM*GasSolbility_vr(NTG,L)
        ENDDO
        VOLWSolute(idg_NH3)  = VOLWSolute(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
        VOLWSolute(idg_NH3B) = VOLWSolute(idg_NH3B)*trcs_VLN_vr(ids_NH4B,L)

        VOLPNH3  = VLsoiAirPMM*trcs_VLN_vr(ids_NH4,L)
        VOLPNH3B = VLsoiAirPMM*trcs_VLN_vr(ids_NH4B,L)
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
              trcg_gcon_loc(NTG)=AZMAX1(trcg_rootml_loc(NTG)/RootPoreVol_pvr(N,L,NZ))
            ENDDO
          ELSE
            trcg_gcon_loc(idg_beg:idg_end-1)=0.0_r8
          ENDIF

          do NTG=idg_beg,idg_end-1
            trc_conc_root_loc(NTG)=AZMAX1(trcs_rootml_loc(NTG)/RootVH2O_pvr(N,L,NZ))
          enddo
          trc_conc_root_loc(idg_O2)=AMIN1(AtmGasc(idg_O2)*GasSolbility_vr(idg_O2,L),trc_conc_root_loc(idg_O2))
          DIFOX                    =DifAqueVolatile(idg_O2)+DIFOP
          do NTG=idg_beg,idg_end
            RSolUptkTransp(NTG)=dtPerPlantRootH2OUptake*trcaqu_conc_soi_loc(NTG)
          enddo
          RSolUptkTransp(idg_NH3)  = RSolUptkTransp(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
          RSolUptkTransp(idg_NH3B) = RSolUptkTransp(idg_NH3B)*trcs_VLN_vr(ids_NH4B,L)
!
!     SOLUTION FOR MASS FLOW + DIFFUSION OF O2 IN AQUEOUS PHASES OF
!     SOIL AND ROOT = ACTIVE UPTAKE OF O2 BY ROOT
!
!     DifAqueVolatile(idg_O2)=O2 
!     dtPerPlantRootH2OUptake=water uptake
!     DIFOP=aqueous diffusivity of O2 within root
!     trcaqu_conc_soi_loc(idg_O2),trc_conc_root_loc(idg_O2)=soil,root aqueous O2 concentration
!     RootOxyDemandPerPlant=O2 demand per plant
!     RootOxyUptakePerPlant=root O2 uptake per plant
!     COXYR: 
!
          X=(DifAqueVolatile(idg_O2)+dtPerPlantRootH2OUptake)*trcaqu_conc_soi_loc(idg_O2)+DIFOP*trc_conc_root_loc(idg_O2)

          IF(X.GT.ZERO .AND. trc_solml_loc(idg_O2).GT.ZERO4Groth_pft(NZ))THEN
            !root take up O2
            B                     = -RootOxyDemandPerPlant-DIFOX*OXKM-X
            C                     = X*RootOxyDemandPerPlant
            RootOxyUptakePerPlant = (-B-SQRT(B*B-4.0_r8*C))/2.0_r8
            COXYR                 = (X-RootOxyUptakePerPlant)/DIFOX
            ROxySoil2UptkPerPlant = RSolUptkTransp(idg_O2)+DifAqueVolatile(idg_O2)*(trcaqu_conc_soi_loc(idg_O2)-COXYR)
            ROxyRoot2UptkPerPlant = DIFOP*(trc_conc_root_loc(idg_O2)-COXYR)
          ELSE
            !soil O2 concentration is too low, uptake is from inside root
            X=DIFOP*trc_conc_root_loc(idg_O2)
            IF(X.GT.ZERO .AND. trcs_rootml_loc(idg_O2).GT.ZERO4Groth_pft(NZ))THEN
              B                     = -RootOxyDemandPerPlant-DIFOP*OXKM-X
              C                     = X*RootOxyDemandPerPlant
              RootOxyUptakePerPlant = (-B-SQRT(B*B-4.0_r8*C))/2.0_r8
              COXYR                 = (X-RootOxyUptakePerPlant)/DIFOP
              ROxySoil2UptkPerPlant = 0.0_r8
              ROxyRoot2UptkPerPlant = DIFOP*(trc_conc_root_loc(idg_O2)-COXYR)
            ELSE
              RootOxyUptakePerPlant = 0.0_r8
              COXYR                 = 0.0_r8
              ROxySoil2UptkPerPlant = 0.0_r8
              ROxyRoot2UptkPerPlant = 0.0_r8
            ENDIF
          ENDIF
!          write(113,*)I*1000+J,MX,M,RootOxyUptakePerPlant,ROxySoil2UptkPerPlant,DifAqueVolatile(idg_O2),DIFOP, &
!            trcaqu_conc_soi_loc(idg_O2),trc_conc_root_loc(idg_O2),trc_solml_loc(idg_O2),trc_solml_vr(idg_O2,L),ROXYLX
!
!     MASS FLOW + DIFFUSIVE EXCHANGE OF OTHER GASES
!     BETWEEN ROOT AND SOIL, CONSTRAINED BY COMPETITION
!     WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!     ROxySoil2Uptk,ROxyRoot2Uptk=aqueous O2 uptake from soil,root
!     PP=PFT population
!     RDFSolute(:),RootUptkSoiSolute(:)=aqueous gas soil-root diffusion,root uptake, > 0 into soil
!     RMF*=soil convective solute flux
!     DIF*=aqueous diffusivity from soil to root
!     C*S1=soil aqueous concentration non-band
!     C*B1=soil aqueous concentration band
!     C*P1=root aqueous concentration
!
          ROxySoil2Uptk      = ROxySoil2UptkPerPlant*PlantPopulation_pft(NZ)
          ROxyRoot2Uptk      = ROxyRoot2UptkPerPlant*PlantPopulation_pft(NZ)
          RDFSolute(idg_CO2) = RSolUptkTransp(idg_CO2)+DifAqueVolatile(idg_CO2)*(trcaqu_conc_soi_loc(idg_CO2)-trc_conc_root_loc(idg_CO2))
          RDXSolute(idg_CO2) = (RootVH2O_pvr(N,L,NZ)*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_CO2)) &
            -VLWatMicPMM*AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_loc(idg_CO2)))/VOLWSP
          
          IF(RDFSolute(idg_CO2).GT.0.0_r8)THEN
            RootUptkSoiSolute(idg_CO2)=AMIN1(AZMAX1(RDXSolute(idg_CO2)),RDFSolute(idg_CO2)*PlantPopulation_pft(NZ))
          ELSE
            RootUptkSoiSolute(idg_CO2)=AMAX1(AZMIN1(RDXSolute(idg_CO2)),RDFSolute(idg_CO2)*PlantPopulation_pft(NZ))
          ENDIF

          IF(N.EQ.ipltroot)THEN
            !fluxes only involve roots
            DO NTG=idg_beg,idg_end-1
              if(NTG/=idg_CO2 .and. NTG/=idg_NH3 .and. NTG/=idg_O2)then
                !soil to root fluxes
                RDFSolute(NTG) = RSolUptkTransp(NTG)+DifAqueVolatile(NTG)*(trcaqu_conc_soi_loc(NTG)-trc_conc_root_loc(NTG))
                RDXSolute(NTG) = (RootVH2O_pvr(N,L,NZ)*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(NTG)) &
                  -VLWatMicPMM*AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_loc(NTG)))/VOLWSP

                IF(RDFSolute(NTG).GT.0.0_r8)THEN
                  RootUptkSoiSolute(NTG)=AMIN1(AZMAX1(RDXSolute(NTG)),RDFSolute(NTG)*PlantPopulation_pft(NZ))
                ELSE
                  RootUptkSoiSolute(NTG)=AMAX1(AZMIN1(RDXSolute(NTG)),RDFSolute(NTG)*PlantPopulation_pft(NZ))
                ENDIF
              ENDIF
            ENDDO

            RDFSolute(idg_NH3)=RSolUptkTransp(idg_NH3)+DifAqueVolatile(idg_NH3)*(trcaqu_conc_soi_loc(idg_NH3)-trc_conc_root_loc(idg_NH3))

            IF(VOLWSA.GT.ZERO4Groth_pft(NZ))THEN
              ZH3PA              = trcs_rootml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
              RDXSolute(idg_NH3) = (RTVLWA*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_NH3)) &
                -VLWatMicPMA*AMAX1(ZERO4Groth_pft(NZ),ZH3PA))/VOLWSA
            ELSE
              RDXSolute(idg_NH3)=0.0_r8
            ENDIF

            IF(RDFSolute(idg_NH3).GT.0.0_r8)THEN
              RootUptkSoiSolute(idg_NH3)=AMIN1(AZMAX1(RDXSolute(idg_NH3)),RDFSolute(idg_NH3)*PlantPopulation_pft(NZ))
            ELSE
              RootUptkSoiSolute(idg_NH3)=AMAX1(AZMIN1(RDXSolute(idg_NH3)),RDFSolute(idg_NH3)*PlantPopulation_pft(NZ))
            ENDIF

            RDFSolute(idg_NH3B)=RSolUptkTransp(idg_NH3B)+DifAqueVolatile(idg_NH3B)*(trcaqu_conc_soi_loc(idg_NH3B)-trc_conc_root_loc(idg_NH3))
            IF(VOLWSB.GT.ZERO4Groth_pft(NZ))THEN
              ZH3PB               = trcs_rootml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)
              RDXSolute(idg_NH3B) = (RTVLWB*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_NH3B)) &
                -VLWatMicPMB*AMAX1(ZERO4Groth_pft(NZ),ZH3PB))/VOLWSB
            ELSE
              RDXSolute(idg_NH3B)=0.0_r8
            ENDIF

            IF(RDFSolute(idg_NH3B).GT.0.0_r8)THEN
              RootUptkSoiSolute(idg_NH3B)=AMIN1(AZMAX1(RDXSolute(idg_NH3B)),RDFSolute(idg_NH3B)*PlantPopulation_pft(NZ))
            ELSE
              RootUptkSoiSolute(idg_NH3B)=AMAX1(AZMIN1(RDXSolute(idg_NH3B)),RDFSolute(idg_NH3B)*PlantPopulation_pft(NZ))
            ENDIF
          ELSE
            RootUptkSoiSolute(idg_beg:idg_end)=0.0_r8
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
!     RootUptkSoiSolute=root aqueous gas uptake
!     ROXYLX=soil net O2 aqueous flux

!     VLWatMicPMM,VLsoiAirPMM=soil micropore water,air volume
!     VOLW*=VLWatMicPMM*gas solubility
!         
          IF(THETPM(M,L).GT.THETX)THEN
            DiffusivitySolutEffP = FracPRoot4Uptake(N,L,NZ)*DiffusivitySolutEff(M,L)
            RDFQSolute(idg_CO2)  = DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(idg_CO2))*VOLWSolute(idg_CO2) &
              -(AMAX1(ZEROS,trc_solml_loc(idg_CO2))-RootUptkSoiSolute(idg_CO2))*VLsoiAirPMM) &
              /(VOLWSolute(idg_CO2)+VLsoiAirPMM)

            RUPOST             = ROxySoil2Uptk-ROXYLX
            RDFQSolute(idg_O2) = DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(idg_O2))*VOLWSolute(idg_O2) &
              -(AMAX1(ZEROS,trc_solml_loc(idg_O2))-RUPOST)*VLsoiAirPMM)/(VOLWSolute(idg_O2)+VLsoiAirPMM)

            IF(N.EQ.ipltroot)THEN
              DO NTG=idg_beg,idg_NH3-1
                if(NTG/=idg_CO2 .and. NTG/=idg_O2)then
                  RDFQSolute(NTG)=DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(NTG))*VOLWSolute(NTG) &
                    -(AMAX1(ZEROS,trc_solml_loc(NTG))-RootUptkSoiSolute(NTG))*VLsoiAirPMM)/(VOLWSolute(NTG)+VLsoiAirPMM)
                endif
              ENDDO

              IF(VOLWSolute(idg_NH3)+VOLPNH3.GT.ZERO4Groth_pft(NZ))THEN
                ZH3GA               = trc_gasml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
                RDFQSolute(idg_NH3) = AMIN1(RootUptkSoiSolute(idg_NH3),AMAX1(-RootUptkSoiSolute(idg_NH3) &
                  ,DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),ZH3GA)*VOLWSolute(idg_NH3) &
                  -(AMAX1(ZEROS,trc_solml_loc(idg_NH3))-RootUptkSoiSolute(idg_NH3))*VOLPNH3) &
                  /(VOLWSolute(idg_NH3)+VOLPNH3)))
              ELSE
                RDFQSolute(idg_NH3)=0.0_r8
              ENDIF

              IF(VOLWSolute(idg_NH3B)+VOLPNH3B.GT.ZERO4Groth_pft(NZ))THEN
                ZH3GB                = trc_gasml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)
                RDFQSolute(idg_NH3B) = AMIN1(RootUptkSoiSolute(idg_NH3),AMAX1(-RootUptkSoiSolute(idg_NH3) &
                  ,DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),ZH3GB)*VOLWSolute(idg_NH3B)   &
                  -(AMAX1(ZEROS,trc_solml_loc(idg_NH3B))-RootUptkSoiSolute(idg_NH3B))*VOLPNH3B)  &
                  /(VOLWSolute(idg_NH3B)+VOLPNH3B)))
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
          trc_gasml_loc(idg_O2)  = trc_gasml_loc(idg_O2)-RDFQSolute(idg_O2)+ROXYFX
          trc_gasml_loc(idg_CO2) = trc_gasml_loc(idg_CO2)-RDFQSolute(idg_CO2)+RCO2FX

          !aqueous concentration in soil
          DO NTG=idg_beg,idg_end
            if(NTG/=idg_O2)then
              trc_solml_loc(NTG)=trc_solml_loc(NTG)+RDFQSolute(NTG)-RootUptkSoiSolute(NTG)
            else
              trc_solml_loc(NTG)=trc_solml_loc(NTG)+RDFQSolute(NTG)-ROxySoil2Uptk
            endif
          enddo
!
!     GAS TRANSFER THROUGH ROOTS
!
          IF(N.EQ.ipltroot.AND.RootPoreVol_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
            RUPNTX=RootUptkSoiSolute(idg_NH3)+RootUptkSoiSolute(idg_NH3B)
!
!     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN ROOTS
!     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
!     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENT
!
!     R*DF1=root gas exchange between gaseous-aqueous phases
!     R*FL1=root gas exchange with atmosphere
!     gas code:CO=CO2,OX=O2,CH=CH4,N2=N2O,NH=NH3,H2=H2
!     trcg_rootml_loc(:),trcs_rootml_loc(:)=gaseous,aqueous gas in root
!     RootVH2O_pvr,RootPoreVol_pvr=root aqueous,gaseous volume
!     VOLW*=RTVLW*gas solubility
!     C*E,C*A1=atmosphere,root gas concentration
!     DF*A=root-atmosphere gas conductance
!
            trcs_maxRootml_loc(idg_CO2) = trcs_rootml_loc(idg_CO2)+RCO2PX
            trcs_maxRootml_loc(idg_NH3) = trcs_rootml_loc(idg_NH3)+RUPNTX

            DO NTG=idg_beg,idg_end-1
              if(NTG/=idg_CO2 .and. NTG/=idg_NH3)then
                if(NTG==idg_O2)then
                  trcs_maxRootml_loc(NTG)=trcs_rootml_loc(NTG)-ROxyRoot2Uptk                
                else
                  trcs_maxRootml_loc(NTG)=trcs_rootml_loc(NTG)+RootUptkSoiSolute(NTG)
                endif
              endif
              Root_gas2sol_flx(NTG)=AMAX1(-trcs_maxRootml_loc(NTG),DFGP*(AMAX1(ZERO4Groth_pft(NZ),trcg_rootml_loc(NTG)) &
                *DisolvedGasVolume(NTG)-trcs_maxRootml_loc(NTG)*RootPoreVol_pvr(N,L,NZ)) &
                /(DisolvedGasVolume(NTG)+RootPoreVol_pvr(N,L,NZ)))
              !>0._r8 into root, <0._r8 into atmosphere, assuming specific rate 1/hr
              trcg_air2root_flx(NTG)=AMIN1(DFAGas(NTG),RootPoreVol_pvr(N,L,NZ))*(AtmGasc(NTG)-trcg_gcon_loc(NTG))
            enddo
          ELSE
            Root_gas2sol_flx(idg_beg:idg_end-1)  = 0.0_r8
            trcg_air2root_flx(idg_beg:idg_end-1) = 0.0_r8
          ENDIF
!
!     UPDATE ROOT AQUEOUS, GASEOUS GAS CONTENTS AND CONCENTRATIONS
!     FOR ROOT AQUEOUS-GASEOUS, GASEOUS-ATMOSPHERE EXCHANGES
!
          DO NTG=idg_beg,idg_end-1
            trcg_rootml_loc(NTG) = trcg_rootml_loc(NTG)-Root_gas2sol_flx(NTG)+trcg_air2root_flx(NTG)
            if(NTG==idg_O2)then
              trcs_rootml_loc(NTG) = trcs_rootml_loc(NTG)+Root_gas2sol_flx(NTG)-ROxyRoot2Uptk      !O2 is consumed.     
!              write(114,*)I+J/24.,MX,M,N,trcs_rootml_loc(NTG),trcg_rootml_loc(NTG),trcg_air2root_flx(NTG),Root_gas2sol_flx(NTG),&
!                DFAGas(NTG),RootPoreVol_pvr(N,L,NZ),ROxySoil2Uptk
            else
              trcs_rootml_loc(NTG) = trcs_rootml_loc(NTG)+Root_gas2sol_flx(NTG)+RootUptkSoiSolute(NTG)              
            endif
          ENDDO
          !releas autotrophic respiration CO2 into roots
          trcs_rootml_loc(idg_CO2) = trcs_rootml_loc(idg_CO2)+RCO2PX
          trcs_rootml_loc(idg_NH3) = trcs_rootml_loc(idg_NH3)+RootUptkSoiSolute(idg_NH3B)
!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
!
!     RCO2S=soil-root gas exchange
!
          DO NTG=idg_beg,idg_end
            if(NTG==idg_O2)then
              RootUptkSoiSol_vr(NTG,N,L,NZ)=RootUptkSoiSol_vr(NTG,N,L,NZ)+ROxySoil2Uptk
            else
              RootUptkSoiSol_vr(NTG,N,L,NZ)=RootUptkSoiSol_vr(NTG,N,L,NZ)+RootUptkSoiSolute(NTG)
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
            trcg_Root_gas2aqu_flx_vr(NTG,N,L,NZ) = trcg_Root_gas2aqu_flx_vr(NTG,N,L,NZ)+Root_gas2sol_flx(NTG)
            trcg_air2root_flx_pvr(NTG,N,L,NZ)    = trcg_air2root_flx_pvr(NTG,N,L,NZ)+trcg_air2root_flx(NTG)
          ENDDO
!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE
!
!     RCO2P=root CO2 emission into root
!     RootO2Uptk_pvr=root O2 uptake from root
!     RO2UptkSoilM_vr=total O2 uptake from soil by all microbial,root popns
!
          RootCO2Emis_pvr(N,L,NZ) = RootCO2Emis_pvr(N,L,NZ)+RCO2PX+RootUptkSoiSolute(idg_CO2)
          RootO2Uptk_pvr(N,L,NZ)  = RootO2Uptk_pvr(N,L,NZ)+ROxyRoot2Uptk  !uptake from O2 in roots
          RO2UptkSoilM_vr(M,L)    = RO2UptkSoilM_vr(M,L)+ROxySoil2Uptk    !uptake from soil O2 for respiration
        ENDDO D90
      ENDIF
    ENDDO D99
!
!     O2 CONSTRAINTS TO ROOT RESPIRATION DEPENDS UPON RATIO
!     OF ROOT O2 UPTAKE 'RO2UptkHeterT' TO ROOT O2 DEMAND 'RootO2Dmnd4Resp_pvr'
!
!     RO2UptkHeterT=O2 uptake from soil+root by each root,myco population
!     RootO2Dmnd4Resp_pvr=O2 demand by each root,myco population
!     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all root processes
!     imposed by O2 uptake
!
    PopPlantO2Uptake_vr           = RootO2Uptk_pvr(N,L,NZ)+RootUptkSoiSol_vr(idg_O2,N,L,NZ)
    RAutoRootO2Limter_rpvr(N,L,NZ) = AMIN1(1.0_r8,AZMAX1(PopPlantO2Uptake_vr/RootO2Dmnd4Resp_pvr(N,L,NZ)))
  ELSE
    PopPlantO2Uptake_vr=0.0_r8
    IF(L.GT.NGTopRootLayer_pft(NZ))THEN
      RAutoRootO2Limter_rpvr(N,L,NZ)=RAutoRootO2Limter_rpvr(N,L-1,NZ)
    ELSE
      RAutoRootO2Limter_rpvr(N,L,NZ)=1.0
    ENDIF
  ENDIF
  end associate
  end subroutine RootSoilGasExchange

end module RootGasMod
