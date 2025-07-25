module RootGasMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: safe_adb, vapsat, AZMAX1, AZMIN1,fixEXConsumpFlux
  use DebugToolMod
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
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine RootSoilGasExchange(I,J,N,L,NZ,FineRootRadius,FracPRoot4Uptake,FracSoiLayByPrimRoot,&
    RootAreaDivRadius_vr,dtPerPlantRootH2OUptake,FOXYX,trc_gasml_loc,trc_solml_loc,PopPlantO2Uptake)
  !
  !In the gas exchange between soil and roots, only O2 and CO2 involve active biochemical production/consumption
  !the other gases undergo physical exchange according to different water volumes inside/outside the roots.
  !O2 uptake includes that from inside the root through diffusion, and that directly from the soil through transpiration+diffusion.
  !CO2 is first released into roots, then diffuse into soil and atmosphere 
  implicit none
  integer , intent(in) :: I,J,N,L,NZ
  real(r8), intent(in) :: FineRootRadius(jroots,JZ1),FracPRoot4Uptake(jroots,JZ1,JP1)
  real(r8), intent(in) :: FracSoiLayByPrimRoot(JZ1,JP1)
  real(r8), intent(in) :: RootAreaDivRadius_vr(jroots,JZ1)
  real(r8), intent(in) :: dtPerPlantRootH2OUptake   !root water uptake 
  real(r8), intent(in) :: FOXYX
  real(r8), intent(inout) :: trc_solml_loc(idg_beg:idg_end)    !local copy of aqueous phase of the volatile tracers
  real(r8), intent(inout) :: trc_gasml_loc(idg_beg:idg_NH3)    !local copy of gaesous phase of the volatile tracers  
  real(r8), intent(out):: PopPlantO2Uptake

  character(len=*), parameter :: subname='RootSoilGasExchange'
  !local variables
  integer :: M,MX
  real(r8) :: B,C
  real(r8) :: trcg_rootml_loc(idg_beg:idg_NH3)     !vol gaseous mass volatile in roots
  real(r8) :: trcs_rootml_loc(idg_beg:idg_NH3)     !vol aqueous mass volatile in roots
  real(r8) :: trcaqu_conc_soi_loc(idg_beg:idg_end)   !aqueous volatile tracer concentration in soil
  real(r8) :: trc_conc_root_loc(idg_beg:idg_NH3)   !aqueous tracer concentration in roots  [g m-3]
  real(r8) :: GasDifc_tscaled(idg_beg:idg_NH3)
  real(r8) :: SolDifc_tscaled(idg_beg:idg_NH3)
  real(r8) :: trcg_gcon_loc(idg_beg:idg_NH3)       !gaseous volatile concentration in roots 
  real(r8) :: COXYR                                !aqueous O2 concentration at root surface used for uptake [gO/m3]
  real(r8) :: trcs_maxRootml_loc(idg_beg:idg_NH3)
  real(r8) :: DIFOP                             !aqueous diffusivity of O2 within root
  real(r8) :: DifAqueVolatile(idg_beg:idg_end)  !aqueous diffusivity from soil to root
  real(r8) :: DIFOX,DiffusivitySolutEffP
  real(r8) :: DFAGas(idg_beg:idg_NH3)                  !maximum volatile flux from atmosphere to roots
  real(r8) :: DFGP,ROxySoil2Uptk

  real(r8) :: O2AquaDiffusvityP
  real(r8) :: RTVLWA  !root H2O vol for NH4
  real(r8) :: RTVLWB  !root H2O vol for NH4B
  real(r8) :: RGasTranspFlxPrev(idg_beg:idg_NH3)       !diagnosed increment gas flux
  real(r8) :: ROXYLX
  real(r8) :: RGas_DisolvSoil_flx(idg_beg:idg_end)     !gas dissolution into aqueous concentration in soil
  real(r8) :: RTCR1,RTCR2
  real(r8) :: RTCRA                                    !root conductance scalar for gas transport between atmosphere and root inner space [m]
  real(r8) :: RTARRX
  real(r8) :: RootCO2Prod_tscaled                      !root CO2 gas efflux due to root respiration at time step for gas flux calculations
  real(r8) :: RRADS
  real(r8) :: RSolUptkTransp(idg_beg:idg_end)          !volatile tracer uptake by transpiration
  real(r8) :: RootOxyUptakePerPlant
  real(r8) :: ROxySoil2UptkPerPlant                    !aqueous O2 uptake flux from soil due to diffusion and transpiration-aided-advection
  real(r8) :: ROxyRoot2UptkPerPlant                    !aqueous O2 diffusion uptake flux from inside roots
  real(r8) :: ROxyRoot2Uptk                            !oyxgen uptake by root, [gO/d2]
  real(r8) :: RootUptkSoiSolute(idg_beg:idg_end)       !root uptake of aqueous volatile from soil to inside roots
  real(r8) :: RDXAqueous(idg_beg:idg_end)              !maximum fluxes can be taken by roots
  real(r8) :: RDFAqueous(idg_beg:idg_end)
  real(r8) :: RUPOST   !total oyxgen uptake from soil by roots and other processes
  real(r8) :: RUPNTX  !total uptake of NH4 from soil into roots
  real(r8) :: Root_gas2sol_flx(idg_beg:idg_NH3)    !gas dissolution into aqueous phase of the volatile tracers in roots
  real(r8) :: trcg_air2root_flx_loc(idg_beg:idg_NH3)   !diffusion flux of gas from atmosphere to inside roots
  real(r8) :: THETW1,THETM
  real(r8) :: RootOxyDemandPerPlant                    !O2 demand per plant at time step of flux calculation
  real(r8) :: DisolvedGasVolume(idg_beg:idg_NH3),VLWatMicPMO,VLWatMicPMM,VLsoiAirPMM
  real(r8) :: VOLWSP,VLWatMicPMA,VLWatMicPMB,VOLWSA,VOLWSB
  real(r8) :: VOLWAqueous(idg_beg:idg_end)  !solubility scaled aqueous volume 
  real(r8) :: VOLPNH3,VOLPNH3B
  real(r8) :: trcg_rootml_beg(idg_beg:idg_NH3)
  real(r8) :: trcs_rootml_beg(idg_beg:idg_NH3)
  real(r8) :: dtrc_err(idg_beg:idg_NH3)
  real(r8) :: X,tcopy
  real(r8) :: ZH3PA,ZH3PB,ZH3GA,ZH3GB
  integer  :: idg
  
!     begin_execution
  associate(                                                              &
    RootStrutElms_pft         => plt_biom%RootStrutElms_pft              ,& !input  :plant root structural element mass, [g d-2]
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft                 ,& !input  :threshold zero for plang growth calculation, [-]
    PlantPopulation_pft       => plt_site%PlantPopulation_pft            ,& !input  :plant population, [d-2]
    CumSoilThickMidL_vr       => plt_site%CumSoilThickMidL_vr            ,& !input  :depth to middle of soil layer from surface of grid cell, [m]
    AtmGasc                   => plt_site%AtmGasc                        ,& !input  :atmospheric gas concentrations, [g m-3]
    ZEROS                     => plt_site%ZEROS                          ,& !input  :threshold zero for numerical stability,[-]
    ZERO                      => plt_site%ZERO                           ,& !input  :threshold zero for numerical stability, [-]
    VLWatMicPM_vr             => plt_site%VLWatMicPM_vr                  ,& !input  :soil micropore water content, [m3 d-2]
    VLsoiAirPM_vr             => plt_site%VLsoiAirPM_vr                  ,& !input  :soil air content, [m3 d-2]
    TortMicPM_vr              => plt_site%TortMicPM_vr                   ,& !input  :micropore soil tortuosity, [m3 m-3]
    FILMM_vr                  => plt_site%FILMM_vr                       ,& !input  :soil water film thickness, [m]
    RGasTranspFlxPrev_vr      => plt_bgcr%RGasTranspFlxPrev_vr           ,& !input  :net gaseous flux, [g d-2 h-1]
    RO2AquaSourcePrev_vr      => plt_bgcr%RO2AquaSourcePrev_vr           ,& !input  :net aqueous O2 flux, [g d-2 h-1]
    ZERO4Uptk_pft             => plt_rbgc%ZERO4Uptk_pft                  ,& !input  :threshold zero for uptake calculation, [-]
    RootRespPotent_pvr        => plt_rbgc%RootRespPotent_pvr             ,& !input  :root respiration unconstrained by O2, [g d-2 h-1]
    RootO2Dmnd4Resp_pvr       => plt_rbgc%RootO2Dmnd4Resp_pvr            ,& !input  :root O2 demand from respiration, [g d-2 h-1]
    RootCO2AutorX_pvr         => plt_rbgc%RootCO2AutorX_pvr              ,& !input  :root respiration from previous time step, [g d-2 h-1]
    TScal4Difsvity_vr         => plt_soilchem%TScal4Difsvity_vr          ,& !input  :temperature effect on diffusivity,[-]
    trcs_VLN_vr               => plt_soilchem%trcs_VLN_vr                ,& !input  :effective relative tracer volume, [-]
    GasDifc_vr                => plt_soilchem%GasDifc_vr                 ,& !input  :gaseous diffusivity, [m2 h-1]
    SoluteDifusvty_vr         => plt_soilchem%SoluteDifusvty_vr          ,& !input  :aqueous diffusivity, [m2 h-1]
    GasSolbility_vr           => plt_soilchem%GasSolbility_vr            ,& !input  :gas solubility, [m3 m-3]
    SoilWatAirDry_vr          => plt_soilchem%SoilWatAirDry_vr           ,& !input  :air-dry water content, [m3 m-3]
    VLSoilMicP_vr             => plt_soilchem%VLSoilMicP_vr              ,& !input  :total micropore volume in layer, [m3 d-2]
    FracAirFilledSoilPoreM_vr => plt_soilchem%FracAirFilledSoilPoreM_vr  ,& !input  :soil air-filled porosity, [m3 m-3]
    DiffusivitySolutEffM_vr   => plt_soilchem%DiffusivitySolutEffM_vr    ,& !input  :coefficient for dissolution - volatilization, [-]
    iPlantCalendar_brch       => plt_pheno%iPlantCalendar_brch           ,& !input  :plant growth stage, [-]
    RootPoreTortu4Gas_pft     => plt_morph%RootPoreTortu4Gas_pft         ,& !input  :power function of root porosity used to calculate root gaseous diffusivity, [-]
    Root1stRadius_pvr         => plt_morph%Root1stRadius_pvr             ,& !input  :root layer diameter primary axes, [m]
    Root2ndMeanLens_pvr       => plt_morph%Root2ndMeanLens_pvr           ,& !input  :root layer average length, [m]
    RootPoreVol_pvr           => plt_morph%RootPoreVol_pvr               ,& !input  :root layer volume air, [m2 d-2]
    RootLenPerPlant_pvr       => plt_morph%RootLenPerPlant_pvr           ,& !input  :root layer length per plant, [m p-1]
    Root2ndXNumL_pvr           => plt_morph%Root2ndXNumL_pvr               ,& !input  :root layer number axes, [d-2]
    Root2ndRadius_pvr         => plt_morph%Root2ndRadius_pvr             ,& !input  :root layer diameter secondary axes, [m]
    RootRaidus_rpft           => plt_morph%RootRaidus_rpft               ,& !input  :root internal radius, [m]
    RootVH2O_pvr              => plt_morph%RootVH2O_pvr                  ,& !input  :root layer volume water, [m2 d-2]
    RootPorosity_pft          => plt_morph%RootPorosity_pft              ,& !input  :root porosity, [m3 m-3]
    Root1stXNumL_pvr          => plt_morph%Root1stXNumL_pvr              ,& !input  :root layer number primary axes, [d-2]
    NGTopRootLayer_pft        => plt_morph%NGTopRootLayer_pft            ,& !input  :soil layer at planting depth, [-]
    MainBranchNum_pft         => plt_morph%MainBranchNum_pft             ,& !input  :number of main branch,[-]
    RootO2Uptk_pvr            => plt_rbgc%RootO2Uptk_pvr                 ,& !inoput :aqueous O2 flux from roots to root water, [g d-2 h-1]
    RAutoRootO2Limter_rpvr    => plt_rbgc%RAutoRootO2Limter_rpvr         ,& !inoput :O2 constraint to root respiration (0-1), [-]
    REcoUptkSoilO2M_vr           => plt_rbgc%REcoUptkSoilO2M_vr                ,& !inoput :total O2 sink, [g d-2 t-1]
    RCO2Emis2Root_pvr         => plt_rbgc%RCO2Emis2Root_pvr              ,& !inoput :aqueous CO2 flux from roots to root water, [g d-2 h-1]
    trcg_air2root_flx_pvr     => plt_rbgc%trcg_air2root_flx_pvr          ,& !inoput :gaseous tracer flux through roots, [g d-2 h-1]
    trcg_Root_gas2aqu_flx_vr  => plt_rbgc%trcg_Root_gas2aqu_flx_vr       ,& !inoput :dissolution (+ve) - volatilization (-ve) gas flux in roots, [g d-2 h-1]
    RootUptkSoiSol_pvr        => plt_rbgc%RootUptkSoiSol_pvr             ,& !inoput :aqueous CO2 flux from roots to soil water, [g d-2 h-1]
    trcg_rootml_pvr           => plt_rbgc%trcg_rootml_pvr                ,& !inoput :root gas content, [g d-2]
    trcs_rootml_pvr           => plt_rbgc%trcs_rootml_pvr                ,& !inoput :root aqueous content, [g d-2]
    RootCO2Ar2Soil_pvr        => plt_rbgc%RootCO2Ar2Soil_pvr             ,& !inoput :root respiration released to soil, [gC d-2 h-1]
    RootCO2Ar2RootX_pvr       => plt_rbgc%RootCO2Ar2RootX_pvr            ,& !inoput :root respiration released to root, [gC d-2 h-1]
    RootO2_TotSink_pvr           => plt_bgcr%RootO2_TotSink_pvr                ,& !output :root O2 sink for autotrophic respiraiton, [gC d-2 h-1]
    RootGasConductance_pvr    => plt_rbgc%RootGasConductance_pvr          & !output :Conductance for gas diffusion [m3 d-2 h-1]
  )
  
  call PrintInfo('beg '//subname)
  IF(RootRespPotent_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ).AND.RootVH2O_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ) &
    .AND.FOXYX.GT.ZERO4Uptk_pft(NZ))THEN
    RootCO2Ar2RootX_pvr(L,NZ)=RootCO2Ar2RootX_pvr(L,NZ)-RootCO2AutorX_pvr(N,L,NZ)
!
!     INITIALIZE VARIABLES USED IN ROOT GAS EXCHANGE
!     (CO2, O2, CH4, N2, N2O, NH3, H2, Ar)
!
!     RTVLWA,RTVLWB=root aqueous volume in non-band,band
!     dts_gas=time step of flux calculation (1/(NPH*NPT))
!     RGasTranspFlxPrev(idg_O2)=diagnosed net O2 gas flux at time step of flux calculation
!     RGasTranspFlxPrev(idg_CO2)=diagnosed net CO2 gas flux at time step of flux calculation
!     ROXYLX=diagnosed net O2 aqueous flux at time step of flux calculation
!

    DO idg=idg_beg,idg_NH3
      trcg_rootml_loc(idg) = AMAX1(ZERO4Groth_pft(NZ),trcg_rootml_pvr(idg,N,L,NZ))
      trcs_rootml_loc(idg) = AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_pvr(idg,N,L,NZ))
      trcg_rootml_beg(idg) = trcg_rootml_pvr(idg,N,L,NZ)
      trcs_rootml_beg(idg) = trcs_rootml_pvr(idg,N,L,NZ)
    ENDDO


    RTVLWA = RootVH2O_pvr(N,L,NZ)*trcs_VLN_vr(ids_NH4,L)
    RTVLWB = RootVH2O_pvr(N,L,NZ)*trcs_VLN_vr(ids_NH4B,L)

    DO idg=idg_beg,idg_NH3
      if(idg==idg_CO2 .or. idg==idg_O2)then
        RGasTranspFlxPrev(idg)  = RGasTranspFlxPrev_vr(idg,L)*dts_gas*FOXYX    
      else
        RGasTranspFlxPrev(idg)  = 0._r8
      endif
    ENDDO
    
    ROXYLX                = -RO2AquaSourcePrev_vr(L)*FOXYX*dts_gas   !>0 into dissolved phase
    RootOxyDemandPerPlant = RootO2Dmnd4Resp_pvr(N,L,NZ)*dts_gas/PlantPopulation_pft(NZ)
!
!     GASEOUS AND AQUEOUS DIFFUSIVITIES IN ROOT AND SOIL
!
    do idg=idg_beg,idg_NH3
      GasDifc_tscaled(idg) = GasDifc_vr(idg,L)*dts_gas*RootPoreTortu4Gas_pft(N,NZ)
      SolDifc_tscaled(idg) = SoluteDifusvty_vr(idg,L)*dts_gas*FOXYX
    enddo
    O2AquaDiffusvityP = SoluteDifusvty_vr(idg_O2,L)*dts_gas

    RGas_DisolvSoil_flx(idg_beg:idg_end) = 0.0_r8
    Root_gas2sol_flx(idg_beg:idg_NH3)    = 0.0_r8
!
!     ROOT CONDUCTANCE TO GAS TRANSFER
!
    IF(RootStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ) .AND. FracSoiLayByPrimRoot(L,NZ).GT.ZERO)THEN
      !primary roots conductance scalar
      RTCR1 = AMAX1(PlantPopulation_pft(NZ),Root1stXNumL_pvr(N,L,NZ))*PICON*Root1stRadius_pvr(N,L,NZ)**2/CumSoilThickMidL_vr(L)
      !secondary roots conductance scalar
      RTCR2 = (Root2ndXNumL_pvr(N,L,NZ)*PICON*Root2ndRadius_pvr(N,L,NZ)**2/Root2ndMeanLens_pvr(N,L,NZ))/FracSoiLayByPrimRoot(L,NZ)
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

    IF(N.EQ.ipltroot .AND. iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).GT.0 &
      .AND. RootLenPerPlant_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
      RTARRX = RootAreaDivRadius_vr(N,L)/RootRaidus_rpft(N,NZ)
      DIFOP  = O2AquaDiffusvityP*RTARRX
      DO idg  = idg_beg, idg_NH3
        DisolvedGasVolume(idg) = RootVH2O_pvr(N,L,NZ)*GasSolbility_vr(idg,L)
        DFAGas(idg)            = GasDifc_tscaled(idg)*RTCRA
      ENDDO
    ELSE
      RTARRX                             = 0.0_r8
      DIFOP                              = 0.0_r8
      DisolvedGasVolume(idg_beg:idg_NH3) = 0.0_r8
      DFAGas(idg_beg:idg_NH3)            = 0.0_r8
      DFAGas(idg_O2)                     = AMIN1(1.e-12_r8,RootPoreVol_pvr(N,L,NZ))
    ENDIF

    DO idg=idg_beg,idg_NH3
      RootGasConductance_pvr(idg,N,L,NZ) = AMIN1(DFAGas(idg),RootPoreVol_pvr(N,L,NZ))
    ENDDO

    DFGP                = AMIN1(1.0,XNPD*SQRT(RootPorosity_pft(N,NZ))*TScal4Difsvity_vr(L))
    !The root respiration below is from the previous time step
    RootCO2Prod_tscaled = -RootCO2AutorX_pvr(N,L,NZ)*dts_gas

!
!     SOLVE FOR GAS EXCHANGE IN SOIL AND ROOTS DURING ROOT UPTAKE
!     AT SMALLER TIME STEP NPH
!
    D99: DO M=1,NPH
!
!     AQUEOUS GAS DIFFUSIVITY THROUGH SOIL WATER TO ROOT
!
      VLWatMicPMO = VLWatMicPM_vr(M,L)*FOXYX
      VLWatMicPMM = VLWatMicPM_vr(M,L)*FracPRoot4Uptake(N,L,NZ)
      VLsoiAirPMM = VLsoiAirPM_vr(M,L)*FracPRoot4Uptake(N,L,NZ)
      VOLWSP      = RootVH2O_pvr(N,L,NZ)+VLWatMicPMM
      VLWatMicPMA = VLWatMicPMM*trcs_VLN_vr(ids_NH4,L)
      VLWatMicPMB = VLWatMicPMM*trcs_VLN_vr(ids_NH4B,L)
      VOLWSA      = RTVLWA+VLWatMicPMA
      VOLWSB      = RTVLWB+VLWatMicPMB
      THETW1      = AZMAX1(VLWatMicPM_vr(M,L)/VLSoilMicP_vr(L))

      IF(THETW1.GT.SoilWatAirDry_vr(L) .AND. FracPRoot4Uptake(N,L,NZ).GT.ZERO4Uptk_pft(NZ))THEN
        THETM  = TortMicPM_vr(M,L)*THETW1
        RRADS  = LOG((FILMM_vr(M,L)+FineRootRadius(N,L))/FineRootRadius(N,L))
        RTARRX = RootAreaDivRadius_vr(N,L)/RRADS

        do idg    = idg_beg, idg_NH3
          DifAqueVolatile(idg)=THETM*SolDifc_tscaled(idg)*RTARRX
        enddo
        DifAqueVolatile(idg_NH3)  = DifAqueVolatile(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
        DifAqueVolatile(idg_NH3B) = DifAqueVolatile(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)

        DO idg=idg_beg,idg_end
          VOLWAqueous(idg)=VLWatMicPMM*GasSolbility_vr(idg,L)
        ENDDO
        VOLWAqueous(idg_NH3)  = VOLWAqueous(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
        VOLWAqueous(idg_NH3B) = VOLWAqueous(idg_NH3B)*trcs_VLN_vr(ids_NH4B,L)

        VOLPNH3  = VLsoiAirPMM*trcs_VLN_vr(ids_NH4,L)
        VOLPNH3B = VLsoiAirPMM*trcs_VLN_vr(ids_NH4B,L)
!
!     MASS FLOW OF GAS FROM SOIL TO ROOT AT SHORTER TIME STEP NPT
!
!     ROXYLX=soil net O2 aqueous flux > 0
!     VLWatMicPMM=micropore water volume
!     RootVH2O_pvr,RootPoreVol_pvr=root aqueous,gaseous volume
!     RMF*=soil convective solute flux:COS=CO2,OXS=O2,CHS=CH4,
!     N2S=N2O,NHS=NH3 non-band,NHB=NH3 band,HGS=H2
!
        D90: DO MX=1,NPT
          tcopy=trc_solml_loc(idg_O2)
          call fixEXConsumpFlux(trc_solml_loc(idg_O2),ROXYLX)
          do idg=idg_beg,idg_end
            if(idg==idg_O2)then
              trcaqu_conc_soi_loc(idg_O2)=AMIN1(AtmGasc(idg_O2)*GasSolbility_vr(idg_O2,L),&
                AZMAX1(trc_solml_loc(idg_O2)/VLWatMicPMO))
            else
              trcaqu_conc_soi_loc(idg)=AZMAX1(trc_solml_loc(idg)/VLWatMicPMM)
            endif            
          enddo
          IF(RootPoreVol_pvr(N,L,NZ).GT.ZERO)THEN
            DO idg=idg_beg,idg_NH3
              trcg_gcon_loc(idg)=AZMAX1(trcg_rootml_loc(idg)/RootPoreVol_pvr(N,L,NZ))
            ENDDO
          ELSE
            trcg_gcon_loc(idg_beg:idg_NH3)=0.0_r8
          ENDIF

          do idg=idg_beg,idg_NH3
            trc_conc_root_loc(idg)=AZMAX1(trcs_rootml_loc(idg)/RootVH2O_pvr(N,L,NZ))
          enddo
          trc_conc_root_loc(idg_O2)=AMIN1(AtmGasc(idg_O2)*GasSolbility_vr(idg_O2,L),trc_conc_root_loc(idg_O2))
          DIFOX                    =DifAqueVolatile(idg_O2)+DIFOP
          do idg=idg_beg,idg_end
            RSolUptkTransp(idg)=dtPerPlantRootH2OUptake*trcaqu_conc_soi_loc(idg)
          enddo
          RSolUptkTransp(idg_NH3)  = RSolUptkTransp(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
          RSolUptkTransp(idg_NH3B) = RSolUptkTransp(idg_NH3B)*trcs_VLN_vr(ids_NH4B,L)
          !
          !     SOLUTION FOR MASS FLOW + DIFFUSION OF O2 IN AQUEOUS PHASES OF
          !     SOIL AND ROOT = ACTIVE UPTAKE OF O2 BY ROOT
          !     Oxygen is treated in its special way, assuming fluxes from soil and root are supporting uptake directory
          !     
          !     DifAqueVolatile(idg_O2)=O2 
          !     dtPerPlantRootH2OUptake=water uptake
          !     trcaqu_conc_soi_loc(idg_O2),trc_conc_root_loc(idg_O2)=soil,root aqueous O2 concentration
          !     RootOxyDemandPerPlant=O2 demand per plant
          !     RootOxyUptakePerPlant=root O2 uptake per plant
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
          !
          !     MASS FLOW + DIFFUSIVE EXCHANGE OF OTHER GASES
          !     BETWEEN ROOT AND SOIL, CONSTRAINED BY COMPETITION
          !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
          !
          RootUptkSoiSolute(idg_beg:idg_end)=0.0_r8          
          ROxySoil2Uptk             = ROxySoil2UptkPerPlant*PlantPopulation_pft(NZ)
          ROxyRoot2Uptk             = ROxyRoot2UptkPerPlant*PlantPopulation_pft(NZ)
          RootUptkSoiSolute(idg_O2) = ROxySoil2Uptk
          RDFAqueous(idg_CO2)       = RSolUptkTransp(idg_CO2)+DifAqueVolatile(idg_CO2)*(trcaqu_conc_soi_loc(idg_CO2)-trc_conc_root_loc(idg_CO2))
          RDXAqueous(idg_CO2)       = (RootVH2O_pvr(N,L,NZ)*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_CO2)) &
            -VLWatMicPMM*AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_loc(idg_CO2)))/VOLWSP

          IF(RDFAqueous(idg_CO2).GT.0.0_r8)THEN
            !CO2 from soil into roots          
            RootUptkSoiSolute(idg_CO2)=AMIN1(AZMAX1(RDXAqueous(idg_CO2)),RDFAqueous(idg_CO2)*PlantPopulation_pft(NZ))
          ELSE
            !CO2 from roots into soil  
            RootUptkSoiSolute(idg_CO2)=AMAX1(AZMIN1(RDXAqueous(idg_CO2)),RDFAqueous(idg_CO2)*PlantPopulation_pft(NZ))
          ENDIF

          IF(N.EQ.ipltroot)THEN
            !fluxes only involve roots
            DO idg=idg_beg,idg_NH3
              if(idg/=idg_CO2 .and. idg/=idg_NH3 .and. idg/=idg_O2)then
                !soil to root aqueous fluxes, transpiration + diffusion
                RDFAqueous(idg) = RSolUptkTransp(idg)+DifAqueVolatile(idg)*(trcaqu_conc_soi_loc(idg)-trc_conc_root_loc(idg))

                RDXAqueous(idg) = (RootVH2O_pvr(N,L,NZ)*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg)) &
                  -VLWatMicPMM*AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_loc(idg)))/VOLWSP

                IF(RDFAqueous(idg).GT.0.0_r8)THEN
                  !flux to roots
                  RootUptkSoiSolute(idg)=AMIN1(AZMAX1(RDXAqueous(idg)),RDFAqueous(idg)*PlantPopulation_pft(NZ))
                ELSE
                  !flux leaves roots
                  RootUptkSoiSolute(idg)=AMAX1(AZMIN1(RDXAqueous(idg)),RDFAqueous(idg)*PlantPopulation_pft(NZ))
                ENDIF
              ENDIF
            ENDDO

            RDFAqueous(idg_NH3)=RSolUptkTransp(idg_NH3)+DifAqueVolatile(idg_NH3)*(trcaqu_conc_soi_loc(idg_NH3)-trc_conc_root_loc(idg_NH3))

            IF(VOLWSA.GT.ZERO4Groth_pft(NZ))THEN
              ZH3PA              = trcs_rootml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
              RDXAqueous(idg_NH3) = (RTVLWA*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_NH3)) &
                -VLWatMicPMA*AMAX1(ZERO4Groth_pft(NZ),ZH3PA))/VOLWSA
            ELSE
              RDXAqueous(idg_NH3)=0.0_r8
            ENDIF


            IF(RDFAqueous(idg_NH3).GT.0.0_r8)THEN
              RootUptkSoiSolute(idg_NH3)=AMIN1(AZMAX1(RDXAqueous(idg_NH3)),RDFAqueous(idg_NH3)*PlantPopulation_pft(NZ))
            ELSE
              RootUptkSoiSolute(idg_NH3)=AMAX1(AZMIN1(RDXAqueous(idg_NH3)),RDFAqueous(idg_NH3)*PlantPopulation_pft(NZ))
            ENDIF

            RDFAqueous(idg_NH3B)=RSolUptkTransp(idg_NH3B)+DifAqueVolatile(idg_NH3B)*(trcaqu_conc_soi_loc(idg_NH3B)-trc_conc_root_loc(idg_NH3))
            IF(VOLWSB.GT.ZERO4Groth_pft(NZ))THEN
              ZH3PB               = trcs_rootml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)
              RDXAqueous(idg_NH3B) = (RTVLWB*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_NH3B)) &
                -VLWatMicPMB*AMAX1(ZERO4Groth_pft(NZ),ZH3PB))/VOLWSB
            ELSE
              RDXAqueous(idg_NH3B)=0.0_r8
            ENDIF

            IF(RDFAqueous(idg_NH3B).GT.0.0_r8)THEN
              RootUptkSoiSolute(idg_NH3B)=AMIN1(AZMAX1(RDXAqueous(idg_NH3B)),RDFAqueous(idg_NH3B)*PlantPopulation_pft(NZ))
            ELSE
              RootUptkSoiSolute(idg_NH3B)=AMAX1(AZMIN1(RDXAqueous(idg_NH3B)),RDFAqueous(idg_NH3B)*PlantPopulation_pft(NZ))
            ENDIF
          ENDIF

          !     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN SOIL
          !     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
          !     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENTS
          !     FROM 'WATSUB'
          !     RootUptkSoiSolute=root aqueous gas uptake
          !     ROXYLX=soil net O2 aqueous flux

          !gas disolution into soil
          IF(FracAirFilledSoilPoreM_vr(M,L).GT.AirFillPore_Min)THEN
            DiffusivitySolutEffP         = FracPRoot4Uptake(N,L,NZ)*DiffusivitySolutEffM_vr(M,L)
            RGas_DisolvSoil_flx(idg_CO2) = DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(idg_CO2))*VOLWAqueous(idg_CO2) &
              -(AMAX1(ZEROS,trc_solml_loc(idg_CO2))-RootUptkSoiSolute(idg_CO2))*VLsoiAirPMM) &
              /(VOLWAqueous(idg_CO2)+VLsoiAirPMM)

            RUPOST                      = ROxySoil2Uptk+ROXYLX
            RGas_DisolvSoil_flx(idg_O2) = DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(idg_O2))*VOLWAqueous(idg_O2) &
              -(AMAX1(ZEROS,trc_solml_loc(idg_O2))-RUPOST)*VLsoiAirPMM)/(VOLWAqueous(idg_O2)+VLsoiAirPMM)

            IF(N.EQ.ipltroot)THEN
              DO idg=idg_beg,idg_NH3-1
                if(idg/=idg_CO2 .and. idg/=idg_O2)then
                  RGas_DisolvSoil_flx(idg)=DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(idg))*VOLWAqueous(idg) &
                    -(AMAX1(ZEROS,trc_solml_loc(idg))-RootUptkSoiSolute(idg))*VLsoiAirPMM)/(VOLWAqueous(idg)+VLsoiAirPMM)
                endif
              ENDDO

              IF(VOLWAqueous(idg_NH3)+VOLPNH3.GT.ZERO4Groth_pft(NZ))THEN
                ZH3GA               = trc_gasml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
                RGas_DisolvSoil_flx(idg_NH3) = AMIN1(RootUptkSoiSolute(idg_NH3),AMAX1(-RootUptkSoiSolute(idg_NH3) &
                  ,DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),ZH3GA)*VOLWAqueous(idg_NH3) &
                  -(AMAX1(ZEROS,trc_solml_loc(idg_NH3))-RootUptkSoiSolute(idg_NH3))*VOLPNH3) &
                  /(VOLWAqueous(idg_NH3)+VOLPNH3)))
              ELSE
                RGas_DisolvSoil_flx(idg_NH3)=0.0_r8
              ENDIF

              IF(VOLWAqueous(idg_NH3B)+VOLPNH3B.GT.ZERO4Groth_pft(NZ))THEN
                ZH3GB                     = trc_gasml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)
                RGas_DisolvSoil_flx(idg_NH3B) = AMIN1(RootUptkSoiSolute(idg_NH3B),AMAX1(-RootUptkSoiSolute(idg_NH3B) &
                  ,DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),ZH3GB)*VOLWAqueous(idg_NH3B)   &
                  -(AMAX1(ZEROS,trc_solml_loc(idg_NH3B))-RootUptkSoiSolute(idg_NH3B))*VOLPNH3B)  &
                  /(VOLWAqueous(idg_NH3B)+VOLPNH3B)))
              ELSE
                RGas_DisolvSoil_flx(idg_NH3B)=0.0_r8
              ENDIF
            ELSE
              DO idg=idg_beg,idg_NH3
                if(idg/=idg_CO2 .and. idg/=idg_O2)RGas_DisolvSoil_flx(idg)=0.0_r8
              ENDDO
            ENDIF
          ELSE
            RGas_DisolvSoil_flx(idg_beg:idg_end)=0.0_r8
          ENDIF
          !
          !     UPDATE GASEOUS, AQUEOUS GAS CONTENTS AND CONCENTRATIONS
          !     FROM GASEOUS-AQUEOUS EXCHANGE, SOIL GAS TRANSFERS
          DO idg=idg_beg,idg_NH3
            trc_gasml_loc(idg)  = trc_gasml_loc(idg)-RGas_DisolvSoil_flx(idg)+RGasTranspFlxPrev(idg)
          ENDDO

          trc_gasml_loc(idg_NH3)  = trc_gasml_loc(idg_NH3)-RGas_DisolvSoil_flx(idg_NH3B)
          call fixEXConsumpFlux(trc_gasml_loc(idg_NH3),RGas_DisolvSoil_flx(idg_NH3B))

          !aqueous concentrations in soil
          DO idg=idg_beg,idg_end
            trc_solml_loc(idg)=trc_solml_loc(idg)+RGas_DisolvSoil_flx(idg)-RootUptkSoiSolute(idg)
          enddo
          !
          !     GAS TRANSFER THROUGH ROOTS
          !
          IF(N.EQ.ipltroot .AND. RootPoreVol_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
            RUPNTX=RootUptkSoiSolute(idg_NH3)+RootUptkSoiSolute(idg_NH3B)
            !
            !     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN ROOTS
            !     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
            !     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENT
            !
            trcs_maxRootml_loc(idg_CO2) = AZMAX1(trcs_rootml_loc(idg_CO2)+RootCO2Prod_tscaled)
            trcs_maxRootml_loc(idg_NH3) = AZMAX1(trcs_rootml_loc(idg_NH3)+RUPNTX)
            trcs_maxRootml_loc(idg_O2)  = AZMAX1(trcs_rootml_loc(idg_O2)-ROxyRoot2Uptk)
            DO idg=idg_beg,idg_NH3
              if(idg/=idg_CO2 .and. idg/=idg_NH3 .and. idg/=idg_O2)then
                trcs_maxRootml_loc(idg)=AZMAX1(trcs_rootml_loc(idg)+RootUptkSoiSolute(idg))
              endif
              
              Root_gas2sol_flx(idg)=AMAX1(-trcs_maxRootml_loc(idg),DFGP*(AMAX1(ZERO4Groth_pft(NZ),trcg_rootml_loc(idg)) &
                *DisolvedGasVolume(idg)-trcs_maxRootml_loc(idg)*RootPoreVol_pvr(N,L,NZ)) &
                /(DisolvedGasVolume(idg)+RootPoreVol_pvr(N,L,NZ)))
              !>0._r8 into root, <0._r8 into atmosphere, assuming specific rate 1/hr
              trcg_air2root_flx_loc(idg)=RootGasConductance_pvr(idg,N,L,NZ)*(AtmGasc(idg)-trcg_gcon_loc(idg))
            enddo
          ELSE
            Root_gas2sol_flx(idg_beg:idg_NH3)      = 0.0_r8
            trcg_air2root_flx_loc(idg_beg:idg_NH3) = 0.0_r8
          ENDIF
!
!     UPDATE ROOT AQUEOUS, GASEOUS GAS CONTENTS AND CONCENTRATIONS'
!     FOR ROOT AQUEOUS-GASEOUS, GASEOUS-ATMOSPHERE EXCHANGES
!          
          !releas autotrophic respiration CO2 into roots
          trcs_rootml_loc(idg_CO2) = trcs_rootml_loc(idg_CO2)+RootCO2Prod_tscaled

          DO idg=idg_beg,idg_NH3
            Root_gas2sol_flx(idg)=AMAX1(AMIN1(trcg_rootml_loc(idg),Root_gas2sol_flx(idg)),-trcs_rootml_loc(idg))

            trcg_rootml_loc(idg) = trcg_rootml_loc(idg)-Root_gas2sol_flx(idg)

            call fixEXConsumpFlux(trcg_rootml_loc(idg),trcg_air2root_flx_loc(idg),-1)

            if(idg==idg_O2)then
              !oxygen is consumed inside roots
              trcs_rootml_loc(idg) = trcs_rootml_loc(idg)+Root_gas2sol_flx(idg)
              call fixEXConsumpFlux(trcs_rootml_loc(idg),ROxyRoot2Uptk)   
            else
              !non-O2 gases are added to the root inside.
              trcs_rootml_loc(idg) = trcs_rootml_loc(idg)+Root_gas2sol_flx(idg)
              call fixEXConsumpFlux(trcs_rootml_loc(idg),RootUptkSoiSolute(idg),-1)
            endif
          ENDDO

          !NH3 taken up from banded soil is added to inside root concentration          
          call fixEXConsumpFlux(trcs_rootml_loc(idg_NH3),RootUptkSoiSolute(idg_NH3B),-1)

!
!     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE'
!
          DO idg=idg_beg,idg_end
            RootUptkSoiSol_pvr(idg,N,L,NZ)=RootUptkSoiSol_pvr(idg,N,L,NZ)+RootUptkSoiSolute(idg)
          enddo
!
!     ACCUMULATE ROOT-ATMOSPHERE GAS EXCHANGE TO HOURLY TIME SCALE'
!
          DO idg=idg_beg,idg_NH3
            trcg_Root_gas2aqu_flx_vr(idg,N,L,NZ) = trcg_Root_gas2aqu_flx_vr(idg,N,L,NZ)+Root_gas2sol_flx(idg)
            trcg_air2root_flx_pvr(idg,N,L,NZ)    = trcg_air2root_flx_pvr(idg,N,L,NZ)+trcg_air2root_flx_loc(idg)
          ENDDO
          !
          ! ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE'
          !
          ! RootO2Uptk_pvr=root O2 uptake from root
          ! REcoUptkSoilO2M_vr=total O2 uptake from soil by all microbial,root popns
          ! Root CO2 emission includes actual CO2 production from respiration and flux exchange with soil
          RCO2Emis2Root_pvr(N,L,NZ) = RCO2Emis2Root_pvr(N,L,NZ)+RootCO2Prod_tscaled+RootUptkSoiSolute(idg_CO2)
          RootO2Uptk_pvr(N,L,NZ)    = RootO2Uptk_pvr(N,L,NZ)+ROxyRoot2Uptk    !uptake from O2 inside roots for root respiration
          REcoUptkSoilO2M_vr(M,L)      = REcoUptkSoilO2M_vr(M,L)+ROxySoil2Uptk      !uptake from soil O2 for root respiration
        ENDDO D90        
      ENDIF

    ENDDO D99

    DO idg=idg_beg,idg_NH3
      trcg_rootml_pvr(idg,N,L,NZ) = trcg_rootml_loc(idg)
      trcs_rootml_pvr(idg,N,L,NZ) = trcs_rootml_loc(idg)

      dtrc_err(idg)=trcg_rootml_beg(idg)+trcs_rootml_beg(idg)-trcg_rootml_loc(idg)-trcs_rootml_loc(idg) &
        +trcg_air2root_flx_pvr(idg,N,L,NZ)
      if(idg==idg_O2)then
        dtrc_err(idg)=dtrc_err(idg)-RootO2Uptk_pvr(N,L,NZ)
      elseif(idg==idg_CO2)then
        dtrc_err(idg)=dtrc_err(idg)+RCO2Emis2Root_pvr(N,L,NZ)
      else
        dtrc_err(idg)=dtrc_err(idg)+RootUptkSoiSol_pvr(idg,N,L,NZ)
      endif      
    ENDDO

    !check mass conservation error
    !
    ! O2 CONSTRAINTS TO ROOT RESPIRATION DEPENDS UPON RATIO'
    ! O2 is taken from inside the root and directly from the aqueous soil O2.

    PopPlantO2Uptake        = RootO2Uptk_pvr(N,L,NZ)+RootUptkSoiSol_pvr(idg_O2,N,L,NZ)
    RootO2_TotSink_pvr(N,L,NZ) = PopPlantO2Uptake     !include O2 uptake from soil and from inside the roots
    !to be used in next iteration
    RAutoRootO2Limter_rpvr(N,L,NZ) = AMIN1(1.0_r8,AZMAX1(PopPlantO2Uptake/RootO2Dmnd4Resp_pvr(N,L,NZ)))
  ELSE
    RootCO2Ar2Soil_pvr(L,NZ) = RootCO2Ar2Soil_pvr(L,NZ)-RootCO2AutorX_pvr(N,L,NZ)
    PopPlantO2Uptake         = 0.0_r8
    RootO2_TotSink_pvr(N,L,NZ)  = 0._r8
    IF(L.GT.NGTopRootLayer_pft(NZ))THEN
      RAutoRootO2Limter_rpvr(N,L,NZ)=RAutoRootO2Limter_rpvr(N,L-1,NZ)
    ELSE
      RAutoRootO2Limter_rpvr(N,L,NZ)=1.0
    ENDIF
    
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine RootSoilGasExchange
  ![tail]
end module RootGasMod
