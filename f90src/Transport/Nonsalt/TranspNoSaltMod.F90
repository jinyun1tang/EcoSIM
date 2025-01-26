module TranspNoSaltMod
!!
! Description:
!
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils,    only: destroy
  USE MiniMathMod,   ONLY: AZMAX1, fixnegmass
  use EcoSIMCtrlMod, only: PrintInfo
  use SOMDataType
  use ChemTranspDataType
  use GridConsts
  use EcoSIMSolverPar
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use ClimForcDataType
  use FertilizerDataType
  use EcosimConst
  use SurfLitterDataType
  use SnowDataType
  use SurfSoilDataType
  use LandSurfDataType
  use RootDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use PlantDataRateType
  use GridDataType
  use BoundaryTranspMod
  use InsideTranspMod
  use TranspNoSaltDataMod
  use IrrigationDataType
  use EcoSiMParDataMod, only : micpar
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  real(r8), PARAMETER :: DPN4=5.7E-07

  public :: TranspNoSalt
  public :: InitTranspNoSalt,DestructTranspNoSalt
  contains

  subroutine InitTranspNoSalt()
  implicit none

  call InitTransfrData

  call InitInsTp
  end subroutine InitTranspNoSalt
!------------------------------------------------------------------------------------------
  subroutine DestructTranspNoSalt
  implicit none

  call DestructTransfrData

  call DestructInsTp
  end subroutine DestructTranspNoSalt

!------------------------------------------------------------------------------------------

  SUBROUTINE TranspNoSalt(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES 3-DIMENSIONAL FLUXES OF ALL SOIL
!     NON-SALT SOLUTES AND GASES
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='TranspNoSalt'
  integer :: MX,MM
  integer :: M
  real(r8) :: WaterFlow2Soil(3,JD,JV,JH)
  real(r8) :: trcsol_Irrig_flxM_vr(ids_beg:ids_end,JZ,JY,JX)  !solute tracer flux due to irrigation
  real(r8) :: RGasAtmDisol2LitrM(idg_beg:idg_NH3,JY,JX)
  real(r8) :: RGasAtmDisol2SoilM(idg_beg:idg_end,JY,JX)

!     execution begins here
  call PrintInfo('beg '//subname)
!
!     TIME STEPS FOR SOLUTE AND GAS FLUX CALCULATIONS
!

  WaterFlow2Soil(:,:,:,:)=0._r8

  call InitFluxStateVars(I,J,NHW,NHE,NVN,NVS,trcsol_Irrig_flxM_vr)

!
! TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
! NPH=NPX,NPT=NPY
! NPH=no. of cycles per hour for water, heat and solute flux calculations
! NPG=NPG=NPH*NPT,number of cycles per hour for gas flux calculations, 
! dt_GasCyc=1.0_r8/NPT, number of cycles for gas flux calculations for 1 iteration in NPH,
  MX=0
  DO  MM=1,NPG
    !compute the ordinal number of the intermediate forcing used for 
    !current iteration
    M=MIN(NPH,INT((MM-1)*dt_GasCyc)+1)

    call ModelTracerHydroFluxMM(I,J,M,MX,NHW, NHE, NVN, NVS,WaterFlow2Soil,RGasAtmDisol2LitrM,RGasAtmDisol2SoilM)

    !   BOUNDARY SOLUTE AND GAS FLUXES
    call XBoundaryFluxMM(M,MX,NHW,NHE,NVN,NVS)

    ! UPDATE STATE VARIABLES FROM TOTAL FLUXES CALCULATED ABOVE
    IF(MM.NE.NPG)call UpdateStateVarMM(I,J,M,MM,MX,NPG,NHW,NHE,NVN,NVS,trcsol_Irrig_flxM_vr,RGasAtmDisol2LitrM,RGasAtmDisol2SoilM)

    !MX records the M location in the moisture-heat calculation
    MX=M
  ENDDO
  call PrintInfo('end '//subname)
  
  END subroutine TranspNoSalt
!------------------------------------------------------------------------------------------

  subroutine InitFluxStateVars(I,J,NHW,NHE,NVN,NVS,trcsol_Irrig_flxM_vr)
  implicit none

  integer, intent(in) :: I, J,NHW,NHE,NVN,NVS
  real(r8), intent(out) :: trcsol_Irrig_flxM_vr(ids_beg:ids_end,JZ,JY,JX)   !solute flux due to irrigation
  character(len=*), parameter :: subname='InitFluxStateVars'
  integer :: NX,NY

  call PrintInfo('beg '//subname)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

!     GAS AND SOLUTE SINKS AND SOURCES IN SURFACE RESIDUE FROM MICROBIAL
!     TRANSFORMATIONS IN 'NITRO' + ROOT EXCHANGE IN 'EXTRACT'
!     + EQUILIBRIA REACTIONS IN 'SOLUTE' AT SUB-HOURLY TIME STEP
!
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
!     dts_gas=1/number of cycles h-1 for gas flux calculations
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!
      call SurfaceSinksandSources(NY,NX)
!
!     INITIALIZE STATE VARIABLES FOR USE IN GAS, SOLUTE FLUX CALCULATIONS
!     PH=pH
!
      call CopyLitterTracer4Transport(NY,NX)
!
!     INITIALIZE SURFACE SOLUTE FLUXES FROM ATMOSPHERE
!
      call ZeroSurfaceSoluteAtmosFlux(NY,NX)
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IN SNOWFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
!     ENTERED IN WEATHER AND IRRIGATION FILES
!
      call PrecipitationSoluteInput(I,J,NY,NX)

!     GAS AND SOLUTE FLUXES AT SUB-HOURLY FLUX TIME STEP
!     ENTERED IN SITE FILE
!
      call SubHourlyFluxesFromSiteFile(NY,NX)
!
!     SOLUTE FLUXES FROM WATSUB.F, NITRO.F, UPTAKE.F, SOLUTE.F
!
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!
      call ImportFluxFromOtherModules(I,NY,NX,trcsol_Irrig_flxM_vr)
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine InitFluxStateVars
!------------------------------------------------------------------------------------------

  subroutine UpdateStateVarMM(I,J,M,MM,MX,NPG,NHW,NHE,NVN,NVS,trcsol_Irrig_flxM_vr,RGasAtmDisol2LitrM,RGasAtmDisol2SoilM)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M, MM
  integer, intent(in) :: MX,NPG, NHW, NHE, NVN, NVS
  real(r8), intent(in) :: trcsol_Irrig_flxM_vr(ids_beg:ids_end,JZ,JY,JX)
  REAL(R8), intent(in) :: RGasAtmDisol2LitrM(idg_beg:idg_NH3,JY,JX)
  real(r8), intent(in) :: RGasAtmDisol2SoilM(idg_beg:idg_end,JY,JX)

  character(len=*),parameter :: subname='UpdateStateVarMM'
  integer :: NY,NX,K,L,idg,ids

  call PrintInfo('beg '//subname)
  !only execute before the last gas iteration
  D9695: DO NX=NHW,NHE
    D9690: DO NY=NVN,NVS
      IF(M.NE.MX)THEN
        ! STATE VARIABLES FOR SOLUTES IN SNOWPACK
        call UpdateSnowTracersM(NY,NX)

        call UpdateSurfTracerM(NY,NX,RGasAtmDisol2LitrM(:,NY,NX),RGasAtmDisol2SoilM(:,NY,NX))
      ENDIF
!
      call UpdateSoilTracersMM(I,J,M,MX,NY,NX,trcsol_Irrig_flxM_vr)
!
!     GAS EXCHANGE IN SURFACE LITTER
!
!     *S2=litter aqueous gas content
!     R*DFG=litter-atmosphere gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
      DO idg=idg_beg,idg_end
        trc_solml2_vr(idg,0,NY,NX)=trc_solml2_vr(idg,0,NY,NX)+RGas_Disol_FlxMM_vr(idg,0,NY,NX)
      ENDDO

      DO L=NU(NY,NX),NL(NY,NX)
        DO ids=ids_beg,ids_end
          trc_solml2_vr(ids,L,NY,NX)=fixnegmass(trc_solml2_vr(ids,L,NY,NX))
          trc_soHml2_vr(ids,L,NY,NX)=fixnegmass(trc_soHml2_vr(ids,L,NY,NX),trc_solml2_vr(ids,L,NY,NX))
        ENDDO  
      ENDDO  
    ENDDO D9690
  ENDDO D9695
  call PrintInfo('end '//subname)
  end subroutine UpdateStateVarMM
!------------------------------------------------------------------------------------------

  subroutine UpdateSnowTracersM(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L,idg,idn
!
!     *W2=solute content of snowpack
!     TQS*=net overland solute flux in snow
!     T*BLS=net solute flux in snowpack
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :N4=NH4,N3=NH3,NO=NO3,1P=HPO4,HP=H2PO4
!
  DO idg=idg_beg,idg_end-1
    trcg_solsml2_snvr(idg,1,NY,NX)=trcg_solsml2_snvr(idg,1,NY,NX)+trcg_SnowDrift(idg,NY,NX)
  ENDDO

  DO idn=ids_nut_beg,ids_nuts_end
    trcn_solsml2_snvr(idn,1,NY,NX)=trcn_solsml2_snvr(idn,1,NY,NX)+trcn_SnowDrift(idn,NY,NX)
  ENDDO

  DO  L=1,JS
    DO idg=idg_beg,idg_end-1
      trcg_solsml2_snvr(idg,L,NY,NX)=trcg_solsml2_snvr(idg,L,NY,NX)+trcg_TBLS_snvr(idg,L,NY,NX)
    ENDDO

    DO idn=ids_nut_beg,ids_nuts_end
      trcn_solsml2_snvr(idn,L,NY,NX)=trcn_solsml2_snvr(idn,L,NY,NX)+trcn_TBLS(idn,L,NY,NX)
    ENDDO
  ENDDO
  end subroutine UpdateSnowTracersM

!------------------------------------------------------------------------------------------

  subroutine UpdateSurfTracerM(NY,NX,RGasAtmDisol2LitrM,RGasAtmDisol2SoilM)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8),intent(in) :: RGasAtmDisol2LitrM(idg_beg:idg_NH3)
  real(r8),intent(in) :: RGasAtmDisol2SoilM(idg_beg:idg_end)

  integer :: K,idg,ids,idom
!     STATE VARIABLES FOR SOLUTES IN SURFACE RESIDUE AND IN
!     MICROPORES AND MACROPORES IN SOIL SURFACE LAYER FROM OVERLAND
!     FLOW AND SURFACE VOLATILIZATION-DISSOLUTION
!
!     *S2=litter solute content
!     R*DFR=gas exchange between atmosphere and surface litter water
!     R*DFS=gas exchange between atmosphere and soil surface water
!     R*FLS=convective + diffusive solute flux into litter,soil surface
!     R*FLW,R*FLB=convective + diffusive solute flux into litter from non-band,band
!     TQR*=net overland solute flux
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
  DO  K=1,jcplx
    DO idom=idom_beg,idom_end
      DOM_MicP2(idom,K,0,NY,NX)=DOM_MicP2(idom,K,0,NY,NX)+DOM_MicpTranspFlxM_3D(idom,K,3,0,NY,NX)
    ENDDO
  ENDDO

  !exclude NH3B
  DO idg=idg_beg,idg_NH3
    trc_solml2_vr(idg,0,NY,NX)=trc_solml2_vr(idg,0,NY,NX)+RGasAtmDisol2LitrM(idg)
  ENDDO

! band does not exist in litter layer
  DO ids=ids_beg,ids_end
    trc_solml2_vr(ids,0,NY,NX)=trc_solml2_vr(ids,0,NY,NX)+R3PoreSolFlx_3D(ids,3,0,NY,NX)
  ENDDO

! include NH3B
  DO idg=idg_beg,idg_end
    trc_solml2_vr(idg,NU(NY,NX),NY,NX)=trc_solml2_vr(idg,NU(NY,NX),NY,NX)+RGasAtmDisol2SoilM(idg)
  ENDDO

  D9680: DO K=1,jcplx
    DO idom=idom_beg,idom_end
      DOM_MicP2(idom,K,0,NY,NX)=DOM_MicP2(idom,K,0,NY,NX)+dom_TFloXSurRunoff(idom,K,NY,NX)
    ENDDO
  ENDDO D9680

!exclude NH3B
  DO idg=idg_beg,idg_end-1
    trc_solml2_vr(idg,0,NY,NX)=trc_solml2_vr(idg,0,NY,NX)+trcg_TFloXSurRunoff(idg,NY,NX)
  ENDDO

  DO ids=ids_nut_beg,ids_nuts_end
    trc_solml2_vr(ids,0,NY,NX)=trc_solml2_vr(ids,0,NY,NX)+trcn_TFloXSurRunoff_2D(ids,NY,NX)
  ENDDO

  end subroutine UpdateSurfTracerM

!------------------------------------------------------------------------------------------

  subroutine UpdateSoilTracersMM(I,J,M,MX,NY,NX,trcsol_Irrig_flxM_vr)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NY,NX,MX
  real(r8), intent(in) :: trcsol_Irrig_flxM_vr(ids_beg:ids_end,JZ,JY,JX)
  integer :: L,K,idg,ids,idom

!     STATE VARIABLES FOR SOLUTES IN MICROPORES AND
!     MACROPORES IN SOIL LAYERS FROM SUBSURFACE FLOW, MICROBIAL
!     AND ROOT EXCHANGE IN 'NITRO' AND 'UPTAKE', AND EQUILIBRIUM
!     REACTIONS IN 'SOLUTE'
!
!     *S2,*B2=micropore solute content in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     T*FLS=net convective + diffusive solute flux through micropores
!     T*FHS=net convective + diffusive solute flux through macropores
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     R*FLZ,R*FBZ=subsurface solute flux in non-band,band
!     R*BBL=bubble flux
!
  D9685: DO L=NU(NY,NX),NL(NY,NX)

    IF(M.NE.MX)THEN
      IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN        
        DO ids=ids_beg,ids_end
          trc_solml2_vr(ids,L,NY,NX) = trc_solml2_vr(ids,L,NY,NX)+trcsol_Irrig_flxM_vr(ids,L,NY,NX)
        ENDDO

        DO idg=idg_beg,idg_end        
          trc_solml2_vr(idg,L,NY,NX)=trc_solml2_vr(idg,L,NY,NX)+trcg_Ebu_vr(idg,L,NY,NX)
        ENDDO

        DO  K=1,jcplx
          DO idom=idom_beg,idom_end
            DOM_MicP2(idom,K,L,NY,NX) = DOM_MicP2(idom,K,L,NY,NX)+DOM_Transp2Micp_vr(idom,K,L,NY,NX)+DOM_XPoreTransp_flx(idom,K,L,NY,NX)
            DOM_MacP2(idom,K,L,NY,NX) = DOM_MacP2(idom,K,L,NY,NX)+DOM_Transp2Macp_flx(idom,K,L,NY,NX)-DOM_XPoreTransp_flx(idom,K,L,NY,NX)
          ENDDO
        ENDDO

        DO ids=ids_beg,ids_end
          trc_solml2_vr(ids,L,NY,NX) = trc_solml2_vr(ids,L,NY,NX)+TR3MicPoreSolFlx_vr(ids,L,NY,NX)+RMac2MicSolFlx_vr(ids,L,NY,NX)
          trc_solml2_vr(ids,L,NY,NX) = fixnegmass(trc_solml2_vr(ids,L,NY,NX))
          
          trc_soHml2_vr(ids,L,NY,NX) = trc_soHml2_vr(ids,L,NY,NX)+TR3MacPoreSolFlx_vr(ids,L,NY,NX)-RMac2MicSolFlx_vr(ids,L,NY,NX)
          trc_soHml2_vr(ids,L,NY,NX) = fixnegmass(trc_soHml2_vr(ids,L,NY,NX),trc_solml2_vr(ids,L,NY,NX))
 
        ENDDO
      ENDIF
    ENDIF
!
!     STATE VARIABLES FOR GASES IN SOIL LAYERS FROM SUBSURFACE FLOW,
!     MICROBIAL AND ROOT EXCHANGE IN 'NITRO' AND 'UPTAKE'
!
!     *G2,*S2=soil gas, solute content
!     R*DFG=water-air gas flux
!     T*FLG=net convective+diffusive gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
    IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DO idg=idg_beg,idg_end
        trc_solml2_vr(idg,L,NY,NX)=trc_solml2_vr(idg,L,NY,NX)+RGas_Disol_FlxMM_vr(idg,L,NY,NX)
      ENDDO

      DO idg=idg_beg,idg_end-1
        trc_gasml2_vr(idg,L,NY,NX)=trc_gasml2_vr(idg,L,NY,NX)+Gas_AdvDif_Flx_vr(idg,L,NY,NX)-RGas_Disol_FlxMM_vr(idg,L,NY,NX)
      ENDDO
      trc_gasml2_vr(idg_NH3,L,NY,NX)=trc_gasml2_vr(idg_NH3,L,NY,NX)-RGas_Disol_FlxMM_vr(idg_NH3B,L,NY,NX)
    ENDIF
  ENDDO D9685
  end subroutine UpdateSoilTracersMM

!------------------------------------------------------------------------------------------

  subroutine SurfaceSinksandSources(NY, NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K,idom
  !dts_gas=1/NPG with NPG=NPH*NPT
  !dts_HeatWatTP=1/NPT  
  !no O2 below 

  RBGCSinkGasMM_vr(idg_CO2,0,NY,NX) = trcs_RMicbUptake_vr(idg_CO2,0,NY,NX)*dts_gas
  RBGCSinkGasMM_vr(idg_CH4,0,NY,NX) = trcs_RMicbUptake_vr(idg_CH4,0,NY,NX)*dts_gas
  RBGCSinkGasMM_vr(idg_N2,0,NY,NX)  = (trcs_RMicbUptake_vr(idg_N2,0,NY,NX)+Micb_N2Fixation_vr(0,NY,NX))*dts_gas
  RBGCSinkGasMM_vr(idg_N2O,0,NY,NX) = trcs_RMicbUptake_vr(idg_N2O,0,NY,NX)*dts_gas
  RBGCSinkGasMM_vr(idg_H2,0,NY,NX)  = trcs_RMicbUptake_vr(idg_H2,0,NY,NX)*dts_gas
  RBGCSinkGasMM_vr(idg_NH3,0,NY,NX) = 0.0_r8

  RBGCSinkSoluteM_vr(idg_NH3,0,NY,NX)   = -TR_NH3_soil_vr(0,NY,NX)*dts_HeatWatTP
  RBGCSinkSoluteM_vr(ids_NH4,0,NY,NX)   = -(RNutMicbTransf_vr(ids_NH4,0,NY,NX)+trcn_RChem_soil_vr(ids_NH4,0,NY,NX))*dts_HeatWatTP
  RBGCSinkSoluteM_vr(ids_NO3,0,NY,NX)   = -(RNutMicbTransf_vr(ids_NO3,0,NY,NX)+trcn_RChem_soil_vr(ids_NO3,0,NY,NX))*dts_HeatWatTP
  RBGCSinkSoluteM_vr(ids_NO2,0,NY,NX)   = -(RNutMicbTransf_vr(ids_NO2,0,NY,NX)+trcn_RChem_soil_vr(ids_NO2,0,NY,NX))*dts_HeatWatTP
  RBGCSinkSoluteM_vr(ids_H2PO4,0,NY,NX) = -(RNutMicbTransf_vr(ids_H2PO4,0,NY,NX)+trcn_RChem_soil_vr(ids_H2PO4,0,NY,NX))*dts_HeatWatTP
  RBGCSinkSoluteM_vr(ids_H1PO4,0,NY,NX) = -(RNutMicbTransf_vr(ids_H1PO4,0,NY,NX)+trcn_RChem_soil_vr(ids_H1PO4,0,NY,NX))*dts_HeatWatTP

  end subroutine SurfaceSinksandSources
!------------------------------------------------------------------------------------------

  subroutine CopyLitterTracer4Transport(NY,NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K,idg,ids,idom

! exclude banded nutrient
  DO idg=idg_beg,idg_NH3
    trc_solml2_vr(idg,0,NY,NX)=trc_solml_vr(idg,0,NY,NX)
  ENDDO

! exclude banded nutrient
  DO ids=ids_nut_beg,ids_nuts_end
    trc_solml2_vr(ids,0,NY,NX)=trc_solml_vr(ids,0,NY,NX)
  ENDDO

  trc_solml2_vr(ids_H1PO4B,0,NY,NX) = trc_solml2_vr(ids_H1PO4,0,NY,NX)
  trc_solml2_vr(ids_H2PO4B,0,NY,NX) = trc_solml2_vr(ids_H2PO4,0,NY,NX)
  trc_solml2_vr(ids_NO3B,0,NY,NX)   = trc_solml2_vr(ids_NO3,0,NY,NX)
  trc_solml2_vr(ids_NO2B,0,NY,NX)   = trc_solml2_vr(ids_NO2,0,NY,NX)
  trc_solml2_vr(ids_NH4B,0,NY,NX)   = trc_solml2_vr(ids_NH4,0,NY,NX)
  trc_solml2_vr(idg_NH3B,0,NY,NX)   = trc_solml2_vr(idg_NH3,0,NY,NX)

  !reset DOM to value before the iteration
  D9979: DO K=1,jcplx
    DO idom=idom_beg,idom_end
      DOM_MicP2(idom,K,0,NY,NX)=DOM_vr(idom,K,0,NY,NX)
    ENDDO
  ENDDO D9979

  end subroutine CopyLitterTracer4Transport
!------------------------------------------------------------------------------------------

  subroutine ZeroSurfaceSoluteAtmosFlux(NY,NX)
  implicit none
  integer, intent(in) :: NY, NX
  integer :: K

  D8855: DO K=1,jcplx
    IF(K.LE.micpar%NumOfPlantLitrCmplxs)THEN
      DOM_MicpTransp_3D(idom_beg:idom_end,K,3,0,NY,NX)=0.0_r8
    ENDIF
    DOM_MicpTransp_3D(idom_beg:idom_end,K,3,NU(NY,NX),NY,NX)=0.0_r8
    DOM_3DMacp_Transp_flx(idom_beg:idom_end,K,3,NU(NY,NX),NY,NX)=0.0_r8
  ENDDO D8855

  end subroutine ZeroSurfaceSoluteAtmosFlux
!------------------------------------------------------------------------------------------

  subroutine PrecipitationSoluteInput(I,J,NY,NX)
  implicit none

  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  integer :: idg

  IF(SnoFalPrec_col(NY,NX).GT.0.0_r8 .OR. &   !snow fall
    (RainFalPrec_col(NY,NX).GT.0.0_r8 .AND. VLSnowHeatCapM_snvr(1,1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX)))THEN !rainfall with snowpack

    DO idg=idg_beg,idg_end-1
      trcVolatile_Xbndl_flx_snvr(idg,1,NY,NX)=Rain2SoilSurf_col(NY,NX)*trcVolatile_rain_conc(idg,NY,NX) &
        +Irrig2SoilSurf(NY,NX)*trcVolatile_irrig_conc(idg,NY,NX)
    ENDDO

    trcn_Xbndl_flx(ids_NH4,1,NY,NX)=(Rain2SoilSurf_col(NY,NX)*NH4_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw
    trcn_Xbndl_flx(ids_NO3,1,NY,NX)=(Rain2SoilSurf_col(NY,NX)*NO3_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw
    trcn_Xbndl_flx(ids_H1PO4,1,NY,NX)=(Rain2SoilSurf_col(NY,NX)*HPO4_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw
    trcn_Xbndl_flx(ids_H2PO4,1,NY,NX)=(Rain2SoilSurf_col(NY,NX)*H2PO4_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IF RAINFALL AND IRRIGATION IS ZERO IF SNOWPACK IS PRESENT
!
!     X*FLS,X*FLB=hourly solute flux to micropores in non-band,band
!
    trcs_TransptMicP_3D(ids_beg:ids_end,3,0,NY,NX)         = 0.0_r8
    trcs_TransptMicP_3D(ids_beg:ids_end,3,NU(NY,NX),NY,NX) = 0.0_r8
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IN RAINFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
!     ENTERED IN WEATHER AND IRRIGATION FILES
!
  ELSEIF((PrecAtm_col(NY,NX).GT.0.0_r8 .OR. IrrigSurface_col(NY,NX).GT.0.0_r8) & !Irrigation or rainfall
    .AND. VLSnowHeatCapM_snvr(1,1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX))THEN     !no snowpack
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IF SNOWFALL AND IRRIGATION IS ZERO AND SNOWPACK IS ABSENT
!
    trcVolatile_Xbndl_flx_snvr(idg_beg:idg_end-1,1,NY,NX) = 0.0_r8
    trcn_Xbndl_flx(ids_nut_beg:ids_nuts_end,1,NY,NX)      = 0.0_r8
    do idg=idg_beg,idg_NH3
      trcs_TransptMicP_3D(idg,3,0,NY,NX) = Rain2LitRSurf_col(NY,NX)*trcVolatile_rain_conc(idg,NY,NX) &
        +Irrig2LitRSurf(NY,NX)*trcVolatile_irrig_conc(idg,NY,NX)
      trcs_TransptMicP_3D(idg,3,NU(NY,NX),NY,NX)=Rain2SoilSurf_col(NY,NX)*trcVolatile_rain_conc(idg,NY,NX) &
        +Irrig2SoilSurf(NY,NX)*trcVolatile_irrig_conc(idg,NY,NX)
    ENDDO

    trcs_TransptMicP_3D(ids_NH4,3,0,NY,NX)=(Rain2LitRSurf_col(NY,NX)*NH4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw
    trcs_TransptMicP_3D(idg_NH3,3,0,NY,NX)=(Rain2LitRSurf_col(NY,NX)*trcVolatile_rain_conc(idg_NH3,NY,NX)+Irrig2LitRSurf(NY,NX)*trcVolatile_irrig_conc(idg_NH3,NY,NX))*natomw
    trcs_TransptMicP_3D(ids_NO3,3,0,NY,NX)=(Rain2LitRSurf_col(NY,NX)*NO3_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw
    trcs_TransptMicP_3D(ids_NO2,3,0,NY,NX)=0.0_r8
    trcs_TransptMicP_3D(ids_H1PO4,3,0,NY,NX)=(Rain2LitRSurf_col(NY,NX)*HPO4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw
    trcs_TransptMicP_3D(ids_H2PO4,3,0,NY,NX)=(Rain2LitRSurf_col(NY,NX)*H2PO4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw

    trcs_TransptMicP_3D(ids_NH4,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*NH4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw)*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
    trcs_TransptMicP_3D(idg_NH3,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*trcVolatile_rain_conc(idg_NH3,NY,NX) &
      +Irrig2SoilSurf(NY,NX)*trcVolatile_irrig_conc(idg_NH3,NY,NX))*natomw)*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
    trcs_TransptMicP_3D(ids_NO3,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*NO3_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw)*trcs_VLN_vr(ids_NO3,NU(NY,NX),NY,NX)

    trcs_TransptMicP_3D(ids_NO2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_TransptMicP_3D(ids_H1PO4,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*HPO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    trcs_TransptMicP_3D(ids_H2PO4,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*H2PO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    trcs_TransptMicP_3D(ids_NH4B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*NH4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw)*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)
    trcs_TransptMicP_3D(idg_NH3B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*trcVolatile_rain_conc(idg_NH3,NY,NX) &
      +Irrig2SoilSurf(NY,NX)*trcVolatile_irrig_conc(idg_NH3,NY,NX))*natomw)*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)
    trcs_TransptMicP_3D(ids_NO3B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*NO3_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw)*trcs_VLN_vr(ids_NO3B,NU(NY,NX),NY,NX)
    trcs_TransptMicP_3D(ids_NO2B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_TransptMicP_3D(ids_H1PO4B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*HPO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
    trcs_TransptMicP_3D(ids_H2PO4B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*H2PO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
!
!     NO SOLUTE FLUXES FROM ATMOSPHERE
!
  ELSE
    trcVolatile_Xbndl_flx_snvr(idg_beg:idg_end-1,1,NY,NX) = 0.0_r8
    trcn_Xbndl_flx(ids_nut_beg:ids_nuts_end,1,NY,NX)      = 0.0_r8

    trcs_TransptMicP_3D(idg_beg:idg_end-1,3,0,NY,NX)        = 0.0_r8
    trcs_TransptMicP_3D(ids_nut_beg:ids_nuts_end,3,0,NY,NX) = 0.0_r8
    trcs_TransptMicP_3D(ids_beg:ids_end,3,NU(NY,NX),NY,NX)  = 0.0_r8

  ENDIF

  trcs_TransptMacP_3D(ids_beg:ids_end,3,NU(NY,NX),NY,NX)=0.0_r8

  end subroutine PrecipitationSoluteInput
!------------------------------------------------------------------------------------------

  subroutine SubHourlyFluxesFromSiteFile(NY,NX)
  use EcoSiMParDataMod, only : micpar
  implicit none

  integer, intent(in) :: NY, NX
  INTEGER :: K,L,ids,idg,idn,idom

  DO  K=1,micpar%NumOfLitrCmplxs
    DO idom=idom_beg,idom_end
      DOM_Flo2LitrM(idom,K,NY,NX)=DOM_MicpTransp_3D(idom,K,3,0,NY,NX)*dts_HeatWatTP    
      DOM_Flo2TopSoilM(idom,K,NY,NX)=DOM_MicpTransp_3D(idom,K,3,NU(NY,NX),NY,NX)*dts_HeatWatTP
    ENDDO
  enddo

  !gas flux into snow
  DO idg=idg_beg,idg_end-1
    trcg_advW_snvr(idg_CO2,1,NY,NX)=trcVolatile_Xbndl_flx_snvr(idg_CO2,1,NY,NX)*dts_HeatWatTP
  ENDDO

  !nutrient flux
  DO idn=ids_nut_beg,ids_nuts_end
    trcn_RBLS(idn,1,NY,NX)=trcn_Xbndl_flx(idn,1,NY,NX)*dts_HeatWatTP
  ENDDO

  DO idg=idg_beg,idg_NH3
    trcg_RFL0(idg,NY,NX)=trcs_TransptMicP_3D(idg,3,0,NY,NX)*dts_HeatWatTP
  ENDDO

  do ids=ids_nut_beg,ids_nuts_end
    trcn_RFL0(ids,NY,NX)=trcs_TransptMicP_3D(ids,3,0,NY,NX)*dts_HeatWatTP
  enddo

  DO ids=ids_beg,ids_end
    trcs_RFL1(ids,NY,NX)=trcs_TransptMicP_3D(ids,3,NU(NY,NX),NY,NX)*dts_HeatWatTP
  ENDDO
!
!     GAS AND SOLUTE SINKS AND SOURCES IN SOIL LAYERS FROM MICROBIAL
!     TRANSFORMATIONS IN 'NITRO' + ROOT EXCHANGE IN 'EXTRACT'
!     + EQUILIBRIA REACTIONS IN 'SOLUTE' AT SUB-HOURLY TIME STEP
!
!     dt_GasCyc=1/number of cycles NPH-1 for gas flux calculations
!
  DO ids=ids_beg,ids_end
    SoluteDifusivitytscaled_vr(ids,0,NY,NX)=SoluteDifusvty_vr(ids,0,NY,NX)*dts_HeatWatTP
  ENDDO

  DO idom=idom_beg,idom_end
    DOMdiffusivity2_vr(idom,0,NY,NX)=DOMdiffusivity_vr(idom,0,NY,NX)*dts_HeatWatTP
  enddo
!
!     INITIAL SOLUTES IN SNOWPACK
!
  DO L=1,JS
    trcg_solsml2_snvr(idg_beg:idg_NH3,L,NY,NX)          = trcg_solsml_snvr(idg_beg:idg_NH3,L,NY,NX)
    trcn_solsml2_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX) = trcn_solsml_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX)
  enddo
  end subroutine SubHourlyFluxesFromSiteFile
!------------------------------------------------------------------------------------------

  subroutine ImportFluxFromOtherModules(I,NY,NX,trcsol_Irrig_flxM_vr)

  implicit none
  integer, intent(in) ::  I,NY, NX
  real(r8), intent(out) :: trcsol_Irrig_flxM_vr(ids_beg:ids_end,JZ,JY,JX)

  integer :: L,K,idg,ids,idom

  DO L=NU(NY,NX),NL(NY,NX)

    RBGCSinkGasMM_vr(idg_CO2,L,NY,NX) = (trcs_RMicbUptake_vr(idg_CO2,L,NY,NX) +trcs_plant_uptake_vr(idg_CO2,L,NY,NX)-TR_CO2_gchem_soil_vr(L,NY,NX))*dts_gas
    RBGCSinkGasMM_vr(idg_CH4,L,NY,NX) = (trcs_RMicbUptake_vr(idg_CH4,L,NY,NX) +trcs_plant_uptake_vr(idg_CH4,L,NY,NX))*dts_gas
    RBGCSinkGasMM_vr(idg_N2,L,NY,NX)  = (trcs_RMicbUptake_vr(idg_N2,L,NY,NX)  +trcs_plant_uptake_vr(idg_N2,L,NY,NX)+Micb_N2Fixation_vr(L,NY,NX))*dts_gas
    RBGCSinkGasMM_vr(idg_N2O,L,NY,NX) = (trcs_RMicbUptake_vr(idg_N2O,L,NY,NX) +trcs_plant_uptake_vr(idg_N2O,L,NY,NX))*dts_gas
    RBGCSinkGasMM_vr(idg_H2,L,NY,NX)  = (trcs_RMicbUptake_vr(idg_H2,L,NY,NX)  +trcs_plant_uptake_vr(idg_H2,L,NY,NX))*dts_gas
    RBGCSinkGasMM_vr(idg_NH3,L,NY,NX) = -TR_NH3_geochem_vr(L,NY,NX)*dts_gas

    RBGCSinkSoluteM_vr(ids_NH4,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NH4,L,NY,NX)-trcn_RChem_soil_vr(ids_NH4,L,NY,NX)+trcs_plant_uptake_vr(ids_NH4,L,NY,NX))*dts_HeatWatTP
    RBGCSinkSoluteM_vr(idg_NH3,L,NY,NX)   = (-TR_NH3_soil_vr(L,NY,NX)+trcs_plant_uptake_vr(idg_NH3,L,NY,NX))*dts_HeatWatTP
    RBGCSinkSoluteM_vr(ids_NO3,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NO3,L,NY,NX)-trcn_RChem_soil_vr(ids_NO3,L,NY,NX)+trcs_plant_uptake_vr(ids_NO3,L,NY,NX))*dts_HeatWatTP
    RBGCSinkSoluteM_vr(ids_NO2,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NO2,L,NY,NX)-trcn_RChem_soil_vr(ids_NO2,L,NY,NX))*dts_HeatWatTP
    RBGCSinkSoluteM_vr(ids_H2PO4,L,NY,NX) = (-RNutMicbTransf_vr(ids_H2PO4,L,NY,NX)-trcn_RChem_soil_vr(ids_H2PO4,L,NY,NX)+trcs_plant_uptake_vr(ids_H2PO4,L,NY,NX))*dts_HeatWatTP
    RBGCSinkSoluteM_vr(ids_H1PO4,L,NY,NX) = (-RNutMicbTransf_vr(ids_H1PO4,L,NY,NX)-trcn_RChem_soil_vr(ids_H1PO4,L,NY,NX)+trcs_plant_uptake_vr(ids_H1PO4,L,NY,NX))*dts_HeatWatTP

    RBGCSinkSoluteM_vr(ids_NH4B,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NH4B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_NH4B,L,NY,NX)+trcs_plant_uptake_vr(ids_NH4B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkSoluteM_vr(idg_NH3B,L,NY,NX)   = (-trcn_RChem_band_soil_vr(idg_NH3B,L,NY,NX)+trcs_plant_uptake_vr(idg_NH3B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkSoluteM_vr(ids_NO3B,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NO3B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_NO3B,L,NY,NX)+trcs_plant_uptake_vr(ids_NO3B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkSoluteM_vr(ids_NO2B,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NO2B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_NO2B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkSoluteM_vr(ids_H2PO4B,L,NY,NX) = (-RNutMicbTransf_vr(ids_H2PO4B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_H2PO4B,L,NY,NX)+trcs_plant_uptake_vr(ids_H2PO4B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkSoluteM_vr(ids_H1PO4B,L,NY,NX) = (-RNutMicbTransf_vr(ids_H1PO4B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_H1PO4B,L,NY,NX)+trcs_plant_uptake_vr(ids_H1PO4B,L,NY,NX))*dts_HeatWatTP
!
!     SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
    trcs_Irrig_vr(idg_CO2,L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*trcVolatile_irrig_conc(idg_CO2,NY,NX)
    trcs_Irrig_vr(idg_CH4,L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*trcVolatile_irrig_conc(idg_CH4,NY,NX)
    trcs_Irrig_vr(idg_O2,L,NY,NX)  = FWatIrrigate2MicP_vr(L,NY,NX)*trcVolatile_irrig_conc(idg_O2,NY,NX)
    trcs_Irrig_vr(idg_N2,L,NY,NX)  = FWatIrrigate2MicP_vr(L,NY,NX)*trcVolatile_irrig_conc(idg_N2,NY,NX)
    trcs_Irrig_vr(idg_N2O,L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*trcVolatile_irrig_conc(idg_N2O,NY,NX)
    trcs_Irrig_vr(idg_H2,L,NY,NX)  = 0.0_r8

    trcs_Irrig_vr(ids_NH4,L,NY,NX)    = FWatIrrigate2MicP_vr(L,NY,NX)*NH4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)*natomw
    trcs_Irrig_vr(idg_NH3,L,NY,NX)    = FWatIrrigate2MicP_vr(L,NY,NX)*trcVolatile_irrig_conc(idg_NH3,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_NO3,L,NY,NX)    = FWatIrrigate2MicP_vr(L,NY,NX)*NO3_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_NO3,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_H1PO4,L,NY,NX)  = FWatIrrigate2MicP_vr(L,NY,NX)*HPO4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)*patomw
    trcs_Irrig_vr(ids_H2PO4,L,NY,NX)  = FWatIrrigate2MicP_vr(L,NY,NX)*H2PO4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)*patomw
    trcs_Irrig_vr(ids_NH4B,L,NY,NX)   = FWatIrrigate2MicP_vr(L,NY,NX)*NH4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)*natomw
    trcs_Irrig_vr(idg_NH3B,L,NY,NX)   = FWatIrrigate2MicP_vr(L,NY,NX)*trcVolatile_irrig_conc(idg_NH3,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_NO3B,L,NY,NX)   = FWatIrrigate2MicP_vr(L,NY,NX)*NO3_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_NO3B,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_H1PO4B,L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*HPO4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)*patomw
    trcs_Irrig_vr(ids_H2PO4B,L,NY,NX) = FWatIrrigate2MicP_vr(L,NY,NX)*H2PO4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)*patomw
!
!     SUB-HOURLY SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     R*FLZ,R*FBZ=subsurface solute flux in non-band,band
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!
    DO ids=ids_beg,ids_end
      trcsol_Irrig_flxM_vr(ids,L,NY,NX)=trcs_Irrig_vr(ids,L,NY,NX)*dts_HeatWatTP
    ENDDO
!
!     GAS AND SOLUTE DIFFUSIVITIES AT SUB-HOURLY TIME STEP
!
    DO idom=idom_beg,idom_end
      DOMdiffusivity2_vr(idom,L,NY,NX)=DOMdiffusivity_vr(idom,L,NY,NX)*dts_HeatWatTP
    enddo

    DO ids=ids_beg,ids_end
      SoluteDifusivitytscaled_vr(ids,L,NY,NX)=SoluteDifusvty_vr(ids,L,NY,NX)*dts_HeatWatTP
    ENDDO

    DO idg=idg_beg,idg_end
      GasDifctScaled_vr(idg,L,NY,NX)=GasDifc_vr(idg,L,NY,NX)*dts_gas
    ENDDO
!
!     STATE VARIABLES FOR GASES AND SOLUTES USED IN 'TranspNoSalt'
!     TO STORE SUB-HOURLY CHANGES DURING FLUX CALCULATIONS
!     INCLUDING TRANSFORMATIONS FROM NITRO.F, UPTAKE.F AND SOLUTE.F
    ! exclude NH3B
    DO idg=idg_beg,idg_NH3
      trc_gasml2_vr(idg,L,NY,NX)=trc_gasml_vr(idg,L,NY,NX)
    ENDDO

    !reset DOM to value before the iteration
    !Concern: RDOMMicProd_vr includes microbial modification of DOM
    !it does not include exchange among different complexes (as one kind of microbial priming). 

    DO  K=1,jcplx
      DO idom=idom_beg,idom_end
        DOM_MicP2(idom,K,L,NY,NX)=DOM_vr(idom,K,L,NY,NX)
        DOM_MacP2(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)
      ENDDO
    enddo
    DO ids=ids_beg,ids_end
      trc_solml2_vr(ids,L,NY,NX) = fixnegmass(trc_solml_vr(ids,L,NY,NX))
      trc_soHml2_vr(ids,L,NY,NX) = fixnegmass(trc_soHml_vr(ids,L,NY,NX))

      if(trc_solml2_vr(ids,L,NY,NX)<0._r8)then
      write(*,*)'micv2 ',trcs_names(ids),L,trc_solml2_vr(ids,L,NY,NX),trc_solml_vr(ids,L,NY,NX)
      stop
      endif
      if(trc_soHml_vr(ids,L,NY,NX)<0._r8)then
      write(*,*)'macv2 ',trcs_names(ids),L,trc_soHml_vr(ids,L,NY,NX)
      endif
    ENDDO

  enddo
  end subroutine ImportFluxFromOtherModules
!

end module TranspNoSaltMod
