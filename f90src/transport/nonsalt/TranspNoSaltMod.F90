module TranspNoSaltMod
!!
! Description:
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  USE MiniMathMod, ONLY : AZMAX1,fixnegmass  
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
  integer :: MX,MM,M
  real(r8) :: WaterFlow2Soil(3,JD,JV,JH)
  real(r8) :: RFLZ_sol_vr(ids_beg:ids_end,JZ,JY,JX)

!     execution begins here
!
!     TIME STEPS FOR SOLUTE AND GAS FLUX CALCULATIONS
!

  WaterFlow2Soil(:,:,:,:)=0._r8
  call InitFluxandStateVariables(I,J,NHW,NHE,NVN,NVS,RFLZ_sol_vr)

!
! TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
! NPH=NPX,NPT=NPY
! NPH=no. of cycles per hour for water, heat and solute flux calculations
! NPG=number of cycles per hour for gas flux calculations, NPG=NPH*NPT
! dt_GasCyc=1/number of cycles for gas flux calculations for 1 iteration in NPH, =1.0_r8/NPT
  MX=0
  DO  MM=1,NPG
    !compute the ordinal number of the intermediate forcing used for 
    !current iteration
    M=MIN(NPH,INT((MM-1)*dt_GasCyc)+1)

    call ModelTracerHydroFlux(I,J,M,MX,NHW, NHE, NVN, NVS,WaterFlow2Soil)

!
!     BOUNDARY SOLUTE AND GAS FLUXES
!
    call BoundaryFlux(M,MX,NHW,NHE,NVN,NVS)
!
!     UPDATE STATE VARIABLES FROM TOTAL FLUXES CALCULATED ABOVE
!
    call UpdateStateVar(I,J,M,MM,MX,NPG,NHW,NHE,NVN,NVS,RFLZ_sol_vr)

    !MX records the M location of the last NPT iterations
    MX=M
  ENDDO
  RETURN
  END subroutine TranspNoSalt
!------------------------------------------------------------------------------------------

  subroutine InitFluxandStateVariables(I,J,NHW,NHE,NVN,NVS,RFLZ_sol_vr)
  implicit none

  integer, intent(in) :: I, J,NHW,NHE,NVN,NVS
  real(r8), intent(out) :: RFLZ_sol_vr(ids_beg:ids_end,JZ,JY,JX)

  integer :: NX,NY
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

!
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
!     R*SK2=total sink from nitro.f, uptake.f, solute.f
!     RCO2NetUptkMicb=net soil CO2 uptake from nitro.f
!     RCH4UptkAutor=net soil CH4 uptake from nitro.f
!     RN2NetUptkMicb=total soil N2 production from nitro.f
!     MicrbN2Fix=total N2 fixation from nitro.f
!     RN2O=net soil N2O uptake from nitro.f
!     RH2NetUptkMicb=net H2 uptake from nitro.f
!     RNH4MicbTransfSoil_vr=net change in NH4 from nitro.f
!     TR_NH4_soil,TR_NH3_soil_vr=NH4,NH3 dissolution from solute.f
!     RNO3MicbTransfSoil_vr=net change in NO3 from nitro.f
!     TRNO3=NO3 dissolution from solute.f
!     RNO2MicbTransfSoil_vr=net change in NO2 from nitro.f
!     TRNO2=NO2 dissolution from solute.f
!     RH2PO4MicbTransfSoil_vr=net change in H2PO4 from nitro.f
!     TR_H2PO4_soil=H2PO4 dissolution from solute.f
!     RH1PO4MicbTransfSoil_vr=net change in HPO4 from nitro.f
!     TR_H1PO4_soil=HPO4 dissolution from solute.f
!
      call SurfaceSinksandSources(NY,NX)
!
!     INITIALIZE STATE VARIABLES FOR USE IN GAS, SOLUTE FLUX CALCULATIONS
!
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 content

!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4
!     PH=pH
!
      call StateVarforGasandSolute(NY,NX)
!
!     INITIALIZE SURFACE SOLUTE FLUXES FROM ATMOSPHERE
!
!     X*FLS,X*FHS=hourly solute flux in macropores,micropores
!
      call SurfaceSolutefromAtms(NY,NX)
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IN SNOWFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
!     ENTERED IN WEATHER AND IRRIGATION FILES
!
!     SnoFalPrec,RainFalPrec=snow,rain
!     VLSnowHeatCapM,VLHeatCapSnowMin_col=current,minimum volumetric heat capacity of snowpack
!     X*BLS=hourly solute flux to snowpack
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation
!     C*R,C*Q=precipitation,irrigation solute concentrations
!     gas code: *CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*H2*=H2
!
      call HourlySoluteFluxes(I,J,NY,NX)
!     GAS AND SOLUTE FLUXES AT SUB-HOURLY FLUX TIME STEP
!     ENTERED IN SITE FILE
!
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     R*BLS,R*FL0,R*FL1,R*FL2=solute flux to snowpack,surface litter,soil surface non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     gas code: *CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*H2*=H2
!
      call SubHourlyFluxesFromSiteFile(NY,NX)
!
!     SOLUTE FLUXES FROM WATSUB.F, NITRO.F, UPTAKE.F, SOLUTE.F
!
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     PH=pH
!     FLWU,TUPWTR=total root water uptake from extract.f
!     R*SK2=total sink from nitro.f, uptake.f, solute.f
!     RCO2NetUptkMicb=net soil CO2 uptake from nitro.f
!     RCH4UptkAutor=net soil CH4 uptake from nitro.f
!     RN2NetUptkMicb=total soil N2 production from nitro.f
!     MicrbN2Fix=total N2 fixation from nitro.f
!     RN2O=net soil N2O uptake from nitro.f
!     RH2NetUptkMicb=net H2 uptake from nitro.f
!     RNH4MicbTransfSoil_vr=net change in NH4 from nitro.f
!     TR_NH4_soil,TR_NH3_soil_vr=NH4,NH3 dissolution from solute.f
!     RNO3MicbTransfSoil_vr=net change in NO3 from nitro.f
!     TRNO3=NO3 dissolution from solute.f
!     RNO2MicbTransfSoil_vr=net change in NO2 from nitro.f
!     TRNO2=NO2 dissolution from solute.f
!     RH2PO4MicbTransfSoil_vr=net change in H2PO4 from nitro.f
!     TR_H2PO4_soil=H2PO4 dissolution from solute.f
!     RH1PO4MicbTransfSoil_vr=net change in HPO4 from nitro.f
!     TR_H1PO4_soil=HPO4 dissolution from solute.f
!     TRootNH4Uptake_pft,TUPNHB=root NH4 uptake in non-band,band from extract.f
!     TRootNO3Uptake_pft,TUPNOB=root NO3 uptake in non-band,band from extract.f
!     TRootH2PO4Uptake_pft,TUPH2B=root H2PO4 uptake in non-band,band from extract.f
!     TRootHPO4Uptake_pft,TUPH1B=root HPO4 uptake in non-band,band from extract.f
!
      call ImportFluxFromOutsideModules(I,NY,NX,RFLZ_sol_vr)
    ENDDO
  ENDDO
  end subroutine InitFluxandStateVariables
!------------------------------------------------------------------------------------------

  subroutine UpdateStateVar(I,J,M,MM,MX,NPG,NHW,NHE,NVN,NVS,RFLZ_sol_vr)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) ::M, MM, MX,NPG, NHW, NHE, NVN, NVS
  real(r8), intent(in) :: RFLZ_sol_vr(ids_beg:ids_end,JZ,JY,JX)

  integer :: NY,NX,K,L,NGT,NTS

  IF(MM.NE.NPG)THEN
    D9695: DO NX=NHW,NHE
      D9690: DO NY=NVN,NVS
        IF(M.NE.MX)THEN
!
!     STATE VARIABLES FOR SOLUTES IN SNOWPACK
          call UpdateStateVarInSnowpack(NY,NX)
!
          call UpdateStateVarInResidue(NY,NX)

        ENDIF
!
        call UpdateSolutesInSoilLayers(I,J,M,MX,NY,NX,RFLZ_sol_vr)
!
!     GAS EXCHANGE IN SURFACE LITTER
!
!     *S2=litter aqueous gas content
!     R*DFG=litter-atmosphere gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
        DO NGT=idg_beg,idg_end
          trc_solml2_vr(NGT,0,NY,NX)=trc_solml2_vr(NGT,0,NY,NX)+RGasDSFlx_vr(NGT,0,NY,NX)
        ENDDO
      ENDDO D9690
    ENDDO D9695
  ENDIF

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO L=NU(NY,NX),NL(NY,NX)
        DO NTS=ids_beg,ids_end
          trc_solml2_vr(NTS,L,NY,NX)=fixnegmass(trc_solml2_vr(NTS,L,NY,NX))
          trc_soHml2_vr(NTS,L,NY,NX)=fixnegmass(trc_soHml2_vr(NTS,L,NY,NX),trc_solml2_vr(NTS,L,NY,NX))
        ENDDO  
      ENDDO  
    ENDDO
  ENDDO
  end subroutine UpdateStateVar
!------------------------------------------------------------------------------------------

  subroutine UpdateStateVarInSnowpack(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L,NTG,NTN
!
!     *W2=solute content of snowpack
!     TQS*=net overland solute flux in snow
!     T*BLS=net solute flux in snowpack
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :N4=NH4,N3=NH3,NO=NO3,1P=HPO4,HP=H2PO4
!
  DO NTG=idg_beg,idg_end-1
    trcg_solsml2_snvr(NTG,1,NY,NX)=trcg_solsml2_snvr(NTG,1,NY,NX)+trcg_SnowDrift(NTG,NY,NX)
  ENDDO

  DO NTN=ids_nut_beg,ids_nuts_end
    trcn_solsml2(NTN,1,NY,NX)=trcn_solsml2(NTN,1,NY,NX)+trcn_SnowDrift(NTN,NY,NX)
  ENDDO

  DO  L=1,JS
    DO NTG=idg_beg,idg_end-1
      trcg_solsml2_snvr(NTG,L,NY,NX)=trcg_solsml2_snvr(NTG,L,NY,NX)+trcg_TBLS_snvr(NTG,L,NY,NX)
    ENDDO

    DO NTN=ids_nut_beg,ids_nuts_end
      trcn_solsml2(NTN,L,NY,NX)=trcn_solsml2(NTN,L,NY,NX)+trcn_TBLS(NTN,L,NY,NX)
    ENDDO
  ENDDO
  end subroutine UpdateStateVarInSnowpack

!------------------------------------------------------------------------------------------

  subroutine UpdateStateVarInResidue(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer :: K,NTG,NTS,idom
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
  DO NTG=idg_beg,idg_end-1
    trc_solml2_vr(NTG,0,NY,NX)=trc_solml2_vr(NTG,0,NY,NX)+RDFR_gas(NTG,NY,NX)
  ENDDO

! band does not exist in litter layer
  DO NTS=ids_beg,ids_end
    trc_solml2_vr(NTS,0,NY,NX)=trc_solml2_vr(NTS,0,NY,NX)+R3PoreSolFlx_3D(NTS,3,0,NY,NX)
  ENDDO

! include NH3B
  DO NTG=idg_beg,idg_end
    trc_solml2_vr(NTG,NU(NY,NX),NY,NX)=trc_solml2_vr(NTG,NU(NY,NX),NY,NX)+RGasSSVol(NTG,NY,NX)
  ENDDO

  D9680: DO K=1,jcplx
    DO idom=idom_beg,idom_end
      DOM_MicP2(idom,K,0,NY,NX)=DOM_MicP2(idom,K,0,NY,NX)+dom_TFloXSurRunoff(idom,K,NY,NX)
    ENDDO
  ENDDO D9680

!exclude NH3B
  DO NTG=idg_beg,idg_end-1
    trc_solml2_vr(NTG,0,NY,NX)=trc_solml2_vr(NTG,0,NY,NX)+trcg_TFloXSurRunoff(NTG,NY,NX)
  ENDDO

  DO NTS=ids_nut_beg,ids_nuts_end
    trc_solml2_vr(NTS,0,NY,NX)=trc_solml2_vr(NTS,0,NY,NX)+trcn_TFloXSurRunoff_2D(NTS,NY,NX)
  ENDDO

  end subroutine UpdateStateVarInResidue

!------------------------------------------------------------------------------------------

  subroutine UpdateSolutesInSoilLayers(I,J,M,MX,NY,NX,RFLZ_sol_vr)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NY,NX,MX
  real(r8), intent(in) :: RFLZ_sol_vr(ids_beg:ids_end,JZ,JY,JX)
  real(r8) :: dval
  integer :: L,K,NTG,NTS,idom

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
        
        DO NTS=ids_beg,ids_end
          dval=trc_solml2_vr(NTS,L,NY,NX)
          trc_solml2_vr(NTS,L,NY,NX)=trc_solml2_vr(NTS,L,NY,NX)+RFLZ_sol_vr(NTS,L,NY,NX)          
        ENDDO

        DO NTG=idg_beg,idg_end        
          dval=trc_solml2_vr(NTG,L,NY,NX)
          trc_solml2_vr(NTG,L,NY,NX)=trc_solml2_vr(NTG,L,NY,NX)+trcg_Ebu_vr(NTG,L,NY,NX)
        ENDDO

        DO  K=1,jcplx
          DO idom=idom_beg,idom_end
            DOM_MicP2(idom,K,L,NY,NX)=DOM_MicP2(idom,K,L,NY,NX)+ &
              DOM_Transp2Micp_vr(idom,K,L,NY,NX)+DOM_XPoreTransp_flx(idom,K,L,NY,NX)
            DOM_MacP2(idom,K,L,NY,NX)=DOM_MacP2(idom,K,L,NY,NX)+ &
              DOM_Transp2Macp_flx(idom,K,L,NY,NX)-DOM_XPoreTransp_flx(idom,K,L,NY,NX)
          ENDDO
        ENDDO

        DO NTS=ids_beg,ids_end
          dval=trc_solml2_vr(NTS,L,NY,NX)
          trc_solml2_vr(NTS,L,NY,NX)=trc_solml2_vr(NTS,L,NY,NX) &
            +TR3MicPoreSolFlx_vr(NTS,L,NY,NX)+RMac2MicSolFlx_vr(NTS,L,NY,NX)
          trc_solml2_vr(NTS,L,NY,NX)=fixnegmass(trc_solml2_vr(NTS,L,NY,NX))
          
          trc_soHml2_vr(NTS,L,NY,NX)=trc_soHml2_vr(NTS,L,NY,NX) &
            +TR3MacPoreSolFlx_vr(NTS,L,NY,NX)-RMac2MicSolFlx_vr(NTS,L,NY,NX)

          trc_soHml2_vr(NTS,L,NY,NX)=fixnegmass(trc_soHml2_vr(NTS,L,NY,NX),trc_solml2_vr(NTS,L,NY,NX))
 
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
      DO NTG=idg_beg,idg_end
        trc_solml2_vr(NTG,L,NY,NX)=trc_solml2_vr(NTG,L,NY,NX)+RGasDSFlx_vr(NTG,L,NY,NX)
      ENDDO

      DO NTG=idg_beg,idg_end-1
        trc_gasml2_vr(NTG,L,NY,NX)=trc_gasml2_vr(NTG,L,NY,NX) &
          +Gas_AdvDif_Flx_vr(NTG,L,NY,NX)-RGasDSFlx_vr(NTG,L,NY,NX)
      ENDDO
      trc_gasml2_vr(idg_NH3,L,NY,NX)=trc_gasml2_vr(idg_NH3,L,NY,NX)-RGasDSFlx_vr(idg_NH3B,L,NY,NX)
    ENDIF
  ENDDO D9685
  end subroutine UpdateSolutesInSoilLayers

!------------------------------------------------------------------------------------------

  subroutine SurfaceSinksandSources(NY, NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K,idom

  RBGCSinkG_vr(idg_CO2,0,NY,NX)=trcg_RMicbTransf_vr(idg_CO2,0,NY,NX)*dts_gas
  RBGCSinkG_vr(idg_CH4,0,NY,NX)=trcg_RMicbTransf_vr(idg_CH4,0,NY,NX)*dts_gas
  RBGCSinkG_vr(idg_N2,0,NY,NX)=(trcg_RMicbTransf_vr(idg_N2,0,NY,NX)+Micb_N2Fixation_vr(0,NY,NX))*dts_gas
  RBGCSinkG_vr(idg_N2O,0,NY,NX)=trcg_RMicbTransf_vr(idg_N2O,0,NY,NX)*dts_gas
  RBGCSinkG_vr(idg_H2,0,NY,NX)=trcg_RMicbTransf_vr(idg_H2,0,NY,NX)*dts_gas
  RBGCSinkG_vr(idg_NH3,0,NY,NX)=0.0_r8

!  DO  K=1,jcplx
!    DO idom=idom_beg,idom_end
!      RDOM_CumEcoProd_vr(idom,K,0,NY,NX)=-REcoDOMProd_vr(idom,K,0,NY,NX)*dts_HeatWatTP
!    ENDDO
!  ENDDO
  RBGCSinkS_vr(ids_NH4,0,NY,NX)   = (-RNutMicbTransf_vr(ids_NH4,0,NY,NX)-trcn_RChem_soil_vr(ids_NH4,0,NY,NX))*dts_HeatWatTP
  RBGCSinkS_vr(idg_NH3,0,NY,NX)   = -TR_NH3_soil_vr(0,NY,NX)*dts_HeatWatTP
  RBGCSinkS_vr(ids_NO3,0,NY,NX)   = (-RNutMicbTransf_vr(ids_NO3,0,NY,NX)-trcn_RChem_soil_vr(ids_NO3,0,NY,NX))*dts_HeatWatTP
  RBGCSinkS_vr(ids_NO2,0,NY,NX)   = (-RNutMicbTransf_vr(ids_NO2,0,NY,NX)-trcn_RChem_soil_vr(ids_NO2,0,NY,NX))*dts_HeatWatTP
  RBGCSinkS_vr(ids_H2PO4,0,NY,NX) = (-RNutMicbTransf_vr(ids_H2PO4,0,NY,NX)-trcn_RChem_soil_vr(ids_H2PO4,0,NY,NX))*dts_HeatWatTP
  RBGCSinkS_vr(ids_H1PO4,0,NY,NX) = (-RNutMicbTransf_vr(ids_H1PO4,0,NY,NX)-trcn_RChem_soil_vr(ids_H1PO4,0,NY,NX))*dts_HeatWatTP

  end subroutine SurfaceSinksandSources
!------------------------------------------------------------------------------------------

  subroutine StateVarforGasandSolute(NY,NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K,NTG,NTS,idom

! exclude banded nutrient
  DO NTG=idg_beg,idg_end-1
    trc_solml2_vr(NTG,0,NY,NX)=trc_solml_vr(NTG,0,NY,NX)
  ENDDO
  !reset DOM to value before the iteration
  D9979: DO K=1,jcplx
    DO idom=idom_beg,idom_end
!      DOM_MicP2(idom,K,0,NY,NX)=AZMAX1(DOM_vr(idom,K,0,NY,NX)-RDOMMicProd_vr(idom,K,0,NY,NX))      
      DOM_MicP2(idom,K,0,NY,NX)=DOM_vr(idom,K,0,NY,NX)
    ENDDO
  ENDDO D9979

! exclude banded nutrient
  DO NTS=ids_nut_beg,ids_nuts_end
    trc_solml2_vr(NTS,0,NY,NX)=trc_solml_vr(NTS,0,NY,NX)
  ENDDO

  trc_solml2_vr(ids_H1PO4B,0,NY,NX)=trc_solml2_vr(ids_H1PO4,0,NY,NX)
  trc_solml2_vr(ids_H2PO4B,0,NY,NX)=trc_solml2_vr(ids_H2PO4,0,NY,NX)
  trc_solml2_vr(ids_NO3B,0,NY,NX)=trc_solml2_vr(ids_NO3,0,NY,NX)
  trc_solml2_vr(ids_NO2B,0,NY,NX)=trc_solml2_vr(ids_NO2,0,NY,NX)
  trc_solml2_vr(ids_NH4B,0,NY,NX)=trc_solml2_vr(ids_NH4,0,NY,NX)
  trc_solml2_vr(idg_NH3B,0,NY,NX)=trc_solml2_vr(idg_NH3,0,NY,NX)

  end subroutine StateVarforGasandSolute
!------------------------------------------------------------------------------------------

  subroutine SurfaceSolutefromAtms(NY,NX)
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
  end subroutine SurfaceSolutefromAtms
!------------------------------------------------------------------------------------------

  subroutine HourlySoluteFluxes(I,J,NY,NX)
  implicit none

  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX

  IF(SnoFalPrec_col(NY,NX).GT.0.0_r8 .OR. (RainFalPrec(NY,NX).GT.0.0_r8 .AND. VLSnowHeatCapM_snvr(1,1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX)))THEN
    trcg_Xbndl_flx(idg_CO2,1,NY,NX)=Rain2SoilSurf_col(NY,NX)*CO2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*CO2_irrig_conc(NY,NX)
    trcg_Xbndl_flx(idg_CH4,1,NY,NX)=Rain2SoilSurf_col(NY,NX)*CH4_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*CH4_irrig_conc(NY,NX)
    trcg_Xbndl_flx(idg_O2,1,NY,NX)=Rain2SoilSurf_col(NY,NX)*O2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*O2_irrig_conc(NY,NX)
    trcg_Xbndl_flx(idg_N2,1,NY,NX)=Rain2SoilSurf_col(NY,NX)*N2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*N2_irrig_conc(NY,NX)
    trcg_Xbndl_flx(idg_N2O,1,NY,NX)=Rain2SoilSurf_col(NY,NX)*N2O_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*N2O_irrig_conc(NY,NX)
    trcg_Xbndl_flx(idg_NH3,1,NY,NX)=(Rain2SoilSurf_col(NY,NX)*NH3_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*NH3_irrig_conc(I,NY,NX))*natomw

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
    trcs_Transp2MicP_3D(ids_beg:ids_end,3,0,NY,NX)=0.0_r8
    trcs_Transp2MicP_3D(ids_beg:ids_end,3,NU(NY,NX),NY,NX)=0.0_r8
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IN RAINFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
!     ENTERED IN WEATHER AND IRRIGATION FILES
!
  ELSEIF((PrecAtm_col(NY,NX).GT.0.0_r8.OR.IrrigSurface_col(NY,NX).GT.0.0_r8) &
    .AND.VLSnowHeatCapM_snvr(1,1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX))THEN
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IF SNOWFALL AND IRRIGATION IS ZERO AND SNOWPACK IS ABSENT
!
!     PrecAtm_col,PRECI=snow+rain,irrigation
!     X*BLS=hourly solute flux to snowpack
!     X*FLS,X*FLB=hourly solute flux to surface litter,soil surface micropore non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     Rain2LitRSurf_col,Irrig2LitRSurf=water flux to surface litter from rain,irrigation
!     FLQGQ,FLQGI=water flux to soil surface from rain,irrigation
!     C*R,C*Q=precipitation,irrigation solute concentrations
!     gas code: *CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*H2*=H2
!
    trcg_Xbndl_flx(idg_beg:idg_end-1,1,NY,NX)=0.0_r8

    trcn_Xbndl_flx(ids_nut_beg:ids_nuts_end,1,NY,NX)=0.0_r8

    trcs_Transp2MicP_3D(idg_CO2,3,0,NY,NX)=Rain2LitRSurf_col(NY,NX)*CO2_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*CO2_irrig_conc(NY,NX)
    trcs_Transp2MicP_3D(idg_CH4,3,0,NY,NX)=Rain2LitRSurf_col(NY,NX)*CH4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*CH4_irrig_conc(NY,NX)
    trcs_Transp2MicP_3D(idg_O2,3,0,NY,NX)=Rain2LitRSurf_col(NY,NX)*O2_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*O2_irrig_conc(NY,NX)
    trcs_Transp2MicP_3D(idg_N2,3,0,NY,NX)=Rain2LitRSurf_col(NY,NX)*N2_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*N2_irrig_conc(NY,NX)
    trcs_Transp2MicP_3D(idg_N2O,3,0,NY,NX)=Rain2LitRSurf_col(NY,NX)*N2O_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*N2O_irrig_conc(NY,NX)
    trcs_Transp2MicP_3D(idg_H2,3,0,NY,NX)=0.0_r8

!    write(115,*)I+J/24.,'tp1',trcs_Transp2MicP_3D(ids_NO3,3,0,NY,NX)
    trcs_Transp2MicP_3D(ids_NH4,3,0,NY,NX)=(Rain2LitRSurf_col(NY,NX)*NH4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw
    trcs_Transp2MicP_3D(idg_NH3,3,0,NY,NX)=(Rain2LitRSurf_col(NY,NX)*NH3_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*NH3_irrig_conc(I,NY,NX))*natomw
    trcs_Transp2MicP_3D(ids_NO3,3,0,NY,NX)=(Rain2LitRSurf_col(NY,NX)*NO3_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw
    trcs_Transp2MicP_3D(ids_NO2,3,0,NY,NX)=0.0_r8
    trcs_Transp2MicP_3D(ids_H1PO4,3,0,NY,NX)=(Rain2LitRSurf_col(NY,NX)*HPO4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw
    trcs_Transp2MicP_3D(ids_H2PO4,3,0,NY,NX)=(Rain2LitRSurf_col(NY,NX)*H2PO4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw
!    write(115,*)I+J/24.,'tp2',trcs_Transp2MicP_3D(ids_NO3,3,0,NY,NX)

    trcs_Transp2MicP_3D(idg_CO2,3,NU(NY,NX),NY,NX)=Rain2SoilSurf_col(NY,NX)*CO2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*CO2_irrig_conc(NY,NX)
    trcs_Transp2MicP_3D(idg_CH4,3,NU(NY,NX),NY,NX)=Rain2SoilSurf_col(NY,NX)*CH4_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*CH4_irrig_conc(NY,NX)
    trcs_Transp2MicP_3D(idg_O2,3,NU(NY,NX),NY,NX)=Rain2SoilSurf_col(NY,NX)*O2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*O2_irrig_conc(NY,NX)
    trcs_Transp2MicP_3D(idg_N2,3,NU(NY,NX),NY,NX)=Rain2SoilSurf_col(NY,NX)*N2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*N2_irrig_conc(NY,NX)
    trcs_Transp2MicP_3D(idg_N2O,3,NU(NY,NX),NY,NX)=Rain2SoilSurf_col(NY,NX)*N2O_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*N2O_irrig_conc(NY,NX)
    trcs_Transp2MicP_3D(idg_H2,3,NU(NY,NX),NY,NX)=0.0_r8

    trcs_Transp2MicP_3D(ids_NH4,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*NH4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw)*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
    trcs_Transp2MicP_3D(idg_NH3,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*NH3_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NH3_irrig_conc(I,NY,NX))*natomw)*trcs_VLN_vr(ids_NH4,NU(NY,NX),NY,NX)
    trcs_Transp2MicP_3D(ids_NO3,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*NO3_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw)*trcs_VLN_vr(ids_NO3,NU(NY,NX),NY,NX)

    trcs_Transp2MicP_3D(ids_NO2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_Transp2MicP_3D(ids_H1PO4,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*HPO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    trcs_Transp2MicP_3D(ids_H2PO4,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*H2PO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN_vr(ids_H1PO4,NU(NY,NX),NY,NX)
    trcs_Transp2MicP_3D(ids_NH4B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*NH4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw)*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)
    trcs_Transp2MicP_3D(idg_NH3B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*NH3_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NH3_irrig_conc(I,NY,NX))*natomw)*trcs_VLN_vr(ids_NH4B,NU(NY,NX),NY,NX)
    trcs_Transp2MicP_3D(ids_NO3B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*NO3_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw)*trcs_VLN_vr(ids_NO3B,NU(NY,NX),NY,NX)
    trcs_Transp2MicP_3D(ids_NO2B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_Transp2MicP_3D(ids_H1PO4B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*HPO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
    trcs_Transp2MicP_3D(ids_H2PO4B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf_col(NY,NX)*H2PO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN_vr(ids_H1PO4B,NU(NY,NX),NY,NX)
!
!     NO SOLUTE FLUXES FROM ATMOSPHERE
!
  ELSE
    trcg_Xbndl_flx(idg_beg:idg_end-1,1,NY,NX)=0.0_r8
    trcn_Xbndl_flx(ids_nut_beg:ids_nuts_end,1,NY,NX)=0.0_r8

    trcs_Transp2MicP_3D(idg_beg:idg_end-1,3,0,NY,NX)=0.0_r8
    trcs_Transp2MicP_3D(ids_nut_beg:ids_nuts_end,3,0,NY,NX)=0.0_r8
    trcs_Transp2MicP_3D(ids_beg:ids_end,3,NU(NY,NX),NY,NX)=0.0_r8

  ENDIF

  trcs_Transp2MacP_3D(ids_beg:ids_end,3,NU(NY,NX),NY,NX)=0.0_r8

  end subroutine HourlySoluteFluxes
!------------------------------------------------------------------------------------------

  subroutine SubHourlyFluxesFromSiteFile(NY,NX)
  use EcoSiMParDataMod, only : micpar
  implicit none

  integer, intent(in) :: NY, NX
  INTEGER :: K,L,NTS,NTG,NTN,idom

  DO  K=1,micpar%NumOfLitrCmplxs
    DO idom=idom_beg,idom_end
      RDOMFL0(idom,K,NY,NX)=DOM_MicpTransp_3D(idom,K,3,0,NY,NX)*dts_HeatWatTP    
      RDOMFL1(idom,K,NY,NX)=DOM_MicpTransp_3D(idom,K,3,NU(NY,NX),NY,NX)*dts_HeatWatTP
    ENDDO
  enddo

  DO NTG=idg_beg,idg_end-1
    trcg_advW_snvr(idg_CO2,1,NY,NX)=trcg_Xbndl_flx(idg_CO2,1,NY,NX)*dts_HeatWatTP
  ENDDO

  DO NTN=ids_nut_beg,ids_nuts_end
    trcn_RBLS(NTN,1,NY,NX)=trcn_Xbndl_flx(NTN,1,NY,NX)*dts_HeatWatTP
  ENDDO

  DO NTG=idg_beg,idg_end-1
    trcg_RFL0(NTG,NY,NX)=trcs_Transp2MicP_3D(NTG,3,0,NY,NX)*dts_HeatWatTP
  ENDDO

  do nts=ids_nut_beg,ids_nuts_end
    trcn_RFL0(nts,NY,NX)=trcs_Transp2MicP_3D(nts,3,0,NY,NX)*dts_HeatWatTP
  enddo

  DO NTS=ids_beg,ids_end
    trcs_RFL1(NTS,NY,NX)=trcs_Transp2MicP_3D(NTS,3,NU(NY,NX),NY,NX)*dts_HeatWatTP
  ENDDO
!
!     GAS AND SOLUTE SINKS AND SOURCES IN SOIL LAYERS FROM MICROBIAL
!     TRANSFORMATIONS IN 'NITRO' + ROOT EXCHANGE IN 'EXTRACT'
!     + EQUILIBRIA REACTIONS IN 'SOLUTE' AT SUB-HOURLY TIME STEP
!
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *SGL*=solute diffusivity from hour1.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     PARR=boundary layer conductance above litter from watsub.f
!     dt_GasCyc=1/number of cycles NPH-1 for gas flux calculations
!
  DO NTS=ids_beg,ids_end
    SoluteDifusvty_vrc(NTS,0,NY,NX)=SoluteDifusvty_vr(NTS,0,NY,NX)*dts_HeatWatTP
  ENDDO

  DO idom=idom_beg,idom_end
    DOMdiffusivity2_vr(idom,0,NY,NX)=DOMdiffusivity_vr(idom,0,NY,NX)*dts_HeatWatTP
  enddo
!
!     INITIAL SOLUTES IN SNOWPACK
!
!     CO2W,CH4W,OXYW,ZNGW,ZN2W,ZN4W,ZN3W,ZNOW,Z1PW,ZHPW=CO2,CH4,O2,N2,N2O,H2 content in snowpack
!
  DO L=1,JS

    trcg_solsml2_snvr(idg_beg:idg_end-1,L,NY,NX)=trcg_solsml_snvr(idg_beg:idg_end-1,L,NY,NX)

    trcn_solsml2(ids_nut_beg:ids_nuts_end,L,NY,NX)=trcn_solsml(ids_nut_beg:ids_nuts_end,L,NY,NX)
  enddo
  end subroutine SubHourlyFluxesFromSiteFile
!------------------------------------------------------------------------------------------

  subroutine ImportFluxFromOutsideModules(I,NY,NX,RFLZ_sol_vr)

  implicit none
  integer, intent(in) ::  I,NY, NX
  real(r8), intent(out) :: RFLZ_sol_vr(ids_beg:ids_end,JZ,JY,JX)

  real(r8) :: FLWU(JZ,JY,JX)
  integer :: L,K,NTG,NTS,idom

  DO L=NU(NY,NX),NL(NY,NX)
    FLWU(L,NY,NX)=TPlantRootH2OUptake_vr(L,NY,NX)*dts_HeatWatTP
    RBGCSinkG_vr(idg_CO2,L,NY,NX) = (trcg_RMicbTransf_vr(idg_CO2,L,NY,NX) +trcs_plant_uptake_vr(idg_CO2,L,NY,NX)-TR_CO2_aqu_soil_vr(L,NY,NX))*dts_gas
    RBGCSinkG_vr(idg_CH4,L,NY,NX) = (trcg_RMicbTransf_vr(idg_CH4,L,NY,NX) +trcs_plant_uptake_vr(idg_CH4,L,NY,NX))*dts_gas
    RBGCSinkG_vr(idg_N2,L,NY,NX)  = (trcg_RMicbTransf_vr(idg_N2,L,NY,NX)  +trcs_plant_uptake_vr(idg_N2,L,NY,NX)+Micb_N2Fixation_vr(L,NY,NX))*dts_gas
    RBGCSinkG_vr(idg_N2O,L,NY,NX) = (trcg_RMicbTransf_vr(idg_N2O,L,NY,NX) +trcs_plant_uptake_vr(idg_N2O,L,NY,NX))*dts_gas
    RBGCSinkG_vr(idg_H2,L,NY,NX)  = (trcg_RMicbTransf_vr(idg_H2,L,NY,NX)  +trcs_plant_uptake_vr(idg_H2,L,NY,NX))*dts_gas
    RBGCSinkG_vr(idg_NH3,L,NY,NX) = -TRN3G(L,NY,NX)*dts_gas

!    DO  K=1,jcplx
!      DO idom=idom_beg,idom_end
!        RDOM_CumEcoProd_vr(idom,K,L,NY,NX)=-REcoDOMProd_vr(idom,K,L,NY,NX)*dts_HeatWatTP
!      enddo
!    ENDDO

    RBGCSinkS_vr(ids_NH4,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NH4,L,NY,NX)-trcn_RChem_soil_vr(ids_NH4,L,NY,NX)+trcs_plant_uptake_vr(ids_NH4,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS_vr(idg_NH3,L,NY,NX)   = (-TR_NH3_soil_vr(L,NY,NX)+trcs_plant_uptake_vr(idg_NH3,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS_vr(ids_NO3,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NO3,L,NY,NX)-trcn_RChem_soil_vr(ids_NO3,L,NY,NX)+trcs_plant_uptake_vr(ids_NO3,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS_vr(ids_NO2,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NO2,L,NY,NX)-trcn_RChem_soil_vr(ids_NO2,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS_vr(ids_H2PO4,L,NY,NX) = (-RNutMicbTransf_vr(ids_H2PO4,L,NY,NX)-trcn_RChem_soil_vr(ids_H2PO4,L,NY,NX)+trcs_plant_uptake_vr(ids_H2PO4,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS_vr(ids_H1PO4,L,NY,NX) = (-RNutMicbTransf_vr(ids_H1PO4,L,NY,NX)-trcn_RChem_soil_vr(ids_H1PO4,L,NY,NX)+trcs_plant_uptake_vr(ids_H1PO4,L,NY,NX))*dts_HeatWatTP

    RBGCSinkS_vr(ids_NH4B,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NH4B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_NH4B,L,NY,NX)+trcs_plant_uptake_vr(ids_NH4B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS_vr(idg_NH3B,L,NY,NX)   = (-trcn_RChem_band_soil_vr(idg_NH3B,L,NY,NX)+trcs_plant_uptake_vr(idg_NH3B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS_vr(ids_NO3B,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NO3B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_NO3B,L,NY,NX)+trcs_plant_uptake_vr(ids_NO3B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS_vr(ids_NO2B,L,NY,NX)   = (-RNutMicbTransf_vr(ids_NO2B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_NO2B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS_vr(ids_H2PO4B,L,NY,NX) = (-RNutMicbTransf_vr(ids_H2PO4B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_H2PO4B,L,NY,NX)+trcs_plant_uptake_vr(ids_H2PO4B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS_vr(ids_H1PO4B,L,NY,NX) = (-RNutMicbTransf_vr(ids_H1PO4B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_H1PO4B,L,NY,NX)+trcs_plant_uptake_vr(ids_H1PO4B,L,NY,NX))*dts_HeatWatTP
!
!     SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     FLU=subsurface water flux from watsub.f
!     R*FLU,R*FBU=subsurface solute flux in non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     C*Q=irrigation solute concentrations
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    trcs_Irrig_vr(idg_CO2,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*CO2_irrig_conc(NY,NX)
    trcs_Irrig_vr(idg_CH4,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*CH4_irrig_conc(NY,NX)
    trcs_Irrig_vr(idg_O2,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*O2_irrig_conc(NY,NX)
    trcs_Irrig_vr(idg_N2,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*N2_irrig_conc(NY,NX)
    trcs_Irrig_vr(idg_N2O,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*N2O_irrig_conc(NY,NX)
    trcs_Irrig_vr(idg_H2,L,NY,NX)=0.0_r8

    trcs_Irrig_vr(ids_NH4,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*NH4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)*natomw
    trcs_Irrig_vr(idg_NH3,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*NH3_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_NO3,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*NO3_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_NO3,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_H1PO4,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*HPO4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)*patomw
    trcs_Irrig_vr(ids_H2PO4,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*H2PO4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)*patomw
    trcs_Irrig_vr(ids_NH4B,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*NH4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)*natomw
    trcs_Irrig_vr(idg_NH3B,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*NH3_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_NO3B,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*NO3_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_NO3B,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_H1PO4B,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*HPO4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)*patomw
    trcs_Irrig_vr(ids_H2PO4B,L,NY,NX)=FWatIrrigate2MicP_vr(L,NY,NX)*H2PO4_irrig_conc(I,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)*patomw
!
!     SUB-HOURLY SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     R*FLZ,R*FBZ=subsurface solute flux in non-band,band
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!
    DO NTS=ids_beg,ids_end
      RFLZ_sol_vr(NTS,L,NY,NX)=trcs_Irrig_vr(NTS,L,NY,NX)*dts_HeatWatTP
    ENDDO
!
!     GAS AND SOLUTE DIFFUSIVITIES AT SUB-HOURLY TIME STEP
!
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!     dts_gas=1/number of cycles h-1 for gas flux calculations
!     *SGL*=solute diffusivity from hour1.f
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!
    DO idom=idom_beg,idom_end
      DOMdiffusivity2_vr(idom,L,NY,NX)=DOMdiffusivity_vr(idom,L,NY,NX)*dts_HeatWatTP
    enddo

    DO NTS=ids_beg,ids_end
      SoluteDifusvty_vrc(NTS,L,NY,NX)=SoluteDifusvty_vr(NTS,L,NY,NX)*dts_HeatWatTP
    ENDDO

    DO NTG=idg_beg,idg_end
      GasDifc_vrc(NTG,L,NY,NX)=GasDifc_vr(NTG,L,NY,NX)*dts_gas
    ENDDO
!
!     STATE VARIABLES FOR GASES AND SOLUTES USED IN 'TranspNoSalt'
!     TO STORE SUB-HOURLY CHANGES DURING FLUX CALCULATIONS
!     INCLUDING TRANSFORMATIONS FROM NITRO.F, UPTAKE.F AND SOLUTE.F
!
!     CO2G,CH4G,OXYG,ZN3G,Z2GG,Z2OG,H2GG=gaseous CO2,CH4,O2,NH3,N2,N2O,H2
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 in micropores

!     OQCH,OQNH,OQPH,OQAH=DOC,DON,DOP,acetate in macropores
!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in non-band micropores
!     ZNH4B,ZNH3B,ZNO3B,ZNO2B,H1POB,H2POB=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in band micropores
!     CO2SH,CH4SH,OXYSH,Z2GSH,Z2OSH,H2GSH=aqueous CO2,CH4,O2,N2,N2O,H2 content in macropores
!     ZNH4SH,ZNH3SH,ZNO3SH,ZNO2SH,H1PO4H,H2PO4H=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in non-band macropores
!     ZNH4BH,ZNH3BH,ZNO3BH,ZNO2BH,H1POBH,H2POBH=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4 in band macropores
! exclude NH3B
    DO NTG=idg_beg,idg_end-1
      trc_gasml2_vr(NTG,L,NY,NX)=trc_gasml_vr(NTG,L,NY,NX)
    ENDDO

    !reset DOM to value before the iteration
    !Concern: RDOMMicProd_vr includes microbial modification of DOM
    !it does not include exchange among different complexes (as one kind of microbial priming). 
    DO  K=1,jcplx
      DO idom=idom_beg,idom_end
        DOM_MicP2(idom,K,L,NY,NX)=DOM_vr(idom,K,L,NY,NX)
!        if(DOM_MicP2(idom,K,L,NY,NX)<0._r8)then
!        print*,'nosaltmode',DOM_MicP2(idom,K,L,NY,NX),RDOMMicProd_vr(idom,K,L,NY,NX)
!        endif
        DOM_MacP2(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)
      ENDDO
    enddo
    DO NTS=ids_beg,ids_end
      trc_solml2_vr(NTS,L,NY,NX)=fixnegmass(trc_solml_vr(NTS,L,NY,NX))
      trc_soHml2_vr(NTS,L,NY,NX)=fixnegmass(trc_soHml_vr(NTS,L,NY,NX))

      if(trc_solml2_vr(NTS,L,NY,NX)<0._r8)then
      write(*,*)'micv2 ',trcs_names(NTS),L,trc_solml2_vr(NTS,L,NY,NX),trc_solml_vr(NTS,L,NY,NX)
      stop
      endif
      if(trc_soHml_vr(NTS,L,NY,NX)<0._r8)then
      write(*,*)'macv2 ',trcs_names(NTS),L,trc_soHml_vr(NTS,L,NY,NX)
      endif
    ENDDO

  enddo
  end subroutine ImportFluxFromOutsideModules
!

end module TranspNoSaltMod
