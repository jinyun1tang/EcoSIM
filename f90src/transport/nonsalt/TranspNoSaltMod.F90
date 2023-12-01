module TranspNoSaltMod
!!
! Description:
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
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
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  real(r8), allocatable :: CHY0(:,:,:)
  real(r8), allocatable :: RFLZ_sol(:,:,:,:)

  real(r8), PARAMETER :: DPN4=5.7E-07
  REAL(r8) :: CCO2SQ,CCH4SQ,COXYSQ,CZ2GSQ,CZ2OSQ,CNH3SQ,CNH3BQ,CH2GSQ

  public :: TranspNoSalt
  public :: InitTranspNoSalt,DestructTranspNoSalt
  contains

  subroutine InitTranspNoSalt()
  implicit none

  allocate(CHY0(0:JZ,JY,JX))
  allocate(RFLZ_sol(ids_beg:ids_end,JZ,JY,JX))

  call InitTransfrData

  call InitInsTp
  end subroutine InitTranspNoSalt
!------------------------------------------------------------------------------------------
  subroutine DestructTranspNoSalt
  implicit none


  call destroy(CHY0)

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
!     execution begins here
!
!     TIME STEPS FOR SOLUTE AND GAS FLUX CALCULATIONS
!

  WaterFlow2Soil(:,:,:,:)=0._r8
  call InitFluxandStateVariables(I,NHW,NHE,NVN,NVS)

!
! TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
! NPH=NPX,NPT=NPY
! NPH=no. of cycles per hour for water, heat and solute flux calculations
! NPG=number of cycles per hour for gas flux calculations, NPG=NPH*NPT
! dt_GasCyc=1/number of cycles NPH-1 for gas flux calculations, =1.0_r8/NPT
  MX=0
  DO  MM=1,NPG
    !compute the ordinal number of the intermediate forcing used for 
    !current iteration
    M=MIN(NPH,INT((MM-1)*dt_GasCyc)+1)

    call ModelTracerHydroFlux(M,MX,NHW, NHE, NVN, NVS,WaterFlow2Soil)
!
!     BOUNDARY SOLUTE AND GAS FLUXES
!
    call BoundaryFlux(M,MX,NHW,NHE,NVN,NVS)
!
!     UPDATE STATE VARIABLES FROM TOTAL FLUXES CALCULATED ABOVE
!
    call UpdateStateVar(M,MM,MX,NPG,NHW,NHE,NVN,NVS)
    !MX records the M location of the last NPT iterations
    MX=M
  ENDDO
  RETURN
  END subroutine TranspNoSalt
!------------------------------------------------------------------------------------------

  subroutine InitFluxandStateVariables(I,NHW,NHE,NVN,NVS)
  implicit none

  integer, intent(in) :: I, NHW,NHE,NVN,NVS

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
!     RCO2O=net soil CO2 uptake from nitro.f
!     RCH4O=net soil CH4 uptake from nitro.f
!     RN2G=total soil N2 production from nitro.f
!     XN2GS=total N2 fixation from nitro.f
!     RN2O=net soil N2O uptake from nitro.f
!     RH2GO=net H2 uptake from nitro.f
!     RNH4MicbTransf_vr=net change in NH4 from nitro.f
!     TR_NH4_soil,TR_NH3_soil_vr=NH4,NH3 dissolution from solute.f
!     RNO3MicbTransf_vr=net change in NO3 from nitro.f
!     TRNO3=NO3 dissolution from solute.f
!     RNO2MicbTransf_vr=net change in NO2 from nitro.f
!     TRNO2=NO2 dissolution from solute.f
!     RH2PO4MicbTransf_vr=net change in H2PO4 from nitro.f
!     TR_H2PO4_soil=H2PO4 dissolution from solute.f
!     RH1PO4MicbTransf_vr=net change in HPO4 from nitro.f
!     TR_H1PO4_soil=HPO4 dissolution from solute.f
!
      call SurfaceSinksandSources(NY,NX)
!
!     INITIALIZE STATE VARIABLES FOR USE IN GAS, SOLUTE FLUX CALCULATIONS
!
!     CO2S,CH4S,OXYS,Z2GS,Z2OS,H2GS=aqueous CO2,CH4,O2,N2,N2O,H2 content

!     ZNH4S,ZNH3S,ZNO3S,ZNO2S,H1PO4,H2PO4=aqueous NH4,NH3,NO3,NO2,HPO4,H2PO4
!     CHY0=H concentration
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
!     VLSnowHeatCapM,VLHeatCapSnowMin=current,minimum volumetric heat capacity of snowpack
!     X*BLS=hourly solute flux to snowpack
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     FLQGQ,FLQGI=water flux to snowpack from rain,irrigation
!     C*R,C*Q=precipitation,irrigation solute concentrations
!     gas code: *CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*H2*=H2
!
      call HourlySoluteFluxes(I,NY,NX)
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
!     CHY0=H concentration
!     PH=pH
!     FLWU,TUPWTR=total root water uptake from extract.f
!     R*SK2=total sink from nitro.f, uptake.f, solute.f
!     RCO2O=net soil CO2 uptake from nitro.f
!     RCH4O=net soil CH4 uptake from nitro.f
!     RN2G=total soil N2 production from nitro.f
!     XN2GS=total N2 fixation from nitro.f
!     RN2O=net soil N2O uptake from nitro.f
!     RH2GO=net H2 uptake from nitro.f
!     RNH4MicbTransf_vr=net change in NH4 from nitro.f
!     TR_NH4_soil,TR_NH3_soil_vr=NH4,NH3 dissolution from solute.f
!     RNO3MicbTransf_vr=net change in NO3 from nitro.f
!     TRNO3=NO3 dissolution from solute.f
!     RNO2MicbTransf_vr=net change in NO2 from nitro.f
!     TRNO2=NO2 dissolution from solute.f
!     RH2PO4MicbTransf_vr=net change in H2PO4 from nitro.f
!     TR_H2PO4_soil=H2PO4 dissolution from solute.f
!     RH1PO4MicbTransf_vr=net change in HPO4 from nitro.f
!     TR_H1PO4_soil=HPO4 dissolution from solute.f
!     TRootNH4Uptake_pft,TUPNHB=root NH4 uptake in non-band,band from extract.f
!     TRootNO3Uptake_pft,TUPNOB=root NO3 uptake in non-band,band from extract.f
!     TRootH2PO4Uptake_pft,TUPH2B=root H2PO4 uptake in non-band,band from extract.f
!     TRootHPO4Uptake_pft,TUPH1B=root HPO4 uptake in non-band,band from extract.f
!
      call ImportFluxFromOutsideModules(I,NY,NX)
    ENDDO
  ENDDO
  end subroutine InitFluxandStateVariables
!------------------------------------------------------------------------------------------

  subroutine UpdateStateVar(M,MM,MX,NPG,NHW,NHE,NVN,NVS)
  implicit none

  integer, intent(in) ::M, MM, MX,NPG, NHW, NHE, NVN, NVS
  integer :: NY,NX,K,L,NGT

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
        call UpdateSolutesInSoilLayers(M,MX,NY,NX)
!
!     GAS EXCHANGE IN SURFACE LITTER
!
!     *S2=litter aqueous gas content
!     R*DFG=litter-atmosphere gas flux
!     gas code:*CO2*=CO2,*OXY*=O2,*CH4*=CH4,*Z2G*=N2,*Z2O*=N2O
!             :*ZN3*=NH3,*H2G*=H2
!
        DO NGT=idg_beg,idg_end
          trc_solml2(NGT,0,NY,NX)=trc_solml2(NGT,0,NY,NX)+RGasDSFlx(NGT,0,NY,NX)
        ENDDO
      ENDDO D9690
    ENDDO D9695
  ENDIF
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
    trcg_solsml2(NTG,1,NY,NX)=trcg_solsml2(NTG,1,NY,NX)+trcg_SnowDrift(NTG,NY,NX)
  ENDDO

  DO NTN=ids_nut_beg,ids_nuts_end
    trcn_solsml2(NTN,1,NY,NX)=trcn_solsml2(NTN,1,NY,NX)+trcn_SnowDrift(NTN,NY,NX)
  ENDDO

  DO  L=1,JS
    DO NTG=idg_beg,idg_end-1
      trcg_solsml2(NTG,L,NY,NX)=trcg_solsml2(NTG,L,NY,NX)+trcg_TBLS(NTG,L,NY,NX)
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

  integer :: K,NTG,NTS
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
    DOM_MicP2(idom_doc,K,0,NY,NX)=DOM_MicP2(idom_doc,K,0,NY,NX)+DOM_3DMicp_Transp_flxM(idom_doc,K,3,0,NY,NX)
    DOM_MicP2(idom_don,K,0,NY,NX)=DOM_MicP2(idom_don,K,0,NY,NX)+DOM_3DMicp_Transp_flxM(idom_don,K,3,0,NY,NX)
    DOM_MicP2(idom_dop,K,0,NY,NX)=DOM_MicP2(idom_dop,K,0,NY,NX)+DOM_3DMicp_Transp_flxM(idom_dop,K,3,0,NY,NX)
    DOM_MicP2(idom_acetate,K,0,NY,NX)=DOM_MicP2(idom_acetate,K,0,NY,NX)+DOM_3DMicp_Transp_flxM(idom_acetate,K,3,0,NY,NX)
  ENDDO
!exclude NH3B
  DO NTG=idg_beg,idg_end-1
    trc_solml2(NTG,0,NY,NX)=trc_solml2(NTG,0,NY,NX)+RDFR_gas(NTG,NY,NX)
  ENDDO

! band does not exist in litter layer
  DO NTS=ids_beg,ids_end
    trc_solml2(NTS,0,NY,NX)=trc_solml2(NTS,0,NY,NX)+R3PoreSolFlx(NTS,3,0,NY,NX)
  ENDDO
! include NH3B
  DO NTG=idg_beg,idg_end
    trc_solml2(NTG,NU(NY,NX),NY,NX)=trc_solml2(NTG,NU(NY,NX),NY,NX)+RGasSSVol(NTG,NY,NX)
  ENDDO

  D9680: DO K=1,jcplx
    DOM_MicP2(idom_doc,K,0,NY,NX)=DOM_MicP2(idom_doc,K,0,NY,NX)+dom_TFloXSurRunoff(idom_doc,K,NY,NX)
    DOM_MicP2(idom_don,K,0,NY,NX)=DOM_MicP2(idom_don,K,0,NY,NX)+dom_TFloXSurRunoff(idom_don,K,NY,NX)
    DOM_MicP2(idom_dop,K,0,NY,NX)=DOM_MicP2(idom_dop,K,0,NY,NX)+dom_TFloXSurRunoff(idom_dop,K,NY,NX)
    DOM_MicP2(idom_acetate,K,0,NY,NX)=DOM_MicP2(idom_acetate,K,0,NY,NX)+dom_TFloXSurRunoff(idom_acetate,K,NY,NX)
  ENDDO D9680
!exclude NH3B
  DO NTG=idg_beg,idg_end-1
    trc_solml2(NTG,0,NY,NX)=trc_solml2(NTG,0,NY,NX)+trcg_TFloXSurRunoff(NTG,NY,NX)
  ENDDO

  DO NTS=ids_nut_beg,ids_nuts_end
    trc_solml2(NTS,0,NY,NX)=trc_solml2(NTS,0,NY,NX)+trcn_TFloXSurRunoff(NTS,NY,NX)
  ENDDO

  end subroutine UpdateStateVarInResidue

!------------------------------------------------------------------------------------------

  subroutine UpdateSolutesInSoilLayers(M,MX,NY,NX)
  implicit none
  integer, intent(in) :: M,NY,NX,MX

  integer :: L,K,NTG,NTS

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
      IF(VLSoilPoreMicP(L,NY,NX).GT.ZEROS2(NY,NX))THEN

        DO NTS=ids_beg,ids_end
          trc_solml2(NTS,L,NY,NX)=trc_solml2(NTS,L,NY,NX)+RFLZ_sol(NTS,L,NY,NX)
        ENDDO

        trc_solml2(idg_CO2,L,NY,NX)=trc_solml2(idg_CO2,L,NY,NX)+trcg_BBL(idg_CO2,L,NY,NX)
        trc_solml2(idg_CH4,L,NY,NX)=trc_solml2(idg_CH4,L,NY,NX)+trcg_BBL(idg_CH4,L,NY,NX)
        trc_solml2(idg_O2,L,NY,NX)=trc_solml2(idg_O2,L,NY,NX)+trcg_BBL(idg_O2,L,NY,NX)
        trc_solml2(idg_N2,L,NY,NX)=trc_solml2(idg_N2,L,NY,NX)+trcg_BBL(idg_N2,L,NY,NX)
        trc_solml2(idg_N2O,L,NY,NX)=trc_solml2(idg_N2O,L,NY,NX)+trcg_BBL(idg_N2O,L,NY,NX)
        trc_solml2(idg_NH3,L,NY,NX)=trc_solml2(idg_NH3,L,NY,NX)+trcg_BBL(idg_NH3,L,NY,NX)
        trc_solml2(idg_NH3B,L,NY,NX)=trc_solml2(idg_NH3B,L,NY,NX)+trcg_BBL(idg_NH3B,L,NY,NX)
        trc_solml2(idg_H2,L,NY,NX)=trc_solml2(idg_H2,L,NY,NX)+trcg_BBL(idg_H2,L,NY,NX)

        DO  K=1,jcplx
          DOM_MicP2(idom_doc,K,L,NY,NX)=DOM_MicP2(idom_doc,K,L,NY,NX)+DOM_Transp2Micp_flx(idom_doc,K,L,NY,NX)+DOM_XPoreTransp_flx(idom_doc,K,L,NY,NX)
          DOM_MicP2(idom_don,K,L,NY,NX)=DOM_MicP2(idom_don,K,L,NY,NX)+DOM_Transp2Micp_flx(idom_don,K,L,NY,NX)+DOM_XPoreTransp_flx(idom_don,K,L,NY,NX)
          DOM_MicP2(idom_dop,K,L,NY,NX)=DOM_MicP2(idom_dop,K,L,NY,NX)+DOM_Transp2Micp_flx(idom_dop,K,L,NY,NX)+DOM_XPoreTransp_flx(idom_dop,K,L,NY,NX)
          DOM_MicP2(idom_acetate,K,L,NY,NX)=DOM_MicP2(idom_acetate,K,L,NY,NX)+DOM_Transp2Micp_flx(idom_acetate,K,L,NY,NX)+DOM_XPoreTransp_flx(idom_acetate,K,L,NY,NX)
          DOM_MacP2(idom_doc,K,L,NY,NX)=DOM_MacP2(idom_doc,K,L,NY,NX)+DOM_Transp2Macp_flx(idom_doc,K,L,NY,NX)-DOM_XPoreTransp_flx(idom_doc,K,L,NY,NX)
          DOM_MacP2(idom_don,K,L,NY,NX)=DOM_MacP2(idom_don,K,L,NY,NX)+DOM_Transp2Macp_flx(idom_don,K,L,NY,NX)-DOM_XPoreTransp_flx(idom_don,K,L,NY,NX)
          DOM_MacP2(idom_dop,K,L,NY,NX)=DOM_MacP2(idom_dop,K,L,NY,NX)+DOM_Transp2Macp_flx(idom_dop,K,L,NY,NX)-DOM_XPoreTransp_flx(idom_dop,K,L,NY,NX)
          DOM_MacP2(idom_acetate,K,L,NY,NX)=DOM_MacP2(idom_acetate,K,L,NY,NX)+DOM_Transp2Macp_flx(idom_acetate,K,L,NY,NX)-DOM_XPoreTransp_flx(idom_acetate,K,L,NY,NX)
        ENDDO
        DO NTS=ids_beg,ids_end
          trc_solml2(NTS,L,NY,NX)=trc_solml2(NTS,L,NY,NX) &
            +R3PorTSolFlx(NTS,L,NY,NX)+RporeSoXFlx(NTS,L,NY,NX)
          trc_soHml2(NTS,L,NY,NX)=trc_soHml2(NTS,L,NY,NX) &
            +R3PorTSoHFlx(NTS,L,NY,NX)-RporeSoXFlx(NTS,L,NY,NX)
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
    IF(VLSoilPoreMicP(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      DO NTG=idg_beg,idg_end
        trc_solml2(NTG,L,NY,NX)=trc_solml2(NTG,L,NY,NX)+RGasDSFlx(NTG,L,NY,NX)
      ENDDO

      DO NTG=idg_beg,idg_end-1
        trc_gasml2(NTG,L,NY,NX)=trc_gasml2(NTG,L,NY,NX) &
          +Gas_AdvDif_Flx_vr(NTG,L,NY,NX)-RGasDSFlx(NTG,L,NY,NX)
      ENDDO
      trc_gasml2(idg_NH3,L,NY,NX)=trc_gasml2(idg_NH3,L,NY,NX)-RGasDSFlx(idg_NH3B,L,NY,NX)
    ENDIF

  ENDDO D9685
  end subroutine UpdateSolutesInSoilLayers

!------------------------------------------------------------------------------------------

  subroutine SurfaceSinksandSources(NY, NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K

  RBGCSinkG(idg_CO2,0,NY,NX)=trcg_RMicbTransf_vr(idg_CO2,0,NY,NX)*dts_gas
  RBGCSinkG(idg_CH4,0,NY,NX)=trcg_RMicbTransf_vr(idg_CH4,0,NY,NX)*dts_gas
  RBGCSinkG(idg_N2,0,NY,NX)=(trcg_RMicbTransf_vr(idg_N2,0,NY,NX)+Micb_N2Fixation_vr(0,NY,NX))*dts_gas
  RBGCSinkG(idg_N2O,0,NY,NX)=trcg_RMicbTransf_vr(idg_N2O,0,NY,NX)*dts_gas
  RBGCSinkG(idg_NH3,0,NY,NX)=0.0_r8
  RBGCSinkG(idg_H2,0,NY,NX)=trcg_RMicbTransf_vr(idg_H2,0,NY,NX)*dts_gas
  DO  K=1,jcplx
    RDOM_micb_cumflx(idom_doc,K,0,NY,NX)=-RDOM_micb_flx(idom_doc,K,0,NY,NX)*dts_HeatWatTP
    RDOM_micb_cumflx(idom_don,K,0,NY,NX)=-RDOM_micb_flx(idom_don,K,0,NY,NX)*dts_HeatWatTP
    RDOM_micb_cumflx(idom_dop,K,0,NY,NX)=-RDOM_micb_flx(idom_dop,K,0,NY,NX)*dts_HeatWatTP
    RDOM_micb_cumflx(idom_acetate,K,0,NY,NX)=-RDOM_micb_flx(idom_acetate,K,0,NY,NX)*dts_HeatWatTP
  ENDDO
  RBGCSinkS(ids_NH4,0,NY,NX)=(-RNutMicbTransf_vr(ids_NH4,0,NY,NX)-trcn_RChem_soil_vr(ids_NH4,0,NY,NX))*dts_HeatWatTP
  RBGCSinkS(idg_NH3,0,NY,NX)=-TR_NH3_soil_vr(0,NY,NX)*dts_HeatWatTP
  RBGCSinkS(ids_NO3,0,NY,NX)=(-RNutMicbTransf_vr(ids_NO3,0,NY,NX)-trcn_RChem_soil_vr(ids_NO3,0,NY,NX))*dts_HeatWatTP
  RBGCSinkS(ids_NO2,0,NY,NX)=(-RNutMicbTransf_vr(ids_NO2,0,NY,NX)-trcn_RChem_soil_vr(ids_NO2,0,NY,NX))*dts_HeatWatTP
  RBGCSinkS(ids_H2PO4,0,NY,NX)=(-RNutMicbTransf_vr(ids_H2PO4,0,NY,NX)-trcn_RChem_soil_vr(ids_H2PO4,0,NY,NX))*dts_HeatWatTP
  RBGCSinkS(ids_H1PO4,0,NY,NX)=(-RNutMicbTransf_vr(ids_H1PO4,0,NY,NX)-trcn_RChem_soil_vr(ids_H1PO4,0,NY,NX))*dts_HeatWatTP

  end subroutine SurfaceSinksandSources
!------------------------------------------------------------------------------------------

  subroutine StateVarforGasandSolute(NY,NX)
  implicit none

  integer, intent(in) :: NY, NX

  integer :: K,NTG,NTS

! exclude banded nutrient
  DO NTG=idg_beg,idg_end-1
    trc_solml2(NTG,0,NY,NX)=trc_solml(NTG,0,NY,NX)
  ENDDO

  D9979: DO K=1,jcplx
    DOM_MicP2(idom_doc,K,0,NY,NX)=DOM(idom_doc,K,0,NY,NX)-RDOM_micb_flx(idom_doc,K,0,NY,NX)
    DOM_MicP2(idom_don,K,0,NY,NX)=DOM(idom_don,K,0,NY,NX)-RDOM_micb_flx(idom_don,K,0,NY,NX)
    DOM_MicP2(idom_dop,K,0,NY,NX)=DOM(idom_dop,K,0,NY,NX)-RDOM_micb_flx(idom_dop,K,0,NY,NX)
    DOM_MicP2(idom_acetate,K,0,NY,NX)=DOM(idom_acetate,K,0,NY,NX)-RDOM_micb_flx(idom_acetate,K,0,NY,NX)
  ENDDO D9979

! exclude banded nutrient
  DO NTS=ids_nut_beg,ids_nuts_end
    trc_solml2(NTS,0,NY,NX)=trc_solml(NTS,0,NY,NX)
  ENDDO
  trc_solml2(ids_H1PO4B,0,NY,NX)=trc_solml2(ids_H1PO4,0,NY,NX)
  trc_solml2(ids_H2PO4B,0,NY,NX)=trc_solml2(ids_H2PO4,0,NY,NX)
  trc_solml2(ids_NO3B,0,NY,NX)=trc_solml2(ids_NO3,0,NY,NX)
  trc_solml2(ids_NO2B,0,NY,NX)=trc_solml2(ids_NO2,0,NY,NX)
  trc_solml2(ids_NH4B,0,NY,NX)=trc_solml2(ids_NH4,0,NY,NX)
  trc_solml2(idg_NH3B,0,NY,NX)=trc_solml2(idg_NH3,0,NY,NX)

  CHY0(0,NY,NX)=10.0_r8**(-(PH(0,NY,NX)-3.0_r8))
  end subroutine StateVarforGasandSolute
!------------------------------------------------------------------------------------------

  subroutine SurfaceSolutefromAtms(NY,NX)
  implicit none
  integer, intent(in) :: NY, NX
  integer :: K

  D8855: DO K=1,jcplx
    IF(K.LE.2)THEN
      DOM_3DMicp_Transp_flx(idom_beg:idom_end,K,3,0,NY,NX)=0.0
    ENDIF
    DOM_3DMicp_Transp_flx(idom_beg:idom_end,K,3,NU(NY,NX),NY,NX)=0.0
    DOM_3DMacp_Transp_flx(idom_beg:idom_end,K,3,NU(NY,NX),NY,NX)=0.0
  ENDDO D8855
  end subroutine SurfaceSolutefromAtms
!------------------------------------------------------------------------------------------

  subroutine HourlySoluteFluxes(I,NY,NX)
  implicit none

  integer, intent(in) :: I
  integer, intent(in) :: NY,NX
  IF(SnoFalPrec(NY,NX).GT.0.0_r8.OR.(RainFalPrec(NY,NX).GT.0.0_r8.AND.VLSnowHeatCapM(1,1,NY,NX).GT.VLHeatCapSnowMin(NY,NX)))THEN
    trcg_XBLS(idg_CO2,1,NY,NX)=Rain2SoilSurf(NY,NX)*CO2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*CO2_irrig_conc(NY,NX)
    trcg_XBLS(idg_CH4,1,NY,NX)=Rain2SoilSurf(NY,NX)*CH4_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*CH4_irrig_conc(NY,NX)
    trcg_XBLS(idg_O2,1,NY,NX)=Rain2SoilSurf(NY,NX)*O2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*O2_irrig_conc(NY,NX)
    trcg_XBLS(idg_N2,1,NY,NX)=Rain2SoilSurf(NY,NX)*N2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*N2_irrig_conc(NY,NX)
    trcg_XBLS(idg_N2O,1,NY,NX)=Rain2SoilSurf(NY,NX)*N2O_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*N2O_irrig_conc(NY,NX)
    trcg_XBLS(idg_NH3,1,NY,NX)=(Rain2SoilSurf(NY,NX)*NH3_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*NH3_irrig_conc(I,NY,NX))*natomw

    trcn_XBLS(ids_NH4,1,NY,NX)=(Rain2SoilSurf(NY,NX)*NH4_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw
    trcn_XBLS(ids_NO3,1,NY,NX)=(Rain2SoilSurf(NY,NX)*NO3_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw
    trcn_XBLS(ids_H1PO4,1,NY,NX)=(Rain2SoilSurf(NY,NX)*HPO4_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw
    trcn_XBLS(ids_H2PO4,1,NY,NX)=(Rain2SoilSurf(NY,NX)*H2PO4_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IF RAINFALL AND IRRIGATION IS ZERO IF SNOWPACK IS PRESENT
!
!     X*FLS,X*FLB=hourly solute flux to micropores in non-band,band
!
    trcs_3DTransp2MicP(ids_beg:ids_end,3,0,NY,NX)=0.0_r8
    trcs_3DTransp2MicP(ids_beg:ids_end,3,NU(NY,NX),NY,NX)=0.0_r8
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SOIL SURFACE
!     IN RAINFALL AND IRRIGATION ACCORDING TO CONCENTRATIONS
!     ENTERED IN WEATHER AND IRRIGATION FILES
!
  ELSEIF((PrecAtm(NY,NX).GT.0.0_r8.OR.IrrigSurface(NY,NX).GT.0.0_r8) &
    .AND.VLSnowHeatCapM(1,1,NY,NX).LE.VLHeatCapSnowMin(NY,NX))THEN
!
!     HOURLY SOLUTE FLUXES FROM ATMOSPHERE TO SNOWPACK
!     IF SNOWFALL AND IRRIGATION IS ZERO AND SNOWPACK IS ABSENT
!
!     PrecAtm,PRECI=snow+rain,irrigation
!     X*BLS=hourly solute flux to snowpack
!     X*FLS,X*FLB=hourly solute flux to surface litter,soil surface micropore non-band,band
!     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
!             :OC=DOC,ON=DON,OP=DOP,OA=acetate
!             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
!             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
!     Rain2LitRSurf,Irrig2LitRSurf=water flux to surface litter from rain,irrigation
!     FLQGQ,FLQGI=water flux to soil surface from rain,irrigation
!     C*R,C*Q=precipitation,irrigation solute concentrations
!     gas code: *CO*=CO2,*OX*=O2,*CH*=CH4,*NG*=N2,*N2*=N2O,*NH*=NH3,*H2*=H2
!
    trcg_XBLS(idg_beg:idg_end-1,1,NY,NX)=0.0_r8

    trcn_XBLS(ids_nut_beg:ids_nuts_end,1,NY,NX)=0.0_r8

    trcs_3DTransp2MicP(idg_CO2,3,0,NY,NX)=Rain2LitRSurf(NY,NX)*CO2_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*CO2_irrig_conc(NY,NX)
    trcs_3DTransp2MicP(idg_CH4,3,0,NY,NX)=Rain2LitRSurf(NY,NX)*CH4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*CH4_irrig_conc(NY,NX)
    trcs_3DTransp2MicP(idg_O2,3,0,NY,NX)=Rain2LitRSurf(NY,NX)*O2_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*O2_irrig_conc(NY,NX)
    trcs_3DTransp2MicP(idg_N2,3,0,NY,NX)=Rain2LitRSurf(NY,NX)*N2_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*N2_irrig_conc(NY,NX)
    trcs_3DTransp2MicP(idg_N2O,3,0,NY,NX)=Rain2LitRSurf(NY,NX)*N2O_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*N2O_irrig_conc(NY,NX)
    trcs_3DTransp2MicP(idg_H2,3,0,NY,NX)=0.0_r8

    trcs_3DTransp2MicP(ids_NH4,3,0,NY,NX)=(Rain2LitRSurf(NY,NX)*NH4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw
    trcs_3DTransp2MicP(idg_NH3,3,0,NY,NX)=(Rain2LitRSurf(NY,NX)*NH3_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*NH3_irrig_conc(I,NY,NX))*natomw
    trcs_3DTransp2MicP(ids_NO3,3,0,NY,NX)=(Rain2LitRSurf(NY,NX)*NO3_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw
    trcs_3DTransp2MicP(ids_NO2,3,0,NY,NX)=0.0_r8
    trcs_3DTransp2MicP(ids_H1PO4,3,0,NY,NX)=(Rain2LitRSurf(NY,NX)*HPO4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw
    trcs_3DTransp2MicP(ids_H2PO4,3,0,NY,NX)=(Rain2LitRSurf(NY,NX)*H2PO4_rain_conc(NY,NX)+Irrig2LitRSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw
    trcs_3DTransp2MicP(idg_CO2,3,NU(NY,NX),NY,NX)=Rain2SoilSurf(NY,NX)*CO2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*CO2_irrig_conc(NY,NX)
    trcs_3DTransp2MicP(idg_CH4,3,NU(NY,NX),NY,NX)=Rain2SoilSurf(NY,NX)*CH4_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*CH4_irrig_conc(NY,NX)
    trcs_3DTransp2MicP(idg_O2,3,NU(NY,NX),NY,NX)=Rain2SoilSurf(NY,NX)*O2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*O2_irrig_conc(NY,NX)
    trcs_3DTransp2MicP(idg_N2,3,NU(NY,NX),NY,NX)=Rain2SoilSurf(NY,NX)*N2_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*N2_irrig_conc(NY,NX)
    trcs_3DTransp2MicP(idg_N2O,3,NU(NY,NX),NY,NX)=Rain2SoilSurf(NY,NX)*N2O_rain_conc(NY,NX)+Irrig2SoilSurf(NY,NX)*N2O_irrig_conc(NY,NX)
    trcs_3DTransp2MicP(idg_H2,3,NU(NY,NX),NY,NX)=0.0_r8

    trcs_3DTransp2MicP(ids_NH4,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf(NY,NX)*NH4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw)*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
    trcs_3DTransp2MicP(idg_NH3,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf(NY,NX)*NH3_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NH3_irrig_conc(I,NY,NX))*natomw)*trcs_VLN(ids_NH4,NU(NY,NX),NY,NX)
    trcs_3DTransp2MicP(ids_NO3,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf(NY,NX)*NO3_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw)*trcs_VLN(ids_NO3,NU(NY,NX),NY,NX)
    trcs_3DTransp2MicP(ids_NO2,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_3DTransp2MicP(ids_H1PO4,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf(NY,NX)*HPO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    trcs_3DTransp2MicP(ids_H2PO4,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf(NY,NX)*H2PO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
    trcs_3DTransp2MicP(ids_NH4B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf(NY,NX)*NH4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NH4_irrig_conc(I,NY,NX))*natomw)*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
    trcs_3DTransp2MicP(idg_NH3B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf(NY,NX)*NH3_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NH3_irrig_conc(I,NY,NX))*natomw)*trcs_VLN(ids_NH4B,NU(NY,NX),NY,NX)
    trcs_3DTransp2MicP(ids_NO3B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf(NY,NX)*NO3_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*NO3_irrig_conc(I,NY,NX))*natomw)*trcs_VLN(ids_NO3B,NU(NY,NX),NY,NX)
    trcs_3DTransp2MicP(ids_NO2B,3,NU(NY,NX),NY,NX)=0.0_r8
    trcs_3DTransp2MicP(ids_H1PO4B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf(NY,NX)*HPO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*HPO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
    trcs_3DTransp2MicP(ids_H2PO4B,3,NU(NY,NX),NY,NX)=((Rain2SoilSurf(NY,NX)*H2PO4_rain_conc(NY,NX) &
      +Irrig2SoilSurf(NY,NX)*H2PO4_irrig_conc(I,NY,NX))*patomw)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
!
!     NO SOLUTE FLUXES FROM ATMOSPHERE
!
  ELSE
    trcg_XBLS(idg_beg:idg_end-1,1,NY,NX)=0.0_r8
    trcn_XBLS(ids_nut_beg:ids_nuts_end,1,NY,NX)=0.0_r8

    trcs_3DTransp2MicP(idg_beg:idg_end-1,3,0,NY,NX)=0.0_r8
    trcs_3DTransp2MicP(ids_nut_beg:ids_nuts_end,3,0,NY,NX)=0.0_r8
    trcs_3DTransp2MicP(ids_beg:ids_end,3,NU(NY,NX),NY,NX)=0.0_r8

  ENDIF

  trcs_3DTransp2MacP(ids_beg:ids_end,3,NU(NY,NX),NY,NX)=0.0_r8

  end subroutine HourlySoluteFluxes
!------------------------------------------------------------------------------------------

  subroutine SubHourlyFluxesFromSiteFile(NY,NX)
  use EcoSiMParDataMod, only : micpar
  implicit none

  integer, intent(in) :: NY, NX
  INTEGER :: K,L,NTS,NTG,NTN

  DO  K=1,micpar%NumOfLitrCmplxs
    ROCFL0(K,NY,NX)=DOM_3DMicp_Transp_flx(idom_doc,K,3,0,NY,NX)*dts_HeatWatTP
    RONFL0(K,NY,NX)=DOM_3DMicp_Transp_flx(idom_don,K,3,0,NY,NX)*dts_HeatWatTP
    ROPFL0(K,NY,NX)=DOM_3DMicp_Transp_flx(idom_dop,K,3,0,NY,NX)*dts_HeatWatTP
    ROAFL0(K,NY,NX)=DOM_3DMicp_Transp_flx(idom_acetate,K,3,0,NY,NX)*dts_HeatWatTP
    ROCFL1(K,NY,NX)=DOM_3DMicp_Transp_flx(idom_doc,K,3,NU(NY,NX),NY,NX)*dts_HeatWatTP
    RONFL1(K,NY,NX)=DOM_3DMicp_Transp_flx(idom_don,K,3,NU(NY,NX),NY,NX)*dts_HeatWatTP
    ROPFL1(K,NY,NX)=DOM_3DMicp_Transp_flx(idom_dop,K,3,NU(NY,NX),NY,NX)*dts_HeatWatTP
    ROAFL1(K,NY,NX)=DOM_3DMicp_Transp_flx(idom_acetate,K,3,NU(NY,NX),NY,NX)*dts_HeatWatTP
  enddo

  DO NTG=idg_beg,idg_end-1
    trcg_RBLS(idg_CO2,1,NY,NX)=trcg_XBLS(idg_CO2,1,NY,NX)*dts_HeatWatTP
  ENDDO

  DO NTN=ids_nut_beg,ids_nuts_end
    trcn_RBLS(NTN,1,NY,NX)=trcn_XBLS(NTN,1,NY,NX)*dts_HeatWatTP
  ENDDO

  DO NTG=idg_beg,idg_end-1
    trcg_RFL0(NTG,NY,NX)=trcs_3DTransp2MicP(NTG,3,0,NY,NX)*dts_HeatWatTP
  ENDDO

  trcn_RFL0(ids_NH4,NY,NX)=trcs_3DTransp2MicP(ids_NH4,3,0,NY,NX)*dts_HeatWatTP
  trcn_RFL0(ids_NO3,NY,NX)=trcs_3DTransp2MicP(ids_NO3,3,0,NY,NX)*dts_HeatWatTP
  trcn_RFL0(ids_NO2,NY,NX)=trcs_3DTransp2MicP(ids_NO2,3,0,NY,NX)*dts_HeatWatTP
  trcn_RFL0(ids_H1PO4,NY,NX)=trcs_3DTransp2MicP(ids_H1PO4,3,0,NY,NX)*dts_HeatWatTP
  trcn_RFL0(ids_H2PO4,NY,NX)=trcs_3DTransp2MicP(ids_H2PO4,3,0,NY,NX)*dts_HeatWatTP

  DO NTS=ids_beg,ids_end
    trcs_RFL1(NTS,NY,NX)=trcs_3DTransp2MicP(NTS,3,NU(NY,NX),NY,NX)*dts_HeatWatTP
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
    SolDifcc(NTS,0,NY,NX)=SolDifc(NTS,0,NY,NX)*dts_HeatWatTP
  ENDDO

  OCSGL2(0,NY,NX)=OCSGL(0,NY,NX)*dts_HeatWatTP
  ONSGL2(0,NY,NX)=ONSGL(0,NY,NX)*dts_HeatWatTP
  OPSGL2(0,NY,NX)=OPSGL(0,NY,NX)*dts_HeatWatTP
  OASGL2(0,NY,NX)=OASGL(0,NY,NX)*dts_HeatWatTP
!
!     INITIAL SOLUTES IN SNOWPACK
!
!     CO2W,CH4W,OXYW,ZNGW,ZN2W,ZN4W,ZN3W,ZNOW,Z1PW,ZHPW=CO2,CH4,O2,N2,N2O,H2 content in snowpack
!
  DO L=1,JS

    trcg_solsml2(idg_beg:idg_end-1,L,NY,NX)=trcg_solsml(idg_beg:idg_end-1,L,NY,NX)

    trcn_solsml2(ids_nut_beg:ids_nuts_end,L,NY,NX)=trcn_solsml(ids_nut_beg:ids_nuts_end,L,NY,NX)
  enddo
  end subroutine SubHourlyFluxesFromSiteFile
!------------------------------------------------------------------------------------------

  subroutine ImportFluxFromOutsideModules(I,NY,NX)
  implicit none

  integer, intent(in) ::  I,NY, NX

  real(r8) :: FLWU(JZ,JY,JX)
  integer :: L,K,NTG,NTS

  DO L=NU(NY,NX),NL(NY,NX)
    CHY0(L,NY,NX)=10.0_r8**(-(PH(L,NY,NX)-3.0_r8))
    FLWU(L,NY,NX)=GridPlantRootH2OUptake_vr(L,NY,NX)*dts_HeatWatTP
    RBGCSinkG(idg_CO2,L,NY,NX)=(trcg_RMicbTransf_vr(idg_CO2,L,NY,NX)+trcs_plant_uptake_vr(idg_CO2,L,NY,NX)-TR_CO2_aqu_soil_vr(L,NY,NX))*dts_gas
    RBGCSinkG(idg_CH4,L,NY,NX)=(trcg_RMicbTransf_vr(idg_CH4,L,NY,NX)+trcs_plant_uptake_vr(idg_CH4,L,NY,NX))*dts_gas
    RBGCSinkG(idg_N2,L,NY,NX)=(trcg_RMicbTransf_vr(idg_N2,L,NY,NX)+Micb_N2Fixation_vr(L,NY,NX)+trcs_plant_uptake_vr(idg_N2,L,NY,NX))*dts_gas
    RBGCSinkG(idg_N2O,L,NY,NX)=(trcg_RMicbTransf_vr(idg_N2O,L,NY,NX)+trcs_plant_uptake_vr(idg_N2O,L,NY,NX))*dts_gas
    RBGCSinkG(idg_NH3,L,NY,NX)=-TRN3G(L,NY,NX)*dts_gas
    RBGCSinkG(idg_H2,L,NY,NX)=(trcg_RMicbTransf_vr(idg_H2,L,NY,NX)+trcs_plant_uptake_vr(idg_H2,L,NY,NX))*dts_gas
    DO  K=1,jcplx
      RDOM_micb_cumflx(idom_doc,K,L,NY,NX)=-RDOM_micb_flx(idom_doc,K,L,NY,NX)*dts_HeatWatTP
      RDOM_micb_cumflx(idom_don,K,L,NY,NX)=-RDOM_micb_flx(idom_don,K,L,NY,NX)*dts_HeatWatTP
      RDOM_micb_cumflx(idom_dop,K,L,NY,NX)=-RDOM_micb_flx(idom_dop,K,L,NY,NX)*dts_HeatWatTP
      RDOM_micb_cumflx(idom_acetate,K,L,NY,NX)=-RDOM_micb_flx(idom_acetate,K,L,NY,NX)*dts_HeatWatTP
    ENDDO
    RBGCSinkS(ids_NH4,L,NY,NX)=(-RNutMicbTransf_vr(ids_NH4,L,NY,NX)-trcn_RChem_soil_vr(ids_NH4,L,NY,NX)+trcs_plant_uptake_vr(ids_NH4,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(idg_NH3,L,NY,NX)=(-TR_NH3_soil_vr(L,NY,NX)+trcs_plant_uptake_vr(idg_NH3,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(ids_NO3,L,NY,NX)=(-RNutMicbTransf_vr(ids_NO3,L,NY,NX)-trcn_RChem_soil_vr(ids_NO3,L,NY,NX)+trcs_plant_uptake_vr(ids_NO3,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(ids_NO2,L,NY,NX)=(-RNutMicbTransf_vr(ids_NO2,L,NY,NX)-trcn_RChem_soil_vr(ids_NO2,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(ids_H2PO4,L,NY,NX)=(-RNutMicbTransf_vr(ids_H2PO4,L,NY,NX)-trcn_RChem_soil_vr(ids_H2PO4,L,NY,NX)+trcs_plant_uptake_vr(ids_H2PO4,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(ids_H1PO4,L,NY,NX)=(-RNutMicbTransf_vr(ids_H1PO4,L,NY,NX)-trcn_RChem_soil_vr(ids_H1PO4,L,NY,NX)+trcs_plant_uptake_vr(ids_H1PO4,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(ids_NH4B,L,NY,NX)=(-RNutMicbTransf_vr(ids_NH4B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_NH4B,L,NY,NX)+trcs_plant_uptake_vr(ids_NH4B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(idg_NH3B,L,NY,NX)=(-trcn_RChem_band_soil_vr(idg_NH3B,L,NY,NX)+trcs_plant_uptake_vr(idg_NH3B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(ids_NO3B,L,NY,NX)=(-RNutMicbTransf_vr(ids_NO3B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_NO3B,L,NY,NX)+trcs_plant_uptake_vr(ids_NO3B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(ids_NO2B,L,NY,NX)=(-RNutMicbTransf_vr(ids_NO2B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_NO2B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(ids_H2PO4B,L,NY,NX)=(-RNutMicbTransf_vr(ids_H2PO4B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_H2PO4B,L,NY,NX)+trcs_plant_uptake_vr(ids_H2PO4B,L,NY,NX))*dts_HeatWatTP
    RBGCSinkS(ids_H1PO4B,L,NY,NX)=(-RNutMicbTransf_vr(ids_H1PO4B,L,NY,NX)-trcn_RChem_band_soil_vr(ids_H1PO4B,L,NY,NX)+trcs_plant_uptake_vr(ids_H1PO4B,L,NY,NX))*dts_HeatWatTP
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
    trcs_Irrig_vr(idg_CO2,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*CO2_irrig_conc(NY,NX)
    trcs_Irrig_vr(idg_CH4,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*CH4_irrig_conc(NY,NX)
    trcs_Irrig_vr(idg_O2,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*O2_irrig_conc(NY,NX)
    trcs_Irrig_vr(idg_N2,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*N2_irrig_conc(NY,NX)
    trcs_Irrig_vr(idg_N2O,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*N2O_irrig_conc(NY,NX)
    trcs_Irrig_vr(idg_H2,L,NY,NX)=0.0_r8

    trcs_Irrig_vr(ids_NH4,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*NH4_irrig_conc(I,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)*natomw
    trcs_Irrig_vr(idg_NH3,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*NH3_irrig_conc(I,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_NO3,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*NO3_irrig_conc(I,NY,NX)*trcs_VLN(ids_NO3,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_H1PO4,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*HPO4_irrig_conc(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)*patomw
    trcs_Irrig_vr(ids_H2PO4,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*H2PO4_irrig_conc(I,NY,NX)*trcs_VLN(ids_H1PO4,L,NY,NX)*patomw
    trcs_Irrig_vr(ids_NH4B,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*NH4_irrig_conc(I,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)*natomw
    trcs_Irrig_vr(idg_NH3B,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*NH3_irrig_conc(I,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_NO3B,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*NO3_irrig_conc(I,NY,NX)*trcs_VLN(ids_NO3B,L,NY,NX)*natomw
    trcs_Irrig_vr(ids_H1PO4B,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*HPO4_irrig_conc(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)*patomw
    trcs_Irrig_vr(ids_H2PO4B,L,NY,NX)=FWatIrrigate2MicP(L,NY,NX)*H2PO4_irrig_conc(I,NY,NX)*trcs_VLN(ids_H1PO4B,L,NY,NX)*patomw
!
!     SUB-HOURLY SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
!
!     R*FLZ,R*FBZ=subsurface solute flux in non-band,band
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations
!
    DO NTS=ids_beg,ids_end
      RFLZ_sol(NTS,L,NY,NX)=trcs_Irrig_vr(NTS,L,NY,NX)*dts_HeatWatTP
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
    OCSGL2(L,NY,NX)=OCSGL(L,NY,NX)*dts_HeatWatTP
    ONSGL2(L,NY,NX)=ONSGL(L,NY,NX)*dts_HeatWatTP
    OPSGL2(L,NY,NX)=OPSGL(L,NY,NX)*dts_HeatWatTP
    OASGL2(L,NY,NX)=OASGL(L,NY,NX)*dts_HeatWatTP

    DO NTS=ids_beg,ids_end
      SolDifcc(NTS,L,NY,NX)=SolDifc(NTS,L,NY,NX)*dts_HeatWatTP
    ENDDO

    DO NTG=idg_beg,idg_end
      GasDifcc(NTG,L,NY,NX)=GasDifc(NTG,L,NY,NX)*dts_gas
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
      trc_gasml2(NTG,L,NY,NX)=trc_gasml(NTG,L,NY,NX)
    ENDDO

    DO  K=1,jcplx
      DOM_MicP2(idom_doc,K,L,NY,NX)=DOM(idom_doc,K,L,NY,NX)-RDOM_micb_flx(idom_doc,K,L,NY,NX)
      DOM_MicP2(idom_don,K,L,NY,NX)=DOM(idom_don,K,L,NY,NX)-RDOM_micb_flx(idom_don,K,L,NY,NX)
      DOM_MicP2(idom_dop,K,L,NY,NX)=DOM(idom_dop,K,L,NY,NX)-RDOM_micb_flx(idom_dop,K,L,NY,NX)
      DOM_MicP2(idom_acetate,K,L,NY,NX)=DOM(idom_acetate,K,L,NY,NX)-RDOM_micb_flx(idom_acetate,K,L,NY,NX)

      DOM_MacP2(idom_doc,K,L,NY,NX)=DOM_Macp(idom_doc,K,L,NY,NX)
      DOM_MacP2(idom_don,K,L,NY,NX)=DOM_Macp(idom_don,K,L,NY,NX)
      DOM_MacP2(idom_dop,K,L,NY,NX)=DOM_Macp(idom_dop,K,L,NY,NX)
      DOM_MacP2(idom_acetate,K,L,NY,NX)=DOM_Macp(idom_acetate,K,L,NY,NX)
    enddo
    DO NTS=ids_beg,ids_end
      trc_solml2(NTS,L,NY,NX)=trc_solml(NTS,L,NY,NX)
      trc_soHml2(NTS,L,NY,NX)=trc_soHml(NTS,L,NY,NX)
    ENDDO

  enddo
  end subroutine ImportFluxFromOutsideModules
!

end module TranspNoSaltMod
