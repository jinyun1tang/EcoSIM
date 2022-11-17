module IngridTranspMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use TrnsfrsDataMod
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use GridDataType
  use SurfLitterDataType
  use SoilBGCDataType
  use SurfSoilDataType
  use SnowDataType
  use ChemTranspDataType
  use EcoSIMSolverPar
  use minimathmod, only : AZMAX1,AZMIN1
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=__FILE__

  real(r8) :: RFEFLS0,RHYFLS0,RCAFLS0,RMGFLS0,RNAFLS0,RKAFLS0
  real(r8) :: ROHFLS0,RSOFLS0,RCLFLS0,RC3FLS0,RHCFLS0,RAL1FS0
  real(r8) :: RAL2FS0,RAL3FS0,RAL4FS0,RALSFS0,RFE1FS0,RFE2FS0
  real(r8) :: RFE3FS0,RFE4FS0,RFESFS0,RCAOFS0,RCACFS0,RCAHFS0
  real(r8) :: RCASFS0,RMGOFS0,RMGCFS0,RMGHFS0,RMGSFS0,RNACFS0
  real(r8) :: RNASFS0,RKASFS0,RH0PFS0,RH3PFS0,RF1PFS0,RF2PFS0
  real(r8) :: RC0PFS0,RC1PFS0,RC2PFS0,RM1PFS0,RALFLS1,RFEFLS1
  real(r8) :: RHYFLS1,RCAFLS1,RMGFLS1,RNAFLS1,RKAFLS1,ROHFLS1
  real(r8) :: RSOFLS1,RCLFLS1,RC3FLS1,RHCFLS1,RAL1FS1,RAL2FS1
  real(r8) :: RAL3FS1,RAL4FS1,RALSFS1,RFE1FS1,RFE2FS1,RFE3FS1
  real(r8) :: RFE4FS1,RFESFS1,RCAOFS1,RCACFS1,RCAHFS1,RCASFS1
  real(r8) :: RMGOFS1,RMGCFS1,RMGHFS1,RMGSFS1,RNACFS1,RNASFS1
  real(r8) :: RKASFS1,RH0PFS1,RH3PFS1,RF1PFS1,RF2PFS1,RC0PFS1
  real(r8) :: RC1PFS1,RC2PFS1,RM1PFS1,RH0BFB1,RH3BFB1,RF1BFB1
  real(r8) :: RF2BFB1,RC0BFB1,RC1BFB1,RC2BFB1,RM1BFB1
  real(r8) :: RALFLS0

  real(r8) :: C0PO41,C3PO41,CFE1P1,CFE2P1,CCA0P1,CCA1P1,CCA2P1
  real(r8) :: CAL1,CFE1,CHY1,CCA1,CMG1,CNA1,CKA1,COH1,CSO41
  real(r8) :: CKAS2,C0PO42,C3PO42,CFE1P2,CFE2P2,CCA0P2,CCA1P2
  real(r8) :: CALS2,CFE12,CFE22,CFE32,CFE42,CFES2,CCAO2,CCAC2
  real(r8) :: CCL1,CCO31,CHCO31,CAL11,CAL21,CAL31,CAL41,CALS1
  real(r8) :: CFE11,CFE21,CFE31,CFE41,CFES1,CCAO1,CCAC1,CCAH1
  real(r8) :: CCAS1,CMGO1,CMGC1,CMGH1,CMGS1,CNAC1,CNAS1,CKAS1
  real(r8) :: CMG1P1,CAL2,CFE2,CHY2,CCA2,CMG2,CNA2,CKA2,COH2
  real(r8) :: CSO42,CCL2,CCO32,CHCO32,CAL12,CAL22,CAL32,CAL42
  real(r8) :: CCAH2,CCAS2,CMGO2,CMGC2,CMGH2,CMGS2,CNAC2,CNAS2
  real(r8) :: CCA2P2,CMG1P2,C0POB2,C3POB2,CF1PB2,CF2PB2,CC0PB2
  real(r8) :: CC1PB2,CC2PB2,CM1PB2,C0POB1,C3POB1,CF1PB1,CF2PB1
  real(r8) :: CC0PB1,CC1PB1,CC2PB1,CM1PB1
  real(r8) :: DFHC0B,DFHC1B,DFHC2B,DFHM1B,DLYR2,DFHAL,DFHFE
  real(r8) :: DFHHY,DFHCA,DFHMG,DFHNA,DFHKA,DFHOH,DFHSO,DFHCL
  real(r8) :: DFHC3,DFHHC,DFHAL1,DFHAL2,DFHAL3,DFHAL4,DFHALS
  real(r8) :: DFHFE1,DFHFE2,DFHFE3,DFHFE4,DFHFES,DFHCAO,DFHCAC
  real(r8) :: DFHCAH,DFHCAS,DFHMGO,DFHMGC,DFHMGH,DFHMGS,DFHNAC
  real(r8) :: DFHNAS,DFHKAS,DFHH0P,DFHH3P,DFHF1P,DFHF2P,DFHC0P
  real(r8) :: DFHC1P,DFHC2P,DFHM1P,DFHH0B,DFHH3B,DFHF1B,DFHF2B
  real(r8) :: DIFPO,DIFAL,DIFFE,DIFHY,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: DIFOH,DIFSO,DIFCL,DIFC3,DIFHC
  real(r8) :: DLYR0,DLYR1,DISPN

  real(r8) :: RFHAL,RFHFE
  real(r8) :: RFHHY,RFHCA,RFHMG,RFHNA,RFHKA,RFHOH,RFHSO,RFHCL
  real(r8) :: RFHC3,RFHHC,RFHAL1,RFHAL2,RFHAL3,RFHAL4,RFHALS
  real(r8) :: RFHFE1,RFHFE2,RFHFE3,RFHFE4,RFHFES,RFHCAO,RFHCAC
  real(r8) :: RFHCAH,RFHCAS,RFHMGO,RFHMGC,RFHMGH,RFHMGS,RFHNAC
  real(r8) :: RFHNAS,RFHKAS,RFHH0P,RFHH3P,RFHF1P,RFHF2P,RFHC0P
  real(r8) :: RFHC1P,RFHC2P,RFHM1P,RFHH0B,RFHH3B,RFHF1B,RFHF2B
  real(r8) :: RFHC0B,RFHC1B,RFHC2B,RFHM1B
  real(r8) :: TORT0,TORT1,TORTL

  real(r8), PARAMETER :: XFRS=0.05_r8
  public :: SaltModelSoluteHydroFlux
  contains


!------------------------------------------------------------------------------------------

  subroutine SaltModelSoluteHydroFlux(I,M,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,M,NHW,NHE,NVN,NVS

  integer :: NY,NX,N1,N2
  real(r8) :: FLWRM1
!     begin_execution

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     ADD SINKS FROM SOLUTE.F
      call SoluteSinksInSoil(NY,NX)
!
!     SOLUTE FLUXES FROM MELTING SNOWPACK TO
!     SOIL SURFACE FROM SNOWMELT IN 'WATSUB' AND
!     CONCENTRATIONS IN SNOWPACK
      call SoluteFluxInSnowpack(M,NY,NX)
!
!     CONVECTIVE SOLUTE EXCHANGE BETWEEN RESIDUE AND SOIL SURFACE
!
      FLWRM1=FLWRM(M,NY,NX)
!
!     FLWRM=litter-soil water flux from watsub.f

      IF(FLWRM1.GT.0.0_r8)THEN
        call Residue2TopsoilSoluteAdvExch(M,NY,NX,FLWRM1)
!
      ELSE
        call Topsoil2ResidueSoluteAdvExch(M,NY,NX,FLWRM1)
      ENDIF
!
!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN RESIDUE AND
!     SOIL SURFACE FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VOLT,DLYR,AREA=soil surface volume, thickness, area
!     VOLWM=micropore water-filled porosity from watsub.f
!
      IF((VOLT(0,NY,NX).GT.ZEROS(NY,NX).AND.VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
        .AND.(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN

        call TopsoilResidueSolutedifusExch(M,NY,NX,FLWRM1)
      ELSE
        call ZeroDifusTopsoilResidueExch
      ENDIF
!
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT FLUXES BETWEEN
!     ADJACENT GRID CELLS = CONVECTIVE + DIFFUSIVE FLUXES

      call TopsoilResidueFluxAdvPlusDifus(NY,NX)
!
      call AccumHourlyTopsoilReisdueFlux(NY,NX)
!
!     MACROPORE-MICROPORE SOLUTE EXCHANGE IN SOIL
!     SURFACE LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!     FINHM=macro-micropore water transfer from watsub.f
!     VOLWM,VOLWHM=micropore,macropore water volume
!     RFL*=convective macropore-micropore solute transfer
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *H2,*2=macropore,micropore solute content
!
!     MACROPORE TO MICROPORE TRANSFER
!
      IF(FINHM(M,NU(NY,NX),NY,NX).GT.0.0)THEN
        call MacToMicPoreSoluteAdvExchange(M,NY,NX)
!
!     MICROPORE TO MACROPORE TRANSFER
!
      ELSEIF(FINHM(M,NU(NY,NX),NY,NX).LT.0.0)THEN
        call MicToMacPoreSoluteAdvExchange(M,NY,NX)
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
      ELSE
        call ZeroMicMacPoreExchange
      ENDIF
!
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
      IF(VOLWHM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
        call MacMicPoreSoluteDifusExchange(M,NY,NX)
      ELSE
        call ZeroMacMicPoreSoluteDifusExch
      ENDIF
!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!      call MacMicPoreFluxAdvPlusDifus(NY,NX)
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
      call AccumHourlyMicMacPoreFlux(NY,NX)
!
!     SOLUTE TRANSPORT FROM WATER OVERLAND FLOW
!     IN 'WATSUB' AND FROM SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE LAYER
!
!     QRM=runoff from watsub.f
!     RQR*0=solute in runoff
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     VOLWM=litter water volume from watsub.f
!     *S2=litter solute content
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!     X*QRS=accumulated hourly solute in runoff
!
!
      N1=NX
      N2=NY
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
        call SoluteFluxBySurfaceOutflow(M,N1,N2)
      ELSE
        call ZeroSoluteFluxbySurfaceflow(N1,N2)
      ENDIF
!
      call UpdateSoluteInSurfNeighbors(M,N1,N2,NY,NX,NHW,NHE,NVN,NVS)
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS
      call UpdateSoluteInSubsurfNeighbors(M,NY,NX,NHW,NHE,NVN,NVS)

    ENDDO
  ENDDO
  end subroutine SaltModelSoluteHydroFlux

!------------------------------------------------------------------------------------------

  subroutine SoluteSinksInSoil(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX

  integer :: L
!     begin_execution
!
!     RZ*2=solute flux at time step for flux calculations
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  DO L=NU(NY,NX),NL(NY,NX)
    ZAL2(L,NY,NX)=ZAL2(L,NY,NX)-RZAL2(L,NY,NX)
    ZFE2(L,NY,NX)=ZFE2(L,NY,NX)-RZFE2(L,NY,NX)
    ZHY2(L,NY,NX)=ZHY2(L,NY,NX)-RZHY2(L,NY,NX)
    ZCA2(L,NY,NX)=ZCA2(L,NY,NX)-RZCA2(L,NY,NX)
    ZMG2(L,NY,NX)=ZMG2(L,NY,NX)-RZMG2(L,NY,NX)
    ZNA2(L,NY,NX)=ZNA2(L,NY,NX)-RZNA2(L,NY,NX)
    ZKA2(L,NY,NX)=ZKA2(L,NY,NX)-RZKA2(L,NY,NX)
    ZOH2(L,NY,NX)=ZOH2(L,NY,NX)-RZOH2(L,NY,NX)
    ZSO42(L,NY,NX)=ZSO42(L,NY,NX)-RZSO42(L,NY,NX)
    ZCL2(L,NY,NX)=ZCL2(L,NY,NX)-RZCL2(L,NY,NX)
    ZCO32(L,NY,NX)=ZCO32(L,NY,NX)-RZCO32(L,NY,NX)
    ZHCO32(L,NY,NX)=ZHCO32(L,NY,NX)-RZHCO32(L,NY,NX)
    ZAL12(L,NY,NX)=ZAL12(L,NY,NX)-RZAL12(L,NY,NX)
    ZAL22(L,NY,NX)=ZAL22(L,NY,NX)-RZAL22(L,NY,NX)
    ZAL32(L,NY,NX)=ZAL32(L,NY,NX)-RZAL32(L,NY,NX)
    ZAL42(L,NY,NX)=ZAL42(L,NY,NX)-RZAL42(L,NY,NX)
    ZALS2(L,NY,NX)=ZALS2(L,NY,NX)-RZALS2(L,NY,NX)
    ZFE12(L,NY,NX)=ZFE12(L,NY,NX)-RZFE12(L,NY,NX)
    ZFE22(L,NY,NX)=ZFE22(L,NY,NX)-RZFE22(L,NY,NX)
    ZFE32(L,NY,NX)=ZFE32(L,NY,NX)-RZFE32(L,NY,NX)
    ZFE42(L,NY,NX)=ZFE42(L,NY,NX)-RZFE42(L,NY,NX)
    ZFES2(L,NY,NX)=ZFES2(L,NY,NX)-RZFES2(L,NY,NX)
    ZCAO2(L,NY,NX)=ZCAO2(L,NY,NX)-RZCAO2(L,NY,NX)
    ZCAC2(L,NY,NX)=ZCAC2(L,NY,NX)-RZCAC2(L,NY,NX)
    ZCAH2(L,NY,NX)=ZCAH2(L,NY,NX)-RZCAH2(L,NY,NX)
    ZCAS2(L,NY,NX)=ZCAS2(L,NY,NX)-RZCAS2(L,NY,NX)
    ZMGO2(L,NY,NX)=ZMGO2(L,NY,NX)-RZMGO2(L,NY,NX)
    ZMGC2(L,NY,NX)=ZMGC2(L,NY,NX)-RZMGC2(L,NY,NX)
    ZMGH2(L,NY,NX)=ZMGH2(L,NY,NX)-RZMGH2(L,NY,NX)
    ZMGS2(L,NY,NX)=ZMGS2(L,NY,NX)-RZMGS2(L,NY,NX)
    ZNAC2(L,NY,NX)=ZNAC2(L,NY,NX)-RZNAC2(L,NY,NX)
    ZNAS2(L,NY,NX)=ZNAS2(L,NY,NX)-RZNAS2(L,NY,NX)
    ZKAS2(L,NY,NX)=ZKAS2(L,NY,NX)-RZKAS2(L,NY,NX)
    H0PO42(L,NY,NX)=H0PO42(L,NY,NX)-RH0PO42(L,NY,NX)
    H3PO42(L,NY,NX)=H3PO42(L,NY,NX)-RH3PO42(L,NY,NX)
    ZFE1P2(L,NY,NX)=ZFE1P2(L,NY,NX)-RZFE1P2(L,NY,NX)
    ZFE2P2(L,NY,NX)=ZFE2P2(L,NY,NX)-RZFE2P2(L,NY,NX)
    ZCA0P2(L,NY,NX)=ZCA0P2(L,NY,NX)-RZCA0P2(L,NY,NX)
    ZCA1P2(L,NY,NX)=ZCA1P2(L,NY,NX)-RZCA1P2(L,NY,NX)
    ZCA2P2(L,NY,NX)=ZCA2P2(L,NY,NX)-RZCA2P2(L,NY,NX)
    ZMG1P2(L,NY,NX)=ZMG1P2(L,NY,NX)-RZMG1P2(L,NY,NX)
    H0POB2(L,NY,NX)=H0POB2(L,NY,NX)-RH0POB2(L,NY,NX)
    H3POB2(L,NY,NX)=H3POB2(L,NY,NX)-RH3POB2(L,NY,NX)
    ZF1PB2(L,NY,NX)=ZF1PB2(L,NY,NX)-RZF1PB2(L,NY,NX)
    ZF2PB2(L,NY,NX)=ZF2PB2(L,NY,NX)-RZF2PB2(L,NY,NX)
    ZC0PB2(L,NY,NX)=ZC0PB2(L,NY,NX)-RZC0PB2(L,NY,NX)
    ZC1PB2(L,NY,NX)=ZC1PB2(L,NY,NX)-RZC1PB2(L,NY,NX)
    ZC2PB2(L,NY,NX)=ZC2PB2(L,NY,NX)-RZC2PB2(L,NY,NX)
    ZM1PB2(L,NY,NX)=ZM1PB2(L,NY,NX)-RZM1PB2(L,NY,NX)
  ENDDO
  end subroutine SoluteSinksInSoil
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxInSnowpack(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  integer :: L,ICHKL,L2
  real(r8) :: VFLWW,VFLWR,VFLWS,VFLWPO4,VFLWPOB

!     begin_execution
!
!     VHCPWM,VHCPWX=current,minimum volumetric heat capacity of snowpack
!     VOLWSL=snowpack water content
!     FLQWM=snowpack water flux
!     R*BLS=solute flux in snowpack
!     X*BLS=hourly solute flux in snowpack
!
  ICHKL=0
  DO L=1,JS
    IF(VHCPWM(M,L,NY,NX).GT.VHCPWX(NY,NX))THEN
      L2=MIN(JS,L+1)
      IF(L.LT.JS.AND.VHCPWM(M,L2,NY,NX).GT.VHCPWX(NY,NX))THEN
        IF(VOLWSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VFLWW=AZMAX1(AMIN1(1.0_r8,FLQWM(M,L2,NY,NX)/VOLWSL(L,NY,NX)))
        ELSE
          VFLWW=1.0_r8
        ENDIF
        RALBLS(L2,NY,NX)=trcs_solsml2(idsa_Al,L,NY,NX)*VFLWW
        RFEBLS(L2,NY,NX)=trcs_solsml2(idsa_Fe,L,NY,NX)*VFLWW
        RHYBLS(L2,NY,NX)=trcs_solsml2(idsa_Hp,L,NY,NX)*VFLWW
        RCABLS(L2,NY,NX)=trcs_solsml2(idsa_Ca,L,NY,NX)*VFLWW
        RMGBLS(L2,NY,NX)=trcs_solsml2(idsa_Mg,L,NY,NX)*VFLWW
        RNABLS(L2,NY,NX)=trcs_solsml2(idsa_Na,L,NY,NX)*VFLWW
        RKABLS(L2,NY,NX)=trcs_solsml2(idsa_K,L,NY,NX)*VFLWW
        ROHBLS(L2,NY,NX)=trcs_solsml2(idsa_OH,L,NY,NX)*VFLWW
        RSOBLS(L2,NY,NX)=trcs_solsml2(idsa_SO4,L,NY,NX)*VFLWW
        RCLBLS(L2,NY,NX)=trcs_solsml2(idsa_Cl,L,NY,NX)*VFLWW
        RC3BLS(L2,NY,NX)=trcs_solsml2(idsa_CO3,L,NY,NX)*VFLWW
        RHCBLS(L2,NY,NX)=trcs_solsml2(idsa_HCO3,L,NY,NX)*VFLWW
        RAL1BS(L2,NY,NX)=trcs_solsml2(idsa_AlOH,L,NY,NX)*VFLWW
        RAL2BS(L2,NY,NX)=trcs_solsml2(idsa_AlOH2,L,NY,NX)*VFLWW
        RAL3BS(L2,NY,NX)=trcs_solsml2(idsa_AlOH3,L,NY,NX)*VFLWW
        RAL4BS(L2,NY,NX)=trcs_solsml2(idsa_AlOH4,L,NY,NX)*VFLWW
        RALSBS(L2,NY,NX)=trcs_solsml2(idsa_AlSO4,L,NY,NX)*VFLWW
        RFE1BS(L2,NY,NX)=trcs_solsml2(idsa_FeOH,L,NY,NX)*VFLWW
        RFE2BS(L2,NY,NX)=trcs_solsml2(idsa_FeOH2,L,NY,NX)*VFLWW
        RFE3BS(L2,NY,NX)=trcs_solsml2(idsa_FeOH3,L,NY,NX)*VFLWW
        RFE4BS(L2,NY,NX)=trcs_solsml2(idsa_FeOH4,L,NY,NX)*VFLWW
        RFESBS(L2,NY,NX)=trcs_solsml2(idsa_FeSO4,L,NY,NX)*VFLWW
        RCAOBS(L2,NY,NX)=trcs_solsml2(idsa_CaOH2,L,NY,NX)*VFLWW
        RCACBS(L2,NY,NX)=trcs_solsml2(idsa_CaCO3,L,NY,NX)*VFLWW
        RCAHBS(L2,NY,NX)=trcs_solsml2(idsa_CaHCO3,L,NY,NX)*VFLWW
        RCASBS(L2,NY,NX)=trcs_solsml2(idsa_CaSO4,L,NY,NX)*VFLWW
        RMGOBS(L2,NY,NX)=trcs_solsml2(idsa_MgOH2,L,NY,NX)*VFLWW
        RMGCBS(L2,NY,NX)=trcs_solsml2(idsa_MgCO3,L,NY,NX)*VFLWW
        RMGHBS(L2,NY,NX)=trcs_solsml2(idsa_MgHCO3,L,NY,NX)*VFLWW
        RMGSBS(L2,NY,NX)=trcs_solsml2(idsa_MgSO4,L,NY,NX)*VFLWW
        RNACBS(L2,NY,NX)=trcs_solsml2(idsa_NaCO3,L,NY,NX)*VFLWW
        RNASBS(L2,NY,NX)=trcs_solsml2(idsa_NaSO4,L,NY,NX)*VFLWW
        RKASBS(L2,NY,NX)=trcs_solsml2(idsa_KSO4,L,NY,NX)*VFLWW
        RH0PBS(L2,NY,NX)=trcs_solsml2(idsa_H0PO4,L,NY,NX)*VFLWW
        RH3PBS(L2,NY,NX)=trcs_solsml2(idsa_H3PO4,L,NY,NX)*VFLWW
        RF1PBS(L2,NY,NX)=trcs_solsml2(idsa_FeHPO4,L,NY,NX)*VFLWW
        RF2PBS(L2,NY,NX)=trcs_solsml2(idsa_FeH2PO4,L,NY,NX)*VFLWW
        RC0PBS(L2,NY,NX)=trcs_solsml2(idsa_CaPO4,L,NY,NX)*VFLWW
        RC1PBS(L2,NY,NX)=trcs_solsml2(idsa_CaHPO4,L,NY,NX)*VFLWW
        RC2PBS(L2,NY,NX)=trcs_solsml2(idsa_CaH2PO4,L,NY,NX)*VFLWW
        RM1PBS(L2,NY,NX)=trcs_solsml2(idsa_MgHPO4,L,NY,NX)*VFLWW

        trcsa_XBLS(idsa_Al,L2,NY,NX)=trcsa_XBLS(idsa_Al,L2,NY,NX)+RALBLS(L2,NY,NX)
        trcsa_XBLS(idsa_Fe,L2,NY,NX)=trcsa_XBLS(idsa_Fe,L2,NY,NX)+RFEBLS(L2,NY,NX)
        trcsa_XBLS(idsa_Hp,L2,NY,NX)=trcsa_XBLS(idsa_Hp,L2,NY,NX)+RHYBLS(L2,NY,NX)
        trcsa_XBLS(idsa_Ca,L2,NY,NX)=trcsa_XBLS(idsa_Ca,L2,NY,NX)+RCABLS(L2,NY,NX)
        trcsa_XBLS(idsa_Mg,L2,NY,NX)=trcsa_XBLS(idsa_Mg,L2,NY,NX)+RMGBLS(L2,NY,NX)
        trcsa_XBLS(idsa_Na,L2,NY,NX)=trcsa_XBLS(idsa_Na,L2,NY,NX)+RNABLS(L2,NY,NX)
        trcsa_XBLS(idsa_K,L2,NY,NX)=trcsa_XBLS(idsa_K,L2,NY,NX)+RKABLS(L2,NY,NX)
        trcsa_XBLS(idsa_OH,L2,NY,NX)=trcsa_XBLS(idsa_OH,L2,NY,NX)+ROHBLS(L2,NY,NX)
        trcsa_XBLS(idsa_SO4,L2,NY,NX)=trcsa_XBLS(idsa_SO4,L2,NY,NX)+RSOBLS(L2,NY,NX)
        trcsa_XBLS(idsa_Cl,L2,NY,NX)=trcsa_XBLS(idsa_Cl,L2,NY,NX)+RCLBLS(L2,NY,NX)
        trcsa_XBLS(idsa_CO3,L2,NY,NX)=trcsa_XBLS(idsa_CO3,L2,NY,NX)+RC3BLS(L2,NY,NX)
        trcsa_XBLS(idsa_HCO3,L2,NY,NX)=trcsa_XBLS(idsa_HCO3,L2,NY,NX)+RHCBLS(L2,NY,NX)
        trcsa_XBLS(idsa_AlOH,L2,NY,NX)=trcsa_XBLS(idsa_AlOH,L2,NY,NX)+RAL1BS(L2,NY,NX)
        trcsa_XBLS(idsa_AlOH2,L2,NY,NX)=trcsa_XBLS(idsa_AlOH2,L2,NY,NX)+RAL2BS(L2,NY,NX)
        trcsa_XBLS(idsa_AlOH3,L2,NY,NX)=trcsa_XBLS(idsa_AlOH3,L2,NY,NX)+RAL3BS(L2,NY,NX)
        trcsa_XBLS(idsa_AlOH4,L2,NY,NX)=trcsa_XBLS(idsa_AlOH4,L2,NY,NX)+RAL4BS(L2,NY,NX)
        trcsa_XBLS(idsa_AlSO4,L2,NY,NX)=trcsa_XBLS(idsa_AlSO4,L2,NY,NX)+RALSBS(L2,NY,NX)
        trcsa_XBLS(idsa_FeOH,L2,NY,NX)=trcsa_XBLS(idsa_FeOH,L2,NY,NX)+RFE1BS(L2,NY,NX)
        trcsa_XBLS(idsa_FeOH2,L2,NY,NX)=trcsa_XBLS(idsa_FeOH2,L2,NY,NX)+RFE2BS(L2,NY,NX)
        trcsa_XBLS(idsa_FeOH3,L2,NY,NX)=trcsa_XBLS(idsa_FeOH3,L2,NY,NX)+RFE3BS(L2,NY,NX)
        trcsa_XBLS(idsa_FeOH4,L2,NY,NX)=trcsa_XBLS(idsa_FeOH4,L2,NY,NX)+RFE4BS(L2,NY,NX)
        trcsa_XBLS(idsa_FeSO4,L2,NY,NX)=trcsa_XBLS(idsa_FeSO4,L2,NY,NX)+RFESBS(L2,NY,NX)
        trcsa_XBLS(idsa_CaOH2,L2,NY,NX)=trcsa_XBLS(idsa_CaOH2,L2,NY,NX)+RCAOBS(L2,NY,NX)
        trcsa_XBLS(idsa_CaCO3,L2,NY,NX)=trcsa_XBLS(idsa_CaCO3,L2,NY,NX)+RCACBS(L2,NY,NX)
        trcsa_XBLS(idsa_CaHCO3,L2,NY,NX)=trcsa_XBLS(idsa_CaHCO3,L2,NY,NX)+RCAHBS(L2,NY,NX)
        trcsa_XBLS(idsa_CaSO4,L2,NY,NX)=trcsa_XBLS(idsa_CaSO4,L2,NY,NX)+RCASBS(L2,NY,NX)
        trcsa_XBLS(idsa_MgOH2,L2,NY,NX)=trcsa_XBLS(idsa_MgOH2,L2,NY,NX)+RMGOBS(L2,NY,NX)
        trcsa_XBLS(idsa_MgCO3,L2,NY,NX)=trcsa_XBLS(idsa_MgCO3,L2,NY,NX)+RMGCBS(L2,NY,NX)
        trcsa_XBLS(idsa_MgHCO3,L2,NY,NX)=trcsa_XBLS(idsa_MgHCO3,L2,NY,NX)+RMGHBS(L2,NY,NX)
        trcsa_XBLS(idsa_MgSO4,L2,NY,NX)=trcsa_XBLS(idsa_MgSO4,L2,NY,NX)+RMGSBS(L2,NY,NX)
        trcsa_XBLS(idsa_NaCO3,L2,NY,NX)=trcsa_XBLS(idsa_NaCO3,L2,NY,NX)+RNACBS(L2,NY,NX)
        trcsa_XBLS(idsa_NaSO4,L2,NY,NX)=trcsa_XBLS(idsa_NaSO4,L2,NY,NX)+RNASBS(L2,NY,NX)
        trcsa_XBLS(idsa_KSO4,L2,NY,NX)=trcsa_XBLS(idsa_KSO4,L2,NY,NX)+RKASBS(L2,NY,NX)
        trcsa_XBLS(idsa_H0PO4,L2,NY,NX)=trcsa_XBLS(idsa_H0PO4,L2,NY,NX)+RH0PBS(L2,NY,NX)
        trcsa_XBLS(idsa_H3PO4,L2,NY,NX)=trcsa_XBLS(idsa_H3PO4,L2,NY,NX)+RH3PBS(L2,NY,NX)
        trcsa_XBLS(idsa_FeHPO4,L2,NY,NX)=trcsa_XBLS(idsa_FeHPO4,L2,NY,NX)+RF1PBS(L2,NY,NX)
        trcsa_XBLS(idsa_FeH2PO4,L2,NY,NX)=trcsa_XBLS(idsa_FeH2PO4,L2,NY,NX)+RF2PBS(L2,NY,NX)
        trcsa_XBLS(idsa_CaPO4,L2,NY,NX)=trcsa_XBLS(idsa_CaPO4,L2,NY,NX)+RC0PBS(L2,NY,NX)
        trcsa_XBLS(idsa_CaHPO4,L2,NY,NX)=trcsa_XBLS(idsa_CaHPO4,L2,NY,NX)+RC1PBS(L2,NY,NX)
        trcsa_XBLS(idsa_CaH2PO4,L2,NY,NX)=trcsa_XBLS(idsa_CaH2PO4,L2,NY,NX)+RC2PBS(L2,NY,NX)
        trcsa_XBLS(idsa_MgHPO4,L2,NY,NX)=trcsa_XBLS(idsa_MgHPO4,L2,NY,NX)+RM1PBS(L2,NY,NX)
      ELSE
        IF(L.LT.JS)THEN
          RALBLS(L2,NY,NX)=0.0_r8
          RFEBLS(L2,NY,NX)=0.0_r8
          RHYBLS(L2,NY,NX)=0.0_r8
          RCABLS(L2,NY,NX)=0.0_r8
          RMGBLS(L2,NY,NX)=0.0_r8
          RNABLS(L2,NY,NX)=0.0_r8
          RKABLS(L2,NY,NX)=0.0_r8
          ROHBLS(L2,NY,NX)=0.0_r8
          RSOBLS(L2,NY,NX)=0.0_r8
          RCLBLS(L2,NY,NX)=0.0_r8
          RC3BLS(L2,NY,NX)=0.0_r8
          RHCBLS(L2,NY,NX)=0.0_r8
          RAL1BS(L2,NY,NX)=0.0_r8
          RAL2BS(L2,NY,NX)=0.0_r8
          RAL3BS(L2,NY,NX)=0.0_r8
          RAL4BS(L2,NY,NX)=0.0_r8
          RALSBS(L2,NY,NX)=0.0_r8
          RFE1BS(L2,NY,NX)=0.0_r8
          RFE2BS(L2,NY,NX)=0.0_r8
          RFE3BS(L2,NY,NX)=0.0_r8
          RFE4BS(L2,NY,NX)=0.0_r8
          RFESBS(L2,NY,NX)=0.0_r8
          RCAOBS(L2,NY,NX)=0.0_r8
          RCACBS(L2,NY,NX)=0.0_r8
          RCAHBS(L2,NY,NX)=0.0_r8
          RCASBS(L2,NY,NX)=0.0_r8
          RMGOBS(L2,NY,NX)=0.0_r8
          RMGCBS(L2,NY,NX)=0.0_r8
          RMGHBS(L2,NY,NX)=0.0_r8
          RMGSBS(L2,NY,NX)=0.0_r8
          RNACBS(L2,NY,NX)=0.0_r8
          RNASBS(L2,NY,NX)=0.0_r8
          RKASBS(L2,NY,NX)=0.0_r8
          RH0PBS(L2,NY,NX)=0.0_r8
          RH3PBS(L2,NY,NX)=0.0_r8
          RF1PBS(L2,NY,NX)=0.0_r8
          RF2PBS(L2,NY,NX)=0.0_r8
          RC0PBS(L2,NY,NX)=0.0_r8
          RC1PBS(L2,NY,NX)=0.0_r8
          RC2PBS(L2,NY,NX)=0.0_r8
          RM1PBS(L2,NY,NX)=0.0_r8
        ENDIF
!
!     SNOWPACK SOLUTE DISCHARGE TO SURFACE LITTER, SOIL SURFACE
!
!     VOLWSL=snowpack water content
!     FLQRM,FLQSM,FLQHM=total water flux to litter,soil micropore,macropore
!     CVRD,BARE=litter cover fraction,1-CVRD
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     R*S0=solute flux to surface litter
!     R*S1,R*B1=solute flux to soil surface non-band,band
!
        IF(ICHKL.EQ.0)THEN
          IF(VOLWSL(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            VFLWR=AZMAX1(AMIN1(1.0_r8,FLQRM(M,NY,NX)/VOLWSL(L,NY,NX)))
            VFLWS=AZMAX1(AMIN1(1.0_r8,(FLQSM(M,NY,NX)+FLQHM(M,NY,NX))/VOLWSL(L,NY,NX)))
          ELSE
            VFLWR=CVRD(NY,NX)
            VFLWS=BARE(NY,NX)
          ENDIF
          VFLWPO4=VFLWS*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
          VFLWPOB=VFLWS*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
          RALFLS0=trcs_solsml2(idsa_Al,L,NY,NX)*VFLWR
          RFEFLS0=trcs_solsml2(idsa_Fe,L,NY,NX)*VFLWR
          RHYFLS0=trcs_solsml2(idsa_Hp,L,NY,NX)*VFLWR
          RCAFLS0=trcs_solsml2(idsa_Ca,L,NY,NX)*VFLWR
          RMGFLS0=trcs_solsml2(idsa_Mg,L,NY,NX)*VFLWR
          RNAFLS0=trcs_solsml2(idsa_Na,L,NY,NX)*VFLWR
          RKAFLS0=trcs_solsml2(idsa_K,L,NY,NX)*VFLWR
          ROHFLS0=trcs_solsml2(idsa_OH,L,NY,NX)*VFLWR
          RSOFLS0=trcs_solsml2(idsa_SO4,L,NY,NX)*VFLWR
          RCLFLS0=trcs_solsml2(idsa_Cl,L,NY,NX)*VFLWR
          RC3FLS0=trcs_solsml2(idsa_CO3,L,NY,NX)*VFLWR
          RHCFLS0=trcs_solsml2(idsa_HCO3,L,NY,NX)*VFLWR
          RAL1FS0=trcs_solsml2(idsa_AlOH,L,NY,NX)*VFLWR
          RAL2FS0=trcs_solsml2(idsa_AlOH2,L,NY,NX)*VFLWR
          RAL3FS0=trcs_solsml2(idsa_AlOH3,L,NY,NX)*VFLWR
          RAL4FS0=trcs_solsml2(idsa_AlOH4,L,NY,NX)*VFLWR
          RALSFS0=trcs_solsml2(idsa_AlSO4,L,NY,NX)*VFLWR
          RFE1FS0=trcs_solsml2(idsa_FeOH,L,NY,NX)*VFLWR
          RFE2FS0=trcs_solsml2(idsa_FeOH2,L,NY,NX)*VFLWR
          RFE3FS0=trcs_solsml2(idsa_FeOH3,L,NY,NX)*VFLWR
          RFE4FS0=trcs_solsml2(idsa_FeOH4,L,NY,NX)*VFLWR
          RFESFS0=trcs_solsml2(idsa_FeSO4,L,NY,NX)*VFLWR
          RCAOFS0=trcs_solsml2(idsa_CaOH2,L,NY,NX)*VFLWR
          RCACFS0=trcs_solsml2(idsa_CaCO3,L,NY,NX)*VFLWR
          RCAHFS0=trcs_solsml2(idsa_CaHCO3,L,NY,NX)*VFLWR
          RCASFS0=trcs_solsml2(idsa_CaSO4,L,NY,NX)*VFLWR
          RMGOFS0=trcs_solsml2(idsa_MgOH2,L,NY,NX)*VFLWR
          RMGCFS0=trcs_solsml2(idsa_MgCO3,L,NY,NX)*VFLWR
          RMGHFS0=trcs_solsml2(idsa_MgHCO3,L,NY,NX)*VFLWR
          RMGSFS0=trcs_solsml2(idsa_MgSO4,L,NY,NX)*VFLWR
          RNACFS0=trcs_solsml2(idsa_NaCO3,L,NY,NX)*VFLWR
          RNASFS0=trcs_solsml2(idsa_NaSO4,L,NY,NX)*VFLWR
          RKASFS0=trcs_solsml2(idsa_KSO4,L,NY,NX)*VFLWR
          RH0PFS0=trcs_solsml2(idsa_H0PO4,L,NY,NX)*VFLWR
          RH3PFS0=trcs_solsml2(idsa_H3PO4,L,NY,NX)*VFLWR
          RF1PFS0=trcs_solsml2(idsa_FeHPO4,L,NY,NX)*VFLWR
          RF2PFS0=trcs_solsml2(idsa_FeH2PO4,L,NY,NX)*VFLWR
          RC0PFS0=trcs_solsml2(idsa_CaPO4,L,NY,NX)*VFLWR
          RC1PFS0=trcs_solsml2(idsa_CaHPO4,L,NY,NX)*VFLWR
          RC2PFS0=trcs_solsml2(idsa_CaH2PO4,L,NY,NX)*VFLWR
          RM1PFS0=trcs_solsml2(idsa_MgHPO4,L,NY,NX)*VFLWR
          RALFLS1=trcs_solsml2(idsa_Al,L,NY,NX)*VFLWS
          RFEFLS1=trcs_solsml2(idsa_Fe,L,NY,NX)*VFLWS
          RHYFLS1=trcs_solsml2(idsa_Hp,L,NY,NX)*VFLWS
          RCAFLS1=trcs_solsml2(idsa_Ca,L,NY,NX)*VFLWS
          RMGFLS1=trcs_solsml2(idsa_Mg,L,NY,NX)*VFLWS
          RNAFLS1=trcs_solsml2(idsa_Na,L,NY,NX)*VFLWS
          RKAFLS1=trcs_solsml2(idsa_K,L,NY,NX)*VFLWS
          ROHFLS1=trcs_solsml2(idsa_OH,L,NY,NX)*VFLWS
          RSOFLS1=trcs_solsml2(idsa_SO4,L,NY,NX)*VFLWS
          RCLFLS1=trcs_solsml2(idsa_Cl,L,NY,NX)*VFLWS
          RC3FLS1=trcs_solsml2(idsa_CO3,L,NY,NX)*VFLWS
          RHCFLS1=trcs_solsml2(idsa_HCO3,L,NY,NX)*VFLWS
          RAL1FS1=trcs_solsml2(idsa_AlOH,L,NY,NX)*VFLWS
          RAL2FS1=trcs_solsml2(idsa_AlOH2,L,NY,NX)*VFLWS
          RAL3FS1=trcs_solsml2(idsa_AlOH3,L,NY,NX)*VFLWS
          RAL4FS1=trcs_solsml2(idsa_AlOH4,L,NY,NX)*VFLWS
          RALSFS1=trcs_solsml2(idsa_AlSO4,L,NY,NX)*VFLWS
          RFE1FS1=trcs_solsml2(idsa_FeOH,L,NY,NX)*VFLWS
          RFE2FS1=trcs_solsml2(idsa_FeOH2,L,NY,NX)*VFLWS
          RFE3FS1=trcs_solsml2(idsa_FeOH3,L,NY,NX)*VFLWS
          RFE4FS1=trcs_solsml2(idsa_FeOH4,L,NY,NX)*VFLWS
          RFESFS1=trcs_solsml2(idsa_FeSO4,L,NY,NX)*VFLWS
          RCAOFS1=trcs_solsml2(idsa_CaOH2,L,NY,NX)*VFLWS
          RCACFS1=trcs_solsml2(idsa_CaCO3,L,NY,NX)*VFLWS
          RCAHFS1=trcs_solsml2(idsa_CaHCO3,L,NY,NX)*VFLWS
          RCASFS1=trcs_solsml2(idsa_CaSO4,L,NY,NX)*VFLWS
          RMGOFS1=trcs_solsml2(idsa_MgOH2,L,NY,NX)*VFLWS
          RMGCFS1=trcs_solsml2(idsa_MgCO3,L,NY,NX)*VFLWS
          RMGHFS1=trcs_solsml2(idsa_MgHCO3,L,NY,NX)*VFLWS
          RMGSFS1=trcs_solsml2(idsa_MgSO4,L,NY,NX)*VFLWS
          RNACFS1=trcs_solsml2(idsa_NaCO3,L,NY,NX)*VFLWS
          RNASFS1=trcs_solsml2(idsa_NaSO4,L,NY,NX)*VFLWS
          RKASFS1=trcs_solsml2(idsa_KSO4,L,NY,NX)*VFLWS
          RH0PFS1=trcs_solsml2(idsa_H0PO4,L,NY,NX)*VFLWPO4
          RH3PFS1=trcs_solsml2(idsa_H3PO4,L,NY,NX)*VFLWPO4
          RF1PFS1=trcs_solsml2(idsa_FeHPO4,L,NY,NX)*VFLWPO4
          RF2PFS1=trcs_solsml2(idsa_FeH2PO4,L,NY,NX)*VFLWPO4
          RC0PFS1=trcs_solsml2(idsa_CaPO4,L,NY,NX)*VFLWPO4
          RC1PFS1=trcs_solsml2(idsa_CaHPO4,L,NY,NX)*VFLWPO4
          RC2PFS1=trcs_solsml2(idsa_CaH2PO4,L,NY,NX)*VFLWPO4
          RM1PFS1=trcs_solsml2(idsa_MgHPO4,L,NY,NX)*VFLWPO4
          RH0BFB1=trcs_solsml2(idsa_H0PO4,L,NY,NX)*VFLWPOB
          RH3BFB1=trcs_solsml2(idsa_H3PO4,L,NY,NX)*VFLWPOB
          RF1BFB1=trcs_solsml2(idsa_FeHPO4,L,NY,NX)*VFLWPOB
          RF2BFB1=trcs_solsml2(idsa_FeH2PO4,L,NY,NX)*VFLWPOB
          RC0BFB1=trcs_solsml2(idsa_CaPO4,L,NY,NX)*VFLWPOB
          RC1BFB1=trcs_solsml2(idsa_CaHPO4,L,NY,NX)*VFLWPOB
          RC2BFB1=trcs_solsml2(idsa_CaH2PO4,L,NY,NX)*VFLWPOB
          RM1BFB1=trcs_solsml2(idsa_MgHPO4,L,NY,NX)*VFLWPOB
          ICHKL=1
        ENDIF
      ENDIF
    ELSE
      RALFLS0=0.0_r8
      RFEFLS0=0.0_r8
      RHYFLS0=0.0_r8
      RCAFLS0=0.0_r8
      RMGFLS0=0.0_r8
      RNAFLS0=0.0_r8
      RKAFLS0=0.0_r8
      ROHFLS0=0.0_r8
      RSOFLS0=0.0_r8
      RCLFLS0=0.0_r8
      RC3FLS0=0.0_r8
      RHCFLS0=0.0_r8
      RAL1FS0=0.0_r8
      RAL2FS0=0.0_r8
      RAL3FS0=0.0_r8
      RAL4FS0=0.0_r8
      RALSFS0=0.0_r8
      RFE1FS0=0.0_r8
      RFE2FS0=0.0_r8
      RFE3FS0=0.0_r8
      RFE4FS0=0.0_r8
      RFESFS0=0.0_r8
      RCAOFS0=0.0_r8
      RCACFS0=0.0_r8
      RCAHFS0=0.0_r8
      RCASFS0=0.0_r8
      RMGOFS0=0.0_r8
      RMGCFS0=0.0_r8
      RMGHFS0=0.0_r8
      RMGSFS0=0.0_r8
      RNACFS0=0.0_r8
      RNASFS0=0.0_r8
      RKASFS0=0.0_r8
      RH0PFS0=0.0_r8
      RH3PFS0=0.0_r8
      RF1PFS0=0.0_r8
      RF2PFS0=0.0_r8
      RC0PFS0=0.0_r8
      RC1PFS0=0.0_r8
      RC2PFS0=0.0_r8
      RM1PFS0=0.0_r8
      RALFLS1=0.0_r8
      RFEFLS1=0.0_r8
      RHYFLS1=0.0_r8
      RCAFLS1=0.0_r8
      RMGFLS1=0.0_r8
      RNAFLS1=0.0_r8
      RKAFLS1=0.0_r8
      ROHFLS1=0.0_r8
      RSOFLS1=0.0_r8
      RCLFLS1=0.0_r8
      RC3FLS1=0.0_r8
      RHCFLS1=0.0_r8
      RAL1FS1=0.0_r8
      RAL2FS1=0.0_r8
      RAL3FS1=0.0_r8
      RAL4FS1=0.0_r8
      RALSFS1=0.0_r8
      RFE1FS1=0.0_r8
      RFE2FS1=0.0_r8
      RFE3FS1=0.0_r8
      RFE4FS1=0.0_r8
      RFESFS1=0.0_r8
      RCAOFS1=0.0_r8
      RCACFS1=0.0_r8
      RCAHFS1=0.0_r8
      RCASFS1=0.0_r8
      RMGOFS1=0.0_r8
      RMGCFS1=0.0_r8
      RMGHFS1=0.0_r8
      RMGSFS1=0.0_r8
      RNACFS1=0.0_r8
      RNASFS1=0.0_r8
      RKASFS1=0.0_r8
      RH0PFS1=0.0_r8
      RH3PFS1=0.0_r8
      RF1PFS1=0.0_r8
      RF2PFS1=0.0_r8
      RC0PFS1=0.0_r8
      RC1PFS1=0.0_r8
      RC2PFS1=0.0_r8
      RM1PFS1=0.0_r8
      RH0BFB1=0.0_r8
      RH3BFB1=0.0_r8
      RF1BFB1=0.0_r8
      RF2BFB1=0.0_r8
      RC0BFB1=0.0_r8
      RC1BFB1=0.0_r8
      RC2BFB1=0.0_r8
      RM1BFB1=0.0_r8
    ENDIF
  ENDDO
  end subroutine SoluteFluxInSnowpack
!------------------------------------------------------------------------------------------

  subroutine Residue2TopsoilSoluteAdvExch(M,NY,NX,FLWRM1)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8) :: VFLW

!     begin_execution
!
!     IF WATER FLUX FROM 'WATSUB' IS FROM RESIDUE TO
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN RESIDUE
!
!     VOLWM=litter water volume
!     RFL*=soil-litter convective solute flux
!     Z*2=litter solute content
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
  IF(VOLWM(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMAX1(AMIN1(VFLWX,FLWRM1/VOLWM(M,0,NY,NX)))
  ELSE
    VFLW=VFLWX
  ENDIF
  RFLAL=VFLW*AZMAX1(ZAL2(0,NY,NX))
  RFLFE=VFLW*AZMAX1(ZFE2(0,NY,NX))
  RFLHY=VFLW*AZMAX1(ZHY2(0,NY,NX))
  RFLCA=VFLW*AZMAX1(ZCA2(0,NY,NX))
  RFLMG=VFLW*AZMAX1(ZMG2(0,NY,NX))
  RFLNA=VFLW*AZMAX1(ZNA2(0,NY,NX))
  RFLKA=VFLW*AZMAX1(ZKA2(0,NY,NX))
  RFLOH=VFLW*AZMAX1(ZOH2(0,NY,NX))
  RFLSO=VFLW*AZMAX1(ZSO42(0,NY,NX))
  RFLCL=VFLW*AZMAX1(ZCL2(0,NY,NX))
  RFLC3=VFLW*AZMAX1(ZCO32(0,NY,NX))
  RFLHC=VFLW*AZMAX1(ZHCO32(0,NY,NX))
  RFLAL1=VFLW*AZMAX1(ZAL12(0,NY,NX))
  RFLAL2=VFLW*AZMAX1(ZAL22(0,NY,NX))
  RFLAL3=VFLW*AZMAX1(ZAL32(0,NY,NX))
  RFLAL4=VFLW*AZMAX1(ZAL42(0,NY,NX))
  RFLALS=VFLW*AZMAX1(ZALS2(0,NY,NX))
  RFLFE1=VFLW*AZMAX1(ZFE12(0,NY,NX))
  RFLFE2=VFLW*AZMAX1(ZFE22(0,NY,NX))
  RFLFE3=VFLW*AZMAX1(ZFE32(0,NY,NX))
  RFLFE4=VFLW*AZMAX1(ZFE42(0,NY,NX))
  RFLFES=VFLW*AZMAX1(ZFES2(0,NY,NX))
  RFLCAO=VFLW*AZMAX1(ZCAO2(0,NY,NX))
  RFLCAC=VFLW*AZMAX1(ZCAC2(0,NY,NX))
  RFLCAH=VFLW*AZMAX1(ZCAH2(0,NY,NX))
  RFLCAS=VFLW*AZMAX1(ZCAS2(0,NY,NX))
  RFLMGO=VFLW*AZMAX1(ZMGO2(0,NY,NX))
  RFLMGC=VFLW*AZMAX1(ZMGC2(0,NY,NX))
  RFLMGH=VFLW*AZMAX1(ZMGH2(0,NY,NX))
  RFLMGS=VFLW*AZMAX1(ZMGS2(0,NY,NX))
  RFLNAC=VFLW*AZMAX1(ZNAC2(0,NY,NX))
  RFLNAS=VFLW*AZMAX1(ZNAS2(0,NY,NX))
  RFLKAS=VFLW*AZMAX1(ZKAS2(0,NY,NX))
  RFLH0P=VFLW*AZMAX1(H0PO42(0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH3P=VFLW*AZMAX1(H3PO42(0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF1P=VFLW*AZMAX1(ZFE1P2(0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF2P=VFLW*AZMAX1(ZFE2P2(0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC0P=VFLW*AZMAX1(ZCA0P2(0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC1P=VFLW*AZMAX1(ZCA1P2(0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC2P=VFLW*AZMAX1(ZCA2P2(0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLM1P=VFLW*AZMAX1(ZMG1P2(0,NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH0B=VFLW*AZMAX1(H0PO42(0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLH3B=VFLW*AZMAX1(H3PO42(0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF1B=VFLW*AZMAX1(ZFE1P2(0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF2B=VFLW*AZMAX1(ZFE2P2(0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC0B=VFLW*AZMAX1(ZCA0P2(0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC1B=VFLW*AZMAX1(ZCA1P2(0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC2B=VFLW*AZMAX1(ZCA2P2(0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLM1B=VFLW*AZMAX1(ZMG1P2(0,NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  end subroutine Residue2TopsoilSoluteAdvExch
!------------------------------------------------------------------------------------------

  subroutine Topsoil2ResidueSoluteAdvExch(M,NY,NX,FLWRM1)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8) :: VFLW

!     begin_execution
!
!     IF WATER FLUX FROM 'WATSUB' IS TO RESIDUE FROM
!     SOIL SURFACE THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN SOIL SURFACE
!
!     VOLWM=litter water volume
!     RFL*=soil-litter convective solute flux
!     Z*2=soil solute content
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
  IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMIN1(AMAX1(-VFLWX,FLWRM1/VOLWM(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=-VFLWX
  ENDIF
  RFLAL=VFLW*AZMAX1(ZAL2(NU(NY,NX),NY,NX))
  RFLFE=VFLW*AZMAX1(ZFE2(NU(NY,NX),NY,NX))
  RFLHY=VFLW*AZMAX1(ZHY2(NU(NY,NX),NY,NX))
  RFLCA=VFLW*AZMAX1(ZCA2(NU(NY,NX),NY,NX))
  RFLMG=VFLW*AZMAX1(ZMG2(NU(NY,NX),NY,NX))
  RFLNA=VFLW*AZMAX1(ZNA2(NU(NY,NX),NY,NX))
  RFLKA=VFLW*AZMAX1(ZKA2(NU(NY,NX),NY,NX))
  RFLOH=VFLW*AZMAX1(ZOH2(NU(NY,NX),NY,NX))
  RFLSO=VFLW*AZMAX1(ZSO42(NU(NY,NX),NY,NX))
  RFLCL=VFLW*AZMAX1(ZCL2(NU(NY,NX),NY,NX))
  RFLC3=VFLW*AZMAX1(ZCO32(NU(NY,NX),NY,NX))
  RFLHC=VFLW*AZMAX1(ZHCO32(NU(NY,NX),NY,NX))
  RFLAL1=VFLW*AZMAX1(ZAL12(NU(NY,NX),NY,NX))
  RFLAL2=VFLW*AZMAX1(ZAL22(NU(NY,NX),NY,NX))
  RFLAL3=VFLW*AZMAX1(ZAL32(NU(NY,NX),NY,NX))
  RFLAL4=VFLW*AZMAX1(ZAL42(NU(NY,NX),NY,NX))
  RFLALS=VFLW*AZMAX1(ZALS2(NU(NY,NX),NY,NX))
  RFLFE1=VFLW*AZMAX1(ZFE12(NU(NY,NX),NY,NX))
  RFLFE2=VFLW*AZMAX1(ZFE22(NU(NY,NX),NY,NX))
  RFLFE3=VFLW*AZMAX1(ZFE32(NU(NY,NX),NY,NX))
  RFLFE4=VFLW*AZMAX1(ZFE42(NU(NY,NX),NY,NX))
  RFLFES=VFLW*AZMAX1(ZFES2(NU(NY,NX),NY,NX))
  RFLCAO=VFLW*AZMAX1(ZCAO2(NU(NY,NX),NY,NX))
  RFLCAC=VFLW*AZMAX1(ZCAC2(NU(NY,NX),NY,NX))
  RFLCAH=VFLW*AZMAX1(ZCAH2(NU(NY,NX),NY,NX))
  RFLCAS=VFLW*AZMAX1(ZCAS2(NU(NY,NX),NY,NX))
  RFLMGO=VFLW*AZMAX1(ZMGO2(NU(NY,NX),NY,NX))
  RFLMGC=VFLW*AZMAX1(ZMGC2(NU(NY,NX),NY,NX))
  RFLMGH=VFLW*AZMAX1(ZMGH2(NU(NY,NX),NY,NX))
  RFLMGS=VFLW*AZMAX1(ZMGS2(NU(NY,NX),NY,NX))
  RFLNAC=VFLW*AZMAX1(ZNAC2(NU(NY,NX),NY,NX))
  RFLNAS=VFLW*AZMAX1(ZNAS2(NU(NY,NX),NY,NX))
  RFLKAS=VFLW*AZMAX1(ZKAS2(NU(NY,NX),NY,NX))
  RFLH0P=VFLW*AZMAX1(H0PO42(NU(NY,NX),NY,NX))
  RFLH3P=VFLW*AZMAX1(H3PO42(NU(NY,NX),NY,NX))
  RFLF1P=VFLW*AZMAX1(ZFE1P2(NU(NY,NX),NY,NX))
  RFLF2P=VFLW*AZMAX1(ZFE2P2(NU(NY,NX),NY,NX))
  RFLC0P=VFLW*AZMAX1(ZCA0P2(NU(NY,NX),NY,NX))
  RFLC1P=VFLW*AZMAX1(ZCA1P2(NU(NY,NX),NY,NX))
  RFLC2P=VFLW*AZMAX1(ZCA2P2(NU(NY,NX),NY,NX))
  RFLM1P=VFLW*AZMAX1(ZMG1P2(NU(NY,NX),NY,NX))
  RFLH0B=VFLW*AZMAX1(H0POB2(NU(NY,NX),NY,NX))
  RFLH3B=VFLW*AZMAX1(H3POB2(NU(NY,NX),NY,NX))
  RFLF1B=VFLW*AZMAX1(ZF1PB2(NU(NY,NX),NY,NX))
  RFLF2B=VFLW*AZMAX1(ZF2PB2(NU(NY,NX),NY,NX))
  RFLC0B=VFLW*AZMAX1(ZC0PB2(NU(NY,NX),NY,NX))
  RFLC1B=VFLW*AZMAX1(ZC1PB2(NU(NY,NX),NY,NX))
  RFLC2B=VFLW*AZMAX1(ZC2PB2(NU(NY,NX),NY,NX))
  RFLM1B=VFLW*AZMAX1(ZM1PB2(NU(NY,NX),NY,NX))
  end subroutine Topsoil2ResidueSoluteAdvExch
!------------------------------------------------------------------------------------------

  subroutine TopsoilResidueSolutedifusExch(M,NY,NX,FLWRM1)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(in) :: FLWRM1
  real(r8) :: VOLWPB,VOLWPA
!     begin_execution
!
!     MICROPORE CONCENTRATIONS FROM WATER IN RESIDUE AND SOIL SURFACE
!
!     C*1,C*2=solute concentration in litter, soil
!     Z*1,Z*2=solute content in litter, soil
!
  VOLWPA=VOLWM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  VOLWPB=VOLWM(M,NU(NY,NX),NY,NX)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  CAL1=AZMAX1(ZAL2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE1=AZMAX1(ZFE2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CHY1=AZMAX1(ZHY2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CCA1=AZMAX1(ZCA2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CMG1=AZMAX1(ZMG2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CNA1=AZMAX1(ZNA2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CKA1=AZMAX1(ZKA2(0,NY,NX)/VOLWM(M,0,NY,NX))
  COH1=AZMAX1(ZOH2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CSO41=AZMAX1(ZSO42(0,NY,NX)/VOLWM(M,0,NY,NX))
  CCL1=AZMAX1(ZCL2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CCO31=AZMAX1(ZCO32(0,NY,NX)/VOLWM(M,0,NY,NX))
  CHCO31=AZMAX1(ZHCO32(0,NY,NX)/VOLWM(M,0,NY,NX))
  CAL11=AZMAX1(ZAL12(0,NY,NX)/VOLWM(M,0,NY,NX))
  CAL21=AZMAX1(ZAL22(0,NY,NX)/VOLWM(M,0,NY,NX))
  CAL31=AZMAX1(ZAL32(0,NY,NX)/VOLWM(M,0,NY,NX))
  CAL41=AZMAX1(ZAL42(0,NY,NX)/VOLWM(M,0,NY,NX))
  CALS1=AZMAX1(ZALS2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE11=AZMAX1(ZFE12(0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE21=AZMAX1(ZFE22(0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE31=AZMAX1(ZFE32(0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE41=AZMAX1(ZFE42(0,NY,NX)/VOLWM(M,0,NY,NX))
  CFES1=AZMAX1(ZFES2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CCAO1=AZMAX1(ZCAO2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CCAC1=AZMAX1(ZCAC2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CCAH1=AZMAX1(ZCAH2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CCAS1=AZMAX1(ZCAS2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CMGO1=AZMAX1(ZMGO2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CMGC1=AZMAX1(ZMGC2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CMGH1=AZMAX1(ZMGH2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CMGS1=AZMAX1(ZMGS2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CNAC1=AZMAX1(ZNAC2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CNAS1=AZMAX1(ZNAS2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CKAS1=AZMAX1(ZKAS2(0,NY,NX)/VOLWM(M,0,NY,NX))
  C0PO41=AZMAX1(H0PO42(0,NY,NX)/VOLWM(M,0,NY,NX))
  C3PO41=AZMAX1(H3PO42(0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE1P1=AZMAX1(ZFE1P2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CFE2P1=AZMAX1(ZFE2P2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CCA0P1=AZMAX1(ZCA0P2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CCA1P1=AZMAX1(ZCA1P2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CCA2P1=AZMAX1(ZCA2P2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CMG1P1=AZMAX1(ZMG1P2(0,NY,NX)/VOLWM(M,0,NY,NX))
  CAL2=AZMAX1(ZAL2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFE2=AZMAX1(ZFE2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CHY2=AZMAX1(ZHY2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCA2=AZMAX1(ZCA2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CMG2=AZMAX1(ZMG2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CNA2=AZMAX1(ZNA2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CKA2=AZMAX1(ZKA2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  COH2=AZMAX1(ZOH2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CSO42=AZMAX1(ZSO42(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCL2=AZMAX1(ZCL2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCO32=AZMAX1(ZCO32(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CHCO32=AZMAX1(ZHCO32(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CAL12=AZMAX1(ZAL12(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CAL22=AZMAX1(ZAL22(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CAL32=AZMAX1(ZAL32(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CAL42=AZMAX1(ZAL42(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CALS2=AZMAX1(ZALS2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFE12=AZMAX1(ZFE12(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFE22=AZMAX1(ZFE22(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFE32=AZMAX1(ZFE32(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFE42=AZMAX1(ZFE42(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CFES2=AZMAX1(ZFES2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCAO2=AZMAX1(ZCAO2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCAC2=AZMAX1(ZCAC2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCAH2=AZMAX1(ZCAH2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CCAS2=AZMAX1(ZCAS2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CMGO2=AZMAX1(ZMGO2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CMGC2=AZMAX1(ZMGC2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CMGH2=AZMAX1(ZMGH2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CMGS2=AZMAX1(ZMGS2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CNAC2=AZMAX1(ZNAC2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CNAS2=AZMAX1(ZNAS2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  CKAS2=AZMAX1(ZKAS2(NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX))
  IF(VOLWPA.GT.ZEROS2(NY,NX))THEN
    C0PO42=AZMAX1(H0PO42(NU(NY,NX),NY,NX)/VOLWPA)
    C3PO42=AZMAX1(H3PO42(NU(NY,NX),NY,NX)/VOLWPA)
    CFE1P2=AZMAX1(ZFE1P2(NU(NY,NX),NY,NX)/VOLWPA)
    CFE2P2=AZMAX1(ZFE2P2(NU(NY,NX),NY,NX)/VOLWPA)
    CCA0P2=AZMAX1(ZCA0P2(NU(NY,NX),NY,NX)/VOLWPA)
    CCA1P2=AZMAX1(ZCA1P2(NU(NY,NX),NY,NX)/VOLWPA)
    CCA2P2=AZMAX1(ZCA2P2(NU(NY,NX),NY,NX)/VOLWPA)
    CMG1P2=AZMAX1(ZMG1P2(NU(NY,NX),NY,NX)/VOLWPA)
  ELSE
    C0PO42=0.0_r8
    C3PO42=0.0_r8
    CFE1P2=0.0_r8
    CFE2P2=0.0_r8
    CCA0P2=0.0_r8
    CCA1P2=0.0_r8
    CCA2P2=0.0_r8
    CMG1P2=0.0_r8
  ENDIF
  IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
    C0POB2=AZMAX1(H0POB2(NU(NY,NX),NY,NX)/VOLWPB)
    C3POB2=AZMAX1(H3POB2(NU(NY,NX),NY,NX)/VOLWPB)
    CF1PB2=AZMAX1(ZF1PB2(NU(NY,NX),NY,NX)/VOLWPB)
    CF2PB2=AZMAX1(ZF2PB2(NU(NY,NX),NY,NX)/VOLWPB)
    CC0PB2=AZMAX1(ZC0PB2(NU(NY,NX),NY,NX)/VOLWPB)
    CC1PB2=AZMAX1(ZC1PB2(NU(NY,NX),NY,NX)/VOLWPB)
    CC2PB2=AZMAX1(ZC2PB2(NU(NY,NX),NY,NX)/VOLWPB)
    CM1PB2=AZMAX1(ZM1PB2(NU(NY,NX),NY,NX)/VOLWPB)
  ELSE
    C0POB2=C0PO42
    C3POB2=C3PO42
    CF1PB2=CFE1P2
    CF2PB2=CFE2P2
    CC0PB2=CCA0P2
    CC1PB2=CCA1P2
    CC2PB2=CCA2P2
    CM1PB2=CMG1P2
  ENDIF
!
!     DIFFUSIVITIES IN RESIDUE AND SOIL SURFACE
!
!     DLYR0,DLYR1=litter, soil surface thickness
!     TORT=tortuosity from hour1.f
!     CVRD,BARE=litter cover fraction,1-CVRD
!     DISP=dispersivity parameter
!     FLWRM=litter-soil water flux from watsub.f
!     DIF*=aqueous diffusivity-dispersivity in soil surface
!     *SGL*=solute diffusivity from hour1.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     DIF*=aqueous diffusivity-dispersivity between litter and soil surface
!
  DLYR0=AMAX1(ZERO2,DLYR(3,0,NY,NX))
  TORT0=TORT(M,0,NY,NX)*CVRD(NY,NX)
  DLYR1=AMAX1(ZERO2,DLYR(3,NU(NY,NX),NY,NX))
  TORT1=TORT(M,NU(NY,NX),NY,NX)
  TORTL=AMIN1(1.0_r8,(TORT0+TORT1)/(DLYR0+DLYR1))
  DISPN=DISP(3,NU(NY,NX),NY,NX)*AMIN1(VFLWX,ABS(FLWRM1/AREA(3,NU(NY,NX),NY,NX)))
  DIFPO=(POSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFAL=(ALSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFFE=(FESGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFHY=(HYSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFCA=(CASGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFMG=(GMSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFNA=(ANSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFKA=(AKSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFOH=(OHSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFSO=(SOSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFCL=(CLSXL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFC3=(C3SGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
  DIFHC=(HCSGL2(NU(NY,NX),NY,NX)*TORTL+DISPN)*AREA(3,NU(NY,NX),NY,NX)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
!     DFV*S,DFV*B=diffusive solute flux between litter and soil surface in non-band,band
!     DIF*=aqueous diffusivity-dispersivity between litter and soil surface
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     C*1,C*2=solute concentration in litter,soil surface
!
  DFVAL=DIFAL*(CAL1-CAL2)
  DFVFE=DIFFE*(CFE1-CFE2)
  DFVHY=DIFHY*(CHY1-CHY2)
  DFVCA=DIFCA*(CCA1-CCA2)
  DFVMG=DIFMG*(CMG1-CMG2)
  DFVNA=DIFNA*(CNA1-CNA2)
  DFVKA=DIFKA*(CKA1-CKA2)
  DFVOH=DIFOH*(COH1-COH2)
  DFVSO=DIFSO*(CSO41-CSO42)
  DFVCL=DIFCL*(CCL1-CCL2)
  DFVC3=DIFC3*(CCO31-CCO32)
  DFVHC=DIFHC*(CHCO31-CHCO32)
  DFVAL1=DIFAL*(CAL11-CAL12)
  DFVAL2=DIFAL*(CAL21-CAL22)
  DFVAL3=DIFAL*(CAL31-CAL32)
  DFVAL4=DIFAL*(CAL41-CAL42)
  DFVALS=DIFAL*(CALS1-CALS2)
  DFVFE1=DIFFE*(CFE11-CFE12)
  DFVFE2=DIFFE*(CFE21-CFE22)
  DFVFE3=DIFFE*(CFE31-CFE32)
  DFVFE4=DIFFE*(CFE41-CFE42)
  DFVFES=DIFFE*(CFES1-CFES2)
  DFVCAO=DIFCA*(CCAO1-CCAO2)
  DFVCAC=DIFCA*(CCAC1-CCAC2)
  DFVCAH=DIFCA*(CCAH1-CCAH2)
  DFVCAS=DIFCA*(CCAS1-CCAS2)
  DFVMGO=DIFMG*(CMGO1-CMGO2)
  DFVMGC=DIFMG*(CMGC1-CMGC2)
  DFVMGH=DIFMG*(CMGH1-CMGH2)
  DFVMGS=DIFMG*(CMGS1-CMGS2)
  DFVNAC=DIFNA*(CNAC1-CNAC2)
  DFVNAS=DIFNA*(CNAS1-CNAS2)
  DFVKAS=DIFKA*(CKAS1-CKAS2)
  DFVH0P=DIFPO*(C0PO41-C0PO42)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVH3P=DIFPO*(C3PO41-C3PO42)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVF1P=DIFPO*(CFE1P1-CFE1P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVF2P=DIFPO*(CFE2P1-CFE2P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC0P=DIFPO*(CCA0P1-CCA0P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC1P=DIFPO*(CCA1P1-CCA1P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC2P=DIFPO*(CCA2P1-CCA2P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVM1P=DIFPO*(CMG1P1-CMG1P2)*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVH0B=DIFPO*(C0PO41-C0POB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVH3B=DIFPO*(C3PO41-C3POB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVF1B=DIFPO*(CFE1P1-CF1PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVF2B=DIFPO*(CFE2P1-CF2PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC0B=DIFPO*(CCA0P1-CC0PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC1B=DIFPO*(CCA1P1-CC1PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC2B=DIFPO*(CCA2P1-CC2PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVM1B=DIFPO*(CMG1P1-CM1PB2)*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  end subroutine TopsoilResidueSolutedifusExch
!------------------------------------------------------------------------------------------

  subroutine ZeroDifusTopsoilResidueExch
!
!     Description:
!
  implicit none
!     begin_execution
  DFVAL=0.0_r8
  DFVFE=0.0_r8
  DFVHY=0.0_r8
  DFVCA=0.0_r8
  DFVMG=0.0_r8
  DFVNA=0.0_r8
  DFVKA=0.0_r8
  DFVOH=0.0_r8
  DFVSO=0.0_r8
  DFVCL=0.0_r8
  DFVC3=0.0_r8
  DFVHC=0.0_r8
  DFVAL1=0.0_r8
  DFVAL2=0.0_r8
  DFVAL3=0.0_r8
  DFVAL4=0.0_r8
  DFVALS=0.0_r8
  DFVFE1=0.0_r8
  DFVFE2=0.0_r8
  DFVFE3=0.0_r8
  DFVFE4=0.0_r8
  DFVFES=0.0_r8
  DFVCAO=0.0_r8
  DFVCAC=0.0_r8
  DFVCAH=0.0_r8
  DFVCAS=0.0_r8
  DFVMGO=0.0_r8
  DFVMGC=0.0_r8
  DFVMGH=0.0_r8
  DFVMGS=0.0_r8
  DFVNAC=0.0_r8
  DFVNAS=0.0_r8
  DFVKAS=0.0_r8
  DFVH0P=0.0_r8
  DFVH3P=0.0_r8
  DFVF1P=0.0_r8
  DFVF2P=0.0_r8
  DFVC0P=0.0_r8
  DFVC1P=0.0_r8
  DFVC2P=0.0_r8
  DFVM1P=0.0_r8
  DFVH0B=0.0_r8
  DFVH3B=0.0_r8
  DFVF1B=0.0_r8
  DFVF2B=0.0_r8
  DFVC0B=0.0_r8
  DFVC1B=0.0_r8
  DFVC2B=0.0_r8
  DFVM1B=0.0_r8
  end subroutine ZeroDifusTopsoilResidueExch
!------------------------------------------------------------------------------------------

  subroutine TopsoilResidueFluxAdvPlusDifus(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
!     begin_execution
!
!     R*FLS=convective + diffusive solute flux between litter, soil surface
!     R*FLW,R*FLB=convective + diffusive solute flux into soil in non-band,band
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     R*FL0,R*FL1=convective flux into surface litter, soil surface
!     RFL*=convective flux between surface litter and soil surface
!     DFV*=diffusive solute flux between litter and soil surface
!
  RALFLS(3,0,NY,NX)=RALFL0(NY,NX)+RALFLS0-RFLAL-DFVAL
  RFEFLS(3,0,NY,NX)=RFEFL0(NY,NX)+RFEFLS0-RFLFE-DFVFE
  RHYFLS(3,0,NY,NX)=RHYFL0(NY,NX)+RHYFLS0-RFLHY-DFVHY
  RCAFLS(3,0,NY,NX)=RCAFL0(NY,NX)+RCAFLS0-RFLCA-DFVCA
  RMGFLS(3,0,NY,NX)=RMGFL0(NY,NX)+RMGFLS0-RFLMG-DFVMG
  RNAFLS(3,0,NY,NX)=RNAFL0(NY,NX)+RNAFLS0-RFLNA-DFVNA
  RKAFLS(3,0,NY,NX)=RKAFL0(NY,NX)+RKAFLS0-RFLKA-DFVKA
  ROHFLS(3,0,NY,NX)=ROHFL0(NY,NX)+ROHFLS0-RFLOH-DFVOH
  RSOFLS(3,0,NY,NX)=RSOFL0(NY,NX)+RSOFLS0-RFLSO-DFVSO
  RCLFLS(3,0,NY,NX)=RCLFL0(NY,NX)+RCLFLS0-RFLCL-DFVCL
  RC3FLS(3,0,NY,NX)=RC3FL0(NY,NX)+RC3FLS0-RFLC3-DFVC3
  RHCFLS(3,0,NY,NX)=RHCFL0(NY,NX)+RHCFLS0-RFLHC-DFVHC
  RAL1FS(3,0,NY,NX)=RAL1F0(NY,NX)+RAL1FS0-RFLAL1-DFVAL1
  RAL2FS(3,0,NY,NX)=RAL2F0(NY,NX)+RAL2FS0-RFLAL2-DFVAL2
  RAL3FS(3,0,NY,NX)=RAL3F0(NY,NX)+RAL3FS0-RFLAL3-DFVAL3
  RAL4FS(3,0,NY,NX)=RAL4F0(NY,NX)+RAL4FS0-RFLAL4-DFVAL4
  RALSFS(3,0,NY,NX)=RALSF0(NY,NX)+RALSFS0-RFLALS-DFVALS
  RFE1FS(3,0,NY,NX)=RFE1F0(NY,NX)+RFE1FS0-RFLFE1-DFVFE1
  RFE2FS(3,0,NY,NX)=RFE2F0(NY,NX)+RFE2FS0-RFLFE2-DFVFE2
  RFE3FS(3,0,NY,NX)=RFE3F0(NY,NX)+RFE3FS0-RFLFE3-DFVFE3
  RFE4FS(3,0,NY,NX)=RFE4F0(NY,NX)+RFE4FS0-RFLFE4-DFVFE4
  RFESFS(3,0,NY,NX)=RFESF0(NY,NX)+RFESFS0-RFLFES-DFVFES
  RCAOFS(3,0,NY,NX)=RCAOF0(NY,NX)+RCAOFS0-RFLCAO-DFVCAO
  RCACFS(3,0,NY,NX)=RCACF0(NY,NX)+RCACFS0-RFLCAC-DFVCAC
  RCAHFS(3,0,NY,NX)=RCAHF0(NY,NX)+RCAHFS0-RFLCAH-DFVCAH
  RCASFS(3,0,NY,NX)=RCASF0(NY,NX)+RCASFS0-RFLCAS-DFVCAS
  RMGOFS(3,0,NY,NX)=RMGOF0(NY,NX)+RMGOFS0-RFLMGO-DFVMGO
  RMGCFS(3,0,NY,NX)=RMGCF0(NY,NX)+RMGCFS0-RFLMGC-DFVMGC
  RMGHFS(3,0,NY,NX)=RMGHF0(NY,NX)+RMGHFS0-RFLMGH-DFVMGH
  RMGSFS(3,0,NY,NX)=RMGSF0(NY,NX)+RMGSFS0-RFLMGS-DFVMGS
  RNACFS(3,0,NY,NX)=RNACF0(NY,NX)+RNACFS0-RFLNAC-DFVNAC
  RNASFS(3,0,NY,NX)=RNASF0(NY,NX)+RNASFS0-RFLNAS-DFVNAS
  RKASFS(3,0,NY,NX)=RKASF0(NY,NX)+RKASFS0-RFLKAS-DFVKAS
  RH0PFS(3,0,NY,NX)=RH0PF0(NY,NX)+RH0PFS0-RFLH0P-DFVH0P-RFLH0B-DFVH0B
  RH3PFS(3,0,NY,NX)=RH3PF0(NY,NX)+RH3PFS0-RFLH3P-DFVH3P-RFLH3B-DFVH3B
  RF1PFS(3,0,NY,NX)=RF1PF0(NY,NX)+RF1PFS0-RFLF1P-DFVF1P-RFLF1B-DFVF1B
  RF2PFS(3,0,NY,NX)=RF2PF0(NY,NX)+RF2PFS0-RFLF2P-DFVF2P-RFLF2B-DFVF2B
  RC0PFS(3,0,NY,NX)=RC0PF0(NY,NX)+RC0PFS0-RFLC0P-DFVC0P-RFLC0B-DFVC0B
  RC1PFS(3,0,NY,NX)=RC1PF0(NY,NX)+RC1PFS0-RFLC1P-DFVC1P-RFLC1B-DFVC1B
  RC2PFS(3,0,NY,NX)=RC2PF0(NY,NX)+RC2PFS0-RFLC2P-DFVC2P-RFLC2B-DFVC2B
  RM1PFS(3,0,NY,NX)=RM1PF0(NY,NX)+RM1PFS0-RFLM1P-DFVM1P-RFLM1B-DFVM1B

  RALFLS(3,NU(NY,NX),NY,NX)=RALFL1(NY,NX)+RALFLS1+RFLAL+DFVAL
  RFEFLS(3,NU(NY,NX),NY,NX)=RFEFL1(NY,NX)+RFEFLS1+RFLFE+DFVFE
  RHYFLS(3,NU(NY,NX),NY,NX)=RHYFL1(NY,NX)+RHYFLS1+RFLHY+DFVHY
  RCAFLS(3,NU(NY,NX),NY,NX)=RCAFL1(NY,NX)+RCAFLS1+RFLCA+DFVCA
  RMGFLS(3,NU(NY,NX),NY,NX)=RMGFL1(NY,NX)+RMGFLS1+RFLMG+DFVMG
  RNAFLS(3,NU(NY,NX),NY,NX)=RNAFL1(NY,NX)+RNAFLS1+RFLNA+DFVNA
  RKAFLS(3,NU(NY,NX),NY,NX)=RKAFL1(NY,NX)+RKAFLS1+RFLKA+DFVKA
  ROHFLS(3,NU(NY,NX),NY,NX)=ROHFL1(NY,NX)+ROHFLS1+RFLOH+DFVOH
  RSOFLS(3,NU(NY,NX),NY,NX)=RSOFL1(NY,NX)+RSOFLS1+RFLSO+DFVSO
  RCLFLS(3,NU(NY,NX),NY,NX)=RCLFL1(NY,NX)+RCLFLS1+RFLCL+DFVCL
  RC3FLS(3,NU(NY,NX),NY,NX)=RC3FL1(NY,NX)+RC3FLS1+RFLC3+DFVC3
  RHCFLS(3,NU(NY,NX),NY,NX)=RHCFL1(NY,NX)+RHCFLS1+RFLHC+DFVHC
  RAL1FS(3,NU(NY,NX),NY,NX)=RAL1F1(NY,NX)+RAL1FS1+RFLAL1+DFVAL1
  RAL2FS(3,NU(NY,NX),NY,NX)=RAL2F1(NY,NX)+RAL2FS1+RFLAL2+DFVAL2
  RAL3FS(3,NU(NY,NX),NY,NX)=RAL3F1(NY,NX)+RAL3FS1+RFLAL3+DFVAL3
  RAL4FS(3,NU(NY,NX),NY,NX)=RAL4F1(NY,NX)+RAL4FS1+RFLAL4+DFVAL4
  RALSFS(3,NU(NY,NX),NY,NX)=RALSF1(NY,NX)+RALSFS1+RFLALS+DFVALS
  RFE1FS(3,NU(NY,NX),NY,NX)=RFE1F1(NY,NX)+RFE1FS1+RFLFE1+DFVFE1
  RFE2FS(3,NU(NY,NX),NY,NX)=RFE2F1(NY,NX)+RFE2FS1+RFLFE2+DFVFE2
  RFE3FS(3,NU(NY,NX),NY,NX)=RFE3F1(NY,NX)+RFE3FS1+RFLFE3+DFVFE3
  RFE4FS(3,NU(NY,NX),NY,NX)=RFE4F1(NY,NX)+RFE4FS1+RFLFE4+DFVFE4
  RFESFS(3,NU(NY,NX),NY,NX)=RFESF1(NY,NX)+RFESFS1+RFLFES+DFVFES
  RCAOFS(3,NU(NY,NX),NY,NX)=RCAOF1(NY,NX)+RCAOFS1+RFLCAO+DFVCAO
  RCACFS(3,NU(NY,NX),NY,NX)=RCACF1(NY,NX)+RCACFS1+RFLCAC+DFVCAC
  RCAHFS(3,NU(NY,NX),NY,NX)=RCAHF1(NY,NX)+RCAHFS1+RFLCAH+DFVCAH
  RCASFS(3,NU(NY,NX),NY,NX)=RCASF1(NY,NX)+RCASFS1+RFLCAS+DFVCAS
  RMGOFS(3,NU(NY,NX),NY,NX)=RMGOF1(NY,NX)+RMGOFS1+RFLMGO+DFVMGO
  RMGCFS(3,NU(NY,NX),NY,NX)=RMGCF1(NY,NX)+RMGCFS1+RFLMGC+DFVMGC
  RMGHFS(3,NU(NY,NX),NY,NX)=RMGHF1(NY,NX)+RMGHFS1+RFLMGH+DFVMGH
  RMGSFS(3,NU(NY,NX),NY,NX)=RMGSF1(NY,NX)+RMGSFS1+RFLMGS+DFVMGS
  RNACFS(3,NU(NY,NX),NY,NX)=RNACF1(NY,NX)+RNACFS1+RFLNAC+DFVNAC
  RNASFS(3,NU(NY,NX),NY,NX)=RNASF1(NY,NX)+RNASFS1+RFLNAS+DFVNAS
  RKASFS(3,NU(NY,NX),NY,NX)=RKASF1(NY,NX)+RKASFS1+RFLKAS+DFVKAS
  RH0PFS(3,NU(NY,NX),NY,NX)=RH0PF1(NY,NX)+RH0PFS1+RFLH0P+DFVH0P
  RH3PFS(3,NU(NY,NX),NY,NX)=RH3PF1(NY,NX)+RH3PFS1+RFLH3P+DFVH3P
  RF1PFS(3,NU(NY,NX),NY,NX)=RF1PF1(NY,NX)+RF1PFS1+RFLF1P+DFVF1P
  RF2PFS(3,NU(NY,NX),NY,NX)=RF2PF1(NY,NX)+RF2PFS1+RFLF2P+DFVF2P
  RC0PFS(3,NU(NY,NX),NY,NX)=RC0PF1(NY,NX)+RC0PFS1+RFLC0P+DFVC0P
  RC1PFS(3,NU(NY,NX),NY,NX)=RC1PF1(NY,NX)+RC1PFS1+RFLC1P+DFVC1P
  RC2PFS(3,NU(NY,NX),NY,NX)=RC2PF1(NY,NX)+RC2PFS1+RFLC2P+DFVC2P
  RM1PFS(3,NU(NY,NX),NY,NX)=RM1PF1(NY,NX)+RM1PFS1+RFLM1P+DFVM1P
  RH0BFB(3,NU(NY,NX),NY,NX)=RH0BF2(NY,NX)+RH0BFB1+RFLH0B+DFVH0B
  RH3BFB(3,NU(NY,NX),NY,NX)=RH3BF2(NY,NX)+RH3BFB1+RFLH3B+DFVH3B
  RF1BFB(3,NU(NY,NX),NY,NX)=RF1BF2(NY,NX)+RF1BFB1+RFLF1B+DFVF1B
  RF2BFB(3,NU(NY,NX),NY,NX)=RF2BF2(NY,NX)+RF2BFB1+RFLF2B+DFVF2B
  RC0BFB(3,NU(NY,NX),NY,NX)=RC0BF2(NY,NX)+RC0BFB1+RFLC0B+DFVC0B
  RC1BFB(3,NU(NY,NX),NY,NX)=RC1BF2(NY,NX)+RC1BFB1+RFLC1B+DFVC1B
  RC2BFB(3,NU(NY,NX),NY,NX)=RC2BF2(NY,NX)+RC2BFB1+RFLC2B+DFVC2B
  RM1BFB(3,NU(NY,NX),NY,NX)=RM1BF2(NY,NX)+RM1BFB1+RFLM1B+DFVM1B
  end subroutine TopsoilResidueFluxAdvPlusDifus
!------------------------------------------------------------------------------------------

  subroutine AccumHourlyTopsoilReisdueFlux(NY,NX)
!
!     Description:
!
      implicit none
      integer, intent(in) :: NY,NX
!     begin_execution
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux
!     X*FLW,X*FLB= hourly convective + diffusive solute flux in non-band,band
!
  trcsa_XFLS(idsa_Al,3,0,NY,NX)=trcsa_XFLS(idsa_Al,3,0,NY,NX)+RALFLS0-RFLAL-DFVAL
  trcsa_XFLS(idsa_Fe,3,0,NY,NX)=trcsa_XFLS(idsa_Fe,3,0,NY,NX)+RFEFLS0-RFLFE-DFVFE
  trcsa_XFLS(idsa_Hp,3,0,NY,NX)=trcsa_XFLS(idsa_Hp,3,0,NY,NX)+RHYFLS0-RFLHY-DFVHY
  trcsa_XFLS(idsa_Ca,3,0,NY,NX)=trcsa_XFLS(idsa_Ca,3,0,NY,NX)+RCAFLS0-RFLCA-DFVCA
  trcsa_XFLS(idsa_Mg,3,0,NY,NX)=trcsa_XFLS(idsa_Mg,3,0,NY,NX)+RMGFLS0-RFLMG-DFVMG
  trcsa_XFLS(idsa_Na,3,0,NY,NX)=trcsa_XFLS(idsa_Na,3,0,NY,NX)+RNAFLS0-RFLNA-DFVNA
  trcsa_XFLS(idsa_K,3,0,NY,NX)=trcsa_XFLS(idsa_K,3,0,NY,NX)+RKAFLS0-RFLKA-DFVKA
  trcsa_XFLS(idsa_OH,3,0,NY,NX)=trcsa_XFLS(idsa_OH,3,0,NY,NX)+ROHFLS0-RFLOH-DFVOH
  trcsa_XFLS(idsa_SO4,3,0,NY,NX)=trcsa_XFLS(idsa_SO4,3,0,NY,NX)+RSOFLS0-RFLSO-DFVSO
  trcsa_XFLS(idsa_Cl,3,0,NY,NX)=trcsa_XFLS(idsa_Cl,3,0,NY,NX)+RCLFLS0-RFLCL-DFVCL
  trcsa_XFLS(idsa_CO3,3,0,NY,NX)=trcsa_XFLS(idsa_CO3,3,0,NY,NX)+RC3FLS0-RFLC3-DFVC3
  trcsa_XFLS(idsa_HCO3,3,0,NY,NX)=trcsa_XFLS(idsa_HCO3,3,0,NY,NX)+RHCFLS0-RFLHC-DFVHC
  trcsa_XFLS(idsa_AlOH,3,0,NY,NX)=trcsa_XFLS(idsa_AlOH,3,0,NY,NX)+RAL1FS0-RFLAL1-DFVAL1
  trcsa_XFLS(idsa_AlOH2,3,0,NY,NX)=trcsa_XFLS(idsa_AlOH2,3,0,NY,NX)+RAL2FS0-RFLAL2-DFVAL2
  trcsa_XFLS(idsa_AlOH3,3,0,NY,NX)=trcsa_XFLS(idsa_AlOH3,3,0,NY,NX)+RAL3FS0-RFLAL3-DFVAL3
  trcsa_XFLS(idsa_AlOH4,3,0,NY,NX)=trcsa_XFLS(idsa_AlOH4,3,0,NY,NX)+RAL4FS0-RFLAL4-DFVAL4
  trcsa_XFLS(idsa_AlSO4,3,0,NY,NX)=trcsa_XFLS(idsa_AlSO4,3,0,NY,NX)+RALSFS0-RFLALS-DFVALS
  trcsa_XFLS(idsa_FeOH,3,0,NY,NX)=trcsa_XFLS(idsa_FeOH,3,0,NY,NX)+RFE1FS0-RFLFE1-DFVFE1
  trcsa_XFLS(idsa_FeOH2,3,0,NY,NX)=trcsa_XFLS(idsa_FeOH2,3,0,NY,NX)+RFE2FS0-RFLFE2-DFVFE2
  trcsa_XFLS(idsa_FeOH3,3,0,NY,NX)=trcsa_XFLS(idsa_FeOH3,3,0,NY,NX)+RFE3FS0-RFLFE3-DFVFE3
  trcsa_XFLS(idsa_FeOH4,3,0,NY,NX)=trcsa_XFLS(idsa_FeOH4,3,0,NY,NX)+RFE4FS0-RFLFE4-DFVFE4
  trcsa_XFLS(idsa_FeSO4,3,0,NY,NX)=trcsa_XFLS(idsa_FeSO4,3,0,NY,NX)+RFESFS0-RFLFES-DFVFES
  trcsa_XFLS(idsa_CaOH2,3,0,NY,NX)=trcsa_XFLS(idsa_CaOH2,3,0,NY,NX)+RCAOFS0-RFLCAO-DFVCAO
  trcsa_XFLS(idsa_CaCO3,3,0,NY,NX)=trcsa_XFLS(idsa_CaCO3,3,0,NY,NX)+RCACFS0-RFLCAC-DFVCAC
  trcsa_XFLS(idsa_CaHCO3,3,0,NY,NX)=trcsa_XFLS(idsa_CaHCO3,3,0,NY,NX)+RCAHFS0-RFLCAH-DFVCAH
  trcsa_XFLS(idsa_CaSO4,3,0,NY,NX)=trcsa_XFLS(idsa_CaSO4,3,0,NY,NX)+RCASFS0-RFLCAS-DFVCAS
  trcsa_XFLS(idsa_MgOH2,3,0,NY,NX)=trcsa_XFLS(idsa_MgOH2,3,0,NY,NX)+RMGOFS0-RFLMGO-DFVMGO
  trcsa_XFLS(idsa_MgCO3,3,0,NY,NX)=trcsa_XFLS(idsa_MgCO3,3,0,NY,NX)+RMGCFS0-RFLMGC-DFVMGC
  trcsa_XFLS(idsa_MgHCO3,3,0,NY,NX)=trcsa_XFLS(idsa_MgHCO3,3,0,NY,NX)+RMGHFS0-RFLMGH-DFVMGH
  trcsa_XFLS(idsa_MgSO4,3,0,NY,NX)=trcsa_XFLS(idsa_MgSO4,3,0,NY,NX)+RMGSFS0-RFLMGS-DFVMGS
  trcsa_XFLS(idsa_NaCO3,3,0,NY,NX)=trcsa_XFLS(idsa_NaCO3,3,0,NY,NX)+RNACFS0-RFLNAC-DFVNAC
  trcsa_XFLS(idsa_NaSO4,3,0,NY,NX)=trcsa_XFLS(idsa_NaSO4,3,0,NY,NX)+RNASFS0-RFLNAS-DFVNAS
  trcsa_XFLS(idsa_KSO4,3,0,NY,NX)=trcsa_XFLS(idsa_KSO4,3,0,NY,NX)+RKASFS0-RFLKAS-DFVKAS
  trcsa_XFLS(idsa_H0PO4,3,0,NY,NX)=trcsa_XFLS(idsa_H0PO4,3,0,NY,NX)+RH0PFS0-RFLH0P-DFVH0P-RFLH0B-DFVH0B
  trcsa_XFLS(idsa_H3PO4,3,0,NY,NX)=trcsa_XFLS(idsa_H3PO4,3,0,NY,NX)+RH3PFS0-RFLH3P-DFVH3P-RFLH3B-DFVH3B
  trcsa_XFLS(idsa_FeHPO4,3,0,NY,NX)=trcsa_XFLS(idsa_FeHPO4,3,0,NY,NX)+RF1PFS0-RFLF1P-DFVF1P-RFLF1B-DFVF1B
  trcsa_XFLS(idsa_FeH2PO4,3,0,NY,NX)=trcsa_XFLS(idsa_FeH2PO4,3,0,NY,NX)+RF2PFS0-RFLF2P-DFVF2P-RFLF2B-DFVF2B
  trcsa_XFLS(idsa_CaPO4,3,0,NY,NX)=trcsa_XFLS(idsa_CaPO4,3,0,NY,NX)+RC0PFS0-RFLC0P-DFVC0P-RFLC0B-DFVC0B
  trcsa_XFLS(idsa_CaHPO4,3,0,NY,NX)=trcsa_XFLS(idsa_CaHPO4,3,0,NY,NX)+RC1PFS0-RFLC1P-DFVC1P-RFLC1B-DFVC1B
  trcsa_XFLS(idsa_CaH2PO4,3,0,NY,NX)=trcsa_XFLS(idsa_CaH2PO4,3,0,NY,NX)+RC2PFS0-RFLC2P-DFVC2P-RFLC2B-DFVC2B
  trcsa_XFLS(idsa_MgHPO4,3,0,NY,NX)=trcsa_XFLS(idsa_MgHPO4,3,0,NY,NX)+RM1PFS0-RFLM1P-DFVM1P-RFLM1B-DFVM1B

  trcsa_XFLS(idsa_Al,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Al,3,NU(NY,NX),NY,NX)+RALFLS1+RFLAL+DFVAL
  trcsa_XFLS(idsa_Fe,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Fe,3,NU(NY,NX),NY,NX)+RFEFLS1+RFLFE+DFVFE
  trcsa_XFLS(idsa_Hp,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Hp,3,NU(NY,NX),NY,NX)+RHYFLS1+RFLHY+DFVHY
  trcsa_XFLS(idsa_Ca,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Ca,3,NU(NY,NX),NY,NX)+RCAFLS1+RFLCA+DFVCA
  trcsa_XFLS(idsa_Mg,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Mg,3,NU(NY,NX),NY,NX)+RMGFLS1+RFLMG+DFVMG
  trcsa_XFLS(idsa_Na,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Na,3,NU(NY,NX),NY,NX)+RNAFLS1+RFLNA+DFVNA
  trcsa_XFLS(idsa_K,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_K,3,NU(NY,NX),NY,NX)+RKAFLS1+RFLKA+DFVKA
  trcsa_XFLS(idsa_OH,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_OH,3,NU(NY,NX),NY,NX)+ROHFLS1+RFLOH+DFVOH
  trcsa_XFLS(idsa_SO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_SO4,3,NU(NY,NX),NY,NX)+RSOFLS1+RFLSO+DFVSO
  trcsa_XFLS(idsa_Cl,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_Cl,3,NU(NY,NX),NY,NX)+RCLFLS1+RFLCL+DFVCL
  trcsa_XFLS(idsa_CO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CO3,3,NU(NY,NX),NY,NX)+RC3FLS1+RFLC3+DFVC3
  trcsa_XFLS(idsa_HCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_HCO3,3,NU(NY,NX),NY,NX)+RHCFLS1+RFLHC+DFVHC
  trcsa_XFLS(idsa_AlOH,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_AlOH,3,NU(NY,NX),NY,NX)+RAL1FS1+RFLAL1+DFVAL1
  trcsa_XFLS(idsa_AlOH2,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_AlOH2,3,NU(NY,NX),NY,NX)+RAL2FS1+RFLAL2+DFVAL2
  trcsa_XFLS(idsa_AlOH3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_AlOH3,3,NU(NY,NX),NY,NX)+RAL3FS1+RFLAL3+DFVAL3
  trcsa_XFLS(idsa_AlOH4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_AlOH4,3,NU(NY,NX),NY,NX)+RAL4FS1+RFLAL4+DFVAL4
  trcsa_XFLS(idsa_AlSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_AlSO4,3,NU(NY,NX),NY,NX)+RALSFS1+RFLALS+DFVALS
  trcsa_XFLS(idsa_FeOH,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeOH,3,NU(NY,NX),NY,NX)+RFE1FS1+RFLFE1+DFVFE1
  trcsa_XFLS(idsa_FeOH2,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeOH2,3,NU(NY,NX),NY,NX)+RFE2FS1+RFLFE2+DFVFE2
  trcsa_XFLS(idsa_FeOH3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeOH3,3,NU(NY,NX),NY,NX)+RFE3FS1+RFLFE3+DFVFE3

  trcsa_XFLS(idsa_FeOH4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeOH4,3,NU(NY,NX),NY,NX)+RFE4FS1+RFLFE4+DFVFE4
  trcsa_XFLS(idsa_FeSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeSO4,3,NU(NY,NX),NY,NX)+RFESFS1+RFLFES+DFVFES
  trcsa_XFLS(idsa_CaOH2,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaOH2,3,NU(NY,NX),NY,NX)+RCAOFS1+RFLCAO+DFVCAO
  trcsa_XFLS(idsa_CaCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaCO3,3,NU(NY,NX),NY,NX)+RCACFS1+RFLCAC+DFVCAC
  trcsa_XFLS(idsa_CaHCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaHCO3,3,NU(NY,NX),NY,NX)+RCAHFS1+RFLCAH+DFVCAH
  trcsa_XFLS(idsa_CaSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaSO4,3,NU(NY,NX),NY,NX)+RCASFS1+RFLCAS+DFVCAS
  trcsa_XFLS(idsa_MgOH2,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgOH2,3,NU(NY,NX),NY,NX)+RMGOFS1+RFLMGO+DFVMGO
  trcsa_XFLS(idsa_MgCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgCO3,3,NU(NY,NX),NY,NX)+RMGCFS1+RFLMGC+DFVMGC
  trcsa_XFLS(idsa_MgHCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgHCO3,3,NU(NY,NX),NY,NX)+RMGHFS1+RFLMGH+DFVMGH
  trcsa_XFLS(idsa_MgSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgSO4,3,NU(NY,NX),NY,NX)+RMGSFS1+RFLMGS+DFVMGS
  trcsa_XFLS(idsa_NaCO3,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_NaCO3,3,NU(NY,NX),NY,NX)+RNACFS1+RFLNAC+DFVNAC
  trcsa_XFLS(idsa_NaSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_NaSO4,3,NU(NY,NX),NY,NX)+RNASFS1+RFLNAS+DFVNAS
  trcsa_XFLS(idsa_KSO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_KSO4,3,NU(NY,NX),NY,NX)+RKASFS1+RFLKAS+DFVKAS
  trcsa_XFLS(idsa_H0PO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_H0PO4,3,NU(NY,NX),NY,NX)+RH0PFS1+RFLH0P+DFVH0P
  trcsa_XFLS(idsa_H3PO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_H3PO4,3,NU(NY,NX),NY,NX)+RH3PFS1+RFLH3P+DFVH3P
  trcsa_XFLS(idsa_FeHPO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeHPO4,3,NU(NY,NX),NY,NX)+RF1PFS1+RFLF1P+DFVF1P
  trcsa_XFLS(idsa_FeH2PO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_FeH2PO4,3,NU(NY,NX),NY,NX)+RF2PFS1+RFLF2P+DFVF2P
  trcsa_XFLS(idsa_CaPO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaPO4,3,NU(NY,NX),NY,NX)+RC0PFS1+RFLC0P+DFVC0P
  trcsa_XFLS(idsa_CaHPO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaHPO4,3,NU(NY,NX),NY,NX)+RC1PFS1+RFLC1P+DFVC1P
  trcsa_XFLS(idsa_CaH2PO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_CaH2PO4,3,NU(NY,NX),NY,NX)+RC2PFS1+RFLC2P+DFVC2P
  trcsa_XFLS(idsa_MgHPO4,3,NU(NY,NX),NY,NX)=trcsa_XFLS(idsa_MgHPO4,3,NU(NY,NX),NY,NX)+RM1PFS1+RFLM1P+DFVM1P
  XH0BFB(3,NU(NY,NX),NY,NX)=XH0BFB(3,NU(NY,NX),NY,NX)+RH0BFB1+RFLH0B+DFVH0B
  XH3BFB(3,NU(NY,NX),NY,NX)=XH3BFB(3,NU(NY,NX),NY,NX)+RH3BFB1+RFLH3B+DFVH3B
  XF1BFB(3,NU(NY,NX),NY,NX)=XF1BFB(3,NU(NY,NX),NY,NX)+RF1BFB1+RFLF1B+DFVF1B
  XF2BFB(3,NU(NY,NX),NY,NX)=XF2BFB(3,NU(NY,NX),NY,NX)+RF2BFB1+RFLF2B+DFVF2B
  XC0BFB(3,NU(NY,NX),NY,NX)=XC0BFB(3,NU(NY,NX),NY,NX)+RC0BFB1+RFLC0B+DFVC0B
  XC1BFB(3,NU(NY,NX),NY,NX)=XC1BFB(3,NU(NY,NX),NY,NX)+RC1BFB1+RFLC1B+DFVC1B
  XC2BFB(3,NU(NY,NX),NY,NX)=XC2BFB(3,NU(NY,NX),NY,NX)+RC2BFB1+RFLC2B+DFVC2B
  XM1BFB(3,NU(NY,NX),NY,NX)=XM1BFB(3,NU(NY,NX),NY,NX)+RM1BFB1+RFLM1B+DFVM1B
  end subroutine AccumHourlyTopsoilReisdueFlux
!------------------------------------------------------------------------------------------

  subroutine MacToMicPoreSoluteAdvExchange(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VFLW
!     begin_execution

  IF(VOLWHM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMAX1(AMIN1(VFLWX,FINHM(M,NU(NY,NX),NY,NX)/VOLWHM(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=VFLWX
  ENDIF
  RFLAL=VFLW*AZMAX1(ZALH2(NU(NY,NX),NY,NX))
  RFLFE=VFLW*AZMAX1(ZFEH2(NU(NY,NX),NY,NX))
  RFLHY=VFLW*AZMAX1(ZHYH2(NU(NY,NX),NY,NX))
  RFLCA=VFLW*AZMAX1(ZCCH2(NU(NY,NX),NY,NX))
  RFLMG=VFLW*AZMAX1(ZMAH2(NU(NY,NX),NY,NX))
  RFLNA=VFLW*AZMAX1(ZNAH2(NU(NY,NX),NY,NX))
  RFLKA=VFLW*AZMAX1(ZKAH2(NU(NY,NX),NY,NX))
  RFLOH=VFLW*AZMAX1(ZOHH2(NU(NY,NX),NY,NX))
  RFLSO=VFLW*AZMAX1(ZSO4H2(NU(NY,NX),NY,NX))
  RFLCL=VFLW*AZMAX1(ZCLH2(NU(NY,NX),NY,NX))
  RFLC3=VFLW*AZMAX1(ZCO3H2(NU(NY,NX),NY,NX))
  RFLHC=VFLW*AZMAX1(ZHCOH2(NU(NY,NX),NY,NX))
  RFLAL1=VFLW*AZMAX1(ZAL1H2(NU(NY,NX),NY,NX))
  RFLAL2=VFLW*AZMAX1(ZAL2H2(NU(NY,NX),NY,NX))
  RFLAL3=VFLW*AZMAX1(ZAL3H2(NU(NY,NX),NY,NX))
  RFLAL4=VFLW*AZMAX1(ZAL4H2(NU(NY,NX),NY,NX))
  RFLALS=VFLW*AZMAX1(ZALSH2(NU(NY,NX),NY,NX))
  RFLFE1=VFLW*AZMAX1(ZFE1H2(NU(NY,NX),NY,NX))
  RFLFE2=VFLW*AZMAX1(ZFE2H2(NU(NY,NX),NY,NX))
  RFLFE3=VFLW*AZMAX1(ZFE3H2(NU(NY,NX),NY,NX))
  RFLFE4=VFLW*AZMAX1(ZFE4H2(NU(NY,NX),NY,NX))
  RFLFES=VFLW*AZMAX1(ZFESH2(NU(NY,NX),NY,NX))
  RFLCAO=VFLW*AZMAX1(ZCAOH2(NU(NY,NX),NY,NX))
  RFLCAC=VFLW*AZMAX1(ZCACH2(NU(NY,NX),NY,NX))
  RFLCAH=VFLW*AZMAX1(ZCAHH2(NU(NY,NX),NY,NX))
  RFLCAS=VFLW*AZMAX1(ZCASH2(NU(NY,NX),NY,NX))
  RFLMGO=VFLW*AZMAX1(ZMGOH2(NU(NY,NX),NY,NX))
  RFLMGC=VFLW*AZMAX1(ZMGCH2(NU(NY,NX),NY,NX))
  RFLMGH=VFLW*AZMAX1(ZMGHH2(NU(NY,NX),NY,NX))
  RFLMGS=VFLW*AZMAX1(ZMGSH2(NU(NY,NX),NY,NX))
  RFLNAC=VFLW*AZMAX1(ZNACH2(NU(NY,NX),NY,NX))
  RFLNAS=VFLW*AZMAX1(ZNASH2(NU(NY,NX),NY,NX))
  RFLKAS=VFLW*AZMAX1(ZKASH2(NU(NY,NX),NY,NX))
  RFLH0P=VFLW*AZMAX1(H0P4H2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH3P=VFLW*AZMAX1(H3P4H2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF1P=VFLW*AZMAX1(ZF1PH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF2P=VFLW*AZMAX1(ZF2PH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC0P=VFLW*AZMAX1(ZC0PH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC1P=VFLW*AZMAX1(ZC1PH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC2P=VFLW*AZMAX1(ZC2PH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLM1P=VFLW*AZMAX1(ZM1PH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH0B=VFLW*AZMAX1(H0PBH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLH3B=VFLW*AZMAX1(H3PBH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF1B=VFLW*AZMAX1(ZF1BH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF2B=VFLW*AZMAX1(ZF2BH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC0B=VFLW*AZMAX1(ZC0BH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC1B=VFLW*AZMAX1(ZC1BH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC2B=VFLW*AZMAX1(ZC2BH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLM1B=VFLW*AZMAX1(ZM1BH2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)

  RALFXS(NU(NY,NX),NY,NX)=RFLAL
  RFEFXS(NU(NY,NX),NY,NX)=RFLFE
  RHYFXS(NU(NY,NX),NY,NX)=RFLHY
  RCAFXS(NU(NY,NX),NY,NX)=RFLCA
  RMGFXS(NU(NY,NX),NY,NX)=RFLMG
  RNAFXS(NU(NY,NX),NY,NX)=RFLNA
  RKAFXS(NU(NY,NX),NY,NX)=RFLKA
  ROHFXS(NU(NY,NX),NY,NX)=RFLOH
  RSOFXS(NU(NY,NX),NY,NX)=RFLSO
  RCLFXS(NU(NY,NX),NY,NX)=RFLCL
  RC3FXS(NU(NY,NX),NY,NX)=RFLC3
  RHCFXS(NU(NY,NX),NY,NX)=RFLHC
  RAL1XS(NU(NY,NX),NY,NX)=RFLAL1
  RAL2XS(NU(NY,NX),NY,NX)=RFLAL2
  RAL3XS(NU(NY,NX),NY,NX)=RFLAL3
  RAL4XS(NU(NY,NX),NY,NX)=RFLAL4
  RALSXS(NU(NY,NX),NY,NX)=RFLALS
  RFE1XS(NU(NY,NX),NY,NX)=RFLFE1
  RFE2XS(NU(NY,NX),NY,NX)=RFLFE2
  RFE3XS(NU(NY,NX),NY,NX)=RFLFE3
  RFE4XS(NU(NY,NX),NY,NX)=RFLFE4
  RFESXS(NU(NY,NX),NY,NX)=RFLFES
  RCAOXS(NU(NY,NX),NY,NX)=RFLCAO
  RCACXS(NU(NY,NX),NY,NX)=RFLCAC
  RCAHXS(NU(NY,NX),NY,NX)=RFLCAH
  RCASXS(NU(NY,NX),NY,NX)=RFLCAS
  RMGOXS(NU(NY,NX),NY,NX)=RFLMGO
  RMGCXS(NU(NY,NX),NY,NX)=RFLMGC
  RMGHXS(NU(NY,NX),NY,NX)=RFLMGH
  RMGSXS(NU(NY,NX),NY,NX)=RFLMGS
  RNACXS(NU(NY,NX),NY,NX)=RFLNAC
  RNASXS(NU(NY,NX),NY,NX)=RFLNAS
  RKASXS(NU(NY,NX),NY,NX)=RFLKAS
  RH0PXS(NU(NY,NX),NY,NX)=RFLH0P
  RH3PXS(NU(NY,NX),NY,NX)=RFLH3P
  RF1PXS(NU(NY,NX),NY,NX)=RFLF1P
  RF2PXS(NU(NY,NX),NY,NX)=RFLF2P
  RC0PXS(NU(NY,NX),NY,NX)=RFLC0P
  RC1PXS(NU(NY,NX),NY,NX)=RFLC1P
  RC2PXS(NU(NY,NX),NY,NX)=RFLC2P
  RM1PXS(NU(NY,NX),NY,NX)=RFLM1P
  RH0BXB(NU(NY,NX),NY,NX)=RFLH0B
  RH3BXB(NU(NY,NX),NY,NX)=RFLH3B
  RF1BXB(NU(NY,NX),NY,NX)=RFLF1B
  RF2BXB(NU(NY,NX),NY,NX)=RFLF2B
  RC0BXB(NU(NY,NX),NY,NX)=RFLC0B
  RC1BXB(NU(NY,NX),NY,NX)=RFLC1B
  RC2BXB(NU(NY,NX),NY,NX)=RFLC2B
  RM1BXB(NU(NY,NX),NY,NX)=RFLM1B
  end subroutine MacToMicPoreSoluteAdvExchange
!------------------------------------------------------------------------------------------

  subroutine MicToMacPoreSoluteAdvExchange(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VFLW
!     begin_execution
  IF(VOLWM(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMIN1(AMAX1(-VFLWX,FINHM(M,NU(NY,NX),NY,NX)/VOLWM(M,NU(NY,NX),NY,NX)))
  ELSE
    VFLW=-VFLWX
  ENDIF
  RFLAL=VFLW*AZMAX1(ZAL2(NU(NY,NX),NY,NX))
  RFLFE=VFLW*AZMAX1(ZFE2(NU(NY,NX),NY,NX))
  RFLHY=VFLW*AZMAX1(ZHY2(NU(NY,NX),NY,NX))
  RFLCA=VFLW*AZMAX1(ZCA2(NU(NY,NX),NY,NX))
  RFLMG=VFLW*AZMAX1(ZMG2(NU(NY,NX),NY,NX))
  RFLNA=VFLW*AZMAX1(ZNA2(NU(NY,NX),NY,NX))
  RFLKA=VFLW*AZMAX1(ZKA2(NU(NY,NX),NY,NX))
  RFLOH=VFLW*AZMAX1(ZOH2(NU(NY,NX),NY,NX))
  RFLSO=VFLW*AZMAX1(ZSO42(NU(NY,NX),NY,NX))
  RFLCL=VFLW*AZMAX1(ZCL2(NU(NY,NX),NY,NX))
  RFLC3=VFLW*AZMAX1(ZCO32(NU(NY,NX),NY,NX))
  RFLHC=VFLW*AZMAX1(ZHCO32(NU(NY,NX),NY,NX))
  RFLAL1=VFLW*AZMAX1(ZAL12(NU(NY,NX),NY,NX))
  RFLAL2=VFLW*AZMAX1(ZAL22(NU(NY,NX),NY,NX))
  RFLAL3=VFLW*AZMAX1(ZAL32(NU(NY,NX),NY,NX))
  RFLAL4=VFLW*AZMAX1(ZAL42(NU(NY,NX),NY,NX))
  RFLALS=VFLW*AZMAX1(ZALS2(NU(NY,NX),NY,NX))
  RFLFE1=VFLW*AZMAX1(ZFE12(NU(NY,NX),NY,NX))
  RFLFE2=VFLW*AZMAX1(ZFE22(NU(NY,NX),NY,NX))
  RFLFE3=VFLW*AZMAX1(ZFE32(NU(NY,NX),NY,NX))
  RFLFE4=VFLW*AZMAX1(ZFE42(NU(NY,NX),NY,NX))
  RFLFES=VFLW*AZMAX1(ZFES2(NU(NY,NX),NY,NX))
  RFLCAO=VFLW*AZMAX1(ZCAO2(NU(NY,NX),NY,NX))
  RFLCAC=VFLW*AZMAX1(ZCAC2(NU(NY,NX),NY,NX))
  RFLCAH=VFLW*AZMAX1(ZCAH2(NU(NY,NX),NY,NX))
  RFLCAS=VFLW*AZMAX1(ZCAS2(NU(NY,NX),NY,NX))
  RFLMGO=VFLW*AZMAX1(ZMGO2(NU(NY,NX),NY,NX))
  RFLMGC=VFLW*AZMAX1(ZMGC2(NU(NY,NX),NY,NX))
  RFLMGH=VFLW*AZMAX1(ZMGH2(NU(NY,NX),NY,NX))
  RFLMGS=VFLW*AZMAX1(ZMGS2(NU(NY,NX),NY,NX))
  RFLNAC=VFLW*AZMAX1(ZNAC2(NU(NY,NX),NY,NX))
  RFLNAS=VFLW*AZMAX1(ZNAS2(NU(NY,NX),NY,NX))
  RFLKAS=VFLW*AZMAX1(ZKAS2(NU(NY,NX),NY,NX))
  RFLH0P=VFLW*AZMAX1(H0PO42(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH3P=VFLW*AZMAX1(H3PO42(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF1P=VFLW*AZMAX1(ZFE1P2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLF2P=VFLW*AZMAX1(ZFE2P2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC0P=VFLW*AZMAX1(ZCA0P2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC1P=VFLW*AZMAX1(ZCA1P2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLC2P=VFLW*AZMAX1(ZCA2P2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLM1P=VFLW*AZMAX1(ZMG1P2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  RFLH0B=VFLW*AZMAX1(H0POB2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLH3B=VFLW*AZMAX1(H3POB2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF1B=VFLW*AZMAX1(ZF1PB2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLF2B=VFLW*AZMAX1(ZF2PB2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC0B=VFLW*AZMAX1(ZC0PB2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC1B=VFLW*AZMAX1(ZC1PB2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLC2B=VFLW*AZMAX1(ZC2PB2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  RFLM1B=VFLW*AZMAX1(ZM1PB2(NU(NY,NX),NY,NX))*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)


  RALFXS(NU(NY,NX),NY,NX)=RFLAL
  RFEFXS(NU(NY,NX),NY,NX)=RFLFE
  RHYFXS(NU(NY,NX),NY,NX)=RFLHY
  RCAFXS(NU(NY,NX),NY,NX)=RFLCA
  RMGFXS(NU(NY,NX),NY,NX)=RFLMG
  RNAFXS(NU(NY,NX),NY,NX)=RFLNA
  RKAFXS(NU(NY,NX),NY,NX)=RFLKA
  ROHFXS(NU(NY,NX),NY,NX)=RFLOH
  RSOFXS(NU(NY,NX),NY,NX)=RFLSO
  RCLFXS(NU(NY,NX),NY,NX)=RFLCL
  RC3FXS(NU(NY,NX),NY,NX)=RFLC3
  RHCFXS(NU(NY,NX),NY,NX)=RFLHC
  RAL1XS(NU(NY,NX),NY,NX)=RFLAL1
  RAL2XS(NU(NY,NX),NY,NX)=RFLAL2
  RAL3XS(NU(NY,NX),NY,NX)=RFLAL3
  RAL4XS(NU(NY,NX),NY,NX)=RFLAL4
  RALSXS(NU(NY,NX),NY,NX)=RFLALS
  RFE1XS(NU(NY,NX),NY,NX)=RFLFE1
  RFE2XS(NU(NY,NX),NY,NX)=RFLFE2
  RFE3XS(NU(NY,NX),NY,NX)=RFLFE3
  RFE4XS(NU(NY,NX),NY,NX)=RFLFE4
  RFESXS(NU(NY,NX),NY,NX)=RFLFES
  RCAOXS(NU(NY,NX),NY,NX)=RFLCAO
  RCACXS(NU(NY,NX),NY,NX)=RFLCAC
  RCAHXS(NU(NY,NX),NY,NX)=RFLCAH
  RCASXS(NU(NY,NX),NY,NX)=RFLCAS
  RMGOXS(NU(NY,NX),NY,NX)=RFLMGO
  RMGCXS(NU(NY,NX),NY,NX)=RFLMGC
  RMGHXS(NU(NY,NX),NY,NX)=RFLMGH
  RMGSXS(NU(NY,NX),NY,NX)=RFLMGS
  RNACXS(NU(NY,NX),NY,NX)=RFLNAC
  RNASXS(NU(NY,NX),NY,NX)=RFLNAS
  RKASXS(NU(NY,NX),NY,NX)=RFLKAS
  RH0PXS(NU(NY,NX),NY,NX)=RFLH0P
  RH3PXS(NU(NY,NX),NY,NX)=RFLH3P
  RF1PXS(NU(NY,NX),NY,NX)=RFLF1P
  RF2PXS(NU(NY,NX),NY,NX)=RFLF2P
  RC0PXS(NU(NY,NX),NY,NX)=RFLC0P
  RC1PXS(NU(NY,NX),NY,NX)=RFLC1P
  RC2PXS(NU(NY,NX),NY,NX)=RFLC2P
  RM1PXS(NU(NY,NX),NY,NX)=RFLM1P
  RH0BXB(NU(NY,NX),NY,NX)=RFLH0B
  RH3BXB(NU(NY,NX),NY,NX)=RFLH3B
  RF1BXB(NU(NY,NX),NY,NX)=RFLF1B
  RF2BXB(NU(NY,NX),NY,NX)=RFLF2B
  RC0BXB(NU(NY,NX),NY,NX)=RFLC0B
  RC1BXB(NU(NY,NX),NY,NX)=RFLC1B
  RC2BXB(NU(NY,NX),NY,NX)=RFLC2B
  RM1BXB(NU(NY,NX),NY,NX)=RFLM1B
  end subroutine MicToMacPoreSoluteAdvExchange
!------------------------------------------------------------------------------------------

  subroutine ZeroMicMacPoreExchange
!
!     Description:
!
  implicit none

  RFLAL=0.0_r8
  RFLFE=0.0_r8
  RFLHY=0.0_r8
  RFLCA=0.0_r8
  RFLMG=0.0_r8
  RFLNA=0.0_r8
  RFLKA=0.0_r8
  RFLOH=0.0_r8
  RFLSO=0.0_r8
  RFLCL=0.0_r8
  RFLC3=0.0_r8
  RFLHC=0.0_r8
  RFLAL1=0.0_r8
  RFLAL2=0.0_r8
  RFLAL3=0.0_r8
  RFLAL4=0.0_r8
  RFLALS=0.0_r8
  RFLFE1=0.0_r8
  RFLFE2=0.0_r8
  RFLFE3=0.0_r8
  RFLFE4=0.0_r8
  RFLFES=0.0_r8
  RFLCAO=0.0_r8
  RFLCAC=0.0_r8
  RFLCAH=0.0_r8
  RFLCAS=0.0_r8
  RFLMGO=0.0_r8
  RFLMGC=0.0_r8
  RFLMGH=0.0_r8
  RFLMGS=0.0_r8
  RFLNAC=0.0_r8
  RFLNAS=0.0_r8
  RFLKAS=0.0_r8
  RFLH0P=0.0_r8
  RFLH3P=0.0_r8
  RFLF1P=0.0_r8
  RFLF2P=0.0_r8
  RFLC0P=0.0_r8
  RFLC1P=0.0_r8
  RFLC2P=0.0_r8
  RFLM1P=0.0_r8
  RFLH0B=0.0_r8
  RFLH3B=0.0_r8
  RFLF1B=0.0_r8
  RFLF2B=0.0_r8
  RFLC0B=0.0_r8
  RFLC1B=0.0_r8
  RFLC2B=0.0_r8
  RFLM1B=0.0_r8
  end subroutine ZeroMicMacPoreExchange
!------------------------------------------------------------------------------------------

  subroutine MacMicPoreSoluteDifusExchange(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VOLWHS,VOLWT,VOLWMNU
!     begin_execution
!
!     VOLWM,VOLWHM=micropore,macropore water volume
!     XFRS*VOLT=maximum macropore volume for solute transfer
!     DFV*=diffusive macropore-micropore solute transfer
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     *H2,*2=macropore,micropore solute content
  VOLWMNU=VOLWM(M,NU(NY,NX),NY,NX)
  VOLWHS=AMIN1(XFRS*VOLT(NU(NY,NX),NY,NX),VOLWHM(M,NU(NY,NX),NY,NX))
  VOLWT=VOLWM(M,NU(NY,NX),NY,NX)+VOLWHS
  DFVAL=XNPH*(AZMAX1(ZALH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZAL2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVFE=XNPH*(AZMAX1(ZFEH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZFE2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVHY=XNPH*(AZMAX1(ZHYH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZHY2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCA=XNPH*(AZMAX1(ZCCH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZCA2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVMG=XNPH*(AZMAX1(ZMAH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZMG2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVNA=XNPH*(AZMAX1(ZNAH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZNA2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVKA=XNPH*(AZMAX1(ZKAH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZKA2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVOH=XNPH*(AZMAX1(ZOHH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZOH2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVSO=XNPH*(AZMAX1(ZSO4H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZSO42(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCL=XNPH*(AZMAX1(ZCLH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZCL2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVC3=XNPH*(AZMAX1(ZCO3H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZCO32(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVHC=XNPH*(AZMAX1(ZHCOH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZHCO32(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVAL1=XNPH*(AZMAX1(ZAL1H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZAL12(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVAL2=XNPH*(AZMAX1(ZAL2H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZAL22(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVAL3=XNPH*(AZMAX1(ZAL3H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZAL32(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVAL4=XNPH*(AZMAX1(ZAL4H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZAL42(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVALS=XNPH*(AZMAX1(ZALSH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZALS2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVFE1=XNPH*(AZMAX1(ZFE1H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZFE12(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVF22=XNPH*(AZMAX1(ZFE2H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZFE22(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVFE3=XNPH*(AZMAX1(ZFE3H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZFE32(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVFE4=XNPH*(AZMAX1(ZFE4H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZFE42(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVFES=XNPH*(AZMAX1(ZFESH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZFES2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCAO=XNPH*(AZMAX1(ZCAOH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZCAO2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCAC=XNPH*(AZMAX1(ZCACH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZCAC2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCAH=XNPH*(AZMAX1(ZCAHH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZCAH2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVCAS=XNPH*(AZMAX1(ZCASH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZCAS2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVMGO=XNPH*(AZMAX1(ZMGOH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZMGO2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVMGC=XNPH*(AZMAX1(ZMGCH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZMGC2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVMGH=XNPH*(AZMAX1(ZMGHH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZMGH2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVMGS=XNPH*(AZMAX1(ZMGSH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZMGS2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVNAC=XNPH*(AZMAX1(ZNACH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZNAC2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVNAS=XNPH*(AZMAX1(ZNASH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZNAS2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVKAS=XNPH*(AZMAX1(ZKASH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZKAS2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVNAC=XNPH*(AZMAX1(ZNACH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZNAC2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT
  DFVH0P=XNPH*(AZMAX1(H0P4H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(H0PO42(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVH3P=XNPH*(AZMAX1(H3P4H2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(H3PO42(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVF1P=XNPH*(AZMAX1(ZF1PH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZFE1P2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVF2P=XNPH*(AZMAX1(ZF2PH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZFE2P2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC0P=XNPH*(AZMAX1(ZC0PH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZCA0P2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC1P=XNPH*(AZMAX1(ZC1PH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZCA1P2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVC2P=XNPH*(AZMAX1(ZC2PH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZCA2P2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVM1P=XNPH*(AZMAX1(ZM1PH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZMG1P2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,NU(NY,NX),NY,NX)
  DFVH0B=XNPH*(AZMAX1(H0PBH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(H0POB2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVH3B=XNPH*(AZMAX1(H3PBH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(H3POB2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVF1B=XNPH*(AZMAX1(ZF1BH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZF1PB2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVF2B=XNPH*(AZMAX1(ZF2BH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZF2PB2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC0B=XNPH*(AZMAX1(ZC0BH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZC0PB2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC1B=XNPH*(AZMAX1(ZC1BH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZC1PB2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVC2B=XNPH*(AZMAX1(ZC2BH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZC2PB2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)
  DFVM1B=XNPH*(AZMAX1(ZM1BH2(NU(NY,NX),NY,NX))*VOLWMNU-AZMAX1(ZM1PB2(NU(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,NU(NY,NX),NY,NX)


!
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     RFL*=convective flux between macropores and micropores
!     DFV*=diffusive solute flux between macropores and micropores
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  RALFXS(NU(NY,NX),NY,NX)=RALFXS(NU(NY,NX),NY,NX)+DFVAL
  RFEFXS(NU(NY,NX),NY,NX)=RFEFXS(NU(NY,NX),NY,NX)+DFVFE
  RHYFXS(NU(NY,NX),NY,NX)=RHYFXS(NU(NY,NX),NY,NX)+DFVHY
  RCAFXS(NU(NY,NX),NY,NX)=RCAFXS(NU(NY,NX),NY,NX)+DFVCA
  RMGFXS(NU(NY,NX),NY,NX)=RMGFXS(NU(NY,NX),NY,NX)+DFVMG
  RNAFXS(NU(NY,NX),NY,NX)=RNAFXS(NU(NY,NX),NY,NX)+DFVNA
  RKAFXS(NU(NY,NX),NY,NX)=RKAFXS(NU(NY,NX),NY,NX)+DFVKA
  ROHFXS(NU(NY,NX),NY,NX)=ROHFXS(NU(NY,NX),NY,NX)+DFVOH
  RSOFXS(NU(NY,NX),NY,NX)=RSOFXS(NU(NY,NX),NY,NX)+DFVSO
  RCLFXS(NU(NY,NX),NY,NX)=RCLFXS(NU(NY,NX),NY,NX)+DFVCL
  RC3FXS(NU(NY,NX),NY,NX)=RC3FXS(NU(NY,NX),NY,NX)+DFVC3
  RHCFXS(NU(NY,NX),NY,NX)=RHCFXS(NU(NY,NX),NY,NX)+DFVHC
  RAL1XS(NU(NY,NX),NY,NX)=RAL1XS(NU(NY,NX),NY,NX)+DFVAL1
  RAL2XS(NU(NY,NX),NY,NX)=RAL2XS(NU(NY,NX),NY,NX)+DFVAL2
  RAL3XS(NU(NY,NX),NY,NX)=RAL3XS(NU(NY,NX),NY,NX)+DFVAL3
  RAL4XS(NU(NY,NX),NY,NX)=RAL4XS(NU(NY,NX),NY,NX)+DFVAL4
  RALSXS(NU(NY,NX),NY,NX)=RALSXS(NU(NY,NX),NY,NX)+DFVALS
  RFE1XS(NU(NY,NX),NY,NX)=RFE1XS(NU(NY,NX),NY,NX)+DFVFE1
  RFE2XS(NU(NY,NX),NY,NX)=RFE2XS(NU(NY,NX),NY,NX)+DFVFE2
  RFE3XS(NU(NY,NX),NY,NX)=RFE3XS(NU(NY,NX),NY,NX)+DFVFE3
  RFE4XS(NU(NY,NX),NY,NX)=RFE4XS(NU(NY,NX),NY,NX)+DFVFE4
  RFESXS(NU(NY,NX),NY,NX)=RFESXS(NU(NY,NX),NY,NX)+DFVFES
  RCAOXS(NU(NY,NX),NY,NX)=RCAOXS(NU(NY,NX),NY,NX)+DFVCAO
  RCACXS(NU(NY,NX),NY,NX)=RCACXS(NU(NY,NX),NY,NX)+DFVCAC
  RCAHXS(NU(NY,NX),NY,NX)=RCAHXS(NU(NY,NX),NY,NX)+DFVCAH
  RCASXS(NU(NY,NX),NY,NX)=RCASXS(NU(NY,NX),NY,NX)+DFVCAS
  RMGOXS(NU(NY,NX),NY,NX)=RMGOXS(NU(NY,NX),NY,NX)+DFVMGO
  RMGCXS(NU(NY,NX),NY,NX)=RMGCXS(NU(NY,NX),NY,NX)+DFVMGC
  RMGHXS(NU(NY,NX),NY,NX)=RMGHXS(NU(NY,NX),NY,NX)+DFVMGH
  RMGSXS(NU(NY,NX),NY,NX)=RMGSXS(NU(NY,NX),NY,NX)+DFVMGS
  RNACXS(NU(NY,NX),NY,NX)=RNACXS(NU(NY,NX),NY,NX)+DFVNAC
  RNASXS(NU(NY,NX),NY,NX)=RNASXS(NU(NY,NX),NY,NX)+DFVNAS
  RKASXS(NU(NY,NX),NY,NX)=RKASXS(NU(NY,NX),NY,NX)+DFVKAS
  RH0PXS(NU(NY,NX),NY,NX)=RH0PXS(NU(NY,NX),NY,NX)+DFVH0P
  RH3PXS(NU(NY,NX),NY,NX)=RH3PXS(NU(NY,NX),NY,NX)+DFVH3P
  RF1PXS(NU(NY,NX),NY,NX)=RF1PXS(NU(NY,NX),NY,NX)+DFVF1P
  RF2PXS(NU(NY,NX),NY,NX)=RF2PXS(NU(NY,NX),NY,NX)+DFVF2P
  RC0PXS(NU(NY,NX),NY,NX)=RC0PXS(NU(NY,NX),NY,NX)+DFVC0P
  RC1PXS(NU(NY,NX),NY,NX)=RC1PXS(NU(NY,NX),NY,NX)+DFVC1P
  RC2PXS(NU(NY,NX),NY,NX)=RC2PXS(NU(NY,NX),NY,NX)+DFVC2P
  RM1PXS(NU(NY,NX),NY,NX)=RM1PXS(NU(NY,NX),NY,NX)+DFVM1P
  RH0BXB(NU(NY,NX),NY,NX)=RH0BXB(NU(NY,NX),NY,NX)+DFVH0B
  RH3BXB(NU(NY,NX),NY,NX)=RH3BXB(NU(NY,NX),NY,NX)+DFVH3B
  RF1BXB(NU(NY,NX),NY,NX)=RF1BXB(NU(NY,NX),NY,NX)+DFVF1B
  RF2BXB(NU(NY,NX),NY,NX)=RF2BXB(NU(NY,NX),NY,NX)+DFVF2B
  RC0BXB(NU(NY,NX),NY,NX)=RC0BXB(NU(NY,NX),NY,NX)+DFVC0B
  RC1BXB(NU(NY,NX),NY,NX)=RC1BXB(NU(NY,NX),NY,NX)+DFVC1B
  RC2BXB(NU(NY,NX),NY,NX)=RC2BXB(NU(NY,NX),NY,NX)+DFVC2B
  RM1BXB(NU(NY,NX),NY,NX)=RM1BXB(NU(NY,NX),NY,NX)+DFVM1B
  end subroutine MacMicPoreSoluteDifusExchange
!------------------------------------------------------------------------------------------

  subroutine ZeroMacMicPoreSoluteDifusExch
!
!     Description:
!
  implicit none
!     begin_execution
  DFVAL=0.0_r8
  DFVFE=0.0_r8
  DFVHY=0.0_r8
  DFVCA=0.0_r8
  DFVMG=0.0_r8
  DFVNA=0.0_r8
  DFVKA=0.0_r8
  DFVOH=0.0_r8
  DFVSO=0.0_r8
  DFVCL=0.0_r8
  DFVC3=0.0_r8
  DFVHC=0.0_r8
  DFVAL1=0.0_r8
  DFVAL2=0.0_r8
  DFVAL3=0.0_r8
  DFVAL4=0.0_r8
  DFVALS=0.0_r8
  DFVFE1=0.0_r8
  DFVFE2=0.0_r8
  DFVFE3=0.0_r8
  DFVFE4=0.0_r8
  DFVFES=0.0_r8
  DFVCAO=0.0_r8
  DFVCAC=0.0_r8
  DFVCAH=0.0_r8
  DFVCAS=0.0_r8
  DFVMGO=0.0_r8
  DFVMGC=0.0_r8
  DFVMGH=0.0_r8
  DFVMGS=0.0_r8
  DFVNAC=0.0_r8
  DFVNAS=0.0_r8
  DFVKAS=0.0_r8
  DFVH0P=0.0_r8
  DFVH3P=0.0_r8
  DFVF1P=0.0_r8
  DFVF2P=0.0_r8
  DFVC0P=0.0_r8
  DFVC1P=0.0_r8
  DFVC2P=0.0_r8
  DFVM1P=0.0_r8
  DFVH0B=0.0_r8
  DFVH3B=0.0_r8
  DFVF1B=0.0_r8
  DFVF2B=0.0_r8
  DFVC0B=0.0_r8
  DFVC1B=0.0_r8
  DFVC2B=0.0_r8
  DFVM1B=0.0_r8
  end subroutine ZeroMacMicPoreSoluteDifusExch
!------------------------------------------------------------------------------------------

  subroutine AccumHourlyMicMacPoreFlux(NY,NX)
!
!     Description:
!
      implicit none
      integer, intent(in) :: NY,NX
!     begin_execution
!
!     X*FXS=hourly convective + diffusive solute flux between macropores and micropores
!     R*FXS=total convective + diffusive solute flux between macropores and micropores
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  XALFXS(NU(NY,NX),NY,NX)=XALFXS(NU(NY,NX),NY,NX)+RALFXS(NU(NY,NX),NY,NX)
  XFEFXS(NU(NY,NX),NY,NX)=XFEFXS(NU(NY,NX),NY,NX)+RFEFXS(NU(NY,NX),NY,NX)
  XHYFXS(NU(NY,NX),NY,NX)=XHYFXS(NU(NY,NX),NY,NX)+RHYFXS(NU(NY,NX),NY,NX)
  XCAFXS(NU(NY,NX),NY,NX)=XCAFXS(NU(NY,NX),NY,NX)+RCAFXS(NU(NY,NX),NY,NX)
  XMGFXS(NU(NY,NX),NY,NX)=XMGFXS(NU(NY,NX),NY,NX)+RMGFXS(NU(NY,NX),NY,NX)
  XNAFXS(NU(NY,NX),NY,NX)=XNAFXS(NU(NY,NX),NY,NX)+RNAFXS(NU(NY,NX),NY,NX)
  XKAFXS(NU(NY,NX),NY,NX)=XKAFXS(NU(NY,NX),NY,NX)+RKAFXS(NU(NY,NX),NY,NX)
  XOHFXS(NU(NY,NX),NY,NX)=XOHFXS(NU(NY,NX),NY,NX)+ROHFXS(NU(NY,NX),NY,NX)
  XSOFXS(NU(NY,NX),NY,NX)=XSOFXS(NU(NY,NX),NY,NX)+RSOFXS(NU(NY,NX),NY,NX)
  XCLFXS(NU(NY,NX),NY,NX)=XCLFXS(NU(NY,NX),NY,NX)+RCLFXS(NU(NY,NX),NY,NX)
  XC3FXS(NU(NY,NX),NY,NX)=XC3FXS(NU(NY,NX),NY,NX)+RC3FXS(NU(NY,NX),NY,NX)
  XHCFXS(NU(NY,NX),NY,NX)=XHCFXS(NU(NY,NX),NY,NX)+RHCFXS(NU(NY,NX),NY,NX)
  XAL1XS(NU(NY,NX),NY,NX)=XAL1XS(NU(NY,NX),NY,NX)+RAL1XS(NU(NY,NX),NY,NX)
  XAL2XS(NU(NY,NX),NY,NX)=XAL2XS(NU(NY,NX),NY,NX)+RAL2XS(NU(NY,NX),NY,NX)
  XAL3XS(NU(NY,NX),NY,NX)=XAL3XS(NU(NY,NX),NY,NX)+RAL3XS(NU(NY,NX),NY,NX)
  XAL4XS(NU(NY,NX),NY,NX)=XAL4XS(NU(NY,NX),NY,NX)+RAL4XS(NU(NY,NX),NY,NX)
  XALSXS(NU(NY,NX),NY,NX)=XALSXS(NU(NY,NX),NY,NX)+RALSXS(NU(NY,NX),NY,NX)
  XFE1XS(NU(NY,NX),NY,NX)=XFE1XS(NU(NY,NX),NY,NX)+RFE1XS(NU(NY,NX),NY,NX)
  XFE2XS(NU(NY,NX),NY,NX)=XFE2XS(NU(NY,NX),NY,NX)+RFE2XS(NU(NY,NX),NY,NX)
  XFE3XS(NU(NY,NX),NY,NX)=XFE3XS(NU(NY,NX),NY,NX)+RFE3XS(NU(NY,NX),NY,NX)
  XFE4XS(NU(NY,NX),NY,NX)=XFE4XS(NU(NY,NX),NY,NX)+RFE4XS(NU(NY,NX),NY,NX)
  XFESXS(NU(NY,NX),NY,NX)=XFESXS(NU(NY,NX),NY,NX)+RFESXS(NU(NY,NX),NY,NX)
  XCAOXS(NU(NY,NX),NY,NX)=XCAOXS(NU(NY,NX),NY,NX)+RCAOXS(NU(NY,NX),NY,NX)
  XCACXS(NU(NY,NX),NY,NX)=XCACXS(NU(NY,NX),NY,NX)+RCACXS(NU(NY,NX),NY,NX)
  XCAHXS(NU(NY,NX),NY,NX)=XCAHXS(NU(NY,NX),NY,NX)+RCAHXS(NU(NY,NX),NY,NX)
  XCASXS(NU(NY,NX),NY,NX)=XCASXS(NU(NY,NX),NY,NX)+RCASXS(NU(NY,NX),NY,NX)
  XMGOXS(NU(NY,NX),NY,NX)=XMGOXS(NU(NY,NX),NY,NX)+RMGOXS(NU(NY,NX),NY,NX)
  XMGCXS(NU(NY,NX),NY,NX)=XMGCXS(NU(NY,NX),NY,NX)+RMGCXS(NU(NY,NX),NY,NX)
  XMGHXS(NU(NY,NX),NY,NX)=XMGHXS(NU(NY,NX),NY,NX)+RMGHXS(NU(NY,NX),NY,NX)
  XMGSXS(NU(NY,NX),NY,NX)=XMGSXS(NU(NY,NX),NY,NX)+RMGSXS(NU(NY,NX),NY,NX)
  XNACXS(NU(NY,NX),NY,NX)=XNACXS(NU(NY,NX),NY,NX)+RNACXS(NU(NY,NX),NY,NX)
  XNASXS(NU(NY,NX),NY,NX)=XNASXS(NU(NY,NX),NY,NX)+RNASXS(NU(NY,NX),NY,NX)
  XKASXS(NU(NY,NX),NY,NX)=XKASXS(NU(NY,NX),NY,NX)+RKASXS(NU(NY,NX),NY,NX)
  XH0PXS(NU(NY,NX),NY,NX)=XH0PXS(NU(NY,NX),NY,NX)+RH0PXS(NU(NY,NX),NY,NX)
  XH3PXS(NU(NY,NX),NY,NX)=XH3PXS(NU(NY,NX),NY,NX)+RH3PXS(NU(NY,NX),NY,NX)
  XF1PXS(NU(NY,NX),NY,NX)=XF1PXS(NU(NY,NX),NY,NX)+RF1PXS(NU(NY,NX),NY,NX)
  XF2PXS(NU(NY,NX),NY,NX)=XF2PXS(NU(NY,NX),NY,NX)+RF2PXS(NU(NY,NX),NY,NX)
  XC0PXS(NU(NY,NX),NY,NX)=XC0PXS(NU(NY,NX),NY,NX)+RC0PXS(NU(NY,NX),NY,NX)
  XC1PXS(NU(NY,NX),NY,NX)=XC1PXS(NU(NY,NX),NY,NX)+RC1PXS(NU(NY,NX),NY,NX)
  XC2PXS(NU(NY,NX),NY,NX)=XC2PXS(NU(NY,NX),NY,NX)+RC2PXS(NU(NY,NX),NY,NX)
  XM1PXS(NU(NY,NX),NY,NX)=XM1PXS(NU(NY,NX),NY,NX)+RM1PXS(NU(NY,NX),NY,NX)
  XH0BXB(NU(NY,NX),NY,NX)=XH0BXB(NU(NY,NX),NY,NX)+RH0BXB(NU(NY,NX),NY,NX)
  XH3BXB(NU(NY,NX),NY,NX)=XH3BXB(NU(NY,NX),NY,NX)+RH3BXB(NU(NY,NX),NY,NX)
  XF1BXB(NU(NY,NX),NY,NX)=XF1BXB(NU(NY,NX),NY,NX)+RF1BXB(NU(NY,NX),NY,NX)
  XF2BXB(NU(NY,NX),NY,NX)=XF2BXB(NU(NY,NX),NY,NX)+RF2BXB(NU(NY,NX),NY,NX)
  XC0BXB(NU(NY,NX),NY,NX)=XC0BXB(NU(NY,NX),NY,NX)+RC0BXB(NU(NY,NX),NY,NX)
  XC1BXB(NU(NY,NX),NY,NX)=XC1BXB(NU(NY,NX),NY,NX)+RC1BXB(NU(NY,NX),NY,NX)
  XC2BXB(NU(NY,NX),NY,NX)=XC2BXB(NU(NY,NX),NY,NX)+RC2BXB(NU(NY,NX),NY,NX)
  XM1BXB(NU(NY,NX),NY,NX)=XM1BXB(NU(NY,NX),NY,NX)+RM1BXB(NU(NY,NX),NY,NX)
  end subroutine AccumHourlyMicMacPoreFlux
!------------------------------------------------------------------------------------------

  subroutine SoluteFluxBySurfaceOutflow(M,N1,N2)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N1,N2
  real(r8) :: VFLW

!     begin_execution

  IF(VOLWM(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
    VFLW=AMIN1(VFLWX,QRM(M,N2,N1)/VOLWM(M,0,N2,N1))
  ELSE
    VFLW=VFLWX
  ENDIF
  RQRAL0(N2,N1)=VFLW*AZMAX1(ZAL2(0,N2,N1))
  RQRFE0(N2,N1)=VFLW*AZMAX1(ZFE2(0,N2,N1))
  RQRHY0(N2,N1)=VFLW*AZMAX1(ZHY2(0,N2,N1))
  RQRCA0(N2,N1)=VFLW*AZMAX1(ZCA2(0,N2,N1))
  RQRMG0(N2,N1)=VFLW*AZMAX1(ZMG2(0,N2,N1))
  RQRNA0(N2,N1)=VFLW*AZMAX1(ZNA2(0,N2,N1))
  RQRKA0(N2,N1)=VFLW*AZMAX1(ZKA2(0,N2,N1))
  RQROH0(N2,N1)=VFLW*AZMAX1(ZOH2(0,N2,N1))
  RQRSO0(N2,N1)=VFLW*AZMAX1(ZSO42(0,N2,N1))
  RQRCL0(N2,N1)=VFLW*AZMAX1(ZCL2(0,N2,N1))
  RQRC30(N2,N1)=VFLW*AZMAX1(ZCO32(0,N2,N1))
  RQRHC0(N2,N1)=VFLW*AZMAX1(ZHCO32(0,N2,N1))
  RQRAL10(N2,N1)=VFLW*AZMAX1(ZAL12(0,N2,N1))
  RQRAL20(N2,N1)=VFLW*AZMAX1(ZAL22(0,N2,N1))
  RQRAL30(N2,N1)=VFLW*AZMAX1(ZAL32(0,N2,N1))
  RQRAL40(N2,N1)=VFLW*AZMAX1(ZAL42(0,N2,N1))
  RQRALS0(N2,N1)=VFLW*AZMAX1(ZALS2(0,N2,N1))
  RQRFE10(N2,N1)=VFLW*AZMAX1(ZFE12(0,N2,N1))
  RQRFE20(N2,N1)=VFLW*AZMAX1(ZFE22(0,N2,N1))
  RQRFE30(N2,N1)=VFLW*AZMAX1(ZFE32(0,N2,N1))
  RQRFE40(N2,N1)=VFLW*AZMAX1(ZFE42(0,N2,N1))
  RQRFES0(N2,N1)=VFLW*AZMAX1(ZFES2(0,N2,N1))
  RQRCAO0(N2,N1)=VFLW*AZMAX1(ZCAO2(0,N2,N1))
  RQRCAC0(N2,N1)=VFLW*AZMAX1(ZCAC2(0,N2,N1))
  RQRCAH0(N2,N1)=VFLW*AZMAX1(ZCAH2(0,N2,N1))
  RQRCAS0(N2,N1)=VFLW*AZMAX1(ZCAS2(0,N2,N1))
  RQRMGO0(N2,N1)=VFLW*AZMAX1(ZMGO2(0,N2,N1))
  RQRMGC0(N2,N1)=VFLW*AZMAX1(ZMGC2(0,N2,N1))
  RQRMGH0(N2,N1)=VFLW*AZMAX1(ZMGH2(0,N2,N1))
  RQRMGS0(N2,N1)=VFLW*AZMAX1(ZMGS2(0,N2,N1))
  RQRNAC0(N2,N1)=VFLW*AZMAX1(ZNAC2(0,N2,N1))
  RQRNAS0(N2,N1)=VFLW*AZMAX1(ZNAS2(0,N2,N1))
  RQRKAS0(N2,N1)=VFLW*AZMAX1(ZKAS2(0,N2,N1))
  RQRH0P0(N2,N1)=VFLW*AZMAX1(H0PO42(0,N2,N1))
  RQRH3P0(N2,N1)=VFLW*AZMAX1(H3PO42(0,N2,N1))
  RQRF1P0(N2,N1)=VFLW*AZMAX1(ZFE1P2(0,N2,N1))
  RQRF2P0(N2,N1)=VFLW*AZMAX1(ZFE2P2(0,N2,N1))
  RQRC0P0(N2,N1)=VFLW*AZMAX1(ZCA0P2(0,N2,N1))
  RQRC1P0(N2,N1)=VFLW*AZMAX1(ZCA1P2(0,N2,N1))
  RQRC2P0(N2,N1)=VFLW*AZMAX1(ZCA2P2(0,N2,N1))
  RQRM1P0(N2,N1)=VFLW*AZMAX1(ZMG1P2(0,N2,N1))
  end subroutine SoluteFluxBySurfaceOutflow
!------------------------------------------------------------------------------------------

  subroutine ZeroSoluteFluxbySurfaceflow(N1,N2)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N1,N2
!     begin_execution

  RQRAL0(N2,N1)=0.0_r8
  RQRFE0(N2,N1)=0.0_r8
  RQRHY0(N2,N1)=0.0_r8
  RQRCA0(N2,N1)=0.0_r8
  RQRMG0(N2,N1)=0.0_r8
  RQRNA0(N2,N1)=0.0_r8
  RQRKA0(N2,N1)=0.0_r8
  RQROH0(N2,N1)=0.0_r8
  RQRSO0(N2,N1)=0.0_r8
  RQRCL0(N2,N1)=0.0_r8
  RQRC30(N2,N1)=0.0_r8
  RQRHC0(N2,N1)=0.0_r8
  RQRAL10(N2,N1)=0.0_r8
  RQRAL20(N2,N1)=0.0_r8
  RQRAL30(N2,N1)=0.0_r8
  RQRAL40(N2,N1)=0.0_r8
  RQRALS0(N2,N1)=0.0_r8
  RQRFE10(N2,N1)=0.0_r8
  RQRFE20(N2,N1)=0.0_r8
  RQRFE30(N2,N1)=0.0_r8
  RQRFE40(N2,N1)=0.0_r8
  RQRFES0(N2,N1)=0.0_r8
  RQRCAO0(N2,N1)=0.0_r8
  RQRCAC0(N2,N1)=0.0_r8
  RQRCAH0(N2,N1)=0.0_r8
  RQRCAS0(N2,N1)=0.0_r8
  RQRMGO0(N2,N1)=0.0_r8
  RQRMGC0(N2,N1)=0.0_r8
  RQRMGH0(N2,N1)=0.0_r8
  RQRMGS0(N2,N1)=0.0_r8
  RQRNAC0(N2,N1)=0.0_r8
  RQRNAS0(N2,N1)=0.0_r8
  RQRKAS0(N2,N1)=0.0_r8
  RQRH0P0(N2,N1)=0.0_r8
  RQRH3P0(N2,N1)=0.0_r8
  RQRF1P0(N2,N1)=0.0_r8
  RQRF2P0(N2,N1)=0.0_r8
  RQRC0P0(N2,N1)=0.0_r8
  RQRC1P0(N2,N1)=0.0_r8
  RQRC2P0(N2,N1)=0.0_r8
  RQRM1P0(N2,N1)=0.0_r8
  end subroutine ZeroSoluteFluxbySurfaceflow
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInSurfNeighbors(M,N1,N2,NY,NX,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N1,N2,NY,NX,NHW,NHE,NVN,NVS
  integer :: N,NN,N4,N5,N4B,N5B
  real(r8) :: FQRM,VFLW
!     begin_execution
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
  DO N=1,2
    DO  NN=1,2
      IF(N.EQ.1)THEN
        IF(NX.EQ.NHE.AND.NN.EQ.1.OR.NX.EQ.NHW.AND.NN.EQ.2)THEN
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N4B=NX-1
          N5B=NY
        ENDIF
      ELSEIF(N.EQ.2)THEN
        IF(NY.EQ.NVS.AND.NN.EQ.1.OR.NY.EQ.NVN.AND.NN.EQ.2)THEN
          cycle
        ELSE
          N4=NX
          N5=NY+1
          N4B=NX
          N5B=NY-1
        ENDIF
      ENDIF
!
!     QRM=runoff from watsub.f
!     RQR*=solute in runoff
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
      IF(QRM(M,N2,N1).GT.ZEROS(N2,N1))THEN
        IF(NN.EQ.1)THEN
          FQRM=QRMN(M,N,2,N5,N4)/QRM(M,N2,N1)
          RQRAL(N,2,N5,N4)=RQRAL0(N2,N1)*FQRM
          RQRFE(N,2,N5,N4)=RQRFE0(N2,N1)*FQRM
          RQRHY(N,2,N5,N4)=RQRHY0(N2,N1)*FQRM
          RQRCA(N,2,N5,N4)=RQRCA0(N2,N1)*FQRM
          RQRMG(N,2,N5,N4)=RQRMG0(N2,N1)*FQRM
          RQRNA(N,2,N5,N4)=RQRNA0(N2,N1)*FQRM
          RQRKA(N,2,N5,N4)=RQRKA0(N2,N1)*FQRM
          RQROH(N,2,N5,N4)=RQROH0(N2,N1)*FQRM
          RQRSO(N,2,N5,N4)=RQRSO0(N2,N1)*FQRM
          RQRCL(N,2,N5,N4)=RQRCL0(N2,N1)*FQRM
          RQRC3(N,2,N5,N4)=RQRC30(N2,N1)*FQRM
          RQRHC(N,2,N5,N4)=RQRHC0(N2,N1)*FQRM
          RQRAL1(N,2,N5,N4)=RQRAL10(N2,N1)*FQRM
          RQRAL2(N,2,N5,N4)=RQRAL20(N2,N1)*FQRM
          RQRAL3(N,2,N5,N4)=RQRAL30(N2,N1)*FQRM
          RQRAL4(N,2,N5,N4)=RQRAL40(N2,N1)*FQRM
          RQRALS(N,2,N5,N4)=RQRALS0(N2,N1)*FQRM
          RQRFE1(N,2,N5,N4)=RQRFE10(N2,N1)*FQRM
          RQRFE2(N,2,N5,N4)=RQRFE20(N2,N1)*FQRM
          RQRFE3(N,2,N5,N4)=RQRFE30(N2,N1)*FQRM
          RQRFE4(N,2,N5,N4)=RQRFE40(N2,N1)*FQRM
          RQRFES(N,2,N5,N4)=RQRFES0(N2,N1)*FQRM
          RQRCAO(N,2,N5,N4)=RQRCAO0(N2,N1)*FQRM
          RQRCAC(N,2,N5,N4)=RQRCAC0(N2,N1)*FQRM
          RQRCAH(N,2,N5,N4)=RQRCAH0(N2,N1)*FQRM
          RQRCAS(N,2,N5,N4)=RQRCAS0(N2,N1)*FQRM
          RQRMGO(N,2,N5,N4)=RQRMGO0(N2,N1)*FQRM
          RQRMGC(N,2,N5,N4)=RQRMGC0(N2,N1)*FQRM
          RQRMGH(N,2,N5,N4)=RQRMGH0(N2,N1)*FQRM
          RQRMGS(N,2,N5,N4)=RQRMGS0(N2,N1)*FQRM
          RQRNAC(N,2,N5,N4)=RQRNAC0(N2,N1)*FQRM
          RQRNAS(N,2,N5,N4)=RQRNAS0(N2,N1)*FQRM
          RQRKAS(N,2,N5,N4)=RQRKAS0(N2,N1)*FQRM
          RQRH0P(N,2,N5,N4)=RQRH0P0(N2,N1)*FQRM
          RQRH3P(N,2,N5,N4)=RQRH3P0(N2,N1)*FQRM
          RQRF1P(N,2,N5,N4)=RQRF1P0(N2,N1)*FQRM
          RQRF2P(N,2,N5,N4)=RQRF2P0(N2,N1)*FQRM
          RQRC0P(N,2,N5,N4)=RQRC0P0(N2,N1)*FQRM
          RQRC1P(N,2,N5,N4)=RQRC1P0(N2,N1)*FQRM
          RQRC2P(N,2,N5,N4)=RQRC2P0(N2,N1)*FQRM
          RQRM1P(N,2,N5,N4)=RQRM1P0(N2,N1)*FQRM
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQR*=hourly solute in runoff
!     RQR*=solute in runoff
!
          XQRAL(N,2,N5,N4)=XQRAL(N,2,N5,N4)+RQRAL(N,2,N5,N4)
          XQRFE(N,2,N5,N4)=XQRFE(N,2,N5,N4)+RQRFE(N,2,N5,N4)
          XQRHY(N,2,N5,N4)=XQRHY(N,2,N5,N4)+RQRHY(N,2,N5,N4)
          XQRCA(N,2,N5,N4)=XQRCA(N,2,N5,N4)+RQRCA(N,2,N5,N4)
          XQRMG(N,2,N5,N4)=XQRMG(N,2,N5,N4)+RQRMG(N,2,N5,N4)
          XQRNA(N,2,N5,N4)=XQRNA(N,2,N5,N4)+RQRNA(N,2,N5,N4)
          XQRKA(N,2,N5,N4)=XQRKA(N,2,N5,N4)+RQRKA(N,2,N5,N4)
          XQROH(N,2,N5,N4)=XQROH(N,2,N5,N4)+RQROH(N,2,N5,N4)
          XQRSO(N,2,N5,N4)=XQRSO(N,2,N5,N4)+RQRSO(N,2,N5,N4)
          XQRCL(N,2,N5,N4)=XQRCL(N,2,N5,N4)+RQRCL(N,2,N5,N4)
          XQRC3(N,2,N5,N4)=XQRC3(N,2,N5,N4)+RQRC3(N,2,N5,N4)
          XQRHC(N,2,N5,N4)=XQRHC(N,2,N5,N4)+RQRHC(N,2,N5,N4)
          XQRAL1(N,2,N5,N4)=XQRAL1(N,2,N5,N4)+RQRAL1(N,2,N5,N4)
          XQRAL2(N,2,N5,N4)=XQRAL2(N,2,N5,N4)+RQRAL2(N,2,N5,N4)
          XQRAL3(N,2,N5,N4)=XQRAL3(N,2,N5,N4)+RQRAL3(N,2,N5,N4)
          XQRAL4(N,2,N5,N4)=XQRAL4(N,2,N5,N4)+RQRAL4(N,2,N5,N4)
          XQRALS(N,2,N5,N4)=XQRALS(N,2,N5,N4)+RQRALS(N,2,N5,N4)
          XQRFE1(N,2,N5,N4)=XQRFE1(N,2,N5,N4)+RQRFE1(N,2,N5,N4)
          XQRFE2(N,2,N5,N4)=XQRFE2(N,2,N5,N4)+RQRFE2(N,2,N5,N4)
          XQRFE3(N,2,N5,N4)=XQRFE3(N,2,N5,N4)+RQRFE3(N,2,N5,N4)
          XQRFE4(N,2,N5,N4)=XQRFE4(N,2,N5,N4)+RQRFE4(N,2,N5,N4)
          XQRFES(N,2,N5,N4)=XQRFES(N,2,N5,N4)+RQRFES(N,2,N5,N4)
          XQRCAO(N,2,N5,N4)=XQRCAO(N,2,N5,N4)+RQRCAO(N,2,N5,N4)
          XQRCAC(N,2,N5,N4)=XQRCAC(N,2,N5,N4)+RQRCAC(N,2,N5,N4)
          XQRCAH(N,2,N5,N4)=XQRCAH(N,2,N5,N4)+RQRCAH(N,2,N5,N4)
          XQRCAS(N,2,N5,N4)=XQRCAS(N,2,N5,N4)+RQRCAS(N,2,N5,N4)
          XQRMGO(N,2,N5,N4)=XQRMGO(N,2,N5,N4)+RQRMGO(N,2,N5,N4)
          XQRMGC(N,2,N5,N4)=XQRMGC(N,2,N5,N4)+RQRMGC(N,2,N5,N4)
          XQRMGH(N,2,N5,N4)=XQRMGH(N,2,N5,N4)+RQRMGH(N,2,N5,N4)
          XQRMGS(N,2,N5,N4)=XQRMGS(N,2,N5,N4)+RQRMGS(N,2,N5,N4)
          XQRNAC(N,2,N5,N4)=XQRNAC(N,2,N5,N4)+RQRNAC(N,2,N5,N4)
          XQRNAS(N,2,N5,N4)=XQRNAS(N,2,N5,N4)+RQRNAS(N,2,N5,N4)
          XQRKAS(N,2,N5,N4)=XQRKAS(N,2,N5,N4)+RQRKAS(N,2,N5,N4)
          XQRH0P(N,2,N5,N4)=XQRH0P(N,2,N5,N4)+RQRH0P(N,2,N5,N4)
          XQRH3P(N,2,N5,N4)=XQRH3P(N,2,N5,N4)+RQRH3P(N,2,N5,N4)
          XQRF1P(N,2,N5,N4)=XQRF1P(N,2,N5,N4)+RQRF1P(N,2,N5,N4)
          XQRF2P(N,2,N5,N4)=XQRF2P(N,2,N5,N4)+RQRF2P(N,2,N5,N4)
          XQRC0P(N,2,N5,N4)=XQRC0P(N,2,N5,N4)+RQRC0P(N,2,N5,N4)
          XQRC1P(N,2,N5,N4)=XQRC1P(N,2,N5,N4)+RQRC1P(N,2,N5,N4)
          XQRC2P(N,2,N5,N4)=XQRC2P(N,2,N5,N4)+RQRC2P(N,2,N5,N4)
          XQRM1P(N,2,N5,N4)=XQRM1P(N,2,N5,N4)+RQRM1P(N,2,N5,N4)
        ELSE
          RQRAL(N,2,N5,N4)=0.0_r8
          RQRFE(N,2,N5,N4)=0.0_r8
          RQRHY(N,2,N5,N4)=0.0_r8
          RQRCA(N,2,N5,N4)=0.0_r8
          RQRMG(N,2,N5,N4)=0.0_r8
          RQRNA(N,2,N5,N4)=0.0_r8
          RQRKA(N,2,N5,N4)=0.0_r8
          RQROH(N,2,N5,N4)=0.0_r8
          RQRSO(N,2,N5,N4)=0.0_r8
          RQRCL(N,2,N5,N4)=0.0_r8
          RQRC3(N,2,N5,N4)=0.0_r8
          RQRHC(N,2,N5,N4)=0.0_r8
          RQRAL1(N,2,N5,N4)=0.0_r8
          RQRAL2(N,2,N5,N4)=0.0_r8
          RQRAL3(N,2,N5,N4)=0.0_r8
          RQRAL4(N,2,N5,N4)=0.0_r8
          RQRALS(N,2,N5,N4)=0.0_r8
          RQRFE1(N,2,N5,N4)=0.0_r8
          RQRFE2(N,2,N5,N4)=0.0_r8
          RQRFE3(N,2,N5,N4)=0.0_r8
          RQRFE4(N,2,N5,N4)=0.0_r8
          RQRFES(N,2,N5,N4)=0.0_r8
          RQRCAO(N,2,N5,N4)=0.0_r8
          RQRCAC(N,2,N5,N4)=0.0_r8
          RQRCAH(N,2,N5,N4)=0.0_r8
          RQRCAS(N,2,N5,N4)=0.0_r8
          RQRMGO(N,2,N5,N4)=0.0_r8
          RQRMGC(N,2,N5,N4)=0.0_r8
          RQRMGH(N,2,N5,N4)=0.0_r8
          RQRMGS(N,2,N5,N4)=0.0_r8
          RQRNAC(N,2,N5,N4)=0.0_r8
          RQRNAS(N,2,N5,N4)=0.0_r8
          RQRKAS(N,2,N5,N4)=0.0_r8
          RQRH0P(N,2,N5,N4)=0.0_r8
          RQRH3P(N,2,N5,N4)=0.0_r8
          RQRF1P(N,2,N5,N4)=0.0_r8
          RQRF2P(N,2,N5,N4)=0.0_r8
          RQRC0P(N,2,N5,N4)=0.0_r8
          RQRC1P(N,2,N5,N4)=0.0_r8
          RQRC2P(N,2,N5,N4)=0.0_r8
          RQRM1P(N,2,N5,N4)=0.0_r8
        ENDIF
!
!     IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
!
        IF(NN.EQ.2)THEN
          IF(N4B.GT.0.AND.N5B.GT.0)THEN
            FQRM=QRMN(M,N,1,N5B,N4B)/QRM(M,N2,N1)
            RQRAL(N,1,N5B,N4B)=RQRAL0(N2,N1)*FQRM
            RQRFE(N,1,N5B,N4B)=RQRFE0(N2,N1)*FQRM
            RQRHY(N,1,N5B,N4B)=RQRHY0(N2,N1)*FQRM
            RQRCA(N,1,N5B,N4B)=RQRCA0(N2,N1)*FQRM
            RQRMG(N,1,N5B,N4B)=RQRMG0(N2,N1)*FQRM
            RQRNA(N,1,N5B,N4B)=RQRNA0(N2,N1)*FQRM
            RQRKA(N,1,N5B,N4B)=RQRKA0(N2,N1)*FQRM
            RQROH(N,1,N5B,N4B)=RQROH0(N2,N1)*FQRM
            RQRSO(N,1,N5B,N4B)=RQRSO0(N2,N1)*FQRM
            RQRCL(N,1,N5B,N4B)=RQRCL0(N2,N1)*FQRM
            RQRC3(N,1,N5B,N4B)=RQRC30(N2,N1)*FQRM
            RQRHC(N,1,N5B,N4B)=RQRHC0(N2,N1)*FQRM
            RQRAL1(N,1,N5B,N4B)=RQRAL10(N2,N1)*FQRM
            RQRAL2(N,1,N5B,N4B)=RQRAL20(N2,N1)*FQRM
            RQRAL3(N,1,N5B,N4B)=RQRAL30(N2,N1)*FQRM
            RQRAL4(N,1,N5B,N4B)=RQRAL40(N2,N1)*FQRM
            RQRALS(N,1,N5B,N4B)=RQRALS0(N2,N1)*FQRM
            RQRFE1(N,1,N5B,N4B)=RQRFE10(N2,N1)*FQRM
            RQRFE2(N,1,N5B,N4B)=RQRFE20(N2,N1)*FQRM
            RQRFE3(N,1,N5B,N4B)=RQRFE30(N2,N1)*FQRM
            RQRFE4(N,1,N5B,N4B)=RQRFE40(N2,N1)*FQRM
            RQRFES(N,1,N5B,N4B)=RQRFES0(N2,N1)*FQRM
            RQRCAO(N,1,N5B,N4B)=RQRCAO0(N2,N1)*FQRM
            RQRCAC(N,1,N5B,N4B)=RQRCAC0(N2,N1)*FQRM
            RQRCAH(N,1,N5B,N4B)=RQRCAH0(N2,N1)*FQRM
            RQRCAS(N,1,N5B,N4B)=RQRCAS0(N2,N1)*FQRM
            RQRMGO(N,1,N5B,N4B)=RQRMGO0(N2,N1)*FQRM
            RQRMGC(N,1,N5B,N4B)=RQRMGC0(N2,N1)*FQRM
            RQRMGH(N,1,N5B,N4B)=RQRMGH0(N2,N1)*FQRM
            RQRMGS(N,1,N5B,N4B)=RQRMGS0(N2,N1)*FQRM
            RQRNAC(N,1,N5B,N4B)=RQRNAC0(N2,N1)*FQRM
            RQRNAS(N,1,N5B,N4B)=RQRNAS0(N2,N1)*FQRM
            RQRKAS(N,1,N5B,N4B)=RQRKAS0(N2,N1)*FQRM
            RQRH0P(N,1,N5B,N4B)=RQRH0P0(N2,N1)*FQRM
            RQRH3P(N,1,N5B,N4B)=RQRH3P0(N2,N1)*FQRM
            RQRF1P(N,1,N5B,N4B)=RQRF1P0(N2,N1)*FQRM
            RQRF2P(N,1,N5B,N4B)=RQRF2P0(N2,N1)*FQRM
            RQRC0P(N,1,N5B,N4B)=RQRC0P0(N2,N1)*FQRM
            RQRC1P(N,1,N5B,N4B)=RQRC1P0(N2,N1)*FQRM
            RQRC2P(N,1,N5B,N4B)=RQRC2P0(N2,N1)*FQRM
            RQRM1P(N,1,N5B,N4B)=RQRM1P0(N2,N1)*FQRM
      !
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQR*=hourly solute in runoff
!     RQR*=solute in runoff
!
            XQRAL(N,1,N5B,N4B)=XQRAL(N,1,N5B,N4B)+RQRAL(N,1,N5B,N4B)
            XQRFE(N,1,N5B,N4B)=XQRFE(N,1,N5B,N4B)+RQRFE(N,1,N5B,N4B)
            XQRHY(N,1,N5B,N4B)=XQRHY(N,1,N5B,N4B)+RQRHY(N,1,N5B,N4B)
            XQRCA(N,1,N5B,N4B)=XQRCA(N,1,N5B,N4B)+RQRCA(N,1,N5B,N4B)
            XQRMG(N,1,N5B,N4B)=XQRMG(N,1,N5B,N4B)+RQRMG(N,1,N5B,N4B)
            XQRNA(N,1,N5B,N4B)=XQRNA(N,1,N5B,N4B)+RQRNA(N,1,N5B,N4B)
            XQRKA(N,1,N5B,N4B)=XQRKA(N,1,N5B,N4B)+RQRKA(N,1,N5B,N4B)
            XQROH(N,1,N5B,N4B)=XQROH(N,1,N5B,N4B)+RQROH(N,1,N5B,N4B)
            XQRSO(N,1,N5B,N4B)=XQRSO(N,1,N5B,N4B)+RQRSO(N,1,N5B,N4B)
            XQRCL(N,1,N5B,N4B)=XQRCL(N,1,N5B,N4B)+RQRCL(N,1,N5B,N4B)
            XQRC3(N,1,N5B,N4B)=XQRC3(N,1,N5B,N4B)+RQRC3(N,1,N5B,N4B)
            XQRHC(N,1,N5B,N4B)=XQRHC(N,1,N5B,N4B)+RQRHC(N,1,N5B,N4B)
            XQRAL1(N,1,N5B,N4B)=XQRAL1(N,1,N5B,N4B)+RQRAL1(N,1,N5B,N4B)
            XQRAL2(N,1,N5B,N4B)=XQRAL2(N,1,N5B,N4B)+RQRAL2(N,1,N5B,N4B)
            XQRAL3(N,1,N5B,N4B)=XQRAL3(N,1,N5B,N4B)+RQRAL3(N,1,N5B,N4B)
            XQRAL4(N,1,N5B,N4B)=XQRAL4(N,1,N5B,N4B)+RQRAL4(N,1,N5B,N4B)
            XQRALS(N,1,N5B,N4B)=XQRALS(N,1,N5B,N4B)+RQRALS(N,1,N5B,N4B)
            XQRFE1(N,1,N5B,N4B)=XQRFE1(N,1,N5B,N4B)+RQRFE1(N,1,N5B,N4B)
            XQRFE2(N,1,N5B,N4B)=XQRFE2(N,1,N5B,N4B)+RQRFE2(N,1,N5B,N4B)
            XQRFE3(N,1,N5B,N4B)=XQRFE3(N,1,N5B,N4B)+RQRFE3(N,1,N5B,N4B)
            XQRFE4(N,1,N5B,N4B)=XQRFE4(N,1,N5B,N4B)+RQRFE4(N,1,N5B,N4B)
            XQRFES(N,1,N5B,N4B)=XQRFES(N,1,N5B,N4B)+RQRFES(N,1,N5B,N4B)
            XQRCAO(N,1,N5B,N4B)=XQRCAO(N,1,N5B,N4B)+RQRCAO(N,1,N5B,N4B)
            XQRCAC(N,1,N5B,N4B)=XQRCAC(N,1,N5B,N4B)+RQRCAC(N,1,N5B,N4B)
            XQRCAH(N,1,N5B,N4B)=XQRCAH(N,1,N5B,N4B)+RQRCAH(N,1,N5B,N4B)
            XQRCAS(N,1,N5B,N4B)=XQRCAS(N,1,N5B,N4B)+RQRCAS(N,1,N5B,N4B)
            XQRMGO(N,1,N5B,N4B)=XQRMGO(N,1,N5B,N4B)+RQRMGO(N,1,N5B,N4B)
            XQRMGC(N,1,N5B,N4B)=XQRMGC(N,1,N5B,N4B)+RQRMGC(N,1,N5B,N4B)
            XQRMGH(N,1,N5B,N4B)=XQRMGH(N,1,N5B,N4B)+RQRMGH(N,1,N5B,N4B)
            XQRMGS(N,1,N5B,N4B)=XQRMGS(N,1,N5B,N4B)+RQRMGS(N,1,N5B,N4B)
            XQRNAC(N,1,N5B,N4B)=XQRNAC(N,1,N5B,N4B)+RQRNAC(N,1,N5B,N4B)
            XQRNAS(N,1,N5B,N4B)=XQRNAS(N,1,N5B,N4B)+RQRNAS(N,1,N5B,N4B)
            XQRKAS(N,1,N5B,N4B)=XQRKAS(N,1,N5B,N4B)+RQRKAS(N,1,N5B,N4B)
            XQRH0P(N,1,N5B,N4B)=XQRH0P(N,1,N5B,N4B)+RQRH0P(N,1,N5B,N4B)
            XQRH3P(N,1,N5B,N4B)=XQRH3P(N,1,N5B,N4B)+RQRH3P(N,1,N5B,N4B)
            XQRF1P(N,1,N5B,N4B)=XQRF1P(N,1,N5B,N4B)+RQRF1P(N,1,N5B,N4B)
            XQRF2P(N,1,N5B,N4B)=XQRF2P(N,1,N5B,N4B)+RQRF2P(N,1,N5B,N4B)
            XQRC0P(N,1,N5B,N4B)=XQRC0P(N,1,N5B,N4B)+RQRC0P(N,1,N5B,N4B)
            XQRC1P(N,1,N5B,N4B)=XQRC1P(N,1,N5B,N4B)+RQRC1P(N,1,N5B,N4B)
            XQRC2P(N,1,N5B,N4B)=XQRC2P(N,1,N5B,N4B)+RQRC2P(N,1,N5B,N4B)
            XQRM1P(N,1,N5B,N4B)=XQRM1P(N,1,N5B,N4B)+RQRM1P(N,1,N5B,N4B)
          ELSE
            RQRAL(N,1,N5B,N4B)=0.0_r8
            RQRFE(N,1,N5B,N4B)=0.0_r8
            RQRHY(N,1,N5B,N4B)=0.0_r8
            RQRCA(N,1,N5B,N4B)=0.0_r8
            RQRMG(N,1,N5B,N4B)=0.0_r8
            RQRNA(N,1,N5B,N4B)=0.0_r8
            RQRKA(N,1,N5B,N4B)=0.0_r8
            RQROH(N,1,N5B,N4B)=0.0_r8
            RQRSO(N,1,N5B,N4B)=0.0_r8
            RQRCL(N,1,N5B,N4B)=0.0_r8
            RQRC3(N,1,N5B,N4B)=0.0_r8
            RQRHC(N,1,N5B,N4B)=0.0_r8
            RQRAL1(N,1,N5B,N4B)=0.0_r8
            RQRAL2(N,1,N5B,N4B)=0.0_r8
            RQRAL3(N,1,N5B,N4B)=0.0_r8
            RQRAL4(N,1,N5B,N4B)=0.0_r8
            RQRALS(N,1,N5B,N4B)=0.0_r8
            RQRFE1(N,1,N5B,N4B)=0.0_r8
            RQRFE2(N,1,N5B,N4B)=0.0_r8
            RQRFE3(N,1,N5B,N4B)=0.0_r8
            RQRFE4(N,1,N5B,N4B)=0.0_r8
            RQRFES(N,1,N5B,N4B)=0.0_r8
            RQRCAO(N,1,N5B,N4B)=0.0_r8
            RQRCAC(N,1,N5B,N4B)=0.0_r8
            RQRCAH(N,1,N5B,N4B)=0.0_r8
            RQRCAS(N,1,N5B,N4B)=0.0_r8
            RQRMGO(N,1,N5B,N4B)=0.0_r8
            RQRMGC(N,1,N5B,N4B)=0.0_r8
            RQRMGH(N,1,N5B,N4B)=0.0_r8
            RQRMGS(N,1,N5B,N4B)=0.0_r8
            RQRNAC(N,1,N5B,N4B)=0.0_r8
            RQRNAS(N,1,N5B,N4B)=0.0_r8
            RQRKAS(N,1,N5B,N4B)=0.0_r8
            RQRH0P(N,1,N5B,N4B)=0.0_r8
            RQRH3P(N,1,N5B,N4B)=0.0_r8
            RQRF1P(N,1,N5B,N4B)=0.0_r8
            RQRF2P(N,1,N5B,N4B)=0.0_r8
            RQRC0P(N,1,N5B,N4B)=0.0_r8
            RQRC1P(N,1,N5B,N4B)=0.0_r8
            RQRC2P(N,1,N5B,N4B)=0.0_r8
            RQRM1P(N,1,N5B,N4B)=0.0_r8
          ENDIF
        ENDIF
      ELSE
        RQRAL(N,2,N5,N4)=0.0_r8
        RQRFE(N,2,N5,N4)=0.0_r8
        RQRHY(N,2,N5,N4)=0.0_r8
        RQRCA(N,2,N5,N4)=0.0_r8
        RQRMG(N,2,N5,N4)=0.0_r8
        RQRNA(N,2,N5,N4)=0.0_r8
        RQRKA(N,2,N5,N4)=0.0_r8
        RQROH(N,2,N5,N4)=0.0_r8
        RQRSO(N,2,N5,N4)=0.0_r8
        RQRCL(N,2,N5,N4)=0.0_r8
        RQRC3(N,2,N5,N4)=0.0_r8
        RQRHC(N,2,N5,N4)=0.0_r8
        RQRAL1(N,2,N5,N4)=0.0_r8
        RQRAL2(N,2,N5,N4)=0.0_r8
        RQRAL3(N,2,N5,N4)=0.0_r8
        RQRAL4(N,2,N5,N4)=0.0_r8
        RQRALS(N,2,N5,N4)=0.0_r8
        RQRFE1(N,2,N5,N4)=0.0_r8
        RQRFE2(N,2,N5,N4)=0.0_r8
        RQRFE3(N,2,N5,N4)=0.0_r8
        RQRFE4(N,2,N5,N4)=0.0_r8
        RQRFES(N,2,N5,N4)=0.0_r8
        RQRCAO(N,2,N5,N4)=0.0_r8
        RQRCAC(N,2,N5,N4)=0.0_r8
        RQRCAH(N,2,N5,N4)=0.0_r8
        RQRCAS(N,2,N5,N4)=0.0_r8
        RQRMGO(N,2,N5,N4)=0.0_r8
        RQRMGC(N,2,N5,N4)=0.0_r8
        RQRMGH(N,2,N5,N4)=0.0_r8
        RQRMGS(N,2,N5,N4)=0.0_r8
        RQRNAC(N,2,N5,N4)=0.0_r8
        RQRNAS(N,2,N5,N4)=0.0_r8
        RQRKAS(N,2,N5,N4)=0.0_r8
        RQRH0P(N,2,N5,N4)=0.0_r8
        RQRH3P(N,2,N5,N4)=0.0_r8
        RQRF1P(N,2,N5,N4)=0.0_r8
        RQRF2P(N,2,N5,N4)=0.0_r8
        RQRC0P(N,2,N5,N4)=0.0_r8
        RQRC1P(N,2,N5,N4)=0.0_r8
        RQRC2P(N,2,N5,N4)=0.0_r8
        RQRM1P(N,2,N5,N4)=0.0_r8
        IF(N4B.GT.0.AND.N5B.GT.0)THEN
          RQRAL(N,1,N5B,N4B)=0.0_r8
          RQRFE(N,1,N5B,N4B)=0.0_r8
          RQRHY(N,1,N5B,N4B)=0.0_r8
          RQRCA(N,1,N5B,N4B)=0.0_r8
          RQRMG(N,1,N5B,N4B)=0.0_r8
          RQRNA(N,1,N5B,N4B)=0.0_r8
          RQRKA(N,1,N5B,N4B)=0.0_r8
          RQROH(N,1,N5B,N4B)=0.0_r8
          RQRSO(N,1,N5B,N4B)=0.0_r8
          RQRCL(N,1,N5B,N4B)=0.0_r8
          RQRC3(N,1,N5B,N4B)=0.0_r8
          RQRHC(N,1,N5B,N4B)=0.0_r8
          RQRAL1(N,1,N5B,N4B)=0.0_r8
          RQRAL2(N,1,N5B,N4B)=0.0_r8
          RQRAL3(N,1,N5B,N4B)=0.0_r8
          RQRAL4(N,1,N5B,N4B)=0.0_r8
          RQRALS(N,1,N5B,N4B)=0.0_r8
          RQRFE1(N,1,N5B,N4B)=0.0_r8
          RQRFE2(N,1,N5B,N4B)=0.0_r8
          RQRFE3(N,1,N5B,N4B)=0.0_r8
          RQRFE4(N,1,N5B,N4B)=0.0_r8
          RQRFES(N,1,N5B,N4B)=0.0_r8
          RQRCAO(N,1,N5B,N4B)=0.0_r8
          RQRCAC(N,1,N5B,N4B)=0.0_r8
          RQRCAH(N,1,N5B,N4B)=0.0_r8
          RQRCAS(N,1,N5B,N4B)=0.0_r8
          RQRMGO(N,1,N5B,N4B)=0.0_r8
          RQRMGC(N,1,N5B,N4B)=0.0_r8
          RQRMGH(N,1,N5B,N4B)=0.0_r8
          RQRMGS(N,1,N5B,N4B)=0.0_r8
          RQRNAC(N,1,N5B,N4B)=0.0_r8
          RQRNAS(N,1,N5B,N4B)=0.0_r8
          RQRKAS(N,1,N5B,N4B)=0.0_r8
          RQRH0P(N,1,N5B,N4B)=0.0_r8
          RQRH3P(N,1,N5B,N4B)=0.0_r8
          RQRF1P(N,1,N5B,N4B)=0.0_r8
          RQRF2P(N,1,N5B,N4B)=0.0_r8
          RQRC0P(N,1,N5B,N4B)=0.0_r8
          RQRC1P(N,1,N5B,N4B)=0.0_r8
          RQRC2P(N,1,N5B,N4B)=0.0_r8
          RQRM1P(N,1,N5B,N4B)=0.0_r8
        ENDIF
      ENDIF
!
!     SNOW DRIFT
!
!     QSM=snow transfer from watsub.f
!     RQS*=solute flux in snow transfer
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     VOLS=volume of snowpack from watsub.f
!     *W2=solute content of snowpack
!
      IF(NN.EQ.1)THEN
!
!     IF NO SNOW DRIFT THEN NO TRANSPORT
!
        IF(ABS(QSM(M,N,N5,N4)).LE.ZEROS2(N2,N1))THEN
          RQSAL(N,N5,N4)=0.0_r8
          RQSFE(N,N5,N4)=0.0_r8
          RQSHY(N,N5,N4)=0.0_r8
          RQSCA(N,N5,N4)=0.0_r8
          RQSMG(N,N5,N4)=0.0_r8
          RQSNA(N,N5,N4)=0.0_r8
          RQSKA(N,N5,N4)=0.0_r8
          RQSOH(N,N5,N4)=0.0_r8
          RQSSO(N,N5,N4)=0.0_r8
          RQSCL(N,N5,N4)=0.0_r8
          RQSC3(N,N5,N4)=0.0_r8
          RQSHC(N,N5,N4)=0.0_r8
          RQSAL1(N,N5,N4)=0.0_r8
          RQSAL2(N,N5,N4)=0.0_r8
          RQSAL3(N,N5,N4)=0.0_r8
          RQSAL4(N,N5,N4)=0.0_r8
          RQSALS(N,N5,N4)=0.0_r8
          RQSFE1(N,N5,N4)=0.0_r8
          RQSFE2(N,N5,N4)=0.0_r8
          RQSFE3(N,N5,N4)=0.0_r8
          RQSFE4(N,N5,N4)=0.0_r8
          RQSFES(N,N5,N4)=0.0_r8
          RQSCAO(N,N5,N4)=0.0_r8
          RQSCAC(N,N5,N4)=0.0_r8
          RQSCAH(N,N5,N4)=0.0_r8
          RQSCAS(N,N5,N4)=0.0_r8
          RQSMGO(N,N5,N4)=0.0_r8
          RQSMGC(N,N5,N4)=0.0_r8
          RQSMGH(N,N5,N4)=0.0_r8
          RQSMGS(N,N5,N4)=0.0_r8
          RQSNAC(N,N5,N4)=0.0_r8
          RQSNAS(N,N5,N4)=0.0_r8
          RQSKAS(N,N5,N4)=0.0_r8
          RQSH0P(N,N5,N4)=0.0_r8
          RQSH3P(N,N5,N4)=0.0_r8
          RQSF1P(N,N5,N4)=0.0_r8
          RQSF2P(N,N5,N4)=0.0_r8
          RQSC0P(N,N5,N4)=0.0_r8
          RQSC1P(N,N5,N4)=0.0_r8
          RQSC2P(N,N5,N4)=0.0_r8
          RQSM1P(N,N5,N4)=0.0_r8
!
!     IF DRIFT IS FROM CURRENT TO ADJACENT GRID CELL
!
        ELSEIF(QSM(M,N,N5,N4).GT.0.0)THEN
          IF(VOLS(N2,N1).GT.ZEROS2(NY,NX))THEN
            VFLW=AZMAX1(AMIN1(VFLWX,QSM(M,N,N5,N4)/VOLS(N2,N1)))
          ELSE
            VFLW=VFLWX
          ENDIF
          RQSAL(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Al,1,N2,N1))
          RQSFE(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Fe,1,N2,N1))
          RQSHY(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Hp,1,N2,N1))
          RQSCA(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Ca,1,N2,N1))
          RQSMG(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Mg,1,N2,N1))
          RQSNA(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Na,1,N2,N1))
          RQSKA(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_K,1,N2,N1))
          RQSOH(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_OH,1,N2,N1))
          RQSSO(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_SO4,1,N2,N1))
          RQSCL(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Cl,1,N2,N1))
          RQSC3(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CO3,1,N2,N1))
          RQSHC(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_HCO3,1,N2,N1))
          RQSAL1(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_AlOH,1,N2,N1))
          RQSAL2(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_AlOH2,1,N2,N1))
          RQSAL3(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_AlOH3,1,N2,N1))
          RQSAL4(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_AlOH4,1,N2,N1))
          RQSALS(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_AlSO4,1,N2,N1))
          RQSFE1(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeOH,1,N2,N1))
          RQSFE2(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeOH2,1,N2,N1))
          RQSFE3(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeOH3,1,N2,N1))
          RQSFE4(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeOH4,1,N2,N1))
          RQSFES(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeSO4,1,N2,N1))
          RQSCAO(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaOH2,1,N2,N1))
          RQSCAC(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaCO3,1,N2,N1))
          RQSCAH(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaHCO3,1,N2,N1))
          RQSCAS(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaSO4,1,N2,N1))
          RQSMGO(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_MgOH2,1,N2,N1))
          RQSMGC(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_MgCO3,1,N2,N1))
          RQSMGH(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_MgHCO3,1,N2,N1))
          RQSMGS(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_MgSO4,1,N2,N1))
          RQSNAC(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_NaCO3,1,N2,N1))
          RQSNAS(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_NaSO4,1,N2,N1))
          RQSKAS(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_KSO4,1,N2,N1))
          RQSH0P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_H0PO4,1,N2,N1))
          RQSH3P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_H3PO4,1,N2,N1))
          RQSF1P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeHPO4,1,N2,N1))
          RQSF2P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeH2PO4,1,N2,N1))
          RQSC0P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaPO4,1,N2,N1))
          RQSC1P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaHPO4,1,N2,N1))
          RQSC2P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaH2PO4,1,N2,N1))
          RQSM1P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_MgHPO4,1,N2,N1))
!
!     IF DRIFT IS TO CURRENT FROM ADJACENT GRID CELL
!
        ELSEIF(QSM(M,N,N5,N4).LT.-ZEROS2(N2,N1))THEN
          IF(VOLS(N5,N4).GT.ZEROS2(N5,N4))THEN
            VFLW=AZMIN1(AMAX1(-VFLWX,QSM(M,N,N5,N4)/VOLS(N5,N4)))
          ELSE
            VFLW=-VFLWX
          ENDIF
          RQSAL(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Al,1,N5,N4))
          RQSFE(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Fe,1,N5,N4))
          RQSHY(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Hp,1,N5,N4))
          RQSCA(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Ca,1,N5,N4))
          RQSMG(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Mg,1,N5,N4))
          RQSNA(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Na,1,N5,N4))
          RQSKA(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_K,1,N5,N4))
          RQSOH(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_OH,1,N5,N4))
          RQSSO(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_SO4,1,N5,N4))
          RQSCL(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_Cl,1,N5,N4))
          RQSC3(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CO3,1,N5,N4))
          RQSHC(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_HCO3,1,N5,N4))
          RQSAL1(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_AlOH,1,N5,N4))
          RQSAL2(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_AlOH2,1,N5,N4))
          RQSAL3(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_AlOH3,1,N5,N4))
          RQSAL4(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_AlOH4,1,N5,N4))
          RQSALS(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_AlSO4,1,N5,N4))
          RQSFE1(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeOH,1,N5,N4))
          RQSFE2(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeOH2,1,N5,N4))
          RQSFE3(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeOH3,1,N5,N4))
          RQSFE4(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeOH4,1,N5,N4))
          RQSFES(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeSO4,1,N5,N4))
          RQSCAO(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaOH2,1,N5,N4))
          RQSCAC(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaCO3,1,N5,N4))
          RQSCAH(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaHCO3,1,N5,N4))
          RQSCAS(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaSO4,1,N5,N4))
          RQSMGO(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_MgOH2,1,N5,N4))
          RQSMGC(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_MgCO3,1,N5,N4))
          RQSMGH(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_MgHCO3,1,N5,N4))
          RQSMGS(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_MgSO4,1,N5,N4))
          RQSNAC(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_NaCO3,1,N5,N4))
          RQSNAS(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_NaSO4,1,N5,N4))
          RQSKAS(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_KSO4,1,N5,N4))
          RQSH0P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_H0PO4,1,N5,N4))
          RQSH3P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_H3PO4,1,N5,N4))
          RQSF1P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeHPO4,1,N5,N4))
          RQSF2P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_FeH2PO4,1,N5,N4))
          RQSC0P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaPO4,1,N5,N4))
          RQSC1P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaHPO4,1,N5,N4))
          RQSC2P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_CaH2PO4,1,N5,N4))
          RQSM1P(N,N5,N4)=VFLW*AZMAX1(trcs_solsml2(idsa_MgHPO4,1,N5,N4))
        ENDIF
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     XQS*=hourly solute in snow transfer
!     RQS*=solute in snow transfer
!
        XQSAL(N,N5,N4)=XQSAL(N,N5,N4)+RQSAL(N,N5,N4)
        XQSFE(N,N5,N4)=XQSFE(N,N5,N4)+RQSFE(N,N5,N4)
        XQSHY(N,N5,N4)=XQSHY(N,N5,N4)+RQSHY(N,N5,N4)
        XQSCA(N,N5,N4)=XQSCA(N,N5,N4)+RQSCA(N,N5,N4)
        XQSMG(N,N5,N4)=XQSMG(N,N5,N4)+RQSMG(N,N5,N4)
        XQSNA(N,N5,N4)=XQSNA(N,N5,N4)+RQSNA(N,N5,N4)
        XQSKA(N,N5,N4)=XQSKA(N,N5,N4)+RQSKA(N,N5,N4)
        XQSOH(N,N5,N4)=XQSOH(N,N5,N4)+RQSOH(N,N5,N4)
        XQSSO(N,N5,N4)=XQSSO(N,N5,N4)+RQSSO(N,N5,N4)
        XQSCL(N,N5,N4)=XQSCL(N,N5,N4)+RQSCL(N,N5,N4)
        XQSC3(N,N5,N4)=XQSC3(N,N5,N4)+RQSC3(N,N5,N4)
        XQSHC(N,N5,N4)=XQSHC(N,N5,N4)+RQSHC(N,N5,N4)
        XQSAL1(N,N5,N4)=XQSAL1(N,N5,N4)+RQSAL1(N,N5,N4)
        XQSAL2(N,N5,N4)=XQSAL2(N,N5,N4)+RQSAL2(N,N5,N4)
        XQSAL3(N,N5,N4)=XQSAL3(N,N5,N4)+RQSAL3(N,N5,N4)
        XQSAL4(N,N5,N4)=XQSAL4(N,N5,N4)+RQSAL4(N,N5,N4)
        XQSALS(N,N5,N4)=XQSALS(N,N5,N4)+RQSALS(N,N5,N4)
        XQSFE1(N,N5,N4)=XQSFE1(N,N5,N4)+RQSFE1(N,N5,N4)
        XQSFE2(N,N5,N4)=XQSFE2(N,N5,N4)+RQSFE2(N,N5,N4)
        XQSFE3(N,N5,N4)=XQSFE3(N,N5,N4)+RQSFE3(N,N5,N4)
        XQSFE4(N,N5,N4)=XQSFE4(N,N5,N4)+RQSFE4(N,N5,N4)
        XQSFES(N,N5,N4)=XQSFES(N,N5,N4)+RQSFES(N,N5,N4)
        XQSCAO(N,N5,N4)=XQSCAO(N,N5,N4)+RQSCAO(N,N5,N4)
        XQSCAC(N,N5,N4)=XQSCAC(N,N5,N4)+RQSCAC(N,N5,N4)
        XQSCAH(N,N5,N4)=XQSCAH(N,N5,N4)+RQSCAH(N,N5,N4)
        XQSCAS(N,N5,N4)=XQSCAS(N,N5,N4)+RQSCAS(N,N5,N4)
        XQSMGO(N,N5,N4)=XQSMGO(N,N5,N4)+RQSMGO(N,N5,N4)
        XQSMGC(N,N5,N4)=XQSMGC(N,N5,N4)+RQSMGC(N,N5,N4)
        XQSMGH(N,N5,N4)=XQSMGH(N,N5,N4)+RQSMGH(N,N5,N4)
        XQSMGS(N,N5,N4)=XQSMGS(N,N5,N4)+RQSMGS(N,N5,N4)
        XQSNAC(N,N5,N4)=XQSNAC(N,N5,N4)+RQSNAC(N,N5,N4)
        XQSNAS(N,N5,N4)=XQSNAS(N,N5,N4)+RQSNAS(N,N5,N4)
        XQSKAS(N,N5,N4)=XQSKAS(N,N5,N4)+RQSKAS(N,N5,N4)
        XQSH0P(N,N5,N4)=XQSH0P(N,N5,N4)+RQSH0P(N,N5,N4)
        XQSH3P(N,N5,N4)=XQSH3P(N,N5,N4)+RQSH3P(N,N5,N4)
        XQSF1P(N,N5,N4)=XQSF1P(N,N5,N4)+RQSF1P(N,N5,N4)
        XQSF2P(N,N5,N4)=XQSF2P(N,N5,N4)+RQSF2P(N,N5,N4)
        XQSC0P(N,N5,N4)=XQSC0P(N,N5,N4)+RQSC0P(N,N5,N4)
        XQSC1P(N,N5,N4)=XQSC1P(N,N5,N4)+RQSC1P(N,N5,N4)
        XQSC2P(N,N5,N4)=XQSC2P(N,N5,N4)+RQSC2P(N,N5,N4)
        XQSM1P(N,N5,N4)=XQSM1P(N,N5,N4)+RQSM1P(N,N5,N4)
      ENDIF
    enddo
  enddo
  end subroutine UpdateSoluteInSurfNeighbors
!------------------------------------------------------------------------------------------
   subroutine SoluteAdvMicropore(M,N,N1,N2,N3,N4,N5,N6)
   implicit none
   integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
   real(r8) :: VFLW
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
!     FLWM=water flux through soil micropore from watsub.f
!     VOLWM=micropore water-filled porosity from watsub.f
!     RFL*S=solute diffusive flux through micropore
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *S2,*B2=micropore solute content in non-band,band
!
  IF(FLWM(M,N,N6,N5,N4).GT.0.0)THEN
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
    IF(VOLWM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FLWM(M,N,N6,N5,N4)/VOLWM(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    RFLAL=VFLW*AZMAX1(ZAL2(N3,N2,N1))
    RFLFE=VFLW*AZMAX1(ZFE2(N3,N2,N1))
    RFLHY=VFLW*AZMAX1(ZHY2(N3,N2,N1))
    RFLCA=VFLW*AZMAX1(ZCA2(N3,N2,N1))
    RFLMG=VFLW*AZMAX1(ZMG2(N3,N2,N1))
    RFLNA=VFLW*AZMAX1(ZNA2(N3,N2,N1))
    RFLKA=VFLW*AZMAX1(ZKA2(N3,N2,N1))
    RFLOH=VFLW*AZMAX1(ZOH2(N3,N2,N1))
    RFLSO=VFLW*AZMAX1(ZSO42(N3,N2,N1))
    RFLCL=VFLW*AZMAX1(ZCL2(N3,N2,N1))
    RFLC3=VFLW*AZMAX1(ZCO32(N3,N2,N1))
    RFLHC=VFLW*AZMAX1(ZHCO32(N3,N2,N1))
    RFLAL1=VFLW*AZMAX1(ZAL12(N3,N2,N1))
    RFLAL2=VFLW*AZMAX1(ZAL22(N3,N2,N1))
    RFLAL3=VFLW*AZMAX1(ZAL32(N3,N2,N1))
    RFLAL4=VFLW*AZMAX1(ZAL42(N3,N2,N1))
    RFLALS=VFLW*AZMAX1(ZALS2(N3,N2,N1))
    RFLFE1=VFLW*AZMAX1(ZFE12(N3,N2,N1))
    RFLFE2=VFLW*AZMAX1(ZFE22(N3,N2,N1))
    RFLFE3=VFLW*AZMAX1(ZFE32(N3,N2,N1))
    RFLFE4=VFLW*AZMAX1(ZFE42(N3,N2,N1))
    RFLFES=VFLW*AZMAX1(ZFES2(N3,N2,N1))
    RFLCAO=VFLW*AZMAX1(ZCAO2(N3,N2,N1))
    RFLCAC=VFLW*AZMAX1(ZCAC2(N3,N2,N1))
    RFLCAH=VFLW*AZMAX1(ZCAH2(N3,N2,N1))
    RFLCAS=VFLW*AZMAX1(ZCAS2(N3,N2,N1))
    RFLMGO=VFLW*AZMAX1(ZMGO2(N3,N2,N1))
    RFLMGC=VFLW*AZMAX1(ZMGC2(N3,N2,N1))
    RFLMGH=VFLW*AZMAX1(ZMGH2(N3,N2,N1))
    RFLMGS=VFLW*AZMAX1(ZMGS2(N3,N2,N1))
    RFLNAC=VFLW*AZMAX1(ZNAC2(N3,N2,N1))
    RFLNAS=VFLW*AZMAX1(ZNAS2(N3,N2,N1))
    RFLKAS=VFLW*AZMAX1(ZKAS2(N3,N2,N1))
    RFLH0P=VFLW*AZMAX1(H0PO42(N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLH3P=VFLW*AZMAX1(H3PO42(N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLF1P=VFLW*AZMAX1(ZFE1P2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLF2P=VFLW*AZMAX1(ZFE2P2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLC0P=VFLW*AZMAX1(ZCA0P2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLC1P=VFLW*AZMAX1(ZCA1P2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLC2P=VFLW*AZMAX1(ZCA2P2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLM1P=VFLW*AZMAX1(ZMG1P2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N3,N2,N1)
    RFLH0B=VFLW*AZMAX1(H0POB2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLH3B=VFLW*AZMAX1(H3POB2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLF1B=VFLW*AZMAX1(ZF1PB2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLF2B=VFLW*AZMAX1(ZF2PB2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLC0B=VFLW*AZMAX1(ZC0PB2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLC1B=VFLW*AZMAX1(ZC1PB2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLC2B=VFLW*AZMAX1(ZC2PB2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    RFLM1B=VFLW*AZMAX1(ZM1PB2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
!
!     IF MICROPORE WATER FLUX FROM 'WATSUB' IS TO CURRENT FROM
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MICROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN ADJACENT GRID CELL
!
  ELSE
    IF(VOLWM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FLWM(M,N,N6,N5,N4)/VOLWM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    RFLAL=VFLW*AZMAX1(ZAL2(N6,N5,N4))
    RFLFE=VFLW*AZMAX1(ZFE2(N6,N5,N4))
    RFLHY=VFLW*AZMAX1(ZHY2(N6,N5,N4))
    RFLCA=VFLW*AZMAX1(ZCA2(N6,N5,N4))
    RFLMG=VFLW*AZMAX1(ZMG2(N6,N5,N4))
    RFLNA=VFLW*AZMAX1(ZNA2(N6,N5,N4))
    RFLKA=VFLW*AZMAX1(ZKA2(N6,N5,N4))
    RFLOH=VFLW*AZMAX1(ZOH2(N6,N5,N4))
    RFLSO=VFLW*AZMAX1(ZSO42(N6,N5,N4))
    RFLCL=VFLW*AZMAX1(ZCL2(N6,N5,N4))
    RFLC3=VFLW*AZMAX1(ZCO32(N6,N5,N4))
    RFLHC=VFLW*AZMAX1(ZHCO32(N6,N5,N4))
    RFLAL1=VFLW*AZMAX1(ZAL12(N6,N5,N4))
    RFLAL2=VFLW*AZMAX1(ZAL22(N6,N5,N4))
    RFLAL3=VFLW*AZMAX1(ZAL32(N6,N5,N4))
    RFLAL4=VFLW*AZMAX1(ZAL42(N6,N5,N4))
    RFLALS=VFLW*AZMAX1(ZALS2(N6,N5,N4))
    RFLFE1=VFLW*AZMAX1(ZFE12(N6,N5,N4))
    RFLFE2=VFLW*AZMAX1(ZFE22(N6,N5,N4))
    RFLFE3=VFLW*AZMAX1(ZFE32(N6,N5,N4))
    RFLFE4=VFLW*AZMAX1(ZFE42(N6,N5,N4))
    RFLFES=VFLW*AZMAX1(ZFES2(N6,N5,N4))
    RFLCAO=VFLW*AZMAX1(ZCAO2(N6,N5,N4))
    RFLCAC=VFLW*AZMAX1(ZCAC2(N6,N5,N4))
    RFLCAH=VFLW*AZMAX1(ZCAH2(N6,N5,N4))
    RFLCAS=VFLW*AZMAX1(ZCAS2(N6,N5,N4))
    RFLMGO=VFLW*AZMAX1(ZMGO2(N6,N5,N4))
    RFLMGC=VFLW*AZMAX1(ZMGC2(N6,N5,N4))
    RFLMGH=VFLW*AZMAX1(ZMGH2(N6,N5,N4))
    RFLMGS=VFLW*AZMAX1(ZMGS2(N6,N5,N4))
    RFLNAC=VFLW*AZMAX1(ZNAC2(N6,N5,N4))
    RFLNAS=VFLW*AZMAX1(ZNAS2(N6,N5,N4))
    RFLKAS=VFLW*AZMAX1(ZKAS2(N6,N5,N4))
    RFLH0P=VFLW*AZMAX1(H0PO42(N6,N5,N4))
    RFLH3P=VFLW*AZMAX1(H3PO42(N6,N5,N4))
    RFLF1P=VFLW*AZMAX1(ZFE1P2(N6,N5,N4))
    RFLF2P=VFLW*AZMAX1(ZFE2P2(N6,N5,N4))
    RFLC0P=VFLW*AZMAX1(ZCA0P2(N6,N5,N4))
    RFLC1P=VFLW*AZMAX1(ZCA1P2(N6,N5,N4))
    RFLC2P=VFLW*AZMAX1(ZCA2P2(N6,N5,N4))
    RFLM1P=VFLW*AZMAX1(ZMG1P2(N6,N5,N4))
    RFLH0B=VFLW*AZMAX1(H0POB2(N6,N5,N4))
    RFLH3B=VFLW*AZMAX1(H3POB2(N6,N5,N4))
    RFLF1B=VFLW*AZMAX1(ZF1PB2(N6,N5,N4))
    RFLF2B=VFLW*AZMAX1(ZF2PB2(N6,N5,N4))
    RFLC0B=VFLW*AZMAX1(ZC0PB2(N6,N5,N4))
    RFLC1B=VFLW*AZMAX1(ZC1PB2(N6,N5,N4))
    RFLC2B=VFLW*AZMAX1(ZC2PB2(N6,N5,N4))
    RFLM1B=VFLW*AZMAX1(ZM1PB2(N6,N5,N4))
  ENDIF

!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores

  RALFLS(N,N6,N5,N4)=RFLAL
  RFEFLS(N,N6,N5,N4)=RFLFE
  RHYFLS(N,N6,N5,N4)=RFLHY
  RCAFLS(N,N6,N5,N4)=RFLCA
  RMGFLS(N,N6,N5,N4)=RFLMG
  RNAFLS(N,N6,N5,N4)=RFLNA
  RKAFLS(N,N6,N5,N4)=RFLKA
  ROHFLS(N,N6,N5,N4)=RFLOH
  RSOFLS(N,N6,N5,N4)=RFLSO
  RCLFLS(N,N6,N5,N4)=RFLCL
  RC3FLS(N,N6,N5,N4)=RFLC3
  RHCFLS(N,N6,N5,N4)=RFLHC
  RAL1FS(N,N6,N5,N4)=RFLAL1
  RAL2FS(N,N6,N5,N4)=RFLAL2
  RAL3FS(N,N6,N5,N4)=RFLAL3
  RAL4FS(N,N6,N5,N4)=RFLAL4
  RALSFS(N,N6,N5,N4)=RFLALS
  RFE1FS(N,N6,N5,N4)=RFLFE1
  RFE2FS(N,N6,N5,N4)=RFLFE2
  RFE3FS(N,N6,N5,N4)=RFLFE3
  RFE4FS(N,N6,N5,N4)=RFLFE4
  RFESFS(N,N6,N5,N4)=RFLFES
  RCAOFS(N,N6,N5,N4)=RFLCAO
  RCACFS(N,N6,N5,N4)=RFLCAC
  RCAHFS(N,N6,N5,N4)=RFLCAH
  RCASFS(N,N6,N5,N4)=RFLCAS
  RMGOFS(N,N6,N5,N4)=RFLMGO
  RMGCFS(N,N6,N5,N4)=RFLMGC
  RMGHFS(N,N6,N5,N4)=RFLMGH
  RMGSFS(N,N6,N5,N4)=RFLMGS
  RNACFS(N,N6,N5,N4)=RFLNAC
  RNASFS(N,N6,N5,N4)=RFLNAS
  RKASFS(N,N6,N5,N4)=RFLKAS
  RH0PFS(N,N6,N5,N4)=RFLH0P
  RH3PFS(N,N6,N5,N4)=RFLH3P
  RF1PFS(N,N6,N5,N4)=RFLF1P
  RF2PFS(N,N6,N5,N4)=RFLF2P
  RC0PFS(N,N6,N5,N4)=RFLC0P
  RC1PFS(N,N6,N5,N4)=RFLC1P
  RC2PFS(N,N6,N5,N4)=RFLC2P
  RM1PFS(N,N6,N5,N4)=RFLM1P
  RH0BFB(N,N6,N5,N4)=RFLH0B
  RH3BFB(N,N6,N5,N4)=RFLH3B
  RF1BFB(N,N6,N5,N4)=RFLF1B
  RF2BFB(N,N6,N5,N4)=RFLF2B
  RC0BFB(N,N6,N5,N4)=RFLC0B
  RC1BFB(N,N6,N5,N4)=RFLC1B
  RC2BFB(N,N6,N5,N4)=RFLC2B
  RM1BFB(N,N6,N5,N4)=RFLM1B
  end subroutine SoluteAdvMicropore
!------------------------------------------------------------------------------------------
  subroutine SoluteDifsMicropore(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)
  REAL(R8) :: VLWPA1,VLWPB1,VLWPA2,VLWPB2
!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MICROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
  IF(THETW1(N3,N2,N1).GT.THETY(N3,N2,N1) &
    .AND.THETW1(N6,N5,N4).GT.THETY(N6,N5,N4) &
    .AND.VOLWM(M,N3,N2,N1).GT.ZEROS2(N2,N1) &
    .AND.VOLWM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
!
!     VOLWM=micropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!
!     MICROPORE CONCENTRATIONS FROM WATER-FILLED POROSITY
!     IN CURRENT AND ADJACENT GRID CELLS
!
!     C*1,C*2=solute concentration in source,destination layer
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *S2,*B2=soil solute content in non-band,band
!
    VLWPA1=VOLWM(M,N3,N2,N1)*trcs_VLN(ids_H1PO4,N3,N2,N1)
    VLWPB1=VOLWM(M,N3,N2,N1)*trcs_VLN(ids_H1PO4B,N3,N2,N1)
    VLWPA2=VOLWM(M,N6,N5,N4)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    VLWPB2=VOLWM(M,N6,N5,N4)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    CAL1=AZMAX1(ZAL2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFE1=AZMAX1(ZFE2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CHY1=AZMAX1(ZHY2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCA1=AZMAX1(ZCA2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CMG1=AZMAX1(ZMG2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CNA1=AZMAX1(ZNA2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CKA1=AZMAX1(ZKA2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    COH1=AZMAX1(ZOH2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CSO41=AZMAX1(ZSO42(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCL1=AZMAX1(ZCL2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCO31=AZMAX1(ZCO32(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CHCO31=AZMAX1(ZHCO32(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CAL11=AZMAX1(ZAL12(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CAL21=AZMAX1(ZAL22(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CAL31=AZMAX1(ZAL32(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CAL41=AZMAX1(ZAL42(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CALS1=AZMAX1(ZALS2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFE11=AZMAX1(ZFE12(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFE21=AZMAX1(ZFE22(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFE31=AZMAX1(ZFE32(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFE41=AZMAX1(ZFE42(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CFES1=AZMAX1(ZFES2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCAO1=AZMAX1(ZCAO2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCAC1=AZMAX1(ZCAC2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCAH1=AZMAX1(ZCAH2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CCAS1=AZMAX1(ZCAS2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CMGO1=AZMAX1(ZMGO2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CMGC1=AZMAX1(ZMGC2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CMGH1=AZMAX1(ZMGH2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CMGS1=AZMAX1(ZMGS2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CNAC1=AZMAX1(ZNAC2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CNAS1=AZMAX1(ZNAS2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    CKAS1=AZMAX1(ZKAS2(N3,N2,N1)/VOLWM(M,N3,N2,N1))
    IF(VLWPA1.GT.ZEROS(N2,N1))THEN
      C0PO41=AZMAX1(H0PO42(N3,N2,N1)/VLWPA1)
      C3PO41=AZMAX1(H3PO42(N3,N2,N1)/VLWPA1)
      CFE1P1=AZMAX1(ZFE1P2(N3,N2,N1)/VLWPA1)
      CFE2P1=AZMAX1(ZFE2P2(N3,N2,N1)/VLWPA1)
      CCA0P1=AZMAX1(ZCA0P2(N3,N2,N1)/VLWPA1)
      CCA1P1=AZMAX1(ZCA1P2(N3,N2,N1)/VLWPA1)
      CCA2P1=AZMAX1(ZCA2P2(N3,N2,N1)/VLWPA1)
      CMG1P1=AZMAX1(ZMG1P2(N3,N2,N1)/VLWPA1)
    ELSE
      C0PO41=0.0_r8
      C3PO41=0.0_r8
      CFE1P1=0.0_r8
      CFE2P1=0.0_r8
      CCA0P1=0.0_r8
      CCA1P1=0.0_r8
      CCA2P1=0.0_r8
      CMG1P1=0.0_r8
    ENDIF
    IF(VLWPB1.GT.ZEROS(N2,N1))THEN
      C0POB1=AZMAX1(H0POB2(N3,N2,N1)/VLWPB1)
      C3POB1=AZMAX1(H3POB2(N3,N2,N1)/VLWPB1)
      CF1PB1=AZMAX1(ZF1PB2(N3,N2,N1)/VLWPB1)
      CF2PB1=AZMAX1(ZF2PB2(N3,N2,N1)/VLWPB1)
      CC0PB1=AZMAX1(ZC0PB2(N3,N2,N1)/VLWPB1)
      CC1PB1=AZMAX1(ZC1PB2(N3,N2,N1)/VLWPB1)
      CC2PB1=AZMAX1(ZC2PB2(N3,N2,N1)/VLWPB1)
      CM1PB1=AZMAX1(ZM1PB2(N3,N2,N1)/VLWPB1)
    ELSE
      C0POB1=C0PO41
      C3POB1=C3PO41
      CF1PB1=CFE1P1
      CF2PB1=CFE2P1
      CC0PB1=CCA0P1
      CC1PB1=CCA1P1
      CC2PB1=CCA2P1
      CM1PB1=CMG1P1
    ENDIF
    CAL2=AZMAX1(ZAL2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFE2=AZMAX1(ZFE2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CHY2=AZMAX1(ZHY2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCA2=AZMAX1(ZCA2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CMG2=AZMAX1(ZMG2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CNA2=AZMAX1(ZNA2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CKA2=AZMAX1(ZKA2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    COH2=AZMAX1(ZOH2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CSO42=AZMAX1(ZSO42(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCL2=AZMAX1(ZCL2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCO32=AZMAX1(ZCO32(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CHCO32=AZMAX1(ZHCO32(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CAL12=AZMAX1(ZAL12(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CAL22=AZMAX1(ZAL22(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CAL32=AZMAX1(ZAL32(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CAL42=AZMAX1(ZAL42(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CALS2=AZMAX1(ZALS2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFE12=AZMAX1(ZFE12(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFE22=AZMAX1(ZFE22(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFE32=AZMAX1(ZFE32(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFE42=AZMAX1(ZFE42(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CFES2=AZMAX1(ZFES2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCAO2=AZMAX1(ZCAO2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCAC2=AZMAX1(ZCAC2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCAH2=AZMAX1(ZCAH2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CCAS2=AZMAX1(ZCAS2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CMGO2=AZMAX1(ZMGO2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CMGC2=AZMAX1(ZMGC2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CMGH2=AZMAX1(ZMGH2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CMGS2=AZMAX1(ZMGS2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CNAC2=AZMAX1(ZNAC2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CNAS2=AZMAX1(ZNAS2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    CKAS2=AZMAX1(ZKAS2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    IF(VLWPA2.GT.ZEROS(N5,N4))THEN
      C0PO42=AZMAX1(H0PO42(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      C3PO42=AZMAX1(H3PO42(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CFE1P2=AZMAX1(ZFE1P2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CFE2P2=AZMAX1(ZFE2P2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CCA0P2=AZMAX1(ZCA0P2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CCA1P2=AZMAX1(ZCA1P2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CCA2P2=AZMAX1(ZCA2P2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CMG1P2=AZMAX1(ZMG1P2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    ELSE
      C0PO42=0.0_r8
      C3PO42=0.0_r8
      CFE1P2=0.0_r8
      CFE2P2=0.0_r8
      CCA0P2=0.0_r8
      CCA1P2=0.0_r8
      CCA2P2=0.0_r8
      CMG1P2=0.0_r8
    ENDIF
    IF(VLWPB2.GT.ZEROS(N5,N4))THEN
      C0POB2=AZMAX1(H0POB2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      C3POB2=AZMAX1(H3POB2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CF1PB2=AZMAX1(ZF1PB2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CF2PB2=AZMAX1(ZF2PB2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CC0PB2=AZMAX1(ZC0PB2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CC1PB2=AZMAX1(ZC1PB2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CC2PB2=AZMAX1(ZC2PB2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
      CM1PB2=AZMAX1(ZM1PB2(N6,N5,N4)/VOLWM(M,N6,N5,N4))
    ELSE
      C0POB2=C0PO42
      C3POB2=C3PO42
      CF1PB2=CFE1P2
      CF2PB2=CFE2P2
      CC0PB2=CCA0P2
      CC1PB2=CCA1P2
      CC2PB2=CCA2P2
      CM1PB2=CMG1P2
    ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MICROPORES
!
!     DLYR=soil layer thickness
!     TORT=micropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     FLWM=water flux through soil micropore from watsub.f
!     DIF*=aqueous diffusivity-dispersivity through micropore
!     *SGL2=solute diffusivity from hour1.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     XDPTH=cross-sectional area/distance between layers
!     C*1,C*2=micropore solute concentration in source,destination layer
!     DFV*=diffusive solute transfer through soil micropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    DLYR1=AMAX1(ZERO2,DLYR(N,N3,N2,N1))
    DLYR2=AMAX1(ZERO2,DLYR(N,N6,N5,N4))
    TORTL=(TORT(M,N3,N2,N1)*DLYR1+TORT(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN=DISP(N,N6,N5,N4)*AMIN1(VFLWX,ABS(FLWM(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))
    DIFPO=(POSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFAL=(ALSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFFE=(FESGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFHY=(HYSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFCA=(CASGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFMG=(GMSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFNA=(ANSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFKA=(AKSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFOH=(OHSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFSO=(SOSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFCL=(CLSXL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFC3=(C3SGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFHC=(HCSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MICROPORES
!
    DFVAL=DIFAL*(CAL1-CAL2)
    DFVFE=DIFFE*(CFE1-CFE2)
    DFVHY=DIFHY*(CHY1-CHY2)
    DFVCA=DIFCA*(CCA1-CCA2)
    DFVMG=DIFMG*(CMG1-CMG2)
    DFVNA=DIFNA*(CNA1-CNA2)
    DFVKA=DIFKA*(CKA1-CKA2)
    DFVOH=DIFOH*(COH1-COH2)
    DFVSO=DIFSO*(CSO41-CSO42)
    DFVCL=DIFCL*(CCL1-CCL2)
    DFVC3=DIFC3*(CCO31-CCO32)
    DFVHC=DIFHC*(CHCO31-CHCO32)
    DFVAL1=DIFAL*(CAL11-CAL12)
    DFVAL2=DIFAL*(CAL21-CAL22)
    DFVAL3=DIFAL*(CAL31-CAL32)
    DFVAL4=DIFAL*(CAL41-CAL42)
    DFVALS=DIFAL*(CALS1-CALS2)
    DFVFE1=DIFFE*(CFE11-CFE12)
    DFVFE2=DIFFE*(CFE21-CFE22)
    DFVFE3=DIFFE*(CFE31-CFE32)
    DFVFE4=DIFFE*(CFE41-CFE42)
    DFVFES=DIFFE*(CFES1-CFES2)
    DFVCAO=DIFCA*(CCAO1-CCAO2)
    DFVCAC=DIFCA*(CCAC1-CCAC2)
    DFVCAH=DIFCA*(CCAH1-CCAH2)
    DFVCAS=DIFCA*(CCAS1-CCAS2)
    DFVMGO=DIFMG*(CMGO1-CMGO2)
    DFVMGC=DIFMG*(CMGC1-CMGC2)
    DFVMGH=DIFMG*(CMGH1-CMGH2)
    DFVMGS=DIFMG*(CMGS1-CMGS2)
    DFVNAC=DIFNA*(CNAC1-CNAC2)
    DFVNAS=DIFNA*(CNAS1-CNAS2)
    DFVKAS=DIFKA*(CKAS1-CKAS2)
    DFVH0P=DIFPO*(C0PO41-C0PO42)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVH3P=DIFPO*(C3PO41-C3PO42)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVF1P=DIFPO*(CFE1P1-CFE1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVF2P=DIFPO*(CFE2P1-CFE2P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC0P=DIFPO*(CCA0P1-CCA0P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC1P=DIFPO*(CCA1P1-CCA1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC2P=DIFPO*(CCA2P1-CCA2P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVM1P=DIFPO*(CMG1P1-CMG1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVH0B=DIFPO*(C0POB1-C0POB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVH3B=DIFPO*(C3POB1-C3POB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVF1B=DIFPO*(CF1PB1-CF1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVF2B=DIFPO*(CF2PB1-CF2PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC0B=DIFPO*(CC0PB1-CC0PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC1B=DIFPO*(CC1PB1-CC1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC2B=DIFPO*(CC2PB1-CC2PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVM1B=DIFPO*(CM1PB1-CM1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
  ELSE
    DFVAL=0.0_r8
    DFVFE=0.0_r8
    DFVHY=0.0_r8
    DFVCA=0.0_r8
    DFVMG=0.0_r8
    DFVNA=0.0_r8
    DFVKA=0.0_r8
    DFVOH=0.0_r8
    DFVSO=0.0_r8
    DFVCL=0.0_r8
    DFVC3=0.0_r8
    DFVHC=0.0_r8
    DFVAL1=0.0_r8
    DFVAL2=0.0_r8
    DFVAL3=0.0_r8
    DFVAL4=0.0_r8
    DFVALS=0.0_r8
    DFVFE1=0.0_r8
    DFVFE2=0.0_r8
    DFVFE3=0.0_r8
    DFVFE4=0.0_r8
    DFVFES=0.0_r8
    DFVCAO=0.0_r8
    DFVCAC=0.0_r8
    DFVCAH=0.0_r8
    DFVCAS=0.0_r8
    DFVMGO=0.0_r8
    DFVMGC=0.0_r8
    DFVMGH=0.0_r8
    DFVMGS=0.0_r8
    DFVNAC=0.0_r8
    DFVNAS=0.0_r8
    DFVKAS=0.0_r8
    DFVH0P=0.0_r8
    DFVH3P=0.0_r8
    DFVF1P=0.0_r8
    DFVF2P=0.0_r8
    DFVC0P=0.0_r8
    DFVC1P=0.0_r8
    DFVC2P=0.0_r8
    DFVM1P=0.0_r8
    DFVH0B=0.0_r8
    DFVH3B=0.0_r8
    DFVF1B=0.0_r8
    DFVF2B=0.0_r8
    DFVC0B=0.0_r8
    DFVC1B=0.0_r8
    DFVC2B=0.0_r8
    DFVM1B=0.0_r8
  ENDIF
!     RFL*=convective flux through micropores
!     DFV*=diffusive solute flux through micropores

  RALFLS(N,N6,N5,N4)=RALFLS(N,N6,N5,N4)+DFVAL
  RFEFLS(N,N6,N5,N4)=RFEFLS(N,N6,N5,N4)+DFVFE
  RHYFLS(N,N6,N5,N4)=RHYFLS(N,N6,N5,N4)+DFVHY
  RCAFLS(N,N6,N5,N4)=RCAFLS(N,N6,N5,N4)+DFVCA
  RMGFLS(N,N6,N5,N4)=RMGFLS(N,N6,N5,N4)+DFVMG
  RNAFLS(N,N6,N5,N4)=RNAFLS(N,N6,N5,N4)+DFVNA
  RKAFLS(N,N6,N5,N4)=RKAFLS(N,N6,N5,N4)+DFVKA
  ROHFLS(N,N6,N5,N4)=ROHFLS(N,N6,N5,N4)+DFVOH
  RSOFLS(N,N6,N5,N4)=RSOFLS(N,N6,N5,N4)+DFVSO
  RCLFLS(N,N6,N5,N4)=RCLFLS(N,N6,N5,N4)+DFVCL
  RC3FLS(N,N6,N5,N4)=RC3FLS(N,N6,N5,N4)+DFVC3
  RHCFLS(N,N6,N5,N4)=RHCFLS(N,N6,N5,N4)+DFVHC
  RAL1FS(N,N6,N5,N4)=RAL1FS(N,N6,N5,N4)+DFVAL1
  RAL2FS(N,N6,N5,N4)=RAL2FS(N,N6,N5,N4)+DFVAL2
  RAL3FS(N,N6,N5,N4)=RAL3FS(N,N6,N5,N4)+DFVAL3
  RAL4FS(N,N6,N5,N4)=RAL4FS(N,N6,N5,N4)+DFVAL4
  RALSFS(N,N6,N5,N4)=RALSFS(N,N6,N5,N4)+DFVALS
  RFE1FS(N,N6,N5,N4)=RFE1FS(N,N6,N5,N4)+DFVFE1
  RFE2FS(N,N6,N5,N4)=RFE2FS(N,N6,N5,N4)+DFVFE2
  RFE3FS(N,N6,N5,N4)=RFE3FS(N,N6,N5,N4)+DFVFE3
  RFE4FS(N,N6,N5,N4)=RFE4FS(N,N6,N5,N4)+DFVFE4
  RFESFS(N,N6,N5,N4)=RFESFS(N,N6,N5,N4)+DFVFES
  RCAOFS(N,N6,N5,N4)=RCAOFS(N,N6,N5,N4)+DFVCAO
  RCACFS(N,N6,N5,N4)=RCACFS(N,N6,N5,N4)+DFVCAC
  RCAHFS(N,N6,N5,N4)=RCAHFS(N,N6,N5,N4)+DFVCAH
  RCASFS(N,N6,N5,N4)=RCASFS(N,N6,N5,N4)+DFVCAS
  RMGOFS(N,N6,N5,N4)=RMGOFS(N,N6,N5,N4)+DFVMGO
  RMGCFS(N,N6,N5,N4)=RMGCFS(N,N6,N5,N4)+DFVMGC
  RMGHFS(N,N6,N5,N4)=RMGHFS(N,N6,N5,N4)+DFVMGH
  RMGSFS(N,N6,N5,N4)=RMGSFS(N,N6,N5,N4)+DFVMGS
  RNACFS(N,N6,N5,N4)=RNACFS(N,N6,N5,N4)+DFVNAC
  RNASFS(N,N6,N5,N4)=RNASFS(N,N6,N5,N4)+DFVNAS
  RKASFS(N,N6,N5,N4)=RKASFS(N,N6,N5,N4)+DFVKAS
  RH0PFS(N,N6,N5,N4)=RH0PFS(N,N6,N5,N4)+DFVH0P
  RH3PFS(N,N6,N5,N4)=RH3PFS(N,N6,N5,N4)+DFVH3P
  RF1PFS(N,N6,N5,N4)=RF1PFS(N,N6,N5,N4)+DFVF1P
  RF2PFS(N,N6,N5,N4)=RF2PFS(N,N6,N5,N4)+DFVF2P
  RC0PFS(N,N6,N5,N4)=RC0PFS(N,N6,N5,N4)+DFVC0P
  RC1PFS(N,N6,N5,N4)=RC1PFS(N,N6,N5,N4)+DFVC1P
  RC2PFS(N,N6,N5,N4)=RC2PFS(N,N6,N5,N4)+DFVC2P
  RM1PFS(N,N6,N5,N4)=RM1PFS(N,N6,N5,N4)+DFVM1P
  RH0BFB(N,N6,N5,N4)=RH0BFB(N,N6,N5,N4)+DFVH0B
  RH3BFB(N,N6,N5,N4)=RH3BFB(N,N6,N5,N4)+DFVH3B
  RF1BFB(N,N6,N5,N4)=RF1BFB(N,N6,N5,N4)+DFVF1B
  RF2BFB(N,N6,N5,N4)=RF2BFB(N,N6,N5,N4)+DFVF2B
  RC0BFB(N,N6,N5,N4)=RC0BFB(N,N6,N5,N4)+DFVC0B
  RC1BFB(N,N6,N5,N4)=RC1BFB(N,N6,N5,N4)+DFVC1B
  RC2BFB(N,N6,N5,N4)=RC2BFB(N,N6,N5,N4)+DFVC2B
  RM1BFB(N,N6,N5,N4)=RM1BFB(N,N6,N5,N4)+DFVM1B

  end subroutine SoluteDifsMicropore
!----------------------------------------------------------------------
  subroutine SoluteAdvMacropore(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  real(r8) :: VFLW
!     FLWHM=water flux through soil macropore from watsub.f
!
  IF(FLWHM(M,N,N6,N5,N4).GT.0.0)THEN
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
!     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN CURRENT GRID CELL
!
!     VOLWHM=macropore water-filled porosity from watsub.f
!     VOLWAH=macropore porosity
!     RFH*=solute diffusive flux through macropore
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *SH2,*BH2=macropore solute content in non-band,band
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    IF(VOLWHM(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FLWHM(M,N,N6,N5,N4)/VOLWHM(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
!
!     ACCOUNT FOR MACROPORE-MICROPORE EXCHANGE IN VERTICAL FLUX
!
    IF(N.EQ.3.AND.VOLAH(N6,N5,N4).GT.VOLWHM(M,N6,N5,N4))THEN
      RFHAL=VFLW*AZMAX1((ZALH2(N3,N2,N1)-AZMIN1(RALFXS(NU(N2,N1),N2,N1))))
      RFHFE=VFLW*AZMAX1((ZFEH2(N3,N2,N1)-AZMIN1(RFEFXS(NU(N2,N1),N2,N1))))
      RFHHY=VFLW*AZMAX1((ZHYH2(N3,N2,N1)-AZMIN1(RHYFXS(NU(N2,N1),N2,N1))))
      RFHCA=VFLW*AZMAX1((ZCCH2(N3,N2,N1)-AZMIN1(RCAFXS(NU(N2,N1),N2,N1))))
      RFHMG=VFLW*AZMAX1((ZMAH2(N3,N2,N1)-AZMIN1(RMGFXS(NU(N2,N1),N2,N1))))
      RFHNA=VFLW*AZMAX1((ZNAH2(N3,N2,N1)-AZMIN1(RNAFXS(NU(N2,N1),N2,N1))))
      RFHKA=VFLW*AZMAX1((ZKAH2(N3,N2,N1)-AZMIN1(RKAFXS(NU(N2,N1),N2,N1))))
      RFHOH=VFLW*AZMAX1((ZOHH2(N3,N2,N1)-AZMIN1(ROHFXS(NU(N2,N1),N2,N1))))
      RFHSO=VFLW*AZMAX1((ZSO4H2(N3,N2,N1)-AZMIN1(RSOFXS(NU(N2,N1),N2,N1))))
      RFHCL=VFLW*AZMAX1((ZCLH2(N3,N2,N1)-AZMIN1(RCLFXS(NU(N2,N1),N2,N1))))
      RFHC3=VFLW*AZMAX1((ZCO3H2(N3,N2,N1)-AZMIN1(RC3FXS(NU(N2,N1),N2,N1))))
      RFHHC=VFLW*AZMAX1((ZHCOH2(N3,N2,N1)-AZMIN1(RHCFXS(NU(N2,N1),N2,N1))))
      RFHAL1=VFLW*AZMAX1((ZAL1H2(N3,N2,N1)-AZMIN1(RAL1XS(NU(N2,N1),N2,N1))))
      RFHAL2=VFLW*AZMAX1((ZAL2H2(N3,N2,N1)-AZMIN1(RAL2XS(NU(N2,N1),N2,N1))))
      RFHAL3=VFLW*AZMAX1((ZAL3H2(N3,N2,N1)-AZMIN1(RAL3XS(NU(N2,N1),N2,N1))))
      RFHAL4=VFLW*AZMAX1((ZAL4H2(N3,N2,N1)-AZMIN1(RAL4XS(NU(N2,N1),N2,N1))))
      RFHALS=VFLW*AZMAX1((ZALSH2(N3,N2,N1)-AZMIN1(RALSXS(NU(N2,N1),N2,N1))))
      RFHFE1=VFLW*AZMAX1((ZFE1H2(N3,N2,N1)-AZMIN1(RFE1XS(NU(N2,N1),N2,N1))))
      RFHFE2=VFLW*AZMAX1((ZFE2H2(N3,N2,N1)-AZMIN1(RFE2XS(NU(N2,N1),N2,N1))))
      RFHFE3=VFLW*AZMAX1((ZFE3H2(N3,N2,N1)-AZMIN1(RFE3XS(NU(N2,N1),N2,N1))))
      RFHFE4=VFLW*AZMAX1((ZFE4H2(N3,N2,N1)-AZMIN1(RFE4XS(NU(N2,N1),N2,N1))))
      RFHFES=VFLW*AZMAX1((ZFESH2(N3,N2,N1)-AZMIN1(RFESXS(NU(N2,N1),N2,N1))))
      RFHCAO=VFLW*AZMAX1((ZCAOH2(N3,N2,N1)-AZMIN1(RCAOXS(NU(N2,N1),N2,N1))))
      RFHCAC=VFLW*AZMAX1((ZCACH2(N3,N2,N1)-AZMIN1(RCACXS(NU(N2,N1),N2,N1))))
      RFHCAH=VFLW*AZMAX1((ZCAHH2(N3,N2,N1)-AZMIN1(RCAHXS(NU(N2,N1),N2,N1))))
      RFHCAS=VFLW*AZMAX1((ZCASH2(N3,N2,N1)-AZMIN1(RCASXS(NU(N2,N1),N2,N1))))
      RFHMGO=VFLW*AZMAX1((ZMGOH2(N3,N2,N1)-AZMIN1(RMGOXS(NU(N2,N1),N2,N1))))
      RFHMGC=VFLW*AZMAX1((ZMGCH2(N3,N2,N1)-AZMIN1(RMGCXS(NU(N2,N1),N2,N1))))
      RFHMGH=VFLW*AZMAX1((ZMGHH2(N3,N2,N1)-AZMIN1(RMGHXS(NU(N2,N1),N2,N1))))
      RFHMGS=VFLW*AZMAX1((ZMGSH2(N3,N2,N1)-AZMIN1(RMGSXS(NU(N2,N1),N2,N1))))
      RFHNAC=VFLW*AZMAX1((ZNACH2(N3,N2,N1)-AZMIN1(RNACXS(NU(N2,N1),N2,N1))))
      RFHNAS=VFLW*AZMAX1((ZNASH2(N3,N2,N1)-AZMIN1(RNASXS(NU(N2,N1),N2,N1))))
      RFHKAS=VFLW*AZMAX1((ZKASH2(N3,N2,N1)-AZMIN1(RKASXS(NU(N2,N1),N2,N1))))
      RFHH0P=VFLW*AZMAX1((H0P4H2(N3,N2,N1)-AZMIN1(RH0PXS(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHH3P=VFLW*AZMAX1((H3P4H2(N3,N2,N1)-AZMIN1(RH3PXS(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHF1P=VFLW*AZMAX1((ZF1PH2(N3,N2,N1)-AZMIN1(RF1PXS(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHF2P=VFLW*AZMAX1((ZF2PH2(N3,N2,N1)-AZMIN1(RF2PXS(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHC0P=VFLW*AZMAX1((ZC0PH2(N3,N2,N1)-AZMIN1(RC0PXS(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHC1P=VFLW*AZMAX1((ZC1PH2(N3,N2,N1)-AZMIN1(RC1PXS(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHC2P=VFLW*AZMAX1((ZC2PH2(N3,N2,N1)-AZMIN1(RC2PXS(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHM1P=VFLW*AZMAX1((ZM1PH2(N3,N2,N1)-AZMIN1(RM1PXS(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4,N3,N2,N1)
      RFHH0B=VFLW*AZMAX1((H0PBH2(N3,N2,N1)-AZMIN1(RH0BXB(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHH3B=VFLW*AZMAX1((H3PBH2(N3,N2,N1)-AZMIN1(RH3BXB(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHF1B=VFLW*AZMAX1((ZF1BH2(N3,N2,N1)-AZMIN1(RF1BXB(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHF2B=VFLW*AZMAX1((ZF2BH2(N3,N2,N1)-AZMIN1(RF2BXB(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHC0B=VFLW*AZMAX1((ZC0BH2(N3,N2,N1)-AZMIN1(RC0BXB(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHC1B=VFLW*AZMAX1((ZC1BH2(N3,N2,N1)-AZMIN1(RC1BXB(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHC2B=VFLW*AZMAX1((ZC2BH2(N3,N2,N1)-AZMIN1(RC2BXB(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
      RFHM1B=VFLW*AZMAX1((ZM1BH2(N3,N2,N1)-AZMIN1(RM1BXB(NU(N2,N1),N2,N1))))*trcs_VLN(ids_H1PO4B,N3,N2,N1)
!
!     OTHERWISE
!
    ELSE
      RFHAL=VFLW*AZMAX1(ZALH2(N3,N2,N1))
      RFHFE=VFLW*AZMAX1(ZFEH2(N3,N2,N1))
      RFHHY=VFLW*AZMAX1(ZHYH2(N3,N2,N1))
      RFHCA=VFLW*AZMAX1(ZCCH2(N3,N2,N1))
      RFHMG=VFLW*AZMAX1(ZMAH2(N3,N2,N1))
      RFHNA=VFLW*AZMAX1(ZNAH2(N3,N2,N1))
      RFHKA=VFLW*AZMAX1(ZKAH2(N3,N2,N1))
      RFHOH=VFLW*AZMAX1(ZOHH2(N3,N2,N1))
      RFHSO=VFLW*AZMAX1(ZSO4H2(N3,N2,N1))
      RFHCL=VFLW*AZMAX1(ZCLH2(N3,N2,N1))
      RFHC3=VFLW*AZMAX1(ZCO3H2(N3,N2,N1))
      RFHHC=VFLW*AZMAX1(ZHCOH2(N3,N2,N1))
      RFHAL1=VFLW*AZMAX1(ZAL1H2(N3,N2,N1))
      RFHAL2=VFLW*AZMAX1(ZAL2H2(N3,N2,N1))
      RFHAL3=VFLW*AZMAX1(ZAL3H2(N3,N2,N1))
      RFHAL4=VFLW*AZMAX1(ZAL4H2(N3,N2,N1))
      RFHALS=VFLW*AZMAX1(ZALSH2(N3,N2,N1))
      RFHFE1=VFLW*AZMAX1(ZFE1H2(N3,N2,N1))
      RFHFE2=VFLW*AZMAX1(ZFE2H2(N3,N2,N1))
      RFHFE3=VFLW*AZMAX1(ZFE3H2(N3,N2,N1))
      RFHFE4=VFLW*AZMAX1(ZFE4H2(N3,N2,N1))
      RFHFES=VFLW*AZMAX1(ZFESH2(N3,N2,N1))
      RFHCAO=VFLW*AZMAX1(ZCAOH2(N3,N2,N1))
      RFHCAC=VFLW*AZMAX1(ZCACH2(N3,N2,N1))
      RFHCAH=VFLW*AZMAX1(ZCAHH2(N3,N2,N1))
      RFHCAS=VFLW*AZMAX1(ZCASH2(N3,N2,N1))
      RFHMGO=VFLW*AZMAX1(ZMGOH2(N3,N2,N1))
      RFHMGC=VFLW*AZMAX1(ZMGCH2(N3,N2,N1))
      RFHMGH=VFLW*AZMAX1(ZMGHH2(N3,N2,N1))
      RFHMGS=VFLW*AZMAX1(ZMGSH2(N3,N2,N1))
      RFHNAC=VFLW*AZMAX1(ZNACH2(N3,N2,N1))
      RFHNAS=VFLW*AZMAX1(ZNASH2(N3,N2,N1))
      RFHKAS=VFLW*AZMAX1(ZKASH2(N3,N2,N1))
      RFHH0P=VFLW*AZMAX1(H0P4H2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHH3P=VFLW*AZMAX1(H3P4H2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHF1P=VFLW*AZMAX1(ZF1PH2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHF2P=VFLW*AZMAX1(ZF2PH2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHC0P=VFLW*AZMAX1(ZC0PH2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHC1P=VFLW*AZMAX1(ZC1PH2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHC2P=VFLW*AZMAX1(ZC2PH2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHM1P=VFLW*AZMAX1(ZM1PH2(N3,N2,N1))*trcs_VLN(ids_H1PO4,N6,N5,N4)
      RFHH0B=VFLW*AZMAX1(H0PBH2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHH3B=VFLW*AZMAX1(H3PBH2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHF1B=VFLW*AZMAX1(ZF1BH2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHF2B=VFLW*AZMAX1(ZF2BH2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHC0B=VFLW*AZMAX1(ZC0BH2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHC1B=VFLW*AZMAX1(ZC1BH2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHC2B=VFLW*AZMAX1(ZC2BH2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
      RFHM1B=VFLW*AZMAX1(ZM1BH2(N3,N2,N1))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    ENDIF
  ELSEIF(FLWHM(M,N,N6,N5,N4).LT.0.0)THEN
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM ADJACENT TO
!     CURRENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE GAS OR SOLUTE CONCENTRATIONS
!     IN ADJACENT GRID CELL
!
    IF(VOLWHM(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FLWHM(M,N,N6,N5,N4)/VOLWHM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    RFHAL=VFLW*AZMAX1(ZALH2(N6,N5,N4))
    RFHFE=VFLW*AZMAX1(ZFEH2(N6,N5,N4))
    RFHHY=VFLW*AZMAX1(ZHYH2(N6,N5,N4))
    RFHCA=VFLW*AZMAX1(ZCCH2(N6,N5,N4))
    RFHMG=VFLW*AZMAX1(ZMAH2(N6,N5,N4))
    RFHNA=VFLW*AZMAX1(ZNAH2(N6,N5,N4))
    RFHKA=VFLW*AZMAX1(ZKAH2(N6,N5,N4))
    RFHOH=VFLW*AZMAX1(ZOHH2(N6,N5,N4))
    RFHSO=VFLW*AZMAX1(ZSO4H2(N6,N5,N4))
    RFHCL=VFLW*AZMAX1(ZCLH2(N6,N5,N4))
    RFHC3=VFLW*AZMAX1(ZCO3H2(N6,N5,N4))
    RFHHC=VFLW*AZMAX1(ZHCOH2(N6,N5,N4))
    RFHAL1=VFLW*AZMAX1(ZAL1H2(N6,N5,N4))
    RFHAL2=VFLW*AZMAX1(ZAL2H2(N6,N5,N4))
    RFHAL3=VFLW*AZMAX1(ZAL3H2(N6,N5,N4))
    RFHAL4=VFLW*AZMAX1(ZAL4H2(N6,N5,N4))
    RFHALS=VFLW*AZMAX1(ZALSH2(N6,N5,N4))
    RFHFE1=VFLW*AZMAX1(ZFE1H2(N6,N5,N4))
    RFHFE2=VFLW*AZMAX1(ZFE2H2(N6,N5,N4))
    RFHFE3=VFLW*AZMAX1(ZFE3H2(N6,N5,N4))
    RFHFE4=VFLW*AZMAX1(ZFE4H2(N6,N5,N4))
    RFHFES=VFLW*AZMAX1(ZFESH2(N6,N5,N4))
    RFHCAO=VFLW*AZMAX1(ZCAOH2(N6,N5,N4))
    RFHCAC=VFLW*AZMAX1(ZCACH2(N6,N5,N4))
    RFHCAH=VFLW*AZMAX1(ZCAHH2(N6,N5,N4))
    RFHCAS=VFLW*AZMAX1(ZCASH2(N6,N5,N4))
    RFHMGO=VFLW*AZMAX1(ZMGOH2(N6,N5,N4))
    RFHMGC=VFLW*AZMAX1(ZMGCH2(N6,N5,N4))
    RFHMGH=VFLW*AZMAX1(ZMGHH2(N6,N5,N4))
    RFHMGS=VFLW*AZMAX1(ZMGSH2(N6,N5,N4))
    RFHNAC=VFLW*AZMAX1(ZNACH2(N6,N5,N4))
    RFHNAS=VFLW*AZMAX1(ZNASH2(N6,N5,N4))
    RFHKAS=VFLW*AZMAX1(ZKASH2(N6,N5,N4))
    RFHH0P=VFLW*AZMAX1(H0P4H2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHH3P=VFLW*AZMAX1(H3P4H2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHF1P=VFLW*AZMAX1(ZF1PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHF2P=VFLW*AZMAX1(ZF2PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHC0P=VFLW*AZMAX1(ZC0PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHC1P=VFLW*AZMAX1(ZC1PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHC2P=VFLW*AZMAX1(ZC2PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHM1P=VFLW*AZMAX1(ZM1PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFHH0B=VFLW*AZMAX1(H0PBH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHH3B=VFLW*AZMAX1(H3PBH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHF1B=VFLW*AZMAX1(ZF1BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHF2B=VFLW*AZMAX1(ZF2BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHC0B=VFLW*AZMAX1(ZC0BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHC1B=VFLW*AZMAX1(ZC1BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHC2B=VFLW*AZMAX1(ZC2BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFHM1B=VFLW*AZMAX1(ZM1BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
!
!     NO MACROPORE FLUX
!
  ELSE
    RFHAL=0.0_r8
    RFHFE=0.0_r8
    RFHHY=0.0_r8
    RFHCA=0.0_r8
    RFHMG=0.0_r8
    RFHNA=0.0_r8
    RFHKA=0.0_r8
    RFHOH=0.0_r8
    RFHSO=0.0_r8
    RFHCL=0.0_r8
    RFHC3=0.0_r8
    RFHHC=0.0_r8
    RFHAL1=0.0_r8
    RFHAL2=0.0_r8
    RFHAL3=0.0_r8
    RFHAL4=0.0_r8
    RFHALS=0.0_r8
    RFHFE1=0.0_r8
    RFHFE2=0.0_r8
    RFHFE3=0.0_r8
    RFHFE4=0.0_r8
    RFHFES=0.0_r8
    RFHCAO=0.0_r8
    RFHCAC=0.0_r8
    RFHCAH=0.0_r8
    RFHCAS=0.0_r8
    RFHMGO=0.0_r8
    RFHMGC=0.0_r8
    RFHMGH=0.0_r8
    RFHMGS=0.0_r8
    RFHNAC=0.0_r8
    RFHNAS=0.0_r8
    RFHKAS=0.0_r8
    RFHH0P=0.0_r8
    RFHH3P=0.0_r8
    RFHF1P=0.0_r8
    RFHF2P=0.0_r8
    RFHC0P=0.0_r8
    RFHC1P=0.0_r8
    RFHC2P=0.0_r8
    RFHM1P=0.0_r8
    RFHH0B=0.0_r8
    RFHH3B=0.0_r8
    RFHF1B=0.0_r8
    RFHF2B=0.0_r8
    RFHC0B=0.0_r8
    RFHC1B=0.0_r8
    RFHC2B=0.0_r8
    RFHM1B=0.0_r8
  ENDIF
!     RFH*=convective flux through macropores
!     DFH*=diffusive solute flux through macropores
!

  RALFHS(N,N6,N5,N4)=RFHAL
  RFEFHS(N,N6,N5,N4)=RFHFE
  RHYFHS(N,N6,N5,N4)=RFHHY
  RCAFHS(N,N6,N5,N4)=RFHCA
  RMGFHS(N,N6,N5,N4)=RFHMG
  RNAFHS(N,N6,N5,N4)=RFHNA
  RKAFHS(N,N6,N5,N4)=RFHKA
  ROHFHS(N,N6,N5,N4)=RFHOH
  RSOFHS(N,N6,N5,N4)=RFHSO
  RCLFHS(N,N6,N5,N4)=RFHCL
  RC3FHS(N,N6,N5,N4)=RFHC3
  RHCFHS(N,N6,N5,N4)=RFHHC
  RAL1HS(N,N6,N5,N4)=RFHAL1
  RAL2HS(N,N6,N5,N4)=RFHAL2
  RAL3HS(N,N6,N5,N4)=RFHAL3
  RAL4HS(N,N6,N5,N4)=RFHAL4
  RALSHS(N,N6,N5,N4)=RFHALS
  RFE1HS(N,N6,N5,N4)=RFHFE1
  RFE2HS(N,N6,N5,N4)=RFHFE2
  RFE3HS(N,N6,N5,N4)=RFHFE3
  RFE4HS(N,N6,N5,N4)=RFHFE4
  RFESHS(N,N6,N5,N4)=RFHFES
  RCAOHS(N,N6,N5,N4)=RFHCAO
  RCACHS(N,N6,N5,N4)=RFHCAC
  RCAHHS(N,N6,N5,N4)=RFHCAH
  RCASHS(N,N6,N5,N4)=RFHCAS
  RMGOHS(N,N6,N5,N4)=RFHMGO
  RMGCHS(N,N6,N5,N4)=RFHMGC
  RMGHHS(N,N6,N5,N4)=RFHMGH
  RMGSHS(N,N6,N5,N4)=RFHMGS
  RNACHS(N,N6,N5,N4)=RFHNAC
  RNASHS(N,N6,N5,N4)=RFHNAS
  RKASHS(N,N6,N5,N4)=RFHKAS
  RH0PHS(N,N6,N5,N4)=RFHH0P
  RH3PHS(N,N6,N5,N4)=RFHH3P
  RF1PHS(N,N6,N5,N4)=RFHF1P
  RF2PHS(N,N6,N5,N4)=RFHF2P
  RC0PHS(N,N6,N5,N4)=RFHC0P
  RC1PHS(N,N6,N5,N4)=RFHC1P
  RC2PHS(N,N6,N5,N4)=RFHC2P
  RM1PHS(N,N6,N5,N4)=RFHM1P
  RH0BHB(N,N6,N5,N4)=RFHH0B
  RH3BHB(N,N6,N5,N4)=RFHH3B
  RF1BHB(N,N6,N5,N4)=RFHF1B
  RF2BHB(N,N6,N5,N4)=RFHF2B
  RC0BHB(N,N6,N5,N4)=RFHC0B
  RC1BHB(N,N6,N5,N4)=RFHC1B
  RC2BHB(N,N6,N5,N4)=RFHC2B
  RM1BHB(N,N6,N5,N4)=RFHM1B

  end subroutine SoluteAdvMacropore

!----------------------------------------------------------------------
  subroutine SoluteDifsMacropore(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

!     DIFFUSIVE FLUXES OF GASES AND SOLUTES BETWEEN CURRENT AND
!     ADJACENT GRID CELL MACROPORES FROM AQUEOUS DIFFUSIVITIES
!     AND CONCENTRATION DIFFERENCES
!
!     VOLWHM=macropore water-filled porosity from watsub.f
!     THETY=hygroscopic water content
!     VOLAH=total macropore volume

  IF(VOLWHM(M,N3,N2,N1).GT.THETY(N3,N2,N1)*VOLAH(N3,N2,N1) &
    .AND.VOLWHM(M,N6,N5,N4).GT.THETY(N6,N5,N4)*VOLAH(N6,N5,N4))THEN
!
!     MACROPORE CONCENTRATIONS IN CURRENT AND ADJACENT GRID CELLS
!
!     C*H1,C*H2=macropore solute concentration in source,destination layer
!     *H2=macropore solute content
!     VOLWHM=macropore water content
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
    CAL1=AZMAX1(ZALH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE1=AZMAX1(ZFEH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CHY1=AZMAX1(ZHYH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCA1=AZMAX1(ZCCH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMG1=AZMAX1(ZMAH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CNA1=AZMAX1(ZNAH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CKA1=AZMAX1(ZKAH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    COH1=AZMAX1(ZOHH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CSO41=AZMAX1(ZSO4H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCL1=AZMAX1(ZCLH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCO31=AZMAX1(ZCO3H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CHCO31=AZMAX1(ZHCOH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CAL11=AZMAX1(ZAL1H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CAL21=AZMAX1(ZAL2H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CAL31=AZMAX1(ZAL3H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CAL41=AZMAX1(ZAL4H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CALS1=AZMAX1(ZALSH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE11=AZMAX1(ZFE1H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE21=AZMAX1(ZFE2H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE31=AZMAX1(ZFE3H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE41=AZMAX1(ZFE4H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFES1=AZMAX1(ZFESH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCAO1=AZMAX1(ZCAOH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCAC1=AZMAX1(ZCACH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCAH1=AZMAX1(ZCAHH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCAS1=AZMAX1(ZCASH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMGO1=AZMAX1(ZMGOH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMGC1=AZMAX1(ZMGCH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMGH1=AZMAX1(ZMGHH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMGS1=AZMAX1(ZMGSH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CNAC1=AZMAX1(ZNACH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CNAS1=AZMAX1(ZNASH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CKAS1=AZMAX1(ZKASH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    C0PO41=AZMAX1(H0P4H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    C3PO41=AZMAX1(H3P4H2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE1P1=AZMAX1(ZF1PH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CFE2P1=AZMAX1(ZF2PH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCA0P1=AZMAX1(ZC0PH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCA1P1=AZMAX1(ZC1PH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CCA2P1=AZMAX1(ZC2PH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CMG1P1=AZMAX1(ZM1PH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    C0POB1=AZMAX1(H0PBH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    C3POB1=AZMAX1(H3PBH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CF1PB1=AZMAX1(ZF1BH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CF2PB1=AZMAX1(ZF2BH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CC0PB1=AZMAX1(ZC0BH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CC1PB1=AZMAX1(ZC1BH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CC2PB1=AZMAX1(ZC2BH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CM1PB1=AZMAX1(ZM1BH2(N3,N2,N1)/VOLWHM(M,N3,N2,N1))
    CAL2=AZMAX1(ZALH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE2=AZMAX1(ZFEH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CHY2=AZMAX1(ZHYH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCA2=AZMAX1(ZCCH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMG2=AZMAX1(ZMAH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CNA2=AZMAX1(ZNAH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CKA2=AZMAX1(ZKAH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    COH2=AZMAX1(ZOHH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CSO42=AZMAX1(ZSO4H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCL2=AZMAX1(ZCLH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCO32=AZMAX1(ZCO3H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CHCO32=AZMAX1(ZHCOH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CAL12=AZMAX1(ZAL1H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CAL22=AZMAX1(ZAL2H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CAL32=AZMAX1(ZAL3H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CAL42=AZMAX1(ZAL4H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CALS2=AZMAX1(ZALSH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE12=AZMAX1(ZFE1H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE22=AZMAX1(ZFE2H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE32=AZMAX1(ZFE3H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE42=AZMAX1(ZFE4H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFES2=AZMAX1(ZFESH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCAO2=AZMAX1(ZCAOH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCAC2=AZMAX1(ZCACH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCAH2=AZMAX1(ZCAHH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCAS2=AZMAX1(ZCASH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMGO2=AZMAX1(ZMGOH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMGC2=AZMAX1(ZMGCH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMGH2=AZMAX1(ZMGHH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMGS2=AZMAX1(ZMGSH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CNAC2=AZMAX1(ZNACH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CNAS2=AZMAX1(ZNASH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CKAS2=AZMAX1(ZKASH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    C0PO42=AZMAX1(H0P4H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    C3PO42=AZMAX1(H3P4H2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE1P2=AZMAX1(ZF1PH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CFE2P2=AZMAX1(ZF2PH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCA0P2=AZMAX1(ZC0PH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCA1P2=AZMAX1(ZC1PH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CCA2P2=AZMAX1(ZC2PH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CMG1P2=AZMAX1(ZM1PH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    C0POB2=AZMAX1(H0PBH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    C3POB2=AZMAX1(H3PBH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CF1PB2=AZMAX1(ZF1BH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CF2PB2=AZMAX1(ZF2BH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CC0PB2=AZMAX1(ZC0BH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CC1PB2=AZMAX1(ZC1BH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CC2PB2=AZMAX1(ZC2BH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
    CM1PB2=AZMAX1(ZM1BH2(N6,N5,N4)/VOLWHM(M,N6,N5,N4))
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MACROPORES
!
!     DLYR=soil layer thickness
!     TORTH=macropore tortuosity from hour1.f
!     DISP=dispersivity parameter
!     FLWHM=water flux through soil macropore from watsub.f
!     DIF*=aqueous diffusivity-dispersivity through macropore
!     *SGL2=solute diffusivity from hour1.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     XDPTH=cross-sectional area/distance between layers
!     C*H1,C*H2=macropore solute concentration in source,destination layer
!     DFH*=diffusive solute transfer through soil macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!
    DLYR1=AMAX1(ZERO2,DLYR(N,N3,N2,N1))
    DLYR2=AMAX1(ZERO2,DLYR(N,N6,N5,N4))
    TORTL=(TORTH(M,N3,N2,N1)*DLYR1+TORTH(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    DISPN=DISP(N,N6,N5,N4)*AMIN1(VFLWX,ABS(FLWHM(M,N,N6,N5,N4)/AREA(N,N6,N5,N4)))
    DIFPO=(POSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFAL=(ALSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFFE=(FESGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFHY=(HYSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFCA=(CASGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFMG=(GMSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFNA=(ANSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFKA=(AKSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFOH=(OHSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFSO=(SOSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFCL=(CLSXL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFC3=(C3SGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
    DIFHC=(HCSGL2(N6,N5,N4)*TORTL+DISPN)*XDPTH(N,N6,N5,N4)
!
!     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
!     MACROPORES
!
    DFHAL=DIFAL*(CAL1-CAL2)
    DFHFE=DIFFE*(CFE1-CFE2)
    DFHHY=DIFHY*(CHY1-CHY2)
    DFHCA=DIFCA*(CCA1-CCA2)
    DFHMG=DIFMG*(CMG1-CMG2)
    DFHNA=DIFNA*(CNA1-CNA2)
    DFHKA=DIFKA*(CKA1-CKA2)
    DFHOH=DIFOH*(COH1-COH2)
    DFHSO=DIFSO*(CSO41-CSO42)
    DFHCL=DIFCL*(CCL1-CCL2)
    DFHC3=DIFC3*(CCO31-CCO32)
    DFHHC=DIFHC*(CHCO31-CHCO32)
    DFHAL1=DIFAL*(CAL11-CAL12)
    DFHAL2=DIFAL*(CAL21-CAL22)
    DFHAL3=DIFAL*(CAL31-CAL32)
    DFHAL4=DIFAL*(CAL41-CAL42)
    DFHALS=DIFAL*(CALS1-CALS2)
    DFHFE1=DIFFE*(CFE11-CFE12)
    DFHFE2=DIFFE*(CFE21-CFE22)
    DFHFE3=DIFFE*(CFE31-CFE32)
    DFHFE4=DIFFE*(CFE41-CFE42)
    DFHFES=DIFFE*(CFES1-CFES2)
    DFHCAO=DIFCA*(CCAO1-CCAO2)
    DFHCAC=DIFCA*(CCAC1-CCAC2)
    DFHCAH=DIFCA*(CCAH1-CCAH2)
    DFHCAS=DIFCA*(CCAS1-CCAS2)
    DFHMGO=DIFMG*(CMGO1-CMGO2)
    DFHMGC=DIFMG*(CMGC1-CMGC2)
    DFHMGH=DIFMG*(CMGH1-CMGH2)
    DFHMGS=DIFMG*(CMGS1-CMGS2)
    DFHNAC=DIFNA*(CNAC1-CNAC2)
    DFHNAS=DIFNA*(CNAS1-CNAS2)
    DFHKAS=DIFKA*(CKAS1-CKAS2)
    DFHH0P=DIFPO*(C0PO41-C0PO42)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHH3P=DIFPO*(C3PO41-C3PO42)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHF1P=DIFPO*(CFE1P1-CFE1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHF2P=DIFPO*(CFE2P1-CFE2P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHC0P=DIFPO*(CCA0P1-CCA0P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHC1P=DIFPO*(CCA1P1-CCA1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHC2P=DIFPO*(CCA2P1-CCA2P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHM1P=DIFPO*(CMG1P1-CMG1P2)*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFHH0B=DIFPO*(C0POB1-C0POB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHH3B=DIFPO*(C3POB1-C3POB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHF1B=DIFPO*(CF1PB1-CF1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHF2B=DIFPO*(CF2PB1-CF2PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHC0B=DIFPO*(CC0PB1-CC0PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHC1B=DIFPO*(CC1PB1-CC1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHC2B=DIFPO*(CC2PB1-CC2PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFHM1B=DIFPO*(CM1PB1-CM1PB2)*trcs_VLN(ids_H1PO4B,N6,N5,N4)
  ELSE
    DFHAL=0.0_r8
    DFHFE=0.0_r8
    DFHHY=0.0_r8
    DFHCA=0.0_r8
    DFHMG=0.0_r8
    DFHNA=0.0_r8
    DFHKA=0.0_r8
    DFHOH=0.0_r8
    DFHSO=0.0_r8
    DFHCL=0.0_r8
    DFHC3=0.0_r8
    DFHHC=0.0_r8
    DFHAL1=0.0_r8
    DFHAL2=0.0_r8
    DFHAL3=0.0_r8
    DFHAL4=0.0_r8
    DFHALS=0.0_r8
    DFHFE1=0.0_r8
    DFHFE2=0.0_r8
    DFHFE3=0.0_r8
    DFHFE4=0.0_r8
    DFHFES=0.0_r8
    DFHCAO=0.0_r8
    DFHCAC=0.0_r8
    DFHCAH=0.0_r8
    DFHCAS=0.0_r8
    DFHMGO=0.0_r8
    DFHMGC=0.0_r8
    DFHMGH=0.0_r8
    DFHMGS=0.0_r8
    DFHNAC=0.0_r8
    DFHNAS=0.0_r8
    DFHKAS=0.0_r8
    DFHH0P=0.0_r8
    DFHH3P=0.0_r8
    DFHF1P=0.0_r8
    DFHF2P=0.0_r8
    DFHC0P=0.0_r8
    DFHC1P=0.0_r8
    DFHC2P=0.0_r8
    DFHM1P=0.0_r8
    DFHH0B=0.0_r8
    DFHH3B=0.0_r8
    DFHF1B=0.0_r8
    DFHF2B=0.0_r8
    DFHC0B=0.0_r8
    DFHC1B=0.0_r8
    DFHC2B=0.0_r8
    DFHM1B=0.0_r8
  ENDIF
!     RFH*=convective flux through macropores
!     DFH*=diffusive solute flux through macropores
!

  RALFHS(N,N6,N5,N4)=RALFHS(N,N6,N5,N4)+DFHAL
  RFEFHS(N,N6,N5,N4)=RFEFHS(N,N6,N5,N4)+DFHFE
  RHYFHS(N,N6,N5,N4)=RHYFHS(N,N6,N5,N4)+DFHHY
  RCAFHS(N,N6,N5,N4)=RCAFHS(N,N6,N5,N4)+DFHCA
  RMGFHS(N,N6,N5,N4)=RMGFHS(N,N6,N5,N4)+DFHMG
  RNAFHS(N,N6,N5,N4)=RNAFHS(N,N6,N5,N4)+DFHNA
  RKAFHS(N,N6,N5,N4)=RKAFHS(N,N6,N5,N4)+DFHKA
  ROHFHS(N,N6,N5,N4)=ROHFHS(N,N6,N5,N4)+DFHOH
  RSOFHS(N,N6,N5,N4)=RSOFHS(N,N6,N5,N4)+DFHSO
  RCLFHS(N,N6,N5,N4)=RCLFHS(N,N6,N5,N4)+DFHCL
  RC3FHS(N,N6,N5,N4)=RC3FHS(N,N6,N5,N4)+DFHC3
  RHCFHS(N,N6,N5,N4)=RHCFHS(N,N6,N5,N4)+DFHHC
  RAL1HS(N,N6,N5,N4)=RAL1HS(N,N6,N5,N4)+DFHAL1
  RAL2HS(N,N6,N5,N4)=RAL2HS(N,N6,N5,N4)+DFHAL2
  RAL3HS(N,N6,N5,N4)=RAL3HS(N,N6,N5,N4)+DFHAL3
  RAL4HS(N,N6,N5,N4)=RAL4HS(N,N6,N5,N4)+DFHAL4
  RALSHS(N,N6,N5,N4)=RALSHS(N,N6,N5,N4)+DFHALS
  RFE1HS(N,N6,N5,N4)=RFE1HS(N,N6,N5,N4)+DFHFE1
  RFE2HS(N,N6,N5,N4)=RFE2HS(N,N6,N5,N4)+DFHFE2
  RFE3HS(N,N6,N5,N4)=RFE3HS(N,N6,N5,N4)+DFHFE3
  RFE4HS(N,N6,N5,N4)=RFE4HS(N,N6,N5,N4)+DFHFE4
  RFESHS(N,N6,N5,N4)=RFESHS(N,N6,N5,N4)+DFHFES
  RCAOHS(N,N6,N5,N4)=RCAOHS(N,N6,N5,N4)+DFHCAO
  RCACHS(N,N6,N5,N4)=RCACHS(N,N6,N5,N4)+DFHCAC
  RCAHHS(N,N6,N5,N4)=RCAHHS(N,N6,N5,N4)+DFHCAH
  RCASHS(N,N6,N5,N4)=RCASHS(N,N6,N5,N4)+DFHCAS
  RMGOHS(N,N6,N5,N4)=RMGOHS(N,N6,N5,N4)+DFHMGO
  RMGCHS(N,N6,N5,N4)=RMGCHS(N,N6,N5,N4)+DFHMGC
  RMGHHS(N,N6,N5,N4)=RMGHHS(N,N6,N5,N4)+DFHMGH
  RMGSHS(N,N6,N5,N4)=RMGSHS(N,N6,N5,N4)+DFHMGS
  RNACHS(N,N6,N5,N4)=RNACHS(N,N6,N5,N4)+DFHNAC
  RNASHS(N,N6,N5,N4)=RNASHS(N,N6,N5,N4)+DFHNAS
  RKASHS(N,N6,N5,N4)=RKASHS(N,N6,N5,N4)+DFHKAS
  RH0PHS(N,N6,N5,N4)=RH0PHS(N,N6,N5,N4)+DFHH0P
  RH3PHS(N,N6,N5,N4)=RH3PHS(N,N6,N5,N4)+DFHH3P
  RF1PHS(N,N6,N5,N4)=RF1PHS(N,N6,N5,N4)+DFHF1P
  RF2PHS(N,N6,N5,N4)=RF2PHS(N,N6,N5,N4)+DFHF2P
  RC0PHS(N,N6,N5,N4)=RC0PHS(N,N6,N5,N4)+DFHC0P
  RC1PHS(N,N6,N5,N4)=RC1PHS(N,N6,N5,N4)+DFHC1P
  RC2PHS(N,N6,N5,N4)=RC2PHS(N,N6,N5,N4)+DFHC2P
  RM1PHS(N,N6,N5,N4)=RM1PHS(N,N6,N5,N4)+DFHM1P
  RH0BHB(N,N6,N5,N4)=RH0BHB(N,N6,N5,N4)+DFHH0B
  RH3BHB(N,N6,N5,N4)=RH3BHB(N,N6,N5,N4)+DFHH3B
  RF1BHB(N,N6,N5,N4)=RF1BHB(N,N6,N5,N4)+DFHF1B
  RF2BHB(N,N6,N5,N4)=RF2BHB(N,N6,N5,N4)+DFHF2B
  RC0BHB(N,N6,N5,N4)=RC0BHB(N,N6,N5,N4)+DFHC0B
  RC1BHB(N,N6,N5,N4)=RC1BHB(N,N6,N5,N4)+DFHC1B
  RC2BHB(N,N6,N5,N4)=RC2BHB(N,N6,N5,N4)+DFHC2B
  RM1BHB(N,N6,N5,N4)=RM1BHB(N,N6,N5,N4)+DFHM1B
  end subroutine SoluteDifsMacropore
!----------------------------------------------------------------------
  subroutine SoluteAdvExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6

  real(r8) :: VFLW
!     MACROPORE-MICROPORE CONVECTIVE SOLUTE EXCHANGE IN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
!     FINHM=macro-micropore water transfer from watsub.f
!     VOLWM,VOLWHM=micropore,macropore water volume
!     RFL*=convective macropore-micropore solute transfer
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *H2,*2=macropore,micropore solute content
!
!     MACROPORE TO MICROPORE TRANSFER
!
  IF(FINHM(M,N6,N5,N4).GT.0.0)THEN
    IF(VOLWHM(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FINHM(M,N6,N5,N4)/VOLWHM(M,N6,N5,N4)))
    ELSE
      VFLW=VFLWX
    ENDIF
    RFLAL=VFLW*AZMAX1(ZALH2(N6,N5,N4))
    RFLFE=VFLW*AZMAX1(ZFEH2(N6,N5,N4))
    RFLHY=VFLW*AZMAX1(ZHYH2(N6,N5,N4))
    RFLCA=VFLW*AZMAX1(ZCCH2(N6,N5,N4))
    RFLMG=VFLW*AZMAX1(ZMAH2(N6,N5,N4))
    RFLNA=VFLW*AZMAX1(ZNAH2(N6,N5,N4))
    RFLKA=VFLW*AZMAX1(ZKAH2(N6,N5,N4))
    RFLOH=VFLW*AZMAX1(ZOHH2(N6,N5,N4))
    RFLSO=VFLW*AZMAX1(ZSO4H2(N6,N5,N4))
    RFLCL=VFLW*AZMAX1(ZCLH2(N6,N5,N4))
    RFLC3=VFLW*AZMAX1(ZCO3H2(N6,N5,N4))
    RFLHC=VFLW*AZMAX1(ZHCOH2(N6,N5,N4))
    RFLAL1=VFLW*AZMAX1(ZAL1H2(N6,N5,N4))
    RFLAL2=VFLW*AZMAX1(ZAL2H2(N6,N5,N4))
    RFLAL3=VFLW*AZMAX1(ZAL3H2(N6,N5,N4))
    RFLAL4=VFLW*AZMAX1(ZAL4H2(N6,N5,N4))
    RFLALS=VFLW*AZMAX1(ZALSH2(N6,N5,N4))
    RFLFE1=VFLW*AZMAX1(ZFE1H2(N6,N5,N4))
    RFLFE2=VFLW*AZMAX1(ZFE2H2(N6,N5,N4))
    RFLFE3=VFLW*AZMAX1(ZFE3H2(N6,N5,N4))
    RFLFE4=VFLW*AZMAX1(ZFE4H2(N6,N5,N4))
    RFLFES=VFLW*AZMAX1(ZFESH2(N6,N5,N4))
    RFLCAO=VFLW*AZMAX1(ZCAOH2(N6,N5,N4))
    RFLCAC=VFLW*AZMAX1(ZCACH2(N6,N5,N4))
    RFLCAH=VFLW*AZMAX1(ZCAHH2(N6,N5,N4))
    RFLCAS=VFLW*AZMAX1(ZCASH2(N6,N5,N4))
    RFLMGO=VFLW*AZMAX1(ZMGOH2(N6,N5,N4))
    RFLMGC=VFLW*AZMAX1(ZMGCH2(N6,N5,N4))
    RFLMGH=VFLW*AZMAX1(ZMGHH2(N6,N5,N4))
    RFLMGS=VFLW*AZMAX1(ZMGSH2(N6,N5,N4))
    RFLNAC=VFLW*AZMAX1(ZNACH2(N6,N5,N4))
    RFLNAS=VFLW*AZMAX1(ZNASH2(N6,N5,N4))
    RFLKAS=VFLW*AZMAX1(ZKASH2(N6,N5,N4))
    RFLH0P=VFLW*AZMAX1(H0P4H2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLH3P=VFLW*AZMAX1(H3P4H2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLF1P=VFLW*AZMAX1(ZF1PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLF2P=VFLW*AZMAX1(ZF2PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC0P=VFLW*AZMAX1(ZC0PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC1P=VFLW*AZMAX1(ZC1PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC2P=VFLW*AZMAX1(ZC2PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLM1P=VFLW*AZMAX1(ZM1PH2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLH0B=VFLW*AZMAX1(H0PBH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLH3B=VFLW*AZMAX1(H3PBH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLF1B=VFLW*AZMAX1(ZF1BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLF2B=VFLW*AZMAX1(ZF2BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC0B=VFLW*AZMAX1(ZC0BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC1B=VFLW*AZMAX1(ZC1BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC2B=VFLW*AZMAX1(ZC2BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLM1B=VFLW*AZMAX1(ZM1BH2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
!
!     MICROPORE TO MACROPORE TRANSFER
!
  ELSEIF(FINHM(M,N6,N5,N4).LT.0.0)THEN
    IF(VOLWM(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FINHM(M,N6,N5,N4)/VOLWM(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    RFLAL=VFLW*AZMAX1(ZAL2(N6,N5,N4))
    RFLFE=VFLW*AZMAX1(ZFE2(N6,N5,N4))
    RFLHY=VFLW*AZMAX1(ZHY2(N6,N5,N4))
    RFLCA=VFLW*AZMAX1(ZCA2(N6,N5,N4))
    RFLMG=VFLW*AZMAX1(ZMG2(N6,N5,N4))
    RFLNA=VFLW*AZMAX1(ZNA2(N6,N5,N4))
    RFLKA=VFLW*AZMAX1(ZKA2(N6,N5,N4))
    RFLOH=VFLW*AZMAX1(ZOH2(N6,N5,N4))
    RFLSO=VFLW*AZMAX1(ZSO42(N6,N5,N4))
    RFLCL=VFLW*AZMAX1(ZCL2(N6,N5,N4))
    RFLC3=VFLW*AZMAX1(ZCO32(N6,N5,N4))
    RFLHC=VFLW*AZMAX1(ZHCO32(N6,N5,N4))
    RFLAL1=VFLW*AZMAX1(ZAL12(N6,N5,N4))
    RFLAL2=VFLW*AZMAX1(ZAL22(N6,N5,N4))
    RFLAL3=VFLW*AZMAX1(ZAL32(N6,N5,N4))
    RFLAL4=VFLW*AZMAX1(ZAL42(N6,N5,N4))
    RFLALS=VFLW*AZMAX1(ZALS2(N6,N5,N4))
    RFLFE1=VFLW*AZMAX1(ZFE12(N6,N5,N4))
    RFLFE2=VFLW*AZMAX1(ZFE22(N6,N5,N4))
    RFLFE3=VFLW*AZMAX1(ZFE32(N6,N5,N4))
    RFLFE4=VFLW*AZMAX1(ZFE42(N6,N5,N4))
    RFLFES=VFLW*AZMAX1(ZFES2(N6,N5,N4))
    RFLCAO=VFLW*AZMAX1(ZCAO2(N6,N5,N4))
    RFLCAC=VFLW*AZMAX1(ZCAC2(N6,N5,N4))
    RFLCAH=VFLW*AZMAX1(ZCAH2(N6,N5,N4))
    RFLCAS=VFLW*AZMAX1(ZCAS2(N6,N5,N4))
    RFLMGO=VFLW*AZMAX1(ZMGO2(N6,N5,N4))
    RFLMGC=VFLW*AZMAX1(ZMGC2(N6,N5,N4))
    RFLMGH=VFLW*AZMAX1(ZMGH2(N6,N5,N4))
    RFLMGS=VFLW*AZMAX1(ZMGS2(N6,N5,N4))
    RFLNAC=VFLW*AZMAX1(ZNAC2(N6,N5,N4))
    RFLNAS=VFLW*AZMAX1(ZNAS2(N6,N5,N4))
    RFLKAS=VFLW*AZMAX1(ZKAS2(N6,N5,N4))
    RFLH0P=VFLW*AZMAX1(H0PO42(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLH3P=VFLW*AZMAX1(H3PO42(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLF1P=VFLW*AZMAX1(ZFE1P2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLF2P=VFLW*AZMAX1(ZFE2P2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC0P=VFLW*AZMAX1(ZCA0P2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC1P=VFLW*AZMAX1(ZCA1P2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLC2P=VFLW*AZMAX1(ZCA2P2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLM1P=VFLW*AZMAX1(ZMG1P2(N6,N5,N4))*trcs_VLN(ids_H1PO4,N6,N5,N4)
    RFLH0B=VFLW*AZMAX1(H0POB2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLH3B=VFLW*AZMAX1(H3POB2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLF1B=VFLW*AZMAX1(ZF1PB2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLF2B=VFLW*AZMAX1(ZF2PB2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC0B=VFLW*AZMAX1(ZC0PB2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC1B=VFLW*AZMAX1(ZC1PB2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLC2B=VFLW*AZMAX1(ZC2PB2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    RFLM1B=VFLW*AZMAX1(ZM1PB2(N6,N5,N4))*trcs_VLN(ids_H1PO4B,N6,N5,N4)
!
!     NO MACROPORE TO MICROPORE TRANSFER
!
  ELSE
    RFLAL=0.0_r8
    RFLFE=0.0_r8
    RFLHY=0.0_r8
    RFLCA=0.0_r8
    RFLMG=0.0_r8
    RFLNA=0.0_r8
    RFLKA=0.0_r8
    RFLOH=0.0_r8
    RFLSO=0.0_r8
    RFLCL=0.0_r8
    RFLC3=0.0_r8
    RFLHC=0.0_r8
    RFLAL1=0.0_r8
    RFLAL2=0.0_r8
    RFLAL3=0.0_r8
    RFLAL4=0.0_r8
    RFLALS=0.0_r8
    RFLFE1=0.0_r8
    RFLFE2=0.0_r8
    RFLFE3=0.0_r8
    RFLFE4=0.0_r8
    RFLFES=0.0_r8
    RFLCAO=0.0_r8
    RFLCAC=0.0_r8
    RFLCAH=0.0_r8
    RFLCAS=0.0_r8
    RFLMGO=0.0_r8
    RFLMGC=0.0_r8
    RFLMGH=0.0_r8
    RFLMGS=0.0_r8
    RFLNAC=0.0_r8
    RFLNAS=0.0_r8
    RFLKAS=0.0_r8
    RFLH0P=0.0_r8
    RFLH3P=0.0_r8
    RFLF1P=0.0_r8
    RFLF2P=0.0_r8
    RFLC0P=0.0_r8
    RFLC1P=0.0_r8
    RFLC2P=0.0_r8
    RFLM1P=0.0_r8
    RFLH0B=0.0_r8
    RFLH3B=0.0_r8
    RFLF1B=0.0_r8
    RFLF2B=0.0_r8
    RFLC0B=0.0_r8
    RFLC1B=0.0_r8
    RFLC2B=0.0_r8
    RFLM1B=0.0_r8
  ENDIF
  RALFXS(N6,N5,N4)=RFLAL
  RFEFXS(N6,N5,N4)=RFLFE
  RHYFXS(N6,N5,N4)=RFLHY
  RCAFXS(N6,N5,N4)=RFLCA
  RMGFXS(N6,N5,N4)=RFLMG
  RNAFXS(N6,N5,N4)=RFLNA
  RKAFXS(N6,N5,N4)=RFLKA
  ROHFXS(N6,N5,N4)=RFLOH
  RSOFXS(N6,N5,N4)=RFLSO
  RCLFXS(N6,N5,N4)=RFLCL
  RC3FXS(N6,N5,N4)=RFLC3
  RHCFXS(N6,N5,N4)=RFLHC
  RAL1XS(N6,N5,N4)=RFLAL1
  RAL2XS(N6,N5,N4)=RFLAL2
  RAL3XS(N6,N5,N4)=RFLAL3
  RAL4XS(N6,N5,N4)=RFLAL4
  RALSXS(N6,N5,N4)=RFLALS
  RFE1XS(N6,N5,N4)=RFLFE1
  RFE2XS(N6,N5,N4)=RFLFE2
  RFE3XS(N6,N5,N4)=RFLFE3
  RFE4XS(N6,N5,N4)=RFLFE4
  RFESXS(N6,N5,N4)=RFLFES
  RCAOXS(N6,N5,N4)=RFLCAO
  RCACXS(N6,N5,N4)=RFLCAC
  RCAHXS(N6,N5,N4)=RFLCAH
  RCASXS(N6,N5,N4)=RFLCAS
  RMGOXS(N6,N5,N4)=RFLMGO
  RMGCXS(N6,N5,N4)=RFLMGC
  RMGHXS(N6,N5,N4)=RFLMGH
  RMGSXS(N6,N5,N4)=RFLMGS
  RNACXS(N6,N5,N4)=RFLNAC
  RNASXS(N6,N5,N4)=RFLNAS
  RKASXS(N6,N5,N4)=RFLKAS
  RH0PXS(N6,N5,N4)=RFLH0P
  RH3PXS(N6,N5,N4)=RFLH3P
  RF1PXS(N6,N5,N4)=RFLF1P
  RF2PXS(N6,N5,N4)=RFLF2P
  RC0PXS(N6,N5,N4)=RFLC0P
  RC1PXS(N6,N5,N4)=RFLC1P
  RC2PXS(N6,N5,N4)=RFLC2P
  RM1PXS(N6,N5,N4)=RFLM1P
  RH0BXB(N6,N5,N4)=RFLH0B
  RH3BXB(N6,N5,N4)=RFLH3B
  RF1BXB(N6,N5,N4)=RFLF1B
  RF2BXB(N6,N5,N4)=RFLF2B
  RC0BXB(N6,N5,N4)=RFLC0B
  RC1BXB(N6,N5,N4)=RFLC1B
  RC2BXB(N6,N5,N4)=RFLC2B
  RM1BXB(N6,N5,N4)=RFLM1B
  end subroutine SoluteAdvExchMicMacpores
!----------------------------------------------------------------------
  subroutine SoluteDifsExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6
  real(r8) :: VOLWHS,VOLWT

!     DIFFUSIVE FLUXES OF SOLUTES BETWEEN MICROPORES AND
!     MACROPORES FROM AQUEOUS DIFFUSIVITIES AND CONCENTRATION DIFFERENCES
!
!     VOLWM,VOLWHM=micropore,macropore water-filled porosity from watsub.f
!     DFV*S,DFV*B=diffusive solute flux between macro- and micropore in non-band,band
!     XNPH=1/no. of cycles h-1 for water, heat and solute flux calculations
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     *2,*H2=solute content of micropores,macropores
!
  IF(VOLWHM(M,N6,N5,N4).GT.ZEROS2(NY,NX))THEN
    VOLWHS=AMIN1(XFRS*VOLT(N6,N5,N4),VOLWHM(M,N6,N5,N4))
    VOLWT=VOLWM(M,N6,N5,N4)+VOLWHS
    DFVAL=XNPH*(AZMAX1(ZALH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZAL2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVFE=XNPH*(AZMAX1(ZFEH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZFE2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVHY=XNPH*(AZMAX1(ZHYH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZHY2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVCA=XNPH*(AZMAX1(ZCCH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZCA2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVMG=XNPH*(AZMAX1(ZMAH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZMG2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVNA=XNPH*(AZMAX1(ZNAH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZNA2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVKA=XNPH*(AZMAX1(ZKAH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZKA2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVOH=XNPH*(AZMAX1(ZOHH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZOH2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVSO=XNPH*(AZMAX1(ZSO4H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZSO42(N6,N5,N4))*VOLWHS)/VOLWT
    DFVCL=XNPH*(AZMAX1(ZCLH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZCL2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVC3=XNPH*(AZMAX1(ZCO3H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZCO32(N6,N5,N4))*VOLWHS)/VOLWT
    DFVHC=XNPH*(AZMAX1(ZHCOH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZHCO32(N6,N5,N4))*VOLWHS)/VOLWT
    DFVAL1=XNPH*(AZMAX1(ZAL1H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZAL12(N6,N5,N4))*VOLWHS)/VOLWT
    DFVAL2=XNPH*(AZMAX1(ZAL2H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZAL22(N6,N5,N4))*VOLWHS)/VOLWT
    DFVAL3=XNPH*(AZMAX1(ZAL3H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZAL32(N6,N5,N4))*VOLWHS)/VOLWT
    DFVAL4=XNPH*(AZMAX1(ZAL4H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZAL42(N6,N5,N4))*VOLWHS)/VOLWT
    DFVALS=XNPH*(AZMAX1(ZALSH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZALS2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVFE1=XNPH*(AZMAX1(ZFE1H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZFE12(N6,N5,N4))*VOLWHS)/VOLWT
    DFVF22=XNPH*(AZMAX1(ZFE2H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZFE22(N6,N5,N4))*VOLWHS)/VOLWT
    DFVFE3=XNPH*(AZMAX1(ZFE3H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZFE32(N6,N5,N4))*VOLWHS)/VOLWT
    DFVFE4=XNPH*(AZMAX1(ZFE4H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZFE42(N6,N5,N4))*VOLWHS)/VOLWT
    DFVFES=XNPH*(AZMAX1(ZFESH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZFES2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVCAO=XNPH*(AZMAX1(ZCAOH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZCAO2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVCAC=XNPH*(AZMAX1(ZCACH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZCAC2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVCAH=XNPH*(AZMAX1(ZCAHH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZCAH2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVCAS=XNPH*(AZMAX1(ZCASH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZCAS2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVMGO=XNPH*(AZMAX1(ZMGOH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZMGO2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVMGC=XNPH*(AZMAX1(ZMGCH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZMGC2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVMGH=XNPH*(AZMAX1(ZMGHH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZMGH2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVMGS=XNPH*(AZMAX1(ZMGSH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZMGS2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVNAC=XNPH*(AZMAX1(ZNACH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZNAC2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVNAS=XNPH*(AZMAX1(ZNASH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZNAS2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVKAS=XNPH*(AZMAX1(ZKASH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZKAS2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVNAC=XNPH*(AZMAX1(ZNACH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZNAC2(N6,N5,N4))*VOLWHS)/VOLWT
    DFVH0P=XNPH*(AZMAX1(H0P4H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(H0PO42(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVH3P=XNPH*(AZMAX1(H3P4H2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(H3PO42(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVF1P=XNPH*(AZMAX1(ZF1PH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZFE1P2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVF2P=XNPH*(AZMAX1(ZF2PH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZFE2P2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC0P=XNPH*(AZMAX1(ZC0PH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZCA0P2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC1P=XNPH*(AZMAX1(ZC1PH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZCA1P2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVC2P=XNPH*(AZMAX1(ZC2PH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZCA2P2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVM1P=XNPH*(AZMAX1(ZM1PH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZMG1P2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4,N6,N5,N4)
    DFVH0B=XNPH*(AZMAX1(H0PBH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(H0POB2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVH3B=XNPH*(AZMAX1(H3PBH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(H3POB2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVF1B=XNPH*(AZMAX1(ZF1BH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZF1PB2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVF2B=XNPH*(AZMAX1(ZF2BH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZF2PB2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC0B=XNPH*(AZMAX1(ZC0BH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZC0PB2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC1B=XNPH*(AZMAX1(ZC1BH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZC1PB2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVC2B=XNPH*(AZMAX1(ZC2BH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZC2PB2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
    DFVM1B=XNPH*(AZMAX1(ZM1BH2(N6,N5,N4))*VOLWM(M,N6,N5,N4)-AZMAX1(ZM1PB2(N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN(ids_H1PO4B,N6,N5,N4)
  ELSE
    DFVAL=0.0_r8
    DFVFE=0.0_r8
    DFVHY=0.0_r8
    DFVCA=0.0_r8
    DFVMG=0.0_r8
    DFVNA=0.0_r8
    DFVKA=0.0_r8
    DFVOH=0.0_r8
    DFVSO=0.0_r8
    DFVCL=0.0_r8
    DFVC3=0.0_r8
    DFVHC=0.0_r8
    DFVAL1=0.0_r8
    DFVAL2=0.0_r8
    DFVAL3=0.0_r8
    DFVAL4=0.0_r8
    DFVALS=0.0_r8
    DFVFE1=0.0_r8
    DFVFE2=0.0_r8
    DFVFE3=0.0_r8
    DFVFE4=0.0_r8
    DFVFES=0.0_r8
    DFVCAO=0.0_r8
    DFVCAC=0.0_r8
    DFVCAH=0.0_r8
    DFVCAS=0.0_r8
    DFVMGO=0.0_r8
    DFVMGC=0.0_r8
    DFVMGH=0.0_r8
    DFVMGS=0.0_r8
    DFVNAC=0.0_r8
    DFVNAS=0.0_r8
    DFVKAS=0.0_r8
    DFVH0P=0.0_r8
    DFVH3P=0.0_r8
    DFVF1P=0.0_r8
    DFVF2P=0.0_r8
    DFVC0P=0.0_r8
    DFVC1P=0.0_r8
    DFVC2P=0.0_r8
    DFVM1P=0.0_r8
    DFVH0B=0.0_r8
    DFVH3B=0.0_r8
    DFVF1B=0.0_r8
    DFVF2B=0.0_r8
    DFVC0B=0.0_r8
    DFVC1B=0.0_r8
    DFVC2B=0.0_r8
    DFVM1B=0.0_r8
  ENDIF
  RALFXS(N6,N5,N4)=RALFXS(N6,N5,N4)+DFVAL
  RFEFXS(N6,N5,N4)=RFEFXS(N6,N5,N4)+DFVFE
  RHYFXS(N6,N5,N4)=RHYFXS(N6,N5,N4)+DFVHY
  RCAFXS(N6,N5,N4)=RCAFXS(N6,N5,N4)+DFVCA
  RMGFXS(N6,N5,N4)=RMGFXS(N6,N5,N4)+DFVMG
  RNAFXS(N6,N5,N4)=RNAFXS(N6,N5,N4)+DFVNA
  RKAFXS(N6,N5,N4)=RKAFXS(N6,N5,N4)+DFVKA
  ROHFXS(N6,N5,N4)=ROHFXS(N6,N5,N4)+DFVOH
  RSOFXS(N6,N5,N4)=RSOFXS(N6,N5,N4)+DFVSO
  RCLFXS(N6,N5,N4)=RCLFXS(N6,N5,N4)+DFVCL
  RC3FXS(N6,N5,N4)=RC3FXS(N6,N5,N4)+DFVC3
  RHCFXS(N6,N5,N4)=RHCFXS(N6,N5,N4)+DFVHC
  RAL1XS(N6,N5,N4)=RAL1XS(N6,N5,N4)+DFVAL1
  RAL2XS(N6,N5,N4)=RAL2XS(N6,N5,N4)+DFVAL2
  RAL3XS(N6,N5,N4)=RAL3XS(N6,N5,N4)+DFVAL3
  RAL4XS(N6,N5,N4)=RAL4XS(N6,N5,N4)+DFVAL4
  RALSXS(N6,N5,N4)=RALSXS(N6,N5,N4)+DFVALS
  RFE1XS(N6,N5,N4)=RFE1XS(N6,N5,N4)+DFVFE1
  RFE2XS(N6,N5,N4)=RFE2XS(N6,N5,N4)+DFVFE2
  RFE3XS(N6,N5,N4)=RFE3XS(N6,N5,N4)+DFVFE3
  RFE4XS(N6,N5,N4)=RFE4XS(N6,N5,N4)+DFVFE4
  RFESXS(N6,N5,N4)=RFESXS(N6,N5,N4)+DFVFES
  RCAOXS(N6,N5,N4)=RCAOXS(N6,N5,N4)+DFVCAO
  RCACXS(N6,N5,N4)=RCACXS(N6,N5,N4)+DFVCAC
  RCAHXS(N6,N5,N4)=RCAHXS(N6,N5,N4)+DFVCAH
  RCASXS(N6,N5,N4)=RCASXS(N6,N5,N4)+DFVCAS
  RMGOXS(N6,N5,N4)=RMGOXS(N6,N5,N4)+DFVMGO
  RMGCXS(N6,N5,N4)=RMGCXS(N6,N5,N4)+DFVMGC
  RMGHXS(N6,N5,N4)=RMGHXS(N6,N5,N4)+DFVMGH
  RMGSXS(N6,N5,N4)=RMGSXS(N6,N5,N4)+DFVMGS
  RNACXS(N6,N5,N4)=RNACXS(N6,N5,N4)+DFVNAC
  RNASXS(N6,N5,N4)=RNASXS(N6,N5,N4)+DFVNAS
  RKASXS(N6,N5,N4)=RKASXS(N6,N5,N4)+DFVKAS
  RH0PXS(N6,N5,N4)=RH0PXS(N6,N5,N4)+DFVH0P
  RH3PXS(N6,N5,N4)=RH3PXS(N6,N5,N4)+DFVH3P
  RF1PXS(N6,N5,N4)=RF1PXS(N6,N5,N4)+DFVF1P
  RF2PXS(N6,N5,N4)=RF2PXS(N6,N5,N4)+DFVF2P
  RC0PXS(N6,N5,N4)=RC0PXS(N6,N5,N4)+DFVC0P
  RC1PXS(N6,N5,N4)=RC1PXS(N6,N5,N4)+DFVC1P
  RC2PXS(N6,N5,N4)=RC2PXS(N6,N5,N4)+DFVC2P
  RM1PXS(N6,N5,N4)=RM1PXS(N6,N5,N4)+DFVM1P
  RH0BXB(N6,N5,N4)=RH0BXB(N6,N5,N4)+DFVH0B
  RH3BXB(N6,N5,N4)=RH3BXB(N6,N5,N4)+DFVH3B
  RF1BXB(N6,N5,N4)=RF1BXB(N6,N5,N4)+DFVF1B
  RF2BXB(N6,N5,N4)=RF2BXB(N6,N5,N4)+DFVF2B
  RC0BXB(N6,N5,N4)=RC0BXB(N6,N5,N4)+DFVC0B
  RC1BXB(N6,N5,N4)=RC1BXB(N6,N5,N4)+DFVC1B
  RC2BXB(N6,N5,N4)=RC2BXB(N6,N5,N4)+DFVC2B
  RM1BXB(N6,N5,N4)=RM1BXB(N6,N5,N4)+DFVM1B
  end subroutine SoluteDifsExchMicMacpores
!----------------------------------------------------------------------
  subroutine SoluteAdvDifsMicMacpore(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)

!     SOLUTE TRANSPORT IN MICROPORES

  call SoluteAdvMicropore(M,N,N1,N2,N3,N4,N5,N6)

!
  call SoluteDifsMicropore(M,N,N1,N2,N3,N4,N5,N6,THETW1)

!
!     SOLUTE TRANSPORT IN MACROPORES
  call SoluteAdvMacropore(M,N,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsMacropore(M,N,N1,N2,N3,N4,N5,N6)
!     TOTAL MICROPORE AND MACROPORE SOLUTE TRANSPORT
!     FLUXES = CONVECTIVE + DIFFUSIVE FLUXES
!
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,R*FLB=convective + diffusive solute flux through micropores in non-band,band
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,R*FHB=convective + diffusive solute flux through macropores in non-band,band
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FLS=hourly convective + diffusive solute flux through micropores
!     X*FLW,X*FLB= hourly convective + diffusive solute flux through micropores in non-band,band
!     X*FHS=hourly convective + diffusive solute flux through macropores
!     X*FHW,X*FHB= hourly convective + diffusive solute flux through macropores in non-band,band
!     R*FLS=convective + diffusive solute flux through micropores
!     R*FLW,X*FLB=convective + diffusive solute flux through micropores in non-band,band
!     R*FHS=convective + diffusive solute flux through macropores
!     R*FHW,X*FHB=convective + diffusive solute flux through macropores in non-band,band
!
  trcsa_XFLS(idsa_Al,N,N6,N5,N4)=trcsa_XFLS(idsa_Al,N,N6,N5,N4)+RALFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_Fe,N,N6,N5,N4)=trcsa_XFLS(idsa_Fe,N,N6,N5,N4)+RFEFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_Hp,N,N6,N5,N4)=trcsa_XFLS(idsa_Hp,N,N6,N5,N4)+RHYFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_Ca,N,N6,N5,N4)=trcsa_XFLS(idsa_Ca,N,N6,N5,N4)+RCAFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_Mg,N,N6,N5,N4)=trcsa_XFLS(idsa_Mg,N,N6,N5,N4)+RMGFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_Na,N,N6,N5,N4)=trcsa_XFLS(idsa_Na,N,N6,N5,N4)+RNAFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_K,N,N6,N5,N4)=trcsa_XFLS(idsa_K,N,N6,N5,N4)+RKAFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_OH,N,N6,N5,N4)=trcsa_XFLS(idsa_OH,N,N6,N5,N4)+ROHFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_SO4,N,N6,N5,N4)=trcsa_XFLS(idsa_SO4,N,N6,N5,N4)+RSOFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_Cl,N,N6,N5,N4)=trcsa_XFLS(idsa_Cl,N,N6,N5,N4)+RCLFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_CO3,N,N6,N5,N4)=trcsa_XFLS(idsa_CO3,N,N6,N5,N4)+RC3FLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_HCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_HCO3,N,N6,N5,N4)+RHCFLS(N,N6,N5,N4)
  trcsa_XFLS(idsa_AlOH,N,N6,N5,N4)=trcsa_XFLS(idsa_AlOH,N,N6,N5,N4)+RAL1FS(N,N6,N5,N4)
  trcsa_XFLS(idsa_AlOH2,N,N6,N5,N4)=trcsa_XFLS(idsa_AlOH2,N,N6,N5,N4)+RAL2FS(N,N6,N5,N4)
  trcsa_XFLS(idsa_AlOH3,N,N6,N5,N4)=trcsa_XFLS(idsa_AlOH3,N,N6,N5,N4)+RAL3FS(N,N6,N5,N4)
  trcsa_XFLS(idsa_AlOH4,N,N6,N5,N4)=trcsa_XFLS(idsa_AlOH4,N,N6,N5,N4)+RAL4FS(N,N6,N5,N4)
  trcsa_XFLS(idsa_AlSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_AlSO4,N,N6,N5,N4)+RALSFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_FeOH,N,N6,N5,N4)=trcsa_XFLS(idsa_FeOH,N,N6,N5,N4)+RFE1FS(N,N6,N5,N4)
  trcsa_XFLS(idsa_FeOH2,N,N6,N5,N4)=trcsa_XFLS(idsa_FeOH2,N,N6,N5,N4)+RFE2FS(N,N6,N5,N4)
  trcsa_XFLS(idsa_FeOH3,N,N6,N5,N4)=trcsa_XFLS(idsa_FeOH3,N,N6,N5,N4)+RFE3FS(N,N6,N5,N4)
  trcsa_XFLS(idsa_FeOH4,N,N6,N5,N4)=trcsa_XFLS(idsa_FeOH4,N,N6,N5,N4)+RFE4FS(N,N6,N5,N4)
  trcsa_XFLS(idsa_FeSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_FeSO4,N,N6,N5,N4)+RFESFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_CaOH2,N,N6,N5,N4)=trcsa_XFLS(idsa_CaOH2,N,N6,N5,N4)+RCAOFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_CaCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_CaCO3,N,N6,N5,N4)+RCACFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_CaHCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_CaHCO3,N,N6,N5,N4)+RCAHFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_CaSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_CaSO4,N,N6,N5,N4)+RCASFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_MgOH2,N,N6,N5,N4)=trcsa_XFLS(idsa_MgOH2,N,N6,N5,N4)+RMGOFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_MgCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_MgCO3,N,N6,N5,N4)+RMGCFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_MgHCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_MgHCO3,N,N6,N5,N4)+RMGHFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_MgSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_MgSO4,N,N6,N5,N4)+RMGSFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_NaCO3,N,N6,N5,N4)=trcsa_XFLS(idsa_NaCO3,N,N6,N5,N4)+RNACFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_NaSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_NaSO4,N,N6,N5,N4)+RNASFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_KSO4,N,N6,N5,N4)=trcsa_XFLS(idsa_KSO4,N,N6,N5,N4)+RKASFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_H0PO4,N,N6,N5,N4)=trcsa_XFLS(idsa_H0PO4,N,N6,N5,N4)+RH0PFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_H3PO4,N,N6,N5,N4)=trcsa_XFLS(idsa_H3PO4,N,N6,N5,N4)+RH3PFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_FeHPO4,N,N6,N5,N4)=trcsa_XFLS(idsa_FeHPO4,N,N6,N5,N4)+RF1PFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_FeH2PO4,N,N6,N5,N4)=trcsa_XFLS(idsa_FeH2PO4,N,N6,N5,N4)+RF2PFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_CaPO4,N,N6,N5,N4)=trcsa_XFLS(idsa_CaPO4,N,N6,N5,N4)+RC0PFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_CaHPO4,N,N6,N5,N4)=trcsa_XFLS(idsa_CaHPO4,N,N6,N5,N4)+RC1PFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_CaH2PO4,N,N6,N5,N4)=trcsa_XFLS(idsa_CaH2PO4,N,N6,N5,N4)+RC2PFS(N,N6,N5,N4)
  trcsa_XFLS(idsa_MgHPO4,N,N6,N5,N4)=trcsa_XFLS(idsa_MgHPO4,N,N6,N5,N4)+RM1PFS(N,N6,N5,N4)

  XH0BFB(N,N6,N5,N4)=XH0BFB(N,N6,N5,N4)+RH0BFB(N,N6,N5,N4)
  XH3BFB(N,N6,N5,N4)=XH3BFB(N,N6,N5,N4)+RH3BFB(N,N6,N5,N4)
  XF1BFB(N,N6,N5,N4)=XF1BFB(N,N6,N5,N4)+RF1BFB(N,N6,N5,N4)
  XF2BFB(N,N6,N5,N4)=XF2BFB(N,N6,N5,N4)+RF2BFB(N,N6,N5,N4)
  XC0BFB(N,N6,N5,N4)=XC0BFB(N,N6,N5,N4)+RC0BFB(N,N6,N5,N4)
  XC1BFB(N,N6,N5,N4)=XC1BFB(N,N6,N5,N4)+RC1BFB(N,N6,N5,N4)
  XC2BFB(N,N6,N5,N4)=XC2BFB(N,N6,N5,N4)+RC2BFB(N,N6,N5,N4)
  XM1BFB(N,N6,N5,N4)=XM1BFB(N,N6,N5,N4)+RM1BFB(N,N6,N5,N4)

  trcsa_XFHS(idsa_Al,N,N6,N5,N4)=trcsa_XFHS(idsa_Al,N,N6,N5,N4)+RALFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Fe,N,N6,N5,N4)=trcsa_XFHS(idsa_Fe,N,N6,N5,N4)+RFEFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Hp,N,N6,N5,N4)=trcsa_XFHS(idsa_Hp,N,N6,N5,N4)+RHYFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Ca,N,N6,N5,N4)=trcsa_XFHS(idsa_Ca,N,N6,N5,N4)+RCAFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Mg,N,N6,N5,N4)=trcsa_XFHS(idsa_Mg,N,N6,N5,N4)+RMGFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Na,N,N6,N5,N4)=trcsa_XFHS(idsa_Na,N,N6,N5,N4)+RNAFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_K,N,N6,N5,N4)=trcsa_XFHS(idsa_K,N,N6,N5,N4)+RKAFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_OH,N,N6,N5,N4)=trcsa_XFHS(idsa_OH,N,N6,N5,N4)+ROHFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_SO4,N,N6,N5,N4)=trcsa_XFHS(idsa_SO4,N,N6,N5,N4)+RSOFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_Cl,N,N6,N5,N4)=trcsa_XFHS(idsa_Cl,N,N6,N5,N4)+RCLFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CO3,N,N6,N5,N4)=trcsa_XFHS(idsa_CO3,N,N6,N5,N4)+RC3FHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_HCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_HCO3,N,N6,N5,N4)+RHCFHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_AlOH,N,N6,N5,N4)=trcsa_XFHS(idsa_AlOH,N,N6,N5,N4)+RAL1HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_AlOH2,N,N6,N5,N4)=trcsa_XFHS(idsa_AlOH2,N,N6,N5,N4)+RAL2HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_AlOH3,N,N6,N5,N4)=trcsa_XFHS(idsa_AlOH3,N,N6,N5,N4)+RAL3HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_AlOH4,N,N6,N5,N4)=trcsa_XFHS(idsa_AlOH4,N,N6,N5,N4)+RAL4HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_AlSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_AlSO4,N,N6,N5,N4)+RALSHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeOH,N,N6,N5,N4)=trcsa_XFHS(idsa_FeOH,N,N6,N5,N4)+RFE1HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeOH2,N,N6,N5,N4)=trcsa_XFHS(idsa_FeOH2,N,N6,N5,N4)+RFE2HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeOH3,N,N6,N5,N4)=trcsa_XFHS(idsa_FeOH3,N,N6,N5,N4)+RFE3HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeOH4,N,N6,N5,N4)=trcsa_XFHS(idsa_FeOH4,N,N6,N5,N4)+RFE4HS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_FeSO4,N,N6,N5,N4)+RFESHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaOH2,N,N6,N5,N4)=trcsa_XFHS(idsa_CaOH2,N,N6,N5,N4)+RCAOHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_CaCO3,N,N6,N5,N4)+RCACHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaHCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_CaHCO3,N,N6,N5,N4)+RCAHHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_CaSO4,N,N6,N5,N4)+RCASHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgOH2,N,N6,N5,N4)=trcsa_XFHS(idsa_MgOH2,N,N6,N5,N4)+RMGOHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_MgCO3,N,N6,N5,N4)+RMGCHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgHCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_MgHCO3,N,N6,N5,N4)+RMGHHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_MgSO4,N,N6,N5,N4)+RMGSHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_NaCO3,N,N6,N5,N4)=trcsa_XFHS(idsa_NaCO3,N,N6,N5,N4)+RNACHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_NaSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_NaSO4,N,N6,N5,N4)+RNASHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_KSO4,N,N6,N5,N4)=trcsa_XFHS(idsa_KSO4,N,N6,N5,N4)+RKASHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_H0PO4,N,N6,N5,N4)=trcsa_XFHS(idsa_H0PO4,N,N6,N5,N4)+RH0PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_H3PO4,N,N6,N5,N4)=trcsa_XFHS(idsa_H3PO4,N,N6,N5,N4)+RH3PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeHPO4,N,N6,N5,N4)=trcsa_XFHS(idsa_FeHPO4,N,N6,N5,N4)+RF1PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_FeH2PO4,N,N6,N5,N4)=trcsa_XFHS(idsa_FeH2PO4,N,N6,N5,N4)+RF2PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaPO4,N,N6,N5,N4)=trcsa_XFHS(idsa_CaPO4,N,N6,N5,N4)+RC0PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaHPO4,N,N6,N5,N4)=trcsa_XFHS(idsa_CaHPO4,N,N6,N5,N4)+RC1PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_CaH2PO4,N,N6,N5,N4)=trcsa_XFHS(idsa_CaH2PO4,N,N6,N5,N4)+RC2PHS(N,N6,N5,N4)
  trcsa_XFHS(idsa_MgHPO4,N,N6,N5,N4)=trcsa_XFHS(idsa_MgHPO4,N,N6,N5,N4)+RM1PHS(N,N6,N5,N4)

  XH0BHB(N,N6,N5,N4)=XH0BHB(N,N6,N5,N4)+RH0BHB(N,N6,N5,N4)
  XH3BHB(N,N6,N5,N4)=XH3BHB(N,N6,N5,N4)+RH3BHB(N,N6,N5,N4)
  XF1BHB(N,N6,N5,N4)=XF1BHB(N,N6,N5,N4)+RF1BHB(N,N6,N5,N4)
  XF2BHB(N,N6,N5,N4)=XF2BHB(N,N6,N5,N4)+RF2BHB(N,N6,N5,N4)
  XC0BHB(N,N6,N5,N4)=XC0BHB(N,N6,N5,N4)+RC0BHB(N,N6,N5,N4)
  XC1BHB(N,N6,N5,N4)=XC1BHB(N,N6,N5,N4)+RC1BHB(N,N6,N5,N4)
  XC2BHB(N,N6,N5,N4)=XC2BHB(N,N6,N5,N4)+RC2BHB(N,N6,N5,N4)
  XM1BHB(N,N6,N5,N4)=XM1BHB(N,N6,N5,N4)+RM1BHB(N,N6,N5,N4)
  end subroutine SoluteAdvDifsMicMacpore
!----------------------------------------------------------------------
  subroutine SoluteAdvDifsExchMicMacpore(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,NY,NX,N1,N2,N3,N4,N5,N6


  call SoluteAdvExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
!
  call SoluteDifsExchMicMacpores(M,N,NY,NX,N1,N2,N3,N4,N5,N6)


!
!     TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
!
!     R*FXS,R*FXB=total convective + diffusive solute flux between macro- and micropore in non-band,band
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!     RFL*=convective flux between macro- and micropore
!     DFV*=diffusive solute flux between macro- and micropore
!

!
!     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
!
!     X*FXS,X*FXB= hourly convective + diffusive solute flux between macro- and micropore in non-band,band
!     R*FXS,R*FXB=convective + diffusive solute flux between macro- and micropore in non-band,band
!
  XALFXS(N6,N5,N4)=XALFXS(N6,N5,N4)+RALFXS(N6,N5,N4)
  XFEFXS(N6,N5,N4)=XFEFXS(N6,N5,N4)+RFEFXS(N6,N5,N4)
  XHYFXS(N6,N5,N4)=XHYFXS(N6,N5,N4)+RHYFXS(N6,N5,N4)
  XCAFXS(N6,N5,N4)=XCAFXS(N6,N5,N4)+RCAFXS(N6,N5,N4)
  XMGFXS(N6,N5,N4)=XMGFXS(N6,N5,N4)+RMGFXS(N6,N5,N4)
  XNAFXS(N6,N5,N4)=XNAFXS(N6,N5,N4)+RNAFXS(N6,N5,N4)
  XKAFXS(N6,N5,N4)=XKAFXS(N6,N5,N4)+RKAFXS(N6,N5,N4)
  XOHFXS(N6,N5,N4)=XOHFXS(N6,N5,N4)+ROHFXS(N6,N5,N4)
  XSOFXS(N6,N5,N4)=XSOFXS(N6,N5,N4)+RSOFXS(N6,N5,N4)
  XCLFXS(N6,N5,N4)=XCLFXS(N6,N5,N4)+RCLFXS(N6,N5,N4)
  XC3FXS(N6,N5,N4)=XC3FXS(N6,N5,N4)+RC3FXS(N6,N5,N4)
  XHCFXS(N6,N5,N4)=XHCFXS(N6,N5,N4)+RHCFXS(N6,N5,N4)
  XAL1XS(N6,N5,N4)=XAL1XS(N6,N5,N4)+RAL1XS(N6,N5,N4)
  XAL2XS(N6,N5,N4)=XAL2XS(N6,N5,N4)+RAL2XS(N6,N5,N4)
  XAL3XS(N6,N5,N4)=XAL3XS(N6,N5,N4)+RAL3XS(N6,N5,N4)
  XAL4XS(N6,N5,N4)=XAL4XS(N6,N5,N4)+RAL4XS(N6,N5,N4)
  XALSXS(N6,N5,N4)=XALSXS(N6,N5,N4)+RALSXS(N6,N5,N4)
  XFE1XS(N6,N5,N4)=XFE1XS(N6,N5,N4)+RFE1XS(N6,N5,N4)
  XFE2XS(N6,N5,N4)=XFE2XS(N6,N5,N4)+RFE2XS(N6,N5,N4)
  XFE3XS(N6,N5,N4)=XFE3XS(N6,N5,N4)+RFE3XS(N6,N5,N4)
  XFE4XS(N6,N5,N4)=XFE4XS(N6,N5,N4)+RFE4XS(N6,N5,N4)
  XFESXS(N6,N5,N4)=XFESXS(N6,N5,N4)+RFESXS(N6,N5,N4)
  XCAOXS(N6,N5,N4)=XCAOXS(N6,N5,N4)+RCAOXS(N6,N5,N4)
  XCACXS(N6,N5,N4)=XCACXS(N6,N5,N4)+RCACXS(N6,N5,N4)
  XCAHXS(N6,N5,N4)=XCAHXS(N6,N5,N4)+RCAHXS(N6,N5,N4)
  XCASXS(N6,N5,N4)=XCASXS(N6,N5,N4)+RCASXS(N6,N5,N4)
  XMGOXS(N6,N5,N4)=XMGOXS(N6,N5,N4)+RMGOXS(N6,N5,N4)
  XMGCXS(N6,N5,N4)=XMGCXS(N6,N5,N4)+RMGCXS(N6,N5,N4)
  XMGHXS(N6,N5,N4)=XMGHXS(N6,N5,N4)+RMGHXS(N6,N5,N4)
  XMGSXS(N6,N5,N4)=XMGSXS(N6,N5,N4)+RMGSXS(N6,N5,N4)
  XNACXS(N6,N5,N4)=XNACXS(N6,N5,N4)+RNACXS(N6,N5,N4)
  XNASXS(N6,N5,N4)=XNASXS(N6,N5,N4)+RNASXS(N6,N5,N4)
  XKASXS(N6,N5,N4)=XKASXS(N6,N5,N4)+RKASXS(N6,N5,N4)
  XH0PXS(N6,N5,N4)=XH0PXS(N6,N5,N4)+RH0PXS(N6,N5,N4)
  XH3PXS(N6,N5,N4)=XH3PXS(N6,N5,N4)+RH3PXS(N6,N5,N4)
  XF1PXS(N6,N5,N4)=XF1PXS(N6,N5,N4)+RF1PXS(N6,N5,N4)
  XF2PXS(N6,N5,N4)=XF2PXS(N6,N5,N4)+RF2PXS(N6,N5,N4)
  XC0PXS(N6,N5,N4)=XC0PXS(N6,N5,N4)+RC0PXS(N6,N5,N4)
  XC1PXS(N6,N5,N4)=XC1PXS(N6,N5,N4)+RC1PXS(N6,N5,N4)
  XC2PXS(N6,N5,N4)=XC2PXS(N6,N5,N4)+RC2PXS(N6,N5,N4)
  XM1PXS(N6,N5,N4)=XM1PXS(N6,N5,N4)+RM1PXS(N6,N5,N4)
  XH0BXB(N6,N5,N4)=XH0BXB(N6,N5,N4)+RH0BXB(N6,N5,N4)
  XH3BXB(N6,N5,N4)=XH3BXB(N6,N5,N4)+RH3BXB(N6,N5,N4)
  XF1BXB(N6,N5,N4)=XF1BXB(N6,N5,N4)+RF1BXB(N6,N5,N4)
  XF2BXB(N6,N5,N4)=XF2BXB(N6,N5,N4)+RF2BXB(N6,N5,N4)
  XC0BXB(N6,N5,N4)=XC0BXB(N6,N5,N4)+RC0BXB(N6,N5,N4)
  XC1BXB(N6,N5,N4)=XC1BXB(N6,N5,N4)+RC1BXB(N6,N5,N4)
  XC2BXB(N6,N5,N4)=XC2BXB(N6,N5,N4)+RC2BXB(N6,N5,N4)
  XM1BXB(N6,N5,N4)=XM1BXB(N6,N5,N4)+RM1BXB(N6,N5,N4)
  end subroutine SoluteAdvDifsExchMicMacpore
!----------------------------------------------------------------------
  subroutine UpdateSoluteInSubsurfNeighbors(M,NY,NX,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX,NHW,NHE,NVN,NVS

  integer :: IFLGB,L,N1,N2,N3,N4,N5,N6
  integer :: LL,N
  real(r8) :: VOLWHS,VOLWT
  real(r8) :: VFLW,VLWPA1,VLWPB1,VLWPA2,VLWPB2
  real(r8) :: THETW1(JZ,JY,JX)

!     begin_execution
!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  IFLGB=0
  DO  L=1,NL(NY,NX)
    N1=NX
    N2=NY
    N3=L
!
!     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
!
    DO  N=NCN(N2,N1),3
      IF(N.EQ.1)THEN
        IF(NX.EQ.NHE)THEN
          cycle
        ELSE
          N4=NX+1
          N5=NY
          N6=L
        ENDIF
      ELSEIF(N.EQ.2)THEN
        IF(NY.EQ.NVS)THEN
          cycle
        ELSE
          N4=NX
          N5=NY+1
          N6=L
        ENDIF
      ELSEIF(N.EQ.3)THEN
        IF(L.EQ.NL(NY,NX))THEN
          cycle
        ELSE
          N4=NX
          N5=NY
          N6=L+1
        ENDIF
      ENDIF
      DO LL=N6,NL(NY,NX)
        IF(VOLX(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
          N6=LL
          exit
        ENDIF
      enddo
!
!     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS FROM
!     WATER CONTENTS AND WATER FLUXES 'FLQM' FROM 'WATSUB'
!
!     VOLX,VOLY=soil volume excluding rock, macropore
!     VLNH4,VLNO3,VLPO4=non-band NH4,NO3,PO4 volume fraction
!     VLNHB,VLNOB,VLPOB=band NH4,NO3,PO4 volume fraction
!     VOLWM,VOLWHM=micropore,macropore water-filled porosity from watsub.f
!     THETW=volumetric water content
!     FLPM=change in air volume
!     XNPT=1/number of cycles NPH-1 for gas flux calculations
!
      IF(VOLX(N3,N2,N1).GT.ZEROS2(NY,NX))THEN
        IF(N3.GE.NUM(N2,N1).AND.N6.GE.NUM(N5,N4).AND.N3.LE.NL(N2,N1).AND.N6.LE.NL(N5,N4))THEN
          THETW1(N3,N2,N1)=AZMAX1(VOLWM(M,N3,N2,N1)/VOLY(N3,N2,N1))
          THETW1(N6,N5,N4)=AZMAX1(VOLWM(M,N6,N5,N4)/VOLY(N6,N5,N4))
!
          call SoluteAdvDifsMicMacpore(M,N,N1,N2,N3,N4,N5,N6,THETW1)

!
!     MACROPORE-MICROPORE SOLUTE EXCHANGE WITHIN SOIL
!     LAYER FROM WATER EXCHANGE IN 'WATSUB' AND
!     FROM MACROPORE OR MICROPORE SOLUTE CONCENTRATIONS
!
          IF(N.EQ.3)THEN
             call SoluteAdvDifsExchMicMacpore(M,N,NY,NX,N1,N2,N3,N4,N5,N6)
          ENDIF
        ELSEIF(N.NE.3)THEN
          THETW1(N3,N2,N1)=0.0_r8
          THETW1(N6,N5,N4)=0.0_r8
          RALFLS(N,N6,N5,N4)=0.0_r8
          RFEFLS(N,N6,N5,N4)=0.0_r8
          RHYFLS(N,N6,N5,N4)=0.0_r8
          RCAFLS(N,N6,N5,N4)=0.0_r8
          RMGFLS(N,N6,N5,N4)=0.0_r8
          RNAFLS(N,N6,N5,N4)=0.0_r8
          RKAFLS(N,N6,N5,N4)=0.0_r8
          ROHFLS(N,N6,N5,N4)=0.0_r8
          RSOFLS(N,N6,N5,N4)=0.0_r8
          RCLFLS(N,N6,N5,N4)=0.0_r8
          RC3FLS(N,N6,N5,N4)=0.0_r8
          RHCFLS(N,N6,N5,N4)=0.0_r8
          RAL1FS(N,N6,N5,N4)=0.0_r8
          RAL2FS(N,N6,N5,N4)=0.0_r8
          RAL3FS(N,N6,N5,N4)=0.0_r8
          RAL4FS(N,N6,N5,N4)=0.0_r8
          RALSFS(N,N6,N5,N4)=0.0_r8
          RFE1FS(N,N6,N5,N4)=0.0_r8
          RFE2FS(N,N6,N5,N4)=0.0_r8
          RFE3FS(N,N6,N5,N4)=0.0_r8
          RFE4FS(N,N6,N5,N4)=0.0_r8
          RFESFS(N,N6,N5,N4)=0.0_r8
          RCAOFS(N,N6,N5,N4)=0.0_r8
          RCACFS(N,N6,N5,N4)=0.0_r8
          RCAHFS(N,N6,N5,N4)=0.0_r8
          RCASFS(N,N6,N5,N4)=0.0_r8
          RMGOFS(N,N6,N5,N4)=0.0_r8
          RMGCFS(N,N6,N5,N4)=0.0_r8
          RMGHFS(N,N6,N5,N4)=0.0_r8
          RMGSFS(N,N6,N5,N4)=0.0_r8
          RNACFS(N,N6,N5,N4)=0.0_r8
          RNASFS(N,N6,N5,N4)=0.0_r8
          RKASFS(N,N6,N5,N4)=0.0_r8
          RH0PFS(N,N6,N5,N4)=0.0_r8
          RH3PFS(N,N6,N5,N4)=0.0_r8
          RF1PFS(N,N6,N5,N4)=0.0_r8
          RF2PFS(N,N6,N5,N4)=0.0_r8
          RC0PFS(N,N6,N5,N4)=0.0_r8
          RC1PFS(N,N6,N5,N4)=0.0_r8
          RC2PFS(N,N6,N5,N4)=0.0_r8
          RM1PFS(N,N6,N5,N4)=0.0_r8
          RH0BFB(N,N6,N5,N4)=0.0_r8
          RH3BFB(N,N6,N5,N4)=0.0_r8
          RF1BFB(N,N6,N5,N4)=0.0_r8
          RF2BFB(N,N6,N5,N4)=0.0_r8
          RC0BFB(N,N6,N5,N4)=0.0_r8
          RC1BFB(N,N6,N5,N4)=0.0_r8
          RC2BFB(N,N6,N5,N4)=0.0_r8
          RM1BFB(N,N6,N5,N4)=0.0_r8
          RALFHS(N,N6,N5,N4)=0.0_r8
          RFEFHS(N,N6,N5,N4)=0.0_r8
          RHYFHS(N,N6,N5,N4)=0.0_r8
          RCAFHS(N,N6,N5,N4)=0.0_r8
          RMGFHS(N,N6,N5,N4)=0.0_r8
          RNAFHS(N,N6,N5,N4)=0.0_r8
          RKAFHS(N,N6,N5,N4)=0.0_r8
          ROHFHS(N,N6,N5,N4)=0.0_r8
          RSOFHS(N,N6,N5,N4)=0.0_r8
          RCLFHS(N,N6,N5,N4)=0.0_r8
          RC3FHS(N,N6,N5,N4)=0.0_r8
          RHCFHS(N,N6,N5,N4)=0.0_r8
          RAL1HS(N,N6,N5,N4)=0.0_r8
          RAL2HS(N,N6,N5,N4)=0.0_r8
          RAL3HS(N,N6,N5,N4)=0.0_r8
          RAL4HS(N,N6,N5,N4)=0.0_r8
          RALSHS(N,N6,N5,N4)=0.0_r8
          RFE1HS(N,N6,N5,N4)=0.0_r8
          RFE2HS(N,N6,N5,N4)=0.0_r8
          RFE3HS(N,N6,N5,N4)=0.0_r8
          RFE4HS(N,N6,N5,N4)=0.0_r8
          RFESHS(N,N6,N5,N4)=0.0_r8
          RCAOHS(N,N6,N5,N4)=0.0_r8
          RCACHS(N,N6,N5,N4)=0.0_r8
          RCAHHS(N,N6,N5,N4)=0.0_r8
          RCASHS(N,N6,N5,N4)=0.0_r8
          RMGOHS(N,N6,N5,N4)=0.0_r8
          RMGCHS(N,N6,N5,N4)=0.0_r8
          RMGHHS(N,N6,N5,N4)=0.0_r8
          RMGSHS(N,N6,N5,N4)=0.0_r8
          RNACHS(N,N6,N5,N4)=0.0_r8
          RNASHS(N,N6,N5,N4)=0.0_r8
          RKASHS(N,N6,N5,N4)=0.0_r8
          RH0PHS(N,N6,N5,N4)=0.0_r8
          RH3PHS(N,N6,N5,N4)=0.0_r8
          RF1PHS(N,N6,N5,N4)=0.0_r8
          RF2PHS(N,N6,N5,N4)=0.0_r8
          RC0PHS(N,N6,N5,N4)=0.0_r8
          RC1PHS(N,N6,N5,N4)=0.0_r8
          RC2PHS(N,N6,N5,N4)=0.0_r8
          RM1PHS(N,N6,N5,N4)=0.0_r8
          RH0BHB(N,N6,N5,N4)=0.0_r8
          RH3BHB(N,N6,N5,N4)=0.0_r8
          RF1BHB(N,N6,N5,N4)=0.0_r8
          RF2BHB(N,N6,N5,N4)=0.0_r8
          RC0BHB(N,N6,N5,N4)=0.0_r8
          RC1BHB(N,N6,N5,N4)=0.0_r8
          RC2BHB(N,N6,N5,N4)=0.0_r8
          RM1BHB(N,N6,N5,N4)=0.0_r8
        ENDIF
      ELSE
        THETW1(N3,N2,N1)=0.0_r8
        THETW1(N6,N5,N4)=0.0_r8
        RALFLS(N,N6,N5,N4)=0.0_r8
        RFEFLS(N,N6,N5,N4)=0.0_r8
        RHYFLS(N,N6,N5,N4)=0.0_r8
        RCAFLS(N,N6,N5,N4)=0.0_r8
        RMGFLS(N,N6,N5,N4)=0.0_r8
        RNAFLS(N,N6,N5,N4)=0.0_r8
        RKAFLS(N,N6,N5,N4)=0.0_r8
        ROHFLS(N,N6,N5,N4)=0.0_r8
        RSOFLS(N,N6,N5,N4)=0.0_r8
        RCLFLS(N,N6,N5,N4)=0.0_r8
        RC3FLS(N,N6,N5,N4)=0.0_r8
        RHCFLS(N,N6,N5,N4)=0.0_r8
        RAL1FS(N,N6,N5,N4)=0.0_r8
        RAL2FS(N,N6,N5,N4)=0.0_r8
        RAL3FS(N,N6,N5,N4)=0.0_r8
        RAL4FS(N,N6,N5,N4)=0.0_r8
        RALSFS(N,N6,N5,N4)=0.0_r8
        RFE1FS(N,N6,N5,N4)=0.0_r8
        RFE2FS(N,N6,N5,N4)=0.0_r8
        RFE3FS(N,N6,N5,N4)=0.0_r8
        RFE4FS(N,N6,N5,N4)=0.0_r8
        RFESFS(N,N6,N5,N4)=0.0_r8
        RCAOFS(N,N6,N5,N4)=0.0_r8
        RCACFS(N,N6,N5,N4)=0.0_r8
        RCAHFS(N,N6,N5,N4)=0.0_r8
        RCASFS(N,N6,N5,N4)=0.0_r8
        RMGOFS(N,N6,N5,N4)=0.0_r8
        RMGCFS(N,N6,N5,N4)=0.0_r8
        RMGHFS(N,N6,N5,N4)=0.0_r8
        RMGSFS(N,N6,N5,N4)=0.0_r8
        RNACFS(N,N6,N5,N4)=0.0_r8
        RNASFS(N,N6,N5,N4)=0.0_r8
        RKASFS(N,N6,N5,N4)=0.0_r8
        RH0PFS(N,N6,N5,N4)=0.0_r8
        RH3PFS(N,N6,N5,N4)=0.0_r8
        RF1PFS(N,N6,N5,N4)=0.0_r8
        RF2PFS(N,N6,N5,N4)=0.0_r8
        RC0PFS(N,N6,N5,N4)=0.0_r8
        RC1PFS(N,N6,N5,N4)=0.0_r8
        RC2PFS(N,N6,N5,N4)=0.0_r8
        RM1PFS(N,N6,N5,N4)=0.0_r8
        RH0BFB(N,N6,N5,N4)=0.0_r8
        RH3BFB(N,N6,N5,N4)=0.0_r8
        RF1BFB(N,N6,N5,N4)=0.0_r8
        RF2BFB(N,N6,N5,N4)=0.0_r8
        RC0BFB(N,N6,N5,N4)=0.0_r8
        RC1BFB(N,N6,N5,N4)=0.0_r8
        RC2BFB(N,N6,N5,N4)=0.0_r8
        RM1BFB(N,N6,N5,N4)=0.0_r8
        RALFHS(N,N6,N5,N4)=0.0_r8
        RFEFHS(N,N6,N5,N4)=0.0_r8
        RHYFHS(N,N6,N5,N4)=0.0_r8
        RCAFHS(N,N6,N5,N4)=0.0_r8
        RMGFHS(N,N6,N5,N4)=0.0_r8
        RNAFHS(N,N6,N5,N4)=0.0_r8
        RKAFHS(N,N6,N5,N4)=0.0_r8
        ROHFHS(N,N6,N5,N4)=0.0_r8
        RSOFHS(N,N6,N5,N4)=0.0_r8
        RCLFHS(N,N6,N5,N4)=0.0_r8
        RC3FHS(N,N6,N5,N4)=0.0_r8
        RHCFHS(N,N6,N5,N4)=0.0_r8
        RAL1HS(N,N6,N5,N4)=0.0_r8
        RAL2HS(N,N6,N5,N4)=0.0_r8
        RAL3HS(N,N6,N5,N4)=0.0_r8
        RAL4HS(N,N6,N5,N4)=0.0_r8
        RALSHS(N,N6,N5,N4)=0.0_r8
        RFE1HS(N,N6,N5,N4)=0.0_r8
        RFE2HS(N,N6,N5,N4)=0.0_r8
        RFE3HS(N,N6,N5,N4)=0.0_r8
        RFE4HS(N,N6,N5,N4)=0.0_r8
        RFESHS(N,N6,N5,N4)=0.0_r8
        RCAOHS(N,N6,N5,N4)=0.0_r8
        RCACHS(N,N6,N5,N4)=0.0_r8
        RCAHHS(N,N6,N5,N4)=0.0_r8
        RCASHS(N,N6,N5,N4)=0.0_r8
        RMGOHS(N,N6,N5,N4)=0.0_r8
        RMGCHS(N,N6,N5,N4)=0.0_r8
        RMGHHS(N,N6,N5,N4)=0.0_r8
        RMGSHS(N,N6,N5,N4)=0.0_r8
        RNACHS(N,N6,N5,N4)=0.0_r8
        RNASHS(N,N6,N5,N4)=0.0_r8
        RKASHS(N,N6,N5,N4)=0.0_r8
        RH0PHS(N,N6,N5,N4)=0.0_r8
        RH3PHS(N,N6,N5,N4)=0.0_r8
        RF1PHS(N,N6,N5,N4)=0.0_r8
        RF2PHS(N,N6,N5,N4)=0.0_r8
        RC0PHS(N,N6,N5,N4)=0.0_r8
        RC1PHS(N,N6,N5,N4)=0.0_r8
        RC2PHS(N,N6,N5,N4)=0.0_r8
        RM1PHS(N,N6,N5,N4)=0.0_r8
        RH0BHB(N,N6,N5,N4)=0.0_r8
        RH3BHB(N,N6,N5,N4)=0.0_r8
        RF1BHB(N,N6,N5,N4)=0.0_r8
        RF2BHB(N,N6,N5,N4)=0.0_r8
        RC0BHB(N,N6,N5,N4)=0.0_r8
        RC1BHB(N,N6,N5,N4)=0.0_r8
        RC2BHB(N,N6,N5,N4)=0.0_r8
        RM1BHB(N,N6,N5,N4)=0.0_r8
      ENDIF
    enddo
  enddo
  end subroutine UpdateSoluteInSubsurfNeighbors
!------------------------------------------------------------------------------------------

  subroutine MacMicPoreFluxAdvPlusDifus(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
!     begin_execution
!
!     R*FXS=convective + diffusive solute flux between macropores and micropores
!     RFL*=convective flux between macropores and micropores
!     DFV*=diffusive solute flux between macropores and micropores
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH2PO4+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  RALFXS(NU(NY,NX),NY,NX)=RFLAL+DFVAL
  RFEFXS(NU(NY,NX),NY,NX)=RFLFE+DFVFE
  RHYFXS(NU(NY,NX),NY,NX)=RFLHY+DFVHY
  RCAFXS(NU(NY,NX),NY,NX)=RFLCA+DFVCA
  RMGFXS(NU(NY,NX),NY,NX)=RFLMG+DFVMG
  RNAFXS(NU(NY,NX),NY,NX)=RFLNA+DFVNA
  RKAFXS(NU(NY,NX),NY,NX)=RFLKA+DFVKA
  ROHFXS(NU(NY,NX),NY,NX)=RFLOH+DFVOH
  RSOFXS(NU(NY,NX),NY,NX)=RFLSO+DFVSO
  RCLFXS(NU(NY,NX),NY,NX)=RFLCL+DFVCL
  RC3FXS(NU(NY,NX),NY,NX)=RFLC3+DFVC3
  RHCFXS(NU(NY,NX),NY,NX)=RFLHC+DFVHC
  RAL1XS(NU(NY,NX),NY,NX)=RFLAL1+DFVAL1
  RAL2XS(NU(NY,NX),NY,NX)=RFLAL2+DFVAL2
  RAL3XS(NU(NY,NX),NY,NX)=RFLAL3+DFVAL3
  RAL4XS(NU(NY,NX),NY,NX)=RFLAL4+DFVAL4
  RALSXS(NU(NY,NX),NY,NX)=RFLALS+DFVALS
  RFE1XS(NU(NY,NX),NY,NX)=RFLFE1+DFVFE1
  RFE2XS(NU(NY,NX),NY,NX)=RFLFE2+DFVFE2
  RFE3XS(NU(NY,NX),NY,NX)=RFLFE3+DFVFE3
  RFE4XS(NU(NY,NX),NY,NX)=RFLFE4+DFVFE4
  RFESXS(NU(NY,NX),NY,NX)=RFLFES+DFVFES
  RCAOXS(NU(NY,NX),NY,NX)=RFLCAO+DFVCAO
  RCACXS(NU(NY,NX),NY,NX)=RFLCAC+DFVCAC
  RCAHXS(NU(NY,NX),NY,NX)=RFLCAH+DFVCAH
  RCASXS(NU(NY,NX),NY,NX)=RFLCAS+DFVCAS
  RMGOXS(NU(NY,NX),NY,NX)=RFLMGO+DFVMGO
  RMGCXS(NU(NY,NX),NY,NX)=RFLMGC+DFVMGC
  RMGHXS(NU(NY,NX),NY,NX)=RFLMGH+DFVMGH
  RMGSXS(NU(NY,NX),NY,NX)=RFLMGS+DFVMGS
  RNACXS(NU(NY,NX),NY,NX)=RFLNAC+DFVNAC
  RNASXS(NU(NY,NX),NY,NX)=RFLNAS+DFVNAS
  RKASXS(NU(NY,NX),NY,NX)=RFLKAS+DFVKAS
  RH0PXS(NU(NY,NX),NY,NX)=RFLH0P+DFVH0P
  RH3PXS(NU(NY,NX),NY,NX)=RFLH3P+DFVH3P
      RF1PXS(NU(NY,NX),NY,NX)=RFLF1P+DFVF1P
      RF2PXS(NU(NY,NX),NY,NX)=RFLF2P+DFVF2P
      RC0PXS(NU(NY,NX),NY,NX)=RFLC0P+DFVC0P
      RC1PXS(NU(NY,NX),NY,NX)=RFLC1P+DFVC1P
      RC2PXS(NU(NY,NX),NY,NX)=RFLC2P+DFVC2P
      RM1PXS(NU(NY,NX),NY,NX)=RFLM1P+DFVM1P
      RH0BXB(NU(NY,NX),NY,NX)=RFLH0B+DFVH0B
      RH3BXB(NU(NY,NX),NY,NX)=RFLH3B+DFVH3B
      RF1BXB(NU(NY,NX),NY,NX)=RFLF1B+DFVF1B
      RF2BXB(NU(NY,NX),NY,NX)=RFLF2B+DFVF2B
      RC0BXB(NU(NY,NX),NY,NX)=RFLC0B+DFVC0B
      RC1BXB(NU(NY,NX),NY,NX)=RFLC1B+DFVC1B
      RC2BXB(NU(NY,NX),NY,NX)=RFLC2B+DFVC2B
      RM1BXB(NU(NY,NX),NY,NX)=RFLM1B+DFVM1B
      end subroutine MacMicPoreFluxAdvPlusDifus
end module IngridTranspMod
