module SoluteMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb
  use SoluteParMod
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use GridConsts
  use FlagDataType
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use SurfLitterDataType
  use RootDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use PlantDataRateType
  use GridDataType
  use EcoSIMSolverPar
  use ChemEquilibriaMod
  use SaltChemEquilibriaMod
  use SoluteChemDataType
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = __FILE__

  real(r8) :: BKVLNH,BKVLNB
  real(r8) :: BKVLNO,BKVLNZ,BKVLPO,BKVLPB,COMA,CNHUA
  real(r8) :: CCO20
  real(r8) :: DWNH4,DPFLW,DNH4S,DNH3S,DXNH4,DWNO3,DNO3S,DNO2S
  real(r8) :: DWPO4,DZH0P,DZH1P,DZH2P,DZH3P,DZF1P,DZF2P,DZC0P
  real(r8) :: DZC1P,DZC2P,DZM1P,DXOH0,DXOH1,DXOH2,DXH1P,DXH2P
  real(r8) :: DPALP,DPFEP,DPCDP,DPCHP,DPCMP,DUKD,DFNSA,DFNSB
  real(r8) :: FLWD,FVLNH4,FVLNO3,FVLPO4
  real(r8) :: P3,RSN4AA,RSN4BA,RSN3AA,RSN3BA,RSNUAA,RSNUBA,RSNOAA
  real(r8) :: RSNOBA,RSN4BB,RSN3BB,RSNUBB,RSNOBB,RN4X,RN3X,RH1PX
  real(r8) :: RH2PX
  real(r8) :: THETWR
  real(r8) :: VOLWMX,VOLWMP,VOLWNX,XVLNH4
!
  real(r8) :: RSNUA, RSNUB


  public :: solute
  contains

  SUBROUTINE solute(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES ALL SOLUTE TRANSFORMATIONS
!     FROM THERMODYNAMIC EQUILIBRIA
!
  implicit none

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW, NHE, NVN, NVS

! declaration of local variables
  integer :: L,NY,NX,NPI
  real(r8) :: RHP1,RHP2,RN3S,RN4S
!     begin_execution
  NPI=INT(NPH/2)
  DO   NX=NHW,NHE
    DO   NY=NVN,NVS
      DO   L=NU(NY,NX),NL(NY,NX)
        IF(VOLX(L,NY,NX).GT.ZEROS2(NY,NX) &
          .AND.VOLWM(NPH,L,NY,NX).GT.ZEROS2(NY,NX))THEN
!
!     WATER VOLUME IN NON-BAND AND BAND SOIL ZONES
!
!     VOLWM=soil water volume
!     VLNH4,VLNHB=fractions of soil volume in NH4 non-band,band
!     VLNO3,VLNOB=fractions of soil volume in N03 non-band,band
!     VLPO4,VLPOB=fractions of soil volume in H2PO4 non-band,band
!     BKVL=soil mass
!
          VOLWNH=VOLWM(NPH,L,NY,NX)*VLNH4(L,NY,NX)
          VOLWNB=VOLWM(NPH,L,NY,NX)*VLNHB(L,NY,NX)
          VOLWNO=VOLWM(NPH,L,NY,NX)*VLNO3(L,NY,NX)
          VOLWNZ=VOLWM(NPH,L,NY,NX)*VLNOB(L,NY,NX)
          VOLWPO=VOLWM(NPH,L,NY,NX)*VLPO4(L,NY,NX)
          VOLWPB=VOLWM(NPH,L,NY,NX)*VLPOB(L,NY,NX)
          IF(BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
            BKVLX=BKVL(L,NY,NX)
            BKVLNH=BKVL(L,NY,NX)*VLNH4(L,NY,NX)
            BKVLNB=BKVL(L,NY,NX)*VLNHB(L,NY,NX)
            BKVLNO=BKVL(L,NY,NX)*VLNO3(L,NY,NX)
            BKVLNZ=BKVL(L,NY,NX)*VLNOB(L,NY,NX)
            BKVLPO=BKVL(L,NY,NX)*VLPO4(L,NY,NX)
            BKVLPB=BKVL(L,NY,NX)*VLPOB(L,NY,NX)
          ELSE
            BKVLX=VOLA(L,NY,NX)
            BKVLNH=VOLWNH
            BKVLNB=VOLWNB
            BKVLNO=VOLWNO
            BKVLNZ=VOLWNZ
            BKVLPO=VOLWPO
            BKVLPB=VOLWPB
          ENDIF

          call UpdateSoilFertlizer(L,NY,NX)
!
!     IF SALT OPTION SELECTED IN SITE FILE
!     THEN SOLVE FULL SET OF EQUILIBRIA REACTIONS
!
!     ISALTG=salt flag from site file
!
          IF(ISALTG.NE.0)THEN

            call SaltChemEquilibria(L,NY,NX)

          ELSE
            call NoSaltChemEquilibria(L,NY,NX)
          ENDIF
!
!     CHANGE IN WIDTHS AND DEPTHS OF FERTILIZER BANDS FROM
!     VERTICAL AND HORIZONTAL DIFFUSION DRIVEN BY CONCENTRATION
!     DIFFERENCES BETWEEN BAND AND NON-BAND SOIL ZONES
!
!     ROWI=band row width
!     FLWD=net vertical flow relative to area
!
!     IF(ROWI(I,NY,NX).GT.0.0)THEN
          FLWD=0.5*(FLW(3,L,NY,NX)+FLW(3,L+1,NY,NX))/AREA(3,L,NY,NX)
!
!     NH4 FERTILIZER BAND
!
          call UpdateNH3FertilizerBandinfo(L,NY,NX)
!
!     NO3 FERTILIZER BAND
!
          call UpdateNO3FertilizerBandinfo(L,NY,NX)
!
!     PO4 FERTILIZER BAND
!
          call UpdatePO4FertilizerBandinfo(L,NY,NX)
!     ENDIF
!
!     SUBTRACT FERTILIZER DISSOLUTION FROM FERTILIZER POOLS
!
!     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=broadcast NH4,NH3,urea,NO3 fertilizer
!     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
!     RSN4AA,RSN4BA=rate of broadcast NH4 fertilizer dissoln in non-band,band
!     RSN3AA,RSN3BA=rate of broadcast NH3 fertilizer dissoln in non-band,band
!     RSNUAA,RSNUBA=rate of broadcast urea fertr dissoln in non-band,band
!     RSNOAA,RSNOBA=rate of broadcast NO3 fertilizer dissoln in non-band,band
!     RSN4BB=rate of banded NH4 fertilizer dissolution in band
!     RSN3BB=rate of banded NH3 fertilizer dissolution in band
!     RSNUBB=rate of banded urea fertilizer dissolution in band
!     RSNOBB=rate of banded NO3 fertilizer dissolution in band
!
          ZNH4FA(L,NY,NX)=ZNH4FA(L,NY,NX)-RSN4AA-RSN4BA
          ZNH3FA(L,NY,NX)=ZNH3FA(L,NY,NX)-RSN3AA-RSN3BA
          ZNHUFA(L,NY,NX)=ZNHUFA(L,NY,NX)-RSNUAA-RSNUBA
          ZNO3FA(L,NY,NX)=ZNO3FA(L,NY,NX)-RSNOAA-RSNOBA
          ZNH4FB(L,NY,NX)=ZNH4FB(L,NY,NX)-RSN4BB
          ZNH3FB(L,NY,NX)=ZNH3FB(L,NY,NX)-RSN3BB
          ZNHUFB(L,NY,NX)=ZNHUFB(L,NY,NX)-RSNUBB
          ZNO3FB(L,NY,NX)=ZNO3FB(L,NY,NX)-RSNOBB
!
!     ADD FERTILIZER DISSOLUTION TO ION FLUXES AND CONVERT TO MASS
!
!     TRN3G=NH3 dissolution from NH3
!     TRN4S,TRN4B=NH4 dissolution in non-band,band
!     TRN3S,TRN3B=NH3 dissolution from urea in non-band,band
!     TRNO3,TRNOB=NO3 dissolution in non-band,band
!
          TRN3G(L,NY,NX)=TRN3G(L,NY,NX)+RSN3AA+RSN3BA+RSN3BB
          TRN4S(L,NY,NX)=TRN4S(L,NY,NX)+RSN4AA
          TRN4B(L,NY,NX)=TRN4B(L,NY,NX)+RSN4BA+RSN4BB
          TRN3S(L,NY,NX)=TRN3S(L,NY,NX)+RSNUAA
          TRN3B(L,NY,NX)=TRN3B(L,NY,NX)+RSNUBA+RSNUBB
          TRNO3(L,NY,NX)=TRNO3(L,NY,NX)+RSNOAA
          TRNOB(L,NY,NX)=TRNOB(L,NY,NX)+RSNOBA+RSNOBB
          TRCO2(L,NY,NX)=TRCO2(L,NY,NX)*12.0
          TRN3G(L,NY,NX)=TRN3G(L,NY,NX)*14.0
          TRN4S(L,NY,NX)=TRN4S(L,NY,NX)*14.0
          TRN4B(L,NY,NX)=TRN4B(L,NY,NX)*14.0
          TRN3S(L,NY,NX)=TRN3S(L,NY,NX)*14.0
          TRN3B(L,NY,NX)=TRN3B(L,NY,NX)*14.0
          TRNO3(L,NY,NX)=TRNO3(L,NY,NX)*14.0
          TRNOB(L,NY,NX)=TRNOB(L,NY,NX)*14.0
          TRNO2(L,NY,NX)=TRNO2(L,NY,NX)*14.0
          TRN2B(L,NY,NX)=TRN2B(L,NY,NX)*14.0
          TRH1P(L,NY,NX)=TRH1P(L,NY,NX)*31.0
          TRH2P(L,NY,NX)=TRH2P(L,NY,NX)*31.0
          TRH1B(L,NY,NX)=TRH1B(L,NY,NX)*31.0
          TRH2B(L,NY,NX)=TRH2B(L,NY,NX)*31.0
        ENDIF
!     IF(I.EQ.116)THEN
!     WRITE(*,9984)'TRN3S',I,J,L,TRN4S(L,NY,NX),TRN3S(L,NY,NX)
!9984  FORMAT(A8,3I4,20F14.7)
!     ENDIF
      ENDDO
!
!     SURFACE RESIDUE
!
      call UpdateSoluteInSurfaceResidue(NX,NY)
!
!     TOTAL ION FLUXES FOR ALL REACTIONS ABOVE
!
!     RN4S=net NH4 flux in litter
!     RN3S=net NH3 flux in litter
!     RHP1,RHP2=net HPO4,H2PO4 flux in litter
!
      RN4S=RNH4-RXN4
      RN3S=-RNH4
      RHP1=-RH2P
      RHP2=RH2P-RPALPX-RPFEPX-RPCADX-2.0*RPCAMX-3.0*RPCAHX
!
!     CONVERT TOTAL ION FLUXES FROM CHANGES IN CONCENTRATION
!     TO CHANGES IN MASS PER UNIT AREA FOR USE IN 'REDIST'
!
!     TRN4S=total NH4 flux
!     TRN3S=total NH3 flux
!     TRH1P,TRH2P=net HPO4,H2PO4 flux
!     TRNX4=total NH4 adsorption
!     TRALPO,TRFEPO,TRCAPD,TRCAPH,TRCAPM
!     =total AlPO4,FePO4,CaHPO4,apatite,Ca(H2PO4)2 precipitation
!
      TRN4S(0,NY,NX)=TRN4S(0,NY,NX)+RN4S*VOLWM(NPH,0,NY,NX)
      TRN3S(0,NY,NX)=TRN3S(0,NY,NX)+RN3S*VOLWM(NPH,0,NY,NX)
      TRH1P(0,NY,NX)=TRH1P(0,NY,NX)+RHP1*VOLWM(NPH,0,NY,NX)
      TRH2P(0,NY,NX)=TRH2P(0,NY,NX)+RHP2*VOLWM(NPH,0,NY,NX)
      TRXN4(0,NY,NX)=TRXN4(0,NY,NX)+RXN4*VOLWM(NPH,0,NY,NX)
      TRALPO(0,NY,NX)=TRALPO(0,NY,NX)+RPALPX*VOLWM(NPH,0,NY,NX)
      TRFEPO(0,NY,NX)=TRFEPO(0,NY,NX)+RPFEPX*VOLWM(NPH,0,NY,NX)
      TRCAPD(0,NY,NX)=TRCAPD(0,NY,NX)+RPCADX*VOLWM(NPH,0,NY,NX)
      TRCAPH(0,NY,NX)=TRCAPH(0,NY,NX)+RPCAHX*VOLWM(NPH,0,NY,NX)
      TRCAPM(0,NY,NX)=TRCAPM(0,NY,NX)+RPCAMX*VOLWM(NPH,0,NY,NX)
      ZNH4FA(0,NY,NX)=ZNH4FA(0,NY,NX)-RSN4AA
      ZNH3FA(0,NY,NX)=ZNH3FA(0,NY,NX)-RSN3AA
      ZNHUFA(0,NY,NX)=ZNHUFA(0,NY,NX)-RSNUAA
      ZNO3FA(0,NY,NX)=ZNO3FA(0,NY,NX)-RSNOAA
      TRN4S(0,NY,NX)=TRN4S(0,NY,NX)+RSN4AA
      TRN3S(0,NY,NX)=TRN3S(0,NY,NX)+RSN3AA+RSNUAA
      TRNO3(0,NY,NX)=TRNO3(0,NY,NX)+RSNOAA
      TRN4S(0,NY,NX)=TRN4S(0,NY,NX)*14.0
      TRN3S(0,NY,NX)=TRN3S(0,NY,NX)*14.0
      TRNO3(0,NY,NX)=TRNO3(0,NY,NX)*14.0
      TRH1P(0,NY,NX)=TRH1P(0,NY,NX)*31.0
      TRH2P(0,NY,NX)=TRH2P(0,NY,NX)*31.0
!     WRITE(*,9989)'TRN4S',I,J,TRN4S(0,NY,NX)
!    2,RN4S,RNH4,RXN4,RSN4AA,VOLWM(NPH,0,NY,NX)
!    3,SPNH4,ZNH4FA(0,NY,NX),THETWR
!9989  FORMAT(A8,2I4,12E12.4)
    ENDDO
  ENDDO
  RETURN
  end subroutine solute
!------------------------------------------------------------------------

  subroutine UreaHydrolysis(L,NY,NX)

  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8) :: CNHUB
!      real(r8), intent(out) :: RSNUA
!      real(r8), intent(out) :: RSNUB
!      real(r8) :: CNHUA, COMA, DUKD
!      real(r8) :: DFNSA
!     begin_execution
!
!     UREA HYDROLYSIS IN BAND AND NON-BAND SOIL ZONES
!
!     VOLQ=biologically active soil water volume from nitro.f
!     COMA=concentration of active biomass
!     TOQCK=total microbial activity from nitro.f
!     DUKD=Km for urea hydrolysis
!
  IF(VOLQ(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    COMA=AMIN1(0.1E+06,TOQCK(L,NY,NX)/VOLQ(L,NY,NX))
  ELSE
    COMA=0.1E+06
  ENDIF
  DUKD=DUKM*(1.0+COMA/DUKI)
!
!     UREA HYDROLYSIS INHIBITION
!
!     ZNHU0,ZNHUI=initial,current inhibition activity
!     RNHUI=rate constants for decline in urea hydrolysis inhibition
!
  IF(ZNHU0(L,NY,NX).GT.ZEROS(NY,NX) &
    .AND.ZNHUI(L,NY,NX).GT.ZEROS(NY,NX))THEN
    ZNHUI(L,NY,NX)=ZNHUI(L,NY,NX) &
      -RNHUI(IUTYP(NY,NX))*ZNHUI(L,NY,NX) &
      *AMAX1(RNHUI(IUTYP(NY,NX)),1.0-ZNHUI(L,NY,NX)/ZNHU0(L,NY,NX))
  ELSE
    ZNHUI(L,NY,NX)=0._r8
  ENDIF
!
!     UREA CONCENTRATION AND HYDROLYSIS IN NON-BAND
!
!     ZNHUFA=urea fertilizer in non-band
!     CNHUA=concentration of urea fertilizer in non-band
!     DFNSA=effect of microbial concentration on urea hydrolysis in non-band
!     RSNUA=rate of urea hydrolysis in non-band
!     SPNHU=specific rate constant for urea hydrolysis
!     TFNQ=temperature effect on microbial activity from nitro.f
!     ZNHUI=current inhibition activity
!
  IF(ZNHUFA(L,NY,NX).GT.ZEROS(NY,NX) &
    .AND.BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
    CNHUA=ZNHUFA(L,NY,NX)/BKVL(L,NY,NX)
  ELSEIF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    CNHUA=ZNHUFA(L,NY,NX)/VOLW(L,NY,NX)
  ELSE
    CNHUA=0._r8
  ENDIF
  DFNSA=CNHUA/(CNHUA+DUKD)
  RSNUA=AMIN1(ZNHUFA(L,NY,NX),SPNHU*TOQCK(L,NY,NX)*DFNSA*TFNQ(L,NY,NX))*(1.0-ZNHUI(L,NY,NX))
!
!     UREA CONCENTRATION AND HYDROLYSIS IN BAND
!
!     ZNHUFB=urea fertilizer in band
!     CNHUB=concentration of urea fertilizer in band
!     DFNSB=effect of microbial concentration on urea hydrolysis in band
!     RSNUB=rate of urea hydrolysis in non-band
!     SPNHU=specific rate constant for urea hydrolysis
!     TFNQ=temperature effect on microbial activity from nitro.f
!
  IF(ZNHUFB(L,NY,NX).GT.ZEROS(NY,NX) &
    .AND.BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
    CNHUB=ZNHUFB(L,NY,NX)/BKVL(L,NY,NX)
  ELSEIF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    CNHUB=ZNHUFB(L,NY,NX)/VOLW(L,NY,NX)
  ELSE
    CNHUB=0._r8
  ENDIF
  DFNSB=CNHUB/(CNHUB+DUKD)
  RSNUB=AMIN1(ZNHUFB(L,NY,NX),SPNHU*TOQCK(L,NY,NX)*DFNSB*TFNQ(L,NY,NX))*(1.0-ZNHUI(L,NY,NX))
!     IF(ZNHUFA(L,NY,NX).GT.ZEROS(NY,NX))THEN
!     WRITE(*,8888)'UREA',I,J,L,IUTYP(NY,NX)
!    2,ZNHUFA(L,NY,NX),ZNHUFB(L,NY,NX),RSNUA,RSNUB
!    2,DFNSA,DFNSB,TFNQ(L,NY,NX),CNHUA,DUKD,DUKM,DUKI,TOQCK(L,NY,NX)
!    3,BKVL(L,NY,NX),SPNHU,ZNHU0(L,NY,NX),ZNHUI(L,NY,NX)
!    4,RNHUI(IUTYP(NY,NX)),VLNH4(L,NY,NX),VLNHB(L,NY,NX)
!    5,THETW(L,NY,NX)
!8888  FORMAT(A8,4I4,40E12.4)
!     ENDIF
  end subroutine UreaHydrolysis
!------------------------------------------------------------------------

  subroutine UpdateSoilFertlizer(L,NY,NX)

  implicit none
  integer, intent(in) :: L,NY,NX

  real(r8) :: RH2BX,RNBX,R3BX,RH1BX
  real(r8) :: VOLWPX
!      real(r8), intent(out) :: XN4B, XN41
!      real(r8) :: RN4X,RN3X
!      real(r8) :: VOLWNX
!      real(r8) :: RNBX
!      real(r8) :: R3BX
!     begin_execution

  call UreaHydrolysis(L,NY,NX)
!
!     NH4, NH3, UREA, NO3 DISSOLUTION IN BAND AND NON-BAND
!     SOIL ZONES FROM FIRST-ORDER FUNCTIONS OF REMAINING
!     FERTILIZER (NOTE: SUPERPHOSPHATE AND ROCK PHOSPHATE
!     ARE REPRESENTED AS MONOCALCIUM PHOSPHATE AND HYDROXYAPATITE
!     MODELLED IN PHOSPHORUS REACTIONS BELOW)
!
!     RSN4AA,RSN4BA=rate of broadcast NH4 fertilizer dissoln in non-band,band
!     RSN3AA,RSN3BA=rate of broadcast NH3 fertilizer dissoln in non-band,band
!     RSNUAA,RSNUBA=rate of broadcast urea fertr dissoln in non-band,band
!     RSNOAA,RSNOBA=rate of broadcast NO3 fertilizer dissoln in non-band,band
!     RSN4BB= rate of banded NH4 fertilizer dissolution in band
!     RSN3BB= rate of banded NH3 fertilizer dissolution in band
!     RSNUBB= rate of banded urea fertilizer dissolution in band
!     RSNOBB= rate of banded NO3 fertilizer dissolution in band
!     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=broadcast NH4,NH3,urea,NO3 fertilizer
!     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
!     VLNH4,VLNHB=fractions of soil volume in NH4 non-band,band
!     VLNO3,VLNOB=fractions of soil volume in N03 non-band,band
!     VLPO4,VLPOB=fractions of soil volume in H2PO4 non-band,band
!     THETW=soil water concentration
!
  RSN4AA=SPNH4*ZNH4FA(L,NY,NX)*VLNH4(L,NY,NX)*THETW(L,NY,NX)
  RSN3AA=SPNH3*ZNH3FA(L,NY,NX)*VLNH4(L,NY,NX)
  RSNUAA=RSNUA*VLNH4(L,NY,NX)*THETW(L,NY,NX)
  RSNOAA=SPNO3*ZNO3FA(L,NY,NX)*VLNO3(L,NY,NX)*THETW(L,NY,NX)
  RSN4BA=SPNH4*ZNH4FA(L,NY,NX)*VLNHB(L,NY,NX)*THETW(L,NY,NX)
  RSN3BA=SPNH3*ZNH3FA(L,NY,NX)*VLNHB(L,NY,NX)
  RSNUBA=RSNUA*VLNHB(L,NY,NX)*THETW(L,NY,NX)
  RSNOBA=SPNO3*ZNO3FA(L,NY,NX)*VLNOB(L,NY,NX)*THETW(L,NY,NX)
  RSN4BB=SPNH4*ZNH4FB(L,NY,NX)*THETW(L,NY,NX)
  RSN3BB=SPNH3*ZNH3FB(L,NY,NX)
  RSNUBB=RSNUB*VLNHB(L,NY,NX)*THETW(L,NY,NX)
  RSNOBB=SPNO3*ZNO3FB(L,NY,NX)*THETW(L,NY,NX)
!
!     SOLUBLE AND EXCHANGEABLE NH4 CONCENTRATIONS
!     IN NON-BAND AND BAND SOIL ZONES
!
!     VOLWNH,VOLWNB=water volume in NH4 non-band,band
!     RN4X,RN3X=NH4,NH3 input from uptake, mineraln, dissoln in non-band
!     RNBX,R3BX=NH4,NH3 input from uptake, mineraln, dissoln in band
!     TUPNH4,TUPNH3=soil-root exchange of NH4,NH3 in non-band from uptake.f
!     TUPNHB,TUPN3B=soil-root exchange of NH4,NH3 in band from uptake.f
!     XNH4S,XNH4B=net change in NH4 in band,non-band from nitro.f
!     CN41,CN31,CN4B,CN3B=total NH4,NH3 concentration in non-band,band
!     XN41,XN4B=adsorbed NH4 concentration in non-band,band
!
  IF(VOLWNH.GT.ZEROS2(NY,NX))THEN
    VOLWNX=14.0*VOLWNH
    RN4X=(-TUPNH4(L,NY,NX)+XNH4S(L,NY,NX)+14.0*RSN4AA)/VOLWNX
    RN3X=(-TUPN3S(L,NY,NX)+14.0*RSNUAA)/VOLWNX
    CN41=AMAX1(0.0,ZNH4S(L,NY,NX)/VOLWNX+RN4X)
    CN31=AMAX1(0.0,ZNH3S(L,NY,NX)/VOLWNX+RN3X)
    XN41=AMAX1(0.0,XN4(L,NY,NX)/BKVLNH)
  ELSE
    RN4X=0._r8
    RN3X=0._r8
    CN41=0._r8
    CN31=0._r8
    XN41=0._r8
  ENDIF
  IF(VOLWNB.GT.ZEROS2(NY,NX))THEN
    VOLWNX=14.0*VOLWNB
    RNBX=(-TUPNHB(L,NY,NX)+XNH4B(L,NY,NX)+14.0*(RSN4BA+RSN4BB))/VOLWNX
    R3BX=(-TUPN3B(L,NY,NX)+14.0*(RSNUBA+RSNUBB))/VOLWNX
    CN4B=AMAX1(0.0,ZNH4B(L,NY,NX)/VOLWNX+RNBX)
    CN3B=AMAX1(0.0,ZNH3B(L,NY,NX)/VOLWNX+R3BX)
    XN4B=AMAX1(0.0,XNB(L,NY,NX)/BKVLNB)
  ELSE
    RNBX=0._r8
    R3BX=0._r8
    CN4B=0._r8
    CN3B=0._r8
    XN4B=0._r8
  ENDIF
!
!     SOLUBLE, EXCHANGEABLE AND PRECIPITATED PO4 CONCENTRATIONS IN
!     NON-BAND AND BAND SOIL ZONES
!
!     VOLWPO,VOLWPB=water volume in H2PO4 non-band,band
!     RH1PX,RH2PX=HPO4,H2PO4 inputs from mineraln, uptake in non-band
!     RH1BX,RH2BX=HPO4,H2PO4 inputs from mineraln, uptake in band
!     XH1PS,XH1BS=net change in HPO4 in band,non-band from nitro.f
!     TUPH1P,TUPH2P=soil-root exch of HPO4,H2PO4 in non-band from uptake.f
!     TUPH1B,TUPH2B=soil-root exch of HPO4,H2PO4 in band from uptake.f
!     CH1P1,CH2P1=HPO4,H2PO4 concentrations in non-band
!     CH1PB,CH2PB=HPO4,H2PO4 concentrations in band
!     XOH11,XOH11,XOH21=concn of adsorption sites R-,R-OH,R-OH2 in non-band
!     XOH1B,XH11B,XH21B=concn of adsorption sites R-,R-OH,R-OH2 in band
!     XH1P1,XH2P1=concentration of adsorbed HPO4,H2PO4 in non-band
!     XH11B,X2P1B=concentration of adsorbed HPO4,H2PO4 in band
!     PALPO1,PFEPO1=concn of precip AlPO4,FEPO4 in non-band
!     PALPOB,PFEPOB=concn of precip AlPO4,FEPO4 in band
!     PCAPM1,PCAPD1,PCAPH1=concn of precip CaH2PO4,CaHPO4,apatite in non-band
!     PCAPMB,PCAPDB,PCAPHB=concn of precip CaH2PO4,CaHPO4,apatite in band
!
  IF(VOLWPO.GT.ZEROS2(NY,NX))THEN
    VOLWPX=31.0*VOLWPO
    RH1PX=(XH1PS(L,NY,NX)-TUPH1P(L,NY,NX))/VOLWPX
    RH2PX=(XH2PS(L,NY,NX)-TUPH2P(L,NY,NX))/VOLWPX
    CH1P1=AMAX1(0.0,H1PO4(L,NY,NX)/VOLWPX+RH1PX)
    CH2P1=AMAX1(0.0,H2PO4(L,NY,NX)/VOLWPX+RH2PX)
    XOH01=AMAX1(0.0,XOH0(L,NY,NX))/BKVLPO
    XOH11=AMAX1(0.0,XOH1(L,NY,NX))/BKVLPO
    XOH21=AMAX1(0.0,XOH2(L,NY,NX))/BKVLPO
    XH1P1=AMAX1(0.0,XH1P(L,NY,NX))/BKVLPO
    XH2P1=AMAX1(0.0,XH2P(L,NY,NX))/BKVLPO
    PALPO1=AMAX1(0.0,PALPO(L,NY,NX))/BKVLPO
    PFEPO1=AMAX1(0.0,PFEPO(L,NY,NX))/BKVLPO
    PCAPM1=AMAX1(0.0,PCAPM(L,NY,NX))/BKVLPO
    PCAPD1=AMAX1(0.0,PCAPD(L,NY,NX))/BKVLPO
    PCAPH1=AMAX1(0.0,PCAPH(L,NY,NX))/BKVLPO
  ELSE
    RH1PX=0._r8
    RH2PX=0._r8
    CH1P1=0._r8
    CH2P1=0._r8
    XOH01=0._r8
    XOH11=0._r8
    XOH21=0._r8
    XH1P1=0._r8
    XH2P1=0._r8
    PALPO1=0._r8
    PFEPO1=0._r8
    PCAPM1=0._r8
    PCAPD1=0._r8
    PCAPH1=0._r8
  ENDIF
  IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
    VOLWPX=31.0*VOLWPB
    RH1BX=(XH1BS(L,NY,NX)-TUPH1B(L,NY,NX))/VOLWPX
    RH2BX=(XH2BS(L,NY,NX)-TUPH2B(L,NY,NX))/VOLWPX
    CH1PB=AMAX1(0.0,H1POB(L,NY,NX)/VOLWPX+RH1BX)
    CH2PB=AMAX1(0.0,H2POB(L,NY,NX)/VOLWPX+RH2BX)
    XH01B=AMAX1(0.0,XOH0B(L,NY,NX))/BKVLPB
    XH11B=AMAX1(0.0,XOH1B(L,NY,NX))/BKVLPB
    XH21B=AMAX1(0.0,XOH2B(L,NY,NX))/BKVLPB
    X1P1B=AMAX1(0.0,XH1PB(L,NY,NX))/BKVLPB
    X2P1B=AMAX1(0.0,XH2PB(L,NY,NX))/BKVLPB
    PALPOB=AMAX1(0.0,PALPB(L,NY,NX))/BKVLPB
    PFEPOB=AMAX1(0.0,PFEPB(L,NY,NX))/BKVLPB
    PCAPMB=AMAX1(0.0,PCPMB(L,NY,NX))/BKVLPB
    PCAPDB=AMAX1(0.0,PCPDB(L,NY,NX))/BKVLPB
    PCAPHB=AMAX1(0.0,PCPHB(L,NY,NX))/BKVLPB
  ELSE
    RH1BX=0._r8
    RH2BX=0._r8
    CH1PB=0._r8
    CH2PB=0._r8
    XH01B=0._r8
    XH11B=0._r8
    XH21B=0._r8
    X1P1B=0._r8
    X2P1B=0._r8
    PALPOB=0._r8
    PFEPOB=0._r8
    PCAPMB=0._r8
    PCAPDB=0._r8
    PCAPHB=0._r8
  ENDIF
  end subroutine UpdateSoilFertlizer

!------------------------------------------------------------------------------------------

  subroutine UpdateNH3FertilizerBandinfo(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX
!     begin_execution
!     IFNHB=banded NH4 fertilizer flag
!     ROWN=NH4 fertilizer band row width
!     DPNH4=NH4 fertilizer band depth
!
  IF(IFNHB(NY,NX).EQ.1.AND.ROWN(NY,NX).GT.0.0)THEN
    IF(L.EQ.NU(NY,NX).OR.CDPTH(L-1,NY,NX).LT.DPNH4(NY,NX))THEN
!
!     NH4 BAND WIDTH
!
!     DWNH4=change in NH4 fertilizer band width
!     WDNHB=layer NH4 fertilizer band width
!     ZNSGL=NH4 diffusivity
!     TORT=tortuosity
!
      DWNH4=0.5*SQRT(ZNSGL(L,NY,NX))*TORT(NPH,L,NY,NX)
      WDNHB(L,NY,NX)=AMIN1(ROWN(NY,NX),AMAX1(0.025,WDNHB(L,NY,NX))+DWNH4)
!
!     NH4 BAND DEPTH
!
!     DPFLW=change in NH4 fertilizer band depth
!     DPNH4,DPNHB=total,layer NH4 fertilizer band depth
!
      IF(CDPTH(L,NY,NX).GE.DPNH4(NY,NX))THEN
        DPFLW=FLWD+DWNH4
        DPNH4(NY,NX)=DPNH4(NY,NX)+DPFLW
        DPNHB(L,NY,NX)=DPNHB(L,NY,NX)+DPFLW
        IF(DPNHB(L,NY,NX).GT.DLYR(3,L,NY,NX))THEN
          DPNHB(L+1,NY,NX)=DPNHB(L+1,NY,NX)+(DPNHB(L,NY,NX)-DLYR(3,L,NY,NX))
          WDNHB(L+1,NY,NX)=WDNHB(L,NY,NX)
          DPNHB(L,NY,NX)=DLYR(3,L,NY,NX)
        ELSEIF(DPNHB(L,NY,NX).LT.0.0)THEN
          DPNHB(L-1,NY,NX)=DPNHB(L-1,NY,NX)+DPNHB(L,NY,NX)
          DPNHB(L,NY,NX)=0._r8
          WDNHB(L,NY,NX)=0._r8
        ENDIF
      ENDIF
!
!     FRACTION OF SOIL LAYER OCCUPIED BY NH4 BAND
!     FROM BAND WIDTH X DEPTH
!
!     VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
!     DLYR=soil layer thickness
!     FVLNH4=relative change in VLNH4
!
      XVLNH4=VLNH4(L,NY,NX)
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
        VLNHB(L,NY,NX)=AMAX1(0.0,AMIN1(0.999,WDNHB(L,NY,NX) &
          /ROWN(NY,NX)*DPNHB(L,NY,NX)/DLYR(3,L,NY,NX)))
      ELSE
        VLNHB(L,NY,NX)=0._r8
      ENDIF
      VLNH4(L,NY,NX)=1._r8-VLNHB(L,NY,NX)
      FVLNH4=AMIN1(0.0,(VLNH4(L,NY,NX)-XVLNH4)/XVLNH4)
!
!     TRANSFER NH4, NH3 FROM NON-BAND TO BAND
!     DURING BAND GROWTH
!
!     DNH4S,DNH3S,DXNH4=transfer of NH4,NH3,exchangeable NH4
!
      DNH4S=FVLNH4*ZNH4S(L,NY,NX)/14.0
      DNH3S=FVLNH4*ZNH3S(L,NY,NX)/14.0
      DXNH4=FVLNH4*XN4(L,NY,NX)
      TRN4S(L,NY,NX)=TRN4S(L,NY,NX)+DNH4S
      TRN4B(L,NY,NX)=TRN4B(L,NY,NX)-DNH4S
      TRN3S(L,NY,NX)=TRN3S(L,NY,NX)+DNH3S
      TRN3B(L,NY,NX)=TRN3B(L,NY,NX)-DNH3S
      TRXN4(L,NY,NX)=TRXN4(L,NY,NX)+DXNH4
      TRXNB(L,NY,NX)=TRXNB(L,NY,NX)-DXNH4
    ELSE
!
!     AMALGAMATE NH4 BAND WITH NON-BAND IF BAND NO LONGER EXISTS
!
      DPNHB(L,NY,NX)=0._r8
      WDNHB(L,NY,NX)=0._r8
      VLNH4(L,NY,NX)=1._r8
      VLNHB(L,NY,NX)=0._r8
      ZNH4S(L,NY,NX)=ZNH4S(L,NY,NX)+ZNH4B(L,NY,NX)
      ZNH3S(L,NY,NX)=ZNH3S(L,NY,NX)+ZNH3B(L,NY,NX)
      ZNH4B(L,NY,NX)=0._r8
      ZNH3B(L,NY,NX)=0._r8
      XN4(L,NY,NX)=XN4(L,NY,NX)+XNB(L,NY,NX)
      XNB(L,NY,NX)=0._r8
    ENDIF
  ENDIF
  end subroutine UpdateNH3FertilizerBandinfo
!------------------------------------------------------------------------------------------

  subroutine UpdatePO4FertilizerBandinfo(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8) :: XVLPO4

!     begin_execution
!     IFPOB=banded H2PO4 fertilizer flag
!     ROWP=H2PO4 fertilizer band row width
!     DPPO4=H2PO4 fertilizer band depth
!
  IF(IFPOB(NY,NX).EQ.1.AND.ROWP(NY,NX).GT.0.0)THEN
    IF(L.EQ.NU(NY,NX).OR.CDPTH(L-1,NY,NX).LT.DPPO4(NY,NX))THEN
!
!     PO4 BAND WIDTH
!     DWPO4=change in H2PO4 fertilizer band width
!     WDPO4=layer H2PO4 fertilizer band width
!     POSGL=H2PO4 diffusivity
!     TORT=tortuosity
!
      DWPO4=0.5*SQRT(POSGL(L,NY,NX))*TORT(NPH,L,NY,NX)
      WDPOB(L,NY,NX)=AMIN1(ROWP(NY,NX),WDPOB(L,NY,NX)+DWPO4)
!
!     PO4 BAND DEPTH
!
!     DPFLW=change in H2PO4 fertilizer band depth
!     DPPO4,DPPOB=total,layer H2PO4 fertilizer band depth
!
      IF(CDPTH(L,NY,NX).GE.DPPO4(NY,NX))THEN
        DPFLW=FLWD+DWPO4
        DPPO4(NY,NX)=DPPO4(NY,NX)+DPFLW
        DPPOB(L,NY,NX)=DPPOB(L,NY,NX)+DPFLW
        IF(DPPOB(L,NY,NX).GT.DLYR(3,L,NY,NX))THEN
          DPPOB(L+1,NY,NX)=DPPOB(L+1,NY,NX)+(DPPOB(L,NY,NX)-DLYR(3,L,NY,NX))
          WDPOB(L+1,NY,NX)=WDPOB(L,NY,NX)
          DPPOB(L,NY,NX)=DLYR(3,L,NY,NX)
        ELSEIF(DPPOB(L,NY,NX).LT.0.0)THEN
          DPPOB(L-1,NY,NX)=DPPOB(L-1,NY,NX)+DPPOB(L,NY,NX)
          DPPOB(L,NY,NX)=0._r8
          WDPOB(L,NY,NX)=0._r8
        ENDIF
      ENDIF
!
!     FRACTION OF SOIL LAYER OCCUPIED BY PO4 BAND
!     FROM BAND WIDTH X DEPTH
!
!     VLPO43,VLPOB=fraction of soil volume in H2PO4 non-band,band
!     DLYR=soil layer thickness
!     FVLPO4=relative change in VLPO4
!
      XVLPO4=VLPO4(L,NY,NX)
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
        VLPOB(L,NY,NX)=AMAX1(0.0,AMIN1(0.999,WDPOB(L,NY,NX) &
          /ROWP(NY,NX)*DPPOB(L,NY,NX)/DLYR(3,L,NY,NX)))
      ELSE
        VLPOB(L,NY,NX)=0._r8
      ENDIF
      VLPO4(L,NY,NX)=1._r8-VLPOB(L,NY,NX)
      FVLPO4=AMIN1(0.0,(VLPO4(L,NY,NX)-XVLPO4)/XVLPO4)
!
!     TRANSFER HPO4,H2PO4 FROM NON-BAND TO BAND
!     DURING BAND GROWTH DEPENDING ON SALT
!     VS. NON-SALT OPTION
!
!     DZ*,DX*,DP*=transfer of solute,adsorbed,precipitated HPO4,H2PO4
!
      IF(ISALTG.NE.0)THEN
        DZH0P=FVLPO4*H0PO4(L,NY,NX)
        DZH1P=FVLPO4*H1PO4(L,NY,NX)/31.0
        DZH2P=FVLPO4*H2PO4(L,NY,NX)/31.0
        DZH3P=FVLPO4*H3PO4(L,NY,NX)
        DZF1P=FVLPO4*ZFE1P(L,NY,NX)
        DZF2P=FVLPO4*ZFE2P(L,NY,NX)
        DZC0P=FVLPO4*ZCA0P(L,NY,NX)
        DZC1P=FVLPO4*ZCA1P(L,NY,NX)
        DZC2P=FVLPO4*ZCA2P(L,NY,NX)
        DZM1P=FVLPO4*ZMG1P(L,NY,NX)
        DXOH0=FVLPO4*XOH0(L,NY,NX)
        DXOH1=FVLPO4*XOH1(L,NY,NX)
        DXOH2=FVLPO4*XOH2(L,NY,NX)
        DXH1P=FVLPO4*XH1P(L,NY,NX)
        DXH2P=FVLPO4*XH2P(L,NY,NX)
        DPALP=FVLPO4*PALPO(L,NY,NX)
        DPFEP=FVLPO4*PFEPO(L,NY,NX)
        DPCDP=FVLPO4*PCAPD(L,NY,NX)
        DPCHP=FVLPO4*PCAPH(L,NY,NX)
        DPCMP=FVLPO4*PCAPM(L,NY,NX)
        TRH0P(L,NY,NX)=TRH0P(L,NY,NX)+DZH0P
        TRH1P(L,NY,NX)=TRH1P(L,NY,NX)+DZH1P
        TRH2P(L,NY,NX)=TRH2P(L,NY,NX)+DZH2P
        TRH3P(L,NY,NX)=TRH3P(L,NY,NX)+DZH3P
        TRF1P(L,NY,NX)=TRF1P(L,NY,NX)+DZF1P
        TRF2P(L,NY,NX)=TRF2P(L,NY,NX)+DZF2P
        TRC0P(L,NY,NX)=TRC0P(L,NY,NX)+DZC0P
        TRC1P(L,NY,NX)=TRC1P(L,NY,NX)+DZC1P
        TRC2P(L,NY,NX)=TRC2P(L,NY,NX)+DZC2P
        TRM1P(L,NY,NX)=TRM1P(L,NY,NX)+DZM1P
        TRH0B(L,NY,NX)=TRH0B(L,NY,NX)-DZH0P
        TRH1B(L,NY,NX)=TRH1B(L,NY,NX)-DZH1P
        TRH2B(L,NY,NX)=TRH2B(L,NY,NX)-DZH2P
        TRH3B(L,NY,NX)=TRH3B(L,NY,NX)-DZH3P
        TRF1B(L,NY,NX)=TRF1B(L,NY,NX)-DZF1P
        TRF2B(L,NY,NX)=TRF2B(L,NY,NX)-DZF2P
        TRC0B(L,NY,NX)=TRC0B(L,NY,NX)-DZC0P
        TRC1B(L,NY,NX)=TRC1B(L,NY,NX)-DZC1P
        TRC2B(L,NY,NX)=TRC2B(L,NY,NX)-DZC2P
        TRM1B(L,NY,NX)=TRM1B(L,NY,NX)-DZM1P
        TRXH0(L,NY,NX)=TRXH0(L,NY,NX)+DXOH0
        TRXH1(L,NY,NX)=TRXH1(L,NY,NX)+DXOH1
        TRXH2(L,NY,NX)=TRXH2(L,NY,NX)+DXOH2
        TRX1P(L,NY,NX)=TRX1P(L,NY,NX)+DXH1P
        TRX2P(L,NY,NX)=TRX2P(L,NY,NX)+DXH2P
        TRBH0(L,NY,NX)=TRBH0(L,NY,NX)-DXOH0
        TRBH1(L,NY,NX)=TRBH1(L,NY,NX)-DXOH1
        TRBH2(L,NY,NX)=TRBH2(L,NY,NX)-DXOH2
        TRB1P(L,NY,NX)=TRB1P(L,NY,NX)-DXH1P
        TRB2P(L,NY,NX)=TRB2P(L,NY,NX)-DXH2P
        TRALPO(L,NY,NX)=TRALPO(L,NY,NX)+DPALP
        TRFEPO(L,NY,NX)=TRFEPO(L,NY,NX)+DPFEP
        TRCAPD(L,NY,NX)=TRCAPD(L,NY,NX)+DPCDP
        TRCAPH(L,NY,NX)=TRCAPH(L,NY,NX)+DPCHP
        TRCAPM(L,NY,NX)=TRCAPM(L,NY,NX)+DPCMP
        TRALPB(L,NY,NX)=TRALPB(L,NY,NX)-DPALP
        TRFEPB(L,NY,NX)=TRFEPB(L,NY,NX)-DPFEP
        TRCPDB(L,NY,NX)=TRCPDB(L,NY,NX)-DPCDP
        TRCPHB(L,NY,NX)=TRCPHB(L,NY,NX)-DPCHP
        TRCPMB(L,NY,NX)=TRCPMB(L,NY,NX)-DPCMP
      ELSE
        DZH1P=FVLPO4*H1PO4(L,NY,NX)/31.0
        DZH2P=FVLPO4*H2PO4(L,NY,NX)/31.0
        DXOH1=FVLPO4*XOH1(L,NY,NX)
        DXOH2=FVLPO4*XOH2(L,NY,NX)
        DXH2P=FVLPO4*XH2P(L,NY,NX)
        DPALP=FVLPO4*PALPO(L,NY,NX)
        DPFEP=FVLPO4*PFEPO(L,NY,NX)
        DPCDP=FVLPO4*PCAPD(L,NY,NX)
        DPCHP=FVLPO4*PCAPH(L,NY,NX)
        DPCMP=FVLPO4*PCAPM(L,NY,NX)
        TRH1P(L,NY,NX)=TRH1P(L,NY,NX)+DZH1P
        TRH2P(L,NY,NX)=TRH2P(L,NY,NX)+DZH2P
        TRXH1(L,NY,NX)=TRXH1(L,NY,NX)+DXOH1
        TRXH2(L,NY,NX)=TRXH2(L,NY,NX)+DXOH2
        TRX2P(L,NY,NX)=TRX2P(L,NY,NX)+DXH2P
        TRH1B(L,NY,NX)=TRH1B(L,NY,NX)-DZH1P
        TRH2B(L,NY,NX)=TRH2B(L,NY,NX)-DZH2P
        TRBH1(L,NY,NX)=TRBH1(L,NY,NX)-DXOH1
        TRBH2(L,NY,NX)=TRBH2(L,NY,NX)-DXOH2
        TRB2P(L,NY,NX)=TRB2P(L,NY,NX)-DXH2P
        TRALPO(L,NY,NX)=TRALPO(L,NY,NX)+DPALP
        TRFEPO(L,NY,NX)=TRFEPO(L,NY,NX)+DPFEP
        TRCAPD(L,NY,NX)=TRCAPD(L,NY,NX)+DPCDP
        TRCAPH(L,NY,NX)=TRCAPH(L,NY,NX)+DPCHP
        TRCAPM(L,NY,NX)=TRCAPM(L,NY,NX)+DPCMP
        TRALPB(L,NY,NX)=TRALPB(L,NY,NX)-DPALP
        TRFEPB(L,NY,NX)=TRFEPB(L,NY,NX)-DPFEP
        TRCPDB(L,NY,NX)=TRCPDB(L,NY,NX)-DPCDP
        TRCPHB(L,NY,NX)=TRCPHB(L,NY,NX)-DPCHP
        TRCPMB(L,NY,NX)=TRCPMB(L,NY,NX)-DPCMP
      ENDIF
    ELSE
!
!     AMALGAMATE PO4 BAND WITH NON-BAND IF BAND NO LONGER EXISTS
!
      DPPOB(L,NY,NX)=0._r8
      WDPOB(L,NY,NX)=0._r8
      VLPOB(L,NY,NX)=0._r8
      VLPO4(L,NY,NX)=1._r8
      H0PO4(L,NY,NX)=H0PO4(L,NY,NX)+H0POB(L,NY,NX)
      H1PO4(L,NY,NX)=H1PO4(L,NY,NX)+H1POB(L,NY,NX)
      H2PO4(L,NY,NX)=H2PO4(L,NY,NX)+H2POB(L,NY,NX)
      H3PO4(L,NY,NX)=H3PO4(L,NY,NX)+H3POB(L,NY,NX)
      ZFE1P(L,NY,NX)=ZFE1P(L,NY,NX)+ZFE1PB(L,NY,NX)
      ZFE2P(L,NY,NX)=ZFE2P(L,NY,NX)+ZFE2PB(L,NY,NX)
      ZCA0P(L,NY,NX)=ZCA0P(L,NY,NX)+ZCA0PB(L,NY,NX)
      ZCA1P(L,NY,NX)=ZCA1P(L,NY,NX)+ZCA1PB(L,NY,NX)
      ZCA2P(L,NY,NX)=ZCA2P(L,NY,NX)+ZCA2PB(L,NY,NX)
      ZMG1P(L,NY,NX)=ZMG1P(L,NY,NX)+ZMG1PB(L,NY,NX)
      H0POB(L,NY,NX)=0._r8
      H1POB(L,NY,NX)=0._r8
      H2POB(L,NY,NX)=0._r8
      H3POB(L,NY,NX)=0._r8
      ZFE1PB(L,NY,NX)=0._r8
      ZFE2PB(L,NY,NX)=0._r8
      ZCA0PB(L,NY,NX)=0._r8
      ZCA1PB(L,NY,NX)=0._r8
      ZCA2PB(L,NY,NX)=0._r8
      ZMG1PB(L,NY,NX)=0._r8
      XOH0(L,NY,NX)=XOH0(L,NY,NX)+XOH0B(L,NY,NX)
      XOH1(L,NY,NX)=XOH1(L,NY,NX)+XOH1B(L,NY,NX)
      XOH2(L,NY,NX)=XOH2(L,NY,NX)+XOH2B(L,NY,NX)
      XH1P(L,NY,NX)=XH1P(L,NY,NX)+XH1PB(L,NY,NX)
      XH2P(L,NY,NX)=XH2P(L,NY,NX)+XH2PB(L,NY,NX)
      XOH0B(L,NY,NX)=0._r8
      XOH1B(L,NY,NX)=0._r8
      XOH2B(L,NY,NX)=0._r8
      XH1PB(L,NY,NX)=0._r8
      XH2PB(L,NY,NX)=0._r8
      PALPO(L,NY,NX)=PALPO(L,NY,NX)+PALPB(L,NY,NX)
      PFEPO(L,NY,NX)=PFEPO(L,NY,NX)+PFEPB(L,NY,NX)
      PCAPD(L,NY,NX)=PCAPD(L,NY,NX)+PCPDB(L,NY,NX)
      PCAPH(L,NY,NX)=PCAPH(L,NY,NX)+PCPHB(L,NY,NX)
      PCAPM(L,NY,NX)=PCAPM(L,NY,NX)+PCPMB(L,NY,NX)
      PALPB(L,NY,NX)=0._r8
      PFEPB(L,NY,NX)=0._r8
      PCPDB(L,NY,NX)=0._r8
      PCPHB(L,NY,NX)=0._r8
      PCPMB(L,NY,NX)=0._r8
    ENDIF
  ENDIF
  end subroutine UpdatePO4FertilizerBandinfo
!------------------------------------------------------------------------------------------

  subroutine UpdateNO3FertilizerBandinfo(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8) :: XVLNO3
!     begin_execution
!     IFNOB=banded NO3 fertilizer flag
!     ROWO=NO3 fertilizer band row width
!     DPNO3=NO3 fertilizer band depth
!
  IF(IFNOB(NY,NX).EQ.1.AND.ROWO(NY,NX).GT.0.0)THEN

    IF(L.EQ.NU(NY,NX).OR.CDPTH(L-1,NY,NX).LT.DPNO3(NY,NX))THEN
!
!     NO3 BAND WIDTH
!
!     DWNO3=change in NO3 fertilizer band width
!     WDNOB=layer NO3 fertilizer band width
!     ZOSGL=NO3 diffusivity
!     TORT=tortuosity
!
      DWNO3=0.5*SQRT(ZOSGL(L,NY,NX))*TORT(NPH,L,NY,NX)
      WDNOB(L,NY,NX)=AMIN1(ROWO(NY,NX),WDNOB(L,NY,NX)+DWNO3)
!
!     NO3 BAND DEPTH
!
!     DPFLW=change in NO3 fertilizer band depth
!     DPNO3,DPNOB=total,layer NO3 fertilizer band depth
!
      IF(CDPTH(L,NY,NX).GE.DPNO3(NY,NX))THEN
        DPFLW=FLWD+DWNO3
        DPNO3(NY,NX)=DPNO3(NY,NX)+DPFLW
        DPNOB(L,NY,NX)=DPNOB(L,NY,NX)+DPFLW
        IF(DPNOB(L,NY,NX).GT.DLYR(3,L,NY,NX))THEN
          DPNOB(L+1,NY,NX)=DPNOB(L+1,NY,NX)+(DPNOB(L,NY,NX)-DLYR(3,L,NY,NX))
          WDNOB(L+1,NY,NX)=WDNOB(L,NY,NX)
          DPNOB(L,NY,NX)=DLYR(3,L,NY,NX)
        ELSEIF(DPNOB(L,NY,NX).LT.0.0)THEN
          DPNOB(L-1,NY,NX)=DPNOB(L-1,NY,NX)+DPNOB(L,NY,NX)
          DPNOB(L,NY,NX)=0._r8
          WDNOB(L,NY,NX)=0._r8
        ENDIF
      ENDIF
!
!     FRACTION OF SOIL LAYER OCCUPIED BY NO3 BAND
!     FROM BAND WIDTH X DEPTH
!
!     VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
!     DLYR=soil layer thickness
!     FVLNO3=relative change in VLNO3
!
      XVLNO3=VLNO3(L,NY,NX)
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
        VLNOB(L,NY,NX)=AMAX1(0.0,AMIN1(0.999,WDNOB(L,NY,NX) &
          /ROWO(NY,NX)*DPNOB(L,NY,NX)/DLYR(3,L,NY,NX)))
      ELSE
        VLNOB(L,NY,NX)=0._r8
      ENDIF
      VLNO3(L,NY,NX)=1._r8-VLNOB(L,NY,NX)
      FVLNO3=AMIN1(0.0,(VLNO3(L,NY,NX)-XVLNO3)/XVLNO3)
!
!     TRANSFER NO3 FROM NON-BAND TO BAND
!     DURING BAND GROWTH
!
!     DNO3S,DNO2S=transfer of NO3,NO2
!
      DNO3S=FVLNO3*ZNO3S(L,NY,NX)/14.0
      DNO2S=FVLNO3*ZNO2S(L,NY,NX)/14.0
      TRNO3(L,NY,NX)=TRNO3(L,NY,NX)+DNO3S
      TRNO2(L,NY,NX)=TRNO2(L,NY,NX)+DNO2S
      TRNOB(L,NY,NX)=TRNOB(L,NY,NX)-DNO3S
      TRN2B(L,NY,NX)=TRN2B(L,NY,NX)-DNO2S
    ELSE
!
!     AMALGAMATE NO3 BAND WITH NON-BAND IF BAND NO LONGER EXISTS
!
      DPNOB(L,NY,NX)=0._r8
      WDNOB(L,NY,NX)=0._r8
      VLNO3(L,NY,NX)=1._r8
      VLNOB(L,NY,NX)=0._r8
      ZNO3S(L,NY,NX)=ZNO3S(L,NY,NX)+ZNO3B(L,NY,NX)
      ZNO2S(L,NY,NX)=ZNO2S(L,NY,NX)+ZNO2B(L,NY,NX)
      ZNO3B(L,NY,NX)=0._r8
      ZNO2B(L,NY,NX)=0._r8
    ENDIF
  ENDIF
  end subroutine UpdateNO3FertilizerBandinfo
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteInSurfaceResidue(NX,NY)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: CAL1,CALX,CCA1,CCAX
  real(r8) :: CCEC,CCEC0,CCO31
  real(r8) :: CFE1,CFEX
  real(r8) :: CH2PA,CH2PD,CH2PF,CH2PH,CH2PM
  real(r8) :: CHY1,COH1,DP
  real(r8) :: FX,S0,S1,XALQ,XCAQ,XCAX
  real(r8) :: XFEQ,XHYQ,XN4Q,XTLQ
!     begin_execution
!     BKVL=litter mass
!
  IF(VOLWM(NPH,0,NY,NX).GT.ZEROS2(NY,NX))THEN
    BKVLX=BKVL(0,NY,NX)
!
!     UREA HYDROLYSIS IN SURFACE RESIDUE
!
!     VOLQ=biologically active litter water volume from nitro.f
!     COMA=concentration of active biomass
!     TOQCK=total microbial activity from nitro.f
!     DUKD=Km for urea hydrolysis
!
    IF(VOLQ(0,NY,NX).GT.ZEROS2(NY,NX))THEN
      COMA=AMIN1(0.1E+06,TOQCK(0,NY,NX)/VOLQ(0,NY,NX))
    ELSE
      COMA=0.1E+06
    ENDIF
    DUKD=DUKM*(1.0+COMA/DUKI)
!
!     UREA HYDROLYSIS INHIBITION
!
!     ZNHU0,ZNHUI=initial,current inhibition activity
!     RNHUI=rate constant for decline in urea hydrolysis inhibition
!
    IF(ZNHU0(0,NY,NX).GT.ZEROS(NY,NX) &
      .AND.ZNHUI(0,NY,NX).GT.ZEROS(NY,NX))THEN
      ZNHUI(0,NY,NX)=ZNHUI(0,NY,NX)-TFNQ(0,NY,NX)**0.25 &
        *RNHUI(IUTYP(NY,NX))*ZNHUI(0,NY,NX) &
        *AMAX1(RNHUI(IUTYP(NY,NX)),1.0-ZNHUI(0,NY,NX)/ZNHU0(0,NY,NX))
    ELSE
      ZNHUI(0,NY,NX)=0._r8
    ENDIF
!
!     UREA CONCENTRATION AND HYDROLYSIS IN SURFACE RESIDUE
!
!     ZNHUFA=urea fertilizer
!     CNHUA=concentration of urea fertilizer
!     DFNSA=effect of microbial concentration on urea hydrolysis
!     RSNUA=rate of urea hydrolysis
!     SPNHU=specific rate constant for urea hydrolysis
!     TFNQ=temperature effect on microbial activity from nitro.f
!
    IF(ZNHUFA(0,NY,NX).GT.ZEROS(NY,NX) &
      .AND.BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
      CNHUA=ZNHUFA(0,NY,NX)/BKVL(0,NY,NX)
      DFNSA=CNHUA/(CNHUA+DUKD)
      RSNUA=AMIN1(ZNHUFA(0,NY,NX),SPNHU*TOQCK(0,NY,NX)*DFNSA*TFNQ(0,NY,NX))*(1.0-ZNHUI(0,NY,NX))
    ELSE
      RSNUA=0._r8
    ENDIF
!     IF(ZNHUFA(0,NY,NX).GT.ZEROS(NY,NX))THEN
!     WRITE(*,8778)'UREA0',I,J,IUTYP(NY,NX)
!    2,ZNHUFA(0,NY,NX),RSNUA
!    2,DFNSA,TFNQ(0,NY,NX),CNHUA,DUKD,DUKM,DUKI,TOQCK(0,NY,NX)
!    3,BKVL(0,NY,NX),TFNQ(0,NY,NX),SPNHU,ZNHU0(0,NY,NX),ZNHUI(0,NY,NX)
!    4,RNHUI(IUTYP(NY,NX))
!8778  FORMAT(A8,3I4,40E12.4)
!     ENDIF
!
!     NH4, NH3, UREA, NO3 DISSOLUTION IN SURFACE RESIDUE
!     FROM FIRST-ORDER FUNCTIONS OF REMAINING
!     FERTILIZER (NOTE: SUPERPHOSPHATE AND ROCK PHOSPHATE
!     ARE REPRESENTED AS MONOCALCIUM PHOSPHATE AND HYDROXYAPATITE
!     MODELLED IN PHOSPHORUS REACTIONS BELOW)
!
!     RSN4AA=rate of broadcast NH4 fertilizer dissoln
!     RSN3AA=rate of broadcast NH3 fertilizer dissoln
!     RSNUAA=rate of broadcast urea fertr dissoln
!     RSNOAA=rate of broadcast NO3 fertilizer dissoln
!
    IF(VOLWRX(NY,NX).GT.ZEROS(NY,NX))THEN
      THETWR=AMIN1(1.0,VOLW(0,NY,NX)/VOLWRX(NY,NX))
    ELSE
      THETWR=1._r8
    ENDIF
    RSN4AA=SPNH4*ZNH4FA(0,NY,NX)*THETWR
    RSN3AA=SPNH3*ZNH3FA(0,NY,NX)
    RSNUAA=RSNUA*THETWR
    RSNOAA=SPNO3*ZNO3FA(0,NY,NX)*THETWR
!
!     SOLUBLE AND EXCHANGEABLE NH4 CONCENTRATIONS
!
!     VOLWNH=water volume
!     RN4X,RN3X=NH4,NH3 input from uptake, mineraln, dissoln
!     XNH4S=net change in NH4 from nitro.f
!     CN41,CN31=total NH4,NH3 concentration
!     XN41=adsorbed NH4 concentration
!
    IF(VOLWM(NPH,0,NY,NX).GT.ZEROS2(NY,NX))THEN
      VOLWMX=14.0*VOLWM(NPH,0,NY,NX)
      RN4X=(XNH4S(0,NY,NX)+14.0*RSN4AA)/VOLWMX
      RN3X=14.0*RSNUAA/VOLWMX
      CN41=AMAX1(ZERO,ZNH4S(0,NY,NX)/VOLWMX+RN4X)
      CN31=AMAX1(ZERO,ZNH3S(0,NY,NX)/VOLWMX+RN3X)
      IF(BKVLX.GT.ZEROS(NY,NX))THEN
        XN41=AMAX1(ZERO,XN4(0,NY,NX)/BKVLX)
      ELSE
        XN41=0._r8
      ENDIF
!
!     SOLUBLE, EXCHANGEABLE AND PRECIPITATED PO4 CONCENTRATIONS
!
!     VOLWMP=water volume
!     RH1PX,RH2PX=HPO4,H2PO4 inputs from mineraln, uptake
!     XH1PS=net change in HPO4 from nitro.f
!     CH1P1,CH2P1=HPO4,H2PO4 concentrations
!
      VOLWMP=31.0*VOLWM(NPH,0,NY,NX)
      RH1PX=XH1PS(0,NY,NX)/VOLWMP
      RH2PX=XH2PS(0,NY,NX)/VOLWMP
      CH1P1=AMAX1(0.0,H1PO4(0,NY,NX)/VOLWMP+RH1PX)
      CH2P1=AMAX1(0.0,H2PO4(0,NY,NX)/VOLWMP+RH2PX)
    ELSE
      RN4X=0._r8
      RN3X=0._r8
      CN41=0._r8
      CN31=0._r8
      XN41=0._r8
      RH1PX=0._r8
      RH2PX=0._r8
      CH1P1=0._r8
      CH2P1=0._r8
    ENDIF
!
!     PHOSPHORUS TRANSFORMATIONS IN SURFACE RESIDUE
!
!     PALPO1,PFEPO1=concn of precip AlPO4,FEPO4
!     PCAPM1,PCAPD1,PCAPH1=concn of precip CaH2PO4,CaHPO4,apatite
!
    IF(BKVLX.GT.ZEROS(NY,NX))THEN
      PALPO1=AMAX1(0.0,PALPO(0,NY,NX)/BKVLX)
      PFEPO1=AMAX1(0.0,PFEPO(0,NY,NX)/BKVLX)
      PCAPM1=AMAX1(0.0,PCAPM(0,NY,NX)/BKVLX)
      PCAPD1=AMAX1(0.0,PCAPD(0,NY,NX)/BKVLX)
      PCAPH1=AMAX1(0.0,PCAPH(0,NY,NX)/BKVLX)
    ELSE
      PALPO1=0._r8
      PFEPO1=0._r8
      PCAPM1=0._r8
      PCAPD1=0._r8
      PCAPH1=0._r8
    ENDIF
!
!     CALCULATE H2PO4 COPRECIPITATES FRPM LITTER PH
!
    CHY1=AMAX1(ZERO,10.0**(-(PH(0,NY,NX)-3.0)))
    COH1=AMAX1(ZERO,DPH2O/CHY1)
    CAL1=AMAX1(ZERO,SPALO/COH1**3)
    CFE1=AMAX1(ZERO,SPFEO/COH1**3)
    CCO20=AMAX1(ZERO,CCO2S(0,NY,NX)/12.0)
    CCO31=AMAX1(ZERO,CCO20*DPCO3/CHY1**2)
    CCA1=AMAX1(ZERO,AMIN1(CCAMX,SPCAC/CCO31))
!
!     ALUMINUM PHOSPHATE (VARISCITE)
!
!     CH2PA,CH2P1=equilibrium,current H2PO4 concentration in litter
!     SYA0P2=solubility product derived from SPALO
!     RPALPX=H2PO4 dissolution from AlPO4 in litter
!
    CH2PA=SYA0P2/(CAL1*COH1**2)
    RPALPX=AMIN1(AMAX1(0.0,4.0E-08*ORGC(0,NY,NX)-PALPO1),AMAX1(-PALPO1,TPD*(CH2P1-CH2PA)))
!
!     IRON PHOSPHATE (STRENGITE)
!
!     CH2PF,CH2P1=equilibrium,current H2PO4 concentration in litter
!     SYF0P2=solubility product derived from SPALO
!     RPFEPX=H2PO4 dissolution from FePO4 in litter
!
    CH2PF=SYF0P2/(CFE1*COH1**2)
    RPFEPX=AMIN1(AMAX1(0.0,2.0E-06*ORGC(0,NY,NX)-PFEPO1),AMAX1(-PFEPO1,TPD*(CH2P1-CH2PF)))
!
!     DICALCIUM PHOSPHATE
!
!     CH2PD,CH2P1=equilibrium,current H2PO4 concentration in litter
!     SYCAD2=solubility product derived from SPALO
!     RPCADX=H2PO4 dissolution from CaHPO4 in litter
!
    CH2PD=SYCAD2/(CCA1*COH1)
    RPCADX=AMIN1(AMAX1(0-.0,5.0E-05*ORGC(0,NY,NX)-PCAPD1),AMAX1(-PCAPD1,TPD*(CH2P1-CH2PD)))
!
!     HYDROXYAPATITE
!
!     CH2PH,CH2P1=equilibrium,current H2PO4 concentration in litter
!     SYCAH2=solubility product derived from SPALO
!     RPCAHX=H2PO4 dissolution from apatite in litter
!
    CH2PH=(SYCAH2/(CCA1**5*COH1**7))**0.333
    RPCAHX=AMIN1(AMAX1(0.0,5.0E-05*ORGC(0,NY,NX)-PCAPH1),AMAX1(-PCAPH1,TPD*(CH2P1-CH2PH)))
!
!     MONOCALCIUM PHOSPHATE
!
!     CH2PM,CH2P1=equilibrium,current H2PO4 concentration in litter
!     SPCAM=solubility product for Ca(H2PO4)2
!     RPCAMX=H2PO4 dissolution from Ca(H2PO4)2 in litter
!
    CH2PM=SQRT(SPCAM/CCA1)
    RPCAMX=AMIN1(AMAX1(0.0,5.0E-05*ORGC(0,NY,NX)-PCAPM1),AMAX1(-PCAPM1*SPPO4,TPD*(CH2P1-CH2PM)))
!     IF(I.GT.315)THEN
!     WRITE(*,2227)'RPPO4',I,J,L,RPCAHX,CH2P1,CH2PA,CH2PH
!    2,SYA0P2,CAL1,COH1,SYCAH2,CCA1,CCO21,CCO31,PCAPH1
!    3,VOLWM(NPH,0,NY,NX),SPCAC/CCO31,H2PO4(0,NY,NX)
!    4,CCO20,DPCO3,CHY1,CCO2S(0,NY,NX)
!2227  FORMAT(A8,3I4,20E12.4)
!     ENDIF
!
!     PHOSPHORUS ANION EXCHANGE IN SURFACE REDISUE
!     CALCULATED FROM EXCHANGE EQUILIBRIA AMONG H2PO4-,
!     HPO4--, H+, OH- AND PROTONATED AND NON-PROTONATED -OH
!     EXCHANGE SITES (NOT CALCULATED)
!
!     H2PO4-H+HPO4
!
!     DPH2P=dissociation constant
!     S1=equilibrium concentration in litter
!     RH2P=H2PO4-H+HPO4 dissociation in litter
!
    DP=DPH2P
    S0=CH1P1+CHY1+DP
    S1=AMAX1(0.0,S0**2-4.0*(CH1P1*CHY1-DP*CH2P1))
    RH2P=TSL*(S0-SQRT(S1))
!
!     EQUILIBRIUM X-CA CONCENTRATION FROM CEC, GAPON COEFFICIENTS
!     AND CATION CONCENTRATIONS
!
!     CCEC,XCEC=cation exchange concentration,capacity
!     XCAX=equilibrium R-Ca concentration
!     GKC4,GKCH,GKCA,GKCM,GKCN,GKCK=Gapon selectivity coefficients for
!     CA-NH4,CA-H,CA-AL,CA-MG,CA-NA,CA-K
!     X*Q=equilibrium exchangeable concentrations
!     XTLQ=total equilibrium exchangeable concentration
!
    IF(BKVLX.GT.ZEROS(NY,NX))THEN
      CCEC0=AMAX1(ZERO,COOH*ORGC(0,NY,NX)/BKVLX)
    ELSE
      CCEC0=ZERO
    ENDIF
    CALX=AMAX1(ZERO,CAL1)**0.333
    CFEX=AMAX1(ZERO,CFE1)**0.333
    CCAX=AMAX1(ZERO,CCA1)**0.500
!
!     EQUILIBRIUM X-CA CONCENTRATION FROM CEC AND CATION
!     CONCENTRATIONS
!
    XCAX=CCEC0/(1.0+GKC4(NU(NY,NX),NY,NX)*CN41/CCAX &
      +GKCH(NU(NY,NX),NY,NX)*CHY1/CCAX &
      +GKCA(NU(NY,NX),NY,NX)*CALX/CCAX &
      +GKCA(NU(NY,NX),NY,NX)*CFEX/CCAX)
    XN4Q=XCAX*CN41*GKC4(NU(NY,NX),NY,NX)
    XHYQ=XCAX*CHY1*GKCH(NU(NY,NX),NY,NX)
    XALQ=XCAX*CALX*GKCA(NU(NY,NX),NY,NX)
    XFEQ=XCAX*CFEX*GKCA(NU(NY,NX),NY,NX)
    XCAQ=XCAX*CCAX
    XTLQ=XN4Q+XHYQ+XALQ+XFEQ+XCAQ
    IF(XTLQ.GT.ZERO)THEN
      FX=CCEC0/XTLQ
    ELSE
      FX=0._r8
    ENDIF
    XN4Q=FX*XN4Q
!
!     NH4 AND NH3 EXCHANGE IN SURFACE RESIDUE
!
!     RXN4=NH4 adsorption in litter
!     TADC0=adsorption rate constant
!     RNH4=NH4-NH3+H dissociation in litter
!     DPN4=NH4 dissociation constant
!
    RXN4=TADC0*(XN4Q-XN41)
    RNH4=(CHY1*CN31-DPN4*CN41)/(DPN4+CHY1)
!     IF(J.EQ.12)THEN
!     WRITE(*,2223)'RXN4',I,J,RXN4,CN41,XN41,TADC0,XN4Q
!    2,CCAX,CCA1,CCO20,CCO31
!    2,XCAQ,CCEC0,FN4X,FCAQ,GKC4(NU(NY,NX),NY,NX)
!    3,PH(0,NY,NX),CHY1,RNH4
!    3,CN31,DPN4,ZNH4S(0,NY,NX),XN4(0,NY,NX),14.0*RSN4AA,RN4X,BKVLX
!    4,BKVL(0,NY,NX),VOLWM(NPH,0,NY,NX)
!2223  FORMAT(A8,2I4,30E12.4)
!     ENDIF
  ELSE
    RSN4AA=0._r8
    RSN3AA=0._r8
    RSNUAA=0._r8
    RSNOAA=0._r8
    RPALPX=0._r8
    RPFEPX=0._r8
    RPCADX=0._r8
    RPCAHX=0._r8
    RPCAMX=0._r8
    RXN4=0._r8
    RNH4=0._r8
    RH2P=0._r8
    RPALPX=0._r8
    RPFEPX=0._r8
    RPCADX=0._r8
    RPCAMX=0._r8
    RPCAHX=0._r8
  ENDIF
  end subroutine UpdateSoluteInSurfaceResidue

END module SoluteMod
