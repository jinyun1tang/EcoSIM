module SoluteMod
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use minimathmod, only : test_aeqb,AZMAX1,AZMIN1
  use SoluteParMod
  use SOMDataType
  use EcosimConst
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

  real(r8) :: RSN4AA,RSN4BA,RSN3AA,RSN3BA,RSNUAA,RSNUBA,RSNOAA
  real(r8) :: RSNOBA,RSN4BB,RSN3BB,RSNUBB,RSNOBB
  real(r8) :: RSNUA, RSNUB

  public :: UreaHydrolysis
  public :: UpdateSoilFertlizer
  public :: UpdateFertilizerBand
  public :: GeoChemEquilibria
  public :: UpdateSoluteinSurfaceResidue
  contains


!------------------------------------------------------------------------

  subroutine GeoChemEquilibria(chemvar, solflx)

  implicit none
  type(chem_var_type), intent(inout) :: chemvar
  type(solute_flx_type), intent(inout) :: solflx

  IF(ISALTG.NE.0)THEN

    call SaltChemEquilibria(chemvar,solflx)

  ELSE
    call NoSaltChemEquilibria(chemvar,solflx)
  ENDIF

  end subroutine GeoChemEquilibria


!------------------------------------------------------------------------

  subroutine UpdateFertilizerBand(L,NY,NX)

  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8) :: FLWD
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
  call UpdateNH3FertilizerBandinfo(L,NY,NX,FLWD)
!
!     NO3 FERTILIZER BAND
!
  call UpdateNO3FertilizerBandinfo(L,NY,NX,FLWD)
!
!     PO4 FERTILIZER BAND
!
  call UpdatePO4FertilizerBandinfo(L,NY,NX,FLWD)
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
  FertN_soil(ifert_nh4,L,NY,NX)=FertN_soil(ifert_nh4,L,NY,NX)-RSN4AA-RSN4BA
  FertN_soil(ifert_nh3,L,NY,NX)=FertN_soil(ifert_nh3,L,NY,NX)-RSN3AA-RSN3BA
  FertN_soil(ifert_urea,L,NY,NX)=FertN_soil(ifert_urea,L,NY,NX)-RSNUAA-RSNUBA
  FertN_soil(ifert_no3,L,NY,NX)=FertN_soil(ifert_no3,L,NY,NX)-RSNOAA-RSNOBA
  FertN_band(ifert_nh4_band,L,NY,NX)=FertN_band(ifert_nh4_band,L,NY,NX)-RSN4BB
  FertN_band(ifert_nh3_band,L,NY,NX)=FertN_band(ifert_nh3_band,L,NY,NX)-RSN3BB
  FertN_band(ifert_urea_band,L,NY,NX)=FertN_band(ifert_urea_band,L,NY,NX)-RSNUBB
  FertN_band(ifert_no3_band,L,NY,NX)=FertN_band(ifert_no3_band,L,NY,NX)-RSNOBB
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
  TRCO2(L,NY,NX)=TRCO2(L,NY,NX)*catomw
  TRN3G(L,NY,NX)=TRN3G(L,NY,NX)*natomw
  TRN4S(L,NY,NX)=TRN4S(L,NY,NX)*natomw
  TRN4B(L,NY,NX)=TRN4B(L,NY,NX)*natomw
  TRN3S(L,NY,NX)=TRN3S(L,NY,NX)*natomw
  TRN3B(L,NY,NX)=TRN3B(L,NY,NX)*natomw
  TRNO3(L,NY,NX)=TRNO3(L,NY,NX)*natomw
  TRNOB(L,NY,NX)=TRNOB(L,NY,NX)*natomw
  TRNO2(L,NY,NX)=TRNO2(L,NY,NX)*natomw
  TRN2B(L,NY,NX)=TRN2B(L,NY,NX)*natomw
  TRH1P(L,NY,NX)=TRH1P(L,NY,NX)*patomw
  TRH2P(L,NY,NX)=TRH2P(L,NY,NX)*patomw
  TRH1B(L,NY,NX)=TRH1B(L,NY,NX)*patomw
  TRH2B(L,NY,NX)=TRH2B(L,NY,NX)*patomw
  end subroutine UpdateFertilizerBand

!------------------------------------------------------------------------

  subroutine UreaHydrolysis(L,NY,NX)

  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8) :: CNHUB,COMA
  real(r8) :: DFNSB,CNHUA,DFNSA,DUKD
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
    COMA=AMIN1(0.1E+06_r8,TOQCK(L,NY,NX)/VOLQ(L,NY,NX))
  ELSE
    COMA=0.1E+06_r8
  ENDIF
  DUKD=DUKM*(1.0_r8+COMA/DUKI)
!
!     UREA HYDROLYSIS INHIBITION
!
!     ZNHU0,ZNHUI=initial,current inhibition activity
!     RNHUI=rate constants for decline in urea hydrolysis inhibition
!
  IF(ZNHU0(L,NY,NX).GT.ZEROS(NY,NX).AND.ZNHUI(L,NY,NX).GT.ZEROS(NY,NX))THEN
    ZNHUI(L,NY,NX)=ZNHUI(L,NY,NX)-RNHUI(IUTYP(NY,NX))*ZNHUI(L,NY,NX) &
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
  IF(FertN_soil(ifert_urea,L,NY,NX).GT.ZEROS(NY,NX).AND.BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
    CNHUA=FertN_soil(ifert_urea,L,NY,NX)/BKVL(L,NY,NX)
  ELSEIF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    CNHUA=FertN_soil(ifert_urea,L,NY,NX)/VOLW(L,NY,NX)
  ELSE
    CNHUA=0._r8
  ENDIF
  DFNSA=CNHUA/(CNHUA+DUKD)
  RSNUA=AMIN1(FertN_soil(ifert_urea,L,NY,NX),SPNHU*TOQCK(L,NY,NX)*DFNSA*TFNQ(L,NY,NX))*(1.0-ZNHUI(L,NY,NX))
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
  IF(FertN_band(ifert_urea_band,L,NY,NX).GT.ZEROS(NY,NX).AND.BKVL(L,NY,NX).GT.ZEROS(NY,NX))THEN
    CNHUB=FertN_band(ifert_urea_band,L,NY,NX)/BKVL(L,NY,NX)
  ELSEIF(VOLW(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    CNHUB=FertN_band(ifert_urea_band,L,NY,NX)/VOLW(L,NY,NX)
  ELSE
    CNHUB=0._r8
  ENDIF
  DFNSB=CNHUB/(CNHUB+DUKD)
  RSNUB=AMIN1(FertN_band(ifert_urea_band,L,NY,NX),SPNHU*TOQCK(L,NY,NX)*DFNSB*TFNQ(L,NY,NX))*(1.0-ZNHUI(L,NY,NX))

  end subroutine UreaHydrolysis
!------------------------------------------------------------------------

  subroutine UpdateSoilFertlizer(L,NY,NX,chemvar)

  implicit none
  integer, intent(in) :: L,NY,NX
  type(chem_var_type), target, intent(inout) :: chemvar
  real(r8) :: RH2BX,RNBX,R3BX,RH1BX
  real(r8) :: VOLWPX,VOLWNX
  real(r8) :: RH2PX,RH1PX,RN3X,RN4X
  real(r8), pointer :: VOLWNH
  real(r8), pointer :: BKVLNH
  real(r8), pointer :: BKVLNB
  real(r8), pointer :: BKVLPO
  real(r8), pointer :: BKVLPB
  real(r8), pointer :: VOLWNB
  real(r8), pointer :: VOLWPO
  real(r8), pointer :: CH1P1
  real(r8), pointer :: CH1PB
  real(r8), pointer :: CH2P1
  real(r8), pointer :: CH2PB
  real(r8), pointer :: CN31
  real(r8), pointer :: CN3B
  real(r8), pointer :: CN41
  real(r8), pointer :: CN4B
  real(r8), pointer :: PALPO1
  real(r8), pointer :: PALPOB
  real(r8), pointer :: PCAPD1
  real(r8), pointer :: PCAPDB
  real(r8), pointer :: PCAPH1
  real(r8), pointer :: PCAPHB
  real(r8), pointer :: PCAPM1
  real(r8), pointer :: PCAPMB
  real(r8), pointer :: PFEPO1
  real(r8), pointer :: PFEPOB
  real(r8), pointer :: X1P1B
  real(r8), pointer :: X2P1B
  real(r8), pointer :: XH11B
  real(r8), pointer :: XH1P1
  real(r8), pointer :: XH21B
  real(r8), pointer :: XH2P1
  real(r8), pointer :: XN41
  real(r8), pointer :: XN4B
  real(r8), pointer :: XOH11
  real(r8), pointer :: XOH21
  real(r8), pointer :: XH01B
  real(r8), pointer :: XOH01
  real(r8), pointer :: VOLWPB

  XOH01  =>  chemvar%XOH01
  XH01B  =>  chemvar%XH01B
  XOH21  =>  chemvar%XOH21
  XOH11  =>  chemvar%XOH11
  XH21B  =>  chemvar%XH21B
  XH1P1  =>  chemvar%XH1P1
  XH11B  =>  chemvar%XH11B
  X1P1B  =>  chemvar%X1P1B
  X2P1B  =>  chemvar%X2P1B
  CH2PB  =>  chemvar%CH2PB
  XN41   =>  chemvar%XN41
  XN4B   =>  chemvar%XN4B
  CN31   =>  chemvar%CN31
  CN3B   =>  chemvar%CN3B
  CN41   =>  chemvar%CN41
  CN4B   =>  chemvar%CN4B
  PALPO1 =>  chemvar%PALPO1
  PALPOB =>  chemvar%PALPOB
  PCAPD1 =>  chemvar%PCAPD1
  PCAPDB =>  chemvar%PCAPDB
  PCAPH1 =>  chemvar%PCAPH1
  PCAPHB =>  chemvar%PCAPHB
  PCAPM1 =>  chemvar%PCAPM1
  PCAPMB =>  chemvar%PCAPMB
  PFEPO1 =>  chemvar%PFEPO1
  PFEPOB =>  chemvar%PFEPOB
  CH2P1  =>  chemvar%CH2P1
  CH1PB  =>  chemvar%CH1PB
  CH1P1  =>  chemvar%CH1P1
  XH2P1  =>  chemvar%XH2P1
  VOLWNH =>  chemvar%VOLWNH
  BKVLNH =>  chemvar%BKVLNH
  BKVLNB =>  chemvar%BKVLNB
  BKVLPO =>  chemvar%BKVLPO
  BKVLPB =>  chemvar%BKVLPB
  VOLWNB =>  chemvar%VOLWNB
  VOLWPB =>  chemvar%VOLWPB
  VOLWPO =>  chemvar%VOLWPO
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
  RSN4AA=SPNH4*FertN_soil(ifert_nh4,L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)*THETW(L,NY,NX)
  RSN3AA=SPNH3*FertN_soil(ifert_nh3,L,NY,NX)*trcs_VLN(ids_NH4,L,NY,NX)
  RSNUAA=RSNUA*trcs_VLN(ids_NH4,L,NY,NX)*THETW(L,NY,NX)
  RSNOAA=SPNO3*FertN_soil(ifert_no3,L,NY,NX)*trcs_VLN(ids_NO3,L,NY,NX)*THETW(L,NY,NX)
  RSN4BA=SPNH4*FertN_soil(ifert_nh4,L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)*THETW(L,NY,NX)
  RSN3BA=SPNH3*FertN_soil(ifert_nh3,L,NY,NX)*trcs_VLN(ids_NH4B,L,NY,NX)
  RSNUBA=RSNUA*trcs_VLN(ids_NH4B,L,NY,NX)*THETW(L,NY,NX)
  RSNOBA=SPNO3*FertN_soil(ifert_no3,L,NY,NX)*trcs_VLN(ids_NO3B,L,NY,NX)*THETW(L,NY,NX)
  RSN4BB=SPNH4*FertN_band(ifert_nh4_band,L,NY,NX)*THETW(L,NY,NX)
  RSN3BB=SPNH3*FertN_band(ifert_nh3_band,L,NY,NX)
  RSNUBB=RSNUB*trcs_VLN(ids_NH4B,L,NY,NX)*THETW(L,NY,NX)
  RSNOBB=SPNO3*FertN_band(ifert_no3_band,L,NY,NX)*THETW(L,NY,NX)
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
    VOLWNX=natomw*VOLWNH
    RN4X=(-TUPNH4(L,NY,NX)+XNH4S(L,NY,NX)+natomw*RSN4AA)/VOLWNX
    RN3X=(-TUPN3S(L,NY,NX)+natomw*RSNUAA)/VOLWNX
    CN41=AZMAX1(trc_solml(ids_NH4,L,NY,NX)/VOLWNX+RN4X)
    CN31=AZMAX1(trc_solml(idg_NH3,L,NY,NX)/VOLWNX+RN3X)
    XN41=AZMAX1(trcx_solml(idx_NH4,L,NY,NX)/BKVLNH)
  ELSE
    RN4X=0._r8
    RN3X=0._r8
    CN41=0._r8
    CN31=0._r8
    XN41=0._r8
  ENDIF
  IF(VOLWNB.GT.ZEROS2(NY,NX))THEN
    VOLWNX=natomw*VOLWNB
    RNBX=(-TUPNHB(L,NY,NX)+XNH4B(L,NY,NX)+natomw*(RSN4BA+RSN4BB))/VOLWNX
    R3BX=(-TUPN3B(L,NY,NX)+natomw*(RSNUBA+RSNUBB))/VOLWNX
    CN4B=AZMAX1(trc_solml(ids_NH4B,L,NY,NX)/VOLWNX+RNBX)
    CN3B=AZMAX1(trc_solml(idg_NH3B,L,NY,NX)/VOLWNX+R3BX)
    XN4B=AZMAX1(trcx_solml(idx_NH4B,L,NY,NX)/BKVLNB)
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
    VOLWPX=patomw*VOLWPO
    RH1PX=(XH1PS(L,NY,NX)-TUPH1P(L,NY,NX))/VOLWPX
    RH2PX=(XH2PS(L,NY,NX)-TUPH2P(L,NY,NX))/VOLWPX
    CH1P1=AZMAX1(trc_solml(ids_H1PO4,L,NY,NX)/VOLWPX+RH1PX)
    CH2P1=AZMAX1(trc_solml(ids_H2PO4,L,NY,NX)/VOLWPX+RH2PX)
    XOH01=AZMAX1(trcx_solml(idx_OHe,L,NY,NX))/BKVLPO
    XOH11=AZMAX1(trcx_solml(idx_OH,L,NY,NX))/BKVLPO
    XOH21=AZMAX1(trcx_solml(idx_OHp,L,NY,NX))/BKVLPO
    XH1P1=AZMAX1(trcx_solml(idx_HPO4,L,NY,NX))/BKVLPO
    XH2P1=AZMAX1(trcx_solml(idx_H2PO4,L,NY,NX))/BKVLPO
    PALPO1=AZMAX1(trcp_salml(idsp_AlPO4,L,NY,NX))/BKVLPO
    PFEPO1=AZMAX1(trcp_salml(idsp_FePO4,L,NY,NX))/BKVLPO
    PCAPM1=AZMAX1(trcp_salml(idsp_CaH2PO4,L,NY,NX))/BKVLPO
    PCAPD1=AZMAX1(trcp_salml(idsp_CaHPO4,L,NY,NX))/BKVLPO
    PCAPH1=AZMAX1(trcp_salml(idsp_HA,L,NY,NX))/BKVLPO
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
    VOLWPX=patomw*VOLWPB
    RH1BX=(XH1BS(L,NY,NX)-TUPH1B(L,NY,NX))/VOLWPX
    RH2BX=(XH2BS(L,NY,NX)-TUPH2B(L,NY,NX))/VOLWPX
    CH1PB=AZMAX1(trc_solml(ids_H1PO4B,L,NY,NX)/VOLWPX+RH1BX)
    CH2PB=AZMAX1(trc_solml(ids_H2PO4B,L,NY,NX)/VOLWPX+RH2BX)
    XH01B=AZMAX1(trcx_solml(idx_OHeB,L,NY,NX))/BKVLPB
    XH11B=AZMAX1(trcx_solml(idx_OHB,L,NY,NX))/BKVLPB
    XH21B=AZMAX1(trcx_solml(idx_OHpB,L,NY,NX))/BKVLPB
    X1P1B=AZMAX1(trcx_solml(idx_HPO4B,L,NY,NX))/BKVLPB
    X2P1B=AZMAX1(trcx_solml(idx_H2PO4B,L,NY,NX))/BKVLPB
    PALPOB=AZMAX1(trcp_salml(idsp_AlPO4B,L,NY,NX))/BKVLPB
    PFEPOB=AZMAX1(trcp_salml(idsp_FePO4B,L,NY,NX))/BKVLPB
    PCAPMB=AZMAX1(trcp_salml(idsp_CaH2PO4B,L,NY,NX))/BKVLPB
    PCAPDB=AZMAX1(trcp_salml(idsp_CaHPO4B,L,NY,NX))/BKVLPB
    PCAPHB=AZMAX1(trcp_salml(idsp_HAB,L,NY,NX))/BKVLPB
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

  subroutine UpdateNH3FertilizerBandinfo(L,NY,NX,FLWD)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: FLWD

  real(r8) :: XVLNH4,DWNH4,DXNH4
  real(r8) :: DPFLW,DNH4S,DNH3S
  real(r8) :: FVLNH4
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
      XVLNH4=trcs_VLN(ids_NH4,L,NY,NX)
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
        trcs_VLN(ids_NH4B,L,NY,NX)=AZMAX1(AMIN1(0.999_r8,WDNHB(L,NY,NX) &
          /ROWN(NY,NX)*DPNHB(L,NY,NX)/DLYR(3,L,NY,NX)))
      ELSE
        trcs_VLN(ids_NH4B,L,NY,NX)=0._r8
      ENDIF
      trcs_VLN(ids_NH4,L,NY,NX)=1._r8-trcs_VLN(ids_NH4B,L,NY,NX)
      trcs_VLN(idg_NH3,L,NY,NX)=trcs_VLN(ids_NH4,L,NY,NX)
      trcs_VLN(idg_NH3B,L,NY,NX)=trcs_VLN(ids_NH4B,L,NY,NX)
      FVLNH4=AZMIN1((trcs_VLN(ids_NH4,L,NY,NX)-XVLNH4)/XVLNH4)
!
!     TRANSFER NH4, NH3 FROM NON-BAND TO BAND
!     DURING BAND GROWTH
!
!     DNH4S,DNH3S,DXNH4=transfer of NH4,NH3,exchangeable NH4
!
      DNH4S=FVLNH4*trc_solml(ids_NH4,L,NY,NX)/natomw
      DNH3S=FVLNH4*trc_solml(idg_NH3,L,NY,NX)/natomw
      DXNH4=FVLNH4*trcx_solml(idx_NH4,L,NY,NX)
      TRN4S(L,NY,NX)=TRN4S(L,NY,NX)+DNH4S
      TRN4B(L,NY,NX)=TRN4B(L,NY,NX)-DNH4S
      TRN3S(L,NY,NX)=TRN3S(L,NY,NX)+DNH3S
      TRN3B(L,NY,NX)=TRN3B(L,NY,NX)-DNH3S
      trcx_TR(idx_NH4,L,NY,NX)=trcx_TR(idx_NH4,L,NY,NX)+DXNH4
      trcx_TR(idx_NH4B,L,NY,NX)=trcx_TR(idx_NH4B,L,NY,NX)-DXNH4
    ELSE
!
!     AMALGAMATE NH4 BAND WITH NON-BAND IF BAND NO LONGER EXISTS
!
      DPNHB(L,NY,NX)=0._r8
      WDNHB(L,NY,NX)=0._r8
      trcs_VLN(ids_NH4,L,NY,NX)=1._r8
      trcs_VLN(ids_NH4B,L,NY,NX)=0._r8
      trcs_VLN(idg_NH3,L,NY,NX)=trcs_VLN(ids_NH4,L,NY,NX)
      trcs_VLN(idg_NH3B,L,NY,NX)=trcs_VLN(ids_NH4B,L,NY,NX)

      trc_solml(ids_NH4,L,NY,NX)=trc_solml(ids_NH4,L,NY,NX)+trc_solml(ids_NH4B,L,NY,NX)
      trc_solml(idg_NH3,L,NY,NX)=trc_solml(idg_NH3,L,NY,NX)+trc_solml(idg_NH3B,L,NY,NX)
      trc_solml(ids_NH4B,L,NY,NX)=0._r8
      trc_solml(idg_NH3B,L,NY,NX)=0._r8
      trcx_solml(idx_NH4,L,NY,NX)=trcx_solml(idx_NH4,L,NY,NX)+trcx_solml(idx_NH4B,L,NY,NX)
      trcx_solml(idx_NH4B,L,NY,NX)=0._r8
    ENDIF
  ENDIF
  end subroutine UpdateNH3FertilizerBandinfo
!------------------------------------------------------------------------------------------

  subroutine UpdatePO4FertilizerBandinfo(L,NY,NX,FLWD)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: FLWD
  real(r8) :: XVLPO4,DPFLW
  real(r8) :: DWPO4,DZH0P,DZH1P
  real(r8) :: DZH2P,DZH3P,DZF1P
  real(r8) :: DZF2P,DZC0P
  real(r8) :: DZC1P,DZC2P,DZM1P,DXOH0,DXOH1,DXOH2,DXH1P,DXH2P
  real(r8) :: DPALP,DPFEP,DPCDP,DPCHP,DPCMP,FVLPO4
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
      XVLPO4=trcs_VLN(ids_H1PO4,L,NY,NX)
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
        trcs_VLN(ids_H1PO4B,L,NY,NX)=AZMAX1(AMIN1(0.999,WDPOB(L,NY,NX) &
          /ROWP(NY,NX)*DPPOB(L,NY,NX)/DLYR(3,L,NY,NX)))
      ELSE
        trcs_VLN(ids_H1PO4B,L,NY,NX)=0._r8
      ENDIF
      trcs_VLN(ids_H1PO4,L,NY,NX)=1._r8-trcs_VLN(ids_H1PO4B,L,NY,NX)

      trcs_VLN(ids_H2PO4B,L,NY,NX)=trcs_VLN(ids_H1PO4B,L,NY,NX)
      trcs_VLN(ids_H2PO4,L,NY,NX)=trcs_VLN(ids_H1PO4,L,NY,NX)
      FVLPO4=AZMIN1((trcs_VLN(ids_H1PO4,L,NY,NX)-XVLPO4)/XVLPO4)
!
!     TRANSFER HPO4,H2PO4 FROM NON-BAND TO BAND
!     DURING BAND GROWTH DEPENDING ON SALT
!     VS. NON-SALT OPTION
!
!     DZ*,DX*,DP*=transfer of solute,adsorbed,precipitated HPO4,H2PO4
!
      IF(ISALTG.NE.0)THEN
        DZH0P=FVLPO4*trcsa_solml(idsa_H0PO4,L,NY,NX)
        DZH1P=FVLPO4*trc_solml(ids_H1PO4,L,NY,NX)/patomw
        DZH2P=FVLPO4*trc_solml(ids_H2PO4,L,NY,NX)/patomw
        DZH3P=FVLPO4*trcsa_solml(idsa_H3PO4,L,NY,NX)
        DZF1P=FVLPO4*trcsa_solml(idsa_FeHPO4,L,NY,NX)
        DZF2P=FVLPO4*trcsa_solml(idsa_FeH2PO4,L,NY,NX)
        DZC0P=FVLPO4*trcsa_solml(idsa_CaPO4,L,NY,NX)
        DZC1P=FVLPO4*trcsa_solml(idsa_CaHPO4,L,NY,NX)
        DZC2P=FVLPO4*trcsa_solml(idsa_CaH2PO4,L,NY,NX)
        DZM1P=FVLPO4*trcsa_solml(idsa_MgHPO4,L,NY,NX)
        DXOH0=FVLPO4*trcx_solml(idx_OHe,L,NY,NX)
        DXOH1=FVLPO4*trcx_solml(idx_OH,L,NY,NX)
        DXOH2=FVLPO4*trcx_solml(idx_OHp,L,NY,NX)
        DXH1P=FVLPO4*trcx_solml(idx_HPO4,L,NY,NX)
        DXH2P=FVLPO4*trcx_solml(idx_H2PO4,L,NY,NX)
        DPALP=FVLPO4*trcp_salml(idsp_AlPO4,L,NY,NX)
        DPFEP=FVLPO4*trcp_salml(idsp_FePO4,L,NY,NX)
        DPCDP=FVLPO4*trcp_salml(idsp_CaHPO4,L,NY,NX)
        DPCHP=FVLPO4*trcp_salml(idsp_HA,L,NY,NX)
        DPCMP=FVLPO4*trcp_salml(idsp_CaH2PO4,L,NY,NX)

        trcsa_TR(idsa_H0PO4,L,NY,NX)=trcsa_TR(idsa_H0PO4,L,NY,NX)+DZH0P
        TRH1P(L,NY,NX)=TRH1P(L,NY,NX)+DZH1P
        TRH2P(L,NY,NX)=TRH2P(L,NY,NX)+DZH2P
        trcsa_TR(idsa_H3PO4,L,NY,NX)=trcsa_TR(idsa_H3PO4,L,NY,NX)+DZH3P
        trcsa_TR(idsa_FeHPO4,L,NY,NX)=trcsa_TR(idsa_FeHPO4,L,NY,NX)+DZF1P
        trcsa_TR(idsa_FeH2PO4,L,NY,NX)=trcsa_TR(idsa_FeH2PO4,L,NY,NX)+DZF2P
        trcsa_TR(idsa_CaPO4,L,NY,NX)=trcsa_TR(idsa_CaPO4,L,NY,NX)+DZC0P
        trcsa_TR(idsa_CaHPO4,L,NY,NX)=trcsa_TR(idsa_CaHPO4,L,NY,NX)+DZC1P
        trcsa_TR(idsa_CaH2PO4,L,NY,NX)=trcsa_TR(idsa_CaH2PO4,L,NY,NX)+DZC2P
        trcsa_TR(idsa_MgHPO4,L,NY,NX)=trcsa_TR(idsa_MgHPO4,L,NY,NX)+DZM1P
        trcsa_TR(idsa_H0PO4B,L,NY,NX)=trcsa_TR(idsa_H0PO4B,L,NY,NX)-DZH0P
        TRH1B(L,NY,NX)=TRH1B(L,NY,NX)-DZH1P
        TRH2B(L,NY,NX)=TRH2B(L,NY,NX)-DZH2P
        trcsa_TR(idsa_H3PO4B,L,NY,NX)=trcsa_TR(idsa_H3PO4B,L,NY,NX)-DZH3P
        trcsa_TR(idsa_FeHPO4B,L,NY,NX)=trcsa_TR(idsa_FeHPO4B,L,NY,NX)-DZF1P
        trcsa_TR(idsa_FeH2PO4B,L,NY,NX)=trcsa_TR(idsa_FeH2PO4B,L,NY,NX)-DZF2P
        trcsa_TR(idsa_CaPO4B,L,NY,NX)=trcsa_TR(idsa_CaPO4B,L,NY,NX)-DZC0P
        trcsa_TR(idsa_CaHPO4B,L,NY,NX)=trcsa_TR(idsa_CaHPO4B,L,NY,NX)-DZC1P
        trcsa_TR(idsa_CaH2PO4B,L,NY,NX)=trcsa_TR(idsa_CaH2PO4B,L,NY,NX)-DZC2P
        trcsa_TR(idsa_MgHPO4B,L,NY,NX)=trcsa_TR(idsa_MgHPO4B,L,NY,NX)-DZM1P
        trcx_TR(idx_OHe,L,NY,NX)=trcx_TR(idx_OHe,L,NY,NX)+DXOH0
        trcx_TR(idx_OH,L,NY,NX)=trcx_TR(idx_OH,L,NY,NX)+DXOH1
        trcx_TR(idx_OHp,L,NY,NX)=trcx_TR(idx_OHp,L,NY,NX)+DXOH2
        trcx_TR(idx_HPO4,L,NY,NX)=trcx_TR(idx_HPO4,L,NY,NX)+DXH1P
        trcx_TR(idx_H2PO4,L,NY,NX)=trcx_TR(idx_H2PO4,L,NY,NX)+DXH2P
        trcx_TR(idx_OHeB,L,NY,NX)=trcx_TR(idx_OHeB,L,NY,NX)-DXOH0
        trcx_TR(idx_OHB,L,NY,NX)=trcx_TR(idx_OHB,L,NY,NX)-DXOH1
        trcx_TR(idx_OHpB,L,NY,NX)=trcx_TR(idx_OHpB,L,NY,NX)-DXOH2
        trcx_TR(idx_HPO4B,L,NY,NX)=trcx_TR(idx_HPO4B,L,NY,NX)-DXH1P
        trcx_TR(idx_H2PO4B,L,NY,NX)=trcx_TR(idx_H2PO4B,L,NY,NX)-DXH2P
        trcp_TR(idsp_AlPO4,L,NY,NX)=trcp_TR(idsp_AlPO4,L,NY,NX)+DPALP
        trcp_TR(idsp_FePO4,L,NY,NX)=trcp_TR(idsp_FePO4,L,NY,NX)+DPFEP
        trcp_TR(idsp_CaHPO4,L,NY,NX)=trcp_TR(idsp_CaHPO4,L,NY,NX)+DPCDP
        trcp_TR(idsp_HA,L,NY,NX)=trcp_TR(idsp_HA,L,NY,NX)+DPCHP
        trcp_TR(idsp_CaH2PO4,L,NY,NX)=trcp_TR(idsp_CaH2PO4,L,NY,NX)+DPCMP
        trcp_TR(idsp_AlPO4B,L,NY,NX)=trcp_TR(idsp_AlPO4B,L,NY,NX)-DPALP
        trcp_TR(idsp_FePO4B,L,NY,NX)=trcp_TR(idsp_FePO4B,L,NY,NX)-DPFEP
        trcp_TR(idsp_CaHPO4B,L,NY,NX)=trcp_TR(idsp_CaHPO4B,L,NY,NX)-DPCDP
        trcp_TR(idsp_HAB,L,NY,NX)=trcp_TR(idsp_HAB,L,NY,NX)-DPCHP
        trcp_TR(idsp_CaH2PO4B,L,NY,NX)=trcp_TR(idsp_CaH2PO4B,L,NY,NX)-DPCMP
      ELSE
        DZH1P=FVLPO4*trc_solml(ids_H1PO4,L,NY,NX)/patomw
        DZH2P=FVLPO4*trc_solml(ids_H2PO4,L,NY,NX)/patomw
        DXOH1=FVLPO4*trcx_solml(idx_OH,L,NY,NX)
        DXOH2=FVLPO4*trcx_solml(idx_OHp,L,NY,NX)
        DXH2P=FVLPO4*trcx_solml(idx_H2PO4,L,NY,NX)
        DPALP=FVLPO4*trcp_salml(idsp_AlPO4,L,NY,NX)
        DPFEP=FVLPO4*trcp_salml(idsp_FePO4,L,NY,NX)
        DPCDP=FVLPO4*trcp_salml(idsp_CaHPO4,L,NY,NX)
        DPCHP=FVLPO4*trcp_salml(idsp_HA,L,NY,NX)
        DPCMP=FVLPO4*trcp_salml(idsp_CaH2PO4,L,NY,NX)
        TRH1P(L,NY,NX)=TRH1P(L,NY,NX)+DZH1P
        TRH2P(L,NY,NX)=TRH2P(L,NY,NX)+DZH2P
        trcx_TR(idx_OH,L,NY,NX)=trcx_TR(idx_OH,L,NY,NX)+DXOH1
        trcx_TR(idx_OHp,L,NY,NX)=trcx_TR(idx_OHp,L,NY,NX)+DXOH2
        trcx_TR(idx_H2PO4,L,NY,NX)=trcx_TR(idx_H2PO4,L,NY,NX)+DXH2P
        TRH1B(L,NY,NX)=TRH1B(L,NY,NX)-DZH1P
        TRH2B(L,NY,NX)=TRH2B(L,NY,NX)-DZH2P
        trcx_TR(idx_OHB,L,NY,NX)=trcx_TR(idx_OHB,L,NY,NX)-DXOH1
        trcx_TR(idx_OHpB,L,NY,NX)=trcx_TR(idx_OHpB,L,NY,NX)-DXOH2
        trcx_TR(idx_H2PO4B,L,NY,NX)=trcx_TR(idx_H2PO4B,L,NY,NX)-DXH2P
        trcp_TR(idsp_AlPO4,L,NY,NX)=trcp_TR(idsp_AlPO4,L,NY,NX)+DPALP
        trcp_TR(idsp_FePO4,L,NY,NX)=trcp_TR(idsp_FePO4,L,NY,NX)+DPFEP
        trcp_TR(idsp_CaHPO4,L,NY,NX)=trcp_TR(idsp_CaHPO4,L,NY,NX)+DPCDP
        trcp_TR(idsp_HA,L,NY,NX)=trcp_TR(idsp_HA,L,NY,NX)+DPCHP
        trcp_TR(idsp_CaH2PO4,L,NY,NX)=trcp_TR(idsp_CaH2PO4,L,NY,NX)+DPCMP
        trcp_TR(idsp_AlPO4B,L,NY,NX)=trcp_TR(idsp_AlPO4B,L,NY,NX)-DPALP
        trcp_TR(idsp_FePO4B,L,NY,NX)=trcp_TR(idsp_FePO4B,L,NY,NX)-DPFEP
        trcp_TR(idsp_CaHPO4B,L,NY,NX)=trcp_TR(idsp_CaHPO4B,L,NY,NX)-DPCDP
        trcp_TR(idsp_HAB,L,NY,NX)=trcp_TR(idsp_HAB,L,NY,NX)-DPCHP
        trcp_TR(idsp_CaH2PO4B,L,NY,NX)=trcp_TR(idsp_CaH2PO4B,L,NY,NX)-DPCMP
      ENDIF
    ELSE
!
!     AMALGAMATE PO4 BAND WITH NON-BAND IF BAND NO LONGER EXISTS
!
      DPPOB(L,NY,NX)=0._r8
      WDPOB(L,NY,NX)=0._r8
      trcs_VLN(ids_H1PO4B,L,NY,NX)=0._r8
      trcs_VLN(ids_H1PO4,L,NY,NX)=1._r8
      trcs_VLN(ids_H2PO4B,L,NY,NX)=trcs_VLN(ids_H1PO4B,L,NY,NX)
      trcs_VLN(ids_H2PO4,L,NY,NX)=trcs_VLN(ids_H1PO4,L,NY,NX)

      trc_solml(ids_H1PO4,L,NY,NX)=trc_solml(ids_H1PO4,L,NY,NX)+trc_solml(ids_H1PO4B,L,NY,NX)
      trc_solml(ids_H2PO4,L,NY,NX)=trc_solml(ids_H2PO4,L,NY,NX)+trc_solml(ids_H2PO4B,L,NY,NX)
      trc_solml(ids_H1PO4B,L,NY,NX)=0._r8
      trc_solml(ids_H2PO4B,L,NY,NX)=0._r8

      trcsa_solml(idsa_H0PO4,L,NY,NX)=trcsa_solml(idsa_H0PO4,L,NY,NX)+trcsa_solml(idsa_H0PO4B,L,NY,NX)
      trcsa_solml(idsa_H3PO4,L,NY,NX)=trcsa_solml(idsa_H3PO4,L,NY,NX)+trcsa_solml(idsa_H3PO4B,L,NY,NX)
      trcsa_solml(idsa_FeHPO4,L,NY,NX)=trcsa_solml(idsa_FeHPO4,L,NY,NX)+trcsa_solml(idsa_FeHPO4B,L,NY,NX)
      trcsa_solml(idsa_FeH2PO4,L,NY,NX)=trcsa_solml(idsa_FeH2PO4,L,NY,NX)+trcsa_solml(idsa_FeH2PO4B,L,NY,NX)
      trcsa_solml(idsa_CaPO4,L,NY,NX)=trcsa_solml(idsa_CaPO4,L,NY,NX)+trcsa_solml(idsa_CaPO4B,L,NY,NX)
      trcsa_solml(idsa_CaHPO4,L,NY,NX)=trcsa_solml(idsa_CaHPO4,L,NY,NX)+trcsa_solml(idsa_CaHPO4B,L,NY,NX)
      trcsa_solml(idsa_CaH2PO4,L,NY,NX)=trcsa_solml(idsa_CaH2PO4,L,NY,NX)+trcsa_solml(idsa_CaH2PO4B,L,NY,NX)
      trcsa_solml(idsa_MgHPO4,L,NY,NX)=trcsa_solml(idsa_MgHPO4,L,NY,NX)+trcsa_solml(idsa_MgHPO4B,L,NY,NX)
      trcsa_solml(idsa_H0PO4B,L,NY,NX)=0._r8
      trcsa_solml(idsa_H3PO4B,L,NY,NX)=0._r8
      trcsa_solml(idsa_FeHPO4B,L,NY,NX)=0._r8
      trcsa_solml(idsa_FeH2PO4B,L,NY,NX)=0._r8
      trcsa_solml(idsa_CaPO4B,L,NY,NX)=0._r8
      trcsa_solml(idsa_CaHPO4B,L,NY,NX)=0._r8
      trcsa_solml(idsa_CaH2PO4B,L,NY,NX)=0._r8
      trcsa_solml(idsa_MgHPO4B,L,NY,NX)=0._r8

      trcx_solml(idx_OHe,L,NY,NX)=trcx_solml(idx_OHe,L,NY,NX)+trcx_solml(idx_OHeB,L,NY,NX)
      trcx_solml(idx_OH,L,NY,NX)=trcx_solml(idx_OH,L,NY,NX)+trcx_solml(idx_OHB,L,NY,NX)
      trcx_solml(idx_OHp,L,NY,NX)=trcx_solml(idx_OHp,L,NY,NX)+trcx_solml(idx_OHpB,L,NY,NX)
      trcx_solml(idx_HPO4,L,NY,NX)=trcx_solml(idx_HPO4,L,NY,NX)+trcx_solml(idx_HPO4B,L,NY,NX)
      trcx_solml(idx_H2PO4,L,NY,NX)=trcx_solml(idx_H2PO4,L,NY,NX)+trcx_solml(idx_H2PO4B,L,NY,NX)
      trcx_solml(idx_OHeB,L,NY,NX)=0._r8
      trcx_solml(idx_OHB,L,NY,NX)=0._r8
      trcx_solml(idx_OHpB,L,NY,NX)=0._r8
      trcx_solml(idx_HPO4B,L,NY,NX)=0._r8
      trcx_solml(idx_H2PO4B,L,NY,NX)=0._r8

      trcp_salml(idsp_AlPO4,L,NY,NX)=trcp_salml(idsp_AlPO4,L,NY,NX)+trcp_salml(idsp_AlPO4B,L,NY,NX)
      trcp_salml(idsp_FePO4,L,NY,NX)=trcp_salml(idsp_FePO4,L,NY,NX)+trcp_salml(idsp_FePO4B,L,NY,NX)
      trcp_salml(idsp_CaHPO4,L,NY,NX)=trcp_salml(idsp_CaHPO4,L,NY,NX)+trcp_salml(idsp_CaHPO4B,L,NY,NX)
      trcp_salml(idsp_HA,L,NY,NX)=trcp_salml(idsp_HA,L,NY,NX)+trcp_salml(idsp_HAB,L,NY,NX)
      trcp_salml(idsp_CaH2PO4,L,NY,NX)=trcp_salml(idsp_CaH2PO4,L,NY,NX)+trcp_salml(idsp_CaH2PO4B,L,NY,NX)
      trcp_salml(idsp_AlPO4B,L,NY,NX)=0._r8
      trcp_salml(idsp_FePO4B,L,NY,NX)=0._r8
      trcp_salml(idsp_CaHPO4B,L,NY,NX)=0._r8
      trcp_salml(idsp_HAB,L,NY,NX)=0._r8
      trcp_salml(idsp_CaH2PO4B,L,NY,NX)=0._r8
    ENDIF
  ENDIF
  end subroutine UpdatePO4FertilizerBandinfo
!------------------------------------------------------------------------------------------

  subroutine UpdateNO3FertilizerBandinfo(L,NY,NX,FLWD)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: FLWD
  real(r8) :: XVLNO3,DWNO3,DNO3S
  real(r8) :: DNO2S,FVLNO3,DPFLW

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
      DWNO3=0.5*SQRT(SolDifc(ids_NO3,L,NY,NX))*TORT(NPH,L,NY,NX)
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
      XVLNO3=trcs_VLN(ids_NO3,L,NY,NX)
      IF(DLYR(3,L,NY,NX).GT.ZERO)THEN
        trcs_VLN(ids_NO3B,L,NY,NX)=AZMAX1(AMIN1(0.999,WDNOB(L,NY,NX) &
          /ROWO(NY,NX)*DPNOB(L,NY,NX)/DLYR(3,L,NY,NX)))
      ELSE
        trcs_VLN(ids_NO3B,L,NY,NX)=0._r8
      ENDIF
      trcs_VLN(ids_NO3,L,NY,NX)=1._r8-trcs_VLN(ids_NO3B,L,NY,NX)

      trcs_VLN(ids_NO2B,L,NY,NX)=trcs_VLN(ids_NO3B,L,NY,NX)
      trcs_VLN(ids_NO2,L,NY,NX)=trcs_VLN(ids_NO3,L,NY,NX)
      FVLNO3=AZMIN1((trcs_VLN(ids_NO3,L,NY,NX)-XVLNO3)/XVLNO3)
!
!     TRANSFER NO3 FROM NON-BAND TO BAND
!     DURING BAND GROWTH
!
!     DNO3S,DNO2S=transfer of NO3,NO2
!
      DNO3S=FVLNO3*trc_solml(ids_NO3,L,NY,NX)/natomw
      DNO2S=FVLNO3*trc_solml(ids_NO2,L,NY,NX)/natomw
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
      trcs_VLN(ids_NO3,L,NY,NX)=1._r8
      trcs_VLN(ids_NO3B,L,NY,NX)=0._r8
      trcs_VLN(ids_NO2,L,NY,NX)=trcs_VLN(ids_NO3,L,NY,NX)
      trcs_VLN(ids_NO2B,L,NY,NX)=trcs_VLN(ids_NO3B,L,NY,NX)

      trc_solml(ids_NO3,L,NY,NX)=trc_solml(ids_NO3,L,NY,NX)+trc_solml(ids_NO3B,L,NY,NX)
      trc_solml(ids_NO2,L,NY,NX)=trc_solml(ids_NO2,L,NY,NX)+trc_solml(ids_NO2B,L,NY,NX)
      trc_solml(ids_NO3B,L,NY,NX)=0._r8
      trc_solml(ids_NO2B,L,NY,NX)=0._r8
    ENDIF
  ENDIF
  end subroutine UpdateNO3FertilizerBandinfo
!------------------------------------------------------------------------------------------

  subroutine UpdateSoluteinSurfaceResidue(NX,NY)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: CAL1,CALX,CCA1,CCAX
  real(r8) :: CCEC,CCEC0,CCO31
  real(r8) :: CFE1,CFEX
  real(r8) :: CH2PA,CH2PD,CH2PF,CH2PH,CH2PM
  real(r8) :: CHY1,COH1,DP
  real(r8) :: FX,S0,S1,XALQ,XCAQ,XCAX
  real(r8) :: XFEQ,XHYQ,XN4Q,XTLQ
  real(r8) :: VOLWMX,VOLWMP,BKVLX
  real(r8) :: THETWR,COMA
  real(r8) :: CNHUA
  real(r8) :: CCO20
  real(r8) :: RH1PX,RH2PX,DUKD,DFNSA,RN3X,RN4X
  real(r8) :: CH1P1,CH2P1,CN31,CN41,PALPO1,PCAPD1
  real(r8) :: PCAPH1,PCAPM1,PFEPO1,RH2P,RNH4,RPALPX
  real(r8) :: RPCADX,RPCAHX,RPCAMX,RPFEPX,RXN4,XN41
  real(r8) :: RHP1,RHP2,RN3S,RN4S
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
      COMA=AMIN1(0.1E+06_r8,TOQCK(0,NY,NX)/VOLQ(0,NY,NX))
    ELSE
      COMA=0.1E+06_r8
    ENDIF
    DUKD=DUKM*(1.0_r8+COMA/DUKI)
!
!     UREA HYDROLYSIS INHIBITION
!
!     ZNHU0,ZNHUI=initial,current inhibition activity
!     RNHUI=rate constant for decline in urea hydrolysis inhibition
!
    IF(ZNHU0(0,NY,NX).GT.ZEROS(NY,NX).AND.ZNHUI(0,NY,NX).GT.ZEROS(NY,NX))THEN
      ZNHUI(0,NY,NX)=ZNHUI(0,NY,NX)-TFNQ(0,NY,NX)**0.25_r8 &
        *RNHUI(IUTYP(NY,NX))*ZNHUI(0,NY,NX) &
        *AMAX1(RNHUI(IUTYP(NY,NX)),1.0_r8-ZNHUI(0,NY,NX)/ZNHU0(0,NY,NX))
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
    IF(FertN_soil(ifert_urea,0,NY,NX).GT.ZEROS(NY,NX) &
      .AND.BKVL(0,NY,NX).GT.ZEROS(NY,NX))THEN
      CNHUA=FertN_soil(ifert_urea,0,NY,NX)/BKVL(0,NY,NX)
      DFNSA=CNHUA/(CNHUA+DUKD)
      RSNUA=AMIN1(FertN_soil(ifert_urea,0,NY,NX),SPNHU*TOQCK(0,NY,NX)*DFNSA*TFNQ(0,NY,NX))*(1.0-ZNHUI(0,NY,NX))
    ELSE
      RSNUA=0._r8
    ENDIF
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
    RSN4AA=SPNH4*FertN_soil(ifert_nh4,0,NY,NX)*THETWR
    RSN3AA=SPNH3*FertN_soil(ifert_nh3,0,NY,NX)
    RSNUAA=RSNUA*THETWR
    RSNOAA=SPNO3*FertN_soil(ifert_no3,0,NY,NX)*THETWR
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
      VOLWMX=natomw*VOLWM(NPH,0,NY,NX)
      RN4X=(XNH4S(0,NY,NX)+natomw*RSN4AA)/VOLWMX
      RN3X=natomw*RSNUAA/VOLWMX
      CN41=AMAX1(ZERO,trc_solml(ids_NH4,0,NY,NX)/VOLWMX+RN4X)
      CN31=AMAX1(ZERO,trc_solml(idg_NH3,0,NY,NX)/VOLWMX+RN3X)
      IF(BKVLX.GT.ZEROS(NY,NX))THEN
        XN41=AMAX1(ZERO,trcx_solml(idx_NH4,0,NY,NX)/BKVLX)
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
      VOLWMP=patomw*VOLWM(NPH,0,NY,NX)
      RH1PX=XH1PS(0,NY,NX)/VOLWMP
      RH2PX=XH2PS(0,NY,NX)/VOLWMP
      CH1P1=AZMAX1(trc_solml(ids_H1PO4,0,NY,NX)/VOLWMP+RH1PX)
      CH2P1=AZMAX1(trc_solml(ids_H2PO4,0,NY,NX)/VOLWMP+RH2PX)
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
      PALPO1=AZMAX1(trcp_salml(idsp_AlPO4,0,NY,NX)/BKVLX)
      PFEPO1=AZMAX1(trcp_salml(idsp_FePO4,0,NY,NX)/BKVLX)
      PCAPM1=AZMAX1(trcp_salml(idsp_CaH2PO4,0,NY,NX)/BKVLX)
      PCAPD1=AZMAX1(trcp_salml(idsp_CaHPO4,0,NY,NX)/BKVLX)
      PCAPH1=AZMAX1(trcp_salml(idsp_HA,0,NY,NX)/BKVLX)
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
    CHY1=AMAX1(ZERO,10.0_r8**(-(PH(0,NY,NX)-3.0_r8)))
    COH1=AMAX1(ZERO,DPH2O/CHY1)
    CAL1=AMAX1(ZERO,SPALO/COH1**3._r8)
    CFE1=AMAX1(ZERO,SPFEO/COH1**3_r8)
    CCO20=AMAX1(ZERO,trc_solcl(idg_CO2,0,NY,NX)/catomw)
    CCO31=AMAX1(ZERO,CCO20*DPCO3/CHY1**2._r8)
    CCA1=AMAX1(ZERO,AMIN1(CCAMX,SPCAC/CCO31))
!
!     ALUMINUM PHOSPHATE (VARISCITE)
!
!     CH2PA,CH2P1=equilibrium,current H2PO4 concentration in litter
!     SYA0P2=solubility product derived from SPALO
!     RPALPX=H2PO4 dissolution from AlPO4 in litter
!
    CH2PA=SYA0P2/(CAL1*COH1**2)
    RPALPX=AMIN1(AZMAX1(4.0E-08*ORGC(0,NY,NX)-PALPO1),AMAX1(-PALPO1,TPD*(CH2P1-CH2PA)))
!
!     IRON PHOSPHATE (STRENGITE)
!
!     CH2PF,CH2P1=equilibrium,current H2PO4 concentration in litter
!     SYF0P2=solubility product derived from SPALO
!     RPFEPX=H2PO4 dissolution from FePO4 in litter
!
    CH2PF=SYF0P2/(CFE1*COH1**2)
    RPFEPX=AMIN1(AZMAX1(2.0E-06*ORGC(0,NY,NX)-PFEPO1),AMAX1(-PFEPO1,TPD*(CH2P1-CH2PF)))
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
    RPCAHX=AMIN1(AZMAX1(5.0E-05*ORGC(0,NY,NX)-PCAPH1),AMAX1(-PCAPH1,TPD*(CH2P1-CH2PH)))
!
!     MONOCALCIUM PHOSPHATE
!
!     CH2PM,CH2P1=equilibrium,current H2PO4 concentration in litter
!     SPCAM=solubility product for Ca(H2PO4)2
!     RPCAMX=H2PO4 dissolution from Ca(H2PO4)2 in litter
!
    CH2PM=SQRT(SPCAM/CCA1)
    RPCAMX=AMIN1(AZMAX1(5.0E-05_r8*ORGC(0,NY,NX)-PCAPM1),AMAX1(-PCAPM1*SPPO4,TPD*(CH2P1-CH2PM)))

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
    S1=AZMAX1(S0**2-4.0*(CH1P1*CHY1-DP*CH2P1))
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
      trcx_TR(idx_NH4,0,NY,NX)=trcx_TR(idx_NH4,0,NY,NX)+RXN4*VOLWM(NPH,0,NY,NX)
      trcp_TR(idsp_AlPO4,0,NY,NX)=trcp_TR(idsp_AlPO4,0,NY,NX)+RPALPX*VOLWM(NPH,0,NY,NX)
      trcp_TR(idsp_FePO4,0,NY,NX)=trcp_TR(idsp_FePO4,0,NY,NX)+RPFEPX*VOLWM(NPH,0,NY,NX)
      trcp_TR(idsp_CaHPO4,0,NY,NX)=trcp_TR(idsp_CaHPO4,0,NY,NX)+RPCADX*VOLWM(NPH,0,NY,NX)
      trcp_TR(idsp_HA,0,NY,NX)=trcp_TR(idsp_HA,0,NY,NX)+RPCAHX*VOLWM(NPH,0,NY,NX)
      trcp_TR(idsp_CaH2PO4,0,NY,NX)=trcp_TR(idsp_CaH2PO4,0,NY,NX)+RPCAMX*VOLWM(NPH,0,NY,NX)
      FertN_soil(ifert_nh4,0,NY,NX)=FertN_soil(ifert_nh4,0,NY,NX)-RSN4AA
      FertN_soil(ifert_nh3,0,NY,NX)=FertN_soil(ifert_nh3,0,NY,NX)-RSN3AA
      FertN_soil(ifert_urea,0,NY,NX)=FertN_soil(ifert_urea,0,NY,NX)-RSNUAA
      FertN_soil(ifert_no3,0,NY,NX)=FertN_soil(ifert_no3,0,NY,NX)-RSNOAA
      TRN4S(0,NY,NX)=TRN4S(0,NY,NX)+RSN4AA
      TRN3S(0,NY,NX)=TRN3S(0,NY,NX)+RSN3AA+RSNUAA
      TRNO3(0,NY,NX)=TRNO3(0,NY,NX)+RSNOAA
      TRN4S(0,NY,NX)=TRN4S(0,NY,NX)*natomw
      TRN3S(0,NY,NX)=TRN3S(0,NY,NX)*natomw
      TRNO3(0,NY,NX)=TRNO3(0,NY,NX)*natomw
      TRH1P(0,NY,NX)=TRH1P(0,NY,NX)*patomw
      TRH2P(0,NY,NX)=TRH2P(0,NY,NX)*patomw
!     WRITE(*,9989)'TRN4S',I,J,TRN4S(0,NY,NX)
!    2,RN4S,RNH4,RXN4,RSN4AA,VOLWM(NPH,0,NY,NX)
!    3,SPNH4,FertN_soil(ifert_nh4,0,NY,NX),THETWR
!9989  FORMAT(A8,2I4,12E12.4)
  end subroutine UpdateSoluteinSurfaceResidue

END module SoluteMod
