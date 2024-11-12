module SnowTransportMod
  use data_kind_mod,  only: r8 => DAT_KIND_R8
  use data_const_mod, only: spval => DAT_CONST_SPVAL
  use EcoSimConst,    only: DENSICE
  use EcoSIMCtrlMod,  only: lverb
  use MiniMathMod,    only: fixEXflux
  use SnowPhysData
  use GridDataType
  use SnowDataType
  use SoilBGCDataType
  use AqueChemDatatype
  use SoilWaterDataType
  use SOMDataType
  use ChemTranspDataType
  use EcoSimSumDataType
implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: VerticalSaltFluxThruSnowpack  
  public :: MassFluxThruSnowRunoff
  public :: DiagSnowChemMass
  public :: OverlandFlowThruSnow
  contains
!------------------------------------------------------------------------------------------

  subroutine VerticalSaltFluxThruSnowpack(I,J,N1,N2,NY,NX)
  !
  !update solute due to snow flux
  implicit none
  integer, intent(in) :: N1,N2,NY,NX,I,J
  integer :: LS, LS2
  integer :: NTG,NTN,NTSA,NTS
!     begin_execution
!     NET WATER AND HEAT FLUXES THROUGH SNOWPACK
!
!     VHCPW,VLHeatCapSnowMin_col=current, minimum snowpack heat capacities
!     CumSno2SnowLay,CumWat2SnowLay,CumIce2SnowLay=net fluxes of snow,water,ice in snowpack
!     CumHeat2SnowLay=convective heat fluxes of snow,water,ice in snowpack
!     XFLWS,WatXfer2SnoLay,IceXfer2SnoLay=snow,water,ice transfer from watsub.f
!     HeatXfer2SnoLay=convective heat flux from snow,water,ice transfer from watsub.f
!     WatConvSno2MicP,WatConvSno2MacP,WatConvSno2LitR=water flux from lowest snow layer to soil macropore,micropore,litter
!     HeatConvSno2Soi,HeatConvSno2LitR=heat flux from lowest snow layer to soil,litter
  
  D1205: DO LS=1,JS

    IF(VLHeatCapSnow_snvr(LS,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
      !id of next snow layer
      LS2=MIN(JS,LS+1)
!
!     IF LOWER LAYER IS IN THE SNOWPACK
!
      IF(LS.LT.JS .AND. VLHeatCapSnow_snvr(LS2,N2,N1).GT.VLHeatCapSnowMin_col(N2,N1))THEN
        !not surface layer, and is heat significant
        !     NET SOLUTE FLUXES THROUGH SNOWPACK
        !
        !     T*BLS=net solute flux in snowpack
        !     X*BLS=solute flux in snowpack from TranspNoSalt.f
        !     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
        !             :OC=DOC,ON=DON,OP=DOP,OA=acetate
        !             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
        !             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
        !
        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1)+trcg_Xbndl_flx(NTG,LS,N2,N1) &
            -trcg_Xbndl_flx(NTG,LS2,N2,N1)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_TBLS(NTN,LS,N2,N1)=trcn_TBLS(NTN,LS,N2,N1)+trcn_Xbndl_flx(NTN,LS,N2,N1) &
            -trcn_Xbndl_flx(NTN,LS2,N2,N1)
        ENDDO
        !
        !     NET SALT FLUXES THROUGH SNOWPACK
        !
        !     T*BLS=net solute flux in snowpack
        !     X*BLS=solute flux in snowpack from TranspSalt.f
        !     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
        !          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
        !          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
        !          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
        !          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
        !          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
        !          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
        !     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
        !          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
        !          :*1=non-band,*B=band
        !
        IF(salt_model)THEN
          DO NTSA=idsalt_beg,idsalt_end
            trcSalt_TBLS(NTSA,LS,N2,N1)=trcSalt_TBLS(NTSA,LS,N2,N1)+trcSaltFlo2SnowLay(NTSA,LS,N2,N1)&
              -trcSaltFlo2SnowLay(NTSA,LS2,N2,N1)
          ENDDO
        ENDIF
        !
        !     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
        !
      ELSE

        ! and NH3B
        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1)+trcg_Xbndl_flx(NTG,LS,N2,N1) &
            -trcs_Transp2MicP_3D(NTG,3,0,N2,N1)-trcs_Transp2MicP_3D(NTG,3,NUM(N2,N1),N2,N1) &
            -trcs_Transp2MacP_3D(NTG,3,NUM(N2,N1),N2,N1)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_TBLS(NTN,LS,N2,N1)=trcn_TBLS(NTN,LS,N2,N1)+trcn_Xbndl_flx(NTN,LS,N2,N1) &
            -trcs_Transp2MicP_3D(NTN,3,0,N2,N1)-trcs_Transp2MicP_3D(NTN,3,NUM(N2,N1),N2,N1) &
            -trcs_Transp2MacP_3D(NTN,3,NUM(N2,N1),N2,N1)
        ENDDO

        !add band flux
        trcg_TBLS(idg_NH3,LS,N2,N1)=trcg_TBLS(idg_NH3,LS,N2,N1) &
          -trcs_Transp2MicP_3D(idg_NH3B,3,NUM(N2,N1),N2,N1)-trcs_Transp2MacP_3D(idg_NH3B,3,NUM(N2,N1),N2,N1)

        DO NTS=0,ids_nuts
          trcn_TBLS(ids_NH4+NTS,LS,N2,N1)=trcn_TBLS(ids_NH4+NTS,LS,N2,N1) &
            -trcs_Transp2MicP_3D(ids_NH4B+NTS,3,NUM(N2,N1),N2,N1)-trcs_Transp2MacP_3D(ids_NH4B+NTS,3,NUM(N2,N1),N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO NTSA=idsalt_beg,idsalt_end
            trcSalt_TBLS(NTSA,LS,NY,NX)=trcSalt_TBLS(NTSA,LS,NY,NX)+trcSaltFlo2SnowLay(NTSA,LS,NY,NX) &
              -trcSalt3DFlo2Cell(NTSA,3,0,N2,N1)-trcSalt3DFlo2Cell(NTSA,3,NUM(N2,N1),N2,N1) &
              -trcSalt_XFHS(NTSA,3,NUM(N2,N1),N2,N1)
          ENDDO

          !add band flux
          DO NTSA=0,idsalt_nuts
            trcSalt_TBLS(idsalt_H0PO4+NTSA,LS,NY,NX)=trcSalt_TBLS(idsalt_H0PO4+NTSA,LS,NY,NX) &
              -trcSalt3DFlo2Cell(idsalt_H0PO4B+NTSA,3,NUM(N2,N1),N2,N1) &
              -trcSalt_XFHS(idsalt_H0PO4B+NTSA,3,NUM(N2,N1),N2,N1)
          ENDDO
        ENDIF
      ENDIF
!
!     WATER,GAS,SOLUTE,SALT FLUXES INTO SNOWPACK SURFACE
!
    ELSEIF(LS.EQ.1)THEN
      IF(abs(SnoXfer2SnoLay_snvr(LS,N2,N1))>0._r8)THEN

        DO NTG=idg_beg,idg_end-1
          trcg_TBLS(NTG,LS,N2,N1)=trcg_TBLS(NTG,LS,N2,N1)+trcg_Xbndl_flx(NTG,LS,N2,N1)
        ENDDO

        DO NTN=ids_nut_beg,ids_nuts_end
          trcn_TBLS(NTN,LS,N2,N1)=trcn_TBLS(NTN,LS,N2,N1)+trcn_Xbndl_flx(NTN,LS,N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO NTSA=idsalt_beg,idsalt_end
            trcSalt_TBLS(NTSA,LS,N2,N1)=trcSalt_TBLS(NTSA,LS,N2,N1)+trcSaltFlo2SnowLay(NTSA,LS,N2,N1)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

  ENDDO D1205
  end subroutine VerticalSaltFluxThruSnowpack

!------------------------------------------------------------------------------------------

  subroutine DiagSnowChemMass(I,J,NY,NX,ENGYSNW)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8),intent(inout):: ENGYSNW   !column heat storage
  real(r8) :: SSW,ENGYW,WS
  integer :: L,nsalts

  if(lverb)write(*,*)'DiagSnowChemMass'
  CALL ChemicalBySnowRedistribution(I,J,NY,NX)
  !
  !     UPDATE STATE VARIABLES WITH TOTAL FLUXES CALCULATED ABOVE
  !
  !     IF(J.EQ.24)THEN
  !
  !     SNOWPACK VARIABLES NEEDED FOR WATER, C, N, P, O, SOLUTE AND
  !     ENERGY BALANCES INCLUDING SUM OF ALL CURRENT STATE VARIABLES,
  !     CUMULATIVE SUMS OF ALL ADDITIONS AND REMOVALS
  !
  DO  L=1,JS
    WS                  = VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE
    if(WS>0._r8)then
      WatMassStore_lnd    = WatMassStore_lnd+WS
      WatMass_col(NY,NX)  = WatMass_col(NY,NX)+WS
      ENGYW               = VLHeatCapSnow_snvr(L,NY,NX)*TKSnow_snvr(L,NY,NX)
      ENGYSNW             = ENGYSNW+ENGYW
      HeatStore_lnd       = HeatStore_lnd+ENGYW
      TGasC_lnd           = TGasC_lnd+trcg_solsml_snvr(idg_CO2,L,NY,NX)+trcg_solsml_snvr(idg_CH4,L,NY,NX)
      DIC_mass_col(NY,NX) = DIC_mass_col(NY,NX)+trcg_solsml_snvr(idg_CO2,L,NY,NX)+trcg_solsml_snvr(idg_CH4,L,NY,NX)
      TSoilO2G_lnd        = TSoilO2G_lnd+trcg_solsml_snvr(idg_O2,L,NY,NX)
      TGasN_lnd           = TGasN_lnd+trcg_solsml_snvr(idg_N2,L,NY,NX)+trcg_solsml_snvr(idg_N2O,L,NY,NX)
      TDisolNH4_lnd       = TDisolNH4_lnd+trcn_solsml(ids_NH4,L,NY,NX)+trcg_solsml_snvr(idg_NH3,L,NY,NX)
      tNO3_lnd            = tNO3_lnd+trcn_solsml(ids_NO3,L,NY,NX)
      TDisolPi_lnd        = TDisolPi_lnd+trcn_solsml(ids_H1PO4,L,NY,NX)+trcn_solsml(ids_H2PO4,L,NY,NX)

      IF(salt_model)THEN
        SSW=0._r8
        do nsalts=idsalt_beg,idsalt_end
          SSW=SSW+trcs_solsml(nsalts,L,NY,NX)*trcSaltIonNumber(nsalts)
        ENDDO  
        TION=TION+SSW
      ENDIF
    ENDIF
  ENDDO
  end subroutine DiagSnowChemMass  
!------------------------------------------------------------------------------------------

  subroutine ChemicalBySnowRedistribution(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX

  integer :: NTA,NTG,NTS
!     begin_execution
!     OVERLAND SNOW REDISTRIBUTION
!
  IF(abs(TDrysnoBySnowRedist(NY,NX))>0._r8)THEN
    DO NTG=idg_beg,idg_end-1
      trcg_solsml_snvr(NTG,1,NY,NX)=trcg_solsml_snvr(NTG,1,NY,NX)+trcg_QSS(NTG,NY,NX)
    ENDDO

    !   write(116,*)I+J/24.,'beg611',trc_solml_vr(ids_NO3,0,NY,NX),trcn_QSS(ids_NO3,NY,NX),trcn_QSS(ids_NO2,NY,NX)

    DO NTS=ids_nut_beg,ids_nuts_end
      trcn_solsml(NTS,1,NY,NX)=trcn_solsml(NTS,1,NY,NX)+trcn_QSS(NTS,NY,NX)
    ENDDO

    !   write(116,*)I+J/24.,'beg612',trc_solml_vr(ids_NO3,0,NY,NX),RNutMicbTransf_vr(ids_NO3,0,NY,NX),trcn_RChem_soil_vr(ids_NO3,0,NY,NX)

    IF(salt_model)THEN
      DO NTA=idsalt_beg,idsalt_end
        trcs_solsml(NTA,1,NY,NX)=trcs_solsml(NTA,1,NY,NX)+trcSalt_TQS(NTA,NY,NX)
      ENDDO
    ENDIF
  ENDIF

  end subroutine ChemicalBySnowRedistribution

!------------------------------------------------------------------------------------------
  subroutine MassFluxThruSnowRunoff(N,N1,N2,N4,N5,N4B,N5B)

  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B
  integer :: NN,NTG,NTN,NTSA

  TDrysnoBySnowRedist(N2,N1)   = TDrysnoBySnowRedist(N2,N1)+DrysnoBySnowRedistrib(N,N2,N1)&
    -DrysnoBySnowRedistrib(N,N5,N4)
  TWatBySnowRedist(N2,N1)      = TWatBySnowRedist(N2,N1)+WatBySnowRedistrib_2DH(N,N2,N1) &
    -WatBySnowRedistrib_2DH(N,N5,N4)
  TIceBySnowRedist(N2,N1)      = TIceBySnowRedist(N2,N1)+IceBySnowRedistrib_2DH(N,N2,N1) &
    -IceBySnowRedistrib_2DH(N,N5,N4)
  THeatBySnowRedist_col(N2,N1) = THeatBySnowRedist_col(N2,N1)+HeatBySnowRedistrib_2DH(N,N2,N1) &
    -HeatBySnowRedistrib_2DH(N,N5,N4)

  D1202: DO NN=1,2
    !gaseous tracers
    DO NTG=idg_beg,idg_end-1
      trcg_TFloXSurRunoff(NTG,N2,N1)=trcg_TFloXSurRunoff(NTG,N2,N1)+trcg_FloXSurRunoff_2D(NTG,N,NN,N2,N1)
    ENDDO

    !nutrient tracres
    DO NTN=ids_nut_beg,ids_nuts_end
      trcn_TFloXSurRunoff_2D(NTN,N2,N1)=trcn_TFloXSurRunoff_2D(NTN,N2,N1)+trcn_FloXSurRunoff_2D(NTN,N,NN,N2,N1)
    ENDDO

    IF(IFLBH(N,NN,N5,N4).EQ.0)THEN    

      DO NTG=idg_beg,idg_end-1
        trcg_TFloXSurRunoff(NTG,N2,N1)=trcg_TFloXSurRunoff(NTG,N2,N1)-trcg_FloXSurRunoff_2D(NTG,N,NN,N5,N4)
      ENDDO
      DO NTN=ids_nut_beg,ids_nuts_end
        trcn_TFloXSurRunoff_2D(NTN,N2,N1)=trcn_TFloXSurRunoff_2D(NTN,N2,N1)-trcn_FloXSurRunoff_2D(NTN,N,NN,N5,N4)
      ENDDO

    ENDIF 
    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      DO NTG=idg_beg,idg_end-1
        trcg_TFloXSurRunoff(NTG,N2,N1)=trcg_TFloXSurRunoff(NTG,N2,N1)-trcg_FloXSurRunoff_2D(NTG,N,NN,N5B,N4B)
      ENDDO
      DO NTN=ids_nut_beg,ids_nuts_end
        trcn_TFloXSurRunoff_2D(NTN,N2,N1)=trcn_TFloXSurRunoff_2D(NTN,N2,N1)-trcn_FloXSurRunoff_2D(NTN,N,NN,N5B,N4B)
      ENDDO

    ENDIF
  ENDDO D1202

  !
  !     NET GAS AND SOLUTE FLUXES FROM RUNOFF AND SNOWPACK
  !
  !     T*QRS=net overland solute flux from runoff
  !     X*QRS=solute in runoff from TranspNoSalt.f
  !     T*QSS=net overland solute flux from snowpack
  !     X*QSS=solute in snowpack flux from TranspNoSalt.f
  !     solute code:CO=CO2,CH=CH4,OX=O2,NG=N2,N2=N2O,HG=H2
  !             :OC=DOC,ON=DON,OP=DOP,OA=acetate
  !             :NH4=NH4,NH3=NH3,NO3=NO3,NO2=NO2,P14=HPO4,PO4=H2PO4 in non-band
  !             :N4B=NH4,N3B=NH3,NOB=NO3,N2B=NO2,P1B=HPO4,POB=H2PO4 in band
  !
  trcg_QSS(idg_CO2,N2,N1)=trcg_QSS(idg_CO2,N2,N1)+trcg_FloXSnow_2DH(idg_CO2,N,N2,N1)-trcg_FloXSnow_2DH(idg_CO2,N,N5,N4)
  trcg_QSS(idg_CH4,N2,N1)=trcg_QSS(idg_CH4,N2,N1)+trcg_FloXSnow_2DH(idg_CH4,N,N2,N1)-trcg_FloXSnow_2DH(idg_CH4,N,N5,N4)
  trcg_QSS(idg_O2,N2,N1)=trcg_QSS(idg_O2,N2,N1)+trcg_FloXSnow_2DH(idg_O2,N,N2,N1)-trcg_FloXSnow_2DH(idg_O2,N,N5,N4)
  trcg_QSS(idg_N2,N2,N1)=trcg_QSS(idg_N2,N2,N1)+trcg_FloXSnow_2DH(idg_N2,N,N2,N1)-trcg_FloXSnow_2DH(idg_N2,N,N5,N4)
  trcg_QSS(idg_N2O,N2,N1)=trcg_QSS(idg_N2O,N2,N1)+trcg_FloXSnow_2DH(idg_N2O,N,N2,N1)-trcg_FloXSnow_2DH(idg_N2O,N,N5,N4)
  trcg_QSS(idg_NH3,N2,N1)=trcg_QSS(idg_NH3,N2,N1)+trcg_FloXSnow_2DH(idg_NH3,N,N2,N1)-trcg_FloXSnow_2DH(idg_NH3,N,N5,N4)

  trcn_QSS(ids_NH4,N2,N1)=trcn_QSS(ids_NH4,N2,N1)+trcn_FloXSnow_2DH(ids_NH4,N,N2,N1)-trcn_FloXSnow_2DH(ids_NH4,N,N5,N4)
  trcn_QSS(ids_NO3,N2,N1)=trcn_QSS(ids_NO3,N2,N1)+trcn_FloXSnow_2DH(ids_NO3,N,N2,N1)-trcn_FloXSnow_2DH(ids_NO3,N,N5,N4)
  trcn_QSS(ids_H1PO4,N2,N1)=trcn_QSS(ids_H1PO4,N2,N1)+trcn_FloXSnow_2DH(ids_H1PO4,N,N2,N1)-trcn_FloXSnow_2DH(ids_H1PO4,N,N5,N4)
  trcn_QSS(ids_H2PO4,N2,N1)=trcn_QSS(ids_H2PO4,N2,N1)+trcn_FloXSnow_2DH(ids_H2PO4,N,N2,N1)-trcn_FloXSnow_2DH(ids_H2PO4,N,N5,N4)

  !     NET SALT FLUXES FROM RUNOFF AND SNOWPACK
  !
!
!     TQR*=net overland solute flux in runoff
!     XQR*=solute in runoff from TranspSalt.f
!     TQS*=net overland solute flux in snow drift
!     XQS*=solute in snow drift from TranspSalt.f
!     salt code: *HY*=H+,*OH*=OH-,*AL*=Al3+,*FE*=Fe3+,*CA*=Ca2+,*MG*=Mg2+
!          :*NA*=Na+,*KA*=K+,*SO4*=SO42-,*CL*=Cl-,*CO3*=CO32-,*HCO3*=HCO3-
!          :*CO2*=CO2,*ALO1*=AlOH2-,*ALOH2=AlOH2-,*ALOH3*=AlOH3
!          :*ALOH4*=AlOH4+,*ALS*=AlSO4+,*FEO1*=FeOH2-,*FEOH2=F3OH2-
!          :*FEOH3*=FeOH3,*FEOH4*=FeOH4+,*FES*=FeSO4+,*CAO*=CaOH
!          :*CAC*=CaCO3,*CAH*=CaHCO3-,*CAS*=CaSO4,*MGO*=MgOH,*MGC*=MgCO3
!          :*MHG*=MgHCO3-,*MGS*=MgSO4,*NAC*=NaCO3-,*NAS*=NaSO4-,*KAS*=KSO4-
!     phosphorus code: *H0P*=PO43-,*H3P*=H3PO4,*F1P*=FeHPO42-,*F2P*=F1H2PO4-
!          :*C0P*=CaPO4-,*C1P*=CaHPO4,*C2P*=CaH4P2O8+,*M1P*=MgHPO4,*COO*=COOH-
!          :*1=non-band,*B=band
!
  IF(salt_model)THEN
    D1203: DO NN=1,2

      DO NTSA=idsalt_beg,idsalt_end
        trcSalt_TQR(NTSA,N2,N1)=trcSalt_TQR(NTSA,N2,N1)+trc_salt_rof_bounds(NTSA,N,NN,N2,N1)
      ENDDO

      IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
! runoff direction
        DO NTSA=idsalt_beg,idsalt_end
          trcSalt_TQR(NTSA,N2,N1)=trcSalt_TQR(NTSA,N2,N1)-trc_salt_rof_bounds(NTSA,N,NN,N5,N4)
        ENDDO
      ENDIF

      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
        DO NTSA=idsalt_beg,idsalt_end
          trcSalt_TQR(NTSA,N2,N1)=trcSalt_TQR(NTSA,N2,N1)-trc_salt_rof_bounds(NTSA,N,NN,N5B,N4B)
        ENDDO
      ENDIF
    ENDDO D1203

    DO NTSA=idsalt_beg,idsalt_end
      trcSalt_TQS(NTSA,N2,N1)=trcSalt_TQS(NTSA,N2,N1)+trcSalt_XQS(NTSA,N,N2,N1)-trcSalt_XQS(NTSA,N,N5,N4)
    ENDDO
  ENDIF

  end subroutine MassFluxThruSnowRunoff
!------------------------------------------------------------------------------------------

  subroutine OverlandFlowThruSnow(I,J,NY,NX)
  implicit none 
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  integer :: NTSA,NTU,NTG
  real(r8) :: dflx

    !    SOLUTES
!  exclude NH3B
!  if(NX==5)write(111,*)I+J/24.,'beg18',trc_solml_vr(ids_NH4,0,1,5),trc_solml_vr(idg_NH3,0,1,5)

  DO NTG=idg_beg,idg_end-1
    dflx=-trcg_TFloXSurRunoff(NTG,NY,NX)
    call fixEXflux(trc_solml_vr(NTG,0,NY,NX),dflx)
    trcg_TFloXSurRunoff(NTG,NY,NX)=-dflx
  ENDDO

!  if(NX==5)write(111,*)I+J/24.,'beg19',trc_solml_vr(ids_NH4,0,1,5),trc_solml_vr(idg_NH3,0,1,5)

!  !   write(116,*)I+J/24.,'beg5xxxx1',trc_solml_vr(ids_NO3,0,NY,NX),trcn_TFloXSurRunoff_2D(ids_NO3,NY,NX)

  DO NTU=ids_nut_beg,ids_nuts_end
    dflx=-trcn_TFloXSurRunoff_2D(NTU,NY,NX)
    call fixEXflux(trc_solml_vr(NTU,0,NY,NX),dflx)
    trcn_TFloXSurRunoff_2D(NTU,NY,NX)=-dflx
  ENDDO

!  !   write(116,*)I+J/24.,'beg5xxxx2',trc_solml_vr(ids_NO3,0,NY,NX),trcn_TFloXSurRunoff_2D(ids_NO3,NY,NX)

!  if(NX==5)write(111,*)I+J/24.,'beg20',trc_solml_vr(ids_NH4,0,1,5),trc_solml_vr(idg_NH3,0,1,5)

  IF(salt_model)THEN
    DO NTSA=idsalt_beg,idsalt_end
      dflx=-trcSalt_TQR(NTSA,NY,NX)
      call fixEXflux(trcSalt_solml_vr(NTSA,0,NY,NX),dflx)
      trcSalt_TQR(NTSA,NY,NX)=-dflx
    ENDDO
  ENDIF
  end subroutine OverlandFlowThruSnow
end module SnowTransportMod