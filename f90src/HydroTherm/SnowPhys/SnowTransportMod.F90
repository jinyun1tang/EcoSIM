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
  integer :: idg,idn,idsalt,NTS
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
        DO idg=idg_beg,idg_NH3
          trcg_TBLS(idg,LS,N2,N1)=trcg_TBLS(idg,LS,N2,N1)+trcg_AquaAdv_flx_snvr(idg,LS,N2,N1) &
            -trcg_AquaAdv_flx_snvr(idg,LS2,N2,N1)
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_TBLS(idn,LS,N2,N1)=trcn_TBLS(idn,LS,N2,N1)+trcn_AquaAdv_flx_snvr(idn,LS,N2,N1) &
            -trcn_AquaAdv_flx_snvr(idn,LS2,N2,N1)
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
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_TBLS(idsalt,LS,N2,N1)=trcSalt_TBLS(idsalt,LS,N2,N1)+trcSalt_AquaAdv_flx_snvr(idsalt,LS,N2,N1)&
              -trcSalt_AquaAdv_flx_snvr(idsalt,LS2,N2,N1)
          ENDDO
        ENDIF
        !
        !     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
        !
      ELSE

        ! and NH3B
        DO idg=idg_beg,idg_NH3
          trcg_TBLS(idg,LS,N2,N1)=trcg_TBLS(idg,LS,N2,N1)+trcg_AquaAdv_flx_snvr(idg,LS,N2,N1) &
            -trcs_TransptMicP_3D(idg,3,0,N2,N1)-trcs_TransptMicP_3D(idg,3,NUM(N2,N1),N2,N1) &
            -trcs_TransptMacP_3D(idg,3,NUM(N2,N1),N2,N1)
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_TBLS(idn,LS,N2,N1)=trcn_TBLS(idn,LS,N2,N1)+trcn_AquaAdv_flx_snvr(idn,LS,N2,N1) &
            -trcs_TransptMicP_3D(idn,3,0,N2,N1)-trcs_TransptMicP_3D(idn,3,NUM(N2,N1),N2,N1) &
            -trcs_TransptMacP_3D(idn,3,NUM(N2,N1),N2,N1)
        ENDDO

        !add band flux
        trcg_TBLS(idg_NH3,LS,N2,N1)=trcg_TBLS(idg_NH3,LS,N2,N1) &
          -trcs_TransptMicP_3D(idg_NH3B,3,NUM(N2,N1),N2,N1)-trcs_TransptMacP_3D(idg_NH3B,3,NUM(N2,N1),N2,N1)

        DO NTS=0,ids_nuts
          trcn_TBLS(ids_NH4+NTS,LS,N2,N1)=trcn_TBLS(ids_NH4+NTS,LS,N2,N1) &
            -trcs_TransptMicP_3D(ids_NH4B+NTS,3,NUM(N2,N1),N2,N1)-trcs_TransptMacP_3D(ids_NH4B+NTS,3,NUM(N2,N1),N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_TBLS(idsalt,LS,NY,NX)=trcSalt_TBLS(idsalt,LS,NY,NX)+trcSalt_AquaAdv_flx_snvr(idsalt,LS,NY,NX) &
              -trcSalt_TransptMicP_3D(idsalt,3,0,N2,N1)-trcSalt_TransptMicP_3D(idsalt,3,NUM(N2,N1),N2,N1) &
              -trcSalt_TransptMacP_3D(idsalt,3,NUM(N2,N1),N2,N1)
          ENDDO

          !add band flux
          DO idsalt=0,idsalt_nuts
            trcSalt_TBLS(idsalt_H0PO4+idsalt,LS,NY,NX)=trcSalt_TBLS(idsalt_H0PO4+idsalt,LS,NY,NX) &
              -trcSalt_TransptMicP_3D(idsalt_H0PO4B+idsalt,3,NUM(N2,N1),N2,N1) &
              -trcSalt_TransptMacP_3D(idsalt_H0PO4B+idsalt,3,NUM(N2,N1),N2,N1)
          ENDDO
        ENDIF
      ENDIF
!
!     WATER,GAS,SOLUTE,SALT FLUXES INTO SNOWPACK SURFACE
!
    ELSEIF(LS.EQ.1)THEN
      IF(abs(SnoXfer2SnoLay_snvr(LS,N2,N1))>0._r8)THEN

        DO idg=idg_beg,idg_NH3
          trcg_TBLS(idg,LS,N2,N1)=trcg_TBLS(idg,LS,N2,N1)+trcg_AquaAdv_flx_snvr(idg,LS,N2,N1)
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_TBLS(idn,LS,N2,N1)=trcn_TBLS(idn,LS,N2,N1)+trcn_AquaAdv_flx_snvr(idn,LS,N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_TBLS(idsalt,LS,N2,N1)=trcSalt_TBLS(idsalt,LS,N2,N1)+trcSalt_AquaAdv_flx_snvr(idsalt,LS,N2,N1)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

  ENDDO D1205
  end subroutine VerticalSaltFluxThruSnowpack

!------------------------------------------------------------------------------------------

  subroutine DiagSnowChemMass(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8) :: SSW,WS
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
      TGasC_lnd           = TGasC_lnd+trcg_solsml_snvr(idg_CO2,L,NY,NX)+trcg_solsml_snvr(idg_CH4,L,NY,NX)
      DIC_mass_col(NY,NX) = DIC_mass_col(NY,NX)+trcg_solsml_snvr(idg_CO2,L,NY,NX)+trcg_solsml_snvr(idg_CH4,L,NY,NX)
      TSoilO2G_lnd        = TSoilO2G_lnd+trcg_solsml_snvr(idg_O2,L,NY,NX)
      TGasN_lnd           = TGasN_lnd+trcg_solsml_snvr(idg_N2,L,NY,NX)+trcg_solsml_snvr(idg_N2O,L,NY,NX)
      TDisolNH4_lnd       = TDisolNH4_lnd+trcn_solsml_snvr(ids_NH4,L,NY,NX)+trcg_solsml_snvr(idg_NH3,L,NY,NX)
      tNO3_lnd            = tNO3_lnd+trcn_solsml_snvr(ids_NO3,L,NY,NX)
      TDisolPi_lnd        = TDisolPi_lnd+trcn_solsml_snvr(ids_H1PO4,L,NY,NX)+trcn_solsml_snvr(ids_H2PO4,L,NY,NX)

      IF(salt_model)THEN
        SSW=0._r8
        do nsalts=idsalt_beg,idsalt_end
          SSW=SSW+trc_Saltml_snvr(nsalts,L,NY,NX)*trcSaltIonNumber(nsalts)
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

  integer :: NTA,idg,NTS
!     begin_execution
!     OVERLAND SNOW REDISTRIBUTION
!
  IF(abs(TDrysnoBySnowRedist(NY,NX))>0._r8)THEN
    !loss of dissolved gases from surface snow
    DO idg=idg_beg,idg_NH3
      trcg_solsml_snvr(idg,1,NY,NX)=trcg_solsml_snvr(idg,1,NY,NX)+trcg_LossXSnowRedist_col(idg,NY,NX)
    ENDDO

    DO NTS=ids_nut_beg,ids_nuts_end
      trcn_solsml_snvr(NTS,1,NY,NX)=trcn_solsml_snvr(NTS,1,NY,NX)+trcn_LossXSnowRedist_col(NTS,NY,NX)
    ENDDO

    IF(salt_model)THEN
      DO NTA=idsalt_beg,idsalt_end
        trc_Saltml_snvr(NTA,1,NY,NX)=trc_Saltml_snvr(NTA,1,NY,NX)+trcSalt_LossXSnowRedist_col(NTA,NY,NX)
      ENDDO
    ENDIF
  ENDIF

  end subroutine ChemicalBySnowRedistribution

!------------------------------------------------------------------------------------------
  subroutine MassFluxThruSnowRunoff(N,N1,N2,N4,N5,N4B,N5B)

  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B
  integer :: NN,idg,idn,idsalt

  TDrysnoBySnowRedist(N2,N1)   = TDrysnoBySnowRedist(N2,N1)+DrySnoBySnoRedistrib_2DH(N,N2,N1)&
    -DrySnoBySnoRedistrib_2DH(N,N5,N4)
  TWatBySnowRedist(N2,N1)      = TWatBySnowRedist(N2,N1)+WatBySnowRedistrib_2DH(N,N2,N1) &
    -WatBySnowRedistrib_2DH(N,N5,N4)
  TIceBySnowRedist(N2,N1)      = TIceBySnowRedist(N2,N1)+IceBySnowRedistrib_2DH(N,N2,N1) &
    -IceBySnowRedistrib_2DH(N,N5,N4)
  THeatBySnowRedist_col(N2,N1) = THeatBySnowRedist_col(N2,N1)+HeatBySnowRedistrib_2DH(N,N2,N1) &
    -HeatBySnowRedistrib_2DH(N,N5,N4)

  D1202: DO NN=1,2
    !gaseous tracers
    DO idg=idg_beg,idg_NH3
      trcg_SurfRunoff_flxM(idg,N2,N1)=trcg_SurfRunoff_flxM(idg,N2,N1)+trcg_FloXSurRunoff_2D(idg,N,NN,N2,N1)
    ENDDO

    !nutrient tracres
    DO idn=ids_nut_beg,ids_nuts_end
      trcn_SurfRunoff_flxM(idn,N2,N1)=trcn_SurfRunoff_flxM(idn,N2,N1)+trcn_FloXSurRunoff_2D(idn,N,NN,N2,N1)
    ENDDO

    IF(IFLBH(N,NN,N5,N4).EQ.0)THEN    

      DO idg=idg_beg,idg_NH3
        trcg_SurfRunoff_flxM(idg,N2,N1)=trcg_SurfRunoff_flxM(idg,N2,N1)-trcg_FloXSurRunoff_2D(idg,N,NN,N5,N4)
      ENDDO
      DO idn=ids_nut_beg,ids_nuts_end
        trcn_SurfRunoff_flxM(idn,N2,N1)=trcn_SurfRunoff_flxM(idn,N2,N1)-trcn_FloXSurRunoff_2D(idn,N,NN,N5,N4)
      ENDDO

    ENDIF 

    IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
      DO idg=idg_beg,idg_NH3
        trcg_SurfRunoff_flxM(idg,N2,N1)=trcg_SurfRunoff_flxM(idg,N2,N1)-trcg_FloXSurRunoff_2D(idg,N,NN,N5B,N4B)
      ENDDO
      DO idn=ids_nut_beg,ids_nuts_end
        trcn_SurfRunoff_flxM(idn,N2,N1)=trcn_SurfRunoff_flxM(idn,N2,N1)-trcn_FloXSurRunoff_2D(idn,N,NN,N5B,N4B)
      ENDDO

    ENDIF
  ENDDO D1202

  !
  !     NET GAS AND SOLUTE FLUXES FROM RUNOFF AND SNOWPACK
  !
  do idg=idg_beg,idg_NH3
    trcg_LossXSnowRedist_col(idg,N2,N1) = trcg_LossXSnowRedist_col(idg,N2,N1)+trcg_FloXSnow_2DH(idg,N,N2,N1)-trcg_FloXSnow_2DH(idg,N,N5,N4)
  ENDDO

  DO idn=ids_nut_beg,ids_nuts_end
    trcn_LossXSnowRedist_col(idn,N2,N1)   = trcn_LossXSnowRedist_col(idn,N2,N1) &
      +trcn_FloXSnow_2DH(idn,N,N2,N1)-trcn_FloXSnow_2DH(idn,N,N5,N4)
  ENDDO  

  !     NET SALT FLUXES FROM RUNOFF AND SNOWPACK
  !
!
!     TQR*=net overland solute flux in runoff
!     XQR*=solute in runoff from TranspSalt.f
!     TQS*=net overland solute flux in snow drift
!     XQS*=solute in snow drift from TranspSalt.f
!
  IF(salt_model)THEN
    D1203: DO NN=1,2

      DO idsalt=idsalt_beg,idsalt_end
        trcSalt_TQR(idsalt,N2,N1)=trcSalt_TQR(idsalt,N2,N1)+trc_salt_rof_bounds(idsalt,N,NN,N2,N1)
      ENDDO

      IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
! runoff direction
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_TQR(idsalt,N2,N1)=trcSalt_TQR(idsalt,N2,N1)-trc_salt_rof_bounds(idsalt,N,NN,N5,N4)
        ENDDO
      ENDIF

      IF(N4B.GT.0.AND.N5B.GT.0.AND.NN.EQ.1)THEN
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_TQR(idsalt,N2,N1)=trcSalt_TQR(idsalt,N2,N1)-trc_salt_rof_bounds(idsalt,N,NN,N5B,N4B)
        ENDDO
      ENDIF
    ENDDO D1203

    DO idsalt=idsalt_beg,idsalt_end
      trcSalt_LossXSnowRedist_col(idsalt,N2,N1)=trcSalt_LossXSnowRedist_col(idsalt,N2,N1) &
        +trcSalt_FloXSnow_2DH(idsalt,N,N2,N1)-trcSalt_FloXSnow_2DH(idsalt,N,N5,N4)
    ENDDO
  ENDIF

  end subroutine MassFluxThruSnowRunoff
!------------------------------------------------------------------------------------------

  subroutine OverlandFlowThruSnow(I,J,NY,NX)
  implicit none 
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  integer :: idsalt,NTU,idg
  real(r8) :: dflx

    !    SOLUTES
!  exclude NH3B
!  if(NX==5)write(111,*)I+J/24.,'beg18',trcs_solml_vr(ids_NH4,0,1,5),trcs_solml_vr(idg_NH3,0,1,5)

  DO idg=idg_beg,idg_NH3
    dflx=-trcg_SurfRunoff_flxM(idg,NY,NX)
    call fixEXflux(trcs_solml_vr(idg,0,NY,NX),dflx)
    trcg_SurfRunoff_flxM(idg,NY,NX)=-dflx
  ENDDO

!  if(NX==5)write(111,*)I+J/24.,'beg19',trcs_solml_vr(ids_NH4,0,1,5),trcs_solml_vr(idg_NH3,0,1,5)

!  !   write(116,*)I+J/24.,'beg5xxxx1',trcs_solml_vr(ids_NO3,0,NY,NX),trcn_SurfRunoff_flxM(ids_NO3,NY,NX)

  DO NTU=ids_nut_beg,ids_nuts_end
    dflx=-trcn_SurfRunoff_flxM(NTU,NY,NX)
    call fixEXflux(trcs_solml_vr(NTU,0,NY,NX),dflx)
    trcn_SurfRunoff_flxM(NTU,NY,NX)=-dflx
  ENDDO

!  !   write(116,*)I+J/24.,'beg5xxxx2',trcs_solml_vr(ids_NO3,0,NY,NX),trcn_SurfRunoff_flxM(ids_NO3,NY,NX)

!  if(NX==5)write(111,*)I+J/24.,'beg20',trcs_solml_vr(ids_NH4,0,1,5),trcs_solml_vr(idg_NH3,0,1,5)

  IF(salt_model)THEN
    DO idsalt=idsalt_beg,idsalt_end
      dflx=-trcSalt_TQR(idsalt,NY,NX)
      call fixEXflux(trcSalt_solml_vr(idsalt,0,NY,NX),dflx)
      trcSalt_TQR(idsalt,NY,NX)=-dflx
    ENDDO
  ENDIF
  end subroutine OverlandFlowThruSnow
end module SnowTransportMod