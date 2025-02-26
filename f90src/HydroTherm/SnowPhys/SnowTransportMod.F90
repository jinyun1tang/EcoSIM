module SnowTransportMod
!
!code to do water and tracer transport in snowpack
!
  use data_kind_mod,  only: r8 => DAT_KIND_R8
  use data_const_mod, only: spval => DAT_CONST_SPVAL
  use EcoSimConst,    only: DENSICE, natomw, patomw
  use EcoSIMCtrlMod,  only: lverb
  use MiniMathMod,    only: fixEXConsumpFlux
  use TracerPropMod,  only: MolecularWeight
  use ClimForcDataType
  use IrrigationDataType
  use SurfSoilDataType
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

  public :: SaltPercolThruSnow  
  public :: MassFluxThruSnowRunoff
  public :: DiagSnowChemMass
  public :: OverlandFlowThruSnow
  public :: TracerThruSnowfall
  public :: SoluteTransportThruSnow
  contains
!------------------------------------------------------------------------------------------

  subroutine SaltPercolThruSnow(I,J,N1,N2,NY,NX)
  !
  !update solute due to snow flux
  implicit none
  integer, intent(in) :: N1,N2,NY,NX,I,J
  integer :: LS, LS2
  integer :: idg,idn,idsalt,ids
!     begin_execution
!     NET WATER AND HEAT FLUXES THROUGH SNOWPACK
!  
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
        DO idg=idg_beg,idg_NH3
          trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)=trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)+trcg_AquaAdv_flx_snvr(idg,LS,N2,N1) &
            -trcg_AquaAdv_flx_snvr(idg,LS2,N2,N1)
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)=trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)+trcn_AquaAdv_flx_snvr(idn,LS,N2,N1) &
            -trcn_AquaAdv_flx_snvr(idn,LS2,N2,N1)
        ENDDO
        !
        !     NET SALT FLUXES THROUGH SNOWPACK
        !
        !
        IF(salt_model)THEN
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,N2,N1)=trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,N2,N1)+trcSalt_AquaAdv_flx_snvr(idsalt,LS,N2,N1)&
              -trcSalt_AquaAdv_flx_snvr(idsalt,LS2,N2,N1)
          ENDDO
        ENDIF
        !
        !     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
        !
      ELSE

        ! and NH3B
        DO idg=idg_beg,idg_NH3
          trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)=trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)+trcg_AquaAdv_flx_snvr(idg,LS,N2,N1) &
            -trcs_TransptMicP_3D(idg,3,0,N2,N1)-trcs_TransptMicP_3D(idg,3,NUM(N2,N1),N2,N1) &
            -trcs_TransptMacP_3D(idg,3,NUM(N2,N1),N2,N1)
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)=trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)+trcn_AquaAdv_flx_snvr(idn,LS,N2,N1) &
            -trcs_TransptMicP_3D(idn,3,0,N2,N1)-trcs_TransptMicP_3D(idn,3,NUM(N2,N1),N2,N1) &
            -trcs_TransptMacP_3D(idn,3,NUM(N2,N1),N2,N1)
        ENDDO

        !add band flux
        trcg_AquaAdv_NetFlx_snvr(idg_NH3,LS,N2,N1)=trcg_AquaAdv_NetFlx_snvr(idg_NH3,LS,N2,N1) &
          -trcs_TransptMicP_3D(idg_NH3B,3,NUM(N2,N1),N2,N1)-trcs_TransptMacP_3D(idg_NH3B,3,NUM(N2,N1),N2,N1)

        DO ids=0,ids_nuts
          trcn_AquaAdv_NetFlx_snvr(ids_NH4+ids,LS,N2,N1)=trcn_AquaAdv_NetFlx_snvr(ids_NH4+ids,LS,N2,N1) &
            -trcs_TransptMicP_3D(ids_NH4B+ids,3,NUM(N2,N1),N2,N1)-trcs_TransptMacP_3D(ids_NH4B+ids,3,NUM(N2,N1),N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,NY,NX)=trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,NY,NX)+trcSalt_AquaAdv_flx_snvr(idsalt,LS,NY,NX) &
              -trcSalt_TransptMicP_3D(idsalt,3,0,N2,N1)-trcSalt_TransptMicP_3D(idsalt,3,NUM(N2,N1),N2,N1) &
              -trcSalt_TransptMacP_3D(idsalt,3,NUM(N2,N1),N2,N1)
          ENDDO

          !add band flux
          DO idsalt=0,idsalt_nuts
            trcSalt_AquaAdv_NetFlx_snvr(idsalt_H0PO4+idsalt,LS,NY,NX)=trcSalt_AquaAdv_NetFlx_snvr(idsalt_H0PO4+idsalt,LS,NY,NX) &
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
          trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)=trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)+trcg_AquaAdv_flx_snvr(idg,LS,N2,N1)
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)=trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)+trcn_AquaAdv_flx_snvr(idn,LS,N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,N2,N1)=trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,N2,N1)+trcSalt_AquaAdv_flx_snvr(idsalt,LS,N2,N1)
          ENDDO
        ENDIF
      ENDIF
    ENDIF

  ENDDO D1205
  end subroutine SaltPercolThruSnow

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

  integer :: idsalt,idg,ids
!     begin_execution
!     OVERLAND SNOW REDISTRIBUTION
!
  IF(abs(TDrysnoBySnowRedist(NY,NX))>0._r8)THEN
    !loss of dissolved gases from surface snow
    DO idg=idg_beg,idg_NH3
      trcg_solsml_snvr(idg,1,NY,NX)=trcg_solsml_snvr(idg,1,NY,NX)+trcg_LossXSnowRedist_col(idg,NY,NX)
      GasHydroLossFlx_col(idg,NY,NX)=GasHydroLossFlx_col(idg,NY,NX)+trcg_LossXSnowRedist_col(idg,NY,NX)
    ENDDO

    DO ids=ids_nut_beg,ids_nuts_end
      trcn_solsml_snvr(ids,1,NY,NX)=trcn_solsml_snvr(ids,1,NY,NX)+trcn_LossXSnowRedist_col(ids,NY,NX)
    ENDDO

    IF(salt_model)THEN
      DO idsalt=idsalt_beg,idsalt_end
        trc_Saltml_snvr(idsalt,1,NY,NX)=trc_Saltml_snvr(idsalt,1,NY,NX)+trcSalt_LossXSnowRedist_col(idsalt,NY,NX)
      ENDDO
    ENDIF
  ENDIF

  end subroutine ChemicalBySnowRedistribution

!------------------------------------------------------------------------------------------
  subroutine MassFluxThruSnowRunoff(N,N1,N2,N4,N5,N4B,N5B)

  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B
  integer :: NN,idg,idn,idsalt

  TDrysnoBySnowRedist(N2,N1)   = TDrysnoBySnowRedist(N2,N1)+DrySnoBySnoRedistrib_2DH(N,N2,N1)-DrySnoBySnoRedistrib_2DH(N,N5,N4)
  TWatBySnowRedist(N2,N1)      = TWatBySnowRedist(N2,N1)+WatBySnowRedistrib_2DH(N,N2,N1)-WatBySnowRedistrib_2DH(N,N5,N4)
  TIceBySnowRedist(N2,N1)      = TIceBySnowRedist(N2,N1)+IceBySnowRedistrib_2DH(N,N2,N1)-IceBySnowRedistrib_2DH(N,N5,N4)
  THeatBySnowRedist_col(N2,N1) = THeatBySnowRedist_col(N2,N1)+HeatBySnowRedistrib_2DH(N,N2,N1)-HeatBySnowRedistrib_2DH(N,N5,N4)

  D1202: DO NN=1,2
    !gaseous tracers
    DO idg=idg_beg,idg_NH3
      trcg_SurfRunoff_flx(idg,N2,N1)=trcg_SurfRunoff_flx(idg,N2,N1)+trcg_FloXSurRunoff_2D(idg,N,NN,N2,N1)
    ENDDO

    !nutrient tracres
    DO idn=ids_nut_beg,ids_nuts_end
      trcn_SurfRunoff_flx(idn,N2,N1)=trcn_SurfRunoff_flx(idn,N2,N1)+trcn_FloXSurRunoff_2D(idn,N,NN,N2,N1)
    ENDDO

    IF(IFLBH(N,NN,N5,N4).EQ.0)THEN    

      DO idg=idg_beg,idg_NH3
        trcg_SurfRunoff_flx(idg,N2,N1)=trcg_SurfRunoff_flx(idg,N2,N1)-trcg_FloXSurRunoff_2D(idg,N,NN,N5,N4)
      ENDDO
      DO idn=ids_nut_beg,ids_nuts_end
        trcn_SurfRunoff_flx(idn,N2,N1)=trcn_SurfRunoff_flx(idn,N2,N1)-trcn_FloXSurRunoff_2D(idn,N,NN,N5,N4)
      ENDDO

    ENDIF 

    IF(N4B.GT.0 .AND. N5B.GT.0 .AND. NN.EQ.iOutflow)THEN
      DO idg=idg_beg,idg_NH3
        trcg_SurfRunoff_flx(idg,N2,N1)=trcg_SurfRunoff_flx(idg,N2,N1)-trcg_FloXSurRunoff_2D(idg,N,NN,N5B,N4B)
      ENDDO
      DO idn=ids_nut_beg,ids_nuts_end
        trcn_SurfRunoff_flx(idn,N2,N1)=trcn_SurfRunoff_flx(idn,N2,N1)-trcn_FloXSurRunoff_2D(idn,N,NN,N5B,N4B)
      ENDDO
    ENDIF
  ENDDO D1202
  !
  !     NET GAS AND SOLUTE FLUXES FROM RUNOFF AND SNOWPACK
  !
  do idg=idg_beg,idg_NH3
    trcg_LossXSnowRedist_col(idg,N2,N1) = trcg_LossXSnowRedist_col(idg,N2,N1)+&
      trcg_FloXSnow_2DH(idg,N,N2,N1)-trcg_FloXSnow_2DH(idg,N,N5,N4)
  ENDDO

  DO idn=ids_nut_beg,ids_nuts_end
    trcn_LossXSnowRedist_col(idn,N2,N1)   = trcn_LossXSnowRedist_col(idn,N2,N1) &
      +trcn_FloXSnow_2DH(idn,N,N2,N1)-trcn_FloXSnow_2DH(idn,N,N5,N4)
  ENDDO  

  !     NET SALT FLUXES FROM RUNOFF AND SNOWPACK
  !
!
  IF(salt_model)THEN
    D1203: DO NN=1,2
      DO idsalt=idsalt_beg,idsalt_end
        trcSalt_SurfRunoff_flx(idsalt,N2,N1)=trcSalt_SurfRunoff_flx(idsalt,N2,N1)+trcSalt_FloXSurRunoff_2D(idsalt,N,NN,N2,N1)
      ENDDO

      IF(IFLBH(N,NN,N5,N4).EQ.0)THEN
! runoff direction
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_SurfRunoff_flx(idsalt,N2,N1)=trcSalt_SurfRunoff_flx(idsalt,N2,N1)-trcSalt_FloXSurRunoff_2D(idsalt,N,NN,N5,N4)
        ENDDO
      ENDIF

      IF(N4B.GT.0 .AND. N5B.GT.0 .AND. NN.EQ.iOutflow)THEN
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_SurfRunoff_flx(idsalt,N2,N1)=trcSalt_SurfRunoff_flx(idsalt,N2,N1)-trcSalt_FloXSurRunoff_2D(idsalt,N,NN,N5B,N4B)
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
  integer :: idsalt,idn,idg
  real(r8) :: dflx

    !    SOLUTES
!  exclude NH3B

  DO idg=idg_beg,idg_NH3
    dflx=-trcg_SurfRunoff_flx(idg,NY,NX)
    call fixEXConsumpFlux(trcs_solml_vr(idg,0,NY,NX),dflx)
    trcg_SurfRunoff_flx(idg,NY,NX)=-dflx
    GasHydroLossFlx_col(idg,NY,NX)=GasHydroLossFlx_col(idg,NY,NX)+trcg_SurfRunoff_flx(idg,NY,NX)
  ENDDO

  DO idn=ids_nut_beg,ids_nuts_end
    dflx=-trcn_SurfRunoff_flx(idn,NY,NX)
    call fixEXConsumpFlux(trcs_solml_vr(idn,0,NY,NX),dflx)
    trcn_SurfRunoff_flx(idn,NY,NX)=-dflx
  ENDDO


  IF(salt_model)THEN
    DO idsalt=idsalt_beg,idsalt_end
      dflx=-trcSalt_SurfRunoff_flx(idsalt,NY,NX)
      call fixEXConsumpFlux(trcSalt_solml_vr(idsalt,0,NY,NX),dflx)
      trcSalt_SurfRunoff_flx(idsalt,NY,NX)=-dflx
    ENDDO
  ENDIF
  end subroutine OverlandFlowThruSnow
!------------------------------------------------------------------------------------------

  subroutine TracerThruSnowfall(I,J,NY,NX)
  implicit none

  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  integer :: idg

  DO idg=idg_beg,idg_NH3
    trcg_AquaAdv_flx_snvr(idg,1,NY,NX)=(Rain2SoilSurf_col(NY,NX)*trcg_rain_mole_conc_col(idg,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcg_irrig_mole_conc_col(idg,NY,NX))*MolecularWeight(idg)
  ENDDO

  trcn_AquaAdv_flx_snvr(ids_NH4,1,NY,NX)   = (Rain2SoilSurf_col(NY,NX)*NH4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*NH4_irrig_mole_conc(I,NY,NX))*natomw
  trcn_AquaAdv_flx_snvr(ids_NO3,1,NY,NX)   = (Rain2SoilSurf_col(NY,NX)*NO3_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*NO3_irrig_mole_conc(I,NY,NX))*natomw
  trcn_AquaAdv_flx_snvr(ids_H1PO4,1,NY,NX) = (Rain2SoilSurf_col(NY,NX)*HPO4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*HPO4_irrig_mole_conc(I,NY,NX))*patomw
  trcn_AquaAdv_flx_snvr(ids_H2PO4,1,NY,NX) = (Rain2SoilSurf_col(NY,NX)*H2PO4_rain_mole_conc(NY,NX) &
    +Irrig2SoilSurf_col(NY,NX)*H2PO4_irrig_mole_conc(I,NY,NX))*patomw
  end subroutine TracerThruSnowfall

!------------------------------------------------------------------------------------------

  subroutine SoluteTransportThruSnow(L,NY,NX)
  implicit none
  integer, intent(in) :: L,NY,NX

  integer :: idg,idn,idsalt
  !     begin_execution
  !     SNOWPACK SOLUTE CONTENT

  DO idg=idg_beg,idg_NH3
    trcg_solsml_snvr(idg,L,NY,NX)=trcg_solsml_snvr(idg,L,NY,NX)+trcg_AquaAdv_NetFlx_snvr(idg,L,NY,NX)
  ENDDO

  DO idn =ids_nut_beg,ids_nuts_end
    trcn_solsml_snvr(idn,L,NY,NX)=trcn_solsml_snvr(idn,L,NY,NX)+trcn_AquaAdv_NetFlx_snvr(idn,L,NY,NX)
  ENDDO
  !
  !
  IF(salt_model)THEN
    DO idsalt=idsalt_beg,idsalt_end
      trc_Saltml_snvr(idsalt,L,NY,NX)=trc_Saltml_snvr(idsalt,L,NY,NX)+trcSalt_AquaAdv_NetFlx_snvr(idsalt,L,NY,NX)
    ENDDO
  ENDIF

  end subroutine SoluteTransportThruSnow

end module SnowTransportMod