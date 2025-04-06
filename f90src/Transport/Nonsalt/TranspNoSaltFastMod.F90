module TranspNoSaltFastMod

  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod
  use DebugToolMod
  use GridDataType  
  use TranspNoSaltDataMod
  use SoilBGCDataType
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SurfLitterDataType
  use AqueChemDatatype
  use SnowDataType
  use EcoSimConst
  use TracerIDMod
  use ChemTranspDataType
  use ClimForcDataType
  use EcoSIMSolverPar
  use LandSurfDataType
  use SurfSoilDataType
  use SoilPropertyDataType
implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: TransptFastNoSaltMM
  contains


  subroutine TransptFastNoSaltMM(I,J,M,NHE,NHW,NVS,NVN)

  !Description:
  !update trcs_solml2_vr and trc_gasml2_vr
  implicit none
  integer, intent(in) :: I,J,M,NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='TransptFastNoSaltMM'
  real(r8) :: pscal(idg_beg:idg_end), dpscal(idg_beg:idg_end)
  real(r8) :: dpscal_max
  integer :: ids

  call PrintInfo('beg '//subname)
  dpscal=1._r8
  dpscal_max=1._r8
  DO while(dpscal_max>1.e-2_r8)
    pscal=1.001_r8

    call ZeroTracerFluxMM(I,J,NHE,NHW,NVS,NVN)

    call LitterGasVolatilDissolMM(I,J,M,NHE,NHW,NVS,NVN)

    call SurfSoilFluxGasDifAdvMM(M,NHE,NHW,NVS,NVN)

    call XBoundaryFluxMM(I,J,M,NHW,NHE,NVN,NVS)

    call TracerFlowXGridsMM(I,J,M,NHE,NHW,NVS,NVN)

    call FastUpdateStateVarsMM(I,J,NHW, NHE, NVN, NVS,dpscal,pscal)

    call AccumSlowFluxesMM(I,J,NHE,NHW,NVS,NVN,dpscal,pscal)

    call copyStateVars(NHE,NHW,NVS,NVN)

    dpscal_max=0._r8
    DO ids=ids_beg,ids_end
      dpscal(ids)=dpscal(ids)*(1._r8-pscal(ids))
      dpscal_max=AMAX1(dpscal_max,dpscal(ids))
    enddo
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine TransptFastNoSaltMM
!------------------------------------------------------------------------------------------
  subroutine copyStateVars(NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: NHE,NHW,NVS,NVN
  integer :: NY,NX,idg,L

  DO NY=NHW,NHE
    DO  NX=NVN,NVS
      DO idg=idg_beg,idg_NH3      
       trcs_solml2_vr(idg,0,NY,NX) = trcs_solml_vr(idg,0,NY,NX)
       DO L=NU(NY,NX),NL(NY,NX)
         trc_gasml2_vr(idg,L,NY,NX)=trcg_gasml_vr(idg,L,NY,NX)
       ENDDO
      enddo

      DO idg=idg_beg,idg_end 
        DO L=NU(NY,NX),NL(NY,NX)    
          trcs_solml2_vr(idg,L,NY,NX) = trcs_solml_vr(idg,L,NY,NX)
        ENDDO
      enddo

    ENDDO
  ENDDO
  end subroutine copyStateVars
!------------------------------------------------------------------------------------------
  subroutine AccumSlowFluxesMM(I,J,NHE,NHW,NVS,NVN,dpscal,pscal)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHE,NHW,NVS,NVN
  real(r8), intent(in) :: dpscal(idg_beg:idg_end)
  real(r8), intent(in):: pscal(idg_beg:idg_end)
  character(len=*), parameter :: subname='AccumSlowFluxesMM'
  real(r8) :: ppscal(idg_beg:idg_end)
  real(r8) :: pscal_max
  integer :: idg,L,NY,NX


  call PrintInfo('beg '//subname)

  pscal_max=0._r8   
  DO idg=idg_beg,idg_end
    ppscal(idg)=dpscal(idg)*pscal(idg)
    pscal_max=AMAX1(pscal(idg),pscal_max)
  ENDDO

  if(pscal_max<0.999_r8)then
    call FastUpdateStateVarsMM(I,J,NHW, NHE, NVN, NVS,ppscal)
  endif

  DO NY=NHW,NHE
    DO  NX=NVN,NVS
      DO idg=idg_beg,idg_NH3          
        GasDiff2Surf_flx_col(idg,NY,NX) = GasDiff2Surf_flx_col(idg,NY,NX)+(RGas_Disol_FlxMM_vr(idg,0,NY,NX) &
          +Gas_AdvDif_FlxMM_3D(idg,3,NU(NY,NX),NY,NX))*ppscal(idg)

        L=0
        RGasTranspFlxPrev_vr(idg,L,NY,NX)  =RGasTranspFlxPrev_vr(idg,L,NY,NX)+ Gas_AdvDif_FlxMM_vr(idg,L,NY,NX)*ppscal(idg)           
        DO L=NU(NY,NX),NL(NY,NX)
          RGasTranspFlxPrev_vr(idg,L,NY,NX)  =RGasTranspFlxPrev_vr(idg,L,NY,NX)+ Gas_AdvDif_FlxMM_vr(idg,L,NY,NX)*ppscal(idg)   
        ENDDO
      ENDDO

    ENDDO
  ENDDO  
  call PrintInfo('end '//subname)
  end subroutine AccumSlowFluxesMM

!------------------------------------------------------------------------------------------
  subroutine ZeroTracerFluxMM(I,J,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J,NHE,NHW,NVS,NVN

  Gas_AdvDif_FlxMM_vr(idg_beg:idg_NH3,:,:,:)   = 0._r8
  Gas_AdvDif_FlxMM_3D(idg_beg:idg_NH3,:,:,:,:) = 0.0_r8

  end subroutine ZeroTracerFluxMM

!------------------------------------------------------------------------------------------

  subroutine LitterGasVolatilDissolMM(I,J,M,NHE,NHW,NVS,NVN)
  !
  !Description:
  !Gas dissolution flux into litter layer
  !This differs from the slow process in subroutine LitterAtmosExchangeM.
  !where the latter used aqueous diffusivity. 
  !These two can be viewed as non-eqiulbrium version of the problem tackled in
  !Tang and Riley (2013), Hydrol. Earth Syst. Sci, 17, 873–893, https://doi.org/10.5194/hess-17-873-2013.
  !
  implicit none

  integer, intent(in) :: I,J, M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  integer :: NY, NX
  real(r8) :: VOLGas
  real(r8) :: trc_gascl
  integer  :: idg
!
  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      IF(VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX) .AND. VLsoiAirPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX) &
        .AND. VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN
        VLWatMicPXA_vr(0,NY,NX)=natomw*VLWatMicPM_vr(M,0,NY,NX)

        do idg=idg_beg,idg_NH3
          trcg_VLWatMicP_vr(idg,0,NY,NX)  = VLWatMicPM_vr(M,0,NY,NX)*GasSolbility_vr(idg,0,NY,NX)
          trc_gascl                       = trcg_gascl_vr(idg,0,NY,NX)*VLsoiAirPM_vr(M,0,NY,NX)
          VOLGas                          = trcg_VLWatMicP_vr(idg,0,NY,NX)+VLsoiAirPM_vr(M,0,NY,NX)
          RGas_Disol_FlxMM_vr(idg,0,NY,NX) = DiffusivitySolutEffM_vr(M,0,NY,NX) &
            *(AMAX1(ZEROS(NY,NX),trc_gascl)*trcg_VLWatMicP_vr(idg,0,NY,NX) &
            - AMAX1(ZEROS(NY,NX),trcs_solml2_vr(idg,0,NY,NX))*VLsoiAirPM_vr(M,0,NY,NX))/VOLGas

        ENDDO
      ENDIF
    ENDDO
  ENDDO
  end subroutine LitterGasVolatilDissolMM

!------------------------------------------------------------------------------------------

  subroutine SurfSoilFluxGasDifAdvMM(M,NHE,NHW,NVS,NVN)
!
! DESCRIPTION:
! surface soil gaseous diffusion, advection, dissolution & volatilization
  implicit none

  integer, intent(in) :: M, NHE,NHW,NVS,NVN

  character(len=*), parameter :: subname='SurfSoilFluxGasDifAdvMM'
  real(r8) :: VFLW
  integer :: idg
  integer :: NY, NX

  call PrintInfo('beg '//subname)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      !Only deal with air-filled grids
      IF(FracAirFilledSoilPoreM_vr(M,NU(NY,NX),NY,NX).GT.AirFillPore_Min &
        .AND. SoilBulkDensity_vr(NU(NY,NX),NY,NX).GT.ZERO)THEN

        call TopSoilGasDifussionMM(M,NY,NX)

        call TopSoilGasAdvectionMM(M,NY,NX)

      ENDIF
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SurfSoilFluxGasDifAdvMM

! ----------------------------------------------------------------------
  subroutine TopSoilGasAdvectionMM(M,NY,NX)
  !
  !Description
  !It assumes that air is not compressible, and the total soil pore volume 
  !is changing slowly (if it does change).
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VFLW
  real(r8) :: RGas_Adv_flxMM(idg_beg:idg_NH3)
  integer :: idg


  IF(WaterFlow2SoilMM_3D(3,NU(NY,NX),NY,NX).GT.0.0_r8)THEN
    !topsoil loses tracer
    IF(VLsoiAirPM_vr(M,NU(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=-AZMAX1(AMIN1(VFLWX,WaterFlow2SoilMM_3D(3,NU(NY,NX),NY,NX)/VLsoiAirPM_vr(M,NU(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    
    DO idg=idg_beg,idg_NH3
      RGas_Adv_flxMM(idg)=VFLW*AZMAX1(trc_gasml2_vr(idg,NU(NY,NX),NY,NX))
    ENDDO
  ELSE
    !topsoil gains tracers
    DO idg=idg_beg,idg_NH3
      RGas_Adv_flxMM(idg)=-WaterFlow2SoilMM_3D(3,NU(NY,NX),NY,NX)*AtmGasCgperm3(idg,NY,NX)
    ENDDO
  ENDIF
  !
  !     TOTAL SOIL GAS FLUX + CONVECTIVE FLUX
  !
  DO idg=idg_beg,idg_NH3
    Gas_AdvDif_FlxMM_3D(idg,3,NU(NY,NX),NY,NX)=Gas_AdvDif_FlxMM_3D(idg,3,NU(NY,NX),NY,NX)+RGas_Adv_flxMM(idg)
  ENDDO
  end subroutine TopSoilGasAdvectionMM

! ----------------------------------------------------------------------
  subroutine GasTransportMM(M,N,N1,N2,N3,N4,N5,N6)

  !Description
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  character(len=*), parameter :: subname='GasTransportMM'
  integer :: idg

  call PrintInfo('beg '//subname)

  IF(FracAirFilledSoilPoreM_vr(M,N3,N2,N1).GT.AirFillPore_Min      &    !grid supports gas
    .AND. FracAirFilledSoilPoreM_vr(M,N6,N5,N4).GT.AirFillPore_Min &    !grid supports gas  
    .AND. VLsoiAirPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1)               &    !source grid has significant air volume
    .AND. VLsoiAirPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN               !dest grid has significant air volume

!     TOTAL SOIL GAS FLUX FROM DIFFUSIVE
    call GasDiffusionMM(M,N,N1,N2,N3,N4,N5,N6)

!     TOTAL SOIL GAS FLUX FROM CONVECTIVE FLUX
    call UpstreamGasAdvectionMM(M,N,N1,N2,N3,N4,N5,N6)

  ENDIF
  call PrintInfo('end '//subname)
  end subroutine GasTransportMM

! ----------------------------------------------------------------------
  subroutine GasDissolutionMM(M,N,N4,N5,N6)
  !
  !Description
  !Compute gas dissolution 
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: N4,N5,N6

  integer :: idg
!
!
  IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
    IF(FracAirFilledSoilPoreM_vr(M,N6,N5,N4).GT.AirFillPore_Min)THEN
      do idg=idg_beg,idg_NH3-1
        RGas_Disol_FlxMM_vr(idg,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
         (AMAX1(ZEROS(N5,N4),trc_gasml2_vr(idg,N6,N5,N4))*trcg_VLWatMicP_vr(idg,N6,N5,N4) &
          -trcs_solml2_vr(idg,N6,N5,N4)*VLsoiAirPM_vr(M,N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg,N6,N5,N4)+VLsoiAirPM_vr(M,N6,N5,N4))
      enddo    

      IF(VLsoiAirPMA_vr(N6,N5,N4).GT.ZEROS2(N5,N4).AND.VLWatMicPXA_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
        RGas_Disol_FlxMM_vr(idg_NH3,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
         (AMAX1(ZEROS(N5,N4),trc_gasml2_vr(idg_NH3,N6,N5,N4))*trcg_VLWatMicP_vr(idg_NH3,N6,N5,N4) &
          -trcs_solml2_vr(idg_NH3,N6,N5,N4)*VLsoiAirPMA_vr(N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg_NH3,N6,N5,N4)+VLsoiAirPMA_vr(N6,N5,N4))
      ENDIF

      IF(VLsoiAirPMB_vr(N6,N5,N4).GT.ZEROS2(N5,N4).AND.VLWatMicPXB_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
        RGas_Disol_FlxMM_vr(idg_NH3B,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
          (AMAX1(ZEROS(N5,N4),trc_gasml2_vr(idg_NH3,N6,N5,N4))*trcg_VLWatMicP_vr(idg_NH3B,N6,N5,N4) &
          -trcs_solml2_vr(idg_NH3B,N6,N5,N4)*VLsoiAirPMB_vr(N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg_NH3B,N6,N5,N4)+VLsoiAirPMB_vr(N6,N5,N4))
      ENDIF

    ENDIF
  ENDIF
  end subroutine GasDissolutionMM

! ----------------------------------------------------------------------
  subroutine UpstreamGasAdvectionMM(M,N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  character(len=*), parameter :: subname='UpstreamGasAdvectionMM'
  real(r8) :: RGasAdv
  real(r8) :: VFLW,FLQW
  integer :: idg
!
!     CONVECTIVE GAS TRANSFER DRIVEN BY SOIL WATER FLUXES
!     FROM 'WATSUB' AND GAS CONCENTRATIONS IN THE ADJACENT GRID CELLS
!     DEPENDING ON WATER FLUX DIRECTION
!
!     by assuming volume conservation, gases and water flow in opposite direction
!
  call PrintInfo('beg '//subname)
  FLQW=WaterFlow2SoilMM_3D(N,N6,N5,N4)  
  
  !flow out of grid (N6,N5,N4)  !water flow into grid, gas out of grid
  IF(FLQW.GT.0.0_r8)THEN
    !dest grid is not saturated
    IF(VLsoiAirPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN   
      VFLW=-AZMAX1(AMIN1(VFLWX,FLQW/VLsoiAirPM_vr(M,N6,N5,N4)))   !negative flow 
    !dest grid is aturated  
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO idg=idg_beg,idg_NH3
      RGasAdv                             = VFLW*AZMAX1(trc_gasml2_vr(idg,N6,N5,N4))
      Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4) = Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4)+RGasAdv
    ENDDO
  !flow out of source grid  
  ELSE
    IF(VLsoiAirPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=-AZMIN1(AMAX1(-VFLWX,FLQW/VLsoiAirPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO idg=idg_beg,idg_NH3
      RGasAdv                             = VFLW*AZMAX1(trc_gasml2_vr(idg,N3,N2,N1))
      Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4) = Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4)+RGasAdv
    ENDDO
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine UpstreamGasAdvectionMM

!------------------------------------------------------------------------------------------

  subroutine XBoundaryFluxMM(I,J,M,NHW,NHE,NVN,NVS)
  !
  !Do transport across the boundaries
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NHW, NHE, NVN, NVS

  character(len=*), parameter :: subname='XBoundaryFluxMM'
  integer :: NY,NX,L
  integer :: NN,N
  integer :: N1,N2,N3,N4,N5,N6,N4B,N5B  !inner grids
  integer :: M1,M2,M3,M4,M5,M6          !boundary exchange
  integer :: LL
  call PrintInfo('beg '//subname)
  !
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      D9585: DO L=NU(NY,NX),NL(NY,NX)
        N1=NX;N2=NY;N3=L
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
        D9580: DO  N=FlowDirIndicator_col(NY,NX),3
          D9575: DO  NN=1,2
            IF(N.EQ.iWestEastDirection)THEN
              !WEST-EAST
              N4 = NX+1; N5 = NY ;N6 = L
              IF(NN.EQ.iFront)THEN   !eastward                
                IF(NX.EQ.NHE)THEN !eastern boundary
                  M1 = NX;M2   = NY;M3 = L
                  M4 = NX+1;M5 = NY;M6 = L
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN  !west                
                IF(NX.EQ.NHW)THEN  !western boundary
                  M1 = NX;M2 = NY;M3 = L
                  M4 = NX;M5 = NY;M6 = L
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN              
              N4  = NX;N5  = NY+1;N6  = L
              IF(NN.EQ.iFront)THEN  
                IF(NY.EQ.NVS)THEN   ! southern boundary
                  M1 = NX;M2 = NY;M3   = L   !source grid
                  M4 = NX;M5 = NY+1;M6 = L   !target grid
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN  !north                
                IF(NY.EQ.NVN)THEN ! northern boundary
                  M1 = NX;M2 = NY;M3 = L !source                   
                  M4 = NX;M5 = NY;M6 = L !target
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iVerticalDirection)THEN !vertical

              N4 = NX;N5 = NY;N6 = L+1  !target
              IF(NN.EQ.iFront)THEN                
                IF(L.EQ.NL(NY,NX))THEN      !lower boundary
                  M1 = NX;M2 = NY;M3 = L    !source grid
                  M4 = NX;M5 = NY;M6 = L+1  !target grid
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN
                !nothing for the upper boundary
                cycle
              ENDIF
            ENDIF            
            !     SURFACE SOLUTE TRANSPORT FROM BOUNDARY SURFACE
            !     RUNOFF IN 'WATSUB' AND CONCENTRATIONS IN THE SURFACE SOIL LAYER
            !           
            call XBoundaryTracerFlowMM(N,NN,M,M1,M2,M3,M4,M5,M6)
          ENDDO D9575
          !
          !     TOTAL SOLUTE FLUX IN MICROPORES AND MACROPORES
          !
          IF(FlowDirIndicator_col(N2,N1).NE.3 .OR. N.EQ.iVerticalDirection)THEN
            DO LL=N6,NL(NY,NX)
              IF(VLSoilPoreMicP_vr(LL,N2,N1).GT.ZEROS2(N2,N1))THEN
                N6=LL
                exit
              ENDIF
            ENDDO
            call NetTracerFlowXSoilPoresMM(NY,NX,N,M,N1,N2,N3,N4,N5,N6)
          ENDIF
        ENDDO D9580
      ENDDO D9585

    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine XBoundaryFluxMM

! ----------------------------------------------------------------------
  subroutine XBoundaryTracerFlowMM(N,NN,M,M1,M2,M3,M4,M5,M6)
  implicit none

  integer, intent(in) :: N, NN, M
  integer, intent(in) :: M1, M2,M3,M4, M5,M6

  character(len=*), parameter :: subname='XBoundaryTracerFlowMM'
  real(r8) :: FLGM,FQRM,VFLW
  integer :: K,idg

! begin_execution
  call PrintInfo('beg '//subname)

!   dt_GasCyc=1/NPT, NPT is number of gas iterations per M

  FLGM=WaterFlow2SoilMM_3D(N,M6,M5,M4)
  
  !make sure out of grid (M3,M2,M1)
  IF(NN.EQ.iFront .AND. FLGM.LT.0.0_r8           &  
    .OR. (NN.EQ.iBehind .AND. FLGM.GT.0.0_r8))THEN  
    IF(VLsoiAirPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=-AMAX1(-VFLWX,AMIN1(VFLWX,FLGM/VLsoiAirPM_vr(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF
    !
    ! HOURLY GAS FLUX FOR USE IN REDIST.F
    !
    DO idg=idg_beg,idg_NH3
      Gas_AdvDif_FlxMM_3D(idg,N,M6,M5,M4) = Gas_AdvDif_FlxMM_3D(idg,N,M6,M5,M4)+VFLW*AZMAX1(trc_gasml2_vr(idg,M3,M2,M1))
    ENDDO
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine XBoundaryTracerFlowMM
!------------------------------------------------------------------------------------------

  subroutine NetTracerFlowXSoilPoresMM(NY,NX,N,M,N1,N2,N3,N4,N5,N6)
  !
  !Description
  !Vertical or lateral flux between grid grids
  implicit none

  integer, intent(in) :: NY,NX,N,M,N1,N2,N3,N4,N5,N6
  integer :: K,LL,ids,idg,idom
!
!     NET GAS FLUX
!
  IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    DO idg=idg_beg,idg_NH3
      Gas_AdvDif_FlxMM_vr(idg,N3,N2,N1)=Gas_AdvDif_FlxMM_vr(idg,N3,N2,N1)+Gas_AdvDif_FlxMM_3D(idg,N,N3,N2,N1) &
        -Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4)
    ENDDO
  ENDIF
  end subroutine NetTracerFlowXSoilPoresMM

! ----------------------------------------------------------------------

  subroutine GasDiffusionMM(M,N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: trc_gasc1(idg_beg:idg_end)
  real(r8) :: trc_gasc2(idg_beg:idg_end)
  character(len=*), parameter :: subname='GasDiffusionMM'
  real(r8) :: CNDC1   !conductance in source grid
  real(r8) :: CNDC2   !conductance in destination grid
  real(r8) :: DFLG2   !air-filled porosity effect on gaseous diffusivity in source
  real(r8) :: DFLGL   !air-filled porosity effect on gaseous diffusivity in destination
  integer :: idg
  real(r8):: GasDifuscoefMM_3D

  call PrintInfo('beg '//subname)
  !
  !     GASEOUS DIFFUSIVITIES
  !
  DFLG2=2.0_r8*AZMAX1(FracAirFilledSoilPoreM_vr(M,N3,N2,N1))*POROQ*FracAirFilledSoilPoreM_vr(M,N3,N2,N1)/POROS_vr(N3,N2,N1) &
    *AREA(N,N3,N2,N1)/DLYR_3D(N,N3,N2,N1)

  DFLGL=2.0_r8*AZMAX1(FracAirFilledSoilPoreM_vr(M,N6,N5,N4))*POROQ*FracAirFilledSoilPoreM_vr(M,N6,N5,N4)/POROS_vr(N6,N5,N4) &
    *AREA(N,N6,N5,N4)/DLYR_3D(N,N6,N5,N4)
  !
  DO idg=idg_beg,idg_NH3
    !     GASOUS CONDUCTANCES
    CNDC1             = DFLG2*GasDifctScaledMM_vr(idg,N3,N2,N1)
    CNDC2             = DFLGL*GasDifctScaledMM_vr(idg,N6,N5,N4)
    GasDifuscoefMM_3D = (CNDC1*CNDC2)/(CNDC1+CNDC2)
    !
    !diffusion flux
    !
    trc_gasc1(idg)                      = AZMAX1(trc_gasml2_vr(idg,N3,N2,N1)/VLsoiAirPM_vr(M,N3,N2,N1))
    trc_gasc2(idg)                      = AZMAX1(trc_gasml2_vr(idg,N6,N5,N4)/VLsoiAirPM_vr(M,N6,N5,N4))
    Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4) = Gas_AdvDif_FlxMM_3D(idg,N,N6,N5,N4)+GasDifuscoefMM_3D*(trc_gasc1(idg)-trc_gasc2(idg))
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine GasDiffusionMM
! ----------------------------------------------------------------------
  subroutine TopSoilGasDifussionMM(M,NY,NX)
  !
  !Description:
  !Gas diffusion between atmosphere and soil. (>0 into soil)
  implicit none
  integer, intent(in) :: M,NY,NX
  character(len=*), parameter :: subname='TopSoilGasDifussionMM'
  real(r8) :: DFLG2
  real(r8) :: trcg_cl2
  real(r8) :: DGQ_cef,GasDifuscoefMM_3D
  integer  :: idg

!     GASEOUS DIFFUSIVITIES
!
  call PrintInfo('beg '//subname)
  DFLG2=AZMAX1(FracAirFilledSoilPoreM_vr(M,NU(NY,NX),NY,NX))*POROQ &
    *FracAirFilledSoilPoreM_vr(M,NU(NY,NX),NY,NX)/POROS_vr(NU(NY,NX),NY,NX) &
    *AREA(3,NU(NY,NX),NY,NX)/AMAX1(ZERO2,DLYR_3D(3,NU(NY,NX),NY,NX))

  DO idg=idg_beg,idg_NH3
    GasDifuscoefMM_3D=DFLG2*GasDifctScaledMM_vr(idg,NU(NY,NX),NY,NX)
!
!     SURFACE GAS CONCENTRATIONS
!
    trcg_cl2=AZMAX1(trc_gasml2_vr(idg,NU(NY,NX),NY,NX)/VLsoiAirPM_vr(M,NU(NY,NX),NY,NX))
!
!     EQUILIBRIUM CONCENTRATIONS AT SOIL SURFACE AT WHICH
!     GASEOUS DIFFUSION THROUGH SOIL SURFACE LAYER = GASEOUS
!     DIFFUSION THROUGH ATMOSPHERE BOUNDARY LAYER CALCULATED
!     FROM GASEOUS DIFFUSIVITY AND BOUNDARY LAYER CONDUCTANCE
!
    DGQ_cef=GasDifuscoefMM_3D*PARGas_CefMM(idg,NY,NX)/(GasDifuscoefMM_3D+PARGas_CefMM(idg,NY,NX))

    Gas_AdvDif_FlxMM_3D(idg,3,NU(NY,NX),NY,NX)=Gas_AdvDif_FlxMM_3D(idg,3,NU(NY,NX),NY,NX)+DGQ_cef*(AtmGasCgperm3(idg,NY,NX)-trcg_cl2)
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine TopSoilGasDifussionMM

!------------------------------------------------------------------------------------------
  subroutine FastUpdateStateVarsMM(I,J,NHW, NHE, NVN, NVS,dpscal,pscal)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: NHW, NHE, NVN, NVS
  real(r8),intent(in) :: dpscal(idg_beg:idg_end)
  real(r8),optional,intent(inout) :: pscal(idg_beg:idg_end)
  
  if(present(pscal))then
    pscal=1.001_r8
    call UpdateStateVarsMM(I,J,NHW, NHE, NVN, NVS,dpscal,pscal)
    pscal=pscal*0.9999_r8
  else
    call UpdateStateVarsMM(I,J,NHW, NHE, NVN, NVS,dpscal)
  endif  
  end subroutine FastUpdateStateVarsMM
!------------------------------------------------------------------------------------------
  subroutine UpdateStateVarsMM(I,J,NHW, NHE, NVN, NVS,dpscal,pscal1)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: NHW, NHE, NVN, NVS
  real(r8),intent(in) :: dpscal(idg_beg:idg_end)
  real(r8), optional, intent(inout) :: pscal1(idg_beg:idg_end)

  character(len=*), parameter :: subname='UpdateStateVarsMM'
  integer :: NY,NX
  integer :: L,ids,idg
  real(r8) :: pscal(idg_beg:idg_end)  
  real(r8) :: flux

  call PrintInfo('beg '//subname)
  if(present(pscal1))then
    pscal=pscal1
  else
    pscal=1.001_r8
  endif  
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      !litter layer
      DO idg=idg_beg,idg_NH3
        if(dpscal(idg)>tiny_p)then
          flux = -RBGCSinkGasMM_vr(idg,0,NY,NX)+RGas_Disol_FlxMM_vr(idg,0,NY,NX)
          flux = flux*dpscal(idg)
          call get_flux_scalar(trcs_solml2_vr(idg,0,NY,NX),flux,trcs_solml_vr(idg,0,NY,NX),pscal(idg))
        endif
      ENDDO
      ! 
      !in the soil      
      !
      DO L=NU(NY,NX),NL(NY,NX)
        !other gases are taken from aqueous phases
        DO idg=idg_beg,idg_NH3-1
          if(dpscal(idg)>tiny_p)then
            flux = -RBGCSinkGasMM_vr(idg,L,NY,NX)+RGas_Disol_FlxMM_vr(idg,L,NY,NX)
            flux=flux*dpscal(idg)
            call get_flux_scalar(trcs_solml2_vr(idg,L,NY,NX),flux, trcs_solml_vr(idg,L,NY,NX),pscal(idg))
          endif
        ENDDO

        DO idg=idg_NH3,idg_NH3B
          if(dpscal(idg)>tiny_p)then
            flux=RGas_Disol_FlxMM_vr(idg,L,NY,NX)
            flux=flux*dpscal(idg)
            call get_flux_scalar(trcs_solml2_vr(idg,L,NY,NX),flux, trcs_solml_vr(idg,L,NY,NX),pscal(idg))
          endif
        ENDDO

        DO idg=idg_beg,idg_NH3-1
          if(dpscal(idg)>tiny_p)then
            flux=-RGas_Disol_FlxMM_vr(idg,L,NY,NX)+Gas_AdvDif_FlxMM_vr(idg,L,NY,NX)
            flux=flux*dpscal(idg)
            call get_flux_scalar(trc_gasml2_vr(idg,L,NY,NX),flux, trcg_gasml_vr(idg,L,NY,NX),pscal(idg))
          endif
        ENDDO

        !gaseous NH3 is consumed through geochemistry
        idg=idg_NH3
        if(dpscal(idg)>tiny_p)then
          flux=-RGas_Disol_FlxMM_vr(idg_NH3B,L,NY,NX)-RBGCSinkGasMM_vr(idg_NH3,L,NY,NX)+Gas_AdvDif_FlxMM_vr(idg,L,NY,NX)   
          flux=flux*dpscal(idg) 
          call get_flux_scalar(trc_gasml2_vr(idg_NH3,L,NY,NX),flux, trcg_gasml_vr(idg_NH3,L,NY,NX),pscal(idg))  
        endif
      ENDDO
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  if(present(pscal1))pscal1=pscal  
  end subroutine UpdateStateVarsMM
!------------------------------------------------------------------------------------------

  subroutine TracerFlowXGridsMM(I,J,M,NHE,NHW,NVS,NVN)
!
! DESCRIPTION
! exchanges tracers within (gaseous vs aqueous phase) and between
! grid cells.
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHE,NHW,NVS,NVN

  character(len=*), parameter :: subname='TracerFlowXGridsMM'
  real(r8) :: VFLW
  real(r8) :: VOLH2A,VOLH2B
  real(r8) :: VLWatMacPS,VOLWT

  integer :: iDisableEbu,N,L,K,LL,NY,NX
  integer :: N1,N2,N3  !source grid index
  integer :: N4,N5,N6  !dest grid index

! begin_execution
  call PrintInfo('beg '//subname)

!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      iDisableEbu=ifalse
      D125: DO L=1,NL(NY,NX)
        !source
        N1=NX;N2=NY;N3=L    
        IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))then
          call GasDissolutionMM(M,N,N3,N2,N1)        
        endif  
        !
        !     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
        !
        D120: DO N=FlowDirIndicator_col(N2,N1),3

          IF(N.EQ.iWestEastDirection)THEN
            IF(NX.EQ.NHE)THEN  !skip eastern boundary
              cycle
            ELSE
              N4 = NX+1;N5 = NY;N6 = L
            ENDIF
          ELSEIF(N.EQ.iNorthSouthDirection)THEN            
            IF(NY.EQ.NVS)THEN  !skip southern boundary
              cycle
            ELSE
              N4 = NX;N5 = NY+1;N6 = L
            ENDIF
          ELSEIF(N.EQ.iVerticalDirection)THEN            
            IF(L.EQ.NL(NY,NX))THEN  !skip bottom boundary
              cycle
            ELSE
              N4 = NX;N5 = NY;N6 = L+1
            ENDIF
          ENDIF

          DO LL=N6,NL(NY,NX)
            IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
              N6=LL
              exit
            ENDIF
          ENDDO
          !
          !     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS FROM
          !     WATER CONTEnsolutes AND WATER FLUXES 'WaterFlow2SoilMM_3D' FROM 'WATSUB'
          !
          IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
            IF(N3.GE.NUM(N2,N1) .AND. N6.GE.NUM(N5,N4) .AND. N3.LE.NL(N2,N1) .AND. N6.LE.NL(N5,N4))THEN
              ! 
              call GasTransportMM(M,N,N1,N2,N3,N4,N5,N6)
            ENDIF           
          ENDIF
        ENDDO D120
      ENDDO D125
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine TracerFlowXGridsMM

end module TranspNoSaltFastMod