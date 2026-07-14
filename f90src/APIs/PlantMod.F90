module PlantMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8,yearIJ_type
  use grosubsMod,        only: GrowPlants
  use PlantPhenolMod,    only: PhenologyUpdate
  use UptakesMod,        only: RootUptakes
  use PlantDisturbMod,   only: PrepLandscapeGrazing
  use PlantMgmtDataType, only: NP_col
  use MiniMathMod,       only: fixEXConsumpFlux
  use EcoSIMCtrlMod,     only: lverb, ldo_sp_mode
  use PlantDebugMod,     only: PrintRootTracer
  use LitterFallMod,     only: ReSeedPlants
  use DebugToolMod,      only: PrintInfo
  use PlantAPI4Uptake
  use TracerIDMod
  use GridDataType
  use PlantDataRateType
  use SoilBGCDataType
  use EcoSimSumDataType
  use PlantAPIData  
  use PlantAPI
  use PlantCanAPI
  use ExtractsMod
  use PlantBalMod
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: PlantModel
  public :: PlantCanopyRadsModel

  contains

  subroutine PlantModel(yearIJ,NHW,NHE,NVN,NVS)
  !
  !run the plant biogeochemistry model
  !
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ  
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: t1,tvegE(NumPlantChemElms)
  integer :: NY,NX,NZ
  character(len=*), parameter :: subname='PlantModel'
333   FORMAT(A8)

  call PrintInfo('beg '//subname)
  if(.not.ldo_sp_mode)then
    call PrepLandscapeGrazing(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)
    plt_site%PlantElemntStoreLandscape(:)=PlantElemntStoreLandscape(:)
  endif

  DO NX=NHW,NHE
    DO NY=NVN,NVS

      if(ldo_sp_mode)then
        !do prescribed phenolgoy
        call PlantUptakeAPISend(yearIJ%I,yearIJ%J,NY,NX)        

        CALL ROOTUPTAKES(yearIJ)
        
        call extracts(yearIJ%I,yearIJ%J)
        
        call PlantUPtakeAPIRecv(yearIJ%I,yearIJ%J,NY,NX)
      else
      
        call  PlantAPISend(yearIJ%I,yearIJ%J,NY,NX)

        call EnterPlantBalance(yearIJ,NP_col(NY,NX))

        !Phenological update, determine living/active branches      
        CALL PhenologyUpdate(yearIJ%I,yearIJ%J)

        !Predict uptake fluxes of nutrients and O2        
        CALL ROOTUPTAKES(yearIJ)

        !Do growth of active branches and roots
        CALL GROWPLANTS(yearIJ)

        !aggregate varaibles
        CALL EXTRACTs(yearIJ%I,yearIJ%J)

        DO NZ=1,NP_col(NY,NX)
          !
          Call ReSeedPlants(yearIJ%I,yearIJ%J,NZ)
        ENDDO

        call ExitPlantBalance(yearIJ,NP_col(NY,NX))

        call PlantAPIRecv(yearIJ%I,yearIJ%J,NY,NX)
        
      endif

    ENDDO
  ENDDO
  
  PlantElemntStoreLandscape(:)=plt_site%PlantElemntStoreLandscape(:)
  call PrintInfo('end '//subname)

  end subroutine PlantModel
!------------------------------------------------------------------------------------------
!  function test_active_plant(NY,NX)result(activeroot)
  !
  !Description:
  !Check the existence of active plants
!  implicit none
!  integer, intent(in) :: NY,NX
!  logical :: activeroot
!  integer :: NZ

!  activeroot=.false.
!  DO NZ=1,NP_col(NY,NX)
!    activeroot=IsPlantActive_pft(NZ,NY,NX).EQ.iTrue .and. PlantPopuLive_pft(NZ,NY,NX)>ZEROS(NY,NX)
!    if(activeroot)return
!  ENDDO
!  end function test_active_plant

!------------------------------------------------------------------------------------------

  subroutine PlantCanopyRadsModel(I,J,NY,NX,DepthSurfWatIce)
  use SurfaceRadiationMod
  use CanopyDataType, only : CanopyLeafArea_lnode
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(in) :: DepthSurfWatIce

  call PlantAPICanMSend(NY,NX)

  call CanopyConditionModel(I,J,DepthSurfWatIce)

  call PlantAPICanMRecv(NY,NX)

  end subroutine PlantCanopyRadsModel

end module PlantMod
