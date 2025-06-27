module PlantMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use grosubsMod,        only: GrowPlants
  use PlantPhenolMod,    only: hfuncs
  use UptakesMod,        only: RootUptakes
  use PlantDisturbMod,   only: PrepLandscapeGrazing
  use PlantMgmtDataType, only: NP_col
  use MiniMathMod,       only: fixEXConsumpFlux
  use EcoSIMCtrlMod,     only: lverb,ldo_sp_mode
  use PlantDebugMod,     only: PrintRootTracer
  use LitterFallMod,     only: ReSeedPlants
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

  subroutine PlantModel(I,J,NHW,NHE,NVN,NVS)
  !
  !run the plant biogeochemistry model
  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: t1
  integer :: NY,NX,NZ

333   FORMAT(A8)

  call PrepLandscapeGrazing(I,J,NHW,NHE,NVN,NVS)

  plt_site%PlantElemntStoreLandscape(:)=PlantElemntStoreLandscape(:)
  DO NX=NHW,NHE
    DO NY=NVN,NVS
  
      call  PlantAPISend(I,J,NY,NX)

      if(ldo_sp_mode)then
        !do prescribed phenolgoy
        CALL ROOTUPTAKES(I,J)
      else
        call EnterPlantBalance(I,J,NP_col(NY,NX))

        !   UPDATE PLANT PHENOLOGY IN 'HFUNC'
        !     zero out plant hourly fluxes
        call ZeroGrosub()  

        if(lverb)WRITE(*,333)'HFUNC'

        DO NZ=1,NP_col(NY,NX)
          call PrintRootTracer(I,J,NZ,'BEEfhfunc')
  !       call SumPlantBiom(I,J,NZ,'bfHFUNCS')
        ENDDO

        !Phenological update, determine living/active branches      
        CALL HFUNCs(I,J)

        DO NZ=1,NP_col(NY,NX)
          call PrintRootTracer(I,J,NZ,'afhfunc')
  !        call SumPlantBiom(I,J,NZ,'bfUPTAKES')
        ENDDO

  !      call SumPlantRootGas(I,J)

        !Predict uptake fluxes of nutrients and O2
        if(lverb)write(*,*)'uptake'
        CALL ROOTUPTAKES(I,J)

  !      call SumPlantRootGas(I,J)

  !      DO NZ=1,NP_col(NY,NX)
  !        call SumPlantBiom(I,J,NZ,'bfGROSUBS')
  !      ENDDO

        !Do growth of active branches and roots
        if(lverb)write(*,*)'grosub'
        CALL GROWPLANTS(I,J)

        if(lverb)write(*,*)'EXTRACT'

        !aggregate varaibles
        CALL EXTRACTs(I,J)
  !      if(I==140 .and. J>=20)write(116,*)'afextract'        

        DO NZ=1,NP_col(NY,NX)
          Call ReSeedPlants(I,J,NZ)
        ENDDO

        call ExitPlantBalance(I,J,NP_col(NY,NX))
      endif
      call PlantAPIRecv(I,J,NY,NX)

    ENDDO
  ENDDO
  PlantElemntStoreLandscape(:)=plt_site%PlantElemntStoreLandscape(:)


  end subroutine PlantModel
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
