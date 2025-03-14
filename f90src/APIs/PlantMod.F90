module PlantMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use grosubsMod,        only: GrowPlants
  use PlantPhenolMod,    only: hfuncs
  use UptakesMod,        only: RootUptakes
  use PlantDisturbMod,   only: PrepLandscapeGrazing
  use PlantMgmtDataType, only: NP
  use MiniMathMod,       only: fixEXConsumpFlux
  use EcoSIMCtrlMod,     only: lverb
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

      call EnterPlantBalance(I,J,NP(NY,NX))

!   UPDATE PLANT PHENOLOGY IN 'HFUNC'
!     zero out plant hourly fluxes
      call ZeroGrosub()  

      if(lverb)WRITE(*,333)'HFUNC'

!     DO NZ=1,NP(NY,NX)
!       call SumPlantBiom(I,J,NZ,'bfHFUNCS')
!     ENDDO
      !Phenological update, determine living/active branches      
      CALL HFUNCs(I,J)

!      DO NZ=1,NP(NY,NX)
!        call SumPlantBiom(I,J,NZ,'bfUPTAKES')
!      ENDDO
!      if(I==140 .and. J>=20)write(116,*)'bfrootupk',I*1000+J        
!      call SumPlantRootGas(I,J)

      !Predict uptake fluxes of nutrients and O2
      if(lverb)write(*,*)'uptake'
      CALL ROOTUPTAKES(I,J)
!      if(I==140 .and. J>=20)write(116,*)'afrootupk',I*1000+J         
!      call SumPlantRootGas(I,J)

!      DO NZ=1,NP(NY,NX)
!        call SumPlantBiom(I,J,NZ,'bfGROSUBS')
!      ENDDO

      !Do growth of active branches and roots
      if(lverb)write(*,*)'grosub'
      CALL GROWPLANTS(I,J)

      if(lverb)write(*,*)'EXTRACT'

      !aggregate varaibles
      CALL EXTRACTs(I,J)
!      if(I==140 .and. J>=20)write(116,*)'afextract'        

      call ExitPlantBalance(I,J,NP(NY,NX))

      call PlantAPIRecv(I,J,NY,NX)

      call ApplyRootUptake2GasTracers(I,J,NY,NX)
    ENDDO
  ENDDO
  PlantElemntStoreLandscape(:)=plt_site%PlantElemntStoreLandscape(:)


  end subroutine PlantModel
!------------------------------------------------------------------------------------------

  subroutine PlantCanopyRadsModel(I,J,NY,NX,DPTH0)
  use CanopyCondsMod
  use CanopyDataType, only : CanopyLeafArea_lpft
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(in) :: DPTH0

  call PlantAPICanMSend(NY,NX)

  call CanopyConditionModel(I,J,DPTH0)

  call PlantAPICanMRecv(NY,NX)

  end subroutine PlantCanopyRadsModel
!------------------------------------------------------------------------------------------

  subroutine ApplyRootUptake2GasTracers(I,J,NY,NX)

  implicit none
  integer,intent(in) :: I,J,NY,NX
  integer :: L,idg
  real(r8):: dmass

  DO L=NU(NY,NX),NL(NY,NX)
    do idg=idg_beg,idg_NH3            
      if (trcs_solml_vr(idg,L,NY,NX)>trcs_plant_uptake_vr(idg,L,NY,NX))then
        trcs_solml_vr(idg,L,NY,NX)=trcs_solml_vr(idg,L,NY,NX)-trcs_plant_uptake_vr(idg,L,NY,NX)
      else
        dmass=trcs_plant_uptake_vr(idg,L,NY,NX)-trcs_solml_vr(idg,L,NY,NX)
        trcs_solml_vr(idg,L,NY,NX)=0._r8        
        if(dmass > trcg_gasml_vr(idg,L,NY,NX))then
          dmass=dmass-trcg_gasml_vr(idg,L,NY,NX)
          trcg_gasml_vr(idg,L,NY,NX)=0._r8
          trcs_plant_uptake_vr(idg,L,NY,NX)=trcs_plant_uptake_vr(idg,L,NY,NX)-dmass
        else
          trcg_gasml_vr(idg,L,NY,NX)=trcg_gasml_vr(idg,L,NY,NX)-dmass
        endif
      endif
    enddo

    idg=idg_NH3B
    if (trcs_solml_vr(idg,L,NY,NX)>trcs_plant_uptake_vr(idg,L,NY,NX))then
      trcs_solml_vr(idg,L,NY,NX)=trcs_solml_vr(idg,L,NY,NX)-trcs_plant_uptake_vr(idg,L,NY,NX)
    else
      dmass=trcs_plant_uptake_vr(idg,L,NY,NX)-trcs_solml_vr(idg,L,NY,NX)
      trcs_solml_vr(idg,L,NY,NX)=0._r8        
      !deficit greater than gas concentration
      if(dmass > trcg_gasml_vr(idg_NH3,L,NY,NX))then
        dmass=dmass-trcg_gasml_vr(idg_NH3,L,NY,NX)
        trcg_gasml_vr(idg_NH3,L,NY,NX)=0._r8
        trcs_plant_uptake_vr(idg,L,NY,NX)=trcs_plant_uptake_vr(idg,L,NY,NX)-dmass
      else
        trcg_gasml_vr(idg_NH3,L,NY,NX)=trcg_gasml_vr(idg_NH3,L,NY,NX)-dmass
      endif
    endif
  ENDDO

  end subroutine ApplyRootUptake2GasTracers
end module PlantMod
