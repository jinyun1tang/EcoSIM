module EcoSIMAPI
  USE EcoSIMCtrlDataType
  use timings      , only : start_timer, end_timer
  use ErosionMod   , only : erosion
  use Hour1Mod     , only : hour1
  use RedistMod    , only : redist
  use GeochemAPI   , only : soluteModel
  use PlantAPI     , only : PlantModel
  use MicBGCAPI    , only : MicrobeModel, MicAPI_Init, MicAPI_cleanup
  use TrnsfrMod    , only : trnsfr
  use TrnsfrsMod   , only : trnsfrs
  use EcoSIMCtrlMod, only : lverb
  use WatsubMod    , only : watsub
implicit none

  character(len=*),private, parameter :: mod_filename = __FILE__
  public :: Run_EcoSIM_one_step

contains

  subroutine Run_EcoSIM_one_step(I,J,NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  real(r8) :: t1

  if(lverb)WRITE(*,334)'HOUR1'
  call start_timer(t1)
  CALL HOUR1(I,J,NHW,NHE,NVN,NVS)
  call end_timer('HOUR1',t1)
  !
  !   CALCULATE SOIL ENERGY BALANCE, WATER AND HEAT FLUXES IN 'WATSUB'
  !
  if(lverb)WRITE(*,334)'WAT'
  call start_timer(t1)
  CALL WATSUB(I,J,NHW,NHE,NVN,NVS)
  call end_timer('WAT',t1)
  !
  !   CALCULATE SOIL BIOLOGICAL TRANSFORMATIONS IN 'NITRO'
  !
  if(lverb)WRITE(*,334)'NIT'
  call start_timer(t1)
  CALL MicrobeModel(I,J,NHW,NHE,NVN,NVS)
  call end_timer('NIT',t1)
  !
  !   UPDATE PLANT biogeochemistry
  !
  if(lverb)WRITE(*,334)'PlantModel'
  call PlantModel(I,J,NHW,NHE,NVN,NVS)
  !
  !
  !   CALCULATE SOLUTE EQUILIBRIA IN 'SOLUTE'
  !
  if(lverb)WRITE(*,334)'SOL'
  call start_timer(t1)
  CALL soluteModel(I,J,NHW,NHE,NVN,NVS)
  call end_timer('SOL',t1)
  !
  !   CALCULATE GAS AND SOLUTE FLUXES IN 'TRNSFR'
  !
  if(lverb)WRITE(*,334)'TRN'
  !    if(I>=170)print*,TKS(0,NVN,NHW)
  call start_timer(t1)
  CALL TRNSFR(I,J,NHW,NHE,NVN,NVS)
  call end_timer('TRN',t1)
  !
  !   CALCULATE ADDITIONAL SOLUTE FLUXES IN 'TRNSFRS' IF SALT OPTION SELECTED
  !
  if(lverb)WRITE(*,334)'TRNS'
  !    if(I>=170)print*,TKS(0,NVN,NHW)
  call start_timer(t1)
  CALL TRNSFRS(I,J,NHW,NHE,NVN,NVS)
  call end_timer('TRNSFRS',t1)
  !
  !   CALCULATE SOIL SEDIMENT TRANSPORT IN 'EROSION'
  !
  if(lverb)WRITE(*,334)'EROSION'
  !    if(I>=170)print*,TKS(0,NVN,NHW)
  call start_timer(t1)
  CALL EROSION(I,J,NHW,NHE,NVN,NVS)
  call end_timer('EROSION',t1)
  !
  !   UPDATE ALL SOIL STATE VARIABLES FOR WATER, HEAT, GAS, SOLUTE
  !   AND SEDIMENT FLUXES IN 'REDIST'
  !
  if(lverb)WRITE(*,334)'RED'
  !    if(I>=170)print*,TKS(0,NVN,NHW)
  call start_timer(t1)
  CALL REDIST(I,J,NHW,NHE,NVN,NVS)
  call end_timer('RED',t1)

334   FORMAT(A8)

  end subroutine Run_EcoSIM_one_step


!------------------------------------------------------------------------

  subroutine SetMesh(NHW,NVN,NHE,NVS)

  use EcoSIMConfig, only : column_mode
  use EcoSIMCtrlMod,only : grid_file_in
  use ncdio_pio
  use netcdf
  USE fileUtil, ONLY : iulog
  use abortutils, only : endrun
  use GridConsts, only : JX,JY,JZ,JH,JV,JD,bounds,JP
!  set up the landscape rectangular mesh
!  beginning(NHW,NVN)
!  o--------------------------x
!  |                          |
!  |                          |
!  |                          |
!  |                          |
!  x--------------------------o
!                             end (NHE,NVS)

  implicit none
  integer, intent(out) :: NHW   !upper corner x index
  integer, intent(out) :: NVN   !upper corner y index
  integer, intent(out) :: NHE   !lower corner x index
  integer, intent(out) :: NVS   !lower corner y index
  integer :: nextra_grid
  INTEGER :: NZ,NY,NX,ic,ip
  type(file_desc_t) :: grid_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar

  call ncd_pio_openfile(grid_nfid, grid_file_in, ncd_nowrite)

  call check_var(grid_nfid, 'NHW', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NHW in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NHW), 'in '//trim(mod_filename))

  call check_var(grid_nfid, 'NHE', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NHE in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NHE), 'in '//trim(mod_filename))

  call check_var(grid_nfid, 'NVN', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NVN in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NVN), 'in '//trim(mod_filename))

  call check_var(grid_nfid, 'NVS', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NVS in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NVS), 'in '//trim(mod_filename))

  call ncd_pio_closefile(grid_nfid)

  bounds%NHW =NHW
  bounds%NVN =NVN
  bounds%NHE =NHE
  bounds%NVS =NVS
  bounds%ncols=((NHE-NHW)+1)*((NVS-NVN)+1)
  bounds%npfts=bounds%ncols*JP
  allocate(bounds%icol((NVS-NVN)+1,(NHE-NHW)+1))
  allocate(bounds%ipft(JP,(NVS-NVN)+1,(NHE-NHW)+1))

  ic=0;ip=0
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      ic=ic+1
      bounds%icol(NY,NX)=ic
      DO NZ=1,JP
        ip=ip+1
        bounds%ipft(NZ,NY,NX)=ip
      ENDDO
    ENDDO
  ENDDO

  nextra_grid=1
  if(column_mode)nextra_grid=0
  JX=(NHE-NHW)+1+nextra_grid
  JY=(NVS-NVN)+1+nextra_grid
  JZ=14
  JH=JX+nextra_grid
  JV=JY+nextra_grid
  JD=JZ+1
  write(iulog,*)'grid size'
  write(iulog,*)'JX=',JX,'JY=',JY,'JZ=',JZ
  end subroutine SetMesh
end module EcoSIMAPI
