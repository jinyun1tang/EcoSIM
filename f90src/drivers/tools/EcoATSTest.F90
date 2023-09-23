program EcoATSTest
  use ATSCPLMod

  implicit none

  ! Initialize ATS and EcoSIM sizes
  type (BGCSizes) :: sizes
  write(*,*) "In driver"
  write(*,*) "Setting Sizes"

  call SetBGCSizes(sizes)

  write(*,*) "Init EcoSIM"
  ! Initialize EcoSIM
  call Init_EcoSIM(sizes)

  write(*,*) "Transferring data"
  ! Transfer data from ATS to EcoSIM
  !integer :: ncol = sizes%ncells_per_col_
  !type (BGCState) :: state
  !type (BGCProperties) :: props
  !call ATS2EcoSIMData(ncol, state, props, sizes)

  ! Run EcoSIM for one time step
  !call Run_EcoSIM_one_step()

  !write(*,*) "Copy back"
  ! Transfer data from EcoSIM to ATS
  !call EcoSIM2ATSData(ncol, state, sizes)

  ! Clean up and finalize
  write(*,*) "EcoATSTest completed successfully."

end program EcoATSTest
