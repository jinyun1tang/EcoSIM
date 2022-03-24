module MicrobialDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__


  integer :: JG

  real(r8),allocatable :: OMC(:,:,:,:,:,:,:)    !microbial biomass C	[g d-2]
  real(r8),allocatable :: OMN(:,:,:,:,:,:,:)    !microbial biomass N	[g d-2]
  real(r8),allocatable :: OMP(:,:,:,:,:,:,:)    !microbial biomass P	[g d-2]
  real(r8),allocatable :: ROXYS(:,:,:,:,:,:)    !aqueous O2 demand	[g d-2 h-1]
  real(r8),allocatable :: ROQCS(:,:,:,:,:,:)    !net microbial DOC flux	[g d-2 h-1]
  real(r8),allocatable :: ROQAS(:,:,:,:,:,:)    !net microbial acetate flux	[g d-2 h-1]
  real(r8),allocatable :: CNOMC(:,:,:,:)        !maximum microbial N:C
  real(r8),allocatable :: CPOMC(:,:,:,:)        !maximum microbial P:C
  real(r8),allocatable :: RINHO(:,:,:,:,:,:)    !microbial NH4 demand in soil	[g d-2 h-1]
  real(r8),allocatable :: RINOO(:,:,:,:,:,:)    !microbial NO3 demand in soil	[g d-2 h-1]
  real(r8),allocatable :: RIPOO(:,:,:,:,:,:)    !microbial PO4 demand in soil	[g d-2 h-1]
  real(r8),allocatable :: RINHOR(:,:,:,:,:)     !microbial NH4 demand in surface litter	[g d-2 h-1]
  real(r8),allocatable :: RIPOOR(:,:,:,:,:)     !microbial PO4 demand in surface litter	[g d-2 h-1]
  real(r8),allocatable :: RINOOR(:,:,:,:,:)     !microbial NO3 demand in surface litter	[g d-2 h-1]
  real(r8),allocatable :: RVMX4(:,:,:,:,:,:)    !total microbial NH4 uptake non-band unconstrained by NH4	[g d-2 h-1]
  real(r8),allocatable :: RVMX3(:,:,:,:,:,:)    !total microbial NO3 uptake non-band unconstrained by NO3	[g d-2 h-1]
  real(r8),allocatable :: RVMX2(:,:,:,:,:,:)    !total microbial NO2 uptake non-band unconstrained by NO2	[g d-2 h-1]
  real(r8),allocatable :: RVMB4(:,:,:,:,:,:)    !total microbial NH4 uptake non-band unconstrained by NH4	[g d-2 h-1]
  real(r8),allocatable :: RVMB3(:,:,:,:,:,:)    !total microbial NO3 uptake band unconstrained by NO3	[g d-2 h-1]
  real(r8),allocatable :: RVMB2(:,:,:,:,:,:)    !total microbial NO2 uptake band unconstrained by NH4	[g d-2 h-1]
  real(r8),allocatable :: RVMX1(:,:,:,:,:,:)    !total microbial N2O uptake unconstrained by N2O	[g d-2 h-1]
  real(r8),allocatable :: RINHB(:,:,:,:,:,:)    !microbial NH4 immobilization (+ve) - mineralization (-ve) band	[g d-2 h-1]
  real(r8),allocatable :: RINOB(:,:,:,:,:,:)    !microbial NO3 immobilization (+ve) - mineralization (-ve) band	[g d-2 h-1]
  real(r8),allocatable :: RIPBO(:,:,:,:,:,:)    !substrate-unlimited H2PO4 mineraln-immobiln
  real(r8),allocatable :: RIPO1(:,:,:,:,:,:)    !substrate-unlimited HPO4 immobilization
  real(r8),allocatable :: RIPB1(:,:,:,:,:,:)    !substrate-unlimited HPO4 mineraln-immobiln
  real(r8),allocatable :: RIPO1R(:,:,:,:,:)     !substrate-unlimited HPO4 immobilization
  real(r8),allocatable :: OMCER(:,:,:,:,:,:,:)  !microbial C  erosion 	[g d-2 h-1]
  real(r8),allocatable :: OMNER(:,:,:,:,:,:,:)  !microbial N  erosion 	[g d-2 h-1]
  real(r8),allocatable :: OMPER(:,:,:,:,:,:,:)  !microbial P  erosion 	[g d-2 h-1]
  real(r8),allocatable :: OMCI(:,:)             !initializes microbial biomass
  real(r8) :: FL(2)                             !allocation to microbial kinetic fractions
  private :: InitAllocate

  contains

  subroutine InitMicrobialData(nguilds)

  implicit none
  integer, intent(in) :: nguilds



  JG=nguilds
  call InitAllocate()

  end subroutine InitMicrobialData

  subroutine InitAllocate

  implicit none


  allocate(OMC(3,JG,7,0:5,0:JZ,JY,JX))
  allocate(OMN(3,JG,7,0:5,0:JZ,JY,JX))
  allocate(OMP(3,JG,7,0:5,0:JZ,JY,JX))
  allocate(ROXYS(JG,7,0:5,0:JZ,JY,JX))
  allocate(ROQCS(JG,7,0:4,0:JZ,JY,JX))
  allocate(ROQAS(JG,7,0:4,0:JZ,JY,JX))
  allocate(CNOMC(3,JG,7,0:5))
  allocate(CPOMC(3,JG,7,0:5))
  allocate(RINHO(JG,7,0:5,0:JZ,JY,JX))
  allocate(RINOO(JG,7,0:5,0:JZ,JY,JX))
  allocate(RIPOO(JG,7,0:5,0:JZ,JY,JX))
  allocate(RINHOR(JG,7,0:5,JY,JX))
  allocate(RIPOOR(JG,7,0:5,JY,JX))
  allocate(RINOOR(JG,7,0:5,JY,JX))
  allocate(RVMX4(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMX3(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMX2(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMB4(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMB3(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMB2(JG,7,0:5,0:JZ,JY,JX))
  allocate(RVMX1(JG,7,0:5,0:JZ,JY,JX))
  allocate(RINHB(JG,7,0:5,0:JZ,JY,JX))
  allocate(RINOB(JG,7,0:5,0:JZ,JY,JX))
  allocate(RIPBO(JG,7,0:5,0:JZ,JY,JX))
  allocate(RIPO1(JG,7,0:5,0:JZ,JY,JX))
  allocate(RIPB1(JG,7,0:5,0:JZ,JY,JX))
  allocate(RIPO1R(JG,7,0:5,JY,JX))
  allocate(OMCER(3*JG,7,0:5,2,2,JV,JH))
  allocate(OMNER(3*JG,7,0:5,2,2,JV,JH))
  allocate(OMPER(3*JG,7,0:5,2,2,JV,JH))
  allocate(OMCI(3*JG,0:4))

  end subroutine InitAllocate
!----------------------------------------------------------------------------------------------

  subroutine DestructMicrobialData
  implicit none

  deallocate(OMC)
  deallocate(OMN)
  deallocate(OMP)
  deallocate(ROXYS)
  deallocate(ROQCS)
  deallocate(ROQAS)
  deallocate(CNOMC)
  deallocate(CPOMC)
  deallocate(RINHO)
  deallocate(RINOO)
  deallocate(RIPOO)
  deallocate(RINHOR)
  deallocate(RIPOOR)
  deallocate(RINOOR)
  deallocate(RVMX4)
  deallocate(RVMX3)
  deallocate(RVMX2)
  deallocate(RVMB4)
  deallocate(RVMB3)
  deallocate(RVMB2)
  deallocate(RVMX1)
  deallocate(RINHB)
  deallocate(RINOB)
  deallocate(RIPBO)
  deallocate(RIPO1)
  deallocate(RIPB1)
  deallocate(RIPO1R)
  deallocate(OMCER)
  deallocate(OMNER)
  deallocate(OMPER)
  deallocate(OMCI)

  end subroutine DestructMicrobialData

end module MicrobialDataType
