module MicrobialDataType
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use EcoSIMConfig, only : jcplx => jcplxc,nlbiomcp=>nlbiomcpc,NFGs=> NFGsc
  use GridConsts
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8),allocatable :: OMC(:,:,:,:,:,:)    !microbial biomass C	[g d-2]
  real(r8),allocatable :: OMN(:,:,:,:,:,:)    !microbial biomass N	[g d-2]
  real(r8),allocatable :: OMP(:,:,:,:,:,:)    !microbial biomass P	[g d-2]
  real(r8),allocatable :: ROXYS(:,:,:,:,:)    !aqueous O2 demand	[g d-2 h-1]
  real(r8),allocatable :: ROQCS(:,:,:,:,:)    !net microbial DOC flux	[g d-2 h-1]
  real(r8),allocatable :: ROQAS(:,:,:,:,:)    !net microbial acetate flux	[g d-2 h-1]
  real(r8),allocatable :: RINHO(:,:,:,:,:)    !microbial NH4 demand in soil	[g d-2 h-1]
  real(r8),allocatable :: RINOO(:,:,:,:,:)    !microbial NO3 demand in soil	[g d-2 h-1]
  real(r8),allocatable :: RIPOO(:,:,:,:,:)    !microbial PO4 demand in soil	[g d-2 h-1]
  real(r8),allocatable :: RINHOR(:,:,:,:)     !microbial NH4 demand in surface litter	[g d-2 h-1]
  real(r8),allocatable :: RIPOOR(:,:,:,:)     !microbial PO4 demand in surface litter	[g d-2 h-1]
  real(r8),allocatable :: RINOOR(:,:,:,:)     !microbial NO3 demand in surface litter	[g d-2 h-1]
  real(r8),allocatable :: RVMX4(:,:,:,:,:)    !total microbial NH4 uptake non-band unconstrained by NH4	[g d-2 h-1]
  real(r8),allocatable :: RVMX3(:,:,:,:,:)    !total microbial NO3 uptake non-band unconstrained by NO3	[g d-2 h-1]
  real(r8),allocatable :: RVMX2(:,:,:,:,:)    !total microbial NO2 uptake non-band unconstrained by NO2	[g d-2 h-1]
  real(r8),allocatable :: RVMB4(:,:,:,:,:)    !total microbial NH4 uptake non-band unconstrained by NH4	[g d-2 h-1]
  real(r8),allocatable :: RVMB3(:,:,:,:,:)    !total microbial NO3 uptake band unconstrained by NO3	[g d-2 h-1]
  real(r8),allocatable :: RVMB2(:,:,:,:,:)    !total microbial NO2 uptake band unconstrained by NH4	[g d-2 h-1]
  real(r8),allocatable :: RVMX1(:,:,:,:,:)    !total microbial N2O uptake unconstrained by N2O	[g d-2 h-1]
  real(r8),allocatable :: RINHB(:,:,:,:,:)    !microbial NH4 immobilization (+ve) - mineralization (-ve) band	[g d-2 h-1]
  real(r8),allocatable :: RINOB(:,:,:,:,:)    !microbial NO3 immobilization (+ve) - mineralization (-ve) band	[g d-2 h-1]
  real(r8),allocatable :: RIPBO(:,:,:,:,:)    !substrate-unlimited H2PO4 mineraln-immobiln
  real(r8),allocatable :: RIPO1(:,:,:,:,:)    !substrate-unlimited HPO4 immobilization
  real(r8),allocatable :: RIPB1(:,:,:,:,:)    !substrate-unlimited HPO4 mineraln-immobiln
  real(r8),allocatable :: RIPO1R(:,:,:,:)     !substrate-unlimited HPO4 immobilization
  real(r8),allocatable :: OMCER(:,:,:,:,:,:)  !microbial C  erosion 	[g d-2 h-1]
  real(r8),allocatable :: OMNER(:,:,:,:,:,:)  !microbial N  erosion 	[g d-2 h-1]
  real(r8),allocatable :: OMPER(:,:,:,:,:,:)  !microbial P  erosion 	[g d-2 h-1]


  real(r8),allocatable :: OMCff(:,:,:,:,:)
  real(r8),allocatable :: OMNff(:,:,:,:,:)
  real(r8),allocatable :: OMPff(:,:,:,:,:)
  real(r8),allocatable :: ROXYSff(:,:,:,:)
  real(r8),allocatable :: RINHOff(:,:,:,:)
  real(r8),allocatable :: RINOOff(:,:,:,:)
  real(r8),allocatable :: RIPOOff(:,:,:,:)
  real(r8),allocatable :: RINHORff(:,:,:)
  real(r8),allocatable :: RIPOORff(:,:,:)
  real(r8),allocatable :: RINOORff(:,:,:)
  real(r8),allocatable :: RVMX4ff(:,:,:,:)
  real(r8),allocatable :: RVMX3ff(:,:,:,:)
  real(r8),allocatable :: RVMX2ff(:,:,:,:)
  real(r8),allocatable :: RVMB4ff(:,:,:,:)
  real(r8),allocatable :: RVMB3ff(:,:,:,:)
  real(r8),allocatable :: RVMB2ff(:,:,:,:)
  real(r8),allocatable :: RVMX1ff(:,:,:,:)
  real(r8),allocatable :: RINHBff(:,:,:,:)
  real(r8),allocatable :: RINOBff(:,:,:,:)
  real(r8),allocatable :: RIPBOff(:,:,:,:)
  real(r8),allocatable :: RIPO1ff(:,:,:,:)
  real(r8),allocatable :: RIPB1ff(:,:,:,:)
  real(r8),allocatable :: RIPO1Rff(:,:,:)
  real(r8),allocatable :: OMCERff(:,:,:,:,:)
  real(r8),allocatable :: OMNERff(:,:,:,:,:)
  real(r8),allocatable :: OMPERff(:,:,:,:,:)

  private :: InitAllocate

  contains

  subroutine InitMicrobialData()

  implicit none

  call InitAllocate

  end subroutine InitMicrobialData

  subroutine InitAllocate

  implicit none
  allocate(OMC(nlbiomcp,NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(OMN(nlbiomcp,NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(OMP(nlbiomcp,NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(ROXYS(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(ROQCS(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(ROQAS(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RINHO(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RINOO(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RIPOO(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RINHOR(NMICBSO,1:jcplx,JY,JX))
  allocate(RIPOOR(NMICBSO,1:jcplx,JY,JX))
  allocate(RINOOR(NMICBSO,1:jcplx,JY,JX))
  allocate(RVMX4(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RVMX3(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RVMX2(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RVMB4(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RVMB3(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RVMB2(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RVMX1(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RINHB(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RINOB(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RIPBO(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RIPO1(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RIPB1(NMICBSO,1:jcplx,0:JZ,JY,JX))
  allocate(RIPO1R(NMICBSO,1:jcplx,JY,JX))
  allocate(OMCER(nlbiomcp*NMICBSO,1:jcplx,2,2,JV,JH))
  allocate(OMNER(nlbiomcp*NMICBSO,1:jcplx,2,2,JV,JH))
  allocate(OMPER(nlbiomcp*NMICBSO,1:jcplx,2,2,JV,JH))

  allocate(OMCff(nlbiomcp,NMICBSA,0:JZ,JY,JX))
  allocate(OMNff(nlbiomcp,NMICBSA,0:JZ,JY,JX))
  allocate(OMPff(nlbiomcp,NMICBSA,0:JZ,JY,JX))
  allocate(ROXYSff(NMICBSA,0:JZ,JY,JX))
  allocate(RINHOff(NMICBSA,0:JZ,JY,JX))
  allocate(RINOOff(NMICBSA,0:JZ,JY,JX))
  allocate(RIPOOff(NMICBSA,0:JZ,JY,JX))
  allocate(RINHORff(NMICBSA,JY,JX))
  allocate(RIPOORff(NMICBSA,JY,JX))
  allocate(RINOORff(NMICBSA,JY,JX))
  allocate(RVMX4ff(NMICBSA,0:JZ,JY,JX))
  allocate(RVMX3ff(NMICBSA,0:JZ,JY,JX))
  allocate(RVMX2ff(NMICBSA,0:JZ,JY,JX))
  allocate(RVMB4ff(NMICBSA,0:JZ,JY,JX))
  allocate(RVMB3ff(NMICBSA,0:JZ,JY,JX))
  allocate(RVMB2ff(NMICBSA,0:JZ,JY,JX))
  allocate(RVMX1ff(NMICBSA,0:JZ,JY,JX))
  allocate(RINHBff(NMICBSA,0:JZ,JY,JX))
  allocate(RINOBff(NMICBSA,0:JZ,JY,JX))
  allocate(RIPBOff(NMICBSA,0:JZ,JY,JX))
  allocate(RIPO1ff(NMICBSA,0:JZ,JY,JX))
  allocate(RIPB1ff(NMICBSA,0:JZ,JY,JX))
  allocate(RIPO1Rff(NMICBSA,JY,JX))
  allocate(OMCERff(nlbiomcp*NMICBSA,2,2,JV,JH))
  allocate(OMNERff(nlbiomcp*NMICBSA,2,2,JV,JH))
  allocate(OMPERff(nlbiomcp*NMICBSA,2,2,JV,JH))
  end subroutine InitAllocate
!----------------------------------------------------------------------------------------------

  subroutine DestructMicrobialData
  implicit none

  if(allocated(OMC))deallocate(OMC)
  if(allocated(OMN))deallocate(OMN)
  if(allocated(OMP))deallocate(OMP)
  if(allocated(ROXYS))deallocate(ROXYS)
  if(allocated(ROQCS))deallocate(ROQCS)
  if(allocated(ROQAS))deallocate(ROQAS)
  if(allocated(RINHO))deallocate(RINHO)
  if(allocated(RINOO))deallocate(RINOO)
  if(allocated(RIPOO))deallocate(RIPOO)
  if(allocated(RINHOR))deallocate(RINHOR)
  if(allocated(RIPOOR))deallocate(RIPOOR)
  if(allocated(RINOOR))deallocate(RINOOR)
  if(allocated(RVMX4))deallocate(RVMX4)
  if(allocated(RVMX3))deallocate(RVMX3)
  if(allocated(RVMX2))deallocate(RVMX2)
  if(allocated(RVMB4))deallocate(RVMB4)
  if(allocated(RVMB3))deallocate(RVMB3)
  if(allocated(RVMB2))deallocate(RVMB2)
  if(allocated(RVMX1))deallocate(RVMX1)
  if(allocated(RINHB))deallocate(RINHB)
  if(allocated(RINOB))deallocate(RINOB)
  if(allocated(RIPBO))deallocate(RIPBO)
  if(allocated(RIPO1))deallocate(RIPO1)
  if(allocated(RIPB1))deallocate(RIPB1)
  if(allocated(RIPO1R))deallocate(RIPO1R)
  if(allocated(OMCER))deallocate(OMCER)
  if(allocated(OMNER))deallocate(OMNER)
  if(allocated(OMPER))deallocate(OMPER)

  if(allocated(OMCff))deallocate(OMCff)
  if(allocated(OMNff))deallocate(OMNff)
  if(allocated(OMPff))deallocate(OMPff)
  if(allocated(ROXYSff))deallocate(ROXYSff)
  if(allocated(RINHOff))deallocate(RINHOff)
  if(allocated(RINOOff))deallocate(RINOOff)
  if(allocated(RIPOOff))deallocate(RIPOOff)
  if(allocated(RINHORff))deallocate(RINHORff)
  if(allocated(RIPOORff))deallocate(RIPOORff)
  if(allocated(RINOORff))deallocate(RINOORff)
  if(allocated(RVMX4ff))deallocate(RVMX4ff)
  if(allocated(RVMX3ff))deallocate(RVMX3ff)
  if(allocated(RVMX2ff))deallocate(RVMX2ff)
  if(allocated(RVMB4ff))deallocate(RVMB4ff)
  if(allocated(RVMB3ff))deallocate(RVMB3ff)
  if(allocated(RVMB2ff))deallocate(RVMB2ff)
  if(allocated(RVMX1ff))deallocate(RVMX1ff)
  if(allocated(RINHBff))deallocate(RINHBff)
  if(allocated(RINOBff))deallocate(RINOBff)
  if(allocated(RIPBOff))deallocate(RIPBOff)
  if(allocated(RIPO1ff))deallocate(RIPO1ff)
  if(allocated(RIPB1ff))deallocate(RIPB1ff)
  if(allocated(RIPO1Rff))deallocate(RIPO1Rff)
  if(allocated(OMCERff))deallocate(OMCERff)
  if(allocated(OMNERff))deallocate(OMNERff)
  if(allocated(OMPERff))deallocate(OMPERff)
  end subroutine DestructMicrobialData

end module MicrobialDataType
