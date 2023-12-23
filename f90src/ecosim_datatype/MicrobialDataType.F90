module MicrobialDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMConfig, only : jcplx => jcplxc,nlbiomcp=>NumLiveMicrbCompts,NumMicbFunGroups=> NumMicbFunGroups
  use GridConsts
  use ElmIDMod, only : NumPlantChemElmnts
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable :: OMEhetr(:,:,:,:,:,:)    !microbial biomass element	[g d-2]
  real(r8),target,allocatable :: ROXYS(:,:,:,:,:)    !aqueous O2 demand	[g d-2 h-1]
  real(r8),target,allocatable :: ROQCS(:,:,:,:,:)    !net microbial DOC flux	[g d-2 h-1]
  real(r8),target,allocatable :: ROQAS(:,:,:,:,:)    !net microbial acetate flux	[g d-2 h-1]
  real(r8),target,allocatable :: RINHO(:,:,:,:,:)    !microbial NH4 demand in soil	[g d-2 h-1]
  real(r8),target,allocatable :: RINOO(:,:,:,:,:)    !microbial NO3 demand in soil	[g d-2 h-1]
  real(r8),target,allocatable :: RIPOO(:,:,:,:,:)    !microbial PO4 demand in soil	[g d-2 h-1]
  real(r8),target,allocatable :: RINHOR(:,:,:,:)     !microbial NH4 demand in surface litter	[g d-2 h-1]
  real(r8),target,allocatable :: RIPOOR(:,:,:,:)     !microbial PO4 demand in surface litter	[g d-2 h-1]
  real(r8),target,allocatable :: RINOOR(:,:,:,:)     !microbial NO3 demand in surface litter	[g d-2 h-1]
  real(r8),target,allocatable :: RVMX4(:,:,:,:,:)    !total microbial NH4 uptake non-band unconstrained by NH4	[g d-2 h-1]
  real(r8),target,allocatable :: RVMX3(:,:,:,:,:)    !total microbial NO3 uptake non-band unconstrained by NO3	[g d-2 h-1]
  real(r8),target,allocatable :: RVMX2(:,:,:,:,:)    !total microbial NO2 uptake non-band unconstrained by NO2	[g d-2 h-1]
  real(r8),target,allocatable :: RVMB4(:,:,:,:,:)    !total microbial NH4 uptake non-band unconstrained by NH4	[g d-2 h-1]
  real(r8),target,allocatable :: RVMB3(:,:,:,:,:)    !total microbial NO3 uptake band unconstrained by NO3	[g d-2 h-1]
  real(r8),target,allocatable :: RVMB2(:,:,:,:,:)    !total microbial NO2 uptake band unconstrained by NH4	[g d-2 h-1]
  real(r8),target,allocatable :: RVMX1(:,:,:,:,:)    !total microbial N2O uptake unconstrained by N2O	[g d-2 h-1]
  real(r8),target,allocatable :: RINHB(:,:,:,:,:)    !microbial NH4 immobilization (+ve) - mineralization (-ve) band	[g d-2 h-1]
  real(r8),target,allocatable :: RINOB(:,:,:,:,:)    !microbial NO3 immobilization (+ve) - mineralization (-ve) band	[g d-2 h-1]
  real(r8),target,allocatable :: RIPBO(:,:,:,:,:)    !substrate-unlimited H2PO4 mineraln-immobiln
  real(r8),target,allocatable :: RIPO1(:,:,:,:,:)    !substrate-unlimited HPO4 immobilization
  real(r8),target,allocatable :: RIPB1(:,:,:,:,:)    !substrate-unlimited HPO4 mineraln-immobiln
  real(r8),target,allocatable :: RIPO1R(:,:,:,:)     !substrate-unlimited HPO4 immobilization
  real(r8),target,allocatable :: OMEERhetr(:,:,:,:,:,:,:)  !microbial C  erosion 	[g d-2 h-1]

  real(r8),target,allocatable :: OMEauto(:,:,:,:,:)
  real(r8),target,allocatable :: ROXYSff(:,:,:,:)
  real(r8),target,allocatable :: RINHOff(:,:,:,:)
  real(r8),target,allocatable :: RINOOff(:,:,:,:)
  real(r8),target,allocatable :: RIPOOff(:,:,:,:)
  real(r8),target,allocatable :: RINHORff(:,:,:)
  real(r8),target,allocatable :: RIPOORff(:,:,:)
  real(r8),target,allocatable :: RINOORff(:,:,:)
  real(r8),target,allocatable :: RVMX4ff(:,:,:,:)
  real(r8),target,allocatable :: RVMX3ff(:,:,:,:)
  real(r8),target,allocatable :: RVMX2ff(:,:,:,:)
  real(r8),target,allocatable :: RVMB4ff(:,:,:,:)
  real(r8),target,allocatable :: RVMB3ff(:,:,:,:)
  real(r8),target,allocatable :: RVMB2ff(:,:,:,:)
  real(r8),target,allocatable :: RVMX1ff(:,:,:,:)
  real(r8),target,allocatable :: RINHBff(:,:,:,:)
  real(r8),target,allocatable :: RINOBff(:,:,:,:)
  real(r8),target,allocatable :: RIPBOff(:,:,:,:)
  real(r8),target,allocatable :: RIPO1ff(:,:,:,:)
  real(r8),target,allocatable :: RIPB1ff(:,:,:,:)
  real(r8),target,allocatable :: RIPO1Rff(:,:,:)
  real(r8),target,allocatable :: OMEERauto(:,:,:,:,:,:)

  private :: InitAllocate

  contains

  subroutine InitMicrobialData()

  implicit none

  call InitAllocate

  end subroutine InitMicrobialData

  subroutine InitAllocate

  implicit none
  allocate(OMEhetr(NumPlantChemElmnts,NumLiveHeterBioms,1:jcplx,0:JZ,JY,JX))
  allocate(ROXYS(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(ROQCS(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(ROQAS(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RINHO(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RINOO(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RIPOO(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RINHOR(NumMicrbHetetrophCmplx,1:jcplx,JY,JX))
  allocate(RIPOOR(NumMicrbHetetrophCmplx,1:jcplx,JY,JX))
  allocate(RINOOR(NumMicrbHetetrophCmplx,1:jcplx,JY,JX))
  allocate(RVMX4(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RVMX3(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RVMX2(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RVMB4(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RVMB3(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RVMB2(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RVMX1(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RINHB(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RINOB(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RIPBO(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RIPO1(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RIPB1(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
  allocate(RIPO1R(NumMicrbHetetrophCmplx,1:jcplx,JY,JX))
  allocate(OMEERhetr(NumPlantChemElmnts,NumLiveHeterBioms,1:jcplx,2,2,JV,JH))

  allocate(OMEauto(NumPlantChemElmnts,NumLiveAutoBioms,0:JZ,JY,JX))
  allocate(ROXYSff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RINHOff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RINOOff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RIPOOff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RINHORff(NumMicrobAutotrophCmplx,JY,JX))
  allocate(RIPOORff(NumMicrobAutotrophCmplx,JY,JX))
  allocate(RINOORff(NumMicrobAutotrophCmplx,JY,JX))
  allocate(RVMX4ff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RVMX3ff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RVMX2ff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RVMB4ff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RVMB3ff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RVMB2ff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RVMX1ff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RINHBff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RINOBff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RIPBOff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RIPO1ff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RIPB1ff(NumMicrobAutotrophCmplx,0:JZ,JY,JX))
  allocate(RIPO1Rff(NumMicrobAutotrophCmplx,JY,JX))
  allocate(OMEERauto(NumPlantChemElmnts,NumLiveAutoBioms,2,2,JV,JH))
  end subroutine InitAllocate
!----------------------------------------------------------------------------------------------

  subroutine DestructMicrobialData
  implicit none

  if(allocated(OMEhetr))deallocate(OMEhetr)
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
  if(allocated(OMEERhetr))deallocate(OMEERhetr)

  if(allocated(OMEauto))deallocate(OMEauto)
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
  if(allocated(OMEERauto))deallocate(OMEERauto)
  end subroutine DestructMicrobialData

end module MicrobialDataType
