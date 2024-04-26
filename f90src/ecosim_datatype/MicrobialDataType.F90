module MicrobialDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMConfig, only : jcplx => jcplxc,nlbiomcp=>NumLiveMicrbCompts,NumMicbFunGroups=> NumMicbFunGroups
  use GridConsts
  use ElmIDMod, only : NumPlantChemElms
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable :: OMEheter(:,:,:,:,:,:)    !microbial biomass element	[g d-2]
  real(r8),target,allocatable :: RO2DmndHetert(:,:,:,:,:)    !aqueous O2 demand	[g d-2 h-1]
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
  real(r8),target,allocatable :: RO2DmndAutort(:,:,:,:)
  real(r8),target,allocatable :: RNH4UptkSoilAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RNO3UptkSoilAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RH2PO4UptkSoilAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RNH4UptkLitrAutor_col(:,:,:)
  real(r8),target,allocatable :: RH2PO4UptkLitrAutor_col(:,:,:)
  real(r8),target,allocatable :: RNO3UptkLitrAutor_col(:,:,:)
  real(r8),target,allocatable :: RNH3OxidAutor(:,:,:,:)
  real(r8),target,allocatable :: RVMX3ff(:,:,:,:)
  real(r8),target,allocatable :: RNO2OxidAutor(:,:,:,:)
  real(r8),target,allocatable :: RNH3OxidAutorBand(:,:,:,:)
  real(r8),target,allocatable :: RVMB3ff(:,:,:,:)
  real(r8),target,allocatable :: RNO2OxidAutorBand(:,:,:,:)
  real(r8),target,allocatable :: RVMX1ff(:,:,:,:)
  real(r8),target,allocatable :: RNH4UptkBandAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RNO3UptkBandAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RH2PO4UptkBandAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RH1PO4UptkSoilAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RH1PO4UptkBandAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RH1PO4UptkLitrAutor_col(:,:,:)
  real(r8),target,allocatable :: OMEERauto(:,:,:,:,:,:)

  private :: InitAllocate

  contains

  subroutine InitMicrobialData()

  implicit none

  call InitAllocate

  end subroutine InitMicrobialData

  subroutine InitAllocate

  implicit none
  allocate(OMEheter(NumPlantChemElms,NumLiveHeterBioms,1:jcplx,0:JZ,JY,JX))
  allocate(RO2DmndHetert(NumMicrbHetetrophCmplx,1:jcplx,0:JZ,JY,JX))
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
  allocate(OMEERhetr(NumPlantChemElms,NumLiveHeterBioms,1:jcplx,2,2,JV,JH))

  allocate(OMEauto(NumPlantChemElms,NumLiveAutoBioms,0:JZ,JY,JX))
  allocate(RO2DmndAutort(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RNH4UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RNO3UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RH2PO4UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RNH4UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX))
  allocate(RH2PO4UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX))
  allocate(RNO3UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX))
  allocate(RNH3OxidAutor(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RVMX3ff(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RNO2OxidAutor(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RNH3OxidAutorBand(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RVMB3ff(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RNO2OxidAutorBand(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RVMX1ff(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RNH4UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RNO3UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RH2PO4UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RH1PO4UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RH1PO4UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX))
  allocate(RH1PO4UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX))
  allocate(OMEERauto(NumPlantChemElms,NumLiveAutoBioms,2,2,JV,JH))
  end subroutine InitAllocate
!----------------------------------------------------------------------------------------------

  subroutine DestructMicrobialData
  implicit none

  if(allocated(OMEheter))deallocate(OMEheter)
  if(allocated(RO2DmndHetert))deallocate(RO2DmndHetert)
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
  if(allocated(RO2DmndAutort))deallocate(RO2DmndAutort)
  if(allocated(RNH4UptkSoilAutor_vr))deallocate(RNH4UptkSoilAutor_vr)
  if(allocated(RNO3UptkSoilAutor_vr))deallocate(RNO3UptkSoilAutor_vr)
  if(allocated(RH2PO4UptkSoilAutor_vr))deallocate(RH2PO4UptkSoilAutor_vr)
  if(allocated(RNH4UptkLitrAutor_col))deallocate(RNH4UptkLitrAutor_col)
  if(allocated(RH2PO4UptkLitrAutor_col))deallocate(RH2PO4UptkLitrAutor_col)
  if(allocated(RNO3UptkLitrAutor_col))deallocate(RNO3UptkLitrAutor_col)
  if(allocated(RNH3OxidAutor))deallocate(RNH3OxidAutor)
  if(allocated(RVMX3ff))deallocate(RVMX3ff)
  if(allocated(RNO2OxidAutor))deallocate(RNO2OxidAutor)
  if(allocated(RNH3OxidAutorBand))deallocate(RNH3OxidAutorBand)
  if(allocated(RVMB3ff))deallocate(RVMB3ff)
  if(allocated(RNO2OxidAutorBand))deallocate(RNO2OxidAutorBand)
  if(allocated(RVMX1ff))deallocate(RVMX1ff)
  if(allocated(RNH4UptkBandAutor_vr))deallocate(RNH4UptkBandAutor_vr)
  if(allocated(RNO3UptkBandAutor_vr))deallocate(RNO3UptkBandAutor_vr)
  if(allocated(RH2PO4UptkBandAutor_vr))deallocate(RH2PO4UptkBandAutor_vr)
  if(allocated(RH1PO4UptkSoilAutor_vr))deallocate(RH1PO4UptkSoilAutor_vr)
  if(allocated(RH1PO4UptkBandAutor_vr))deallocate(RH1PO4UptkBandAutor_vr)
  if(allocated(RH1PO4UptkLitrAutor_col))deallocate(RH1PO4UptkLitrAutor_col)
  if(allocated(OMEERauto))deallocate(OMEERauto)
  end subroutine DestructMicrobialData

end module MicrobialDataType
