module MicrobialDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMConfig, only : jcplx => jcplxc,nlbiomcp=>NumLiveMicrbCompts,NumMicbFunGrupsPerCmplx=> NumMicbFunGrupsPerCmplx
  use GridConsts
  use ElmIDMod, only : NumPlantChemElms
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable :: mBiomeHeter_vr(:,:,:,:,:,:)    !microbial biomass element	[g d-2]
  real(r8),target,allocatable :: RO2DmndHetert(:,:,:,:,:)    !aqueous O2 demand	[g d-2 h-1]
  real(r8),target,allocatable :: RDOCUptkHeter_vr(:,:,:,:,:)    !net microbial DOC flux	[g d-2 h-1]
  real(r8),target,allocatable :: RAcetateUptkHeter_vr(:,:,:,:,:)    !net microbial acetate flux	[g d-2 h-1]
  real(r8),target,allocatable :: RNH4DmndSoilHeter_vr(:,:,:,:,:)    !microbial NH4 demand in soil	[g d-2 h-1]
  real(r8),target,allocatable :: RNO3DmndSoilHeter_vr(:,:,:,:,:)    !microbial NO3 demand in soil	[g d-2 h-1]
  real(r8),target,allocatable :: RH2PO4DmndSoilHeter_vr(:,:,:,:,:)    !microbial PO4 demand in soil	[g d-2 h-1]
  real(r8),target,allocatable :: RNH4DmndLitrHeter_col(:,:,:,:)     !microbial NH4 demand in surface litter	[g d-2 h-1]
  real(r8),target,allocatable :: RH2PO4DmndLitrHeter_col(:,:,:,:)     !microbial PO4 demand in surface litter	[g d-2 h-1]
  real(r8),target,allocatable :: RNO3DmndLitrHeter_col(:,:,:,:)     !microbial NO3 demand in surface litter	[g d-2 h-1]
  real(r8),target,allocatable :: RNO3ReduxDmndSoilHeter_vr(:,:,:,:,:)    !total microbial NO3 uptake non-band unconstrained by NO3	[g d-2 h-1]
  real(r8),target,allocatable :: RNO2DmndReduxSoilHeter_vr(:,:,:,:,:)    !total microbial NO2 uptake non-band unconstrained by NO2	[g d-2 h-1]
  real(r8),target,allocatable :: RNO3ReduxDmndBandHeter_vr(:,:,:,:,:)    !total microbial NO3 uptake band unconstrained by NO3	[g d-2 h-1]
  real(r8),target,allocatable :: RNO2DmndReduxBandHeter_vr(:,:,:,:,:)    !total microbial NO2 uptake band unconstrained by NH4	[g d-2 h-1]
  real(r8),target,allocatable :: RN2ODmndReduxHeter_vr(:,:,:,:,:)    !total microbial N2O uptake unconstrained by N2O	[g d-2 h-1]
  real(r8),target,allocatable :: RNH4DmndBandHeter_vr(:,:,:,:,:)    !microbial NH4 immobilization (+ve) - mineralization (-ve) band	[g d-2 h-1]
  real(r8),target,allocatable :: RNO3DmndBandHeter_vr(:,:,:,:,:)    !microbial NO3 immobilization (+ve) - mineralization (-ve) band	[g d-2 h-1]
  real(r8),target,allocatable :: RH2PO4DmndBandHeter_vr(:,:,:,:,:)    !substrate-unlimited H2PO4 mineraln-immobiln
  real(r8),target,allocatable :: RH1PO4DmndSoilHeter_vr(:,:,:,:,:)    !substrate-unlimited HPO4 immobilization
  real(r8),target,allocatable :: RH1PO4DmndBandHeter_vr(:,:,:,:,:)    !substrate-unlimited HPO4 mineraln-immobiln
  real(r8),target,allocatable :: RH1PO4DmndLitrHeter_col(:,:,:,:)     !substrate-unlimited HPO4 immobilization
  real(r8),target,allocatable :: OMEERhetr(:,:,:,:,:,:,:)  !microbial C  erosion 	[g d-2 h-1]

  real(r8),target,allocatable :: mBiomeAutor_vr(:,:,:,:,:)
  real(r8),target,allocatable :: RO2DmndAutort_vr(:,:,:,:)
  real(r8),target,allocatable :: RNH4UptkSoilAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RNO3UptkSoilAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RH2PO4UptkSoilAutor_vr(:,:,:,:)
  real(r8),target,allocatable :: RNH4UptkLitrAutor_col(:,:,:)
  real(r8),target,allocatable :: RH2PO4UptkLitrAutor_col(:,:,:)
  real(r8),target,allocatable :: RNO3UptkLitrAutor_col(:,:,:)
  real(r8),target,allocatable :: RNH3OxidAutor(:,:,:,:)
  real(r8),target,allocatable :: RNO2OxidAutor(:,:,:,:)
  real(r8),target,allocatable :: RNH3OxidAutorBand(:,:,:,:)
  real(r8),target,allocatable :: RNO2OxidAutorBand(:,:,:,:)
  real(r8),target,allocatable :: RN2ODmndReduxAutor_vr(:,:,:,:)
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
  allocate(mBiomeHeter_vr(NumPlantChemElms,NumLiveHeterBioms,1:jcplx,0:JZ,JY,JX));mBiomeHeter_vr=0._r8
  allocate(RO2DmndHetert(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX)); RO2DmndHetert=0._r8
  allocate(RDOCUptkHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RDOCUptkHeter_vr=0._r8
  allocate(RAcetateUptkHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RAcetateUptkHeter_vr=0._r8
  allocate(RNH4DmndSoilHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RNH4DmndSoilHeter_vr=0._r8
  allocate(RNO3DmndSoilHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RNO3DmndSoilHeter_vr=0._r8
  allocate(RH2PO4DmndSoilHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RH2PO4DmndSoilHeter_vr=0._r8
  allocate(RNH4DmndLitrHeter_col(NumHetetrMicCmplx,1:jcplx,JY,JX));RNH4DmndLitrHeter_col=0._r8
  allocate(RH2PO4DmndLitrHeter_col(NumHetetrMicCmplx,1:jcplx,JY,JX));RH2PO4DmndLitrHeter_col=0._r8
  allocate(RNO3DmndLitrHeter_col(NumHetetrMicCmplx,1:jcplx,JY,JX));RNO3DmndLitrHeter_col=0._r8
  allocate(RNO3ReduxDmndSoilHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RNO3ReduxDmndSoilHeter_vr=0._r8
  allocate(RNO2DmndReduxSoilHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RNO2DmndReduxSoilHeter_vr=0._r8
  allocate(RNO3ReduxDmndBandHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RNO3ReduxDmndBandHeter_vr=0._r8
  allocate(RNO2DmndReduxBandHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RNO2DmndReduxBandHeter_vr=0._r8
  allocate(RN2ODmndReduxHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RN2ODmndReduxHeter_vr=0._r8
  allocate(RNH4DmndBandHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RNH4DmndBandHeter_vr=0._r8
  allocate(RNO3DmndBandHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RNO3DmndBandHeter_vr=0._r8
  allocate(RH2PO4DmndBandHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RH2PO4DmndBandHeter_vr=0._r8
  allocate(RH1PO4DmndSoilHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RH1PO4DmndSoilHeter_vr=0._r8
  allocate(RH1PO4DmndBandHeter_vr(NumHetetrMicCmplx,1:jcplx,0:JZ,JY,JX));RH1PO4DmndBandHeter_vr=0._r8
  allocate(RH1PO4DmndLitrHeter_col(NumHetetrMicCmplx,1:jcplx,JY,JX));RH1PO4DmndLitrHeter_col=0._r8
  allocate(OMEERhetr(NumPlantChemElms,NumLiveHeterBioms,1:jcplx,2,2,JV,JH));OMEERhetr=0._r8

  allocate(mBiomeAutor_vr(NumPlantChemElms,NumLiveAutoBioms,0:JZ,JY,JX));mBiomeAutor_vr=0._r8
  allocate(RO2DmndAutort_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RO2DmndAutort_vr=0._r8
  allocate(RNH4UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNH4UptkSoilAutor_vr=0._r8
  allocate(RNO3UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNO3UptkSoilAutor_vr=0._r8
  allocate(RH2PO4UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RH2PO4UptkSoilAutor_vr=0._r8
  allocate(RNH4UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX));RNH4UptkLitrAutor_col=0._r8
  allocate(RH2PO4UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX));RH2PO4UptkLitrAutor_col=0._r8
  allocate(RNO3UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX));RNO3UptkLitrAutor_col=0._r8
  allocate(RNH3OxidAutor(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNH3OxidAutor=0._r8
  allocate(RNO2OxidAutor(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNO2OxidAutor=0._r8
  allocate(RNH3OxidAutorBand(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNH3OxidAutorBand=0._r8
  allocate(RNO2OxidAutorBand(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNO2OxidAutorBand=0._r8
  allocate(RN2ODmndReduxAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RN2ODmndReduxAutor_vr=0._r8
  allocate(RNH4UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNH4UptkBandAutor_vr=0._r8
  allocate(RNO3UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNO3UptkBandAutor_vr=0._r8
  allocate(RH2PO4UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RH2PO4UptkBandAutor_vr=0._r8
  allocate(RH1PO4UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RH1PO4UptkSoilAutor_vr=0._r8
  allocate(RH1PO4UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RH1PO4UptkBandAutor_vr=0._r8
  allocate(RH1PO4UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX));RH1PO4UptkLitrAutor_col=0._r8
  allocate(OMEERauto(NumPlantChemElms,NumLiveAutoBioms,2,2,JV,JH));OMEERauto=0._r8
  end subroutine InitAllocate
!----------------------------------------------------------------------------------------------

  subroutine DestructMicrobialData
  implicit none

  if(allocated(mBiomeHeter_vr))deallocate(mBiomeHeter_vr)
  if(allocated(RO2DmndHetert))deallocate(RO2DmndHetert)
  if(allocated(RDOCUptkHeter_vr))deallocate(RDOCUptkHeter_vr)
  if(allocated(RAcetateUptkHeter_vr))deallocate(RAcetateUptkHeter_vr)
  if(allocated(RNH4DmndSoilHeter_vr))deallocate(RNH4DmndSoilHeter_vr)
  if(allocated(RNO3DmndSoilHeter_vr))deallocate(RNO3DmndSoilHeter_vr)
  if(allocated(RH2PO4DmndSoilHeter_vr))deallocate(RH2PO4DmndSoilHeter_vr)
  if(allocated(RNH4DmndLitrHeter_col))deallocate(RNH4DmndLitrHeter_col)
  if(allocated(RH2PO4DmndLitrHeter_col))deallocate(RH2PO4DmndLitrHeter_col)
  if(allocated(RNO3DmndLitrHeter_col))deallocate(RNO3DmndLitrHeter_col)
  if(allocated(RNO3ReduxDmndSoilHeter_vr))deallocate(RNO3ReduxDmndSoilHeter_vr)
  if(allocated(RNO2DmndReduxSoilHeter_vr))deallocate(RNO2DmndReduxSoilHeter_vr)
  if(allocated(RNO3ReduxDmndBandHeter_vr))deallocate(RNO3ReduxDmndBandHeter_vr)
  if(allocated(RNO2DmndReduxBandHeter_vr))deallocate(RNO2DmndReduxBandHeter_vr)
  if(allocated(RN2ODmndReduxHeter_vr))deallocate(RN2ODmndReduxHeter_vr)
  if(allocated(RNH4DmndBandHeter_vr))deallocate(RNH4DmndBandHeter_vr)
  if(allocated(RNO3DmndBandHeter_vr))deallocate(RNO3DmndBandHeter_vr)
  if(allocated(RH2PO4DmndBandHeter_vr))deallocate(RH2PO4DmndBandHeter_vr)
  if(allocated(RH1PO4DmndSoilHeter_vr))deallocate(RH1PO4DmndSoilHeter_vr)
  if(allocated(RH1PO4DmndBandHeter_vr))deallocate(RH1PO4DmndBandHeter_vr)
  if(allocated(RH1PO4DmndLitrHeter_col))deallocate(RH1PO4DmndLitrHeter_col)
  if(allocated(OMEERhetr))deallocate(OMEERhetr)

  if(allocated(mBiomeAutor_vr))deallocate(mBiomeAutor_vr)
  if(allocated(RO2DmndAutort_vr))deallocate(RO2DmndAutort_vr)
  if(allocated(RNH4UptkSoilAutor_vr))deallocate(RNH4UptkSoilAutor_vr)
  if(allocated(RNO3UptkSoilAutor_vr))deallocate(RNO3UptkSoilAutor_vr)
  if(allocated(RH2PO4UptkSoilAutor_vr))deallocate(RH2PO4UptkSoilAutor_vr)
  if(allocated(RNH4UptkLitrAutor_col))deallocate(RNH4UptkLitrAutor_col)
  if(allocated(RH2PO4UptkLitrAutor_col))deallocate(RH2PO4UptkLitrAutor_col)
  if(allocated(RNO3UptkLitrAutor_col))deallocate(RNO3UptkLitrAutor_col)
  if(allocated(RNH3OxidAutor))deallocate(RNH3OxidAutor)
  if(allocated(RNO2OxidAutor))deallocate(RNO2OxidAutor)
  if(allocated(RNH3OxidAutorBand))deallocate(RNH3OxidAutorBand)
  if(allocated(RNO2OxidAutorBand))deallocate(RNO2OxidAutorBand)
  if(allocated(RN2ODmndReduxAutor_vr))deallocate(RN2ODmndReduxAutor_vr)
  if(allocated(RNH4UptkBandAutor_vr))deallocate(RNH4UptkBandAutor_vr)
  if(allocated(RNO3UptkBandAutor_vr))deallocate(RNO3UptkBandAutor_vr)
  if(allocated(RH2PO4UptkBandAutor_vr))deallocate(RH2PO4UptkBandAutor_vr)
  if(allocated(RH1PO4UptkSoilAutor_vr))deallocate(RH1PO4UptkSoilAutor_vr)
  if(allocated(RH1PO4UptkBandAutor_vr))deallocate(RH1PO4UptkBandAutor_vr)
  if(allocated(RH1PO4UptkLitrAutor_col))deallocate(RH1PO4UptkLitrAutor_col)
  if(allocated(OMEERauto))deallocate(OMEERauto)
  end subroutine DestructMicrobialData

end module MicrobialDataType
