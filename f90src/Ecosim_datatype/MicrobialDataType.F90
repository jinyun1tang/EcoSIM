module MicrobialDataType
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use EcoSIMConfig,  only: jcplx => jcplxc, nlbiomcp=>NumLiveMicrbCompts, NumMicbFunGrupsPerCmplx=> NumMicbFunGrupsPerCmplx
  use ElmIDMod,      only: NumPlantChemElms
  use abortutils,    only: destroy
  use GridConsts
implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable :: mBiomeHeter_vr(:,:,:,:,:,:)        !microbial biomass chemical element,	[g d-2]
  real(r8),target,allocatable :: RO2DmndHetert(:,:,:,:,:)           !aqueous O2 demand,	[g d-2 h-1]
  real(r8),target,allocatable :: RDOCUptkHeter_vr(:,:,:,:,:)        !net microbial DOC flux,	[g d-2 h-1]
  real(r8),target,allocatable :: RAcetateUptkHeter_vr(:,:,:,:,:)    !net microbial acetate flux,	[g d-2 h-1]
  real(r8),target,allocatable :: RNH4DmndSoilHeter_vr(:,:,:,:,:)    !heterotrophic microbial NH4 demand in soil,	[g d-2 h-1]
  real(r8),target,allocatable :: RNO3DmndSoilHeter_vr(:,:,:,:,:)    !heterotrophic microbial NO3 demand in soil,	[g d-2 h-1]
  real(r8),target,allocatable :: RH2PO4DmndSoilHeter_vr(:,:,:,:,:)  !heterotrophic microbial PO4 demand in soil,	[g d-2 h-1]
  real(r8),target,allocatable :: RNH4DmndLitrHeter_col(:,:,:,:)     !heterotrophic microbial NH4 demand in surface litter,	[g d-2 h-1]
  real(r8),target,allocatable :: RH2PO4DmndLitrHeter_col(:,:,:,:)     !heterotrophic microbial PO4 demand in surface litter,	[g d-2 h-1]
  real(r8),target,allocatable :: RNO3DmndLitrHeter_col(:,:,:,:)     !heterotrophic microbial NO3 demand in surface litter,	[g d-2 h-1]
  real(r8),target,allocatable :: RNO3ReduxDmndSoilHeter_vr(:,:,:,:,:)    !total heterotrophic microbial NO3 uptake non-band unconstrained by NO3,	[g d-2 h-1]
  real(r8),target,allocatable :: RNO2DmndReduxSoilHeter_vr(:,:,:,:,:)    !total heterotrophic microbial NO2 uptake non-band unconstrained by NO2,	[g d-2 h-1]
  real(r8),target,allocatable :: RNO3ReduxDmndBandHeter_vr(:,:,:,:,:)    !total heterotrophic microbial NO3 uptake in band soil unconstrained by NO3,	[g d-2 h-1]
  real(r8),target,allocatable :: RNO2DmndReduxBandHeter_vr(:,:,:,:,:)    !total heterotrophic microbial NO2 uptake in band soil unconstrained by NO2,	[g d-2 h-1]
  real(r8),target,allocatable :: RN2ODmndReduxHeter_vr(:,:,:,:,:)    !total heterotrophic microbial N2O uptake unconstrained by N2O,	[g d-2 h-1]
  real(r8),target,allocatable :: RNH4DmndBandHeter_vr(:,:,:,:,:)    !heterotrophic microbial NH4 immobilization (+ve) - mineralization (-ve) in band soil,	[g d-2 h-1]
  real(r8),target,allocatable :: RNO3DmndBandHeter_vr(:,:,:,:,:)    !heterotrophic microbial NO3 immobilization (+ve) - mineralization (-ve) in band soil,	[g d-2 h-1]
  real(r8),target,allocatable :: RH2PO4DmndBandHeter_vr(:,:,:,:,:)    !heterotrophic substrate-unlimited H2PO4 mineraln-immobiln in band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RH1PO4DmndSoilHeter_vr(:,:,:,:,:)    !heterotrophic substrate-unlimited HPO4 immobilization in non-band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RH1PO4DmndBandHeter_vr(:,:,:,:,:)    !heterotrophic substrate-unlimited HPO4 mineraln-immobiln in band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RH1PO4DmndLitrHeter_col(:,:,:,:)     !heterotrophic substrate-unlimited HPO4 immobilization in surface litter, [g d-2 h-1]
  real(r8),target,allocatable :: OMEERhetr_2D(:,:,:,:,:,:,:)           !heterotrophic microbial C loss through erosion 	[g d-2 h-1]
  real(r8),target,allocatable :: mBiomeAutor_vr(:,:,:,:,:)            !autotrophic microbial biomass chemical element,[g d-2]
  real(r8),target,allocatable :: RO2DmndAutort_vr(:,:,:,:)            !aqueous O2 demand by autotrophic microbes, [g d-2 h-1]
  real(r8),target,allocatable :: RNH4UptkSoilAutor_vr(:,:,:,:)        !autotrophic microbial NH4 demand in soil, [g d-2 h-1]
  real(r8),target,allocatable :: RNO3UptkSoilAutor_vr(:,:,:,:)         !autotrophic microbial NO3 demand in soil, [g d-2 h-1]
  real(r8),target,allocatable :: RH2PO4UptkSoilAutor_vr(:,:,:,:)       !autotrophic microbes H2PO4 demand in soil, [g d-2 h-1]
  real(r8),target,allocatable :: RNH4UptkLitrAutor_col(:,:,:)        !autotrophic microbial NH4 demand in surface litter, [g d-2 h-1]
  real(r8),target,allocatable :: RH2PO4UptkLitrAutor_col(:,:,:)      !autotrophic microbial H2PO4 demand in surface litter, [g d-2 h-1]
  real(r8),target,allocatable :: RNO3UptkLitrAutor_col(:,:,:)        !autotrophic microbial NO3 demand in surface litter, [g d-2 h-1]
  real(r8),target,allocatable :: RNH3OxidAutor_vr(:,:,:,:)             !autotrophic NH3 oxidation in non-band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RNO2OxidAutor_vr(:,:,:,:)            !autotrophic NO2 oxidation in non-band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RNH3OxidAutorBand_vr(:,:,:,:)       !autotrophic NH3 oxidation in band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RNO2OxidAutorBand_vr(:,:,:,:)      !autotrophic NO2 oxidation in band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RN2ODmndReduxAutor_vr(:,:,:,:)     !autotrophic microbial N2O Demand for reduction, [g d-2 h-1]
  real(r8),target,allocatable :: RNH4UptkBandAutor_vr(:,:,:,:)      !autotrophic microbial NH4 demand in band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RNO3UptkBandAutor_vr(:,:,:,:)      !autotrophic microbial NO3 demand in band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RH2PO4UptkBandAutor_vr(:,:,:,:)    !autotrophic microbial H2PO4 demand in band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RH1PO4UptkSoilAutor_vr(:,:,:,:)    !autotrophic microbial H1PO4 demand in non-band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RH1PO4UptkBandAutor_vr(:,:,:,:)    !autotrophic microbial H1PO4 demand in band soil, [g d-2 h-1]
  real(r8),target,allocatable :: RH1PO4UptkLitrAutor_col(:,:,:)     !autotrophic microibal H1pO4 demand in surface litter, [g d-2 h-1]
  real(r8),target,allocatable :: OMEERauto_2D(:,:,:,:,:,:)             !autotrophic microbial biomass loss through erosion, [g d-2 h-1]

  private :: InitAllocate

  contains

  subroutine InitMicrobialData()

  implicit none

  call InitAllocate

  end subroutine InitMicrobialData

  subroutine InitAllocate

  implicit none
  allocate(mBiomeHeter_vr(NumPlantChemElms,NumLiveHeterBioms,1:jcplx,0:JZ,JY,JX));mBiomeHeter_vr=0._r8
  allocate(RO2DmndHetert(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX)); RO2DmndHetert=0._r8
  allocate(RDOCUptkHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RDOCUptkHeter_vr=0._r8
  allocate(RAcetateUptkHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RAcetateUptkHeter_vr=0._r8
  allocate(RNH4DmndSoilHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RNH4DmndSoilHeter_vr=0._r8
  allocate(RNO3DmndSoilHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RNO3DmndSoilHeter_vr=0._r8
  allocate(RH2PO4DmndSoilHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RH2PO4DmndSoilHeter_vr=0._r8
  allocate(RNH4DmndLitrHeter_col(NumHetetr1MicCmplx,1:jcplx,JY,JX));RNH4DmndLitrHeter_col=0._r8
  allocate(RH2PO4DmndLitrHeter_col(NumHetetr1MicCmplx,1:jcplx,JY,JX));RH2PO4DmndLitrHeter_col=0._r8
  allocate(RNO3DmndLitrHeter_col(NumHetetr1MicCmplx,1:jcplx,JY,JX));RNO3DmndLitrHeter_col=0._r8
  allocate(RNO3ReduxDmndSoilHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RNO3ReduxDmndSoilHeter_vr=0._r8
  allocate(RNO2DmndReduxSoilHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RNO2DmndReduxSoilHeter_vr=0._r8
  allocate(RNO3ReduxDmndBandHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RNO3ReduxDmndBandHeter_vr=0._r8
  allocate(RNO2DmndReduxBandHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RNO2DmndReduxBandHeter_vr=0._r8
  allocate(RN2ODmndReduxHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RN2ODmndReduxHeter_vr=0._r8
  allocate(RNH4DmndBandHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RNH4DmndBandHeter_vr=0._r8
  allocate(RNO3DmndBandHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RNO3DmndBandHeter_vr=0._r8
  allocate(RH2PO4DmndBandHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RH2PO4DmndBandHeter_vr=0._r8
  allocate(RH1PO4DmndSoilHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RH1PO4DmndSoilHeter_vr=0._r8
  allocate(RH1PO4DmndBandHeter_vr(NumHetetr1MicCmplx,1:jcplx,0:JZ,JY,JX));RH1PO4DmndBandHeter_vr=0._r8
  allocate(RH1PO4DmndLitrHeter_col(NumHetetr1MicCmplx,1:jcplx,JY,JX));RH1PO4DmndLitrHeter_col=0._r8
  allocate(OMEERhetr_2D(NumPlantChemElms,NumLiveHeterBioms,1:jcplx,2,2,JV,JH));OMEERhetr_2D=0._r8
  allocate(mBiomeAutor_vr(NumPlantChemElms,NumLiveAutoBioms,0:JZ,JY,JX));mBiomeAutor_vr=0._r8
  allocate(RO2DmndAutort_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RO2DmndAutort_vr=0._r8
  allocate(RNH4UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNH4UptkSoilAutor_vr=0._r8
  allocate(RNO3UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNO3UptkSoilAutor_vr=0._r8
  allocate(RH2PO4UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RH2PO4UptkSoilAutor_vr=0._r8
  allocate(RNH4UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX));RNH4UptkLitrAutor_col=0._r8
  allocate(RH2PO4UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX));RH2PO4UptkLitrAutor_col=0._r8
  allocate(RNO3UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX));RNO3UptkLitrAutor_col=0._r8
  allocate(RNH3OxidAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNH3OxidAutor_vr=0._r8
  allocate(RNO2OxidAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNO2OxidAutor_vr=0._r8
  allocate(RNH3OxidAutorBand_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNH3OxidAutorBand_vr=0._r8
  allocate(RNO2OxidAutorBand_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNO2OxidAutorBand_vr=0._r8
  allocate(RN2ODmndReduxAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RN2ODmndReduxAutor_vr=0._r8
  allocate(RNH4UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNH4UptkBandAutor_vr=0._r8
  allocate(RNO3UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RNO3UptkBandAutor_vr=0._r8
  allocate(RH2PO4UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RH2PO4UptkBandAutor_vr=0._r8
  allocate(RH1PO4UptkSoilAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RH1PO4UptkSoilAutor_vr=0._r8
  allocate(RH1PO4UptkBandAutor_vr(NumMicrobAutrophCmplx,0:JZ,JY,JX));RH1PO4UptkBandAutor_vr=0._r8
  allocate(RH1PO4UptkLitrAutor_col(NumMicrobAutrophCmplx,JY,JX));RH1PO4UptkLitrAutor_col=0._r8
  allocate(OMEERauto_2D(NumPlantChemElms,NumLiveAutoBioms,2,2,JV,JH));OMEERauto_2D=0._r8
  end subroutine InitAllocate
!----------------------------------------------------------------------------------------------

  subroutine DestructMicrobialData
  implicit none

  call destroy(mBiomeHeter_vr)
  call destroy(RO2DmndHetert)
  call destroy(RDOCUptkHeter_vr)
  call destroy(RAcetateUptkHeter_vr)
  call destroy(RNH4DmndSoilHeter_vr)
  call destroy(RNO3DmndSoilHeter_vr)
  call destroy(RH2PO4DmndSoilHeter_vr)
  call destroy(RNH4DmndLitrHeter_col)
  call destroy(RH2PO4DmndLitrHeter_col)
  call destroy(RNO3DmndLitrHeter_col)
  call destroy(RNO3ReduxDmndSoilHeter_vr)
  call destroy(RNO2DmndReduxSoilHeter_vr)
  call destroy(RNO3ReduxDmndBandHeter_vr)
  call destroy(RNO2DmndReduxBandHeter_vr)
  call destroy(RN2ODmndReduxHeter_vr)
  call destroy(RNH4DmndBandHeter_vr)
  call destroy(RNO3DmndBandHeter_vr)
  call destroy(RH2PO4DmndBandHeter_vr)
  call destroy(RH1PO4DmndSoilHeter_vr)
  call destroy(RH1PO4DmndBandHeter_vr)
  call destroy(RH1PO4DmndLitrHeter_col)
  call destroy(OMEERhetr_2D)
  call destroy(mBiomeAutor_vr)
  call destroy(RO2DmndAutort_vr)
  call destroy(RNH4UptkSoilAutor_vr)
  call destroy(RNO3UptkSoilAutor_vr)
  call destroy(RH2PO4UptkSoilAutor_vr)
  call destroy(RNH4UptkLitrAutor_col)
  call destroy(RH2PO4UptkLitrAutor_col)
  call destroy(RNO3UptkLitrAutor_col)
  call destroy(RNH3OxidAutor_vr)
  call destroy(RNO2OxidAutor_vr)
  call destroy(RNH3OxidAutorBand_vr)
  call destroy(RNO2OxidAutorBand_vr)
  call destroy(RN2ODmndReduxAutor_vr)
  call destroy(RNH4UptkBandAutor_vr)
  call destroy(RNO3UptkBandAutor_vr)
  call destroy(RH2PO4UptkBandAutor_vr)
  call destroy(RH1PO4UptkSoilAutor_vr)
  call destroy(RH1PO4UptkBandAutor_vr)
  call destroy(RH1PO4UptkLitrAutor_col)
  call destroy(OMEERauto_2D)
  end subroutine DestructMicrobialData

end module MicrobialDataType
