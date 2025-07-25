module MicFluxTypeMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL     
  use EcoSiMParDataMod, only : micpar
  use abortutils, only : destroy
  use TracerIDMod
implicit none

  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  
  type, public :: micfluxtype

  real(r8) :: RCO2NetUptkMicb
  real(r8) :: RCH4UptkAutor
  real(r8) :: RH2NetUptkMicb
  real(r8) :: RO2UptkMicb
  real(r8) :: RN2NetUptkMicb
  real(r8) :: RN2ONetUptkMicb
  real(r8) :: RNH4MicbReliz2Soil
  real(r8) :: RNO3MicbReliz2Soil
  real(r8) :: RNO2MicbReliz2Soil
  real(r8) :: RH2PO4MicbReliz2Soil
  real(r8) :: RH1PO4MicbReliz2Soil
  real(r8) :: RNH4MicbReliz2Band
  real(r8) :: RNO3MicbReliz2Band
  real(r8) :: RNO2MicbReliz2Band
  real(r8) :: RH2PO4MicbReliz2Band
  real(r8) :: RH1PO4MicbReliz2Band
  real(r8) :: MicrbN2Fix                  !>0 N fixation
  real(r8) :: RNO2DmndSoilChemo
  real(r8) :: RNO2DmndBandChemo
  real(r8) :: CDOMuptk1                   !total DOC+acetate uptake for aerobic growth, [gC d-2 h-1]
  real(r8) :: tROMT
  real(r8) :: tRGOXP
  real(r8) :: tRGOZP
  real(r8) :: tGROMO
  real(r8) :: CDOMuptk2                   !total DOC+acetate uptake for denitrifcation, [gC d-2 h-1]
  real(r8) :: GrosAssimhr
  real(r8) :: NetCAssimhr                     !>0 uptake/ASSIMILATION
  real(r8) :: NetNH4Mineralize                !>0 uptake/ASSIMILATION
  real(r8) :: NetPO4Mineralize                !>0 uptake/ASSIMILATION
  real(r8) :: RNiDemand
  real(r8) :: RPiDemand  
  real(r8) :: TRDOE2DIE(1:NumPlantChemElms)              !cumulative conversion of organic element to inorganic element  
  real(r8) :: tRHydlySOM(1:NumPlantChemElms)
  real(r8) :: tRHydlyBioReSOM(1:NumPlantChemElms)
  real(r8) :: tRHydlySoprtOM(1:NumPlantChemElms)
  real(r8) :: RNO2DmndSoilChemoPrev
  real(r8) :: RNO2DmndBandChemoPrev  
  real(r8), allocatable :: RGrowthRespAutor(:)      !growth respiration of autotrophs
  real(r8), allocatable :: REcoDOMProd(:,:)
  real(r8), allocatable :: RO2DmndAutort(:)
  real(r8), allocatable :: RNH3OxidAutor(:)
  real(r8), allocatable :: RNH3OxidAutorBand(:)
  real(r8), allocatable :: RNO2OxidAutor(:)
  real(r8), allocatable :: RNO2OxidAutorBand(:)
  real(r8), allocatable :: RO2DmndHetert(:,:)
  real(r8), allocatable :: RDOCUptkHeter(:,:)             !potential DOC (unlimited) uptake flux, [gC d-2 h-1]
  real(r8), allocatable :: RAcetateUptkHeter(:,:)
  real(r8), allocatable :: RNO3ReduxDmndSoilHeter(:,:)
  real(r8), allocatable :: RNO3ReduxDmndBandHeter(:,:)
  real(r8), allocatable :: RNO2DmndReduxSoilHeter(:,:)
  real(r8), allocatable :: RNO2DmndReduxBandHeter(:,:)
  real(r8), allocatable :: RN2ODmndReduxHeter(:,:)
  real(r8), allocatable :: REcoUptkSoilO2M(:)
  real(r8), allocatable :: RNH4DmndSoilHeter(:,:)
  real(r8), allocatable :: RNH4DmndBandHeter(:,:)
  real(r8), allocatable :: RNO3DmndSoilHeter(:,:)
  real(r8), allocatable :: RNO3DmndBandHeter(:,:)
  real(r8), allocatable :: RH2PO4DmndSoilHeter(:,:)
  real(r8), allocatable :: RH2PO4DmndBandHeter(:,:)
  real(r8), allocatable :: RH1PO4DmndSoilHeter(:,:)
  real(r8), allocatable :: RH1PO4DmndBandHeter(:,:)
  real(r8), allocatable :: RNH4DmndLitrHeter(:,:)
  real(r8), allocatable :: RNO3DmndLitrHeter(:,:)
  real(r8), allocatable :: RH2PO4DmndLitrHeter(:,:)
  real(r8), allocatable :: RH1PO4DmndLitrHeter(:,:)
  real(r8), allocatable :: RNH4UptkSoilAutor(:)
  real(r8), allocatable :: RNH4UptkBandAutor(:)
  real(r8), allocatable :: RNO3UptkSoilAutor(:)
  real(r8), allocatable :: RNO3UptkBandAutor(:)
  real(r8), allocatable :: RH2PO4UptkSoilAutor(:)
  real(r8), allocatable :: RH2PO4UptkBandAutor(:)
  real(r8), allocatable :: RH1PO4UptkSoilAutor(:)
  real(r8), allocatable :: RH1PO4UptkBandAutor(:)
  real(r8), allocatable :: RNH4UptkLitrAutor(:)
  real(r8), allocatable :: RNO3UptkLitrAutor(:)
  real(r8), allocatable :: RH2PO4UptkLitrAutor(:)
  real(r8), allocatable :: RH1PO4UptkLitrAutor(:)
  real(r8), allocatable :: RHydlySOCK(:)            !solid organic carbon hydrolysis for each organic matter complex, [gC d-2]
  real(r8), allocatable :: RNH4UptkSoilAutorPrev(:)   
  real(r8), allocatable :: RNH4UptkBandAutorPrev(:)   
  real(r8), allocatable :: RNO3UptkSoilAutorPrev(:)   
  real(r8), allocatable :: RNO3UptkBandAutorPrev(:)   
  real(r8), allocatable :: RH2PO4UptkSoilAutorPrev(:) 
  real(r8), allocatable :: RH2PO4UptkBandAutorPrev(:) 
  real(r8), allocatable :: RH1PO4UptkSoilAutorPrev(:) 
  real(r8), allocatable :: RH1PO4UptkBandAutorPrev(:) 
  real(r8), allocatable :: RO2DmndAutortPrev(:)       

  real(r8), allocatable :: RNH4DmndSoilHeterPrev(:,:) 
  real(r8), allocatable :: RNH4DmndBandHeterPrev(:,:) 
  real(r8), allocatable :: RNO3DmndSoilHeterPrev(:,:) 
  real(r8), allocatable :: RNO3DmndBandHeterPrev(:,:)  
  real(r8), allocatable :: RH2PO4DmndSoilHeterPrev(:,:)
  real(r8), allocatable :: RH2PO4DmndBandHeterPrev(:,:)
  real(r8), allocatable :: RH1PO4DmndSoilHeterPrev(:,:)
  real(r8), allocatable :: RH1PO4DmndBandHeterPrev(:,:)
  real(r8), allocatable :: RO2DmndHetertPrev(:,:)      
  real(r8), allocatable :: RDOCUptkHeterPrev(:,:)      
  real(r8), allocatable :: RAcetateUptkHeterPrev(:,:)  
  real(r8), allocatable :: RNH3OxidAutorPrev(:) 
  real(r8), allocatable :: RNH3OxidAutorBandPrev(:)
  real(R8), allocatable :: RGrowthRespHeter(:,:)  !growth respiraiton of heterotrophs
  real(r8), allocatable :: RNH4DmndLitrHeterPrev(:,:);
  real(r8), allocatable :: RNO3DmndLitrHeterPrev(:,:)   
  real(r8), allocatable :: RH2PO4DmndLitrHeterPrev(:,:) 
  real(r8), allocatable :: RH1PO4DmndLitrHeterPrev(:,:) 
  real(r8), allocatable :: RNH4UptkLitrAutorPrev(:)     
  real(r8), allocatable :: RNO3UptkLitrAutorPrev(:)     
  real(r8), allocatable :: RH2PO4UptkLitrAutorPrev(:)   
  real(r8), allocatable :: RH1PO4UptkLitrAutorPrev(:)   
  real(r8), allocatable :: RNO2OxidAutorPrev(:)
  real(r8), allocatable :: RNO2OxidAutorBandPrev(:)
  real(r8), allocatable :: RNO3ReduxDmndSoilHeterPrev(:,:)
  real(r8), allocatable :: RNO3ReduxDmndBandHeterPrev(:,:)
  real(r8), allocatable :: RNO2DmndReduxSoilHeterPrev(:,:)
  real(r8), allocatable :: RNO2DmndReduxBandHeterPrev(:,:)
  real(r8), allocatable :: RN2ODmndReduxHeterPrev(:,:)
  contains
   procedure, public :: Init
   procedure, public :: Destroy=> Destruct
   procedure, public :: ZeroOut
  end type micfluxtype

  contains

  subroutine Init(this)

  implicit none
  class(micfluxtype) :: this
  integer :: jcplx,JG,NumMicbFunGrupsPerCmplx
  integer :: NumHetetr1MicCmplx, NumMicrobAutrophCmplx,NPH
  NPH=60
  jcplx                   = micpar%jcplx
  JG                      = micpar%jguilds
  NumMicbFunGrupsPerCmplx = micpar%NumMicbFunGrupsPerCmplx
  NumHetetr1MicCmplx      = micpar%NumHetetr1MicCmplx
  NumMicrobAutrophCmplx   = micpar%NumMicrobAutrophCmplx

  allocate(this%RGrowthRespHeter(1:NumHetetr1MicCmplx,1:jcplx));     this%RGrowthRespHeter=0._r8
  allocate(this%RNH4DmndSoilHeterPrev(1:NumHetetr1MicCmplx,1:jcplx)) ;this%RNH4DmndSoilHeterPrev=0._r8
  allocate(this%RNH4DmndBandHeterPrev(1:NumHetetr1MicCmplx,1:jcplx)) ;this%RNH4DmndBandHeterPrev=0._r8
  allocate(this%RNO3DmndSoilHeterPrev(1:NumHetetr1MicCmplx,1:jcplx)) ;this%RNO3DmndSoilHeterPrev=0._r8
  allocate(this%RNO3DmndBandHeterPrev(1:NumHetetr1MicCmplx,1:jcplx))  ;this%RNO3DmndBandHeterPrev=0._r8
  allocate(this%RH2PO4DmndSoilHeterPrev(1:NumHetetr1MicCmplx,1:jcplx));this%RH2PO4DmndSoilHeterPrev=0._r8
  allocate(this%RH2PO4DmndBandHeterPrev(1:NumHetetr1MicCmplx,1:jcplx));this%RH2PO4DmndBandHeterPrev=0._r8
  allocate(this%RH1PO4DmndSoilHeterPrev(1:NumHetetr1MicCmplx,1:jcplx));this%RH1PO4DmndSoilHeterPrev=0._r8
  allocate(this%RH1PO4DmndBandHeterPrev(1:NumHetetr1MicCmplx,1:jcplx));this%RH1PO4DmndBandHeterPrev=0._r8
  allocate(this%RO2DmndHetertPrev(1:NumHetetr1MicCmplx,1:jcplx))      ;this%RO2DmndHetertPrev=0._r8
  allocate(this%RDOCUptkHeterPrev(1:NumHetetr1MicCmplx,1:jcplx))      ;this%RDOCUptkHeterPrev=0._r8
  allocate(this%RAcetateUptkHeterPrev(1:NumHetetr1MicCmplx,1:jcplx))  ;this%RAcetateUptkHeterPrev=0._r8

  allocate(this%RNH4UptkSoilAutorPrev(1:NumMicrobAutrophCmplx));this%RNH4UptkSoilAutorPrev=0._r8   
  allocate(this%RNH4UptkBandAutorPrev(1:NumMicrobAutrophCmplx));this%RNH4UptkBandAutorPrev=0._r8   
  allocate(this%RNO3UptkSoilAutorPrev(1:NumMicrobAutrophCmplx));  this%RNO3UptkSoilAutorPrev=0._r8
  allocate(this%RNO3UptkBandAutorPrev(1:NumMicrobAutrophCmplx)) ;this%RNO3UptkBandAutorPrev=0._r8  
  allocate(this%RH2PO4UptkSoilAutorPrev(1:NumMicrobAutrophCmplx)) ;this%RH2PO4UptkSoilAutorPrev=0._r8
  allocate(this%RH2PO4UptkBandAutorPrev(1:NumMicrobAutrophCmplx)) ;this%RH2PO4UptkBandAutorPrev=0._r8
  allocate(this%RH1PO4UptkSoilAutorPrev(1:NumMicrobAutrophCmplx)) ;this%RH1PO4UptkSoilAutorPrev=0._r8
  allocate(this%RH1PO4UptkBandAutorPrev(1:NumMicrobAutrophCmplx)) ;this%RH1PO4UptkBandAutorPrev=0._r8
  allocate(this%RO2DmndAutortPrev(1:NumMicrobAutrophCmplx))       ;this%RO2DmndAutortPrev=0._r8
  allocate(this%RGrowthRespAutor(1:NumMicrobAutrophCmplx));this%RGrowthRespAutor=0._r8
  allocate(this%RHydlySOCK(1:jcplx)); this%RHydlySOCK=0._r8
  allocate(this%REcoUptkSoilO2M(NPH));this%REcoUptkSoilO2M = spval
  allocate(this%REcoDOMProd(idom_beg:idom_end,1:jcplx));this%REcoDOMProd=spval
  allocate(this%RO2DmndHetert(NumHetetr1MicCmplx,1:jcplx));this%RO2DmndHetert=spval
  allocate(this%RDOCUptkHeter(NumHetetr1MicCmplx,1:jcplx));this%RDOCUptkHeter=spval
  allocate(this%RAcetateUptkHeter(NumHetetr1MicCmplx,1:jcplx));this%RAcetateUptkHeter=spval
  allocate(this%RNO3ReduxDmndSoilHeter(NumHetetr1MicCmplx,1:jcplx));this%RNO3ReduxDmndSoilHeter=spval
  allocate(this%RNO3ReduxDmndBandHeter(NumHetetr1MicCmplx,1:jcplx));this%RNO3ReduxDmndBandHeter=spval
  allocate(this%RNO2DmndReduxSoilHeter(NumHetetr1MicCmplx,1:jcplx));this%RNO2DmndReduxSoilHeter=spval
  allocate(this%RNO2DmndReduxBandHeter(NumHetetr1MicCmplx,1:jcplx));this%RNO2DmndReduxBandHeter=spval
  allocate(this%RN2ODmndReduxHeter(NumHetetr1MicCmplx,1:jcplx));this%RN2ODmndReduxHeter=spval
  allocate(this%RNH4DmndSoilHeter(NumHetetr1MicCmplx,1:jcplx));this%RNH4DmndSoilHeter=spval
  allocate(this%RNH4DmndBandHeter(NumHetetr1MicCmplx,1:jcplx));this%RNH4DmndBandHeter=spval
  allocate(this%RNO3DmndSoilHeter(NumHetetr1MicCmplx,1:jcplx));this%RNO3DmndSoilHeter=spval
  allocate(this%RNO3DmndBandHeter(NumHetetr1MicCmplx,1:jcplx));this%RNO3DmndBandHeter=spval
  allocate(this%RH2PO4DmndSoilHeter(NumHetetr1MicCmplx,1:jcplx));this%RH2PO4DmndSoilHeter=spval
  allocate(this%RH2PO4DmndBandHeter(NumHetetr1MicCmplx,1:jcplx));this%RH2PO4DmndBandHeter=spval
  allocate(this%RH1PO4DmndSoilHeter(NumHetetr1MicCmplx,1:jcplx));this%RH1PO4DmndSoilHeter=spval
  allocate(this%RH1PO4DmndBandHeter(NumHetetr1MicCmplx,1:jcplx));this%RH1PO4DmndBandHeter=spval
  allocate(this%RNH4DmndLitrHeter(NumHetetr1MicCmplx,1:jcplx));this%RNH4DmndLitrHeter=spval
  allocate(this%RNO3DmndLitrHeter(NumHetetr1MicCmplx,1:jcplx));this%RNO3DmndLitrHeter=spval
  allocate(this%RH2PO4DmndLitrHeter(NumHetetr1MicCmplx,1:jcplx));this%RH2PO4DmndLitrHeter=spval
  allocate(this%RH1PO4DmndLitrHeter(NumHetetr1MicCmplx,1:jcplx));this%RH1PO4DmndLitrHeter=spval
  allocate(this%RO2DmndAutort(NumMicrobAutrophCmplx));this%RO2DmndAutort=spval
  allocate(this%RNH4UptkSoilAutor(NumMicrobAutrophCmplx));this%RNH4UptkSoilAutor=spval
  allocate(this%RNH4UptkBandAutor(NumMicrobAutrophCmplx));this%RNH4UptkBandAutor=spval
  allocate(this%RNO3UptkSoilAutor(NumMicrobAutrophCmplx));this%RNO3UptkSoilAutor=spval
  allocate(this%RNO3UptkBandAutor(NumMicrobAutrophCmplx));this%RNO3UptkBandAutor=spval
  allocate(this%RH2PO4UptkSoilAutor(NumMicrobAutrophCmplx));this%RH2PO4UptkSoilAutor=spval
  allocate(this%RH2PO4UptkBandAutor(NumMicrobAutrophCmplx));this%RH2PO4UptkBandAutor=spval
  allocate(this%RH1PO4UptkSoilAutor(NumMicrobAutrophCmplx));this%RH1PO4UptkSoilAutor=spval
  allocate(this%RH1PO4UptkBandAutor(NumMicrobAutrophCmplx));this%RH1PO4UptkBandAutor=spval
  allocate(this%RNH4UptkLitrAutor(NumMicrobAutrophCmplx));this%RNH4UptkLitrAutor=spval
  allocate(this%RNO3UptkLitrAutor(NumMicrobAutrophCmplx));this%RNO3UptkLitrAutor=spval
  allocate(this%RH2PO4UptkLitrAutor(NumMicrobAutrophCmplx));this%RH2PO4UptkLitrAutor=spval
  allocate(this%RH1PO4UptkLitrAutor(NumMicrobAutrophCmplx));this%RH1PO4UptkLitrAutor=spval
  allocate(this%RNH3OxidAutor(NumMicrobAutrophCmplx));this%RNH3OxidAutor=spval
  allocate(this%RNH3OxidAutorBand(NumMicrobAutrophCmplx));this%RNH3OxidAutorBand=spval
  allocate(this%RNO2OxidAutor(NumMicrobAutrophCmplx));this%RNO2OxidAutor=spval
  allocate(this%RNO2OxidAutorBand(NumMicrobAutrophCmplx));this%RNO2OxidAutorBand=spval

  allocate(this%RNH4DmndLitrHeterPrev(1:NumHetetr1MicCmplx,1:micpar%NumOfLitrCmplxs));this%RNH4DmndLitrHeterPrev=0._r8
  allocate(this%RNO3DmndLitrHeterPrev(1:NumHetetr1MicCmplx,1:micpar%NumOfLitrCmplxs));this%RNO3DmndLitrHeterPrev=0._r8
  allocate(this%RH2PO4DmndLitrHeterPrev(1:NumHetetr1MicCmplx,1:micpar%NumOfLitrCmplxs));this%RH2PO4DmndLitrHeterPrev=0._r8
  allocate(this%RH1PO4DmndLitrHeterPrev(1:NumHetetr1MicCmplx,1:micpar%NumOfLitrCmplxs));this%RH1PO4DmndLitrHeterPrev=0._r8
  allocate(this%RNH4UptkLitrAutorPrev(1:NumMicrobAutrophCmplx));this%RNH4UptkLitrAutorPrev=0._r8
  allocate(this%RNO3UptkLitrAutorPrev(1:NumMicrobAutrophCmplx));this%RNO3UptkLitrAutorPrev=0._r8
  allocate(this%RH2PO4UptkLitrAutorPrev(1:NumMicrobAutrophCmplx));this%RH2PO4UptkLitrAutorPrev=0._r8
  allocate(this%RH1PO4UptkLitrAutorPrev(1:NumMicrobAutrophCmplx));this%RH1PO4UptkLitrAutorPrev=0._r8

  allocate(this%RNO2OxidAutorPrev(1:NumMicrobAutrophCmplx)) ;this%RNO2OxidAutorPrev=0._r8     
  allocate(this%RNO2OxidAutorBandPrev(1:NumMicrobAutrophCmplx));this%RNO2OxidAutorBandPrev=0._r8
  allocate(this%RNH3OxidAutorPrev(1:NumMicrobAutrophCmplx)); this%RNH3OxidAutorPrev=0._r8
  allocate(this%RNH3OxidAutorBandPrev(1:NumMicrobAutrophCmplx));this%RNH3OxidAutorBandPrev=0._r8
  allocate(this%RNO3ReduxDmndSoilHeterPrev(1:NumHetetr1MicCmplx,1:jcplx));this%RNO3ReduxDmndSoilHeterPrev=0._r8
  allocate(this%RNO3ReduxDmndBandHeterPrev(1:NumHetetr1MicCmplx,1:jcplx));this%RNO3ReduxDmndBandHeterPrev=0._r8
  allocate(this%RNO2DmndReduxSoilHeterPrev(1:NumHetetr1MicCmplx,1:jcplx));this%RNO2DmndReduxSoilHeterPrev=0._r8
  allocate(this%RNO2DmndReduxBandHeterPrev(1:NumHetetr1MicCmplx,1:jcplx));this%RNO2DmndReduxBandHeterPrev=0._r8
  allocate(this%RN2ODmndReduxHeterPrev(1:NumHetetr1MicCmplx,1:jcplx));this%RN2ODmndReduxHeterPrev=0._r8
  end subroutine Init
!------------------------------------------------------------------------------------------
  subroutine ZeroOut(this)

  implicit none
  class(micfluxtype) :: this
  integer :: jcplx,JG,NumMicbFunGrupsPerCmplx

  this%RGrowthRespAutor       = 0._r8
  this%RGrowthRespHeter       = 0._r8
  this%TRDOE2DIE              = 0._r8
  this%REcoUptkSoilO2M        = 0._r8
  this%REcoDOMProd            = 0._r8
  this%RO2DmndHetert          = 0._r8
  this%RDOCUptkHeter          = 0._r8
  this%RAcetateUptkHeter      = 0._r8
  this%RNO3ReduxDmndSoilHeter = 0._r8
  this%RNO3ReduxDmndBandHeter = 0._r8
  this%RNO2DmndReduxSoilHeter = 0._r8
  this%RNO2DmndReduxBandHeter = 0._r8
  this%RN2ODmndReduxHeter     = 0._r8
  this%RNH4DmndSoilHeter      = 0._r8
  this%RNH4DmndBandHeter      = 0._r8
  this%RNO3DmndSoilHeter      = 0._r8
  this%RNO3DmndBandHeter      = 0._r8
  this%RH2PO4DmndSoilHeter    = 0._r8
  this%RH2PO4DmndBandHeter    = 0._r8
  this%RH1PO4DmndSoilHeter    = 0._r8
  this%RH1PO4DmndBandHeter    = 0._r8
  this%RNH4DmndLitrHeter      = 0._r8
  this%RNO3DmndLitrHeter      = 0._r8
  this%RH2PO4DmndLitrHeter    = 0._r8
  this%RH1PO4DmndLitrHeter    = 0._r8
  this%RO2DmndAutort          = 0._r8
  this%RNH4UptkSoilAutor      = 0._r8
  this%RNH4UptkBandAutor      = 0._r8
  this%RNO3UptkSoilAutor      = 0._r8
  this%RNO3UptkBandAutor      = 0._r8
  this%RH2PO4UptkSoilAutor    = 0._r8
  this%RH2PO4UptkBandAutor    = 0._r8
  this%RH1PO4UptkSoilAutor    = 0._r8
  this%RH1PO4UptkBandAutor    = 0._r8
  this%RNH4UptkLitrAutor      = 0._r8
  this%RNO3UptkLitrAutor      = 0._r8
  this%RH2PO4UptkLitrAutor    = 0._r8
  this%RH1PO4UptkLitrAutor    = 0._r8
  this%RNH3OxidAutor          = 0._r8
  this%RNH3OxidAutorBand      = 0._r8
  this%RNO2OxidAutor          = 0._r8
  this%RNO2OxidAutorBand      = 0._r8
  this%tRHydlySOM             = 0._r8
  this%tRHydlyBioReSOM        = 0._r8
  this%tRHydlySoprtOM         = 0._r8
  this%NetCAssimhr            = 0._r8
  this%NetNH4Mineralize       = 0._r8
  this%NetPO4Mineralize       = 0._r8
  this%RPiDemand              = 0._r8
  this%RNiDemand              = 0._r8
  this%GrosAssimhr            = 0._r8
  this%CDOMuptk1              = 0._r8
  this%CDOMuptk2              = 0._r8
  this%tROMT                  = 0._r8
  this%tGROMO                 = 0._r8
  this%tRGOXP                 = 0._r8
  this%tRGOZP                 = 0._r8
  this%RHydlySOCK             = 0._r8

  end subroutine ZeroOut
!------------------------------------------------------------------------------------------
  subroutine Destruct(this)
  use abortutils, only : destroy
  implicit none
  class(micfluxtype) :: this

  call destroy(this%RGrowthRespAutor)
  call destroy(this%RGrowthRespHeter)
  call destroy(this%RHydlySOCK)
  call destroy(this%REcoDOMProd)
  call destroy(this%RNH3OxidAutor)
  call destroy(this%RNH3OxidAutorBand)
  call destroy(this%RNO2OxidAutor)
  call destroy(this%RNO2OxidAutorBand)
  call destroy(this%RO2DmndHetert)
  call destroy(this%RDOCUptkHeter)
  call destroy(this%RAcetateUptkHeter)
  call destroy(this%RNO3ReduxDmndSoilHeter)
  call destroy(this%RNO3ReduxDmndBandHeter)
  call destroy(this%RNO2DmndReduxSoilHeter)
  call destroy(this%RNO2DmndReduxBandHeter)
  call destroy(this%RN2ODmndReduxHeter)
  call destroy(this%REcoUptkSoilO2M)
  call destroy(this%RNH4DmndSoilHeter)
  call destroy(this%RNH4DmndBandHeter)
  call destroy(this%RNO3DmndSoilHeter)
  call destroy(this%RNO3DmndBandHeter)
  call destroy(this%RH2PO4DmndSoilHeter)
  call destroy(this%RH2PO4DmndBandHeter)
  call destroy(this%RH1PO4DmndSoilHeter)
  call destroy(this%RH1PO4DmndBandHeter)
  call destroy(this%RNH4DmndLitrHeter)
  call destroy(this%RNO3DmndLitrHeter)
  call destroy(this%RH2PO4DmndLitrHeter)
  call destroy(this%RH1PO4DmndLitrHeter)
  call destroy(this%RNH4UptkSoilAutor)
  call destroy(this%RNH4UptkBandAutor)
  call destroy(this%RNO3UptkSoilAutor)
  call destroy(this%RNO3UptkBandAutor)
  call destroy(this%RH2PO4UptkSoilAutor)
  call destroy(this%RH2PO4UptkBandAutor)
  call destroy(this%RH1PO4UptkSoilAutor)
  call destroy(this%RH1PO4UptkBandAutor)
  call destroy(this%RNH4UptkLitrAutor)
  call destroy(this%RNO3UptkLitrAutor)
  call destroy(this%RH2PO4UptkLitrAutor)
  call destroy(this%RH1PO4UptkLitrAutor)
  call destroy(this%RO2DmndAutort)

  end subroutine Destruct
end module MicFluxTypeMod
