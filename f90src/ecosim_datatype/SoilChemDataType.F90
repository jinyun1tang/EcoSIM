module SoilChemDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none

  save

  real(r8),allocatable :: CZ2GS(:,:,:)
  real(r8),allocatable :: CNH4S(:,:,:)
  real(r8),allocatable :: CNH3S(:,:,:)
  real(r8),allocatable :: CNO3S(:,:,:)
  real(r8),allocatable :: CPO4S(:,:,:)
  real(r8),allocatable :: CNH4B(:,:,:)
  real(r8),allocatable :: CNH3B(:,:,:)
  real(r8),allocatable :: CNO3B(:,:,:)
  real(r8),allocatable :: CPO4B(:,:,:)
  real(r8),allocatable :: CNO2S(:,:,:)
  real(r8),allocatable :: CNH3G(:,:,:)
  real(r8),allocatable :: CZ2GG(:,:,:)
  real(r8),allocatable :: CZ2OG(:,:,:)
  real(r8),allocatable :: CZ2OS(:,:,:)
  real(r8),allocatable :: OXYG(:,:,:)
  real(r8),allocatable :: OXYS(:,:,:)
  real(r8),allocatable :: OXYSH(:,:,:)
  real(r8),allocatable :: CO2G(:,:,:)
  real(r8),allocatable :: CO2S(:,:,:)
  real(r8),allocatable :: CO2SH(:,:,:)
  real(r8),allocatable :: CH4G(:,:,:)
  real(r8),allocatable :: CH4S(:,:,:)
  real(r8),allocatable :: CH4SH(:,:,:)
  real(r8),allocatable :: COXYG(:,:,:)
  real(r8),allocatable :: CCH4G(:,:,:)
  real(r8),allocatable :: COXYS(:,:,:)
  real(r8),allocatable :: CCO2G(:,:,:)
  real(r8),allocatable :: CCO2S(:,:,:)
  real(r8),allocatable :: CCH4S(:,:,:)
  real(r8),allocatable :: CH1P4(:,:,:)
  real(r8),allocatable :: CH1P4B(:,:,:)
  real(r8),allocatable :: CNO2B(:,:,:)
  real(r8),allocatable :: H2GS(:,:,:)
  real(r8),allocatable :: CH2GS(:,:,:)
  real(r8),allocatable :: CH2P4(:,:,:)
  real(r8),allocatable :: CH2P4B(:,:,:)
  real(r8),allocatable :: H2GSH(:,:,:)
  private :: InitAllocate

  contains


  subroutine InitSoilChemData

  implicit none

  call InitAllocate

  end subroutine InitSoilChemData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none

  include "parameters.h"


  allocate(CZ2GS(0:JZ,JY,JX))
  allocate(CNH4S(0:JZ,JY,JX))
  allocate(CNH3S(0:JZ,JY,JX))
  allocate(CNO3S(0:JZ,JY,JX))
  allocate(CPO4S(JZ,JY,JX))
  allocate(CNH4B(0:JZ,JY,JX))
  allocate(CNH3B(0:JZ,JY,JX))
  allocate(CNO3B(0:JZ,JY,JX))
  allocate(CPO4B(0:JZ,JY,JX))
  allocate(CNO2S(0:JZ,JY,JX))
  allocate(CNH3G(0:JZ,JY,JX))
  allocate(CZ2GG(0:JZ,JY,JX))
  allocate(CZ2OG(0:JZ,JY,JX))
  allocate(CZ2OS(0:JZ,JY,JX))
  allocate(OXYG(JZ,JY,JX))
  allocate(OXYS(0:JZ,JY,JX))
  allocate(OXYSH(JZ,JY,JX))
  allocate(CO2G(JZ,JY,JX))
  allocate(CO2S(0:JZ,JY,JX))
  allocate(CO2SH(JZ,JY,JX))
  allocate(CH4G(JZ,JY,JX))
  allocate(CH4S(0:JZ,JY,JX))
  allocate(CH4SH(JZ,JY,JX))
  allocate(COXYG(0:JZ,JY,JX))
  allocate(CCH4G(0:JZ,JY,JX))
  allocate(COXYS(0:JZ,JY,JX))
  allocate(CCO2G(0:JZ,JY,JX))
  allocate(CCO2S(0:JZ,JY,JX))
  allocate(CCH4S(0:JZ,JY,JX))
  allocate(CH1P4(0:JZ,JY,JX))
  allocate(CH1P4B(0:JZ,JY,JX))
  allocate(CNO2B(0:JZ,JY,JX))
  allocate(H2GS(0:JZ,JY,JX))
  allocate(CH2GS(0:JZ,JY,JX))
  allocate(CH2P4(0:JZ,JY,JX))
  allocate(CH2P4B(0:JZ,JY,JX))
  allocate(H2GSH(JZ,JY,JX))

  end subroutine InitAllocate

end module SoilChemDataType
