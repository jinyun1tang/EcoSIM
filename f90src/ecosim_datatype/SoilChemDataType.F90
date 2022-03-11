module SoilChemDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
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
  real(r8),allocatable :: H2GG(:,:,:)
  real(r8),allocatable :: CH2GG(:,:,:)
  private :: InitAllocate

  contains


  subroutine InitSoilChemData

  implicit none

  call InitAllocate

  end subroutine InitSoilChemData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none




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
  allocate(COXYG(0:JZ,JY,JX));COXYG=0._r8
  allocate(CCH4G(0:JZ,JY,JX));CCH4G=0._r8
  allocate(COXYS(0:JZ,JY,JX));COXYS=0._r8
  allocate(CCO2G(0:JZ,JY,JX));CCO2G=0._r8
  allocate(CCO2S(0:JZ,JY,JX));CCO2S=0._r8
  allocate(CCH4S(0:JZ,JY,JX));CCH4S=0._r8
  allocate(CH1P4(0:JZ,JY,JX))
  allocate(CH1P4B(0:JZ,JY,JX))
  allocate(CNO2B(0:JZ,JY,JX))
  allocate(H2GS(0:JZ,JY,JX))
  allocate(CH2GS(0:JZ,JY,JX))
  allocate(CH2P4(0:JZ,JY,JX))
  allocate(CH2P4B(0:JZ,JY,JX))
  allocate(H2GSH(JZ,JY,JX))
  allocate(H2GG(JZ,JY,JX))
  allocate(CH2GG(0:JZ,JY,JX))

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructSoilChemData

  implicit none

    deallocate(CZ2GS)
  deallocate(CNH4S)
  deallocate(CNH3S)
  deallocate(CNO3S)
  deallocate(CPO4S)
  deallocate(CNH4B)
  deallocate(CNH3B)
  deallocate(CNO3B)
  deallocate(CPO4B)
  deallocate(CNO2S)
  deallocate(CNH3G)
  deallocate(CZ2GG)
  deallocate(CZ2OG)
  deallocate(CZ2OS)
  deallocate(OXYG)
  deallocate(OXYS)
  deallocate(OXYSH)
  deallocate(CO2G)
  deallocate(CO2S)
  deallocate(CO2SH)
  deallocate(CH4G)
  deallocate(CH4S)
  deallocate(CH4SH)
  deallocate(COXYG)
  deallocate(CCH4G)
  deallocate(COXYS)
  deallocate(CCO2G)
  deallocate(CCO2S)
  deallocate(CCH4S)
  deallocate(CH1P4)
  deallocate(CH1P4B)
  deallocate(CNO2B)
  deallocate(H2GS)
  deallocate(CH2GS)
  deallocate(CH2P4)
  deallocate(CH2P4B)
  deallocate(H2GSH)
  deallocate(H2GG)
  deallocate(CH2GG)

  end subroutine DestructSoilChemData
end module SoilChemDataType
