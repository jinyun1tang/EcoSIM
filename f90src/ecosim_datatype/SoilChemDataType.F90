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
  real(r8),allocatable :: PH(:,:,:)
  real(r8),allocatable :: CEC(:,:,:)
  real(r8),allocatable :: AEC(:,:,:)

  private :: InitAllocate

  contains


  subroutine InitSoilChemData

  implicit none

  call InitAllocate

  end subroutine InitSoilChemData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate
  implicit none

  allocate(CZ2GS(0:JZ,JY,JX));CZ2GS(0:JZ,JY,JX)=0._r8
  allocate(CNH4S(0:JZ,JY,JX));CNH4S(0:JZ,JY,JX)=0._r8
  allocate(CNH3S(0:JZ,JY,JX));CNH3S(0:JZ,JY,JX)=0._r8
  allocate(CNO3S(0:JZ,JY,JX));CNO3S(0:JZ,JY,JX)=0._r8
  allocate(CPO4S(JZ,JY,JX));CPO4S(JZ,JY,JX)=0._r8
  allocate(CNH4B(0:JZ,JY,JX));CNH4B(0:JZ,JY,JX)=0._r8
  allocate(CNH3B(0:JZ,JY,JX));CNH3B(0:JZ,JY,JX)=0._r8
  allocate(CNO3B(0:JZ,JY,JX));CNO3B(0:JZ,JY,JX)=0._r8
  allocate(CPO4B(0:JZ,JY,JX));CPO4B(0:JZ,JY,JX)=0._r8
  allocate(CNO2S(0:JZ,JY,JX));CNO2S(0:JZ,JY,JX)=0._r8
  allocate(CNH3G(0:JZ,JY,JX));CNH3G(0:JZ,JY,JX)=0._r8
  allocate(CZ2GG(0:JZ,JY,JX));CZ2GG(0:JZ,JY,JX)=0._r8
  allocate(CZ2OG(0:JZ,JY,JX));CZ2OG(0:JZ,JY,JX)=0._r8
  allocate(CZ2OS(0:JZ,JY,JX));CZ2OS(0:JZ,JY,JX)=0._r8
  allocate(OXYG(JZ,JY,JX));OXYG(JZ,JY,JX)=0._r8
  allocate(OXYS(0:JZ,JY,JX));OXYS(0:JZ,JY,JX)=0._r8
  allocate(OXYSH(JZ,JY,JX));OXYSH(JZ,JY,JX)=0._r8
  allocate(CO2G(JZ,JY,JX));CO2G(JZ,JY,JX)=0._r8
  allocate(CO2S(0:JZ,JY,JX));CO2S(0:JZ,JY,JX)=0._r8
  allocate(CO2SH(JZ,JY,JX));CO2SH(JZ,JY,JX)=0._r8
  allocate(CH4G(JZ,JY,JX));CH4G(JZ,JY,JX)=0._r8
  allocate(CH4S(0:JZ,JY,JX));CH4S(0:JZ,JY,JX)=0._r8
  allocate(CH4SH(JZ,JY,JX));CH4SH(JZ,JY,JX)=0._r8
  allocate(COXYG(0:JZ,JY,JX));COXYG(0:JZ,JY,JX)=0._r8
  allocate(CCH4G(0:JZ,JY,JX));CCH4G(0:JZ,JY,JX)=0._r8
  allocate(COXYS(0:JZ,JY,JX));COXYS(0:JZ,JY,JX)=0._r8
  allocate(CCO2G(0:JZ,JY,JX));CCO2G(0:JZ,JY,JX)=0._r8
  allocate(CCO2S(0:JZ,JY,JX));CCO2S(0:JZ,JY,JX)=0._r8
  allocate(CCH4S(0:JZ,JY,JX));CCH4S(0:JZ,JY,JX)=0._r8
  allocate(CH1P4(0:JZ,JY,JX));CH1P4(0:JZ,JY,JX)=0._r8
  allocate(CH1P4B(0:JZ,JY,JX));CH1P4B(0:JZ,JY,JX)=0._r8
  allocate(CNO2B(0:JZ,JY,JX));CNO2B(0:JZ,JY,JX)=0._r8
  allocate(H2GS(0:JZ,JY,JX));H2GS(0:JZ,JY,JX)=0._r8
  allocate(CH2GS(0:JZ,JY,JX));CH2GS(0:JZ,JY,JX)=0._r8
  allocate(CH2P4(0:JZ,JY,JX));CH2P4(0:JZ,JY,JX)=0._r8
  allocate(CH2P4B(0:JZ,JY,JX));CH2P4B(0:JZ,JY,JX)=0._r8
  allocate(H2GSH(JZ,JY,JX));H2GSH(JZ,JY,JX)=0._r8
  allocate(H2GG(JZ,JY,JX));H2GG(JZ,JY,JX)=0._r8
  allocate(CH2GG(0:JZ,JY,JX));CH2GG(0:JZ,JY,JX)=0._r8
  allocate(PH(0:JZ,JY,JX));PH(0:JZ,JY,JX)=0._r8
  allocate(CEC(JZ,JY,JX));CEC(JZ,JY,JX)=0._r8
  allocate(AEC(JZ,JY,JX));AEC(JZ,JY,JX)=0._r8
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructSoilChemData

  implicit none

  if(allocated(CZ2GS))deallocate(CZ2GS)
  if(allocated(CNH4S))deallocate(CNH4S)
  if(allocated(CNH3S))deallocate(CNH3S)
  if(allocated(CNO3S))deallocate(CNO3S)
  if(allocated(CPO4S))deallocate(CPO4S)
  if(allocated(CNH4B))deallocate(CNH4B)
  if(allocated(CNH3B))deallocate(CNH3B)
  if(allocated(CNO3B))deallocate(CNO3B)
  if(allocated(CPO4B))deallocate(CPO4B)
  if(allocated(CNO2S))deallocate(CNO2S)
  if(allocated(CNH3G))deallocate(CNH3G)
  if(allocated(CZ2GG))deallocate(CZ2GG)
  if(allocated(CZ2OG))deallocate(CZ2OG)
  if(allocated(CZ2OS))deallocate(CZ2OS)
  if(allocated(OXYG))deallocate(OXYG)
  if(allocated(OXYS))deallocate(OXYS)
  if(allocated(OXYSH))deallocate(OXYSH)
  if(allocated(CO2G))deallocate(CO2G)
  if(allocated(CO2S))deallocate(CO2S)
  if(allocated(CO2SH))deallocate(CO2SH)
  if(allocated(CH4G))deallocate(CH4G)
  if(allocated(CH4S))deallocate(CH4S)
  if(allocated(CH4SH))deallocate(CH4SH)
  if(allocated(COXYG))deallocate(COXYG)
  if(allocated(CCH4G))deallocate(CCH4G)
  if(allocated(COXYS))deallocate(COXYS)
  if(allocated(CCO2G))deallocate(CCO2G)
  if(allocated(CCO2S))deallocate(CCO2S)
  if(allocated(CCH4S))deallocate(CCH4S)
  if(allocated(CH1P4))deallocate(CH1P4)
  if(allocated(CH1P4B))deallocate(CH1P4B)
  if(allocated(CNO2B))deallocate(CNO2B)
  if(allocated(H2GS))deallocate(H2GS)
  if(allocated(CH2GS))deallocate(CH2GS)
  if(allocated(CH2P4))deallocate(CH2P4)
  if(allocated(CH2P4B))deallocate(CH2P4B)
  if(allocated(H2GSH))deallocate(H2GSH)
  if(allocated(H2GG))deallocate(H2GG)
  if(allocated(CH2GG))deallocate(CH2GG)
  if(allocated(PH))deallocate(PH)
  if(allocated(CEC))deallocate(CEC)
  if(allocated(AEC))deallocate(AEC)
  end subroutine DestructSoilChemData
end module SoilChemDataType
