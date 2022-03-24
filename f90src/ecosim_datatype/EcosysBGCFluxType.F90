module EcosysBGCFluxType

!
! Ecosystm fluxes for C, N, and P budget
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: TCAN(JY,JX)                       !total net CO2 fixation
  real(r8) :: TRN(JY,JX)                        !ecosystem net radiation, [MJ d-2 h-1]
  real(r8) :: TLE(JY,JX)                        !ecosystem latent heat flux, [MJ d-2 h-1]
  real(r8) :: TSH(JY,JX)                        !ecosystem sensible heat flux, [MJ d-2 h-1]
  real(r8) :: TGH(JY,JX)                        !ecosystem storage heat flux, [MJ d-2 h-1]
  real(r8) :: TGPP(JY,JX)                       !ecosystem GPP, [g d-2 h-1]
  real(r8) :: TRAU(JY,JX)                       !ecosystem autotrophic respiration, [g d-2 h-1]
  real(r8) :: TNPP(JY,JX)                       !ecosystem NPP, [g d-2 h-1]
  real(r8) :: THRE(JY,JX)                       !ecosystem heterotrophic respiration, [g d-2 h-1]
  real(r8) :: XHVSTC(JY,JX)                     !ecosystem harvest C, [g d-2]
  real(r8) :: XHVSTN(JY,JX)                     !ecosystem harvest N, [g d-2]
  real(r8) :: XHVSTP(JY,JX)                     !ecosystem harvest P, [g d-2]
  real(r8) :: TRINH4(JY,JX)                     !total NH4 net mineraln (-ve) or immobiln (+ve)
  real(r8) :: TRIPO4(JY,JX)                     !total H2PO4 net mineraln (-ve) or immobiln (+ve)
  real(r8) :: GPP(JY,JX)                        !gross primary productivity, [g d-2 h-1]
  real(r8) :: TCCAN(JY,JX)                      !total net CO2 fixation
  real(r8) :: ZCSNC(JY,JX)                      !total litterfall C, [g d-2 h-1]
  real(r8) :: ZZSNC(JY,JX)                      !total litterfall N, [g d-2 h-1]
  real(r8) :: ZPSNC(JY,JX)                      !total litterfall P, [g d-2 h-1]
  real(r8) :: RECO(JY,JX)                       !ecosystem respiration, [g d-2 h-1]
  real(r8) :: TNBP(JY,JX)                       !total NBP, [g d-2]

  real(r8) :: RP14X(0:JZ,JY,JX)                 !HPO4 demand in non-band by all microbial,root,myco populations
  real(r8) :: RP14Y(0:JZ,JY,JX)                 !HPO4 demand in non-band by all microbial,root,myco populations
  real(r8) :: RP1BX(0:JZ,JY,JX)                 !HPO4 demand in band by all microbial,root,myco populations
  real(r8) :: RP1BY(0:JZ,JY,JX)                 !HPO4 demand in band by all microbial,root,myco populations

end module EcosysBGCFluxType
