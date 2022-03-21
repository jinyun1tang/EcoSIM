module PhenologyDataType

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridDataType
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = __FILE__

  real(r8) :: FVRN(0:5)
  real(r8) :: FWOOD(0:1)
  real(r8) :: FWOODN(0:1)
  real(r8) :: FWOODP(0:1)
  REAL(R8) :: FWODB(0:1)
  real(r8) :: FWODLN(0:1)
  real(r8) :: FWODLP(0:1)
  REAL(R8) :: FWODSN(0:1)
  real(r8) :: FWODSP(0:1)
  real(r8) :: FWODR(0:1)
  real(r8) :: FWODRN(0:1)
  real(r8) :: FWODRP(0:1)

  real(r8) :: PSILZ(JP,JY,JX)                   !minimum daily canopy water potential, [MPa]
  real(r8) :: VRNY(JC,JP,JY,JX)                 !initial heat requirement for spring leafout/dehardening, [h]
  real(r8) :: VRNZ(JC,JP,JY,JX)                 !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8) :: GROUP(JC,JP,JY,JX)                !plant maturity group, [-]
  real(r8) :: FLG4(JC,JP,JY,JX)                 !flag to detect physiological maturity from  grain fill , [-]
  real(r8) :: GSTGI(JC,JP,JY,JX)                !normalized node number during vegetative growth stages , [-]
  real(r8) :: VSTG(JC,JP,JY,JX)                 !leaf number, [-]
  real(r8) :: VSTGX(JC,JP,JY,JX)                !leaf number at floral initiation, [-]
  real(r8) :: DGSTGI(JC,JP,JY,JX)               !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8) :: DGSTGF(JC,JP,JY,JX)               !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8) :: PSTG(JC,JP,JY,JX)                 !node number, [-]
  real(r8) :: PSTGI(JC,JP,JY,JX)                !node number at floral initiation, [-]
  real(r8) :: GSTGF(JC,JP,JY,JX)                !normalized node number during reproductive growth stages, [-]
  real(r8) :: PSTGF(JC,JP,JY,JX)                !node number at anthesis, [-]
  real(r8) :: TGSTGI(JC,JP,JY,JX)               !normalized node number during vegetative growth stages , [-]
  real(r8) :: TGSTGF(JC,JP,JY,JX)               !normalized node number during reproductive growth stages , [-]
  real(r8) :: VRNS(JC,JP,JY,JX)                 !heat requirement for spring leafout/dehardening, [h]
  real(r8) :: VRNF(JC,JP,JY,JX)                 !cold requirement for autumn leafoff/hardening, [h]
  real(r8) :: ATRP(JC,JP,JY,JX)                 !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8) :: FLGZ(JC,JP,JY,JX)                 !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]

  integer :: NBR(JP,JY,JX)                     !branch number, [-]
  integer :: NB1(JP,JY,JX)                     !number of main branch, [-]
  integer :: NRT(JP,JY,JX)                     !root primary axis number, [-]
  integer :: KVSTGN(JC,JP,JY,JX)               !leaf growth stage counter, [-]
  integer :: IDAY(10,JC,JP,JY,JX)              !plant growth stage, [-]
  integer :: NBT(JP,JY,JX)                     !branch number, [-]
  integer :: NBTB(JC,JP,JY,JX)                 !branch number, [-]
  integer :: NINR(JC,JP,JY,JX)                 !maximum soil layer number for roox axes, [-]
  integer :: KLEAF(JC,JP,JY,JX)                !leaf number, [-]
  integer :: IDTHR(JP,JY,JX)                   !flag to detect root system death , [-]
  integer :: IFLGR(JC,JP,JY,JX)                !branch phenology flag, [-]
  integer :: IFLGQ(JC,JP,JY,JX)                !branch phenology flag, [h]
  integer :: IFLGG(JC,JP,JY,JX)                !branch phenology flag, [-]
  integer :: IDTHP(JP,JY,JX)                   !flag to detect canopy death , [-]
  integer :: KVSTG(JC,JP,JY,JX)                !leaf growth stage counter, [-]
  integer :: IDTHB(JC,JP,JY,JX)                !flag to detect branch death , [-]
  integer :: NIX(JP,JY,JX)                     !maximum soil layer number for all root axes, [-]
  integer :: NI(JP,JY,JX)                      !maximum soil layer number for all root axes, [-]
  integer :: NG(JP,JY,JX)                      !soil layer at planting depth, [-]
  integer :: IFLGP(JC,JP,JY,JX)                !branch phenology flag, [-]
  integer :: IFLGA(JC,JP,JY,JX)                !branch phenology flag, [-]
  integer :: IFLGE(JC,JP,JY,JX)                !branch phenology flag, [-]
  integer :: IFLGF(JC,JP,JY,JX)                !branch phenology flag, [-]
  real(r8) :: KLEAFX(JC,JP,JY,JX)              !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION

  real(r8) :: CTC(JP,JY,JX)                     !temperature below which seed set is adversely affected, [oC]
  real(r8) :: XRNI(JP,JY,JX)                    !rate of node initiation, [h-1 at 25 oC]
  real(r8) :: XRLA(JP,JY,JX)                    !rate of leaf initiation, [h-1 at 25 oC]
  real(r8) :: WDLF(JP,JY,JX)                    !leaf length:width ratio, [-]
  real(r8) :: PB(JP,JY,JX)                      !branch nonstructural C content required for new branch, [g g-1]
  real(r8) :: SLA1(JP,JY,JX)                    !leaf area:mass during growth, [m2 g-1]
  real(r8) :: SSL1(JP,JY,JX)                    !petiole length:mass during growth, [m g-1]
  real(r8) :: SNL1(JP,JY,JX)                    !internode length:mass during growth, [m g-1]
  real(r8) :: VRNL(JC,JP,JY,JX)                 !hours above threshold temperature required for spring leafout/dehardening, [-]
  real(r8) :: HTC(JP,JY,JX)                     !temperature above which seed set is adversely affected, [oC]
  real(r8) :: SSTX(JP,JY,JX)                    !sensitivity to HTC (seeds oC-1 above HTC)
  real(r8) :: FNOD(JP,JY,JX)                    !parameter for allocation of growth to nodes, [-]
  real(r8) :: DMLF(JP,JY,JX)                    !leaf growth yield, [g g-1]
  real(r8) :: DMSHE(JP,JY,JX)                   !sheath growth yield, [g g-1]
  real(r8) :: DMSTK(JP,JY,JX)                   !stalk growth yield, [g g-1]
  real(r8) :: DMRSV(JP,JY,JX)                   !reserve growth yield, [g g-1]
  real(r8) :: DMHSK(JP,JY,JX)                   !husk growth yield, [g g-1]
  real(r8) :: DMEAR(JP,JY,JX)                   !ear growth yield, [g g-1]
  real(r8) :: DMGR(JP,JY,JX)                    !grain growth yield, [g g-1]
  real(r8) :: DMRT(JP,JY,JX)                    !root growth yield, [g g-1]
  real(r8) :: DMND(JP,JY,JX)                    !nodule growth yield, [g g-1]
  real(r8) :: XDL(JP,JY,JX)                     !critical daylength for phenological progress, [h]
  real(r8) :: XPPD(JP,JY,JX)                    !difference between current and critical daylengths used to calculate  phenological progress, [h]
  real(r8) :: CFI(JP,JY,JX)                     !initial clumping factor for self-shading in canopy layer, [-]
  real(r8) :: XTLI(JP,JY,JX)                    !number of nodels in seed, [-]
  real(r8) :: STMX(JP,JY,JX)                    !maximum grain node number per branch, [-]
  real(r8) :: SDMX(JP,JY,JX)                    !maximum grain number per node , [-]
  real(r8) :: GRMX(JP,JY,JX)                    !maximum grain size   , [g]
  real(r8) :: PPZ(JP,JY,JX)                     !plant population at seeding, [m-2]
  real(r8) :: GRDM(JP,JY,JX)                    !grain size at seeding, [g]
  real(r8) :: GFILL(JP,JY,JX)                   !maximum rate of fill per grain, [g h-1]
  real(r8) :: SDPTHI(JP,JY,JX)                  !planting depth, [m]
  real(r8) :: VRNX(JC,JP,JY,JX)                 !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8) :: CFX(JP,JY,JX)                     !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8) :: CF(JP,JY,JX)                      !clumping factor for self-shading in canopy layer, [-]
  integer :: NNOD(JP,JY,JX)                     !number of concurrently growing nodes
  real(r8) :: OFFST(JP,JY,JX)                   !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
  real(r8) :: TCZ(JP,JY,JX)                     !threshold temperature for spring leafout/dehardening, [oC]
  real(r8) :: PR(JP,JY,JX)                      !threshold root nonstructural C content for initiating new root axis, [g g-1]
  real(r8) :: TCZD                              !basal value for threshold temperature for spring leafout/dehardening	oC
  real(r8) :: TCXD                              !basal value for threshold temperature for autumn leafoff/hardening	oC

  contains
  subroutine InitPhenologyData

  implicit none


  FVRN =real((/0.75,0.5,0.5,0.5,0.5,0.5/),r8)
  end subroutine InitPhenologyData
end module PhenologyDataType
