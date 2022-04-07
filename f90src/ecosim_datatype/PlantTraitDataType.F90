module PlantTraitDataType

!
!!
! data types of plant trait characteristics that cannot be grouped into canopy
! or roots
  use data_kind_mod, only : r8 => SHR_KIND_R8
  use GridConsts
  implicit none
  character(len=*), private, parameter :: mod_filename = __FILE__


!allocation parameter
  real(r8) :: FVRN(0:5)              !allocation parameter
  real(r8) :: FWOOD(0:1)             !woody C allocation
  real(r8) :: FWOODN(0:1)            !woody N allocation
  real(r8) :: FWOODP(0:1)            !woody P allocation
  REAL(R8) :: FWODB(0:1)             !woody C allocation
  real(r8) :: FWODLN(0:1)            !allocation
  real(r8) :: FWODLP(0:1)
  REAL(R8) :: FWODSN(0:1)            !N woody fraction in petiole
  real(r8) :: FWODSP(0:1)            !P woody fraction in petiole
  real(r8) :: FWODR(0:1)
  real(r8) :: FWODRN(0:1)
  real(r8) :: FWODRP(0:1)

!physical structure, height and area
  real(r8) :: ARSTK(JC,JC,JP,JY,JX)             !stem layer area, [m2 d-2]
  real(r8) :: ARLFP(JP,JY,JX)                   !plant leaf area, [m2 d-2]
  real(r8) :: ARLFS(JP,JY,JX)                   !plant leaf area, [m2 d-2]
  real(r8) :: ARLFX(JY,JX)                      !total canopy leaf area, [m2 d-2]
  real(r8) :: ARSTP(JP,JY,JX)                   !plant stem area, [m2 d-2]
  real(r8) :: ZC(JP,JY,JX)                      !canopy height, [m]
  real(r8) :: ARSTX(JY,JX)                      !total canopy stem area, [m2 d-2]
  real(r8) :: ARLFT(JC,JY,JX)                   !total leaf area, [m2 d-2]
  real(r8) :: ARSTT(JC,JY,JX)                   !total stem area, [m2 d-2]
  real(r8) :: ARLFC(JY,JX)                      !total canopy leaf area, [m2 d-2]
  real(r8) :: ARSTC(JY,JX)                      !total canopy stem area, [m2 d-2]
  real(r8) :: ARLSS(JY,JX)                      !stalk area of combined,each PFT canopy
  integer :: NG(JP,JY,JX)                      !soil layer at planting depth, [-]
  real(r8) :: SDPTHI(JP,JY,JX)                  !planting depth, [m]
  real(r8) :: SDPTH(JP,JY,JX)                   !seeding depth, [m]
  real(r8) :: SDVL(JP,JY,JX)                    !seed volume, [m3 ]
  real(r8) :: SDLG(JP,JY,JX)                    !seed length, [m]
  real(r8) :: SDAR(JP,JY,JX)                    !seed surface area, [m2]
  real(r8) :: HTCTL(JP,JY,JX)                   !cotyledon height, [m]
  real(r8) :: ZT(JY,JX)                         !canopy height , [m]
  real(r8) :: ZL(0:JC,JY,JX)                    !canopy layer height , [m]
  real(r8) :: ANGBR(JP,JY,JX)                   !branching angle, [degree from horizontal]
  real(r8) :: ANGSH(JP,JY,JX)                   !sheath angle, [degree from horizontal]
  real(r8) :: RAZ(JP,JY,JX)                     !canopy roughness height, [m]
  real(r8) :: HTSTZ(JP,JY,JX)                   !canopy height, [m]
  real(r8) :: ARLF(0:JNODS,JC,JP,JY,JX)            !leaf area, [m2 d-2]
  real(r8) :: HTSHE(0:JNODS,JC,JP,JY,JX)           !sheath height, [m]
  real(r8) :: HTNODE(0:JNODS,JC,JP,JY,JX)          !internode height, [m]
  real(r8) :: ARLFB(JC,JP,JY,JX)                !branch leaf area, [m2 d-2]
  real(r8) :: ARLFZ(JC,JP,JY,JX)                !branch leaf area, [m2 d-2]
  real(r8) :: HTSHEX(JC,JP,JY,JX)               !branch height, [m]
  real(r8) :: GRNOB(JC,JP,JY,JX)                !branch grain number, [d-2]
  real(r8) :: GRNXB(JC,JP,JY,JX)                !branch potential grain number, [d-2]
  real(r8) :: GRNO(JP,JY,JX)                    !canopy grain number, [d-2]
  real(r8) :: PP(JP,JY,JX)                      !plant population, [d-2]
  real(r8) :: HTNODX(0:JNODS,JC,JP,JY,JX)          !internode height, [m]

!stoichiometry
  real(r8) :: CNLF(JP,JY,JX)                    !maximum leaf N:C ratio, [g g-1]
  real(r8) :: CPLF(JP,JY,JX)                    !maximum leaf P:C ratio, [g g-1]
  real(r8) :: CNSHE(JP,JY,JX)                   !sheath N:C ratio, [g g-1]
  real(r8) :: CNSTK(JP,JY,JX)                   !stalk N:C ratio, [g g-1]
  real(r8) :: CNRSV(JP,JY,JX)                   !reserve N:C ratio, [g g-1]
  real(r8) :: CNHSK(JP,JY,JX)                   !husk N:C ratio, [g g-1]
  real(r8) :: CNEAR(JP,JY,JX)                   !ear N:C ratio, [g g-1]
  real(r8) :: CNGR(JP,JY,JX)                    !grain N:C ratio, [g g-1]
  real(r8) :: CNND(JP,JY,JX)                    !nodule N:C ratio, [g g-1]
  real(r8) :: CPSHE(JP,JY,JX)                   !sheath P:C ratio, [g g-1]
  real(r8) :: CPSTK(JP,JY,JX)                   !stalk P:C ratio, [g g-1]
  real(r8) :: CPRSV(JP,JY,JX)                   !reserve P:C ratio, [g g-1]
  real(r8) :: CPHSK(JP,JY,JX)                   !husk P:C ratio, [g g-1]
  real(r8) :: CPEAR(JP,JY,JX)                   !ear P:C ratio, [g g-1]
  real(r8) :: CPGR(JP,JY,JX)                    !grain P:C ratio, [g g-1]
  real(r8) :: CPND(JP,JY,JX)                    !nodule P:C ratio, [g g-1]
  real(r8) :: CNWS(JP,JY,JX)                    !C:N ratio in remobilizable nonstructural biomass, [-]
  real(r8) :: CPWS(JP,JY,JX)                    !C:P ratio in remobilizable nonstructural biomass, [-]

  real(r8) :: OSMO(JP,JY,JX)                    !canopy osmotic potential when canopy water potential = 0 MPa, [MPa]
  real(r8) :: TCX(JP,JY,JX)                     !threshold temperature for autumn leafoff/hardening, [oC]
  real(r8) :: ZTYPI(JP,JY,JX)                   !initial plant thermal adaptation zone, [-]
  real(r8) :: ZTYP(JP,JY,JX)                    !plant thermal adaptation zone, [-]
  real(r8) :: GROUP(JC,JP,JY,JX)                !plant maturity group, [-]
  real(r8) :: GROUPI(JP,JY,JX)                  !acclimated plant maturity group, [-]
  real(r8) :: GROUPX(JP,JY,JX)                  !initial plant maturity group, [-]
  real(r8) :: PPI(JP,JY,JX)                     !initial plant population, [m-2]
  real(r8) :: WTSTDI(JP,JY,JX)                  !initial standing dead C, [g C m-2]
  real(r8) :: PPX(JP,JY,JX)                     !plant population, [m-2]
  integer :: IFLGT(JY,JX)                       !number of active PFT
  real(r8) :: PPT(JY,JX)                        !total plant population, [d-2]
  real(r8) :: PPZ(JP,JY,JX)                     !plant population at seeding, [m-2]

  real(r8) :: WSTR(JP,JY,JX)                    !canopy plant water stress indicator, number of hours PSILT < PSILY, []
  real(r8) :: OSTR(JP,JY,JX)                    !plant O2 stress indicator, []
  real(r8) :: TFN3(JP,JY,JX)                    !canopy temperature growth function, [-]
  real(r8) :: TCG(JP,JY,JX)                     !canopy growth temperature, [oC]
  real(r8) :: TKG(JP,JY,JX)                     !canopy growth temperature, [K]

! yield

  real(r8) :: DMSHE(JP,JY,JX)                   !sheath growth yield, [g g-1]
  real(r8) :: DMSTK(JP,JY,JX)                   !stalk growth yield, [g g-1]
  real(r8) :: DMRSV(JP,JY,JX)                   !reserve growth yield, [g g-1]
  real(r8) :: DMHSK(JP,JY,JX)                   !husk growth yield, [g g-1]
  real(r8) :: DMEAR(JP,JY,JX)                   !ear growth yield, [g g-1]
  real(r8) :: DMGR(JP,JY,JX)                    !grain growth yield, [g g-1]
  real(r8) :: DMND(JP,JY,JX)                    !nodule growth yield, [g g-1]
  real(r8) :: DMLF(JP,JY,JX)                    !leaf growth yield, [g g-1]

!leaf/sheath/petiole
  real(r8) :: VRNY(JC,JP,JY,JX)                 !initial heat requirement for spring leafout/dehardening, [h]
  real(r8) :: VRNZ(JC,JP,JY,JX)                 !initial cold requirement for autumn leafoff/hardening, [h]
  real(r8) :: VSTG(JC,JP,JY,JX)                 !leaf number, [-]
  real(r8) :: VSTGX(JC,JP,JY,JX)                !leaf number at floral initiation, [-]
  real(r8) :: VRNS(JC,JP,JY,JX)                 !heat requirement for spring leafout/dehardening, [h]
  real(r8) :: VRNF(JC,JP,JY,JX)                 !cold requirement for autumn leafoff/hardening, [h]
  integer :: KLEAF(JC,JP,JY,JX)                 !leaf number, [-]
  integer :: KVSTGN(JC,JP,JY,JX)                !leaf growth stage counter, [-]
  real(r8) :: KLEAFX(JC,JP,JY,JX)               !NUMBER OF MINIMUM LEAFED NODE USED IN GROWTH ALLOCATION
  integer :: KVSTG(JC,JP,JY,JX)                 !leaf growth stage counter, [-]
  real(r8) :: XRLA(JP,JY,JX)                    !rate of leaf initiation, [h-1 at 25 oC]
  real(r8) :: WDLF(JP,JY,JX)                    !leaf length:width ratio, [-]
  real(r8) :: SLA1(JP,JY,JX)                    !leaf area:mass during growth, [m2 g-1]
  real(r8) :: TCZD                              !basal value for threshold temperature for spring leafout/dehardening	oC
  real(r8) :: TCXD                              !basal value for threshold temperature for autumn leafoff/hardening	oC
  real(r8) :: TCZ(JP,JY,JX)                     !threshold temperature for spring leafout/dehardening, [oC]
  real(r8) :: SSL1(JP,JY,JX)                    !petiole length:mass during growth, [m g-1]
  real(r8) :: VRNL(JC,JP,JY,JX)                 !hours above threshold temperature required for spring leafout/dehardening, [-]

!branch
  integer :: NBR(JP,JY,JX)                      !branch number, [-]
  integer :: NBT(JP,JY,JX)                      !branch number, [-]
  integer :: NBTB(JC,JP,JY,JX)                  !branch number, [-]
  integer :: NB1(JP,JY,JX)                      !number of main branch, [-]
  integer :: IFLGR(JC,JP,JY,JX)                 !branch phenology flag, [-]
  integer :: IFLGQ(JC,JP,JY,JX)                 !branch phenology flag, [h]
  integer :: IFLGG(JC,JP,JY,JX)                 !branch phenology flag, [-]
  integer :: IFLGP(JC,JP,JY,JX)                 !branch phenology flag, [-]
  integer :: IFLGA(JC,JP,JY,JX)                 !branch phenology flag, [-]
  integer :: IFLGE(JC,JP,JY,JX)                 !branch phenology flag, [-]
  integer :: IFLGF(JC,JP,JY,JX)                 !branch phenology flag, [-]
  integer :: IDTHB(JC,JP,JY,JX)                 !flag to detect branch death , [-]
  real(r8) :: PB(JP,JY,JX)                      !branch nonstructural C content required for new branch, [g g-1]

!node/stalk
  real(r8) :: GSTGI(JC,JP,JY,JX)                !normalized node number during vegetative growth stages , [-]
  real(r8) :: DGSTGI(JC,JP,JY,JX)               !gain in normalized node number during vegetative growth stages , [h-1]
  real(r8) :: DGSTGF(JC,JP,JY,JX)               !gain in normalized node number during reproductive growth stages, [h-1]
  real(r8) :: PSTG(JC,JP,JY,JX)                 !node number, [-]
  real(r8) :: PSTGI(JC,JP,JY,JX)                !node number at floral initiation, [-]
  real(r8) :: GSTGF(JC,JP,JY,JX)                !normalized node number during reproductive growth stages, [-]
  real(r8) :: PSTGF(JC,JP,JY,JX)                !node number at anthesis, [-]
  real(r8) :: TGSTGI(JC,JP,JY,JX)               !normalized node number during vegetative growth stages , [-]
  real(r8) :: TGSTGF(JC,JP,JY,JX)               !normalized node number during reproductive growth stages , [-]
  real(r8) :: XRNI(JP,JY,JX)                    !rate of node initiation, [h-1 at 25 oC]
  real(r8) :: SNL1(JP,JY,JX)                    !internode length:mass during growth, [m g-1]
  real(r8) :: FNOD(JP,JY,JX)                    !parameter for allocation of growth to nodes, [-]
  integer :: NNOD(JP,JY,JX)                     !number of concurrently growing nodes

!canopy
  real(r8) :: PSILZ(JP,JY,JX)                   !minimum daily canopy water potential, [MPa]
  real(r8) :: CFX(JP,JY,JX)                     !clumping factor for self-shading in canopy layer at current LAI, [-]
  real(r8) :: CF(JP,JY,JX)                      !clumping factor for self-shading in canopy layer, [-]
  integer :: IDTHP(JP,JY,JX)                    !flag to detect canopy death , [-]

!grain/seed
  real(r8) :: STMX(JP,JY,JX)                    !maximum grain node number per branch, [-]
  real(r8) :: SDMX(JP,JY,JX)                    !maximum grain number per node , [-]
  real(r8) :: GRMX(JP,JY,JX)                    !maximum grain size   , [g]
  real(r8) :: XTLI(JP,JY,JX)                    !number of nodes in seed, [-]
  real(r8) :: GRDM(JP,JY,JX)                    !grain size at seeding, [g]
  real(r8) :: GFILL(JP,JY,JX)                   !maximum rate of fill per grain, [g h-1]
!plant
  real(r8) :: FLG4(JC,JP,JY,JX)                 !flag to detect physiological maturity from  grain fill , [-]
  real(r8) :: ATRP(JC,JP,JY,JX)                 !counter for mobilizing nonstructural C during spring leafout/dehardening, [h]
  real(r8) :: FLGZ(JC,JP,JY,JX)                 !counter for mobilizing nonstructural C during autumn leafoff/hardening, [h]
  integer :: IDAY(10,JC,JP,JY,JX)               !plant growth stage, [-]
  real(r8) :: CTC(JP,JY,JX)                     !temperature below which seed set is adversely affected, [oC]
  real(r8) :: HTC(JP,JY,JX)                     !temperature above which seed set is adversely affected, [oC]
  real(r8) :: SSTX(JP,JY,JX)                    !sensitivity to HTC (seeds oC-1 above HTC)
  real(r8) :: XDL(JP,JY,JX)                     !critical daylength for phenological progress, [h]
  real(r8) :: XPPD(JP,JY,JX)                    !difference between current and critical daylengths used to calculate  phenological progress, [h]
  real(r8) :: CFI(JP,JY,JX)                     !initial clumping factor for self-shading in canopy layer, [-]
  real(r8) :: VRNX(JC,JP,JY,JX)                 !number of hours below set temperature required for autumn leafoff/hardening, [-]
  real(r8) :: OFFST(JP,JY,JX)                   !adjustment of Arhhenius curves for plant thermal acclimation, [oC]
  contains

  subroutine InitPlantTraits
  implicit none

  FVRN =real((/0.75,0.5,0.5,0.5,0.5,0.5/),r8)

  end subroutine InitPlantTraits
end module PlantTraitDataType
