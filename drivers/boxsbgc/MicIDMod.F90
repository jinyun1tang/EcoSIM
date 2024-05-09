module MicIDMod
use ChemIDMod
implicit none
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  integer :: cidg_CO2      !gaseous tracer CO2
  integer :: cidg_CH4      !gaseous tracer CH4
  integer :: cid_OXYG      !gaseous tracer O2
  integer :: cid_Z2GG      !gaseous tracer N2
  integer :: cid_Z2OG      !gaseous tracer N2O
  integer :: cid_ZN3G      !gaseous tracer NH3
  integer :: cid_H2GG      !gaseous tracer H2
  integer :: cid_Z2GS      !aqueous N2
  integer :: cid_ZNH3G     !gaseous NH3 IN soil, [gN d-2]
  integer :: cid_ZNH3B     !NH3 band micropore, [gN d-2]
  integer :: cid_ZNH3S     !NH3 non-band micropore, [gN d-2]
  integer :: cid_ZNH4B     !NH4 band micropore, [g d-2]
  integer :: cid_ZNH4S     !NH4 non-band micropore, [g d-2]
  integer :: cid_H1POB     !soil aqueous HPO4 content micropore band, [g d-2]
  integer :: cid_H1PO4     !soil aqueous HPO4 content micropore non-band, [g d-2]
  integer :: cid_ZNO2B     !NO2(-)  band micropore, [g d-2]
  integer :: cid_ZNO2S     !NO2(-)  non-band micropore, [g d-2]
  integer :: cid_H2POB     !H2PO4 band micropore, [g d-2]
  integer :: cid_H2PO4     !H2PO4 non-band micropore, [g d-2]
  integer :: cid_CCO2S     !aqueous CO2 concentration micropore	[g m-3]
  integer :: cid_CNO2S     !NO2(-) concentration non-band micropore	[g m-3]
  integer :: cid_CNO2B     !NO2 concentration band micropore	[g m-3]
  integer :: cid_CZ2OS     !aqueous N2O concentration micropore	[g m-3]
  integer :: cid_Z2OS      !aqueous N2O mass in micropore, [g d-2]
  integer :: cid_COXYS     !aqueous O2 concentration micropore	[g m-3]
  integer :: cid_OXYS      !aqueous O2 mass in micropore	[g m-2]
  integer :: cid_COXYG     !gaseous O2 concentration	[g m-3]
  integer :: cid_CZ2GG     !gaseous N2 concentration [g m-3]
  integer :: cid_CZ2GS     !aqueous N2 concentration [g m-3]
  integer :: cid_CH2GS     !aqueous H2 concentration [g m-3]
  integer :: cid_H2GS      !aqueous H2 mass [g m-3]
  integer :: cid_CCH4G     !gaseous CH4 concentration in micropore [g m-3]
  integer :: cid_CH4S      !aqueous CH4 mass in  micropore	[g d-2]
  integer :: cid_ZNFN0     !initial nitrification inhibition activity
  integer :: cid_ZNFNI     !current nitrification inhibition activity
  integer :: cid_oqc_b,cid_oqc_e     !dissolved mass organic C micropore	[gC d-2]
  integer :: cid_oqn_b,cid_oqn_e     !dissolved mass organic N micropore	[gN d-2]
  integer :: cid_oqp_b,cid_oqp_e     !dissolved mass organic P micropore	[gP d-2]
  integer :: cid_oqa_b,cid_oqa_e     !dissolved mass acetate micropore [gC d-2]
  integer :: cid_ohc_b,cid_ohc_e     !adsorbed mass soil C	[gC d-2]
  integer :: cid_ohn_b,cid_ohn_e     !adsorbed mass soil N	[gN d-2]
  integer :: cid_ohp_b,cid_ohp_e     !adsorbed mass soil P	[gP d-2]
  integer :: cid_oha_b,cid_oha_e     !adsorbed mass soil acetate	[gC d-2]
  integer :: cid_osc_b,cid_osc_e     !humus mass soil C	[gC d-2]
  integer :: cid_osa_b,cid_osa_e     !humus mass soil acetate	[gC d-2]
  integer :: cid_osn_b,cid_osn_e     !humus mass soil N	[gN d-2]
  integer :: cid_osp_b,cid_osp_e     !humus mass soil P	[gP d-2]
  integer :: cid_orc_b,cid_orc_e     !microbial mass residue C [gC d-2]
  integer :: cid_orn_b,cid_orn_e     !microbial mass residue N [gN d-2]
  integer :: cid_orp_b,cid_orp_e     !microbial mass residue P [gP d-2]
  integer :: cid_mBiomeHeter_b,cid_mBiomeHeter_e     !microbial biomass component	[g d-2]
  integer :: cid_mBiomeAutor_b,cid_mBiomeAutor_e   !autotrophic microbial biomass component	[g d-2]

  integer :: fid_XCODFS             !CO2 dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XCHDFS             !CH4 dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XOXDFS             !O2  dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XNGDFS             !N2  dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XN2DFS             !N2O dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XN3DFS             !NH3 dissolution (+)-volatiziation (-) with respect to atmosphere in non-band soil
  integer :: fid_XNBDFS             !NH3 dissolution (+)-volatiziation (-) with respect to atmosphere in band soil
  integer :: fid_XHGDFS             !H2  dissolution (+)-volatiziation (-) with respect to atmosphere
  integer :: fid_XCODFG             !CO2 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XCHDFG             !CH4 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XOXDFG             !O2 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XNGDFG             !N2 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XN2DFG             !N2O dissolution (+)-volatiziation (-) in soil
  integer :: fid_XN3DFG             !NH3 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XNBDFG             !NH3 dissolution (+)-volatiziation (-) in band soil
  integer :: fid_XHGDFG             !H2 dissolution (+)-volatiziation (-) in soil
  integer :: fid_XCOFLG             !CO2 gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XCHFLG             !CH4 gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XOXFLG             !O2 gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XNGFLG             !N2 gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XN2FLG             !N2O gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XN3FLG             !N3H gaseous exchange with atmosphere (-) into atmosphere
  integer :: fid_XHGFLG             !H2 gaseous exchange with atmosphere (-) into atmosphere


  integer :: fid_RO2GasXchangePrev               !net gaseous O2 flux from previous hour, [g d-2 h-1]
  integer :: fid_RO2EcoDmndPrev               !total root + microbial O2 uptake, [g d-2 h-1]
  integer :: fid_RNH4EcoDmndSoilPrev               !total root + microbial NH4(+) uptake non-band, [gN d-2 h-1]
  integer :: fid_RNO3EcoDmndSoilPrev               !total root + microbial NO3(-) uptake non-band, [gN d-2 h-1]
  integer :: fid_RNO2EcoUptkSoilPrev               !total root + microbial NO2(-) uptake non-band, [gN d-2 h-1]
  integer :: fid_RN2OEcoUptkSoilPrev               !total root + microbial N2O uptake, [g d-2 h-1]
  integer :: fid_RH2PO4EcoDmndSoilPrev               !total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  integer :: fid_RH1PO4EcoDmndSoilPrev               !HPO4 demand in non-band by all microbial,root,myco populations [gP d-2 h-1]
  integer :: fid_RNH4EcoDmndBandPrev               !total root + microbial NH4 uptake band, [gN d-2 h-1]
  integer :: fid_RNO3EcoDmndBandPrev               !total root + microbial NO3(-) uptake band, [gN d-2 h-1]
  integer :: fid_RNO2EcoUptkBandPrev               !total root + microbial NO2(-) uptake band, [gN d-2 h-1]
  integer :: fid_RH2PO4EcoDmndBandPrev               !total root + microbial PO4 uptake band, [gP d-2 h-1]
  integer :: fid_RH1PO4EcoDmndBandPrev               !HPO4 demand in band by all microbial,root,myco populations [gP d-2 h-1]
  integer :: fid_RDOMEcoDmndPrev_b,fid_RDOMEcoDmndPrev_e !total root + microbial DOC uptake, [gC d-2 h-1]
  integer :: fid_RAcetateEcoDmndPrev_b,fid_RAcetateEcoDmndPrev_e !total root + microbial acetate uptake, [gC d-2 h-1]
  integer :: fid_RNH4DmndSoilHeter_b,fid_RNH4DmndSoilHeter_e
  integer :: fid_RNH4DmndBandHeter_b,fid_RNH4DmndBandHeter_e
  integer :: fid_RNO3DmndSoilHeter_b,fid_RNO3DmndSoilHeter_e
  integer :: fid_RNO3DmndBandHeter_b,fid_RNO3DmndBandHeter_e
  integer :: fid_RH2PO4DmndSoilHeter_b,fid_RH2PO4DmndSoilHeter_e
  integer :: fid_RH2PO4DmndBandHeter_b,fid_RH2PO4DmndBandHeter_e
  integer :: fid_RH1PO4DmndSoilHeter_b,fid_RH1PO4DmndSoilHeter_e
  integer :: fid_RH1PO4DmndBandHeter_b,fid_RH1PO4DmndBandHeter_e
  integer :: fid_RO2DmndHetert_b,fid_RO2DmndHetert_e
end module MicIDMod
