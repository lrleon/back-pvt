# coding: utf-8

declare_correlation_subtype("HydrocarbonGasSpecificGravity", "GasCorrelation",
                            "P_{pcHC}");

begin_correlation("YghcWichertAziz", "Sgg", 0.55)
add_title("METHOD OF WICHERT & AZIZ FOR THE CALCULATION OF THE GAS GRAVITY OF THE HYDROCARBON PORTION")
add_parameter("yg", "Sgg", "Gas specific gravity")
add_parameter("n2", "MoleFraction", "N2 concentration")
add_parameter("co2", "MoleFraction", "CO2 concentration")
add_parameter("h2s", "MoleFraction", "H2S concentration")
add_precondition("n2", "co2", "h2s")
add_note("Method of Wichert & Aziz for the calculation of the gas gravity of the hydrocarbon portion.")
add_internal_note("The original reference is not available. The correlation was verified by using secondary references: Standing (1977) and Bánzer (1996). Date: September 12 2016.")
add_author("Wichert & Aziz (Gas Gravity Correction)")
add_ref("wichert:1972")
add_ref("standing:1977")
# add_ref("banzer:1996") Secondary reference
set_hidden_blackoil_grid()
set_hidden_drygas_grid()
set_hidden_wetgas_grid()
end_correlation()

