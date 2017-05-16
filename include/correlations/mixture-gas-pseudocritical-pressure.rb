# coding: utf-8

declare_correlation_subtype("MixtureGasPseudocriticalPressure",
              "GasCorrelation", "P_{pcM}")

################################################################

begin_correlation("PpcmKayMixingRule", "psia")
add_title("KAY'S MIXING RULE, CALCULATION OF PSEUDOCRITICAL PRESSURE OF THE WHOLE GAS MIXTURE")
add_parameter("ppchc", "psia", "Gas pseudocritical pressure of the hydrocarbon portion")
add_parameter("n2", "MoleFraction", "N2 concentration")
add_parameter("co2", "MoleFraction", "CO2 concentration")
add_parameter("h2s", "MoleFraction", "H2S concentration")
add_precondition("n2", "co2", "h2s")
add_note("The value of the pseudocritical pressure of the gas hydrocarbon portion is adjusted for nonhydrocarbon content on the basis of Kay's mixing rule.")
add_author("Kay (Mixing Rule)")
add_ref("kay:1936")
add_ref("whitson:2000")
add_ref("standing:1977")
# add_ref("banzer:1996") Secondary reference
add_internal_note("The original reference is not available. The correlation was verified by using secondary references: Whitson & Brulé (2000), Standing (1977) and Bánzer (1996). Date: September 12 2016.")
set_hidden_blackoil_grid()
set_hidden_drygas_grid()
set_hidden_wetgas_grid()
end_correlation()

################################################################

begin_correlation("AdjustedppcmWichertAziz", "psia")
add_title("WICHERT & AZIZ CORRELATION, CALCULATION OF PSEUDOCRITICAL PRESSURE OF THE WHOLE GAS MIXTURE")
add_parameter("ppcm", "psia", "Gas pseudocritical pressure of the mixture")
add_parameter("tpcm", "Rankine", "Gas pseudocritical temperature of the mixture")
add_parameter("co2", "MoleFraction", "CO2 concentration", "Quantity<MolePercent>(0)", "Quantity<MolePercent>(54.4)") 
add_parameter("h2s", "MoleFraction", "H2S concentration", "Quantity<MolePercent>(0)", "Quantity<MolePercent>(73.8)")
add_precondition("co2", "h2s")
add_note("The pseudocritical pressure of the mixture is adjusted by using the Wichert & Aziz correlation when the gas contains significant fractions of acid gases, specifically carbon dioxide and hydrogen sulfide.")
add_note("The authors indicate that the correlation has an average absolute error of 0.97 % over these ranges: 154 psia < p < 7,026 psia and 40 °F < T < 300 °F.")
add_internal_note("The original reference is not available. The correlation was verified by using a secondary reference: Standing (1977). Date: September 12 2016.")
add_internal_note("The development ranges were taken from PetroWiki (http://petrowiki.org/Real_gases).")
add_author("Wichert & Aziz (Mixture Correction)")
add_ref("wichert:1972")
add_ref("standing:1977")
add_ref("petroWiki:2016:6")
set_hidden_blackoil_grid()
set_hidden_drygas_grid()
set_hidden_wetgas_grid()
end_correlation()

################################################################

begin_correlation("TpcmKayMixingRule", "MixtureGasPseudocriticalTemperature", "Rankine")
add_title("KAY'S MIXING RULE, CALCULATION OF PSEUDOCRITICAL TEMPERATURE OF THE WHOLE GAS MIXTURE")
add_parameter("tpchc", "Rankine", "Gas pseudocritical temperature of the hydrocarbon portion")
add_parameter("n2", "MoleFraction", "N2 concentration")
add_parameter("co2", "MoleFraction", "CO2 concentration")
add_parameter("h2s", "MoleFraction", "H2S concentration")
add_precondition("n2", "co2", "h2s")
add_note("The value of the pseudocritical temperature of the gas hydrocarbon portion is adjusted for nonhydrocarbon content on the basis of Kay's mixing rule.")
add_internal_note("The original reference is not available. The correlation was verified by using secondary references: Whitson & Brulé (2000), Standing (1977) and Bánzer (1996). Date: September 12 2016.")
add_ref("kay:1936")
add_ref("whitson:2000")
add_ref("standing:1977")
# add_ref("banzer:1996") Secondary reference
add_author("Kay (Mixing Rule)")
set_hidden_blackoil_grid()
set_hidden_drygas_grid()
set_hidden_wetgas_grid()
end_correlation()

################################################################

begin_correlation("AdjustedtpcmWichertAziz", "MixtureGasPseudocriticalTemperature",
         "Rankine")
add_title("WICHERT & AZIZ CORRELATION, CALCULATION OF PSEUDOCRITICAL TEMPERATURE OF THE WHOLE GAS MIXTURE")
add_parameter("tpcm", "Rankine", "Gas pseudocritical temperature of the mixture") 
add_parameter("co2", "MoleFraction", "CO2 concentration", "Quantity<MolePercent>(0)", "Quantity<MolePercent>(54.4)") 
add_parameter("h2s", "MoleFraction", "H2S concentration", "Quantity<MolePercent>(0)", "Quantity<MolePercent>(73.8)")
add_precondition("co2", "h2s")
add_note("The pseudocritical temperature of the mixture is adjusted by using the Wichert & Aziz correlation when the gas contains significant fractions of acid gases, specifically carbon dioxide and hydrogen sulfide.")
add_note("The authors indicate that the correlation has an average absolute error of 0.97 % over this range: 40 °F < T < 300 °F.")
add_ref("wichert:1972")
add_ref("standing:1977")
add_ref("petroWiki:2016:6")
add_internal_note("The original reference is not available. The correlation was verified by using a secondary reference: Standing (1977). Date: September 12 2016.")
add_internal_note("The development ranges were taken from PetroWiki (http://petrowiki.org/Real_gases).")
add_author("Wichert & Aziz (Mixture Correction)")
set_hidden_blackoil_grid()
set_hidden_drygas_grid()
set_hidden_wetgas_grid()
end_correlation()
