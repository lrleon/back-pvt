# coding: utf-8

declare_correlation_subtype("GasViscosity", "GasCorrelation", "\\\\mu_{g}")

begin_correlation("UgCarrKB", "CP")
add_title("CARR, KOBAYASHI & BURROWS CORRELATION, CALCULATION OF GAS VISCOSITY")
add_parameter("t", "Fahrenheit", "Temperature", 100, 300)
add_parameter("tpr", "PseudoReducedTemperature", "Gas pseudoreduced temperature", 1.2, 3) 
add_parameter("ppr", "PseudoReducedPressure", "Gas pseudoreduced pressure", 1, 20)
add_parameter("yg", "Sgg", "Gas specific gravity", 0.55, 1.55)
add_parameter("n2", "MoleFraction", "N2 concentration")
add_parameter("co2", "MoleFraction", "CO2 concentration")
add_parameter("h2s", "MoleFraction", "H2S concentration")
add_precondition("n2", "co2", "h2s")
add_author("Carr, Kobayashi & Burrows")
add_db("Carr, Kobayashi & Burrows (1954) developed charts that are the most widely used for calculating the viscosity of natural gas from the pseudoreduced temperature and pressure.")
add_ref("carr:1954")
add_ref("dempsey:1965")
add_ref("standing:1977")
add_ref("petroWiki:2016:2")
# add_ref("banzer:1996") Secondary reference
add_note("Carr et al. (1954) presented the graphs of the correlation. Dempsey (1965) expressed them in a mathematical form.")
add_internal_note("A linear effect of concentration was assumed to apply over the concentration range from 0 to 15 mol percent of non-hydrocarbon components.")
add_internal_note("The equation was verified by using secondary references: Bánzer (1996) and Standing (1977). Date: September 29 2016.")
add_internal_note("The development ranges were presented by Bánzer (1996).")
add_internal_note("The description was obtained from PetroWiki.")
end_correlation()

################################################################

begin_correlation("UgLeeGE", "mP")
add_title("LEE, GONZALEZ & EAKIN CORRELATION, CALCULATION OF GAS VISCOSITY")
add_parameter("t", "Rankine", "Temperature", "Quantity<Fahrenheit>(100)", "Quantity<Fahrenheit>(340)") 
add_parameter("p", "psia", "Pressure", 100, 8000)
add_parameter("yg", "Sgg", "Gas specific gravity")
add_parameter("z", "Zfactor", "Gas compressibility factor")
add_author("Lee, Gonzalez & Eakin")
add_ref("lee:1966")
# add_ref("banzer:1996") Secondary reference
add_db("Based on experimental viscosity and density data of four natural gases. The samples were furnished by the Atlantic Richfield Co., the Continental Oil Co. and the Pan American Petroleum Corp.")
add_internal_note("The correlation was verified by using the original reference and a secondary one: Bánzer (1996). Date: October 03 2016.")
end_correlation()

################################################################

begin_correlation("UgDeanStiel", "CP")
add_title("DEAN & STIEL CORRELATION, CALCULATION OF GAS VISCOSITY")
add_author("Dean & Stiel")
add_parameter("t", "Rankine", "Temperature")
add_parameter("p", "psia", "Pressure")
add_parameter("tpc", "Rankine", "Gas pseudocritical temperature")
add_parameter("ppc", "psia", "Gas pseudocritical pressure")
add_parameter("yg", "Sgg", "Gas specific gravity")
add_parameter("z", "Zfactor", "Gas compressibility factor")
add_ref("dean:1958")
add_ref("gawish:2005")
add_internal_note("The original reference is not available. The correlation was verified by using a secondary reference: Gawish & Al-Homadhi (2005). Date: October 03 2016.")
end_correlation()
