# coding: utf-8

declare_correlation_subtype("SaturatedOilDensity", "OilFunction", "\\rho_{ob}");

# verificada con python!
begin_correlation("PobBradley", "Lb_ft3")
add_title("CALCULATION OF SATURATED OIL DENSITY")
add_parameter("yg", "Sgg", "Gas specific gravity")
add_parameter("rs", "SCF_STB", "Solution GOR")
add_parameter("bob", "RB_STB", "Saturated oil formation volume factor")
add_parameter("yo", "Sg_do", "Oil specific gravity")
add_synonym("yo", "api", "api")
add_author("Standard Equation")
add_ref("beggs:1987")
end_correlation()

################################################################    
