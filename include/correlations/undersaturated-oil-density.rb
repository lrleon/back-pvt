# coding: utf-8

declare_correlation_subtype("UndersaturatedOilDensity", "OilFunction",
                            "\\rho_{oa}")

# verificada con python!
begin_correlation("PoaBradley", "Lb_ft3")
add_title("CALCULATION OF UNDERSATURATED OIL DENSITY")
add_db("Calculation of oil density at pressures above the bubble point.")
add_parameter("pobp", "Lb_ft3", "Oil density at Pb")
add_parameter("p", "psia", "Pressure")
add_parameter("pb", "psia", "Bubble point pressure")
add_parameter("co", "psia_1", "Oil isothermal compressibility")
add_author("Standard Equation")
add_ref("beggs:1987")
add_precondition("p", "pb")
end_correlation()

################################################################    
