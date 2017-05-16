# coding: utf-8

declare_correlation_subtype("PureWaterDensity", "WaterCorrelation", "\\\\rho_{w}");

begin_correlation("PpwSpiveyMN", "Gr_cm3")
add_parameter("t", "Celsius", "Temperature", 0, 275)
add_parameter("p", "mPascal", "Pressure", 0.1, 200)
add_title("SPIVEY, McCAIN & NORTH CORRELATION, CALCULATION OF PURE WATER DENSITY")
add_author("Spivey, McCain & North")
add_ref("spivey:2004")
add_note("Spivey, McCain & North (2004) used the IAPWS-95 international standard equation of state for water to determine the coefficients of the correlation for pure water.")
add_internal_note("The correlation was verified by using the original reference.")
set_hidden_blackoil_grid()
set_hidden_drygas_grid()
set_hidden_wetgas_grid()
end_correlation()
