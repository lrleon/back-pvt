declare_correlation_subtype("VaporPressure", "WaterCorrelation", "P_v") # TODO: Ixhel

begin_correlation("PvSpiveyMN", "mPascal")
add_title("SPIVEY, McCAIN & NORTH, CALCULATION OF VAPOR PRESSURE")
add_db("Vapor pressure of pure water, calculated from the IAWPS-95 formulation.")
add_parameter("t", "Kelvin", "Temperature") 
add_author("Spivey, McCain & North")
add_ref("spivey:2004")
add_ref("mcCain:2011")
add_internal_note("The correlation was verified by using the original reference and McCain et al. (2011) as a secondary reference. Date: October 19 2016.")
set_hidden_blackoil_grid()
set_hidden_drygas_grid()
set_hidden_wetgas_grid()
end_correlation()
