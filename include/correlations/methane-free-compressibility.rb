declare_correlation_subtype("MethaneFreeCompressibility", "WaterCorrelation",
                            "Cg_{fw}") # TODO: Ixhel

begin_correlation("CgfwSpiveyMN", "mPa_1")
add_title("SPIVEY, McCAIN & NORTH, CALCULATION OF COMPRESSIBILITY OF METHANE-FREE BRINE")
add_db("Calculation of compressibility of methane-free brine.")
add_parameter("t", "Celsius", "Temperature") 
add_parameter("p", "mPascal", "Pressure")
add_parameter("nacl", "Molality_NaCl", "Dissolved salt concentration")
add_author("Spivey, McCain & North")
add_ref("spivey:2004")
add_ref("mcCain:2011")
add_internal_note("The correlation was verified by using the original reference and McCain et al. (2011) as a secondary reference. Date: October 19 2016.")
set_hidden()
end_correlation()

