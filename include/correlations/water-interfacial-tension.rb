# coding: utf-8

declare_correlation_subtype("WaterInterfacialTension", "WaterCorrelation",
			    "\\\\sigma_{gw}")

begin_correlation("SgwJenningsNewman", "dynes_cm")
add_title("JENNINGS & NEWMAN CORRELATION, CALCULATION OF GAS-WATER INTERFACIAL TENSION")
add_db("Based on measurements of interfacial tension of water against mixtures of methane and normal decane, under reservoir conditions of temperature and pressure.")
add_db("The interfacial tension data were obtained by the pendent drop method. Seven systems covering the methane-decane composition range from 100 percent decane to 100 percent methane were studied at three temperatures and 11 pressure intervals in the single-phase, hydrocarbon region.")
add_parameter("t", "Fahrenheit", "Temperature", 74, 350)
add_parameter("p", "psia", "Pressure", 14.7, 8000)
add_author("Jennings & Newman")
add_ref("jennings:1971")
add_ref("mcCain:1990")
# add_ref("banzer:1996") Secondary reference
add_internal_note("The data bank and development ranges were obtained from the original reference and a secondary one: McCain (1990).")
add_internal_note("Jennings & Newman (1971) presented a graphical correlation with a pressure range up to 12000 psia. McCain (1990) presented the correlation in a mathematical form, indicating that it should not be used for pressures above 8000 psia.")
add_internal_note("The correlation was verified by using Bánzer (1996) and McCain (1990) as secondary references. Date: October 13 2016.")
end_correlation()
