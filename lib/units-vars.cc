
# include <pvt-units.H>

// the following data is declared in units.H
UnitItemTable PhysicalQuantity::tbl;

UnitItemTable Unit::tbl;
DynSetTree<const Unit * const> Unit::unit_tbl;

UnitHashTbl __unit_name_name_tbl;
UnitHashTbl __unit_name_symbol_tbl;
UnitHashTbl __unit_symbol_name_tbl;
UnitHashTbl __unit_symbol_symbol_tbl;
CompoundUnitTbl __compound_unit_tbl;

// The following global singleton variables are generated by extract-cv script
