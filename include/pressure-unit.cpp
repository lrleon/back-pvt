

# ifndef PRESSURE_UNIT_H
# define PRESSURE_UNIT_H

# include <units.H>


Declare_Physical_Quantity(Pressure, "P",
              "Pressure is the force applied perpendicular to the surface of an object per unit area over which that force is distributed");

Declare_Unit(Bar, "bar",
         "A bar is a metric measurement unit of pressure. One bar is equivalent to ten newtons (N) per square centimeter (cm²)",
         Pressure, -1.0135326860, 2068.43);

Declare_Unit(Pascal, "Pa",
         "The pascal is the SI derived unit of pressure used to quantify internal pressure, stress, Young's modulus and
         ultimate tensile strength. It is defined as one newton per square metre.",
          Pressure, -101353.2686010438, 206843405.31);

Declare_Unit(kPascal, "kPa",
         "a kPascal is a common multiple unit of the pascal 1 kPa = 1000 Pa,",
          Pressure, -101.3532686010, 206843.41);

Declare_Unit(psia, "psia",
         "Pounds per square inch absolute (psia) is used to make it clear that the pressure is relative
         to a vacuum rather than the ambient atmospheric pressure. Since atmospheric pressure 
         at sea level is around 14.7 psi",

          Pressure, 0, 29985.3);

Declare_Unit(psig, "psig",
         "Pound-force per square inch is a unit of pressure or of stress based on avoirdupois units.
         It is the pressure resulting from a force of one pound-force applied to 
         an area of one square inch",
          Pressure, -14.7, 30000);

Declare_Unit(Atmophere, "atm",
         "Pound-force per square inch is a unit of pressure or of stress based on avoirdupois units.
         It is the pressure resulting from a force of one pound-force applied to 
         an area of one square inch",
          Pressure, -1, 2041.39);

// pressure conversions 
  
// To atmophere (atm)
Declare_Conversion(Bar, Atmophere, v) { return v * 0.986923267; }
Declare_Conversion(Pascal, Atmophere, v)  { return v * 9.86923267e-6; }
Declare_Conversion(kPascal, Atmophere, v)  { return v * 9.86923267e-3; }
Declare_Conversion(psia, Atmophere, v)  { return  v  / 14.69594877551; }
Declare_Conversion(psig, Atmophere, v)  { return (v + 14.69594877551) / 14.69594877551; }

// to Bar 
Declare_Conversion(Atmophere, Bar, v)  { return v * 1.0132499658; }
Declare_Conversion(Pascal, Bar, v)  { return v * 1.0e-5; }
Declare_Conversion(kPascal, Bar, v)  { return v * 1.0e-2; }
Declare_Conversion(psia, Bar, v)  { return v * 0.0689478018; }
Declare_Conversion(psig, Bar, v) { return (v + 14.695948775) * 0.0689478018; }

// To Pascal
Declare_Conversion(Atmophere, Pascal, v) { return v * 101324.99658; }
Declare_Conversion(Bar, Pascal, v) { return  v * 100000.0; }
Declare_Conversion(psia, Pascal, v) { return  v * 6894.757293178308 ; }
Declare_Conversion(psig, Pascal, v) { return (v + 14.695948775) * 6894.757293178308; }
Declare_Conversion(kPascal, Pascal, v) { return v * 1000; }

// To kPascal
Declare_Conversion(Atmophere, kPascal, v) { return v * 101.32499658; }
Declare_Conversion(Bar, kPascal, v) { return  v * 100.0000; }
Declare_Conversion(psia, kPascal, v) { return  v * 6.894757293178308; }
Declare_Conversion(psig, kPascal, v) { return (v + 14.695948775) * 6.894757293178308 ; }
Declare_Conversion(Pascal, kPascal, v) { return v / 1000; }

// To Absolute Pound Force per Square Inch (psia)
Declare_Conversion(Atmophere, psia, v) { return v * 14.695948775; }
Declare_Conversion(Bar, psia, v) { return v * 14.503773773; }
Declare_Conversion(Pascal, psia, v) { return v * 1.4503773773e-4; }
Declare_Conversion(kPascal, psia, v) { return v * 1.4503773773e-1; }
Declare_Conversion(psig, psia, v) { return v + 14.695948775; }

// To Pound Force per Square Inch (psig)
Declare_Conversion(psia, psig, v) { return v - 14.695948775; }
Declare_Conversion(Atmophere, psig, v) { return v * 14.695948775 - 14.695948775; }
Declare_Conversion(Bar, psig, v) { return v * 14.503773773 - 14.695948775; }
Declare_Conversion(Pascal, psig, v) { return v * 1.4503773773e-4 - 14.695948775; }
Declare_Conversion(kPascal, psig, v) { return v * 1.4503773773e-1 - 14.695948775; }
  

# endif // PRESSURE_UNIT_H



