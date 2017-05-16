# ifndef WATER_VOLUME_FACTOR_IMPL_H
# define WATER_VOLUME_FACTOR_IMPL_H

# include <correlations/water-volume-factor.H>


void BwSpiveyMNDry::precondition(const Quantity<mPascal> & p,
				 const Quantity<mPascal> & pv) const
{
  if (p < pv)
    {
      ostringstream s;
      s << "BwSpiveyMN::impl: received pressure " << p
	<< " is less than computed vapor pressure " << pv;
      ZENTHROW(WrongValueVaporPressure, s.str());
    }  
}


inline double
BwSpiveyMNDry::impl(const double & t, // comes in Celsius
		    const double & p,
		    const double & pv,
		    const double & nacl) noexcept
{
  const double & m = nacl;

  // Transformation from °C to °K
  const double tk = unit_convert<Celsius, Kelvin>(t);

  // pw: DENSITY OF BRINE WITH DISSOLVED METHANE
  // Density of methane-free brine at the temperature and pressure of
  // evaluation 
  const double pgfw =
    PwSpiveyMNGasFree::get_instance().impl(t, p, nacl);

  // Specific volume of the gas-free brine [cm³/g]
  const double vgfw = 1.0/pgfw;

  // mCH4w: SOLUBILITY OF METHANE IN BRINE [gmol NaCl/kgH2O] 
  // AT THE TEMPERATURE AND PRESSURE OF EVALUATION
  const double mch4w =
    MwSpiveyMN::get_instance().impl(t, p, pv, nacl);

  // VMCH4w: PARTIAL MOLAR VOLUME OF METHANE IN BRINE AT THE
  // TEMPERATURE AND PRESSURE OF EVALUATION
  
  // Derivatives with respect to P
  const double c3 = -8.5658e-2 + 1.31961e-2*log(tk) + 7.338/tk +
    9.05e-2/(680 - tk);
  
  const double c4 = 2.260e-1/tk + 2 * -1.227e-3 * p/tk;

  constexpr double R = 8.314467; // Universal gas constant [MPa
				 // cm³/gmol K] (McCain et al., 2011) 
            
  // Partial molar volume of methane in brine
  const double vmch4w = R*tk*(c3 + 2*m*c4);

  const Quantity<Celsius> Tstd = 15.555556; // °C (60 °F)  
  constexpr double Pstd = 0.101352; // MPa (14.7 psia)

  // Density of methane-free brine [g/cm³] at standard temperature and pressure
  const double pgfwstd =
    PwSpiveyMNGasFree::get_instance().impl(Tstd.raw(), Pstd, nacl);

  // Specific volume of methane-free brine [cm³/g]
  const double vgfwstd = 1.0/pgfwstd;

  // Bw: FORMATION VOLUME FACTOR [res cm³/ST cm³]       
  // Molecular weight MNaCl: 58.4428 g/gmol
  const double n = (1000 + m * 58.4428) * vgfw + mch4w * vmch4w; // Volume
								 // at
								 // reservoir
								 // conditions 
  const double d = (1000 + m * 58.4428) * vgfwstd; // Volume at stock
						   // tank conditions
  const double bw = n/d;

  return double(bw);
}


# endif
