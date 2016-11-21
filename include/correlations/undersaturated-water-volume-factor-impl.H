# ifndef UNDERSATURATED_WATER_VOLUME_FACTOR_IMPL_H
# define UNDERSATURATED_WATER_VOLUME_FACTOR_IMPL_H

# include <correlations/undersaturated-water-volume-factor.H>
# include <correlations/spivey-dry.H>

inline Quantity<Sm3_Sm3>
BwaSpiveyMN::impl(const Quantity<Celsius> & t,
		  const Quantity<mPascal> & p,
		  const Quantity<Molality_NaCl> & nacl_concentration) const
{
  // Transformation from °C to °K
  const double tk = Quantity<Kelvin>(t).raw();

  // Vapor pressure of pure water, calculated from the IAWPS-95 formulation
  const double pv =
    PvSpiveyMN::get_instance().impl(tk).raw();

  if (p < pv)
    {
      ostringstream s;
      s << "BwaSpiveyMN::impl: received pressure " << p
	<< " is less than computed vapor pressure " << pv << " "
	<< p.unit.name;
      throw domain_error(s.str());
    }

  const double bw =
    BwSpiveyMNDry::get_instance().impl(t, p, pv, nacl_concentration).raw();

  return Quantity<Sm3_Sm3>(bw);
}

inline Quantity<RB_STB>
BwaMcCain::impl(const Quantity<psia> & p,
                const Quantity<psia> & pb,
		const Quantity<RB_STB> & bwbp,
		const Quantity<psia_1> & cw) const
{

  const double bw = bwbp.raw() * exp(cw * (pb.raw() - p.raw()));

  return Quantity<RB_STB>(bw);
}






# endif // UNDERSATURATED_WATER_VOLUME_FACTOR_IMPL_H