# ifndef MIXTURE_GAS_PSEUDOCRITICAL_PRESSURE_IMPL_H
# define MIXTURE_GAS_PSEUDOCRITICAL_PRESSURE_IMPL_H

# include <correlations/mixture-gas-pseudocritical-pressure.H>

void
PpcmKayMixingRule::precondition(const Quantity<MoleFraction> & n2,
				const Quantity<MoleFraction> & co2,
				const Quantity<MoleFraction> & h2s)
  const
{
  if (n2.raw() + co2.raw() + h2s.raw() >= 1)
    {
    ostringstream s;
    s << "PpcMKayMixingRule::impl: invalid values for nonhydrocarbon concentrations "
      << n2 << ", " << co2 << ", " << h2s;
    ZENTHROW(WrongNonhydrocarbonsConcentrations, s.str());
    }
}

inline double
PpcmKayMixingRule::impl(const double & ppchc,
			const double & n2,
			const double & co2,
			const double & h2s) noexcept
{
  const double ppcm = (1 - n2 - co2 - h2s)*ppchc + 493*n2 + 1071*co2 + 1306 * h2s;
  return ppcm;
}

void AdjustedppcmWichertAziz::precondition
(const Quantity<MoleFraction> & co2,
 const Quantity<MoleFraction> & h2s) const
{
  if (co2.raw() + h2s.raw() >= 1)
    {
    ostringstream s;
    s << "AdjustedppcMWichertAziz::impl: invalid values for nonhydrocarbon concentrations "
      << co2 << ", " << h2s;
    ZENTHROW(WrongNonhydrocarbonsConcentrations, s.str());
    }
}

inline double
AdjustedppcmWichertAziz::impl(const double & ppcm,
			      const double & tpcm,
			      const double & co2,
			      const double & h2s) noexcept
{
  const double a = co2 + h2s;
  const double b = h2s;
        
  const double e = 120*(pow(a, 0.9) - pow(a, 1.6)) + 15*(sqrt(b) - pow(b, 4));

  const double n = ppcm * (tpcm - e);
  const double d = tpcm + b * (1 - b) * e;
        
  const double adjusted_ppc = n/d;

  return adjusted_ppc;
}


# endif
