# ifndef MIXTURE_GAS_PSEUDOCRITICAL_TEMP_IMPL_H
# define MIXTURE_GAS_PSEUDOCRITICAL_TEMP_IMPL_H

# include <correlations/mixture-gas-pseudocritical-temp.H>

void TpcmKayMixingRule::precondition
(const Quantity<MoleFraction> & n2,
 const Quantity<MoleFraction> & co2,
 const Quantity<MoleFraction> & h2s) const
{
  if (n2.raw() + co2.raw() +
      h2s.raw() >= 1)
    {
    ostringstream s;
    s << "In TpcmKayMixingRule::impl: invalid values for nonhydrocarbon concentrations "
      << n2 << ", " << co2 << ", "
      << h2s;
    ZENTHROW(WrongNonhydrocarbonsConcentrations, s.str());
    }
}

inline double
TpcmKayMixingRule::impl(const double & tpchc,
			const double & n2,
			const double & co2,
			const double & h2s) noexcept
{
  const double tpcm = (1 - n2 - co2 - h2s)*tpchc + 227*n2 + 548*co2 + 672*h2s;

  return tpcm;
}

void AdjustedtpcmWichertAziz::precondition
(const Quantity<MoleFraction> & co2,
 const Quantity<MoleFraction> & h2s) const
{
  if (co2.raw() + h2s.raw() >= 1)
    {
    ostringstream s;
    s << "In AdjustedtpcmWichertAziz::impl: invalid values for nonhydrocarbon concentrations "
      << co2 << ", " << h2s;
    ZENTHROW(WrongNonhydrocarbonsConcentrations, s.str());
    }
}

inline double
AdjustedtpcmWichertAziz::impl(const double & tpcm,
			      const double & co2,
			      const double & h2s) noexcept
{
  const double a = co2 + h2s;
  const double b = h2s;
  const double b4 = b*b*b*b;

  const double e = 120*(pow(a,  0.9) - pow(a, 1.6)) + 15*(sqrt(b) - b4);

  const double ret = tpcm - e;

  return ret;
}


#endif
