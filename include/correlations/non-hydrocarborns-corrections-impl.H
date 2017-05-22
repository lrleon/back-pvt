# ifndef NON_HYDROCARBORNS_CORRECTIONS_IMPL_H
# define NON_HYDROCARBORNS_CORRECTIONS_IMPL_H

# include <correlations/non-hydrocarborns-corrections.H>

void
YghcWichertAziz::precondition(const Quantity<MoleFraction> & n2,
			      const Quantity<MoleFraction> & co2,
			      const Quantity<MoleFraction> & h2s) const
{
  if (n2.raw() + co2.raw() + h2s.raw() >= 1)
    {
      ostringstream s;
      s << "invalid values for nonhydrocarbon concentrations "
	<< n2 << "+ " << co2 << " + " << h2s << " >= " << 1 << "";
      ZENTHROW(WrongNonhydrocarbonsConcentrations, s.str());
    }
}

inline double
YghcWichertAziz::impl(const double & yg,
		      const double & n2,
		      const double & co2,
		      const double & h2s) noexcept
{
  const double yghc = (yg - 0.967*n2 - 1.52*co2 - 1.18*h2s) / (1 - n2 - co2 - h2s);

  return yghc;
}


# endif // NON_HYDROCARBORNS_CORRECTIONS_IMPL_H
