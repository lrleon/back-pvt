# ifndef COMPOUND_CORR_H
# define COMPOUND_CORR_H

# include <pvt-units.H>
# include "correlation.H"

template <class Corr1, class Corr2, template <class, class> class Pars>
class CompoundCorrelation : public Pars<Corr1, Corr2>
{
  double min_val = 0;
  double max_val = numeric_limits<double>::max();

  double c1 = 0, c2 = 0, m1 = 1, m2 = 1;

  const Quantity<psia> pivot;

public:

  DynList<string> names() const
  {
    DynList<string> ret = this->corr1_me()->names.keys();
    ret.append(this->corr2_me()->names.keys());
    return ret;
  }

  auto search_correlation(const Quantity<psia> & p) const
  {
    return static_cast<const Correlation*>(p < pivot ? this->corr1_me() :
					   this->corr2_me());
  }

  bool has_parameter(const Quantity<psia> & p, const string & name) const
  {
    return p <= pivot ? this->corr1_me().names.has(name) :
      this->corr2_me().names.has(name);
  }

  const Unit * result_unit = nullptr; // it will be the result unit of
				      // the first defined correlation

private:

  static double no_convert(double val) { return val; }

  Unit_Convert_Fct_Ptr convert_corr2 = no_convert;

public:

  CompoundCorrelation(Corr1 * ptr1, Corr2 * ptr2, const Quantity<psia> & pivot)
    : Pars<Corr1, Corr2>(ptr1, ptr2), pivot(pivot),
    result_unit(&ptr1->unit)
  {
    if (&ptr1->unit != &ptr2->unit)
      {
	convert_corr2 = search_conversion(ptr2->unit, *result_unit);
	if (convert_corr2 == nullptr)
	  ZENTHROW(UnitConversionNotFound, "Conversion from " + ptr2->unit.name +
		   " to " + result_unit->name + " not found");
      }
  }

  CompoundCorrelation(Corr1 * ptr1, Corr2 * ptr2, const Quantity<psia> & pivot,
		      double c1, double m1, double c2, double m2)
    : CompoundCorrelation(ptr1, ptr2, pivot)
  {
    this->c1 = c1;
    this->m1 = m1;
    this->c2 = c2;
    this->m2 = m2;
  }

  void set_min(double min_val) { this->min_val = min_val; }

  void set_max(double max_val) { this->max_val = max_val; }

  void set_min(const VtlQuantity & min_val)
  {
    set_min(VtlQuantity(psia::get_instance(), min_val).raw());
  }

  void set_max(const VtlQuantity & max_val)
  {
    set_max(VtlQuantity(psia::get_instance(), max_val).raw());
  }

# define DRY_CALL_CORR()						\
  VtlQuantity val;							\
  if (p_val <= pivot)							\
    val = c1 + m1*this->corr1_me()->compute_by_names(pars, check);	\
  else									\
    {									\
      val = c2 + m2*this->corr2_me()->compute_by_names(pars, check);	\
      val = convert_corr2(val.raw());					\
    }									\
  if (val < min_val)							\
    return VtlQuantity(val.unit, min_val);				\
  if (val > max_val)							\
    return VtlQuantity(val.unit, max_val);				\
  return val

# define DRY_COMPUTE(TYPE)						\
  VtlQuantity								\
  compute_by_names(const TYPE & pars, bool check = true) const		\
  {									\
    auto p_par_ptr = pars.find_ptr([] (auto p) { return p.first == "p"; }); \
    if (p_par_ptr == nullptr)						\
      ZENTHROW(ParameterNameNotFound, "parameter p not found");		\
									\
    VtlQuantity p_val = { p_par_ptr->first, p_par_ptr->second };	\
									\
    DRY_CALL_CORR();							\
  }

  DRY_COMPUTE(DynList<Correlation::ParByName>);
  DRY_COMPUTE(DynList<Correlation::NamedPar>);

  VtlQuantity compute_by_names(const ParList & pars, bool check = true) const
  {
    VtlQuantity p_val = pars.search("p");
    DRY_CALL_CORR();
  }

# undef DRY_CALL_CORR
# undef DRY_COMPUTE
  
  VtlQuantity compute(const Quantity<psia> & p, bool check = true)
  {
    this->set_p(p);
    VtlQuantity val(*result_unit);

    if (p <= pivot)
      val = c1 + m1*this->corr1_me()->compute(check);
    else
      val = c2 + m2*this->corr2_me()->compute(check);
    if (val < min_val)
      return VtlQuantity(val.unit, min_val);
    if (val > max_val)
      return VtlQuantity(val.unit, max_val);
    return val;
  }
};


# endif 
