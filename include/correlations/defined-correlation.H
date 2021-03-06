
# ifndef DEFINED_CORRELATION_H
# define DEFINED_CORRELATION_H 1

# include "correlation.H"

class DefinedCorrelation
{
  struct Interval
  {
    const Correlation * correlation_ptr = nullptr;
    double start, end;
    const DynSetTree<string> * par_names = nullptr;

    double c = 0, m = 1; // set if the correlation is tuned
    bool tuned = false;

    const Unit * result_unit = nullptr;
    const Unit * tuned_unit = nullptr;

    Interval() {}

    Interval(const Correlation * ptr, double start, double end)
      : correlation_ptr(ptr), start(start), end(end)
    {
      if (correlation_ptr)
	par_names = &ptr->parameter_names();
      if (start <= end)
	return;
      ostringstream s;
      s << "start val " << start << " is greater than end val " << end;
      throw domain_error(s.str());
    }

    Interval(const Correlation * ptr, double start, double end,
	     double c, double m, const Unit & tuned_unit)
      : Interval(ptr, start, end)
    {
      this->c = c;
      this->m = m;
      this->tuned_unit = &tuned_unit;
      tuned = this->c != 0 and this->m != 1;
    }

    VtlQuantity to_result_unit(const VtlQuantity & val) const
    {
      assert(result_unit != nullptr);
      if (result_unit == &val.unit)
	return val;
      return VtlQuantity(*result_unit, val);
    }

    VtlQuantity
    compute_and_check(const DynList<VtlQuantity> & pars,
		      bool check = true) const
    {
      auto val = tuned ?
	correlation_ptr->tuned_compute_and_check(pars, c, m, *tuned_unit, check) :
	correlation_ptr->compute_and_check(pars, check);
      return to_result_unit(val);
    }

    VtlQuantity compute(const DynList<Correlation::ParByName> & pars,
			bool check = true) const
    {
      auto val = tuned ?
	correlation_ptr->tuned_compute_by_names(pars, c, m, *tuned_unit, check) :
	correlation_ptr->compute_by_names(pars, check);
      return to_result_unit(val);
    }

    VtlQuantity compute(const DynList<Correlation::NamedPar> & pars,
			bool check = true) const
    {
      auto val = tuned ?
	correlation_ptr->tuned_compute_by_names(pars, c, m, *tuned_unit, check) :
	correlation_ptr->compute_by_names(pars, check);
      return to_result_unit(val);
    }

    VtlQuantity compute(const ParList & pars, bool check = true) const
    {
      auto val = tuned ?
	correlation_ptr->tuned_compute_by_names(pars, c, m, *tuned_unit, check) :
	correlation_ptr->compute_by_names(pars, check);
      return to_result_unit(val);
    }
  };

  struct Cmp
  {
    bool operator () (const Interval & i1, const Interval & i2) const noexcept
    {
      return i1.end < i2.start;
    }
  };

  DynSetTree<Interval, Avl_Tree, Cmp> intervals;
  string main_par_name;
      // store all the parameter names for all correlations
  DynSetTree<string> par_names;

  double min_val = 0;
  double max_val = numeric_limits<double>::max();

  const Unit & unit; // unit of pivot parameter

public:

  DynList<const Correlation*> correlations() const
  {
    return intervals.maps<const Correlation*> ([] (const auto & i)
					       { return i.correlation_ptr; });
  }

  const Unit * result_unit = nullptr; // it will be the result unit of
				      // the first defined correlation

  /// Returns a list union of all parameters names. Each item is a
  /// pair with the parameter name and a list of parameter synonyms
  DynList<pair<string, DynList<string>>> parameter_list() const
  {
    using P = pair<string, DynList<string>>;
    DynMapTree<string, DynList<string>> pars;
    intervals.for_each([&pars] (const auto & i)
      {
	i.correlation_ptr->get_preconditions().
	  for_each([&pars] (const auto & par)
          {
	    auto ptr = pars.insert(P(par.name, DynList<string>()));
	    if (ptr)
	      par.get_synonyms().for_each([ptr] (const auto & p)
                {
		  ptr->second.append(p.first);
		});	    
	  });
      });
    return pars.items_ptr().maps<P>([] (auto ptr)
      {
	return P(ptr->first, ptr->second);
      });
  }

    /// Returns a list union of all parameters. Each item is a pair with
  /// the parameter name and a list of parameter synonyms
  static DynList<pair<string, DynList<string>>>
    parameter_list(const DynList<const Correlation*> & l)
  {
    using P = pair<string, DynList<string>>;
    DynMapTree<string, DynList<string>> pars;
    l.for_each([&pars] (const auto & p)
      {
	p->get_preconditions().for_each([&pars] (const auto & par)
          {
	    auto ptr = pars.insert(P(par.name, DynList<string>()));
	    if (ptr)
	      par.get_synonyms().for_each([ptr] (const auto & p)
                {
		  ptr->second.append(p.first);
		});
	  });
      });
    return pars.items_ptr().maps<P>([] (auto ptr)
      {
	return P(ptr->first, ptr->second);
      });
  }

  Interval * search_interval(double val) const noexcept
  {
    Interval i(nullptr, val, val);
    return intervals.search(i);
  }

  Interval * search_interval(const VtlQuantity & v) const
  {
    return search_interval(VtlQuantity(unit, v).raw());
  }

  const Correlation * search_correlation(double val) const 
  {
    auto interval_ptr = search_interval(val);
    if (not interval_ptr)
      return nullptr;
    return interval_ptr->correlation_ptr;
  }

  const Correlation * search_correlation(const VtlQuantity & val) const
  {
    return search_correlation(VtlQuantity(unit, val).raw());
  }

  const DynSetTree<string> & search_parameters(double val) const
  {
    auto interval_ptr = search_interval(val);
    if (not interval_ptr)
      ZENTHROW(InvalidValue, "value " + ::to_string(val) + " not found");
    return *interval_ptr->par_names;
  }

  const DynSetTree<string> & search_parameters(const VtlQuantity & val) const
  {
    return search_parameters(VtlQuantity(unit, val).raw());
  }

  bool contains(double val) const { return search_interval(val); }

  bool contains(const VtlQuantity & val) const
  {
    return search_interval(VtlQuantity(unit, val).raw());
  }

  DynList<Interval> interval_list() const { return intervals.keys(); }

  DefinedCorrelation(const string & main_par_name, const Unit & unit)
    : main_par_name(main_par_name), unit(unit) {}

private:

  void set_result_unit(Interval * interval, const Correlation * corr_ptr,
		       double start, double end)
  {
    if (interval)
      {
	if (result_unit == nullptr)
	  result_unit = &corr_ptr->unit;
	interval->result_unit = result_unit;
	corr_ptr->get_preconditions().for_each([this] (const auto & par)
          {
	    par_names.insert(par.name);
	  });
	return;
      }
    ostringstream s;
    s << "DefinedCorrelation::add_correlation: interval [" << start << ", "
      << end << "] is already inserted";
    throw domain_error(s.str());
  }

public:

  void add_tuned_correlation(const Correlation * corr_ptr,
			     double start, double end,
			     double c, double m, const Unit & tuned_unit)
  {
    auto interval =
      intervals.insert(Interval(corr_ptr, start, end, c, m, tuned_unit));
    set_result_unit(interval, corr_ptr, start, end);
  }

  void add_tuned_correlation(const Correlation * correlation_ptr,
			     const VtlQuantity & start, const VtlQuantity & end,
			     double c, double m, const Unit & tuned_unit)
  {
    add_tuned_correlation(correlation_ptr, VtlQuantity(unit, start).raw(),
			  VtlQuantity(unit, end).raw(), c, m, tuned_unit);
  }

  void add_correlation(const Correlation * corr_ptr,
		       double start, double end)
  {
    auto interval = intervals.insert(Interval(corr_ptr, start, end));
    set_result_unit(interval, corr_ptr, start, end);
  }

  void add_correlation(const Correlation * correlation_ptr,
		       const VtlQuantity & start, const VtlQuantity & end)
  {
    add_correlation(correlation_ptr, VtlQuantity(unit, start).raw(),
		    VtlQuantity(unit, end).raw());
  }

  void set_min(double min_val) { this->min_val = min_val; }

  void set_max(double max_val) { this->max_val = max_val; }

  void set_min(const VtlQuantity & min_val)
  {
    set_min(VtlQuantity(unit, min_val).raw());
  }

  void set_max(const VtlQuantity & max_val)
  {
    set_max(VtlQuantity(unit, max_val).raw());
  }

  VtlQuantity
  compute_by_names(const DynList<Correlation::ParByName> & pars,
		   bool check = true) const
  {
    auto main_par_ptr =
      pars.find_ptr([this] (auto p) { return p.first == main_par_name; });
    if (main_par_ptr == nullptr)
      {
	ostringstream s;
	s << "DefinedCorrelation: main parameter name " << main_par_name
	  << " was not found in parameter list";
	throw domain_error(s.str());
      }

    VtlQuantity main_val = { main_par_ptr->first, main_par_ptr->second };

    Interval * interval_ptr = search_interval(main_val);
    if (interval_ptr == nullptr)
      {
	ostringstream s;
	s << "DefinedCorrelation: value " << main_val
	  << " was not found in any interval";
	throw domain_error(s.str());
      }

    VtlQuantity val = interval_ptr->compute(pars, check);
    if (val < min_val)
      return VtlQuantity(val.unit, min_val);
    if (val > max_val)
      return VtlQuantity(val.unit, max_val);
    return val;
  }

  VtlQuantity compute_by_names(const DynList<Correlation::NamedPar> & pars,
			       bool check = true) const
  {
    auto main_par_ptr =
      pars.find_ptr([this] (auto t) { return get<1>(t) == main_par_name; });
    if (main_par_ptr == nullptr)
      {
	ostringstream s;
	s << "DefinedCorrelation: main parameter name " << main_par_name
	  << " was not found in parameter list";
	throw domain_error(s.str());
      }

    VtlQuantity main_val = { *get<3>(*main_par_ptr), get<2>(*main_par_ptr) };

    Interval * interval_ptr = search_interval(main_val);
    if (interval_ptr == nullptr)
      {
	ostringstream s;
	s << "DefinedCorrelation: value " << main_val
	  << " was not found in any interval";
	throw domain_error(s.str());
      }

    VtlQuantity val = interval_ptr->compute(pars, check);
    if (val < min_val)
      return VtlQuantity(val.unit, min_val);
    if (val > max_val)
      return VtlQuantity(val.unit, max_val);
    return val;
  }

  VtlQuantity compute_by_names(const ParList & pars, bool check = true) const
  {
    VtlQuantity main_val = pars.search(main_par_name);
    Interval * interval_ptr = search_interval(main_val);
    if (interval_ptr == nullptr)
      {
	ostringstream s;
	s << "DefinedCorrelation: value " << main_val
	  << " was not found in any interval";
	throw domain_error(s.str());
      }

    VtlQuantity val = interval_ptr->compute(pars, check);
    if (val < min_val)
      return VtlQuantity(val.unit, min_val);
    if (val > max_val)
      return VtlQuantity(val.unit, max_val);
    return val;
  }
};


# endif // DEFINED_CORRELATION_H
