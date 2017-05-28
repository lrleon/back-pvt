
# include <fstream>

# include <ah-string-utils.H>
# include <ahSort.H>
# include <ah-dispatcher.H>
# include <tclap/CmdLine.h>

# include <correlations/pvt-correlations.H>
# include <metadata/z-calibrate.H>

using namespace TCLAP;


struct PZArg
{
  const Unit * tunit_ptr = nullptr;
  const Unit * punit_ptr = nullptr;

  double t = 0;
  DynList<double> p;
  DynList<double> z;

  PZArg() {}

  static void
  read_and_validate_unit(const PhysicalQuantity & pq,
			 istringstream & iss, const Unit *& unit_ptr)
  {
    string data;
    if (not (iss >> data))
      ZENTHROW(CommandLineError, "cannot read unit for " + pq.name);

    unit_ptr = Unit::search(data);
    if (unit_ptr == nullptr)
      ZENTHROW(CommandLineError, "unit " + data + " not found");

    if (&unit_ptr->physical_quantity != &pq)
      ZENTHROW(CommandLineError, "unit " + data + " is not for " + pq.name);
  }

  static double read_double(istringstream & iss)
  {
    string data;
    iss >> data;
    if (not is_double(data))
      ZENTHROW(CommandLineError, "read value " + data + " is not a double");
    return atof(data);
  }

  PZArg & operator = (const string & str)
  {
    string data;
    istringstream iss(str);

    t = read_double(iss);
    read_and_validate_unit(Temperature::get_instance(), iss, tunit_ptr);
    read_and_validate_unit(Pressure::get_instance(), iss, punit_ptr);

    DynList<double> l;
    size_t n = 0;    
    for (; iss.good(); ++n)
      l.append(read_double(iss));

    if ((n % 2) != 0)
      ZENTHROW(CommandLineError, "Number of values " + to_string(n) +
	       " is not even");

    for (size_t i = 0; i < n/2; ++i)
      p.append(l.remove_first());

    for (size_t i = 0; i < n/2; ++i)
      z.append(l.remove_first());

    if (not is_sorted(p))
      {
	p = p.rev();
	z = z.rev();
      }
      
    if (not is_sorted(p) and not is_inversely_sorted(p))
      ZENTHROW(CommandLineError, "pressure values are not monotone");

    mutable_unit_convert(*punit_ptr, p, psia::get_instance());

    return *this;
  }

  friend ostream & operator << (ostream & out, const PZArg & a)
  {
    out << "t = " << a.t << " " << a.tunit_ptr->name << endl
	<< "pressure unit = " << a.punit_ptr->name << endl
	<< "p =";
    a.p.for_each([&out] (auto v) { out << " " << v; });
    out << endl
	<< "z =";
    a.z.for_each([&out] (auto v) { out << " " << v; });
    return out;
  }
};

static const DynSetTree<string> valid = { "yg", "co2", "n2", "h2s" };

struct ArgUnit
{
  string name;
  string unit_name;

  ArgUnit & operator = (const string & str)
  {
    istringstream iss(str);
    if (not (iss >> name >> unit_name))
      ZENTHROW(CommandLineError, str + " is not a pair par-name unit");

    if (not valid.contains(name))
      ZENTHROW(CommandLineError, name + " is an invalid parameter name");

    return *this;
  }

  ArgUnit() {}

  friend ostream& operator << (ostream &os, const ArgUnit & a) 
  {
    return os << a.name << " " << a.unit_name;
  }
};

struct PlotNumbers
{
  size_t n = 0;
  DynList<size_t> numbers;

  PlotNumbers() {}

  PlotNumbers & operator = (const string & str)
  {
    string data;
    istringstream iss(str);
    for (; iss >> data; ++n)
      if (not is_size_t(data))
	ZENTHROW(CommandLineError, data + " is not a unsigned integer");
      else
	numbers.append(atol(data));

    return *this;
  }
};

namespace TCLAP
{
  template<> struct ArgTraits<ArgUnit> { typedef StringLike ValueCategory; };
  template <> struct ArgTraits<PZArg> { typedef StringLike ValueCategory; };
  template <> struct ArgTraits<PlotNumbers> { typedef StringLike ValueCategory; };
}

CmdLine cmd = { "ztuner", ' ', "0" };		 

# define Declare_Arg(NAME, v)						\
  ValueArg<double> NAME##_arg = { "", #NAME, #NAME, false, v, #NAME, cmd }; \
  const Unit * NAME##_unit = nullptr;

Declare_Arg(yg, 0.6);
Declare_Arg(co2, 0);
Declare_Arg(n2, 0);
Declare_Arg(h2s, 0);

MultiArg<PZArg> zvalues = { "", "z", "z", false,
			    "t tunit punit p-list z-list", cmd };

ValueArg<string> fname = { "f", "file", "file name", false, "",
			   "file name", cmd };

SwitchArg save = { "s", "save", "save json", cmd };

SwitchArg print = { "p", "print", "print data", cmd };

SwitchArg eol = { "n", "eol", "print end of line", cmd };

MultiArg<ArgUnit> unit = { "", "unit", "change unit of input data", false,
			   "unit \"par-name unit\"", cmd };

vector<string> sort_types = { "sumsq", "c", "m", "num" };
ValuesConstraint<string> allowed_sort_types = sort_types;
ValueArg<string> sort = { "", "sort", "sort type", false,
			  "num", &allowed_sort_types, cmd };

vector<string> output_types = { "R", "csv", "mat" };
ValuesConstraint<string> allowed_output_types = output_types;
ValueArg<string> output = { "", "output", "output type", false,
			    "mat", &allowed_output_types, cmd };

SwitchArg solve = { "S", "solve", "solve z", cmd };

SwitchArg check = { "c", "check", "check z application ranges", cmd };

SwitchArg exceptions = { "e", "exceptions", "prints exceptions", cmd };

ValueArg<PlotNumbers> plot = { "P", "plot", "plot", false, PlotNumbers(),
			       "plot", cmd };

// Checks whether the parameter par_name has a change of
// unity. ref_unit is the default unit of the parameter. If there was
// no change specification for par_name, then returns ref_unit
const Unit * test_par_unit_change(const string & par_name,
				  const Unit & unit_ref)
{
  if (not valid.contains(par_name))
    {
      cout << "for option --unit " << par_name << ": invalid parameter name"
	   << endl;
      abort();
    }

  const Unit * unit_ptr = &unit_ref;
  auto & pq = unit_ref.physical_quantity;
  for (const auto & par : unit.getValue()) // traverse list of changes
    if (par.name == par_name)
      {
	unit_ptr = Unit::search_by_name(par.unit_name);
	if (unit_ptr == nullptr)
	  {
	    cout << "In unit change for " << par_name << ": unit name "
		 << par.unit_name << " not found" << endl;
	    abort();
	  }

	if (&pq != &unit_ptr->physical_quantity)
	  {
	    cout << "For " << par_name << " unit: physical quantity "
		 << pq.name << " is invalid" << endl;
	    abort();
	  }
	return unit_ptr;
      }
  return unit_ptr;
}

# define Set_Par(NAME, UNIT)						\
  if (NAME##_arg.isSet())							\
    {									\
      NAME##_unit = test_par_unit_change(#NAME, UNIT::get_instance());	\
      data_ptr->NAME = VtlQuantity(*NAME##_unit, NAME##_arg.getValue()); \
    }

void process_input(unique_ptr<Ztuner> & data_ptr)
{
  if (data_ptr == nullptr)
    data_ptr = unique_ptr<Ztuner>
      (new Ztuner(VtlQuantity(*yg_unit, yg_arg.getValue()),
		  VtlQuantity(*n2_unit, n2_arg.getValue()),
		  VtlQuantity(*co2_unit, co2_arg.getValue()),
		  VtlQuantity(*h2s_unit, h2s_arg.getValue())));

  Set_Par(yg, Sgg);
  Set_Par(n2, MolePercent);
  Set_Par(co2, MolePercent);
  Set_Par(h2s, MolePercent);

  for (auto & z : zvalues.getValue())
    data_ptr->add_z(VtlQuantity(*z.tunit_ptr, z.t), move(z.p), move(z.z));
}

void terminate_app()
{
  if (eol.getValue())
    cout << endl;
  exit(0);
}

unique_ptr<Ztuner> data;

void process_print()
{
  if (not print.isSet())
    return;

  cout << *data << endl;
  terminate_app();
}

# define Define_Cmp(NAME)						\
  static auto cmp_##NAME = [] (const Ztuner::Zcomb & z1, const Ztuner::Zcomb & z2) \
    {									\
      return Ztuner::NAME(z1) < Ztuner::NAME(z2);			\
    }

void process_plot()
{
  assert(plot.isSet());
  assert(not data->zcomb_list.is_empty());
  assert(not plot.getValue().numbers.is_empty());

  const PlotNumbers & numbers = plot.getValue();
  cout << numbers.numbers.size() << endl;
  numbers.numbers.for_each([] (auto i) { cout << i << " "; }); cout << endl;
  if (not numbers.numbers.all([n = data->zcomb_list.size()] (auto i)
			      { return i < n; }))
    ZENTHROW(CommandLineError, "Invalid number in plot list");

  DynList<string> header = data->basic_header();
  DynList<DynList<double>> cols = transpose(data->vals());

  header.for_each([] (auto & s) { cout << s << "  "; }); cout << endl;

  cout << "Not yet implemented" << endl;
}

void process_solve()
{
  Define_Cmp(sumsq);
  Define_Cmp(c);
  Define_Cmp(m);
  Define_Cmp(num);
  static DynMapTree<string, bool (*)(const Ztuner::Zcomb&, const Ztuner::Zcomb&)>
    cmp = { {"sumsq", cmp_sumsq}, {"c", cmp_c}, {"m", cmp_m}, {"num", cmp_num} };
  static auto format_mat = [] (const DynList<Ztuner::Zcomb> & l)
    {
      return to_string(format_string(Ztuner::to_dynlist(l)));
    };
  static auto format_csv = [] (const DynList<Ztuner::Zcomb> & l)
    {
      return to_string(format_string_csv(Ztuner::to_dynlist(l)));
    };
  static auto format_R = [] (const DynList<Ztuner::Zcomb> &) -> string
    {
      return "R format is incompatible for solve option";
    };  
  static AHDispatcher<string, string (*)(const DynList<Ztuner::Zcomb>&)>
    dispatcher("mat", format_mat, "csv", format_csv, "R", format_R);

  if (not solve.isSet())
    return;

  auto l = data->solve(check.getValue());
		       
  if (exceptions.getValue())
    data->exception_list.for_each([] (auto & s) { cout << s << endl; });
  else if (plot.isSet())
    process_plot();
  else
    cout << dispatcher.run(output.getValue(),
			   Aleph::sort(l, cmp.find(::sort.getValue())));
  terminate_app();
}

int main(int argc, char *argv[])
{
  cmd.parse(argc, argv);

  if (fname.isSet())
    {
      const string & file_name = fname.getValue();
      if (not exists_file(file_name) and not save.isSet())
	ZENTHROW(CommandLineError, "file with name " + file_name + " not found");
      else if (not save.isSet())
	data = unique_ptr<Ztuner>(new Ztuner(ifstream(file_name)));
    }

  process_input(data);
  if (save.getValue())
    {
      if (not fname.isSet())
	ZENTHROW(CommandLineError,
		 "for save option file name has not been set");
      ofstream out(fname.getValue());
      out << data->to_json().dump(2) << endl;
    }

  process_print();
  process_solve();
  
}