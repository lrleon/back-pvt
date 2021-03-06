
# include <tuple>
# include <iostream>

# include <tclap/CmdLine.h>

# include <metadata/pvt-analyse.H>

using namespace std;
using namespace TCLAP;

int main(int argc, char *argv[])
{
  CmdLine cmd = { argv[0], ' ', "0" };

  ValueArg<string> json_file = { "f", "json-file", "json file", false,
				 "data.json", "json file name", cmd };

  cmd.parse(argc, argv);

  ifstream input(json_file.getValue());
  if (not input)
    {
      cout << "cannot open " << json_file.getValue() << endl;
      return 0;
    }

  PvtAnalyzer pvt(input);

  auto bob_vals = pvt.get_data().search_variable("Below Pb", "bob");
  if (not get<0>(bob_vals))
    {
      cout << "The data set does not contain the bob variable" << endl;
      abort();
    }

  cout << pvt.get_data().full_desc() << endl;

  auto boa_samples = pvt.get_data().values("Above Pb", "boa").first;
  cout << Rvector("boa", boa_samples) << endl;

  auto boa_corr = Correlation::array().filter([] (auto p)
    { return p->target_name() == "Boa"; });
  cout << "Boa correlations" << endl;
  boa_corr.for_each([] (auto p) { cout << "  " << p->call_string() << endl; });
  
  auto boa_apliable = pvt.boa_correlations();
  if (boa_apliable.is_empty())
    {
      cout << "There is not any boa correlation" << endl;
      return 0;
    }  

  auto boa_valid = pvt.boa_valid_correlations();
  if (boa_valid.is_empty())
    {
      cout << "Valid correlation is is empty" << endl;
      return 0;
    }

  auto c = boa_valid.get_first();
  cout << c->call_string() << endl << endl;

  auto p = pvt.get_data().values("Above Pb", "p").first;

  cout << Rvector("p", p) << endl;

  cout << "Boa valid correlations:" << endl;
  boa_valid.for_each([] (auto p) { cout << " " << p->call_string() << endl; });
  cout << endl;

  auto boa_fits = sort(pvt.boa_correlations_lfits(),
		       [] (const auto & t1, const auto &t2)
		       {
			 return get<4>(t1).sumsq < get<4>(t2).sumsq;
		       });

  cout << get<4>(boa_fits.get_first()).to_string() << endl
       << pvt.boa_lfit(get<2>(boa_fits.get_first())).to_string() << endl;

  assert(get<4>(boa_fits.get_first()) ==
	 pvt.boa_lfit(get<2>(boa_fits.get_first())));

  boa_fits.for_each([&pvt] (auto c) { cout << pvt.to_string(c) << endl; });

  boa_valid.for_each([&pvt] (auto t)
    {
      cout << Rvector(t->name, pvt.get_data().compute(1, t)) << endl;
    });

  auto boa_lfits_list = pvt.boa_lfits_list(boa_fits);

  cout << pvt.to_R("tuned.", boa_lfits_list) << endl;

  auto p_below = pvt.get_data().values("Below Pb", "p").first;
  auto best_bob_correlation = &BobStanding::get_instance();
  auto below_fit = pvt.bob_lfit(best_bob_correlation);

  auto p_above = pvt.get_data().values("Above Pb", "p").first;
  auto best_boa_correlation = get<2>(boa_fits.get_first());
  auto above_fit = pvt.boa_lfit(best_boa_correlation);

  DefinedCorrelation defcorr("p", psia::get_instance());

  defcorr.add_tuned_correlation(best_bob_correlation,
				p_below.get_first(), p_below.get_last(),
				below_fit.c, below_fit.m);
  defcorr.add_tuned_correlation(best_boa_correlation,
				p_above.get_first(), p_above.get_last(),
				above_fit.c, above_fit.m);

  auto l = pvt.get_data().generate_samples_by_name(best_boa_correlation, 0,
						   "p", 0, 15000);
  l.for_each([] (const auto & row)
    {
      row.for_each([] (auto p) { cout << p.first << " = " << p.second << " "; });
      cout << endl;
    });

  auto lc = pvt.get_data().generate_samples_by_name(defcorr, "p", {0, 1});
  cout << endl
       << "All samples:" << endl;
  lc.for_each([] (const auto & p)
    {
      cout << p.first->call_string() << endl;
      p.second.for_each([] (auto l)
         {
	   l.for_each([] (auto p)
		      { cout << "  " << p.first << " = " << p.second << " "; });
	   cout << endl;
	 });
      cout << endl;
    });

  cout << endl
       << endl
       << "Parameter list:" << endl;

  auto samples = pvt.get_data().samples_by_name(defcorr, "p", { 0, 1 });

  samples.for_each([] (const auto & l)
      {
	l.for_each([] (auto p)
		   {
		     cout << p.first << " = " << p.second << " ";
		   });
	cout << endl;
      });
  cout << endl;

  auto values = samples.maps<double>([&defcorr] (const auto & pars)
    {
      return defcorr.compute_by_names(pars);
    });

  auto pressure = pvt.get_data().values(DynList<size_t>({0, 1}), "p").first;
  auto bo = to_dynlist(pvt.get_data().values(0, "bob").first);
  bo.append(to_dynlist(pvt.get_data().values(1, "boa").first));

  cout << Rvector("P", pressure) << endl
       << Rvector("Bo", values) << endl
       << Rvector("Bo.Lab", bo) << endl;

}

