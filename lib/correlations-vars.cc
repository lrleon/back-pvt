
# include <ahSort.H>
# include <ah-stl-utils.H>
# include <correlations/pvt-correlations.H>

# include <json.hpp>

# include <pvt-units.H>

# include "gitversion.H"

using namespace Aleph;
using json = nlohmann::json;
  
size_t Correlation::counter = 0;
DynMapTree<string, const Correlation *> Correlation::tbl;
Array<const Correlation *> Correlation::correlations_tbl;

CorrelationInstantiater __correlations; 

static json to_json(const CorrelationPar & p) 
{
  json j;
  j["name"] = p.name;
  j["desc"] = p.description;
  j["unit"] = p.unit.name;
  j["unit_name"] = p.unit.name;
  j["physicalq"] = p.unit.physical_quantity.name;
  j["minv"] = p.min_val.get_value();
  j["minvu"] = p.min_val.unit.name;
  j["maxv"] = p.max_val.get_value();
  j["maxvu"] = p.min_val.unit.name;
  j["minauthor"] = p.min_from_author;
  j["maxauthor"] = p.max_from_author;
  j["latex"] = p.latex_symbol;
  return j;
}

static json to_json(const Correlation & c) 
{
  json j;
  j["maxv"] = c.max_val;
  j["minv"] = c.min_val;
  j["min_from_author"] = c.min_from_author;
  j["max_from_author"] = c.max_from_author;
  j["unit"] = c.unit.name;
  j["refs"] = to_vector(c.refs.maps<string>([] (const auto & r)
    { return r->to_string(); }));
  j["notes"] = to_vector(c.notes);
  j["db"] = to_vector(c.db);
  j["title"] = c.title;
  j["author"] = c.author;
  j["latex"] = c.latex_symbol;
  j["subtype"] = c.subtype_name;
  j["type"] = c.type_name;
  j["name"] = c.name;
  j["hidden"] = c.hidden;
  j["hidden_blackoil_grid"] = c.hidden_blackoil_grid;
  j["hidden_wetgas_grid"] = c.hidden_wetgas_grid;
  j["hidden_drygas_grid"] = c.hidden_drygas_grid;
  j["hidden_calc"] = c.hidden_calc;
  j["id"] = c.id;

  auto jpars = c.get_preconditions().maps<json>([] (const auto & par)
					       { return ::to_json(par); });
  auto pars = to_vector(jpars);
  j["pars"] = pars;

  return j;
}

string Correlation::to_json() const
{
  return ::to_json(*this).dump(2);
}

static json to_json_concise(const Correlation & c) 
{
  json j;
  j["author"] = c.author;
  j["latex"] = c.latex_symbol;
  j["name"] = c.name;
  j["hidden"] = c.hidden;
  j["hidden_blackoil_grid"] = c.hidden_blackoil_grid;
  j["hidden_wetgas_grid"] = c.hidden_wetgas_grid;
  j["hidden_drygas_grid"] = c.hidden_drygas_grid;
  j["hidden_calc"] = c.hidden_calc;
  j["id"] = c.id;

  return j;
}

string Correlation::to_json(const string & subtype_name)
{
  DynList<json> l;
  for (auto it = correlations_tbl.get_it(); it.has_curr(); it.next())
    {
      auto corr_ptr = it.get_curr();
      if (corr_ptr->subtype_name == subtype_name)
	l.append(to_json_concise(*corr_ptr));
    }

  if (l.is_empty())
    throw domain_error("Invalid subtype " + subtype_name);
  
  json j;
  j["Correlations"] = to_vector(l);
  return j.dump(2);
}

string Correlation::json_of_all_correlations()
{
  DynMapTree<string, DynMapTree<string, DynList<const Correlation * const>>>
    tree;

  // insert all correlations in the tree
  for (auto it = tbl.get_it(); it.has_curr(); it.next())
    {
      auto corr_ptr = it.get_curr().second;
      auto & subtree = tree[corr_ptr->type_name];
      // TODO: si subtype_name está marcado hidden ==> ignorar
      auto & l = subtree[corr_ptr->subtype_name];
      l.append(corr_ptr);
    }

  DynList<json> jl;
  for (auto it = tree.get_it(); it.has_curr(); it.next())
    {
      auto p = it.get_curr();
      auto type_name = p.first;

      json j;
      j["relation"] = type_name;
      vector<json> properties;

      auto subtree = p.second;
      for (auto i = subtree.get_it(); i.has_curr(); i.next())
	{
	  auto p = i.get_curr();
	  auto l = p.second;

	  json physical_property;
	  physical_property["name"] = p.first;

	  vector<json> corrs;

	  for (auto it = l.get_it(); it.has_curr(); it.next())
	    {
	      auto corr_ptr = it.get_curr();
	      //corrs.push_back(::to_json(*corr_ptr));
	      corrs.push_back(::to_json_concise(*corr_ptr));
	    }

	  physical_property["relations"] = corrs;
	  properties.push_back(physical_property);
	}

      j["nproperties"] = properties.size();
      j["Physicalq"] = properties;
      jl.append(j);
    }

  auto properties = to_vector(jl);
  json j;
  j["Physical_properties"] = properties;
  j["version"] = GITVERSION;

  return j.dump(2);
}

