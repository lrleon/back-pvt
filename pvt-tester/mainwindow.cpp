

#include "mainwindow.h"
#include "ui_mainwindow.h"


void MainWindow::set_substype_combo(const string &type_name)
{
  subtypes = Correlation::subtype_list(type_name);

  ui->corr_subtype_combo->clear();
  for (auto it = subtypes.get_it(); it.has_curr(); it.next())
    ui->corr_subtype_combo->addItem(it.get_curr().c_str());

  auto first = subtypes.get_first();
  set_corr_combo(first);
}

void MainWindow::set_corr_combo(const string &subtype_name)
{
  correlations = Correlation::list(subtype_name);

  ui->corr_combo->clear();
  for (auto it = correlations.get_it(); it.has_curr(); it.next())
    ui->corr_combo->addItem(it.get_curr()->name.c_str());

  set_correlation(correlations.get_first());
}

void MainWindow::set_tech_note(const string &note)
{
  ui->tech_note->setText(note.c_str());
}

void MainWindow::set_correlation(const Correlation * __correlation)
{  
  correlation = __correlation;
  set_tech_note(correlation->full_desc(70));

  auto frame = ui->parameters_area;

  // delete previous correlations combo
  for (auto it = pars_vals.get_it(); it.has_curr(); it.next())
    {
      auto t = it.get_curr();
      auto lyt = get<0>(t);
      auto name = get<1>(t);
      auto spin_box = get<2>(t);
      auto units_combo = get<3>(t);

      units_combo->setParent(nullptr);
      spin_box->setParent(nullptr);
      units_combo->setParent(nullptr);

      delete units_combo;
      delete spin_box;
      delete name;
      delete lyt;
    }

  pars_vals.empty();
  auto pars = correlation->get_preconditions();

  size_t i = 0;
  for (auto it = pars.get_it(); it.has_curr(); it.next())
    {
      const CorrelationPar & par = it.get_curr();
      auto par_name = QString::fromStdString(par.name);

      QHBoxLayout * lyt = new QHBoxLayout;

      QLabel * name = new QLabel(par_name + " = ", this);

      QDoubleSpinBox * spin_box = new QDoubleSpinBox(this);
      spin_box->setDecimals(10);
      spin_box->setMinimum(par.min_val.raw());
      spin_box->setMaximum(par.max_val.raw());
      spin_box->setSingleStep(par.step_size(20));

      connect(spin_box, SIGNAL(valueChanged(double)),
              this, SLOT(par_val_changed(double)));

      const string unit_symbol = par.unit.symbol;

      QComboBox * units_combo = new QComboBox(this);
      units_combo->addItem(unit_symbol.c_str());
      par.unit.sibling_units().for_each([units_combo] (auto ptr)
      {
        units_combo->addItem(ptr->symbol.c_str());
      });

      connect(units_combo, SIGNAL(activated(QString)),
              this, SLOT(par_unit_changed(QString)));

      lyt->addWidget(name);
      lyt->addWidget(spin_box);
      lyt->addWidget(units_combo);
      frame->addLayout(lyt);
      pars_vals.append(make_tuple(lyt, name, spin_box, units_combo,
                                  unit_symbol, i));
    }

  auto result_combo = ui->result_unit_combo;
  result_combo->clear();
  result_combo->addItem(correlation->unit.symbol.c_str());
  correlation->unit.sibling_units().for_each([result_combo] (auto ptr)
  {
    result_combo->addItem(ptr->symbol.c_str());
  });

  compute();
}

void MainWindow::set_exception(const string &msg)
{
  const string str = "Status: " + msg;
  ui->status->setText(str.c_str());
  ui->result->setText("Result = exception!");
}

void MainWindow::reset_status()
{
  ui->status->setText("Status:");
}

void MainWindow::show_result()
{
  double val = result.first;
  auto combo_unit =
    Unit::search_by_symbol(ui->result_unit_combo->currentText().toStdString());

  if (combo_unit != result.second)
    val = unit_convert(*result.second, result.first, *combo_unit);

  const string str = string("Result = ") + to_string(val);
  ui->result->setText(str.c_str());
  ui->result->show();
}

void MainWindow::compute()
{
  DynList<VtlQuantity> pars;
  for (auto it = pars_vals.get_it(); it.has_curr(); it.next())
    {
      auto par = it.get_curr();
      auto val = get<2>(par)->value();
      auto symbol = get<3>(par)->currentText().toStdString();
      VtlQuantity v(symbol, val);
      pars.append(v);
    }

  try
  {
    auto r = correlation->compute_and_check(pars);
    result.first = r.raw();
    result.second = &r.unit;
    show_result();
    reset_status();
  }
  catch (exception & e)
  {
    set_exception(e.what());
  }
}

MainWindow::MainWindow(QWidget *parent) :
  QMainWindow(parent),
  ui(new Ui::MainWindow)
{
  ui->setupUi(this);

  types = Correlation::type_list();
  for (auto it = types.get_it(); it.has_curr(); it.next())
    ui->corr_type_combo->addItem(it.get_curr().c_str());

  set_substype_combo(types.get_first());
}

MainWindow::~MainWindow()
{
  delete ui;
}

void MainWindow::on_corr_type_combo_activated(const QString &arg1)
{
  auto type_name = arg1.toStdString();
  set_substype_combo(type_name);
}

void MainWindow::on_corr_subtype_combo_activated(const QString &arg1)
{
  auto subtype_name = arg1.toStdString();
  set_corr_combo(subtype_name);
}

void MainWindow::on_corr_combo_activated(const QString &arg1)
{
  auto corr_name = arg1.toStdString();
  set_correlation(Correlation::search_by_name(corr_name));
}

void MainWindow::on_result_unit_combo_activated(const QString &arg1)
{
  try
  {
    show_result();
  }
  catch (exception & e)
  {
    set_exception(e.what());
  }
}

void MainWindow::par_unit_changed(const QString &arg1)
{
  QObject * owner = sender();

  auto ptr =
      pars_vals.find_ptr([owner] (auto t) { return get<3>(t) == owner; });
  if (ptr == nullptr)
    {
      set_exception("Warning: unit not found");
      return;
    }

  QComboBox * unit_combo = get<3>(*ptr);
  const auto old_unit_symbol = get<4>(*ptr);
  const auto new_unit_symbol = unit_combo->currentText().toStdString();

  auto conversion_fct =
      search_conversion_fct(old_unit_symbol, new_unit_symbol);
  if (conversion_fct == nullptr)
    {
      ostringstream s;
      s << "Not found unit conversion from " << old_unit_symbol << " to "
        << new_unit_symbol;
      throw domain_error(s.str());
    }

  QDoubleSpinBox * spin_box = get<2>(*ptr);

  double old_min = spin_box->minimum();
  double old_val = spin_box->value();
  double old_max = spin_box->maximum();

  double new_min = conversion_fct(old_min);
  double new_val = conversion_fct(old_val);
  double new_max = conversion_fct(old_max);
  if (new_max < new_min)
    swap(new_max, new_min);

  double step_size = (new_max - new_min)/20;
  spin_box->setSingleStep(step_size);
  spin_box->setRange(new_min, new_max);  
  spin_box->setValue(new_val);

  get<4>(*ptr) = new_unit_symbol;
}

void MainWindow::par_val_changed(double val)
{
  QObject * owner = sender();

  auto ptr =
      pars_vals.find_ptr([owner] (auto t) { return get<2>(t) == owner; });
  if (ptr == nullptr)
    {
      set_exception("Panic: widget not found");
      return;
    }

  compute();
}
