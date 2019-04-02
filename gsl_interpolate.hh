#pragma once

#include <vector>
#include <gsl/gsl_interp.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>

template <typename C>
class gsl_lin_interp_pp {
  C xdata;
  C ydata;
  double* x;
  double* y;
  size_t n;

  gsl_interp *interpolation;
  gsl_interp_accel * accelerator;

  public:
  gsl_lin_interp_pp(C xd, C yd):
    xdata(std::move(xd)), ydata(std::move(yd)),
    x(xdata.data()), y(ydata.data()),
    n(xdata.size())
  {
    if (xdata.size() != ydata.size()) throw;
    interpolation = gsl_interp_alloc (gsl_interp_linear,n);
    gsl_interp_init(interpolation, x, y, n);
    accelerator =  gsl_interp_accel_alloc();
  }

  ~gsl_lin_interp_pp() {
    gsl_interp_free(interpolation);
    gsl_interp_accel_free(accelerator);
  }

  double operator()(double xp) {
    return gsl_interp_eval(interpolation, x, y, xp, accelerator);
  }
};

template <typename C>
gsl_lin_interp_pp<C> make_gsl_lin_interp(const C & xd, C const & yd) {
  return gsl_lin_interp_pp<C>(xd, yd);
}

gsl_lin_interp_pp<std::vector<double>> 
make_gsl_lin_interp(std::initializer_list<double> xd, std::initializer_list<double> yd) {
  return make_gsl_lin_interp(std::vector<double>(xd), std::vector<double>(yd));
}

gsl_lin_interp_pp<std::vector<double>>
make_gsl_lin_interp(const std::string & filename, const size_t col1, const size_t col2) {
  std::ifstream file(filename);
  std::string line;
  std::vector<double> vx;
  std::vector<double> vy;
  while (std::getline(file, line)) {
    line = line.substr(0, line.find("#"));
    std::stringstream ss(line);
    std::vector<double> vline(std::istream_iterator<double>(ss), {});
    if (vline.empty()) continue;
    auto check_col = [&](const size_t col) {
      if (col > vline.size()) {
        throw std::out_of_range("Asked for colum "+std::to_string(col)+" but "+filename+" only has "+std::to_string(vline.size())+" columns.");
      }
    };
    check_col(col1);
    check_col(col2);
    vx.push_back(vline[col1-1]);
    vy.push_back(vline[col2-1]);
  }
  if (vx.size() == 0) {
    throw std::runtime_error("No usable data in file " + filename);
  }
  return make_gsl_lin_interp(std::move(vx), std::move(vy));
}

