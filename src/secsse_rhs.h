//  Copyright (c) 2021 - 2023, Thijs Janzen
//  Copyright (c) 2023, Hanno Hildenbrandt
//
//  Distributed under the Boost Software License, Version 1.0. (See
//  accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

#pragma once
#include <Rcpp.h>
#include <RcppParallel.h>
#include <type_traits>
#include <vector>

enum states {DE_0, DE_1, DM3_0, DM3_1, E_0, E_1, D_A};

namespace loglik {

  template <typename T>
  using rvector = RcppParallel::RVector<T>;

  template <typename T>
  using rmatrix = RcppParallel::RMatrix<T>;


  template <typename T>
  class vector_view_t {
  public:
    vector_view_t(T* data, size_t n) : first_(data), n_(n) {};

    size_t size() const noexcept { return n_; }
    T* begin() noexcept { return first_; }
    T* end() noexcept { return first_ + n_; }
    T& operator[](size_t i) { return *(first_ + i); }
    void advance(size_t s) noexcept { first_ += s; }

  private:
    T* first_ = nullptr;
    size_t n_ = 0;
  };

  class ode_rhs {
    const rvector<const double> m_; // extinction rates
    const rvector<const double> l_; // speciation rates
    const rvector<const double> q_; // transition rates

  public:

    // constructor
    ode_rhs(const Rcpp::NumericVector ll,
            const Rcpp::NumericVector& m,
            const Rcpp::NumericMatrix& q)
    : m_(m), q_(q), l_(ll) {
    }

    size_t size() const noexcept { return m_.size(); }

    void mergebranch(const std::vector<double>& N,
                     const std::vector<double>& M,
                     std::vector<double>& out) const {

      // substitute the code below with your own node-merging code

    }

    // this is the dx/dt calculation // true rhs that gets integrated
    // along the branches
    void operator()(const std::vector<double>& x,
                    std::vector<double>& dxdt,
                    const double /* t */) const
    {

      // substitute with your own code below:

     // vector x is:
     // [
     // DE_0, DE_1
     // DM3_0, DM3_1
     // E_0, E_1
     // D_A
     // ]

     dxdt[DE_0] = -(lambda_c_0 + mu_0 + q_01) * x[DE_0] + 2 * lambda_0_c * x[DE_0] * x[E_0] + q_01 * x[E_1];

      // etc












    }
  };
} // namespace secsse
