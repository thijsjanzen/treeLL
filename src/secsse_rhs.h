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

class ode_tree {

  enum states {DE_0, DE_1, DM3_0, DM3_1, E_0, E_1, DA_3, num_states};

  const rvector<const double> lc_; // cladogenesis rates
  const rvector<const double> m_; //  extinction rates
  const rvector<const double> g_; //  colonization rates
  const rvector<const double> la_; // anagenesis rates
  const rvector<const double> q_; // transition rates
  const double p_;

  // mu_0, mu_1
  // lambda_c_0, lambda_c_1
  // lambda_a_0, lambda_a_1
  // gamma_0, gamma_1
  // p
  // q_01, q_10

public:

  // constructor
  ode_tree(const Rcpp::NumericVector& lc,
          const Rcpp::NumericVector& la,
          const Rcpp::NumericVector& m,
          const Rcpp::NumericVector& g,
          const Rcpp::NumericVector& q,
          const double p)
    : lc_(lc), m_(m), g_(g), la_(la), q_(q), p_(p) {
  }

  size_t size() const noexcept { return num_states; }

  void mergebranch(const std::vector<double>& N,
                   const std::vector<double>& M,
                   std::vector<double>& out) const {

    // substitute the code below with your own node-merging code
    size_t s = num_states;
    out.resize(s);

    out[DE_0]  = lc_[0] * N[DE_0] * M[DE_0];
    out[DE_1]  = lc_[1] * N[DE_1] * M[DE_1];

    out[DM3_0] = N[DM3_0];
    out[DM3_1] = N[DM3_1];

    out[E_0]   = N[E_0];
    out[E_1]   = N[E_1];
    out[DA_3]  = N[DA_3];
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
    // DA3
    // ]

    dxdt[DE_0]  = -(lc_[0] + m_[0] + q_[0]) * x[DE_0] + 2 * lc_[0] * x[DE_0] * x[E_0] + q_[0] * x[DE_1];

    dxdt[DE_1]  = -(lc_[1] + m_[1] + q_[1]) * x[DE_1] + 2 * lc_[1] * x[DE_1] * x[E_1] + q_[1] * x[DE_0];

    dxdt[DM3_0] = -(la_[0] + lc_[0] + m_[0] + g_[1] + q_[0]) * x[DM3_0]
                  + (la_[0] * x[E_0]  + lc_[0] * x[E_0] * x[E_0] + m_[0] + p_ * q_[0] * x[E_1]) * x[DA_3]
                  + (1 - p_) * q_[0] * x[DM3_1] + g_[1] * x[DM3_1];

    dxdt[DM3_1] = -(lc_[1] + la_[1] + m_[1] + g_[0] + q_[1]) * x[DM3_1]
                  + (la_[1] * x[E_1]  + lc_[1] * x[E_1] * x[E_1] + m_[1] + p_ * q_[1] * x[E_0]) * x[DA_3]
                  + (1 - p_) * q_[1] * x[DM3_0] + g_[0] * x[DM3_0];


    dxdt[E_0]   = m_[0] - (lc_[0] + m_[0] + q_[0]) * x[E_0] +  lc_[0] * x[E_0] * x[E_0] + q_[0] * x[E_1];


    dxdt[E_1]   = m_[1] - (lc_[1] + m_[1] + q_[1]) * x[E_1] +  lc_[1] * x[E_1] * x[E_1] + q_[1] * x[E_0];

    dxdt[DA_3]  = -(g_[0] + g_[1]) * x[DA_3] + g_[0] * x[DM3_0] + g_[1] * x[DM3_1];
  }
};
} // namespace secsse
