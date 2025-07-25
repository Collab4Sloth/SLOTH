/*
 * Copyright © CEA 2023
 *
 * \brief Analytical potentials and derivatives used by phase-field models
 *
 * \file PhaseFieldPotentials.hpp
 * \author ci230846
 * \date 20/03/2023
 */
// TODO(ci) mettre de l'ordre dans ce qui est autorisé et non pour la manipulation des potentiels

#pragma once
#include <algorithm>
#include <functional>
#include <string>
#include <vector>

#include "Utils/Utils.hpp"

template <int ORDER, ThermodynamicsPotentialDiscretization SCHEME>
struct potential_function {};

template <int ORDER, ThermodynamicsPotentialDiscretization SCHEME,
          ThermodynamicsPotentials POTENTIAL>
class PotentialFunctions {
 private:
  template <class... Args>
  std::function<double(const double&)> getLog(Args... args) {
    potential_function<ORDER, SCHEME> func;
    return func.getLog(args...);
  }
  template <class... Args>
  std::function<double(const double&)> getW(Args... args) {
    potential_function<ORDER, SCHEME> func;
    return func.getW(args...);
  }
  template <class... Args>
  std::function<double(const double&)> getWW(Args... args) {
    potential_function<ORDER, SCHEME> func;
    return func.getWW(args...);
  }

  template <class... Args>
  std::function<double(const double&)> getF(Args... args) {
    potential_function<ORDER, SCHEME> func;
    return func.getF(args...);
  }

  template <class... Args>
  std::function<double(const double&)> getH(Args... args) {
    potential_function<ORDER, SCHEME> func;
    return func.getH(args...);
  }

  template <class... Args>
  std::function<double(const double&)> getX(Args... args) {
    potential_function<ORDER, SCHEME> func;
    return func.getX(args...);
  }

 public:
  PotentialFunctions();

  template <class... Args>
  std::function<double(const double&)> getPotentialFunction(Args... args);

  ~PotentialFunctions();
};

//////////////////////////////////////////////////////
// IMPLICIT
//////////////////////////////////////////////////////
///////////////////////
// ORDER  = 0
///////////////////////
template <>
struct potential_function<0, ThermodynamicsPotentialDiscretization::Implicit> {
  /**
   * @brief Logarithmic potential W(x)=x*log(x) + (1-x)*log(1-x)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLog(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto eps = 1.e-6;
      const auto xx = std::min(std::max(x, eps), 1. - eps);
      const auto pot = xx * log(xx) + (1.0 - xx) * log(1.0 - xx);
      return pot;
    });
  }
  /**
   * @brief Double Well potential W(x)=x² * (1-x)²
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getW(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = x * x * (1.0 - x) * (1.0 - x);
      return pot;
    });
  }
  template <typename... Args>
  std::function<double(const double&)> getWW(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const double a = 0.3;
      const double b = 0.7;
      const auto pot = (x - a) * (x - a) * (b - x) * (b - x);
      return pot;
    });
  }
  /**
   * @brief F potential F(x)=(x^4 -2x^2)/4
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getF(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 0.25 * (x * x - 1.0) * (x * x - 1.0);
      return pot;
    });
  }
  /**
   * @brief Interpolation potential H(x)=x³ * (6x²-15x+10)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getH(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = x * x * x * (6.0 * x * x - 15.0 * x + 10.0);
      return pot;
    });
  }
  /**
   * @brief Identity potential X(x)=x
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getX(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = x;
      return pot;
    });
  }
};
///////////////////////
// ORDER = 1
///////////////////////
template <>
struct potential_function<1, ThermodynamicsPotentialDiscretization::Implicit> {
  /**
   * @brief First derivative of the Logarithmic potential W(x)=x*log(x) + (1-x)*log(1-x)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLog(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto eps = 1.e-6;
      const auto xx = std::min(std::max(x, eps), 1. - eps);
      const auto pot = log(xx) - log(1.0 - xx);
      return pot;
    });
  }

  /**
   * @brief First derivative of the double Well potential W(x)=x² * (1-x)²
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getW(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 2. * x * (1.0 - x) * (1.0 - 2. * x);
      return pot;
    });
  }
  template <typename... Args>
  std::function<double(const double&)> getWW(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const double a = 0.3;
      const double b = 0.7;
      const auto pot = 2. * (x - a) * (b - x) * (a + b - 2. * x);
      return pot;
    });
  }
  /**
   * @brief First derivative of the F potential F(x)=(x^4 -2x^2)/4
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getF(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = x * x * x - x;
      return pot;
    });
  }
  /**
   * @brief First derivative of the interpolation potential H(x)=x³ * (6x²-15x+10)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getH(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 30. * x * x * (1.0 - x) * (1.0 - x);
      return pot;
    });
  }
  /**
   * @brief First derivative of the identity potential X(x)=x
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getX(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 1.;
      return pot;
    });
  }
};
///////////////////////
// ORDER = 2
///////////////////////
template <>
struct potential_function<2, ThermodynamicsPotentialDiscretization::Implicit> {
  /**
   * @brief Second derivative of the Logarithmic potential W(x)=x*log(x) + (1-x)*log(1-x)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLog(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto eps = 1.e-6;
      const auto xx = std::min(std::max(x, eps), 1. - eps);
      const auto pot = 1. / (xx * (1. - xx));
      return pot;
    });
  }
  /**
   * @brief Second derivative of the double Well potential W(x)=x² * (1-x)²
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getW(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 2. * (1. - 6. * x + 6. * x * x);
      return pot;
    });
  }
  template <typename... Args>
  std::function<double(const double&)> getWW(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const double a = 0.3;
      const double b = 0.7;
      const auto pot = 2. * (a * a + b * b + 4. * a * b - 6. * x * a - 6. * x * b + 6. * x * x);
      return pot;
    });
  }
  /**
   * @brief Second derivative of the F potential F(x)=(x^4 -2x^2)/4
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getF(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 3. * x * x - 1.0;
      return pot;
    });
  }
  /**
   * @brief Second derivative of the interpolation potential H(x)=x³ * (6x²-15x+10)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getH(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 60. * x * (1.0 - x) * (1.0 - 2. * x);
      return pot;
    });
  }
  /**
   * @brief Second derivative of the identity potential X(x)=x
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getX(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 0.;
      return pot;
    });
  }
};
//////////////////////////////////////////////////////
// EXPLICIT
//////////////////////////////////////////////////////
///////////////////////
// ORDER  = 0
///////////////////////
template <>
struct potential_function<0, ThermodynamicsPotentialDiscretization::Explicit> {
  /**
   * @brief Logarithmic potential W(x)=x*log(x) + (1-x)*log(1-x)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLog(Args... args) {
    mfem::mfem_error("potential_function::getLog: explicit scheme not available");
  }
  /**
   * @brief Double Well potential W(x)=x² * (1-x)²
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getW(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = xn * xn * (1.0 - xn) * (1.0 - xn);
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getW: only one argument is expected for explicit scheme");
    }
  }
  template <typename... Args>
  std::function<double(const double&)> getWW(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = xn * xn * (1.0 - xn) * (1.0 - xn);
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getW: only one argument is expected for explicit scheme");
    }
  }
  /**
   * @brief F potential F(x)=(x² -1)/4
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getF(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = 0.25 * (xn * xn - 1.0) * (xn * xn - 1.0);
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getF: only one argument is expected for explicit scheme");
    }
  }
  /**
   * @brief Interpolation potential H(x)=x³ * (6x²-15x+10)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getH(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = xn * xn * xn * (6.0 * xn * xn - 15.0 * xn + 10.0);
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getH: only one argument is expected for explicit scheme");
    }
  }
  /**
   * @brief Identity potential X(x)=x
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getX(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = xn;
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getX: only one argument is expected for explicit scheme");
    }
  }
};
///////////////////////
// ORDER = 1
///////////////////////
template <>
struct potential_function<1, ThermodynamicsPotentialDiscretization::Explicit> {
  /**
   * @brief First derivative Logarithmic potential W(x)=x*log(x) + (1-x)*log(1-x)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLog(Args... args) {
    mfem::mfem_error("potential_function::getLog: explicit scheme not available");
  }
  /**
   * @brief First derivative of the double Well potential W(x)=x² * (1-x)²
   *        with explicit scheme (as implicit scheme)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getW(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = 2. * xn * (1.0 - xn) * (1.0 - 2. * xn);
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getW: only one argument is expected for explicit scheme");
    }
  }
  template <typename... Args>
  std::function<double(const double&)> getWW(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = 2. * xn * (1.0 - xn) * (1.0 - 2. * xn);
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getW: only one argument is expected for explicit scheme");
    }
  }

  /**
   * @brief First derivative of the F potential F(x)=(x² -1)/4
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getF(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = xn * xn * xn - xn;
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getF: only one argument is expected for explicit scheme");
    }
  }
  /**
   * @brief First derivative of the interpolation potential H(x)=x³ * (6x²-15x+10)
   *        with explicit scheme (as implicit scheme)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getH(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = 30. * xn * xn * (1.0 - xn) * (1.0 - xn);
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getH: only one argument is expected for explicit scheme");
    }
  }
  /**
   * @brief First derivative of the identity potential X(x)=x
   *        with explicit scheme (as implicit scheme)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getX(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 1.;
      return pot;
    });
  }
};
///////////////////////
// ORDER = 2
///////////////////////
template <>
struct potential_function<2, ThermodynamicsPotentialDiscretization::Explicit> {
  /**
   * @brief Second derivative Logarithmic potential W(x)=x*log(x) + (1-x)*log(1-x)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLog(Args... args) {
    mfem::mfem_error("potential_function::getLog: explicit scheme not available");
  }
  /**
   * @brief Second derivative of the double Well potential W(x)=x² * (1-x)²
   *        with explicit scheme (as implicit scheme)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getW(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 0.;  // 2. * (1. - 6. * x + 6. * x * x);
      return pot;
    });
  }
  template <typename... Args>
  std::function<double(const double&)> getWW(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 0.;  // 2. * (1. - 6. * x + 6. * x * x);
      return pot;
    });
  }
  /**
   * @brief Second derivative of the F potential F(x)=(x² -1)/4
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getF(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 0.;  // 3. * x * x - 1.0;
      return pot;
    });
  }
  /**
   * @brief Second derivative of the interpolation potential H(x)=x³ * (6x²-15x+10)
   *        with explicit scheme (as implicit scheme)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getH(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 0.;  // 60. * x * (1.0 - x) * (1.0 - 2. * x);
      return pot;
    });
  }
  /**
   * @brief Second derivative of the identity potential X(x)=x
   *        with explicit scheme (as implicit scheme)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getX(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 0.;
      return pot;
    });
  }
};
//////////////////////////////////////////////////////
// SEMI-IMPLICIT
//////////////////////////////////////////////////////
///////////////////////
// ORDER  = 0
///////////////////////
template <>
struct potential_function<0, ThermodynamicsPotentialDiscretization::SemiImplicit> {
  /**
   * @brief Logarithmic potential W(x)=x*log(x) + (1-x)*log(1-x)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLog(Args... args) {
    mfem::mfem_error("potential_function::getLog: semi-explicit scheme not available");
  }

  /**
   * @brief Double Well potential W(x)=x² * (1-x)²
   *        with semi-implicit scheme (as implicit/explicit schemes)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getW(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = x * x * (1.0 - x) * (1.0 - x);
      return pot;
    });
  }
  template <typename... Args>
  std::function<double(const double&)> getWW(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = x * x * (1.0 - x) * (1.0 - x);
      return pot;
    });
  }
  /**
   * @brief F potential F(x)=(x² -1)/4
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getF(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 0.25 * (x * x - 1.0) * (x * x - 1.0);
      return pot;
    });
  }
  /**
   * @brief Interpolation potential H(x)=x³ * (6x²-15x+10)
   *        with semi-implicit scheme (as implicit/explicit schemes)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getH(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = x * x * x * (6.0 * x * x - 15.0 * x + 10.0);
      return pot;
    });
  }
  /**
   * @brief Identity potential X(x)=x
   *        with semi-implicit scheme (as implicit/explicit schemes)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getX(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = x;
      return pot;
    });
  }
};
///////////////////////
// ORDER = 1
///////////////////////
template <>
struct potential_function<1, ThermodynamicsPotentialDiscretization::SemiImplicit> {
  /**
   * @brief First derivative of Logarithmic potential W(x)=x*log(x) + (1-x)*log(1-x)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLog(Args... args) {
    mfem::mfem_error("potential_function::getLog: semi-explicit scheme not available");
  }
  /**
   * @brief First derivative of the double Well potential W(x)=x² * (1-x)²
   *        with semi-implicit scheme
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getW(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = (1.0 - x - xn) * (x + xn - x * x - xn * xn);
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getW: only one argument is expected for smei-implicit scheme");
    }
  }
  template <typename... Args>
  std::function<double(const double&)> getWW(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = (1.0 - x - xn) * (x + xn - x * x - xn * xn);
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getW: only one argument is expected for smei-implicit scheme");
    }
  }
  /**
   * @brief First derivative of the F potential F(x)=(x² -1)/4
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getF(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = x * (x * x - 1.0);
      return pot;
    });
  }
  /**
   * @brief First derivative of the interpolation potential H(x)=x³ * (6x²-15x+10)
   *        with semi-implicit scheme (as implicit/explicit schemes)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getH(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 30. * x * x * (1.0 - x) * (1.0 - x);
      return pot;
    });
  }
  /**
   * @brief First derivative of the identity potential X(x)=x
   *        with semi-implicit scheme (as implicit/explicit schemes)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getX(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 1.;
      return pot;
    });
  }
};
///////////////////////
// ORDER = 2
///////////////////////
template <>
struct potential_function<2, ThermodynamicsPotentialDiscretization::SemiImplicit> {
  /**
   * @brief Second derivative of Logarithmic potential W(x)=x*log(x) + (1-x)*log(1-x)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getLog(Args... args) {
    mfem::mfem_error("potential_function::getLog: semi-explicit scheme not available");
  }
  /**
   * @brief Second derivative of the double Well potential W(x)=x² * (1-x)²
   *        with semi-implicit scheme
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getW(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = ((1.0 - x - xn) * (1.0 - 2.0 * x) - (x + xn - x * x - xn * xn));
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getW: only one argument is expected for smei-implicit scheme");
    }
  }
  template <typename... Args>
  std::function<double(const double&)> getWW(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = ((1.0 - x - xn) * (1.0 - 2.0 * x) - (x + xn - x * x - xn * xn));
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getW: only one argument is expected for smei-implicit scheme");
    }
  }
  /**
   * @brief Second derivative of the F potential F(x)=(x² -1)/4
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getF(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 3. * x * x - 1.0;
      return pot;
    });
  }
  /**
   * @brief Second derivative of the interpolation potential H(x)=x³ * (6x²-15x+10)
   *        with semi-implicit scheme
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getH(Args... args) {
    auto v = std::vector<double>{args...};

    if (v.size() == 1) {
      const auto xn = v[0];
      return std::function<double(const double&)>([xn](double x) {
        const auto pot = 30. * (1.0 - x - xn) * (x + xn - x * x - xn * xn);
        return pot;
      });
    } else {
      mfem::mfem_error(
          "potential_function::getH: only one argument is expected for smei-implicit scheme");
    }
  }
  /**
   * @brief Second derivative of the identity potential X(x)=x
   *        with semi-implicit scheme (as implicit/explicit schemes)
   *
   * @tparam Args
   * @param args
   * @return std::function<double(const double&)>
   */
  template <typename... Args>
  std::function<double(const double&)> getX(Args... args) {
    return std::function<double(const double&)>([](double x) {
      const auto pot = 0.;
      return pot;
    });
  }
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Construct a new potential function:: potential function object
 *
 */
template <int ORDER, ThermodynamicsPotentialDiscretization SCHEME,
          ThermodynamicsPotentials POTENTIAL>
PotentialFunctions<ORDER, SCHEME, POTENTIAL>::PotentialFunctions() {}

/**
 * @brief return the function associated with the POTENTIAL, its ORDER of
 * derivative and its time discretization scheme
 *
 * @tparam ORDER
 * @tparam SCHEME
 * @tparam POTENTIAL
 * @tparam Args
 * @param args
 * @return std::function<double(const double&)>
 */
template <int ORDER, ThermodynamicsPotentialDiscretization SCHEME,
          ThermodynamicsPotentials POTENTIAL>
template <class... Args>
std::function<double(const double&)>
PotentialFunctions<ORDER, SCHEME, POTENTIAL>::getPotentialFunction(Args... args) {
  switch (POTENTIAL) {
    case ThermodynamicsPotentials::LOG:
      return this->getLog(args...);
    case ThermodynamicsPotentials::W:
      return this->getW(args...);
    case ThermodynamicsPotentials::WW:
      return this->getWW(args...);
    case ThermodynamicsPotentials::F:
      return this->getF(args...);
    case ThermodynamicsPotentials::H:
      return this->getH(args...);
    case ThermodynamicsPotentials::X:
      return this->getX(args...);
    default:
      mfem::mfem_error(
          "PotentialFunctions::getPotentialFunctions: double well, H interpolation and identity "
          "potential function  are available");
      break;
  }
}

/**
 * @brief Destroy the potential function :: potential function  object
 *
 */
template <int ORDER, ThermodynamicsPotentialDiscretization SCHEME,
          ThermodynamicsPotentials POTENTIAL>
PotentialFunctions<ORDER, SCHEME, POTENTIAL>::~PotentialFunctions() {}
