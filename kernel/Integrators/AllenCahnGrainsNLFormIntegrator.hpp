/**
 * @file AllenCahnGrainNLFormIntegrator.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief VF for Allen-Cahn equation solved in a polycrystal
 * @todo Extend for different interaction coefficient
 * @version 0.1
 * @date 2025-06-22
 *
 * Copyright CEA (c) 2025
 *
 */
#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Coefficients/LambdaCoefficient.hpp"
#include "Coefficients/MobilityCoefficient.hpp"
#include "Coefficients/OmegaCoefficient.hpp"
#include "Coefficients/PhaseFieldPotentials.hpp"
#include "Integrators/SlothNLFormIntegrator.hpp"
#include "Parameters/Parameters.hpp"
#include "Utils/Utils.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]

#pragma once

/**
 * @brief Class dedicated to the FV of the Allen Cahn equation
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
class AllenCahnGrainNLFormIntegrator : public mfem::BlockNonlinearFormIntegrator,
                                       public SlothNLFormIntegrator<VARS> {
 private:
  mfem::DenseMatrix gradPsi;
  mfem::Vector Psi, gradU;

  PotentialFunctions<1, SCHEME, ENERGY> energy_first_derivative_potential_;
  PotentialFunctions<2, SCHEME, ENERGY> energy_second_derivative_potential_;

  template <typename... Args>
  double mobility(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
                  const Parameters& parameters);

  template <typename... Args>
  double omega(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
               const Parameters& parameters);

  template <typename... Args>
  double lambda(mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
                const Parameters& parameters);

  FType double_well_derivative(const int order_derivative, mfem::ElementTransformation& Tr,
                               const mfem::IntegrationPoint& ir);

  void check_variables_consistency();

 protected:
  std::vector<mfem::ParGridFunction> u_old_;
  std::vector<mfem::ParGridFunction> aux_gf_;
  std::vector<mfem::Vector> aux_old_gf_;
  std::vector<std::vector<std::string>> aux_gf_infos_;
  std::vector<mfem::ParGridFunction> temp_gf_;
  bool scale_mobility_by_temperature_{false};

  virtual FType energy_derivatives(const int order_derivative, mfem::ElementTransformation& Tr,
                                   const mfem::IntegrationPoint& ir);

 public:
  AllenCahnGrainNLFormIntegrator(const std::vector<mfem::ParGridFunction>& u_old,
                                 const Parameters& params, std::vector<VARS*> auxvars);
  ~AllenCahnGrainNLFormIntegrator();

  virtual void AssembleElementVector(const mfem::Array<const mfem::FiniteElement*>& el,
                                     mfem::ElementTransformation& Tr,
                                     const mfem::Array<const mfem::Vector*>& elfun,
                                     const mfem::Array<mfem::Vector*>& elvec);

  virtual void AssembleElementGrad(const mfem::Array<const mfem::FiniteElement*>& el,
                                   mfem::ElementTransformation& Tr,
                                   const mfem::Array<const mfem::Vector*>& elfun,
                                   const mfem::Array2D<mfem::DenseMatrix*>& elmats);

  std::unique_ptr<HomogeneousEnergyCoefficient<ENERGY>> get_energy(
      std::vector<mfem::ParGridFunction*> gfu, const double omega);
  std::unique_ptr<GradientEnergyCoefficient> get_grad_energy(
      std::vector<mfem::ParGridFunction*> gfu, const double lambda);
};
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

/**
 * @brief Energy derivatives contribution
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param order_derivative
 * @return FuncType
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FType AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::energy_derivatives(
    const int order_derivative, mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  return [this, order_derivative, &Tr, &ir](const double& u) {
    return this->double_well_derivative(order_derivative, Tr, ir)(u);
  };
}

/**
 * @brief Return the value of the mobility coefficient at integration point
 *
 * @remark Actually, only explicit mobility is available. (See diffusion integrator  to extend
 * towards implicit mobility. )
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param Tr
 * @param ip
 * @param u
 * @param parameters
 * @return double
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
template <typename... Args>
double AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::mobility(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  const int nElement = Tr.ElementNo;
  // TODO(cci) return vector of double
  MobilityCoefficient<0, MOBI> mobi_coeff(&this->u_old_[0], parameters);
  double mob_coeff = mobi_coeff.Eval(Tr, ip);
  if (this->scale_mobility_by_temperature_) {
    mob_coeff /= this->temp_gf_[0].GetValue(nElement, ip);
  }
  return mob_coeff;
}

/**
 * @brief Return the value of the omega coefficient at integration point
 *
 * @remark Actually, only constant omega is available.
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param Tr
 * @param ip
 * @param u
 * @param parameters
 * @return double
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
template <typename... Args>
double AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::omega(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  OmegaCoefficient<0, Omega::Constant> omega_coeff(&this->u_old_[0], parameters);
  return omega_coeff.Eval(Tr, ip);
}
/**
 * @brief Return the value of the lambda coefficient at integration point
 *
 * @remark Actually, only constant lambda is available.
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param Tr
 * @param ip
 * @param u
 * @param parameters
 * @return double
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
template <typename... Args>
double AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::lambda(
    mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ip, const double u,
    const Parameters& parameters) {
  LambdaCoefficient<0, Lambda::Constant> lambda_coeff(&this->u_old_[0], parameters);
  return lambda_coeff.Eval(Tr, ip);
}

/**
 * @brief double_well_derivative contribution
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param order_derivative
 * @return std::function<double(const double&, const double&)>
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
FType AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::double_well_derivative(
    const int order_derivative, mfem::ElementTransformation& Tr, const mfem::IntegrationPoint& ir) {
  return FType([this, order_derivative, &Tr, &ir](const double& u) {
    const auto& un = this->u_old_[0].GetValue(Tr, ir);
    FType W_derivative;
    if (order_derivative == 1) {
      W_derivative = this->energy_first_derivative_potential_.getPotentialFunction(un);
    } else if (order_derivative == 2) {
      W_derivative = this->energy_second_derivative_potential_.getPotentialFunction(un);
    } else {
      std::runtime_error("Error while setting the order of derivative : only 1 and 2 are allowed.");
    }
    const double omega = this->omega(Tr, ir, u, this->params_);

    const auto& w_prime = omega * W_derivative(u);
    return w_prime;
  });
}

/**
 * @brief Construct a new AllenCahnGrainNLFormIntegrator object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param u_old
 * @param omega
 * @param lambda
 * @param mob
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::AllenCahnGrainNLFormIntegrator(
    const std::vector<mfem::ParGridFunction>& u_old, const Parameters& params,
    std::vector<VARS*> auxvars)
    : SlothNLFormIntegrator<VARS>(params, auxvars), u_old_(u_old) {
  this->check_variables_consistency();
}

/**
 * @brief Check variables consistency
 *
 * @tparam VARS
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::check_variables_consistency() {
  this->aux_gf_ = this->get_aux_gf();
  this->aux_old_gf_ = this->get_aux_old_gf();
  this->aux_gf_infos_ = this->get_aux_infos();

  // Temperature scaling for mobility
  bool temperature_found = false;
  for (std::size_t i = 0; i < this->aux_gf_infos_.size(); ++i) {
    const auto& variable_info = this->aux_gf_infos_[i];
    MFEM_VERIFY(!variable_info.empty(), "Empty variable_info encountered.");
    size_t vsize = variable_info.size();

    MFEM_VERIFY(vsize >= 1,
                "AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>: at least "
                "one additionnal information is expected for auxiliary variables associated with "
                "this integrator");
    const std::string& symbol = toUpperCase(variable_info.back());
    if (symbol == "T") {
      this->temp_gf_.emplace_back(std::move(this->aux_gf_[i]));
      temperature_found = true;
      break;
    }
  }
  if (this->params_.has_parameter("ScaleMobilityByTemperature")) {
    this->scale_mobility_by_temperature_ =
        this->params_.template get_param_value<bool>("ScaleMobilityByTemperature");
    if (this->scale_mobility_by_temperature_) {
      MFEM_VERIFY(
          temperature_found,
          "AllenCahnGrainNLFormIntegrator: "
          "Temperature variable required to scale mobility, but not found in auxiliary variables");
    }
  }
}

/**
 * @brief Residual part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param el
 * @param Tr
 * @param elfun
 * @param elvect
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::AssembleElementVector(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun, const mfem::Array<mfem::Vector*>& elvect) {
  int num_blocks = el.Size();
  for (int blk = 0; blk < num_blocks; ++blk) {
    // Catch_Time_Section("AllenCahnGrainNLFormIntegrator:AssembleElementVector");
    int nd = el[blk]->GetDof();
    int dim = el[blk]->GetDim();
    gradPsi.SetSize(nd, dim);
    Psi.SetSize(nd);
    gradU.SetSize(dim);
    // elvect.SetSize(nd);
    elvect[blk]->SetSize(nd);
    *elvect[blk] = 0.;

    const mfem::IntegrationRule* ir =
        &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());
    // elvect = 0.0;
    for (int i = 0; i < ir->GetNPoints(); i++) {
      const mfem::IntegrationPoint& ip = ir->IntPoint(i);
      el[blk]->CalcShape(ip, Psi);  //
      Tr.SetIntPoint(&ip);

      const auto& u = *elfun[blk] * Psi;
      // Interaction term : 2 * u_i (sum_j!=i u__j * u_j)
      double interaction = 0.0;
      for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
        interaction += *elfun[off_blk] * Psi;
      }
      interaction -= u * u;
      interaction *= 2.0 * u;

      // Laplacian : given u, compute (grad(u), grad(psi)), psi is shape function.
      // given u (elfun), compute grad(u)
      el[blk]->CalcPhysDShape(Tr, gradPsi);
      gradPsi.MultTranspose(*elfun[blk], gradU);
      const double coef_mob = this->mobility(Tr, ip, u, this->params_) * ip.weight * Tr.Weight();
      const double lambda = this->lambda(Tr, ip, u, this->params_);
      gradU *= coef_mob * lambda;
      gradPsi.AddMult(gradU, *elvect[blk]);

      // Given u, compute (w'(u) + interaction term , psi), psi is shape function
      const double ww = coef_mob * (this->energy_derivatives(1, Tr, ip)(u) + interaction);
      add(*elvect[blk], ww, Psi, *elvect[blk]);
    }
  }
}

/**
 * @brief Jacobian part of the non linear problem
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param el
 * @param Tr
 * @param elfun
 * @param elmat
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
void AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::AssembleElementGrad(
    const mfem::Array<const mfem::FiniteElement*>& el, mfem::ElementTransformation& Tr,
    const mfem::Array<const mfem::Vector*>& elfun,
    const mfem::Array2D<mfem::DenseMatrix*>& elmats) {
  // Catch_Time_Section("AllenCahnGrainNLFormIntegrator::AssembleElementGrad");
  int num_blocks = el.Size();

  // Diagonal
  for (int blk = 0; blk < num_blocks; ++blk) {
    for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
      // int nd = el.GetDof();
      // int dim = el.GetDim();
      int nd = el[blk]->GetDof();
      int dim = el[blk]->GetDim();

      gradPsi.SetSize(nd, dim);
      Psi.SetSize(nd);
      // elmat.SetSize(nd);
      elmats(blk, off_blk)->SetSize(nd);
      *elmats(blk, off_blk) = 0.0;

      const mfem::IntegrationRule* ir =
          &mfem::IntRules.Get(el[blk]->GetGeomType(), 2 * el[blk]->GetOrder() + Tr.OrderW());

      // elmat = 0.0;
      for (int i = 0; i < ir->GetNPoints(); i++) {
        const mfem::IntegrationPoint& ip = ir->IntPoint(i);
        el[blk]->CalcShape(ip, Psi);
        const auto& u = *elfun[blk] * Psi;
        const auto& v = *elfun[off_blk] * Psi;
        const double coef_mob = this->mobility(Tr, ip, u, this->params_) * ip.weight * Tr.Weight();

        double interaction = 0.0;

        Tr.SetIntPoint(&ip);
        double double_well = 0;
        if (off_blk == blk) {
          double_well = this->energy_derivatives(2, Tr, ip)(u);
          el[blk]->CalcPhysDShape(Tr, gradPsi);
          const double lambda = this->lambda(Tr, ip, u, this->params_);

          AddMult_a_AAt(coef_mob * lambda, gradPsi, *elmats(blk, off_blk));

          for (int off_blk = 0; off_blk < num_blocks; ++off_blk) {
            interaction += *elfun[off_blk] * Psi;
          }
          interaction -= u * u;
          interaction *= 2.0;
        } else {
          interaction = 4.0 * u * v;
        }

        double fun_val = coef_mob * (double_well + interaction);
        AddMult_a_VVt(fun_val, Psi, *elmats(blk, off_blk));  // w'(u)*(du, psi)
      }
    }
  }
}

/**
 * @brief Return the energy coefficient associated with the integrator
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @tparam Args
 * @param gfu
 * @param lambda
 * @param omega
 * @return EnergyCoefficient<SCHEME, ENERGY>
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
std::unique_ptr<HomogeneousEnergyCoefficient<ENERGY>>
AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::get_energy(
    std::vector<mfem::ParGridFunction*> gfu, const double omega) {
  return std::make_unique<HomogeneousEnergyCoefficient<ENERGY>>(gfu[0], omega);
}

/**
 * @brief  Return the gradient energy coefficient associated with the integrator
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 * @param gfu
 * @param lambda
 * @return GradientEnergyCoefficient
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
std::unique_ptr<GradientEnergyCoefficient>
AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::get_grad_energy(
    std::vector<mfem::ParGridFunction*> gfu, const double lambda) {
  return std::make_unique<GradientEnergyCoefficient>(gfu[0], lambda);
}

/**
 * @brief Destroy the AllenCahnGrainNLFormIntegrator  object
 *
 * @tparam SCHEME
 * @tparam ENERGY
 * @tparam MOBI
 */
template <class VARS, ThermodynamicsPotentialDiscretization SCHEME, ThermodynamicsPotentials ENERGY,
          Mobility MOBI>
AllenCahnGrainNLFormIntegrator<VARS, SCHEME, ENERGY, MOBI>::~AllenCahnGrainNLFormIntegrator() {}
