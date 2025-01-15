/**
 * @file Calphad_problem.hpp
 * @author ci230846  (clement.introini@cea.fr)
 * @brief Class used to defined a Calphad_Problem objet
 * @version 0.1
 * @date 2025-01-07
 *
 * Copyright CEA (c) 2025
 *
 */
#include <functional>
#include <list>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "Convergence/PhysicalConvergence.hpp"
#include "Parameters/Parameter.hpp"
#include "Problems/ProblemBase.hpp"
#include "mfem.hpp"  // NOLINT [no include the directory when naming mfem include file]


#pragma once

class NonLinearSys : public mfem::Operator {

    public : 
        mutable mfem::DenseMatrix J;
        NonLinearSys() : Operator(2), J(2,2) {}
    

        double varphi;
        double c;
    double interpolationp(double varphi) const
{
   // return std::pow(varphi,3) * (10 - 15 * varphi + 6 * std::pow(varphi,2)); 
   return varphi;
}
double mu(double c, double ceq) const
{
    return ((c-ceq));
}
    virtual void Mult(const mfem::Vector &x, mfem::Vector &y) const override
    {
        
        double x1 = x(0);
        double x2 = x(1);
                
        y(0) = interpolationp(this->varphi) * x1 + (1 - interpolationp(this->varphi)) * x2 - this->c;
        y(1) = mu(x1,0.2) - mu(x2,0.8);

    }
    
    mfem::Operator& GetGradient(const mfem::Vector &x) const override {
    
        double x1 = x(0);
        double x2 = x(1);
        double eps = 1e-10;

       J(0,0) = interpolationp(this->varphi) ;
       J(0,1) =  1 - interpolationp(this->varphi);
       J(1,0) = 1.;
       J(1,1) = - 1.;

        return J;
            
    }

    void get_infos(const double& varphi, std::vector<double>& c_) {

      this->varphi = varphi;
      this->c = c_[0];
      
    }

}; // end class






template <class VAR, class PST>
class InterfaceProblem : public ProblemBase<VAR, PST> {
 private:

  void check_variables_consistency();

  std::vector<mfem::Vector> get_tp_conditions();
  std::vector<std::tuple<std::string, std::string>> get_chemical_system();
  std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
  get_output_system(const std::vector<std::vector<std::string>>& unks_info,
                    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk);

 public:
  template <class... Args>
  InterfaceProblem(const Parameters& params, VAR& variables, PST& pst,
                  const PhysicalConvergence& convergence, std::list<int> pop_elem,
                  Args&&... auxvariable);

  template <class... Args>
  InterfaceProblem(const std::string& name, const Parameters& params, VAR& variables, PST& pst,
                  const PhysicalConvergence& convergence, std::list<int> pop_elem,
                  Args&&... auxvariable);

  template <class... Args>
  InterfaceProblem(const Parameters& params, VAR& variables, PST& pst,
                  const PhysicalConvergence& convergence, Args&&... auxvariable);

  template <class... Args>
  InterfaceProblem(const std::string& name, VAR& variables, PST& pst,
                  const PhysicalConvergence& convergence, Args&&... auxvariable);

  /////////////////////////////////////////////////////
  void initialize(const double& initial_time) override;

  /////////////////////////////////////////////////////

  void do_time_step(double& next_time, const double& current_time, double current_time_step,
                    const int iter, std::vector<std::unique_ptr<mfem::Vector>>& unks,
                    const std::vector<std::vector<std::string>>& unks_info) override;

  void get_parameters();

  /////////////////////////////////////////////////////
  void finalize() override;
  /////////////////////////////////////////////////////

  ~InterfaceProblem();
};

////////////////////////////////
////////////////////////////////

/**
 * @brief Construct a new  Calphad_Problem object (with auxiliary variables)
 * @warning At least two auxiliary variable. The first auxiliary variable is the temperature and the
 * seccond is the pressure.

 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param params
 * @param variables
 * @param pst
 * @param convergence
 * @param pop_elem
 * @param auxvariables
 */
template <class VAR, class PST>
template <class... Args>
InterfaceProblem< VAR, PST>::InterfaceProblem(const Parameters& params, VAR& variables,
                                                    PST& pst,
                                                    const PhysicalConvergence& convergence,
                                                    std::list<int> pop_elem, Args&&... auxvariables)
    : ProblemBase<VAR, PST>(" Problem", variables, pst, convergence, pop_elem,
                            auxvariables...) {
}

/**
 * @brief Construct a new InterfaceProblem<VAR, PST>::InterfaceProblem object
 *
 * @tparam 
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param name
 * @param params
 * @param variables
 * @param pst
 * @param convergence
 * @param pop_elem
 * @param auxvariables
 */
template <class VAR, class PST>
template <class... Args>
InterfaceProblem<VAR, PST>::InterfaceProblem(const std::string& name,
                                                    const Parameters& params, VAR& variables,
                                                    PST& pst,
                                                    const PhysicalConvergence& convergence,
                                                    std::list<int> pop_elem, Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, pst, convergence, pop_elem, auxvariables...) {
}

/**
 * @brief Construct a new _Problem<VAR, PST>::_Problem object
 *
 * @tparam 
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param params
 * @param variables
 * @param pst
 * @param convergence
 * @param auxvariables
 */
template <class VAR, class PST>
template <class... Args>
InterfaceProblem<VAR, PST>::InterfaceProblem(const Parameters& params, VAR& variables,
                                                    PST& pst,
                                                    const PhysicalConvergence& convergence,
                                                    Args&&... auxvariables)
    : ProblemBase<VAR, PST>(" Problem", variables, pst, convergence, auxvariables...) {
}

/**
 * @brief Construct a new InterfaceProblem<VAR, PST>::InterfaceProblem object
 *
 * @tparam 
 * @tparam VAR
 * @tparam PST
 * @tparam Args
 * @param name
 * @param params
 * @param variables
 * @param pst
 * @param convergence
 * @param auxvariables
 */
template <class VAR, class PST>
template <class... Args>
InterfaceProblem<VAR, PST>::InterfaceProblem(const std::string& name, VAR& variables,
                                                    PST& pst,
                                                    const PhysicalConvergence& convergence,
                                                    Args&&... auxvariables)
    : ProblemBase<VAR, PST>(name, variables, pst, convergence, auxvariables...) {
}

/**
 * @brief Check consistency of variables. Name of a variable must have one of the following prefix:
 * xx_ for molar fraction,
 * nn_ for mole number,
 * mu_ for chemical potential,
 * df_ for driving force
 *
 */
template <class VAR, class PST>
void InterfaceProblem<VAR, PST>::check_variables_consistency() {

}

/**
 * @brief Initialization of the problem depending on the initialization of the Base objet
 *
 * @tparam 
 * @tparam VAR
 * @tparam PST
 * @param initial_time
 */
template <class VAR, class PST>
void InterfaceProblem<VAR, PST>::initialize(const double& initial_time) {
}

/**
 * @brief  Do a time-step by calling Step method of the ODE

 *
 * @tparam VAR
 * @tparam PST
 * @param next_time
 * @param current_time
 * @param current_time_step
 * @param iter
 * @param vect_unk
 * @param unks_info
 */
template <class VAR, class PST>
void InterfaceProblem<VAR, PST>::do_time_step(
    double& next_time, const double& current_time, double current_time_step, const int iter,
    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk,
    const std::vector<std::vector<std::string>>& unks_info) {
  int rank = mfem::Mpi::WorldRank();
  mfem::real_t local_varphi, local_c;
  std::vector<mfem::Vector> tp_gf = this->get_tp_conditions();



        mfem::NewtonSolver ourSolver;

      std::unique_ptr<mfem::GMRESSolver> linearSolver = std::make_unique<mfem::GMRESSolver>();
      std::unique_ptr<NonLinearSys> F = std::make_unique<NonLinearSys>();

  for (unsigned i = 0; i < tp_gf[0].Size(); i++)
  {
    if (tp_gf[0](i) > 0.01 && tp_gf[0](i) < 0.99) {
      std::cout << "------------------------------" << std::endl;
      std::cout << "je suis dans l'interface ici" << std::endl;
      local_varphi = tp_gf[0](i);
      local_c = tp_gf[1](i);
      std::cout << "local c : " << local_c << std::endl;

      std::vector<double> lala = {local_c};
      F->get_infos(local_varphi,lala);
      ourSolver.SetOperator(*F);
      ourSolver.SetSolver(*linearSolver);
      ourSolver.SetRelTol(1e-8);
      ourSolver.SetAbsTol(1e-12);
      ourSolver.SetMaxIter(10000);
      ourSolver.SetPrintLevel(1);

      linearSolver->SetRelTol(1e-8);
      linearSolver->SetAbsTol(1e-12);
      linearSolver->SetMaxIter(10000);
      linearSolver->SetPrintLevel(0);

      mfem::Vector x(2);
      x(0) = local_c;
      x(1) = local_c;

      mfem::Vector b(2);
      b = 0.0;

      ourSolver.Mult(b, x);

      std::cout << "phi : " << local_varphi << " x_alpha : " << x(0) << " x_beta : " << x(1) << std::endl;
      double veri = local_varphi * x(0) + (1-local_varphi) * x(1) - local_c;
      std::cout << "Verification phi x c_a + (1-phi) x c_b - c =  " << veri << std::endl;

    }

  }



  const size_t unk_size = vect_unk.size();
  for (size_t i = 0; i < unk_size; i++) {
    auto& unk_i = *(vect_unk[i]);
    this->unknown_.emplace_back(unk_i);
  }
  next_time = current_time + current_time_step;
}



/**
 * @brief Get all parameters associated with the problem.
 *
 * @tparam 
 * @tparam VAR
 * @tparam PST
 */
template <class VAR, class PST>
void InterfaceProblem<VAR, PST>::get_parameters() {
}

/**
 * @brief Finalization of the problem depending on the finalization of the Base objet
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 */
template <class VAR, class PST>
void InterfaceProblem< VAR, PST>::finalize() {
}

/**
 * @brief Define the temperature and the pressure use to calculate  equilibria
 * @warning By convention temperature, pressure are given in the first Variables objet Composition
 * is  given in the second Variables objet
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @return std::vector<mfem::Vector>
 */
template <class VAR, class PST>
std::vector<mfem::Vector> InterfaceProblem< VAR, PST>::get_tp_conditions() {
  std::vector<mfem::Vector> aux_gf;

  for (const auto& auxvar_vec : this->auxvariables_) {
    for (const auto& auxvar : auxvar_vec->getVariables()) {
      const auto gf = auxvar.get_unknown();
      aux_gf.emplace_back(gf);
    }
  }

  return aux_gf;
}

/**
 * @brief Define the composition use to calculate  equilibria
 * @warning By convention temperature, pressure are given in the first Variables objet Composition
 * is  given in the second Variables objet
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @return std::vector<std::tuple<std::string, std::string>>
 */
template <class VAR, class PST>
std::vector<std::tuple<std::string, std::string>>
InterfaceProblem< VAR, PST>::get_chemical_system() {
  std::vector<std::tuple<std::string, std::string>> chemical_system;


  return chemical_system;
}

/**
 * @brief Define the output system
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 * @param unks_info
 * @param vect_unk
 * @return std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
 */
template <class VAR, class PST>
std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
InterfaceProblem< VAR, PST>::get_output_system(
    const std::vector<std::vector<std::string>>& unks_info,
    std::vector<std::unique_ptr<mfem::Vector>>& vect_unk) {
   std::vector<std::tuple<std::vector<std::string>, std::reference_wrapper<mfem::Vector>>>
      output_system;
  return output_system;
}

/**
 * @brief Destroy the InterfaceProblem< VAR, PST>::InterfaceProblem object
 *
 * @tparam CALPHAD
 * @tparam VAR
 * @tparam PST
 */
template <class VAR, class PST>
InterfaceProblem< VAR, PST>::~InterfaceProblem() {}
