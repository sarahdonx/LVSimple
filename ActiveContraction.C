// Copyright (c) 2011-2013, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of New York University nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// IBAMR INCLUDES
#include <ibtk/libmesh_utilities.h>

// LIBMESH INCLUDES
#include <libmesh/dof_map.h>
#include <libmesh/equation_systems.h>
#include <libmesh/quadrature.h>

// APPLICATION INCLUDES
#include <ActiveContraction.h>
#include <BoundaryConditions.h>
#include <MechanicsModel.h>

// STATIC VARIABLES
int ActiveContraction::act_system_num, ActiveContraction::T_system_num;
namespace  // private namespace
{
// Model constants.
static const double a = 0.35;            // dimensionless
static const double A1 = -29.0;          // dimensionless
static const double A2 = 138.0;          // dimensionless
static const double A3 = 129.0;          // dimensionless
static const double alpha_0 = 8.0;       // sec^-1
static const double alpha_1 = 30.0;      // sec^-1
static const double alpha_2 = 130.0;     // sec^-1
static const double alpha_3 = 625.0;     // sec^-1
static const double alpha_r1 = 2.0;      // sec^-1
static const double alpha_r2 = 1.75;     // sec^-1
static const double beta_0 = 4.9;        // dimensionless
static const double beta_1 = -4.0;       // dimensionless
static const double Ca_50_ref = 1.05;    // uM
static const double Ca_TRPN_max = 70.0;  // uM
static const double gamma_trpn = 2.0;    // dimensionless
static const double k_on = 100.0;        // uM^-1 sec^-1
static const double k_refoff = 200.0;    // sec^-1
static const double K_Z = 0.15;          // dimensionless
static const double n = 3.0;             // dimensionless
static const double n_r = 3.0;           // dimensionless
static const double T_ref = 56.2;        // kPa = 1000 N m^-2 = 10000 dyne cm^-2
static const double z_p = 0.85;          // dimensionless
}

// CLASS IMPLEMENTATION

void
ActiveContraction::update_active_tension_model_state_variables(
    EquationSystems* equation_systems,
    const double time,
    const double dt)
{
    const MeshBase& mesh = equation_systems->get_mesh();
    const int dim = mesh.mesh_dimension();
    libMesh::UniquePtr<QBase> qrule = QBase::build(QGAUSS, NDIM, CONSTANT);

    System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
    const DofMap& X_dof_map = X_system.get_dof_map();
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(X_dof_map.variable_type(d) == X_dof_map.variable_type(0));
    std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);
    libMesh::UniquePtr<FEBase> X_fe(FEBase::build(dim, X_dof_map.variable_type(0)));
    X_fe->attach_quadrature_rule(qrule.get());
    const std::vector<std::vector<VectorValue<double> > >& dphi_X = X_fe->get_dphi();

    System& U_system = equation_systems->get_system<System>(IBFEMethod::VELOCITY_SYSTEM_NAME);
    const DofMap& U_dof_map = U_system.get_dof_map();
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(U_dof_map.variable_type(d) == U_dof_map.variable_type(0));
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(U_dof_map.variable_type(d) == X_dof_map.variable_type(0));

    System& f0_system = equation_systems->get_system<System>(MechanicsModel::f0_system_num);
    const DofMap& f0_dof_map = f0_system.get_dof_map();
    for (unsigned d = 0; d < NDIM; ++d) TBOX_ASSERT(f0_dof_map.variable_type(d) == f0_dof_map.variable_type(0));
    std::vector<std::vector<unsigned int> > f0_dof_indices(NDIM);

    System& T_system = equation_systems->get_system<System>(T_system_num);
    const DofMap& T_dof_map = T_system.get_dof_map();
    std::vector<unsigned int> T_dof_indices;

    System& act_system = equation_systems->get_system<System>(act_system_num);
    const DofMap& act_dof_map = act_system.get_dof_map();
    std::vector<std::vector<unsigned int> > act_dof_indices(NUM_ACT_VARS);

    //X_system.solution->localize(*X_system.current_local_solution);
    //NumericVector<double>& X_data = *(X_system.current_local_solution);
    //X_data.close();
    NumericVector<double>* X_vec = X_system.solution.get();
    NumericVector<double>* X_ghost_vec = X_system.current_local_solution.get();
    X_vec->localize(*X_ghost_vec);

    //U_system.solution->localize(*U_system.current_local_solution);
    //NumericVector<double>& U_data = *(U_system.current_local_solution);
    //U_data.close();
    NumericVector<double>* U_vec = U_system.solution.get();
    NumericVector<double>* U_ghost_vec = U_system.current_local_solution.get();
    U_vec->localize(*U_ghost_vec);

    NumericVector<double>* f0_vec = f0_system.solution.get();
    NumericVector<double>* f0_ghost_vec = f0_system.current_local_solution.get();
    f0_vec->localize(*f0_ghost_vec);
    //f0_system.solution->localize(*f0_system.current_local_solution);
    //NumericVector<double>&  f0_data = *(f0_system.solution);
    T_system.solution->localize(*T_system.current_local_solution);
    NumericVector<double>&   T_data = *(T_system.solution);

    act_system.solution->localize(*act_system.current_local_solution);
    NumericVector<double>& act_data = *(act_system.solution);

    double Ca_i_max = -1.0e300;
    double Ca_i_min = +1.0e300;
    double Ca_b_max = -1.0e300;
    double Ca_b_min = +1.0e300;
    double T_max = -1.0e300;

    TensorValue<double> FF, dFF_dt;
    boost::multi_array<double,2> X_node, U_node;
    const MeshBase::const_element_iterator el_begin = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator el_end   = mesh.active_local_elements_end();
    for (MeshBase::const_element_iterator el_it = el_begin; el_it != el_end; ++el_it)
    {
        Elem* const elem = *el_it;

        X_fe->reinit(elem);
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            X_dof_map.dof_indices(elem, X_dof_indices[d], d);
        }

        for (unsigned int d = 0; d < NDIM; ++d)
        {
            f0_dof_map.dof_indices(elem, f0_dof_indices[d], d);
        }

        T_dof_map.dof_indices(elem, T_dof_indices, 0);

        for (unsigned int d = 0; d < NUM_ACT_VARS; ++d)
        {
            act_dof_map.dof_indices(elem, act_dof_indices[d], d);
        }

        const unsigned int n_qp = qrule->n_points();
        TBOX_ASSERT(n_qp == 1);
        const unsigned int qp = 0;

        get_values_for_interpolation(X_node, *X_ghost_vec, X_dof_indices);
        jacobian(FF,qp,X_node,dphi_X);

        get_values_for_interpolation(U_node, *U_ghost_vec, X_dof_indices);
        jacobian(dFF_dt,qp,U_node,dphi_X);

        VectorValue<double> f0;
        for (unsigned int d = 0; d < NDIM; ++d)
        {
            f0(d) = (*f0_ghost_vec)(f0_dof_indices[d][0]);  // piecewise constant representation
        }
        const VectorValue<double> f = FF*f0;

        const double lambda = f.norm(); // the size() is depcrated
        const double dlambda_dt = (0.5/lambda)*f0*((FF.transpose()*dFF_dt + dFF_dt.transpose()*FF)*f0);

        double Ca_i = act_data(act_dof_indices[CA_I_IDX][0]);
        double Ca_b = act_data(act_dof_indices[CA_B_IDX][0]);
        double Q1   = act_data(act_dof_indices[  Q1_IDX][0]);
        double Q2   = act_data(act_dof_indices[  Q2_IDX][0]);
        double Q3   = act_data(act_dof_indices[  Q3_IDX][0]);
        double z    = act_data(act_dof_indices[   Z_IDX][0]);

        NHS_RK2_step(Ca_i, Ca_b, Q1, Q2, Q3, z, lambda, dlambda_dt, time, dt);

        act_data.set(act_dof_indices[CA_I_IDX][0], Ca_i);
        act_data.set(act_dof_indices[CA_B_IDX][0], Ca_b);
        act_data.set(act_dof_indices[  Q1_IDX][0],   Q1);
        act_data.set(act_dof_indices[  Q2_IDX][0],   Q2);
        act_data.set(act_dof_indices[  Q3_IDX][0],   Q3);
        act_data.set(act_dof_indices[   Z_IDX][0],    z);

        const double zz = max(min(lambda,1.15),0.8);
        assert(n   == 3.0);
        assert(n_r == 3.0);
        const double z_p_n_r = z_p*z_p*z_p;
        const double K_Z_n_r = K_Z*K_Z*K_Z;

        const double Ca_50 = Ca_50_ref*(1.0+beta_1*(zz-1.0));
        const double Ca_TRPN_50 = Ca_TRPN_max*Ca_50/(Ca_50+(k_refoff/k_on)*(1.0-(1.0+beta_0*(zz-1.0))*0.5/gamma_trpn));

        const double Ca_TRPN_50_Ca_TRPN_max_n = (Ca_TRPN_50*Ca_TRPN_50*Ca_TRPN_50)/(Ca_TRPN_max*Ca_TRPN_max*Ca_TRPN_max);

        const double K1 = alpha_r2*(z_p_n_r/z_p)*n_r*K_Z_n_r/((z_p_n_r+K_Z_n_r)*(z_p_n_r+K_Z_n_r));
        const double K2 = alpha_r2*(z_p_n_r/(z_p_n_r+K_Z_n_r))*(1.0-n_r*K_Z_n_r/(z_p_n_r+K_Z_n_r));
        const double z_max = (alpha_0/Ca_TRPN_50_Ca_TRPN_max_n-K2)/(alpha_r1+K1+alpha_0/Ca_TRPN_50_Ca_TRPN_max_n);

        const double T_0_max = T_ref*(1.0+beta_0*(zz-1.0));
        const double T_0 = T_0_max*z/z_max;

        const double Q_sum = Q1+Q2+Q3;

        //updated active model 22/01/2014
        double T = (Q_sum < 0.0 ? T_0*(a*Q_sum+1.0)/(1.0-Q_sum) : T_0*(1.0+(2.0+a)*Q_sum)/(1.0+Q_sum));
        T = max(0.0,T); //make sure T is always greater than 0
		    T = min(T,gamma_trpn*T_ref);

        T_data.set(T_dof_indices[0], T);

        Ca_i_max = max(Ca_i,Ca_i_max);
        Ca_i_min = min(Ca_i,Ca_i_min);
        Ca_b_max = max(Ca_b,Ca_b_max);
        Ca_b_min = min(Ca_b,Ca_b_min);
        T_max = max(T,T_max);
    }
    T_data.close();
    act_data.close();

    //make sure it will update the solution in the local processor, which will be used to calculate the contribution to the stress
    T_system.solution->localize(*T_system.current_local_solution);
    act_system.solution->localize(*act_system.current_local_solution);


    plog << "Ca_i_max = " << SAMRAI_MPI::maxReduction(Ca_i_max) << "\n";
    plog << "Ca_i_min = " << SAMRAI_MPI::minReduction(Ca_i_min) << "\n";
    plog << "Ca_b_max = " << SAMRAI_MPI::maxReduction(Ca_b_max) << "\n";
    plog << "Ca_b_min = " << SAMRAI_MPI::minReduction(Ca_b_min) << "\n";
    plog << "T_max    = " << SAMRAI_MPI::maxReduction(T_max)    << "\n";
    return;
}// update_active_tension_model_state_variables

void
ActiveContraction::NHS_RK2_step(
    double& Ca_i,
    double& Ca_b,
    double& Q1,
    double& Q2,
    double& Q3,
    double& z,
    const double lambda,
    const double dlambda_dt,
    const double time,
    const double dt)
{
    double Ca_i_new = Ca_i;
    double Ca_b_new = Ca_b;
    double   Q1_new = Q1;
    double   Q2_new = Q2;
    double   Q3_new = Q3;
    double    z_new = z;
    NHS_euler_step(Ca_i_new, Ca_b_new, Q1_new, Q2_new, Q3_new, z_new, lambda, dlambda_dt, time   , dt);
    NHS_euler_step(Ca_i_new, Ca_b_new, Q1_new, Q2_new, Q3_new, z_new, lambda, dlambda_dt, time+dt, dt);
    Ca_i = 0.5*(Ca_i+Ca_i_new);
    Ca_b = 0.5*(Ca_b+Ca_b_new);
    Q1   = 0.5*(Q1  +  Q1_new);
    Q2   = 0.5*(Q2  +  Q2_new);
    Q3   = 0.5*(Q3  +  Q3_new);
    z    = 0.5*(z   +   z_new);
    return;
}// NHS_RK2_step

void
ActiveContraction::NHS_euler_step(
    double& Ca_i,
    double& Ca_b,
    double& Q1,
    double& Q2,
    double& Q3,
    double& z,
    const double lambda,
    const double dlambda_dt,
    const double time,
    const double dt)
{
    // Simple Ca transient.
    static const double Ca_max = 1.0;   // uM
    static const double Ca_o = 0.01;    // uM
    static const double tau_Ca = 0.06;  // sec

    //const double t_constant = 0.0; //decide to wait or not
    const double t_shift = time-BoundaryConditions::t_end_diastole;
    //v3
    //const double dCa_i_dt = t_shift > 0.0 ? ((Ca_max-Ca_o)/tau_Ca)*exp(1.0-t_shift/tau_Ca)*(1.0-t_shift/tau_Ca) : 0.0;

    //v4
    double dCa_i_dt = t_shift > 0.0 ? ((Ca_max-Ca_o)/tau_Ca)*exp(1.0-t_shift/tau_Ca)*(1.0-t_shift/tau_Ca) : 0.0;
    if (time>=0.0599 + BoundaryConditions::t_end_diastole)
    {
        dCa_i_dt = 0.0;
    }

    // The model is only valid for 0.8 <= lambda <= 1.15.
    const double zz = max(min(lambda,1.15),0.8);

    // Tropomyosin kinetics.
    assert(n   == 3.0);
    assert(n_r == 3.0);
    const double z_p_n_r = z_p*z_p*z_p;
    const double K_Z_n_r = K_Z*K_Z*K_Z;

    const double z_n_r = z*z*z;

    const double Ca_50 = Ca_50_ref*(1.0+beta_1*(zz-1.0));
    const double Ca_TRPN_50 = Ca_TRPN_max*Ca_50/(Ca_50+(k_refoff/k_on)*(1.0-(1.0+beta_0*(zz-1.0))*0.5/gamma_trpn));

    const double Ca_b_Ca_TRPN_50_n = (Ca_b*Ca_b*Ca_b)/(Ca_TRPN_50*Ca_TRPN_50*Ca_TRPN_50);
    const double Ca_TRPN_50_Ca_TRPN_max_n = (Ca_TRPN_50*Ca_TRPN_50*Ca_TRPN_50)/(Ca_TRPN_max*Ca_TRPN_max*Ca_TRPN_max);

    const double K1 = alpha_r2*(z_p_n_r/z_p)*n_r*K_Z_n_r/((z_p_n_r+K_Z_n_r)*(z_p_n_r+K_Z_n_r));
    const double K2 = alpha_r2*(z_p_n_r/(z_p_n_r+K_Z_n_r))*(1.0-n_r*K_Z_n_r/(z_p_n_r+K_Z_n_r));
    const double z_max = (alpha_0/Ca_TRPN_50_Ca_TRPN_max_n-K2)/(alpha_r1+K1+alpha_0/Ca_TRPN_50_Ca_TRPN_max_n);

    const double dz_dt = alpha_0*Ca_b_Ca_TRPN_50_n*(1.0-z)-alpha_r1*z-alpha_r2*z_n_r/(z_n_r+K_Z_n_r);

    // Tension development and crossbridge dynamics.
    const double T_0_max = T_ref*(1.0+beta_0*(zz-1.0));
    const double T_0 = T_0_max*z/z_max;

    const double Q_sum = Q1+Q2+Q3;
    double T = (Q_sum < 0.0 ? T_0*(1.0+a*Q_sum)/(1.0-Q_sum) : T_0*(1.0+(2.0+a)*Q_sum)/(1.0+Q_sum));
    T = min(T,gamma_trpn*T_ref);

    const double dQ1_dt = A1*dlambda_dt-alpha_1*Q1;
    const double dQ2_dt = A2*dlambda_dt-alpha_2*Q2;
    const double dQ3_dt = A3*dlambda_dt-alpha_3*Q3;

    // Troponin C-Calcium binding.
    //const double k_off = k_refoff*(1.0-T/(gamma_trpn*T_ref));
    //updated by Hao, in order to keep k_off above zero. 05/11/2013
    const double k_off = k_refoff*(1.0-T/(gamma_trpn*T_ref));
    const double dCa_b_dt = k_on*Ca_i*(Ca_TRPN_max-Ca_b)-k_off*Ca_b;

    // Update time-dependent variables.
    Ca_i += dt*dCa_i_dt;  Ca_i = max(0.0,Ca_i);
    Ca_b += dt*dCa_b_dt;  Ca_b = max(0.0,Ca_b);
    Ca_b = min(Ca_b, Ca_TRPN_max); //added by Hao 01/11/2013 to ensure Ca_b wont exceed the maximum value, in fact koff>=0 should be enough
    z    += dt*dz_dt;
    Q1   += dt*dQ1_dt;
    Q2   += dt*dQ2_dt;
    Q3   += dt*dQ3_dt;
    return;
}// NHS_euler_step
