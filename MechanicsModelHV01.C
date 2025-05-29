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

// APPLICATION INCLUDES
#include <ActiveContraction.h>
#include <MechanicsModel.h>

// STATIC VARIABLES
int MechanicsModel::f0_system_num, MechanicsModel::s0_system_num;
bool MechanicsModel::enable_active_tension, MechanicsModel::normalize_stress;
double MechanicsModel::T_scale, MechanicsModel::beta_s;
double MechanicsModel::I1_dev_max;
double MechanicsModel::I1_dev_min;
double MechanicsModel::J_dev_max;
double MechanicsModel::J_dev_min;
double MechanicsModel::I1_dil_max;
double MechanicsModel::I1_dil_min;
double MechanicsModel::J_dil_max;
double MechanicsModel::J_dil_min;
namespace // private namespace
{
#if 0
// Material parameters (from Ogden & Holzapfel).
static const double a   =  5.900e+2;  // dyne/cm^2
static const double b   =  8.023   ;  // nondimensional
static const double af  = 1.8472e+5;  // dyne/cm^2
static const double bf  = 16.026   ;  // nondimensional
static const double as  =  2.481e+4;  // dyne/cm^2
static const double bs  = 11.120   ;  // nondimensional
static const double afs =  2.150e+3;  // dyne/cm^2
static const double bfs = 11.436   ;  // nondimensional
#else
//new optimized value 31/01/2013
//static const double a   = 2.362e+3;  // dyne/cm^2
//static const double b   = 5.081    ;  // nondimensional
//static const double af  = 1.46625e+4;  // dyne/cm^2
//static const double bf  = 4.15   ;  // nondimensional
//static const double as  = 8.625e+3;  // dyne/cm^2
//static const double bs  = 1.6  ;  // nondimensional
//static const double afs = 2.953e+3;  // dyne/cm^2
//static const double bfs = 1.3     ;  // nondimensional
//opt for HV21 using abaqus
//0.224487,	1.621500,	2.426717,	1.826862,	0.556238,	0.774678,	0.390516,	1.695000
static const double a   = 2244.87;  // dyne/cm^2
static const double b   = 1.6215;  // nondimensional
static const double af  = 2.4267e4;  // dyne/cm^2
static const double bf  = 1.8268;  // nondimensional
static const double as  = 5562.38;  // dyne/cm^2
static const double bs  = 0.7746;  // nondimensional
static const double afs = 3905.16;  // dyne/cm^2
static const double bfs = 1.695;  // nondimensional
#endif
}

// CLASS IMPLEMENTATION

void
MechanicsModel::PK1_dev_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /* X */,  // current   location
    const libMesh::Point& /* s */,  // reference location
    Elem* const elem,
    const vector<const vector<double>*>&  var_data,
    const vector<const vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
    double /* data_time */,
    void* /* ctx */)
{
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF,NDIM);

    // Isotropic contribution.
    const double I1 = CC.tr();
    const double J = FF.det();

    I1_dev_max = max(I1_dev_max, I1);
    I1_dev_min = min(I1_dev_min, I1);

    J_dev_max = max(J_dev_max, J);
    J_dev_min = min(J_dev_min, J);

    //PP = a*exp(b*(I1-3.0))*FF;
    PP = a*exp(b*(I1-3.0))*FF;
	//-a*exp(b*(I1-3.0))*FF_inv_trans;

   const std::vector<double>& f0_var = *var_data[0];
   const std::vector<double>& s0_var = *var_data[1];


   const VectorValue<double> f0(f0_var[0], f0_var[1], f0_var[2]);
   const VectorValue<double> s0(s0_var[0], s0_var[1], s0_var[2]);


    // Fiber contribution.
    //NumericVector<double>* f0_vec = system_data[0];
    //VectorValue<double> f0;
    //for (unsigned int d = 0; d < NDIM; ++d)
    //{
    //    f0(d) = (*f0_vec)(elem->dof_number(f0_system_num,d,0));
    //}
    const double I4f = f0*(CC*f0);
    if (I4f > 1.0)
    {
        PP += 2.0*af*(I4f-1.0)*exp(bf*(I4f-1.0)*(I4f-1.0))*FF*outer_product(f0,f0);
    }

    // Sheet contribution.
    //NumericVector<double>* s0_vec = system_data[1];
    //VectorValue<double> s0;
    //for (unsigned int d = 0; d < NDIM; ++d)
    //{
    //    s0(d) = (*s0_vec)(elem->dof_number(s0_system_num,d,0));
    //}
    const double I4s = s0*(CC*s0);
    if (I4s > 1.0)
    {
        PP += 2.0*as*(I4s-1.0)*exp(bs*(I4s-1.0)*(I4s-1.0))*FF*outer_product(s0,s0);
    }

    // Fiber-sheet contribution.
    const double I8fs = f0*(CC*s0);
    PP += afs*I8fs*exp(bfs*I8fs*I8fs)*FF*(outer_product(f0,s0)+outer_product(s0,f0));

    // Active tension contribution.
    //if (enable_active_tension)
    //{
    //    NumericVector<double>* T_vec = system_data[2];
    //    const double T = T_scale*(*T_vec)(elem->dof_number(ActiveContraction::T_system_num,0,0))*1.0e4;  // NOTE: T is stored in kPa; must be converted to dyne/cm^2
    //    if (T > 0.0)
    //    {
    //        PP += J*T*FF*outer_product(f0,f0);
    //    }
    //}
    return;
}// PK1_dev_stress_function

void
MechanicsModel::PK1_dil_stress_function(
    TensorValue<double>& PP,
    const TensorValue<double>& FF,
    const libMesh::Point& /* X */,  // current   location
    const libMesh::Point& /* s */,  // reference location
    Elem* const elem,
    const vector<const vector<double>*>& var_data,
    const vector<const vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
    double /* data_time */,
    void* /* ctx */)
{
    const TensorValue<double> CC = FF.transpose() * FF;
    const TensorValue<double> FF_inv_trans = tensor_inverse_transpose(FF,NDIM);
    const double I1 = CC.tr();
    const double J = FF.det();

    I1_dil_max = max(I1_dil_max, I1);
    I1_dil_min = min(I1_dil_min, I1);

    J_dil_max = max(J_dil_max, J);
    J_dil_min = min(J_dil_min, J);

    PP = -a*exp(b*(I1-3.0))*FF_inv_trans;

    const std::vector<double>& f0_var = *var_data[0];
    const std::vector<double>& s0_var = *var_data[1];
    const std::vector<double>& acT_var = *var_data[2];

    const VectorValue<double> f0(f0_var[0], f0_var[1], f0_var[2]);
    const VectorValue<double> s0(s0_var[0], s0_var[1], s0_var[2]);
    const double acT = acT_var[0];

	if (enable_active_tension)
    {
        //NumericVector<double>* T_vec = system_data[2];
        const double T = T_scale*acT*1.0e4;  // NOTE: T is stored in kPa; must be converted to dyne/cm^2
        if (T > 0.0)
        {
            PP += J*T*FF*outer_product(f0,f0);
        }
    }

    if (!MathUtilities<double>::equalEps(beta_s, 0.0))
    {
        PP += (beta_s*log(CC.det()))*FF_inv_trans; // volume constraint
    }
    return;
}// PK1_dil_stress_function
