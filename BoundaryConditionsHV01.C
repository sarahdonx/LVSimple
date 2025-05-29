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

// APPLICATION INCLUDES
#include <BoundaryConditions.h>

// LIBMESH INCLUDES
#include <libmesh/boundary_info.h>
#include <libmesh/point.h>

// STATIC VARIABLES
double BoundaryConditions::P_load, BoundaryConditions::t_load, BoundaryConditions::kappa;
double BoundaryConditions::P_load_es;
BoundaryInfo* BoundaryConditions::boundary_info;
double BoundaryConditions::P_current_loading;
double BoundaryConditions::t_end_diastole;

double BoundaryConditions::p_loading_bcast_from_root;
double BoundaryConditions::LV_volume;
std::vector<int>  BoundaryConditions::LV_endo_points_list;
int BoundaryConditions::LV_NoOfEndoNode;
std::vector< std::vector<double>  > BoundaryConditions::LV_endo_points;


// CLASS IMPLEMENTATION

double
BoundaryConditions::loading_pressure(
    double time)
{
    double P = 0.0;
    double Pfactor = P_load_es/P_load;
    //double t_constant = 0.0;
    //return P_load*(time < t_load ? time/t_load : 1.0);

    if (time <= t_load)
    {
        P = P_load*time/(t_load);
    }
    else if (time>t_load && time<=t_end_diastole)
    {
        P = P_load;
    }
    else if (time>t_end_diastole && time<t_end_diastole + 0.2)
    {
        P = P_load*(Pfactor-1)*(time<=t_end_diastole + 0.2? (time-t_end_diastole)/(0.2): 1.0 ) + P_load;
    }
    else
    {
        P = Pfactor*P_load;
    }

    P_current_loading = P; //to output current P loading
    return P;
}// loading_pressure

void
BoundaryConditions::loading_force_function(
    double& P,
    const libMesh::VectorValue<double>& n,
  	const libMesh::VectorValue<double>& N,
    const TensorValue<double>& /*FF*/,
    const libMesh::Point& /*x*/,
    const libMesh::Point& /*X*/,
    Elem* const elem,
    unsigned short int side,
    const vector<const vector<double>*>& /*system_data*/,
    const vector<const vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
    double time,
    void* /*ctx*/)
{
    // 4096 ===>  epicardium
    // 4097 ===> endocardium
    // 4098 ===>        base
    // const vector<short int>& bdry_ids = boundary_info->boundary_ids(elem, side);
    //const vector<short int>& bdry_ids = boundary_info->boundary_ids(elem, side);
     vector<short int> bdry_ids;
     boundary_info->boundary_ids(elem, side, bdry_ids);

    if (find(bdry_ids.begin(),bdry_ids.end(),4097) != bdry_ids.end())
    {
        P = loading_pressure(time)*1333.2239;
    }
    else if (find(bdry_ids.begin(),bdry_ids.end(),5090) != bdry_ids.end() ||
             find(bdry_ids.begin(),bdry_ids.end(),5091) != bdry_ids.end())
    {  //for right ventricle
		P = loading_pressure(time)*1333.2239/3.0; //only one third of the left ventricular pressure
	  }
    else
    {
        P = 0.0;
    }
    return;
}// loading_force_function

void
BoundaryConditions::tether_force_function(
    VectorValue<double>& F,
    const libMesh::VectorValue<double>& n,
  	const libMesh::VectorValue<double>& N,
    const TensorValue<double>& /*FF*/,
    const libMesh::Point& x,
    const libMesh::Point& X,
    Elem* const elem,
    const unsigned short int side,
    const vector<const vector<double>*>& /*system_data*/,
    const vector<const vector<libMesh::VectorValue<double> >*>& /*system_grad_var_data*/,
    double /*time*/,
    void* /*ctx*/)
{
    F.zero();

    // 4096 ===>  epicardium
    // 4097 ===> endocardium
    // 4098 ===>        base
    //const vector<short int>& bdry_ids = boundary_info->boundary_ids(elem, side);
    vector<short int> bdry_ids;
    boundary_info->boundary_ids(elem, side, bdry_ids);
    if (find(bdry_ids.begin(),bdry_ids.end(),4098) != bdry_ids.end())
    {
        // approximately constrain motion to be in the radial direction at base
        VectorValue<double> r = x - libMesh::Point(0.0,0.0,x(2));  // radial displacement in current   configuration
        VectorValue<double> R = X - libMesh::Point(0.0,0.0,X(2));  // radial displacement in reference configuration
        VectorValue<double> r_perp = (R*r)*R/R.norm_sq() - r;
        F = kappa*r_perp;
        F(2) += kappa*(X(2)-x(2));
    }
    return;
}// tether_force_function



void
BoundaryConditions::readingPoints(MeshBase& mesh)
{

	   std::ifstream ifsendo("endoList.point");
	   int NoOfEndoNode, nodeEndoID;

	   ifsendo>>NoOfEndoNode;
	   LV_NoOfEndoNode = NoOfEndoNode;



	   LV_endo_points_list.resize(NoOfEndoNode);
	   LV_endo_points.resize(LV_NoOfEndoNode);
	   for (unsigned int i = 0; i< LV_NoOfEndoNode; i++)
	   {
		      LV_endo_points[i].resize(3);
	   }

	   unsigned int IDtemp=1; //reused from the initial defintion
	   unsigned int pIndex = 0;

	   //initialize end_points
	   while (!ifsendo.eof()& IDtemp<=NoOfEndoNode){
			ifsendo>>nodeEndoID;
			IDtemp++;
			nodeEndoID = nodeEndoID - 1 ; //start from 0

		    LV_endo_points_list[pIndex]=nodeEndoID;

			pIndex = pIndex + 1;

		}


		printf("processor %d read %d points\n", mesh.processor_id(), pIndex);

	   return ;
}


void
BoundaryConditions::readingPointsGeneral(MeshBase& mesh,
                                    std::vector<int>& points_list,
                                    std::vector< std::vector<double> >& points,
                                    int& NoOfPoints,
                                    std::string file_name)
{
	std::ifstream ifsendo(file_name.c_str());
	int NoOfEndoNode, nodeEndoID;

	ifsendo>>NoOfEndoNode;
	NoOfPoints = NoOfEndoNode;



	   points_list.resize(NoOfPoints);
	   points.resize(NoOfPoints);
	   for (unsigned int i = 0; i< NoOfPoints; i++)
	   {
		      points[i].resize(3);
	   }

	   unsigned int IDtemp=1; //reused from the initial defintion
	   unsigned int pIndex = 0;

	   //initialize end_points
	   while (!ifsendo.eof()& IDtemp<=NoOfPoints){
			ifsendo>>nodeEndoID;
			IDtemp++;
			nodeEndoID = nodeEndoID - 1 ; //start from 0

		    points_list[pIndex]=nodeEndoID;
            //cout << points_list[pIndex]<< "\n";

			pIndex = pIndex + 1;

		}


		printf("processor %d read %d points\n", mesh.processor_id(), pIndex);
	return ;
}


void
BoundaryConditions::updatePointsPosition(EquationSystems* equation_systems)
{
	       const MeshBase& mesh = equation_systems->get_mesh();
	       const unsigned int dim = mesh.mesh_dimension();
		   System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);



		   //libMesh::Point  node_ref;
		   //libMesh::Node    node_ref;

		   for (unsigned int i = 0; i < LV_NoOfEndoNode; i++)
		   {
			   int nodeEndoID = LV_endo_points_list[i]; //start from 0

			   //node_ref = mesh.point(nodeEndoID);
			   const libMesh::Node& node_ref = mesh.node_ref(nodeEndoID);

			   LV_endo_points[i][0] = node_ref(0);
			   LV_endo_points[i][1] = node_ref(1);
			   LV_endo_points[i][2] = node_ref(2);
		   }

		   return ;
		   //double volume = tetVolumeCalculationByPoints(endo_points, BoundaryConditions::LV_NoOfEndoNode);
}


void BoundaryConditions::updatePointsPositionGeneral(EquationSystems* equation_systems,
                                                    std::vector<int> & points_list,
                                     std::vector< std::vector<double> >& points_coor,
                                                    int NoOfPoints)
{



	       const MeshBase& mesh = equation_systems->get_mesh();
	       const unsigned int dim = mesh.mesh_dimension();

		   System& X_system = equation_systems->get_system<System>(IBFEMethod::COORDS_SYSTEM_NAME);
		   const DofMap& X_dof_map = X_system.get_dof_map();
		   std::vector<std::vector<unsigned int> > X_dof_indices(NDIM);

		   X_system.solution->localize(*X_system.current_local_solution);
           NumericVector<double>& X_data = *(X_system.current_local_solution);
           X_data.close();

           //copy data to one processor, default is 0 processor
           std::vector<double> X_data_vec;
	       X_system.solution->localize_to_one(X_data_vec);
           const unsigned int X_sys_num = X_system.number();

      // print some information for debug
      //equation_systems->print_info();
      //mesh.print_info();

		    //printf("working on the updating, processor: %d\n", X_system.processor_id());
		    MPI_Barrier(X_system.comm().get());
            if (0 == X_system.processor_id())
      {
			   //printf("updating node position for vol cal processor 0\n");
			   for (unsigned int i = 0; i < NoOfPoints; i++)
			   {
				   unsigned int nodeEndoID = points_list[i]; //start from 0
           //printf("update node: %d using system %d\n", nodeEndoID,X_sys_num);
				    //node_ref = mesh.point(nodeEndoID);
				    //libMesh::Point  node_ref;
		        //libMesh::Node&    node_ref;
            //printf("access node : %d\n", nodeEndoID);
		        const libMesh::Node& node_ref = mesh.node_ref(nodeEndoID); //const is needed otherwise it just abort
            //printf("obtain node ref : %d\n", nodeEndoID);


				   const int x_dof_index = node_ref.dof_number(X_sys_num, 0, 0);
				   const int y_dof_index = node_ref.dof_number(X_sys_num, 1, 0);
				   const int z_dof_index = node_ref.dof_number(X_sys_num, 2, 0);

				   points_coor[i][0] = X_data_vec[x_dof_index];
				   points_coor[i][1] = X_data_vec[y_dof_index];
				   points_coor[i][2] = X_data_vec[z_dof_index];
			   }
          printf("updating node position for vol cal processor 0 done\n");
		   }

		   //double volume = tetVolumeCalculationByPoints(endo_points, BoundaryConditions::LV_NoOfEndoNode);
		   return;
}
