// C++ include files that we need
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

// Functions to initialize the library.
#include "libmesh/libmesh.h"
// Basic include files needed for the mesh functionality.
#include "libmesh/mesh.h"
#include "libmesh/point.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/boundary_info.h"

//added by Hao 
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"

#include "libmesh/exodusII_io.h"

using namespace std;
using namespace libMesh;

int main (int argc, char** argv)
{
    LibMeshInit init (argc, argv);

    const unsigned int dim = 3;
    if (dim>3)
	{
		libmesh_error();
	}

    ///reading tetgen mesh manually
	Mesh meshHex(init.comm(), dim);
	meshHex.allow_renumbering(false);
	ifstream ifsNode("heart_real.1.node");
	//ifstream ifsEle("heart_real.1.ele");
   	ifstream ifsEle("heart_real.1.ele.Libmesh_HEX8.dat");
	int NoOfEle, NoOfNode, tempv1, tempv2, tempv3;
	int IDNo, Node1, Node2, Node3, Node4,Node5,Node6,Node7,Node8;
	int side1B,side2B,side3B,side4B,side5B,side6B;
	int side1ID,side2ID,side3ID,side4ID,side5ID,side6ID;
	double CoorX, CoorY, CoorZ;

    ifsNode>>NoOfNode>>tempv1>>tempv2>>tempv3;
	ifsEle>>NoOfEle>>tempv1>>tempv2;

	meshHex.reserve_nodes(NoOfNode);
	meshHex.reserve_elem(NoOfEle);

	///read node inside mesh
    int IDtemp=1;
	while (!ifsNode.eof()& IDtemp<=NoOfNode){
		ifsNode>>IDNo>>CoorX>>CoorY>>CoorZ;
		IDtemp++;
		Point p(CoorX,CoorY,CoorZ);
		meshHex.add_point(p,IDNo-1);
	}
	cout<<IDNo<<'\t'<<CoorX<<'\t'<<CoorY<<'\t'<<CoorZ<<endl;

	///read element inside mesh
	IDtemp = 1;
	while(!ifsEle.eof() & IDtemp<=NoOfEle){
		ifsEle>>IDNo>>Node1>>Node2>>Node3>>Node4>>Node5>>Node6>>Node7>>Node8>>side1B>>side2B>>side3B>>side4B>>side5B>>side6B>>side1ID>>side2ID>>side3ID>>side4ID>>side5ID>>side6ID;
		Elem* elem = meshHex.add_elem(new Hex8);
		elem->set_node(0) = meshHex.node_ptr(Node1-1);
		elem->set_node(1) = meshHex.node_ptr(Node2-1);
		elem->set_node(2) = meshHex.node_ptr(Node3-1);
		elem->set_node(3) = meshHex.node_ptr(Node4-1);
		elem->set_node(4) = meshHex.node_ptr(Node5-1);
		elem->set_node(5) = meshHex.node_ptr(Node6-1);
		elem->set_node(6) = meshHex.node_ptr(Node7-1);
		elem->set_node(7) = meshHex.node_ptr(Node8-1);
		elem->set_id()=IDNo-1;
		elem->subdomain_id()=1;
     		//set up the boundary information 
		if(side1B == 1)
		{   
            meshHex.boundary_info->add_side(elem,0,side1ID);
		}
		if(side2B == 1)
		{	
            meshHex.boundary_info->add_side(elem,1,side2ID);
		}
		if(side3B == 1)
		{	
            meshHex.boundary_info->add_side(elem,2,side3ID);
		}
		if(side4B == 1)
		{	
            meshHex.boundary_info->add_side(elem,3,side4ID);
		}
		if(side5B == 1)
		{   
            meshHex.boundary_info->add_side(elem,4,side5ID);
		}
		if(side6B == 1)
		{	
            meshHex.boundary_info->add_side(elem,5,side6ID);
		}

	
		IDtemp++;
        }
 	
	meshHex.prepare_for_use();
	meshHex.print_info();	
	
///output the mesh file
	#ifdef LIBMESH_HAVE_EXODUS_API
  		ExodusII_IO meshHex_writer(meshHex);
 		meshHex_writer.write("heart_exodusII_hex.e");
  	//	//mesh_writer.write(argv[2]);
	#endif
//	TetGenIO meshTet_writer(meshTet);
//	meshTet.prepare_for_use();
//	meshTet_writer.write("heart_new");

	meshHex.write("heart_real_hex.xda");

    return 0;
}