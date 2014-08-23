//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.85. June 9, 2008.
//LIC// 
//LIC// Copyright (C) 2006-2008 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================
//Driver for 3D rectangular driven cavity

#include <fenv.h>

//Generic includes
#include "generic.h"
#include "navier_stokes.h"

#include "meshes/simple_cubic_mesh.h"
//#include "meshes/simple_cubic_tet_mesh.h"

// Get the mesh
#include "meshes/tetgen_mesh.h"
#include "meshes/brick_from_tet_mesh.h"

// My own header
#include "./../rayheader.h"
using namespace std;

using namespace oomph;

// Alias the namespace for convenience.
//namespace NSPP = NavierStokesProblemParameters;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

namespace Global_Parameters
{
  unsigned Noel = 4;

  double Length = 1.0;

  double Rey = 200.0;

  DocLinearSolverInfo* Doc_linear_solver_info_pt;


  inline double get_prescribed_inflow_for_quarter(const double& y, 
      const double& z)
  {
#ifdef PARANOID
    // Quick check that we are in the correct range of
    // y and z coordinates.
    if( !((y > 0.5)&&(z > 0.5)) )
    {
      std::ostringstream error_message;
      error_message << "Prescribed inflow: incorrect range of y, z.\n"
        << "y = " << y << ", z = " << z << std::endl;

      throw OomphLibError(error_message.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    return (y-0.5)*(1.0-y)*(z-0.5)*(1.0-z);
  }

  inline double get_prescribed_inflow_for_quarter(const double& t,
      const double& y,
      const double& z)
  {
    // For the velocity profile in the x direction.
    // 1) First form the parabolic profile
   
    const double ux_scaling = -cos(MathematicalConstants::Pi*t)/2.0 + 0.51;
    return get_prescribed_inflow_for_quarter(y,z) * ux_scaling;
  } 


}




//==start_of_problem_class============================================
/// Test problem for Fp/PCD preconditioner
//====================================================================
template<class ELEMENT>
class CubeProblem : public Problem
{

public:


 /// Constructor
 CubeProblem();

 
 /// Destructor: Cleanup
 virtual ~CubeProblem()
  {
//    GenericProblemSetup::clean_up_solver_memory();
//    LPH::clean_up_memory();

    delete Bulk_mesh_pt;
//   delete Prec_pt;
//   delete P_matrix_preconditioner_pt;
//   delete F_matrix_preconditioner_pt;
  }

//////////////////////////////

// void actions_before_implicit_timestep()
// {
//   Doc_linear_solver_info_pt->clear_current_time_step();
//
//   {
//     // Inflow in upper half of inflow boundary
//     const unsigned ibound=Inflow_boundary; 
//     const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
//     for (unsigned inod=0;inod<num_nod;inod++)
//     {
//       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
//       //      const double x=nod_pt->x(0);
//       const double y=nod_pt->x(1);
//       const double z=nod_pt->x(2);
//
//       // Only do the inflow for a quarter of the boundary
//       if((y > 0.5) && (z > 0.5))
//       {
//
//         const double time=time_pt()->time();
//
//         double ux = 0.0;
//
//         CL::get_prescribed_inflow(time,y,z,ux);
//
//         nod_pt->set_value(0,ux);
//         nod_pt->set_value(1,0.0);
//         nod_pt->set_value(2,0.0);
//       }
//     }
//   }
// } // end of actions_before_implicit_timestep



////////////////////////////////


 void actions_before_distribute()
 {

//   if(NSPP::Distribute_problem)
//   {
//     GenericProblemSetup::delete_flux_elements(Surface_mesh_pt);
//     rebuild_global_mesh();
//   }
 }

 void actions_after_distribute()
 {
//   if(NSPP::Distribute_problem)
//   {
//     create_parall_outflow_lagrange_elements(Outflow_boundary,
//         Bulk_mesh_pt,
//         Surface_mesh_pt);
//     rebuild_global_mesh();
//   }
 }


 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   
   // Now set the first pressure dof in the first element to 0.0
//   if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt
 

 /// Update the after solve (empty)
 void actions_after_newton_solve()
 {


 }

 /// Update the problem specs before solve. 
 void actions_before_newton_solve()
 {
   // Start a new "time step"
   Doc_linear_solver_info_pt->setup_new_time_step();

    // Inflow in upper half of inflow boundary
    const unsigned ibound=Inflow_boundary; 
    const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      const double y=nod_pt->x(1);
      const double z=nod_pt->x(2);

      if( (y > 0.5) && (z > 0.5) )
      {

      const double ux 
        = Global_Parameters::get_prescribed_inflow_for_quarter(y,z);

      nod_pt->set_value(0,ux);
      nod_pt->set_value(1,0.0);
      nod_pt->set_value(2,0.0);
      }
     }


  // Initialise counter for iterations
//  Global_Variables::Iterations.clear();
//  Global_Variables::Linear_solver_time.clear();
  
 } // end_of_actions_before_newton_solve
 void actions_after_newton_step()
 {
//   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
//   {
//     NSPP::doc_iter_times(this,Doc_linear_solver_info_pt);
//   }
 }
 
 /// Global error norm for adaptive time-stepping
 double global_temporal_error_norm();


 /// Pointer to the "bulk" mesh
 Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}


 void doc_solution(Mesh* bulk_mesh_pt)
 {
  std::ofstream some_file;
  std::stringstream filename;

  filename << "tmp_soln.dat";

  // Number of plot points
  const unsigned npts=5;

  // Output solution
  some_file.open(filename.str().c_str());
  bulk_mesh_pt->output(some_file,npts);
  some_file.close();
 }

//private:

 /// Solver
 IterativeLinearSolver* Solver_pt;

 /// Solver
 Preconditioner* Prec_pt;

 /// Inexact solver for P block
 Preconditioner* P_matrix_preconditioner_pt;

 /// Inexact solver for F block
 Preconditioner* F_matrix_preconditioner_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;



 unsigned Left_boundary;
 unsigned Right_boundary;
 unsigned Front_boundary;
 unsigned Back_boundary;
 unsigned Bottom_boundary;
 unsigned Top_boundary;

 /// ID of inflow boundary
 unsigned Inflow_boundary;

 /// ID of outflow boundary
 unsigned Outflow_boundary;

 /// Pointer to the "bulk" mesh
 Mesh* Bulk_mesh_pt;
 

}; // end_of_problem_class

//==start_of_constructor==================================================
/// Constructor for cube problem
//        
//        y
//        |
//        |
//        |
//        |____________ x
//        /
//       /
//      /
//     /
//    z
//        --------------
//       /|            /      Left = 4 (Inflow)
//      / |           / |     Right = 2 (Outflow)
//     /  |          /  |     Front = 5
//    /   |         /   |     Back = 0
//    --------------    |     Bottom = 1
//    |   |--------|----|     Top = 3
//    |   /        |   /
//    |  /         |  /
//    | /          | /
//    |/           |/
//    --------------
// 
//========================================================================
  template<class ELEMENT>
CubeProblem<ELEMENT>::CubeProblem()
{ 
  // Assign boundary IDs defined in rayheader.h
  Left_boundary = 4;
  Right_boundary = 2;
  Front_boundary = 5;
  Back_boundary = 0;
  Bottom_boundary = 1;
  Top_boundary = 3;

  Inflow_boundary=Left_boundary;
  Outflow_boundary=Right_boundary;
  Doc_linear_solver_info_pt = Global_Parameters::Doc_linear_solver_info_pt;

//  if(NSPP::Time_type == NSPP::Time_type_STEADY)
//  {
//    oomph_info << "Doing steady state.\n"
//      << "No time stepper added." << std::endl;
//  }
//  else if(NSPP::Time_type == NSPP::Time_type_ADAPT)
//  {
//    oomph_info << "Adding adaptive time stepper" << std::endl; 
//    add_time_stepper_pt(new BDF<2>(true));
//  }
//  else if(NSPP::Time_type == NSPP::Time_type_FIXED)
//  {
//    oomph_info << "Adding non-adaptive time stepper" << std::endl; 
//    add_time_stepper_pt(new BDF<2>);
//  }
//  else
//  {
//    std::ostringstream err_msg;
//    err_msg << "Time stepper for Time_type: "
//      << NSPP::Time_type << std::endl;
//
//    throw OomphLibError(err_msg.str(),
//        OOMPH_CURRENT_FUNCTION,
//        OOMPH_EXCEPTION_LOCATION);
//  }


  // Setup mesh
  const unsigned noel = Global_Parameters::Noel;
  const double length = Global_Parameters::Length;

//  if(NSPP::Time_type == NSPP::Time_type_STEADY)
  {
    Bulk_mesh_pt = 
      new SimpleCubicMesh<ELEMENT>(noel,noel,noel,
          length, length, length);
  }
//  else if((NSPP::Time_type == NSPP::Time_type_ADAPT)
//      ||NSPP::Time_type == NSPP::Time_type_FIXED)
//  {
//    Bulk_mesh_pt = 
//      new SimpleCubicMesh<ELEMENT>(noel,noel,noel, 
//          length, length, length,
//          time_stepper_pt());
//  }
//  else
//  {
//    std::ostringstream err_msg;
//    err_msg << "No mesh set up for Time_type: "
//      << NSPP::Time_type << std::endl;
//
//    throw OomphLibError(err_msg.str(),
//        OOMPH_CURRENT_FUNCTION,
//        OOMPH_EXCEPTION_LOCATION);
//  }


  // Add the two sub meshes to the problem
  add_sub_mesh(Bulk_mesh_pt);

  // Combine all submeshes into a single Mesh
  build_global_mesh();


  // Set the boundary conditions for this problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here. 
  unsigned num_bound = Bulk_mesh_pt->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
      Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      // Loop over values (u, v and w velocities)
       nod_pt->pin(0);
       nod_pt->pin(1);
       nod_pt->pin(2);

       nod_pt->set_value(0,0.0);
       nod_pt->set_value(1,0.0);
       nod_pt->set_value(2,0.0);
    }
  } // end loop over boundaries!


  {
    unsigned ibound=Outflow_boundary;
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      // Only free if node is ONLY on a single boundary
      std::set<unsigned>* bnd_pt=0;
      nod_pt->get_boundaries_pt(bnd_pt);
      if (bnd_pt!=0)
      {
        if (bnd_pt->size()<2)
        {
          {
            const double y = nod_pt->x(1);
            const double z = nod_pt->x(2);
            
            if ((y<0.5)&&(z<0.5))
            {
              nod_pt->unpin(0);
            }
          }
        }
      }
    }
  }

  // Complete the build of all elements so they are fully functional

  //Find number of elements in mesh
  unsigned n_element = Bulk_mesh_pt->nelement();

  // Loop over the elements to set up element-specific 
  // things that cannot be handled by constructor
  for(unsigned e=0;e<n_element;e++)
  {
    // Upcast from GeneralisedElement to the present element
    NavierStokesEquations<3>* el_pt = 
      dynamic_cast<NavierStokesEquations<3>*>(Bulk_mesh_pt->element_pt(e));

    //Set the Reynolds number
    el_pt->re_pt() = &Global_Parameters::Rey;
    el_pt->re_st_pt() = &Global_Parameters::Rey;
  } // end loop over elements


  // Now set the first pressure value in element 0 to 0.0
  // if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);

  // Setup equation numbering scheme
  oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 



//  
//  if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
//  {
//    Vector<Mesh*> mesh_pt;
//    if(NSPP::Prob_id == CL::PID_CU_PO_QUARTER)
//    {
//      mesh_pt.resize(2,0);
//      mesh_pt[0] = Bulk_mesh_pt;
//      mesh_pt[1] = Surface_mesh_pt;
//    }
//
//    LPH::Mesh_pt = mesh_pt;
//    LPH::Problem_pt = this;
//    Prec_pt = LPH::get_preconditioner();
//  }
//
// const double solver_tol = 1.0e-6;
// const double newton_tol = 1.0e-6;
// GenericProblemSetup::setup_solver(NSPP::Max_solver_iteration,
//                                   solver_tol,newton_tol,
//                                   NSPP::Solver_type,this,Prec_pt);
} // end_of_constructor


  template<class ELEMENT>
double CubeProblem<ELEMENT>::global_temporal_error_norm()
{
  return GenericProblemSetup::global_temporal_error_norm(this,
      3,Bulk_mesh_pt);
} // end of global_temporal_error_norm




//============start_of_create_traction_elements==========================
/// Create Navier-Stokes traction elements on outflow boundary
//=======================================================================
//template<class ELEMENT>
//void CubeProblem::create_traction_elements()
//{
//
// unsigned b=Outflow_boundary;
//
// // How many bulk elements are adjacent to boundary b?
// unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
//
// // Loop over the bulk elements adjacent to boundary b?
// for(unsigned e=0;e<n_element;e++)
//  {
//   // Get pointer to the bulk element that is adjacent to boundary b
//   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
//    Bulk_mesh_pt->boundary_element_pt(b,e));
//   
//   //What is the index of the face of element e along boundary b
//   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
//   
//   // Build the corresponding prescribed-flux element
//   NavierStokesTractionElement<ELEMENT>* flux_element_pt = new 
//      NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
//
//   //Add the prescribed-flux element to the surface mesh
//   Surface_mesh_pt->add_element_pt(flux_element_pt);
//   
//   // Set the pointer to the prescribed traction function
//   flux_element_pt->traction_fct_pt() = &Global_Variables::prescribed_traction;
//   
//  } //end of loop over bulk elements adjacent to boundary b
//
// // Now rebuild the global mesh
// rebuild_global_mesh();
//
// // Reassign equation numbers
// oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
//
//} // end of create_traction_elements


//Template<class ELEMENT>
//Void CubeProblem<ELEMENT>::create_inflow_traction_elements(
//    const unsigned &b, 
//    Mesh* const &bulk_mesh_pt, 
//    Mesh* const &surface_mesh_pt)
//{
// // How many bulk elements are adjacent to boundary b?
// unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
//
// // Loop over the bulk elements adjacent to boundary b?   // Set the boundary conditions for this problem: All nodes are
//   // free by default -- just pin the ones that have Dirichlet conditions
//   // here. 
//   unsigned num_bound = Bulk_mesh_pt->nboundary();
//   for(unsigned ibound=0;ibound<num_bound;ibound++)
//    {
//     unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
//     for (unsigned inod=0;inod<num_nod;inod++)
//      {
//       // Loop over values (u, v and w velocities)
//       for (unsigned i=0;i<3;i++)
//        {
//         Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(i); 
//        }
//      }
//    } // end loop over boundaries
//
// 
//
//   // OUTFLOW ONLY, for inflow, check out the before solve
//// if (Problem_id==Global_Variables::Through_flow)
//  {
//   unsigned ibound=Outflow_boundary;
//   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
//   for (unsigned inod=0;inod<num_nod;inod++)
//    {
//     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
//     // Only free if node is ONLY on a single boundary
//     std::set<unsigned>* bnd_pt=0;
//     nod_pt->get_boundaries_pt(bnd_pt);
//     if (bnd_pt!=0)
//      {
//       if (bnd_pt->size()<2)
//        {
//         if (!(nod_pt->is_on_boundary(0)))
//          {
//           if ((nod_pt->x(1)<0.5)&&nod_pt->x(2)<0.5) nod_pt->unpin(0);
//          }
//        }
//      }
//    }
//  }
//
// // Complete the build of all elements so they are fully functional
//
// //Find number of elements in mesh
// unsigned n_element = Bulk_mesh_pt->nelement();
//
// // Loop over the elements to set up element-specific 
// // things that cannot be handled by constructor
// for(unsigned e=0;e<n_element;e++)
//  {
//   // Upcast from GeneralisedElement to the present element
//   NavierStokesEquations<3>* el_pt = 
//    dynamic_cast<NavierStokesEquations<3>*>(Bulk_mesh_pt->element_pt(e));
//   
//   //Set the Reynolds number
//   el_pt->re_pt() = &Global_Variables::Re;
//  } // end loop over elements
// 
// 
//
// // Now set the first pressure value in element 0 to 0.0
//// if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);
//
// // Setup equation numbering scheme
// oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
// for(unsigned e=0;e<n_element;e++)
//  {
//   // Get pointer to the bulk element that is adjacent to boundary b
//   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
//    Bulk_mesh_pt->boundary_element_pt(b,e));
//
//   // Loop through all of the nodes on this element and find and out if
//   // all the y and z coordinates with within [0.5]^2
//   const unsigned nbulk_nod = bulk_elem_pt->nnode();
//   bool within_inflow = true;
//
//   for(unsigned nod_i = 0; (nod_i < nbulk_nod) && within_inflow; nod_i++)
//   {
//     Node* bulk_nod_pt = bulk_elem_pt->node_pt(nod_i);
//     const double y = bulk_nod_pt->x(1);
//     const double z = bulk_nod_pt->x(2);
//     if((y <= 0.5) || (z <= 0.5))
//     {
//       within_inflow = false;
//     }
//   }
//   
//   if(within_inflow)
//   {
//     //What is the index of the face of element e along boundary b
//     int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);
//
//     // Build the corresponding prescribed-flux element
//     NavierStokesTractionElement<ELEMENT>* flux_element_pt = new 
//       NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,face_index);
//
//     //Add the prescribed-flux element to the surface mesh
//     Surface_mesh_pt->add_element_pt(flux_element_pt);
//
//     // Set the pointer to the prescribed traction function
//     flux_element_pt->traction_fct_pt() 
//       = &Global_Variables::inflow_prescribed_traction;
//   }
//   
//  } //end of loop over bulk elements adjacent to boundary b
//
// // Now rebuild the global mesh
// rebuild_global_mesh();
//
// // Reassign equation numbers
// oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
//
//} // end of create_traction_elements






std::string create_label()
{
  // We want the unique problem label, then the generic problem label, then
  // preconditioner used.
  //
  // Because the unique problem label contains a label for the problem
  // and parameters such as angle/noel, we want the unique problem 
  // identifier to be first, then the parameters last, with the generic
  // problem stuff in between.
  //
  // i.e.
  // SqPo + NSPP::label + LPH::label Ang Noel.
  
  std::string label = "tmp_label";
//  std::string label = CL::prob_str()
//                      + NSPP::create_label() 
//                      + LPH::create_label() 
//                      + CL::ang_deg_str() + CL::noel_str();
  return label;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_main======================================================
/// Driver for Fp preconditioner
//=====================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Problem dimension.
  const unsigned dim = 3;

  // Set up doc info - used to store information on solver and iteration time.
  DocLinearSolverInfo doc_linear_solver_info;

  Global_Parameters::Doc_linear_solver_info_pt = &doc_linear_solver_info;
  // Again, pass this to the NSPP and LPH
  
//  NSPP::Doc_linear_solver_info_pt = &doc_linear_solver_info;
//  LPH::Doc_linear_solver_info_pt = &doc_linear_solver_info;

  // Set the Label_pt
//  LPH::Label_str_pt = &NSPP::Label_str;
//  LPH::Vis_pt = &NSPP::Vis;
//  CL::Prob_id_pt = &NSPP::Prob_id;

//  NSPP::Time_start = 0.0;
//  NSPP::Time_end = 1.0; 


//  RaySpace::Use_tetgen_mesh = true;

  // Store commandline arguments
//  CommandLineArgs::setup(argc,argv);

//  NSPP::setup_commandline_flags();
//  LPH::setup_commandline_flags();
//  CL::setup_commandline_flags(); 

  // Parse the above flags.
//  CommandLineArgs::parse_and_assign();
//  CommandLineArgs::doc_specified_flags();

  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////

  // dim = 3
//  NSPP::generic_problem_setup(dim);
//  LPH::generic_setup();
//  CL::generic_setup(); 




  //////////////////////////////////////  

//  if(RaySpace::Use_tetgen_mesh)
//  {
//   // Build the problem 
//  CubeProblem <QTaylorHoodElement<3> >problem;
//
//
//  // Solve the problem 
//  //              problem.newton_solve();
//
//  if(NSPP::Distribute_problem)
//  {
//    problem.distribute();
//  }
//
//  NSPP::Label_str = create_label();
//
//  time_t rawtime;
//  time(&rawtime);
//
//  std::cout << "RAYDOING: "
//    << NSPP::Label_str
//    << " on " << ctime(&rawtime) << std::endl;
//
//  problem.unsteady_run(); 
//  
//  }
//  else
  {
  // Build the problem 
  CubeProblem <QTaylorHoodElement<3> >problem;


  // Solve the problem 
  //              problem.newton_solve();

//  if(NSPP::Distribute_problem)
  {
    problem.distribute();
  }

  std::string label = create_label();
//  NSPP::Label_str = create_label();

  time_t rawtime;
  time(&rawtime);

  std::cout << "RAYDOING: "
    << label
    << " on " << ctime(&rawtime) << std::endl;

  problem.newton_solve();



  problem.doc_solution(problem.bulk_mesh_pt());
//  GenericProblemSetup::doc_solution(problem.bulk_mesh_pt(),0);

//  GenericProblemSetup::unsteady_run(&problem,
//                                    &doc_linear_solver_info,
//                                    problem.bulk_mesh_pt());


  }

  //////////////////////////////////////////////////////////////////////////
  ////////////// Outputting results ////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

//  if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
//  {
//    // Get the global oomph-lib communicator 
//    const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();
//
//    // my rank and number of processors. 
//    // This is used later for putting the data.
//    const unsigned my_rank = comm_pt->my_rank();
//    const unsigned nproc = comm_pt->nproc();
//
//    // Variable to indicate if we want to output to a file or not.
//    bool output_to_file = false;
//
//    // The output file.
//    std::ofstream outfile;
//
//    // If we want to output to a file, we create the outfile.
//    if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir"))
//    {
//      output_to_file = true;
//      std::ostringstream filename_stream;
//      filename_stream << NSPP::Itstime_dir_str<<"/"
//        << NSPP::Label_str
//        <<"NP"<<nproc<<"R"<<my_rank;
//      outfile.open(filename_stream.str().c_str());
//    }
//
//    // Stringstream to hold the results. We do not output the results
//    // (timing/iteration counts) as we get it since it will interlace with the
//    // other processors and becomes hard to read.
//    std::ostringstream results_stream;
//
//    // Get the 3D vector which holds the iteration counts and timing results.
//    Vector<Vector<Vector<double> > > iters_times
//      = NSPP::Doc_linear_solver_info_pt->iterations_and_times();
//
//    // Since this is a steady state problem, there is only
//    // one "time step", thus it is essentially a 2D vector 
//    // (the outer-most vector is of size 1).
//
//    // Loop over the time steps and output the iterations, prec setup time and
//    // linear solver time.
//    unsigned ntimestep = iters_times.size();
//    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
//    {
//      ResultsFormat::format_rayits(intimestep,&iters_times,&results_stream);
//    }
//    
//    ResultsFormat::format_rayavgits(&iters_times,&results_stream);
//    ResultsFormat::format_rayavavgits(&iters_times,&results_stream);
//
//    // Now doing the preconditioner setup time.
//    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
//    {
//      ResultsFormat::format_prectime(intimestep,&iters_times,&results_stream);
//    }
//
//    ResultsFormat::format_avgprectime(&iters_times,&results_stream);
//    ResultsFormat::format_avavgprectime(&iters_times,&results_stream);
//    // Now doing the linear solver time.
//    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
//    {
//      ResultsFormat::format_solvertime(intimestep,&iters_times,&results_stream);
//    }
//    
//    ResultsFormat::format_avgsolvertime(&iters_times,&results_stream);
//    ResultsFormat::format_avavgsolvertime(&iters_times,&results_stream); 
//
//    // Print the result to oomph_info one processor at a time...
//    // This still doesn't seem to always work, since there are other calls
//    // to oomph_info before this one...
//    for (unsigned proc_i = 0; proc_i < nproc; proc_i++) 
//    {
//      if(proc_i == my_rank)
//      {
//        oomph_info << "\n" 
//          << "========================================================\n"
//          << results_stream.str()
//          << "========================================================\n"
//          << "\n" << std::endl;
//      }
//      MPI_Barrier(MPI_COMM_WORLD);
//    }
//
//    if(output_to_file)
//    {
//      outfile << "\n" << results_stream.str();
//      outfile.close();
//    }
//  }



#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif


} // end_of_main


