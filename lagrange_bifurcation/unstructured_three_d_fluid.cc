//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//           Version 0.90. August 3, 2009.
//LIC// 
//LIC// Copyright (C) 2006-2009 Matthias Heil and Andrew Hazel
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
// Driver code for a simple unstructured fluid problem using a mesh
// generated from an input file generated by the 3d mesh generator
// tetgen


//Generic routines
#include "generic.h"
#include "constitutive.h"
#include "navier_stokes.h"

// Get the mesh
#include "meshes/tetgen_mesh.h"
#include "meshes/brick_from_tet_mesh.h" 

// My own header
//#include "./../rayheader.h"
#include "./../ray_preconditioner_creation.h"
#include "./../ray_navier_stokes_parameters.h"
#include "./../ray_general_problem_parameters.h"

using namespace std;
using namespace oomph;


namespace GenProbHelpers = GeneralProblemHelpers;
namespace PrecHelpers = PreconditionerHelpers;
namespace NSHelpers = NavierStokesHelpers;

// Alias the namespace for convenience.
//namespace NSPP = NavierStokesProblemParameters;
//namespace LPH = LagrangianPreconditionerHelpers;
//namespace BL = BifurcationLagrange;

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
namespace ProblemHelpers
{

  int Prob_id = -1;

  std::string Prob_str = "";
  std::string Mesh_folder_str = "";
  double Mesh_area = 0.0;

  inline void specify_command_line_flags()
  {
    CommandLineArgs::specify_command_line_flag("--prob_id", &Prob_id);
    CommandLineArgs::specify_command_line_flag("--mesh_area", &Mesh_area);
  }

  inline void setup_command_line_flags()
  {
        // Check the problem id
    if(!CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
    {
      std::ostringstream err_msg;
      err_msg << "Please set --prob_id" << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION); 
    }


        // Set Noel_str, used for book keeping.
    if(CommandLineArgs::command_line_flag_has_been_set("--mesh_area"))
    {
   // Set the string to load the files.
   // The Mesh_area parameter is a double.
   // The actual mesh files are in tetgen_files/xdyz
   // where d represents the decimal place.
   // So we need to replace the decimal in the RNS::Mesh_area parameter with
   // d.
   std::ostringstream tmp_stringstream;
   tmp_stringstream << Mesh_area;
   Mesh_folder_str = tmp_stringstream.str();
   std::replace(Mesh_folder_str.begin(), Mesh_folder_str.end(),
                '.','d');
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please supply the min element area using --mesh_area\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }


  }


 // Get the prescribed inflow velocity for the steady state problem.
 // This is assumed that x and y are in the range [-1,1].
 inline double get_prescribed_inflow(const double& x,
                                     const double& y)
 {
   return (1 - x) * (x- (-1)) * (1 - y) * (y - (-1));
 }


 // Scale the steady state prescribed velocity inflow above
 // by the time.
 inline double get_prescribed_inflow(const double& t,
                                     const double& x,
                                     const double& y)
 {
   const double scaling = -cos(MathematicalConstants::Pi*t)/2.0 + 0.5;
   return (get_prescribed_inflow(x,y) * scaling);
 } 


 inline std::string prob_str()
 {
   std::string prob_str = "";
   if(Prob_id == 0)
   {
     prob_str = "Bi";
   }
   else
   {
      std::ostringstream error_message;
      error_message << "No other problem ids done." << std::endl;

      throw OomphLibError(error_message.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
   }
   return prob_str;
 }

}

namespace ProbHelpers = ProblemHelpers;

//======start_problem_class===========================================
/// Unstructured fluid problem
//====================================================================
template<class ELEMENT>
class UnstructuredFluidProblem : public Problem
{

public:

 /// Constructor: 
 UnstructuredFluidProblem();

 /// Destructor (empty)
 ~UnstructuredFluidProblem(){}


 void actions_before_implicit_timestep()
 {
   if(GenProbHelpers::Time_type != GenProbHelpers::Time_type_STEADY)
   {
     // NOTE: before an implicit time step, we clear the previous times.
     // (This has no effect if no times were added i.e.
     // we are starting the first implicit time step)
     //
     // The logic for this is that we want to keep only the time step we
     // use, not the rejected ones.
     //
     // In the unsteady solve loop, it is there at we add new storage for
     // a new time step.
     Doc_linear_solver_info_pt->clear_current_time_step();

     const unsigned ibound = Inflow_boundary_id[0];
     const unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod = 0; inod < num_nod; inod++) 
     {
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);
       const double x = nod_pt->x(0);
       const double y = nod_pt->x(1);

       const double time = time_pt()->time();

       const double uz = ProbHelpers::get_prescribed_inflow(time,x,y);

       nod_pt->set_value(0,0.0);
       nod_pt->set_value(1,0.0);
       nod_pt->set_value(2,uz);
     }
   }
 }
 
 /// Update the problem specs before solve. 
 void actions_before_newton_solve()
 {

//   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
//   {
//    // Initialise counters for each newton solve.
//    Doc_linear_solver_info_pt->setup_new_time_step();
//   }
//
//
//
//
//   if(NSPP::Steady_state)
//   {
//     // Set the inflow, this is boundary 0
//     const unsigned ibound = Inflow_boundary_id[0];
//     const unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
//     for (unsigned inod = 0; inod < num_nod; inod++) 
//     {
//       Node*  nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);
//       const double x = nod_pt->x(0);
//       const double y = nod_pt->x(1);
//       //const double z = nod_pt->x(2);
//
//       // x and y are in the range [-1,1]
//       const double uz = BL::get_prescribed_inflow(x,y);
//
//       nod_pt->set_value(0,0.0);
//       nod_pt->set_value(1,0.0);
//       nod_pt->set_value(2,uz);
//     }
//   }
 }

 void actions_after_newton_step()
 {
   if(GenProbHelpers::Solver_type !=
       GenProbHelpers::Solver_type_DIRECT_SOLVE)
   {
     GenProbHelpers::doc_iter_times(this,Doc_linear_solver_info_pt);
   }
 }

 void actions_before_distribute()
 {
     GenProbHelpers::delete_flux_elements(Surface_mesh_pt);
     rebuild_global_mesh();
 }

 void actions_after_distribute()
 {
     const unsigned n_outflow_boundary = Outflow_boundary_id.size();
     for (unsigned ibound = 0; ibound < n_outflow_boundary; ibound++) 
     {
       const unsigned current_bound = Outflow_boundary_id[ibound];
       create_parall_outflow_lagrange_elements(current_bound,
           Tangent_direction,
           Bulk_mesh_pt,
           Surface_mesh_pt);
     }
     rebuild_global_mesh();
 }


 /// Global error norm for adaptive time-stepping
 double global_temporal_error_norm();

 void create_parall_outflow_lagrange_elements(const unsigned &b,
                                              Vector<double> &tangent_direction,
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);
 //private:
 /// Pointer to the "bulk" mesh
 Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}


 /// Bulk fluid mesh
// TetgenMesh<ELEMENT>* Bulk_mesh_pt;
 Mesh* Bulk_mesh_pt;

 Mesh* Surface_mesh_pt;

 Vector<double> Tangent_direction;

 Preconditioner* Prec_pt;

 Preconditioner* NS_matrix_preconditioner_pt;

 Preconditioner* P_matrix_preconditioner_pt;
 Preconditioner* F_matrix_preconditioner_pt;

 IterativeLinearSolver* Solver_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;

 /// \short IDs of fluid mesh boundaries along which inflow boundary conditions
 /// are applied
 Vector<unsigned> Inflow_boundary_id;

// /// \short IDs of fluid mesh boundaries along which inflow boundary conditions
// /// are applied
 Vector<unsigned> Outflow_boundary_id;

};



//==========start_constructor=============================================
/// Constructor for unstructured 3D fluid problem
//========================================================================
  template<class ELEMENT>
UnstructuredFluidProblem<ELEMENT>::UnstructuredFluidProblem()
{ 
Doc_linear_solver_info_pt = GenProbHelpers::Doc_linear_solver_info_pt;

  // Add a new time stepper if not doing steady state.
  if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_STEADY)
  {
    oomph_info 
      << "RAYINFO: Doing steady state, no time stepper added." << std::endl;
  }
  else if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_ADAPT)
  {
    oomph_info 
      << "RAYINFO: Adding adaptive time stepper" << std::endl;
    add_time_stepper_pt(new BDF<2>(true));
  }
  else if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_FIXED)
  {
    oomph_info << "RAYINFO: Adding non-adaptive time stepper" << std::endl;
    add_time_stepper_pt(new BDF<2>);
  }
  else
  {
    std::ostringstream err_msg;
    err_msg << "Time stepper for Time_type: "
      << GenProbHelpers::Time_type << std::endl;

    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }



  //Create fluid bulk mesh, sub-dividing "corner" elements
  string mesh_folder = "tetgen_files/" + ProbHelpers::Mesh_folder_str +"/";

  string node_file_name=mesh_folder+"fsi_bifurcation_fluid.1.node";
  string element_file_name=mesh_folder+"fsi_bifurcation_fluid.1.ele";
  string face_file_name=mesh_folder+"fsi_bifurcation_fluid.1.face";
  bool split_corner_elements=true;

  // Check if we want the tetrahedral or hexahedral
//  if(NSPP::Mesh_type == NSPP::MeshType_TETRAHEDRAL)
//  {
//    if(NSPP::Steady_state)
//    {
//      Bulk_mesh_pt =  new TetgenMesh<ELEMENT>(node_file_name,
//          element_file_name,
//          face_file_name,
//          split_corner_elements);
//    }
//    else
//    {
//      Bulk_mesh_pt =  new TetgenMesh<ELEMENT>(node_file_name,
//          element_file_name,
//          face_file_name,
//          split_corner_elements,
//          time_stepper_pt());
//    }
//  }
//  else if(NSPP::Mesh_type == NSPP::MeshType_HEXAHEDRAL)
  {
    if(GenProbHelpers::Time_type != GenProbHelpers::Time_type_STEADY)
    {
      Bulk_mesh_pt = new BrickFromTetMesh<ELEMENT>(node_file_name,
          element_file_name,
          face_file_name,
          split_corner_elements);
    }
    else
    {
      Bulk_mesh_pt = new BrickFromTetMesh<ELEMENT>(node_file_name,
          element_file_name,
          face_file_name,
          split_corner_elements,
          time_stepper_pt());
    }
  }
//  else
//  {
//       std::ostringstream err_msg;
//        err_msg << "Please set --mesh_type" << std::endl;
//
//        throw OomphLibError(err_msg.str(),
//            OOMPH_CURRENT_FUNCTION,
//            OOMPH_EXCEPTION_LOCATION);
//  }

  oomph_info << "Calling setup_boundary_element_info()" << std::endl; 
  // Find elements next to boundaries
  Bulk_mesh_pt->setup_boundary_element_info();
  oomph_info << "Done setup_boundary_element_info()" << std::endl; 


  // The following corresponds to the boundaries as specified by
  // facets in the tetgen input:

  Inflow_boundary_id.resize(1);
  Inflow_boundary_id[0] = 0;

  Outflow_boundary_id.resize(2);
  Outflow_boundary_id[0] = 1;
  Outflow_boundary_id[1] = 2;


  // First pin all boundary nodes, then we unpin those on the outflow boundary.
  const unsigned num_bound = Bulk_mesh_pt->nboundary();
  for (unsigned ibound = 0; ibound < num_bound; ibound++) 
  {
    const unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod = 0; inod < num_nod; inod++) 
    {
      // Loop over velocity nodes
      const unsigned n_velocity_nodes = 3;
      for (unsigned iv = 0; iv < n_velocity_nodes; iv++) 
      {
        // Locally cache the node, since we use it more than once.
        Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);

        // Pin and (just to be safe!) set the value to zero.
        nod_pt->pin(iv);
        nod_pt->set_value(iv,0.0);
      } // for - loop over the velocity nodes.
    } // for - loop over the nodes
  } // for - loop over all the boundaries

  // Unpin the nodes at the outflow boundary, but only the node which
  // is on a single boundary!
  // Outflow boundaries are one and two
  const unsigned n_outflow_boundary = Outflow_boundary_id.size();

  for (unsigned ibound = 0; ibound < n_outflow_boundary; ibound++)
  {
    const unsigned current_bound = Outflow_boundary_id[ibound];

    const unsigned num_nod = Bulk_mesh_pt->nboundary_node(current_bound);

    for (unsigned inod = 0; inod < num_nod; inod++) 
    {
      Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(current_bound,inod);

      // Only free if node is ONLY on a single boundary
      std::set<unsigned>*bnd_pt=0;
      nod_pt->get_boundaries_pt(bnd_pt);
      if (bnd_pt != 0) 
      {
        if (bnd_pt->size()<2) 
        {
          nod_pt->unpin(0);
          nod_pt->unpin(1);
          nod_pt->unpin(2);
        } // if there is less than two boundaries
      } // if the boundary pointer is not null
    } // for - loop over nodes
  } // for - loop over outflow boundaries


  // Create the surface mesh for the parallel outflow elements
  // on boundary 1 and 2. To be safe, we give a general direction for the
  // tangent vector. Recall that now we have two outflow faces which
  // planes intersect. This means that the automatically calculated tangent
  // vector may switch (even within an element if the face is unfortunately
  // aligned with one of the axis).
  Tangent_direction.resize(3,0);
  Tangent_direction[0] = 0;
  Tangent_direction[1] = 1;
  Tangent_direction[2] = 0;

  Surface_mesh_pt = new Mesh;
  for (unsigned ibound = 0; ibound < n_outflow_boundary; ibound++) 
  {
    const unsigned current_bound = Outflow_boundary_id[ibound];
    create_parall_outflow_lagrange_elements(current_bound,
        Tangent_direction,
        Bulk_mesh_pt,
        Surface_mesh_pt);
  }

  // Combine all the sub meshes.
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);

  // Build the global mesh
  build_global_mesh();

  // Set up equation numbering scheme
  oomph_info << "Number of equations: " << assign_eqn_numbers() << std::endl;

//  // Complete the build of the fluid elements so they are fully functional
//  //----------------------------------------------------------------------
  const unsigned n_element = Bulk_mesh_pt->nelement();
  for(unsigned e=0;e<n_element;e++)
  {
    // Upcast from GeneralisedElement to the present element
    ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

    //Set the Reynolds number
    el_pt->re_pt() = &NSHelpers::Rey;
    el_pt->re_st_pt() = &NSHelpers::Rey;
  } 


  F_matrix_preconditioner_pt
    = PrecHelpers::create_f_p_amg_preconditioner(PrecHelpers::F_amg_param,0);
  P_matrix_preconditioner_pt
    = PrecHelpers::create_f_p_amg_preconditioner(PrecHelpers::P_amg_param,1);

  NS_matrix_preconditioner_pt = PrecHelpers::create_lsc_preconditioner(
      this,
      Bulk_mesh_pt,
      F_matrix_preconditioner_pt,
      P_matrix_preconditioner_pt);

  Vector<Mesh*> mesh_pt(2,0);
  mesh_pt[0] = Bulk_mesh_pt;
  mesh_pt[1] = Surface_mesh_pt;
  Prec_pt = PrecHelpers::create_lgr_precondiitoner(
      mesh_pt,
      PrecHelpers::W_solver,
      NS_matrix_preconditioner_pt);

  const double solver_tol = 1.0e-6;
  const double newton_tol = 1.0e-6;

  Solver_pt = GenProbHelpers::setup_solver(
      GenProbHelpers::Max_solver_iteration,
      solver_tol,newton_tol,
      this,Prec_pt);



} // end constructor



template<class ELEMENT>
double UnstructuredFluidProblem<ELEMENT>::global_temporal_error_norm()
{
  return NSHelpers::global_temporal_error_norm(
      this,
      NSHelpers::Dim,
      Bulk_mesh_pt);

} // end of global_temporal_error_norm



//============start_of_fluid_traction_elements==============================
/// Create fluid traction elements 
//=======================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::create_parall_outflow_lagrange_elements
(const unsigned &b, Vector<double>& tangent_direction,
 Mesh* const &bulk_mesh_pt, Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);
   
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));
     
   //What is the index of the face of the element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   {
    // Build the corresponding impose_impenetrability_element
    ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new
     ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
                                           face_index);

    flux_element_pt->set_tangent_direction(&tangent_direction);
    surface_mesh_pt->add_element_pt(flux_element_pt);

    // Loop over the nodes
    unsigned nnod=flux_element_pt->nnode();
    for (unsigned j=0;j<nnod;j++)
     {
      Node* nod_pt = flux_element_pt->node_pt(j);
           
      // Determine which outflow boundary it is, left or right?
           
      if (  (nod_pt->is_on_boundary(7))||(nod_pt->is_on_boundary(8))
            ||(nod_pt->is_on_boundary(9))||(nod_pt->is_on_boundary(10))
            ||(nod_pt->is_on_boundary(11))||(nod_pt->is_on_boundary(12))
            ||(nod_pt->is_on_boundary(13))||(nod_pt->is_on_boundary(14)))
       {
        // How many nodal values were used by the "bulk" element
        // that originally created this node?
        unsigned n_bulk_value=flux_element_pt->nbulk_value(j);

        // The remaining ones are Lagrange multipliers and we pin them.
        unsigned nval=nod_pt->nvalue();
        for (unsigned j=n_bulk_value;j<nval;j++)
         {
          nod_pt->pin(j);
         }
       }
     }
   }
  }
} // end of create_parall_outflow_lagrange_elements



std::string create_label()
{
  
  std::string label = ProbHelpers::prob_str()
                      + NSHelpers::create_label()
                      + PrecHelpers::Lgr_prec_str
                      + ProbHelpers::Mesh_folder_str;
  return label;
}



//=============start_main=================================================
/// Demonstrate how to solve an unstructured 3D fluids problem
//========================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Set up doc info - used to store information on solver and iteration time.
  DocLinearSolverInfo doc_linear_solver_info;

  // The TWO things:
  GenProbHelpers::Doc_linear_solver_info_pt = &doc_linear_solver_info;
  NSHelpers::Dim = 3;

  // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  // From GeneralProblemHelpers, we need to set:
  // --time_type, 
  //   0 = steady state, 
  //   1 = adaptive, 
  //   2 = fixed.
  // --solver_type, 
  //   0 = exact solver, 
  //   1 = OOMPH-LIB's GMRES, 
  //   2 = trilinos GMRES
  // -- dist_prob (no arguments)
  //
  // --max_solver_iter - set an integer
  //
  // --dt - only set if doing fixed time stepping.
  // --time_start - always set if not steady state
  // --time_end - always set if not steady state
  //
  // -- doc_soln [soln_dir]
  //
  // --itstimedir [results_dir]
  //
  GenProbHelpers::specify_command_line_flags();


  // From NavierStokesHelpers, set:
  // --visc 0 or 1
  // --rey a double
  // --rey_start - only set if looping through reynolds numbers
  // --rey_incre - same as above
  // --rey_end - same as above
  NSHelpers::specify_command_line_flags();

  // From PreconditionerHelpers, set:
  // --f_solver 0 (exact) or 1 (amg)
  // --p_solver 0 (exact) or 1 (amg)
  //
  // if f_solver or p_solver is 1, we NEED to additionally set (for f_solver)
  // --f_amg_iter (usually 1)
  //
  // --f_amg_smiter (usually 2)
  //
  // --f_amg_sim_smoo
  //   0 = Jacobi, IMPORTANT: set --f_amg_damp (to something like 1).
  //   1 = Gauss-Seidel, sequential
  //       (very slow in parallel!)
  //   2 = Gauss-Seidel, interior points in parallel, boundary sequential
  //       (slow in parallel!)
  //   3 = hybrid Gauss-Seidel or SOR, forward solve
  //   4 = hybrid Gauss-Seidel or SOR, backward solve
  //   6 = hybrid symmetric Gauss-Seidel or SSOR 
  //
  // OR set this (NOT BOTH):
  // --f_amg_com_smoo
  //    
  //   6 = Schwarz
  //   7 = Pilut
  //   8 = ParaSails
  //   9 = Euclid
  //
  // --f_amg_str - strength of dependence
  //
  // --f_amg_coarse:
  //    0 = CLJP (parallel coarsening using independent sets)
  //    1 = classical RS with no boundary treatment (not recommended
  //        in parallel)
  //    3 = modified RS with 3rd pass to add C points on the boundaries
  //    6 = Falgout (uses 1 then CLJP using interior coarse points as
  //        first independent set) THIS IS DEFAULT ON DOCUMENTATION
  //    8 = PMIS (parallel coarsening using independent sets - lower
  //        complexities than 0, maybe also slower convergence)
  //    10= HMIS (one pass RS on each processor then PMIS on interior
  //        coarse points as first independent set)
  //    11= One pass RS on each processor (not recommended)
  //
  // --print_f_hypre - to print the hypre parameters to confirm.
  //    The information is extracted from the preconditioner after it has
  //    been created, so it is good to always do this.
  //
  // REPEAT FOR p solver if it is also using AMG
  // THERE ARE 6 things to set!
  //
  PrecHelpers::specify_command_line_flags();

  // --prob_id
  // --ang
  // --noel
  ProbHelpers::specify_command_line_flags();



  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags(); 



  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////

  GenProbHelpers::setup_command_line_flags();
  NSHelpers::setup_command_line_flags();
  PrecHelpers::setup_command_line_flags();
  ProbHelpers::setup_command_line_flags();

  std::string label = "";

//  if(NSPP::Mesh_type == NSPP::MeshType_TETRAHEDRAL)
//  {
//    //Set up the problem
//    UnstructuredFluidProblem<TTaylorHoodElement<3> > problem;
//
//    if(NSPP::Distribute_problem)
//    {
//      problem.distribute();
//    }
//
//    NSPP::Label_str = create_label();
//
//    time_t rawtime;
//    time(&rawtime);
//
//    std::cout << "RAYDOING: "
//      << NSPP::Label_str
//      << " on " << ctime(&rawtime) << std::endl;
//
//
//    // Solve the problem
//    if(NSPP::Steady_state)
//    {
//      problem.newton_solve();
//
//      //Output solution
//      problem.doc_solution(0);
//    }
//    else
//    {
//      problem.unsteady_run();
//    }
//  }
//  else if(NSPP::Mesh_type == NSPP::MeshType_HEXAHEDRAL)
  {
    //Set up the problem
    UnstructuredFluidProblem<QTaylorHoodElement<3> > problem;

  if(GenProbHelpers::Distribute_problem)
  {
    const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();
    const unsigned nproc = comm_pt->nproc();

    if(nproc == 1)
    {
      oomph_info << "RAYINFO: only 1 core, "
                 << "not distributing problem." << std::endl;
    }
    else
    {
      oomph_info << "RAYINFO: I am distributing the problem" << std::endl;

      problem.distribute();
    }
  }
  GenProbHelpers::Distribute_problem = problem.distributed();
  oomph_info << "Problem.distributed() is " 
             << GenProbHelpers::Distribute_problem << std::endl; 

  label = create_label();
    time_t rawtime;
    time(&rawtime);

    std::cout << "RAYDOING: "
      << label
      << " on " << ctime(&rawtime) << std::endl;

  // There are two types of solves, one for steady state, another for
  // time stepping.
  if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_STEADY)
  {
    problem.newton_solve();

    if(GenProbHelpers::Doc_soln_flag)
    {
      GenProbHelpers::doc_solution(problem.bulk_mesh_pt(),
          GenProbHelpers::Soln_dir_str,
          label);
    }
  }
  else
  {
    GenProbHelpers::unsteady_run(&problem,
        problem.bulk_mesh_pt(),
        &doc_linear_solver_info,
        label,
        GenProbHelpers::Soln_dir_str);
  }

  }

  //////////////////////////////////////////////////////////////////////////
  ////////////// Outputting results ////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////


  if(GenProbHelpers::Solver_type != 
     GenProbHelpers::Solver_type_DIRECT_SOLVE)
  {
    // Get the global oomph-lib communicator 
    const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

    // my rank and number of processors. 
    // This is used later for putting the data.
    const unsigned my_rank = comm_pt->my_rank();
    const unsigned nproc = comm_pt->nproc();

    // Variable to indicate if we want to output to a file or not.
    bool output_to_file = false;

    // The output file.
    std::ofstream outfile;

    // If we want to output to a file, we create the outfile.
    if(GenProbHelpers::Doc_time_flag)
    {
      output_to_file = true;
      std::ostringstream filename_stream;
      filename_stream << GenProbHelpers::Itstime_dir_str<<"/"
        << label
        <<"NP"<<nproc<<"R"<<my_rank;
      outfile.open(filename_stream.str().c_str());
    }

    // Stringstream to hold the results. We do not output the results
    // (timing/iteration counts) as we get it since it will interlace with the
    // other processors and becomes hard to read.
    std::ostringstream results_stream;

    // Get the 3D vector which holds the iteration counts and timing results.
    Vector<Vector<Vector<double> > > iters_times
      = GenProbHelpers::Doc_linear_solver_info_pt->iterations_and_times();

    // Since this is a steady state problem, there is only
    // one "time step", thus it is essentially a 2D vector 
    // (the outer-most vector is of size 1).

    // Loop over the time steps and output the iterations, prec setup time and
    // linear solver time.
    unsigned ntimestep = iters_times.size();
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat2::format_rayits(intimestep,&iters_times,&results_stream);
    }
    
    ResultsFormat2::format_rayavgits(&iters_times,&results_stream);
    ResultsFormat2::format_rayavavgits(&iters_times,&results_stream);

    // Now doing the preconditioner setup time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat2::format_prectime(intimestep,&iters_times,&results_stream);
    }

    ResultsFormat2::format_avgprectime(&iters_times,&results_stream);
    ResultsFormat2::format_avavgprectime(&iters_times,&results_stream);
    // Now doing the linear solver time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat2::format_solvertime(intimestep,&iters_times,&results_stream);
    }
    
    ResultsFormat2::format_avgsolvertime(&iters_times,&results_stream);
    ResultsFormat2::format_avavgsolvertime(&iters_times,&results_stream); 

    // Print the result to oomph_info one processor at a time...
    // This still doesn't seem to always work, since there are other calls
    // to oomph_info before this one...
    for (unsigned proc_i = 0; proc_i < nproc; proc_i++) 
    {
      if(proc_i == my_rank)
      {
        oomph_info << "\n" 
          << "========================================================\n"
          << results_stream.str()
          << "========================================================\n"
          << "\n" << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }

    if(output_to_file)
    {
      outfile << "\n" << results_stream.str();
      outfile.close();
    }
  }



#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS); 

} // end_of_main




