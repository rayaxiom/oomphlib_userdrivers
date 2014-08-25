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

//#include "./../rayheader.h"
#include "./../ray_preconditioner_creation.h"
#include "./../ray_navier_stokes_parameters.h"
#include "./../ray_general_problem_parameters.h"

using namespace std;

using namespace oomph;


namespace GenProbHelpers = GeneralProblemHelpers;
namespace PrecHelpers = PreconditionerHelpers;
namespace NSHelpers = NavierStokesHelpers;

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

namespace Global_Parameters
{

  // Problem specific parameters
  unsigned Noel = 0;
  int Prob_id = -1;

  const double Length = 1.0;

  inline std::string prob_str()
  {
    std::string prob_str = "";
    if(Prob_id == 0)
    {
      prob_str = "Cuq";
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



  inline std::string noel_str()
  {
    std::ostringstream noel_stream;
    if(Noel != 0)
    {
      noel_stream << "N" << Noel;
    }
    else
    {
      std::ostringstream error_message;
      error_message << "Noel is zero. Have you called --noel?" << std::endl;

      throw OomphLibError(error_message.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    return noel_stream.str();
  }





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

    // Note: +0.51, so at time = 0, there is still come velocity
    const double ux_scaling = -cos(MathematicalConstants::Pi*t)/2.0 + 0.51;
    return get_prescribed_inflow_for_quarter(y,z) * ux_scaling;
  } 


}

namespace GlobalParam = Global_Parameters;


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

      delete Solver_pt;
      delete Prec_pt;
      delete P_matrix_preconditioner_pt;
      delete F_matrix_preconditioner_pt;
      delete Bulk_mesh_pt;
      //   delete Prec_pt;
      //   delete P_matrix_preconditioner_pt;
      //   delete F_matrix_preconditioner_pt;
    }

    //////////////////////////////

    void actions_before_implicit_timestep()
    {
      if(GenProbHelpers::Time_type != GenProbHelpers::Time_type_FIXED)
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

        // Inflow in upper half of inflow boundary
        const unsigned ibound=Inflow_boundary; 
        const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
        for (unsigned inod=0;inod<num_nod;inod++)
        {
          Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

          std::set<unsigned>* bnd_pt=0;

          nod_pt->get_boundaries_pt(bnd_pt);
          if(bnd_pt !=0)
          {
            if(bnd_pt->size() < 2)
            {
              // Imposing velocity along x axis, depends on y-z plane only.
              const double y=nod_pt->x(1);
              const double z=nod_pt->x(2);

              // Only do the inflow for a quarter of the boundary
              if((y > 0.5) && (z > 0.5))
              {

                const double time=time_pt()->time();

                double ux 
                  = GlobalParam::get_prescribed_inflow_for_quarter(time,y,z);

                nod_pt->set_value(0,ux);
                nod_pt->set_value(1,0.0);
                nod_pt->set_value(2,0.0);
              }
            }
          }
        }
      }
    } // end of actions_before_implicit_timestep



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
      // If steady state, we set the boundary conditions
      // before the newton solve. Otherwise we set it before
      // implicit time step.
      if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_STEADY)
      {
        if(GenProbHelpers::Solver_type 
            != GenProbHelpers::Solver_type_DIRECT_SOLVE)
        {
          // Start a new "time step"
          Doc_linear_solver_info_pt->setup_new_time_step();
        }

        // Inflow in upper half of inflow boundary
        const unsigned ibound=Inflow_boundary; 
        const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
        for (unsigned inod=0;inod<num_nod;inod++)
        {
          Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
          std::set<unsigned>* bnd_pt=0;
          nod_pt->get_boundaries_pt(bnd_pt);
          if(bnd_pt != 0)
          {
            if(bnd_pt->size() < 2)
            {


              const double y=nod_pt->x(1);
              const double z=nod_pt->x(2);

              if( (y > 0.5) && (z > 0.5) )
              {

                const double ux 
                  = GlobalParam::get_prescribed_inflow_for_quarter(y,z);

                nod_pt->set_value(0,ux);
                nod_pt->set_value(1,0.0);
                nod_pt->set_value(2,0.0);
              }
            }
          }
        }
      }
    } // end_of_actions_before_newton_solve


    void actions_after_newton_step()
    {
      if(GenProbHelpers::Solver_type 
          != GenProbHelpers::Solver_type_DIRECT_SOLVE)
      {
        GenProbHelpers::doc_iter_times(this,Doc_linear_solver_info_pt);
      }
    }

    /// Global error norm for adaptive time-stepping
    double global_temporal_error_norm();


    /// Pointer to the "bulk" mesh
    Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}


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
  Doc_linear_solver_info_pt = GenProbHelpers::Doc_linear_solver_info_pt;

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


  // Setup mesh
  const unsigned noel = GlobalParam::Noel;
  const double length = GlobalParam::Length;

  if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_STEADY)
  {
    Bulk_mesh_pt = 
      new SimpleCubicMesh<ELEMENT>(noel,noel,noel,
          length, length, length);
  }
  else if((GenProbHelpers::Time_type == GenProbHelpers::Time_type_ADAPT)
      ||GenProbHelpers::Time_type == GenProbHelpers::Time_type_FIXED)
  {
    Bulk_mesh_pt = 
      new SimpleCubicMesh<ELEMENT>(noel,noel,noel, 
          length, length, length,
          time_stepper_pt());
  }
  else
  {
    std::ostringstream err_msg;
    err_msg << "No mesh set up for Time_type: "
      << GenProbHelpers::Time_type << std::endl;

    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }


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
    el_pt->re_pt() = &NSHelpers::Rey;
//    el_pt->re_st_pt() = &NSHelpers::Rey;
  } // end loop over elements


  // Now set the first pressure value in element 0 to 0.0
  // if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);

  // Setup equation numbering scheme
  oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 


  // Set up preconditioner and solver.
  F_matrix_preconditioner_pt 
    = PrecHelpers::create_f_p_amg_preconditioner(PrecHelpers::F_amg_param,0);
  P_matrix_preconditioner_pt 
    = PrecHelpers::create_f_p_amg_preconditioner(PrecHelpers::P_amg_param,1);

  Prec_pt = PrecHelpers::create_lsc_preconditioner(this,Bulk_mesh_pt,
      F_matrix_preconditioner_pt,
      P_matrix_preconditioner_pt);
  const double solver_tol = 1.0e-6;
  const double newton_tol = 1.0e-6;

  Solver_pt = GenProbHelpers::setup_solver(
      GenProbHelpers::Max_solver_iteration,
      solver_tol,newton_tol,
      this,Prec_pt);

} // end_of_constructor


  template<class ELEMENT>
double CubeProblem<ELEMENT>::global_temporal_error_norm()
{
  return NSHelpers::global_temporal_error_norm(this,
      NSHelpers::Dim,
      Bulk_mesh_pt);
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


  std::string label = Global_Parameters::prob_str()
    + NSHelpers::create_label()
    + PrecHelpers::NS_prec_str
    + Global_Parameters::noel_str();
  //  std::string label = CL::prob_str()
  //                      + NSPP::create_label() 
  //                      + LPH::create_label() 
  //                      + CL::ang_deg_str() + CL::noel_str();
  return label;
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////




inline void problem_specific_setup_commandline_flags()
{
  CommandLineArgs::specify_command_line_flag("--noel", 
      &GlobalParam::Noel);

  CommandLineArgs::specify_command_line_flag("--prob_id", 
      &GlobalParam::Prob_id);
}


inline void problem_specific_generic_setup()
{
  if(!CommandLineArgs::command_line_flag_has_been_set("--noel"))
  {
    std::ostringstream err_msg;
    err_msg << "Please set --noel" << std::endl;
    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }

  if(!CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
  {
    std::ostringstream err_msg;
    err_msg << "Please set --prob_id" << std::endl;
    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }

}


//==start_of_main======================================================
/// Driver for Fp preconditioner
//=====================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Set up doc info - used to store information on solver and iteration time.
  DocLinearSolverInfo doc_linear_solver_info;

  GenProbHelpers::Doc_linear_solver_info_pt = &doc_linear_solver_info;
  NSHelpers::Dim = 3;

  // Store commandline arguments
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
  GenProbHelpers::setup_commandline_flags();

  // From NavierStokesHelpers, set:
  // --visc 0 or 1
  // --rey a double
  // --rey_start - only set if looping through reynolds numbers
  // --rey_incre - same as above
  // --rey_end - same as above
  NSHelpers::setup_commandline_flags();

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
  PrecHelpers::setup_commandline_flags();

  // --noel number of elements in 1D
  // --prob_id - currently, the only prob id is 0
  problem_specific_setup_commandline_flags();

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();


  GenProbHelpers::generic_setup();
  NSHelpers::generic_setup();
  PrecHelpers::generic_setup();
  problem_specific_generic_setup();



  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////
  std::string label = "";

  {
    // Build the problem 
    CubeProblem <QTaylorHoodElement<3> >problem;


    // Solve the problem 
    //              problem.newton_solve();

    if(GenProbHelpers::Distribute_problem)
    {
      oomph_info << "RAYINFO: I am distributing the problem" << std::endl; 
      
      problem.distribute();
    }

    label = create_label();
    //  NSPP::Label_str = create_label();

    time_t rawtime;
    time(&rawtime);

    std::cout << "RAYDOING: "
      << label
      << " on " << ctime(&rawtime) << std::endl;

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

    //  GenericProblemSetup::unsteady_run(&problem,
    //                                    &doc_linear_solver_info,
    //                                    problem.bulk_mesh_pt());


  }

  //////////////////////////////////////////////////////////////////////////
  ////////////// Outputting results ////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

    if(GenProbHelpers::Solver_type != GenProbHelpers::Solver_type_DIRECT_SOLVE)
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
        ResultsFormat::format_rayits(intimestep,&iters_times,&results_stream);
      }
      
      ResultsFormat::format_rayavgits(&iters_times,&results_stream);
      ResultsFormat::format_rayavavgits(&iters_times,&results_stream);
  
      // Now doing the preconditioner setup time.
      for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
      {
        ResultsFormat::format_prectime(intimestep,&iters_times,&results_stream);
      }
  
      ResultsFormat::format_avgprectime(&iters_times,&results_stream);
      ResultsFormat::format_avavgprectime(&iters_times,&results_stream);
      // Now doing the linear solver time.
      for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
      {
        ResultsFormat::format_solvertime(intimestep,&iters_times,&results_stream);
      }
      
      ResultsFormat::format_avgsolvertime(&iters_times,&results_stream);
      ResultsFormat::format_avavgsolvertime(&iters_times,&results_stream); 
  
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
  MPI_Helpers::finalize();
#endif


} // end_of_main


