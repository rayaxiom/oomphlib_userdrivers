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
#include "meshes/simple_cubic_tet_mesh.h"


#include "./../ray_preconditioner_creation.h"
#include "./../ray_navier_stokes_parameters.h"
#include "./../ray_general_problem_parameters.h"


using namespace std;

using namespace oomph;

namespace PH = PreconditionerHelpers;
namespace GPH = GeneralProblemHelpers;
namespace NSH = NavierStokesHelpers;

 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////




namespace Global_Parameters
{

  // Problem specific parameters

  unsigned Noel = 4;
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

    const double ux_scaling = -cos(MathematicalConstants::Pi*t)/2.0 + 0.51;
    return get_prescribed_inflow_for_quarter(y,z) * ux_scaling;
  } 


}

namespace GP = Global_Parameters;


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Test problem for Fp/PCD preconditioner
//====================================================================

class FpTestProblem : public Problem
{

public:


 /// Constructor
 FpTestProblem();

 
 /// Destructor: Cleanup
 ~FpTestProblem()
  {
//   delete Solver_pt;
//   delete Prec_pt;
//   delete P_matrix_preconditioner_pt;
//   delete F_matrix_preconditioner_pt;
  }


 
 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {

  } // end_of_actions_after_adapt
 

 /// Actions after Newton step record number of iterations
 void actions_after_newton_step() 
  {                               

  }  

 /// Update the after solve (empty)
 void actions_after_newton_solve(){}

 /// Update the problem specs before solve. 
 void actions_before_newton_solve()
 {
   {
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
           = GP::get_prescribed_inflow_for_quarter(y,z);

         nod_pt->set_value(0,ux);
         nod_pt->set_value(1,0.0);
         nod_pt->set_value(2,0.0);
       }
     }

   }
  
 } // end_of_actions_before_newton_solve

 /// Doc the solution
 void doc_solution(DocInfo& doc_info);


 /// Pointer to the "bulk" mesh
 Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

private:

 /// Solver
 IterativeLinearSolver* Solver_pt;

 /// Solver
 Preconditioner* Prec_pt;

 /// Inexact solver for P block
 Preconditioner* P_matrix_preconditioner_pt;

 /// Inexact solver for F block
 Preconditioner* F_matrix_preconditioner_pt;

 /// ID of driven boundary
 unsigned Driven_boundary;

 /// ID of inflow boundary
 unsigned Inflow_boundary;

 /// ID of outflow boundary
 unsigned Outflow_boundary;


 /// Pointer to the "bulk" mesh
 Mesh* Bulk_mesh_pt;

}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for DrivenCavity problem 
//========================================================================
FpTestProblem::FpTestProblem()
{ 
 
 
 const unsigned noel = GP::Noel;
 const double length = GP::Length;

 Bulk_mesh_pt = 
   new SimpleCubicMesh<QTaylorHoodElement<3> >(noel,noel,noel,
                                               length,length,length);

 Driven_boundary=0;
 Inflow_boundary=4;
 Outflow_boundary=2;

 // Create "surface mesh" that will contain only the prescribed-traction 
 // elements.

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);

 // Combine all submeshes into a single Mesh
 build_global_mesh();
 

   
 Prec_pt = PH::create_lsc_preconditioner(this,Bulk_mesh_pt,0,0);
 
   
   // Set the boundary conditions for this problem: All nodes are
   // free by default -- just pin the ones that have Dirichlet conditions
   // here. 
   unsigned num_bound = Bulk_mesh_pt->nboundary();
   for(unsigned ibound=0;ibound<num_bound;ibound++)
    {
     unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Loop over values (u, v and w velocities)
       for (unsigned i=0;i<3;i++)
        {
         Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(i); 
        }
      }
    } // end loop over boundaries

 

 // In/outflow bcs
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
         //if (!(nod_pt->is_on_boundary(0)))
          {
           if ((nod_pt->x(1)<0.5)&&nod_pt->x(2)<0.5) nod_pt->unpin(0);
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
   el_pt->re_pt() = &NSH::Rey;
  } // end loop over elements
 
 

 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 


 const double solver_tol = 1.0e-6;
 const double newton_tol = 1.0e-6;

   GPH::setup_solver(GPH::Max_solver_iteration,solver_tol,newton_tol,
                    this,Prec_pt);

  Solver_pt = GPH::Iterative_linear_solver_pt;
} // end_of_constructor



inline void problem_specific_setup_commandline_flags()
{
    CommandLineArgs::specify_command_line_flag("--noel", 
        &GP::Noel);

        CommandLineArgs::specify_command_line_flag("--prob_id", 
        &GP::Prob_id);
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


  DocLinearSolverInfo doc_linear_solver_info;

  GPH::Doc_linear_solver_info_pt = &doc_linear_solver_info;
  NSH::Dim = 3;


  // Store commandline arguments
  CommandLineArgs::setup(argc,argv);

  GPH::setup_commandline_flags();
  NSH::setup_commandline_flags();
  PH::setup_commandline_flags();

  problem_specific_setup_commandline_flags();


  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();


  GPH::generic_setup();
  NSH::generic_setup();
  PH::generic_setup();
  problem_specific_generic_setup();


  // Build the problem 
  FpTestProblem problem;

  // Solve the problem 
  problem.newton_solve();


#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif


} // end_of_main


