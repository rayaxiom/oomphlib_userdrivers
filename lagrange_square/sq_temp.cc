//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
//LIC//=====================================================================
//Driver for 2D tilted rectangle

#include <fenv.h>
#include <sstream>
#include <iomanip>
#include <ios>

// Oomph-lib includes
#include "generic.h"
#include "navier_stokes.h"

// The 2D mesh
#include "meshes/simple_rectangular_quadmesh.h"
#include "meshes/rectangular_quadmesh.h"

// My own header
#include "./../rayheader.h"

using namespace std;

using namespace oomph;

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Variables
{

 // Problem dimension
 static const unsigned Dim = 2;

 // Min and max x value respectively.
 static const double X_min = 0.0;
 static const double X_max = 1.0;

 // Min and max y value respectively.
 static const double Y_min = 0.0;
 static const double Y_max = 1.0;

 // The domain length in the x and y direction respectively.
 static const double Lx = X_max - X_min;
 static const double Ly = Y_max - Y_min;

 /// Reynolds number
 double Re = 100.0;

 // Tilting angle of the domain with the x-axis
 double Ang_deg = 30.0;
 double Ang_rad = -1.0;

 // Number of elements in 1D
 unsigned Noel = 4;

 // Use LSC preconditioner for the Navier-Stokes block?
 bool Use_lsc = false;

 // Use Boomer AMG for the momentum block?
 bool Use_amg_for_f = false;

 // Use Boomer AMG for the pressure block?
 bool Use_amg_for_p = false;

 // Convert degrees to radians
 inline double degtorad(const double& ang_deg)
 {
   return ang_deg * (MathematicalConstants::Pi / 180.0);
 }

 /// Storage for number of iterations during Newton steps 
 Vector<unsigned> Iterations;

} // namespace Global_Variables



namespace oomph
{
//========================================================================
/// \short A Sloping Mesh  class.
///
/// derived from RectangularQuadMesh:
/// the same mesh rotated with an angle phi
//========================================================================
 template<class ELEMENT>
 class SlopingQuadMesh : public RectangularQuadMesh<ELEMENT>
 {
 public:

  /// Constructor.
  SlopingQuadMesh(const unsigned& nx, const unsigned& ny,
                  const double& lx,  const double& ly, const double& phi ) :
   RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly)
   {
    // Find out how many nodes there are
    unsigned n_node=this->nnode();

    // Loop over all nodes
    for (unsigned n=0;n<n_node;n++)
     {
      // Pointer to node:
      Node* nod_pt=this->node_pt(n);

      // Get the x/y coordinates
      double x=nod_pt->x(0);
      double y=nod_pt->x(1);

      // Set new nodal coordinates
      nod_pt->x(0)=x*cos(phi)-y*sin(phi);
      nod_pt->x(1)=x*sin(phi)+y*cos(phi);
     }
   }
 };
} // end of namespace oomph


//===start_of_problem_class=============================================
//======================================================================

template<class ELEMENT>
class TiltedCavityProblem : public Problem
{
public:

 /// \short Constructor: Pass number of elements in x and y directions and
 /// lengths
 TiltedCavityProblem();

 /// Update before solve is empty
 void actions_before_newton_solve()
 {
   // Alias the namespace for convenience
   namespace GV = Global_Variables;
   GV::Iterations.clear();

   namespace NSPP = NavierStokesProblemParameters;
   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
   {
    // Initialise counters for each newton solve.
    Doc_linear_solver_info_pt->setup_new_time_step();
   }
 }

 /// \short Update after solve is empty
 void actions_after_newton_solve()
 {
 }

 void actions_after_newton_step()
 {
    // Alias the namespace for convenience
    namespace GV = Global_Variables;

    unsigned iters = 0;
    // Get the iteration counts
#ifdef PARANOID
    IterativeLinearSolver* iterative_solver_pt
      = dynamic_cast<IterativeLinearSolver*>
      (this->linear_solver_pt());
    if(iterative_solver_pt == 0)
    {
      std::ostringstream error_message;
      error_message << "Cannot cast the solver pointer." << std::endl;

      throw OomphLibError(error_message.str(),
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      iters = iterative_solver_pt->iterations();
      GV::Iterations.push_back(iters);
    }
#else
    iters = static_cast<IterativeLinearSolver*>
       (this->linear_solver_pt())->iterations();
    GV::Iterations.push_back(iters);
#endif


   namespace NSPP = NavierStokesProblemParameters;
   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
   {
     NSPP::doc_iter_times(this,Doc_linear_solver_info_pt);
   }
 }

 void actions_before_distribute()
 {
   namespace NSPP = NavierStokesProblemParameters;
   namespace SL = SquareLagrange;

   if(NSPP::Distribute_problem)
   {
     if(NSPP::Prob_id == SL::PID_SQ_PO)
     {
       GenericProblemSetup::delete_flux_elements(Surface_mesh_P_pt);

      rebuild_global_mesh();
     }
     else
     {
   std::ostringstream err_msg;
   err_msg << "Please set up the distributed bit for problem id: "
           << NSPP::Prob_id << ".\n"
           << std::endl;

   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
     }
   }
 }

 void actions_after_distribute()
 {
   namespace NSPP = NavierStokesProblemParameters;
   namespace SL = SquareLagrange;

   if(NSPP::Distribute_problem)
   {
   if(NSPP::Prob_id == SL::PID_SQ_PO)
   {
     create_parall_outflow_lagrange_elements(1,
         Bulk_mesh_pt,Surface_mesh_P_pt);
     rebuild_global_mesh();
   }
   else
   {
   std::ostringstream err_msg;
   err_msg << "Please set up the distributed bit for problem id: "
           << NSPP::Prob_id << ".\n"
           << std::endl;

   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
   }
   }
 }

 /// Doc the solution
 void doc_solution();

 /// \short Create lagrange elements on boundary b of the Mesh pointed
 /// to by bulk_mesh_pt and add them to the Mesh object pointed to by
 /// surface_mesh_pt
 void create_parall_outflow_lagrange_elements(const unsigned &b,
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);

 void create_impenetrable_lagrange_elements(const unsigned &b,
                                            Mesh* const &bulk_mesh_pt,
                                            Mesh* const &surface_mesh_pt);

 void set_inflow_BC(const unsigned &b,
                    Mesh* const &bulk_mesh_pt);
 void set_nonslip_BC(const unsigned &b,
                     Mesh* const &bulk_mesh_pt);

private:

 void set_mesh_bc_for_SqPo();
 void set_mesh_bc_for_SqTf();
 void set_mesh_bc_for_SqVa();

 /// Pointer to the "bulk" mesh
 //SlopingQuadMesh<ELEMENT>* Bulk_mesh_pt;
 Mesh* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_T_pt;
 Mesh* Surface_mesh_P_pt;

 // Preconditioner
 Preconditioner* Prec_pt;
 // Solver
 IterativeLinearSolver* Solver_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;

 unsigned Right_bound;
 unsigned Left_bound;
 unsigned Top_bound;
 unsigned Bottom_bound;

};




//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT> // rrrback - changed here.
TiltedCavityProblem<ELEMENT>::TiltedCavityProblem()
{
  // Alias the namespace for convenience
  namespace NSPP = NavierStokesProblemParameters;
  namespace LPH = LagrangianPreconditionerHelpers;
  namespace SL = SquareLagrange;

  Bottom_bound = 0;
  Right_bound = 1;
  Top_bound = 2;
  Left_bound = 3;
  
  Doc_linear_solver_info_pt = NSPP::Doc_linear_solver_info_pt;

      /// Setup the mesh
    // # of elements in x-direction
    const unsigned nx=SL::Noel;
    
    // # of elements in y-direction
    const unsigned ny=SL::Noel;
    
    // Domain length in x-direction
    const double lx=SL::Lx;

    // Domain length in y-direction
    const double ly=SL::Ly;
  // First we set the mesh.
  if((NSPP::Prob_id == SL::PID_SQ_TMP) ||
     (NSPP::Prob_id == SL::PID_SQ_PO)  ||
     (NSPP::Prob_id == SL::PID_SQ_TF)  ||
     (NSPP::Prob_id == SL::PID_SQ_TFPO) )
  {
    // This is the tilted cavity mesh.
    Bulk_mesh_pt =
      new SlopingQuadMesh<ELEMENT>(nx,ny,lx,ly,SL::Ang);
  }
  else if (NSPP::Prob_id == SL::PID_SQ_VA)
  {
    Bulk_mesh_pt = new RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly);
  }
  else
  {
   std::ostringstream err_msg;
   err_msg << "There is no mesh for the problem ID: "
           << NSPP::Prob_id << ".\n"
           << std::endl;

   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
  }

  // Set the boundary conditions
  if(NSPP::Prob_id == SL::PID_SQ_PO)
  {
    set_mesh_bc_for_SqPo();
  }
  else if(NSPP::Prob_id == SL::PID_SQ_TF)
  {
    set_mesh_bc_for_SqTf();
  }
  else if(NSPP::Prob_id == SL::PID_SQ_VA)
  {
    set_mesh_bc_for_SqVa();
  }
  else
  {
    std::ostringstream err_msg;
   err_msg << "There are no boundary conditions configured for prob_id: "
           << NSPP::Prob_id << ".\n"
           << std::endl;

   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
  }




// // Top boundary is slip.
// current_bound = 2;
// num_nod= mesh_pt()->nboundary_node(current_bound);
// for(unsigned inod=0;inod<num_nod;inod++)
// {
//   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);
//
//   if(!nod_pt->is_on_boundary(3))
//   {
//     nod_pt->unpin(0);
//     nod_pt->pin(1);
//
//     nod_pt->set_value(1,0.0);
//   }
// }


 
 //set_nonslip_BC(0,Bulk_mesh_pt);
// set_nonslip_BC(2,Bulk_mesh_pt);

// set_inflow_BC(if_b,Bulk_mesh_pt);
 
 //Complete the problem setup to make the elements fully functional

 //Loop over the elements
 unsigned n_el = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &NSPP::Rey;

  } // for(unsigned e=0;e<n_el;e++)

 //Assign equation numbers
 oomph_info << "\n equation numbers : "<< assign_eqn_numbers() << std::endl;
 

 // Only do this bit if we do NOT have a direct solver.
 if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
 {
   // Create the vector of mesh pointers!
   Vector<Mesh*> mesh_pt;
   if(NSPP::Prob_id == SL::PID_SQ_PO)
   {
     mesh_pt.resize(2,0);
     mesh_pt[0] = Bulk_mesh_pt;
     mesh_pt[1] = Surface_mesh_P_pt;
   }
   else if(NSPP::Prob_id == SL::PID_SQ_TF)
   {
     mesh_pt.resize(2,0);
     mesh_pt[0] = Bulk_mesh_pt;
     mesh_pt[1] = Surface_mesh_T_pt;
   }
   else if(NSPP::Prob_id == SL::PID_SQ_VA)
   {
     mesh_pt.resize(1,0);
     mesh_pt[0] = Bulk_mesh_pt;
   }

   // Quick check that the correct preconditioner is chosen.
   if((NSPP::Prob_id == SL::PID_SQ_VA) && 
       !CommandLineArgs::command_line_flag_has_been_set("--lsc_only"))
   {
     std::ostringstream err_msg;
     err_msg << "You have requested Vanilla Navier-Stokes problem,\n"
       << "NSPP::Prob_id is " << NSPP::Prob_id << "\n"
       << "But you have not set the flag --lsc_only.\n" 
       << "Please choose you preconditioner parameters again.\n"
       << std::endl;

     throw OomphLibError(err_msg.str(),
         OOMPH_CURRENT_FUNCTION,
         OOMPH_EXCEPTION_LOCATION);
   }

   LPH::Mesh_pt = mesh_pt;
   LPH::Problem_pt = this;
   Prec_pt = LPH::get_preconditioner();
 }
 const double solver_tol = 1.0e-6;
 const double newton_tol = 1.0e-6;
 GenericProblemSetup::setup_solver(NSPP::Max_solver_iteration,
                                   solver_tol,newton_tol,
                                   NSPP::Solver_type,this,Prec_pt);
}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::doc_solution()
{

  namespace NSPP = NavierStokesProblemParameters;
  
  std::ofstream some_file;
  std::stringstream filename;
  filename << NSPP::Soln_dir_str<<"/"<<NSPP::Label_str<<".dat";

  // Number of plot points
  unsigned npts=5;

  // Output solution
  some_file.open(filename.str().c_str());
  Bulk_mesh_pt->output(some_file,npts);
  some_file.close();
}

//============RAYRAY===========
/// RAYRAY
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::set_mesh_bc_for_SqPo()
{
  // Alias the namespace for convenience
  namespace SL = SquareLagrange;

  // Assign the boundaries:
  //             2 non slip
  //         ----------
  //         |        |
  // 3 Inflow|        |1 P.O.
  //         |        |
  //         ----------
  //             0 non slip
//  unsigned if_b = 3; // inflow
  const unsigned po_b = 1; // parallel outflow

  // Create a "surface mesh" that will contain only
  // ImposeParallelOutflowElements in boundary 1
  // The constructor just creates the mesh without
  // giving it any elements, nodes, etc.
  Surface_mesh_P_pt = new Mesh;

  // Create ImposeParallelOutflowElement from all elements that are
  // adjacent to the Neumann boundary.
  create_parall_outflow_lagrange_elements(po_b,
                                          Bulk_mesh_pt,Surface_mesh_P_pt);

  // Add the two meshes to the problem.
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_P_pt);
  
  // combine all sub-meshes into a single mesh.
  build_global_mesh();
  const unsigned num_bound = mesh_pt()->nboundary();

  // Set the boundary conditions for this problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here.
  for(unsigned ibound=0;ibound<num_bound;ibound++)
  { 
    if(ibound != po_b)
    {
      unsigned num_nod=mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
      {
        // Get node
        Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
 
        nod_pt->pin(0);
        nod_pt->pin(1);
        
        nod_pt->set_value(0,0);
        nod_pt->set_value(1,0);
 
      }
    }
  }

 // Which boundary are we dealing with?
 unsigned current_bound;
 
 // The number of nodes on a boundary.
 unsigned num_nod;

 // Inflow is at boundary 3
 current_bound = 3;
 num_nod= mesh_pt()->nboundary_node(current_bound);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);

   // Pin both velocity components
   nod_pt->pin(0);
   nod_pt->pin(1);

   // Get the x and y cartesian coordinates
   double x0=nod_pt->x(0);
   double x1=nod_pt->x(1);

   // Tilt x1 by -SL::Ang, this will give us the original coordinate.
   double x1_old = x0*sin(-SL::Ang) + x1*cos(-SL::Ang); 

   // Now calculate the parabolic inflow at this point.
   double u0_old = (x1_old - SL::Y_min)*(SL::Y_max - x1_old);
   
   // Now apply the rotation to u0_old, using rotation matrices.
   // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
   // velocity in the x direction only. There is no velocity
   // in the y direction.
   double u0=u0_old*cos(SL::Ang);
   double u1=u0_old*sin(SL::Ang);

   nod_pt->set_value(0,u0);
   nod_pt->set_value(1,u1);
 }
} // set_mesh_bc_for_SqPo

//============RAYRAY===========
/// RAYRAY
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::set_mesh_bc_for_SqTf()
{
  // Alias the namespace for convenience
  namespace SL = SquareLagrange;

  // Assign the boundaries:
  //             2 slip bc
  //         ----------
  //         |        |
  // 3 Inflow|        |1 Imposed outflow
  //         |        |
  //         ----------
  //             0 non slip
//  const unsigned if_b = 3; // inflow
//  const unsigned po_b = 1; // parallel outflow
  const unsigned slip_b = 2;

  // Create a "surface mesh" that will contain only
  // ImposeParallelOutflowElements in boundary 1
  // The constructor just creates the mesh without
  // giving it any elements, nodes, etc.
  Surface_mesh_T_pt = new Mesh;

  // Create ImposeParallelOutflowElement from all elements that are
  // adjacent to the Neumann boundary.
  create_impenetrable_lagrange_elements(slip_b,Bulk_mesh_pt,
                                        Surface_mesh_T_pt);
//  create_parall_outflow_lagrange_elements(slip_b,
//                                          Bulk_mesh_pt,Surface_mesh_T_pt);

  // Add the two meshes to the problem.
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_T_pt);
  
  // combine all sub-meshes into a single mesh.
  build_global_mesh();
//  const unsigned num_bound = mesh_pt()->nboundary();


 // Which boundary are we dealing with?
 unsigned current_bound;
 
 // The number of nodes on a boundary.
 unsigned num_nod;

 // Inflow is at boundary 3
 current_bound = 3;
 num_nod= mesh_pt()->nboundary_node(current_bound);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);

   // Pin both velocity components
   nod_pt->pin(0);
   nod_pt->pin(1);

   // Get the x and y cartesian coordinates
   double x0=nod_pt->x(0);
   double x1=nod_pt->x(1);

   // Tilt x1 by -SL::Ang, this will give us the original coordinate.
   double x1_old = x0*sin(-SL::Ang) + x1*cos(-SL::Ang);

   // Now calculate the parabolic inflow at this point.
   double u0_old = (x1_old - SL::Y_min)*(2*SL::Y_max - x1_old);
   
   // Now apply the rotation to u0_old, using rotation matrices.
   // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
   // velocity in the x direction only. There is no velocity
   // in the y direction.
   double u0=u0_old*cos(SL::Ang);
   double u1=u0_old*sin(SL::Ang);

   nod_pt->set_value(0,u0);
   nod_pt->set_value(1,u1);
 }

 // Now do the outflow.
 // This is on boundary 1.
 current_bound = 1;
 num_nod= mesh_pt()->nboundary_node(current_bound);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);

   // Pin both velocity components
   nod_pt->pin(0);
   nod_pt->pin(1);

   // Get the x and y cartesian coordinates
   double x0=nod_pt->x(0);
   double x1=nod_pt->x(1);

   // Tilt x1 by -SL::Ang, this will give us the original coordinate.
   double x1_old = x0*sin(-SL::Ang) + x1*cos(-SL::Ang);

   // Now calculate the parabolic inflow at this point.
   double u0_old = (x1_old - SL::Y_min)*(2*SL::Y_max - x1_old);
   
   // Now apply the rotation to u0_old, using rotation matrices.
   // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
   // velocity in the x direction only. There is no velocity
   // in the y direction.
   double u0=u0_old*cos(SL::Ang);
   double u1=u0_old*sin(SL::Ang);

//   nod_pt->unpin(0);
   nod_pt->set_value(0,u0);
   nod_pt->set_value(1,u1);
 }

 // Now we do the bottom boundary, this is boundary 0
 current_bound = 0;
 num_nod= mesh_pt()->nboundary_node(current_bound);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);

   // Pin both velocity components
   nod_pt->pin(0);
   nod_pt->pin(1);

   nod_pt->set_value(0,0.0);
   nod_pt->set_value(1,0.0);
 }

} // set_mesh_bc_for_SqPo

template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::set_mesh_bc_for_SqVa()
{
  // Alias the namespace for convenience
  namespace SL = SquareLagrange;

  // Assign the boundaries:
  //             2 non slip
  //         ----------
  //         |        |
  // 3 Inflow|        |1 P.O.
  //         |        |
  //         ----------
  //             0 non slip
//  unsigned if_b = 3; // inflow
  unsigned po_b = 1; // parallel outflow

  // Create a "surface mesh" that will contain only
  // ImposeParallelOutflowElements in boundary 1
  // The constructor just creates the mesh without
  // giving it any elements, nodes, etc.
//  Surface_mesh_P_pt = new Mesh;

  // Create ImposeParallelOutflowElement from all elements that are
  // adjacent to the Neumann boundary.
//  create_parall_outflow_lagrange_elements(po_b,
//                                          Bulk_mesh_pt,Surface_mesh_P_pt);

  // Add the two meshes to the problem.
  add_sub_mesh(Bulk_mesh_pt);
//  add_sub_mesh(Surface_mesh_P_pt);
  
  // combine all sub-meshes into a single mesh.
  build_global_mesh();
  const unsigned num_bound = mesh_pt()->nboundary();

  // Leave the x velocity on po_b to unpinned.
  {
    unsigned num_nod = mesh_pt()->nboundary_node(po_b);
    for (unsigned inod = 0; inod < num_nod; inod++) 
    {
      Node* nod_pt = mesh_pt()->boundary_node_pt(po_b,inod);
      // Unpin x
      nod_pt->unpin(0);

      // Pin y
      nod_pt->pin(1);
      // Set y value to zero.
      nod_pt->set_value(1,0);
    }
  }


  // Set the boundary conditions for this problem: All nodes are
  // free by default -- just pin the ones that have Dirichlet conditions
  // here.
  for(unsigned ibound=0;ibound<num_bound;ibound++)
  { 
    if(ibound != po_b)
    {
      unsigned num_nod=mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
      {
        // Get node
        Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
 
        nod_pt->pin(0);
        nod_pt->pin(1);
        
        nod_pt->set_value(0,0);
        nod_pt->set_value(1,0);
 
      }
    }
  }

 // Which boundary are we dealing with?
 unsigned current_bound;
 
 // The number of nodes on a boundary.
 unsigned num_nod;

 // Inflow is at boundary 3
 current_bound = 3;
 num_nod=mesh_pt()->nboundary_node(current_bound);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);

   // Pin both velocity components
   nod_pt->pin(0);
   nod_pt->pin(1);

   // Get the y cartesian coordinates
   double x1=nod_pt->x(1);

   // Now calculate the parabolic inflow at this point.
   double u0 = (x1 - SL::Y_min)*(SL::Y_max - x1);
   
   double u1=0.0;

   nod_pt->set_value(0,u0);
   nod_pt->set_value(1,u1);
 }
} // set_mesh_bc_for_SqPo

//============RAYRAY===========
/// RAYRAY
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
set_nonslip_BC(const unsigned &b,
               Mesh* const &bulk_mesh_pt)
{
  unsigned num_nod = bulk_mesh_pt->nboundary_node(b);
  unsigned dim = bulk_mesh_pt->finite_element_pt(0)->node_pt(0)->ndim();
   
  for(unsigned inod=0;inod<num_nod;inod++)
   {
    Node* nod_pt=bulk_mesh_pt->boundary_node_pt(b,inod);
    
    // pin all velocity components and set the value to zero.
    for (unsigned velo_i = 0; velo_i < dim; velo_i++) 
    {
      nod_pt->pin(velo_i);
      nod_pt->set_value(velo_i,0);
    }
   }
}

//============RAYRAY===========
/// RAYRAY
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
set_inflow_BC(const unsigned &b,
              Mesh* const &bulk_mesh_pt)
{

 // Alias the namespace for convenience
 namespace SL = SquareLagrange;

 // Check that the dimension is correct.
#ifdef PARANOID
  unsigned dim = bulk_mesh_pt->finite_element_pt(0)->node_pt(0)->ndim();
  if(dim != 2)
   {
     std::ostringstream err_msg;
     err_msg << "Inflow implemented for dim = 2 only." << std::endl;

     throw OomphLibError(err_msg.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
   }
#endif

  unsigned num_nod = bulk_mesh_pt->nboundary_node(b);
   
  for(unsigned inod=0;inod<num_nod;inod++)
   {
    Node* nod_pt=bulk_mesh_pt->boundary_node_pt(b,inod);

    // Pin both velocity components
    nod_pt->pin(0);
    nod_pt->pin(1);

    // Get the x and y cartesian coordinates.
    double x0 = nod_pt->x(0);
    double x1 = nod_pt->x(1);

    // Tilt x1 back the coordinate so we get the original coordinate.
    double x1_old = x0*sin(-SL::Ang) + x1*cos(-SL::Ang);
    
    // Now calculate the parabolic inflow at this point
    //double u0_old = (x1_old - SL::Y_min)*(SL::Y_max - x1_old);
    double u0_old = (x1_old - SL::Y_min)*(2.0 - x1_old);

    // Now apply the rotation to u0_old, using rotation matrices.
    // with x = u0_old and y = 0, i.e. R*[u;0] since we have the
    // velocity in the x direction only. There is no velocity
    // in the y direction.
    double u0=u0_old*cos(SL::Ang);
    double u1=u0_old*sin(SL::Ang);
    
    nod_pt->set_value(0,u0);
    nod_pt->set_value(1,u1); 
   }
}


//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
create_parall_outflow_lagrange_elements(const unsigned &b,
                                        Mesh* const &bulk_mesh_pt,
                                        Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));

   // What is the index of the face of element e along boundary b?
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding impose_impenetrability_element
   ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new
    ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
                                          face_index,0);


   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 0 or 2?
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(2)))
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

//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
create_impenetrable_lagrange_elements(const unsigned &b,
                                        Mesh* const &bulk_mesh_pt,
                                        Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));

   //What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding impose_impenetrability_element
   ImposeImpenetrabilityElement<ELEMENT>* flux_element_pt = new
    ImposeImpenetrabilityElement<ELEMENT>(bulk_elem_pt,
                                          face_index);
//   ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new
//    ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
//                                          face_index);

   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 3 or 1?
     if ((nod_pt->is_on_boundary(3))||(nod_pt->is_on_boundary(1)))
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
  
  namespace NSPP = NavierStokesProblemParameters;
  namespace LPH = LagrangianPreconditionerHelpers;
  namespace SL = SquareLagrange;

  // Create a custom label for the F solve.
  std::string prec_string = "";
  if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
  {
  std::string w_string = "";
  if(LPH::W_solver == 0)
  {
    w_string = "We";
  }
  else if(LPH::W_solver == 1)
  {
    w_string = "Wcg";
  }
  else
  {
   std::ostringstream err_msg;
   err_msg << "Something went wrong setting up the W string."
           << "LPH::W_solver is " << LPH::W_solver << ".\n"
           << std::endl;

   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
  }

  prec_string += w_string;


  std::string ns_string = "";
  if(LPH::NS_solver == 0)
  {
    ns_string = "Ne";
  }
  else if(LPH::NS_solver == 1)
  {
    ns_string = "Nl";
  }
  else
  {
   std::ostringstream err_msg;
   err_msg << "Something went wrong setting up the NS string."
           << "LPH::NS_solver is " << LPH::NS_solver << ".\n"
           << std::endl;

   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
  }

  prec_string += ns_string;

  std::string f_string = "";
  std::string p_string = "";
  if(LPH::NS_solver == 1)
  {
    f_string += "F_";
    p_string += "P_";

    // First do F string. We either have Exact or AMG.
    if(LPH::F_solver == 0)
    {
      f_string += "e";
    }
    else if(LPH::F_solver == 96)
    {
      // Determine the cycle
      std::stringstream f_cycle_stream;
      f_cycle_stream << LPH::f_amg_iterations << "v" 
        << LPH::f_amg_smoother_iterations;

      f_string += f_cycle_stream.str();

      // Determine strength of dependence
      std::stringstream f_strn_stream;
      f_strn_stream << "Strn" << LPH::f_amg_strength;

      f_string += f_strn_stream.str();

      // Determine the coarsening strategy
      const int f_amg_coarsening = LPH::f_amg_coarsening;
      std::string f_amg_coarsening_str = "";
      switch (f_amg_coarsening)
      {
        case 0:
          f_amg_coarsening_str = "CLPJ";
          break;
        case 1:
          f_amg_coarsening_str = "cRS";
          break;
        case 3:
          f_amg_coarsening_str = "mRS";
          break;
        case 6:
          f_amg_coarsening_str = "Falgout";
          break;
        case 8:
          f_amg_coarsening_str = "PMIS";
          break;
        case 10:
          f_amg_coarsening_str = "HMIS";
          break;
        case 11:
          f_amg_coarsening_str = "oRS";
          break;
        default:            
          {
            std::ostringstream err_msg;
            err_msg << "Something wrong setting f coarsening string."
              << "LPH::f_amg_coarsening is " << LPH::f_amg_coarsening << ".\n"
              << std::endl;

            throw OomphLibError(err_msg.str(),
                OOMPH_CURRENT_FUNCTION,
                OOMPH_EXCEPTION_LOCATION);
          }
          break;
      } // Switch f_amg_coarsening.

      f_string += f_amg_coarsening_str;

      // Now do the smoother string.
      const int f_amg_simple_smoother = LPH::f_amg_simple_smoother;
      const int f_amg_complex_smoother = LPH::f_amg_complex_smoother;
      std::string f_amg_smoother_str = "";

      if(f_amg_simple_smoother >= 0)
      {
        if(f_amg_simple_smoother == 0)
        {
          f_amg_smoother_str = "Jac";

          // If it is Jacobi, we have the damping!
          std::ostringstream dampstream;
          dampstream << LPH::f_amg_damping;
          f_amg_smoother_str += dampstream.str();
        }
        else if(f_amg_simple_smoother == 1)
        {
          f_amg_smoother_str = "Gs";
        }
        else if(f_amg_simple_smoother == 2)
        {
          f_amg_smoother_str = "Gspinter";
        }
        else if(f_amg_simple_smoother == 3)
        {
          f_amg_smoother_str = "SORfs";
        }
        else if(f_amg_simple_smoother == 4)
        {
          f_amg_smoother_str = "SORbs";
        }
        else if(f_amg_simple_smoother == 6)
        {
          f_amg_smoother_str = "SSOR";
        }
        else
        {
          std::ostringstream err_msg;
          err_msg << "Something wrong setting f smoother string."
            << "LPH::f_amg_simple_smoother is " 
            << LPH::f_amg_simple_smoother << ".\n"
            << std::endl;

          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
      }
      else if(f_amg_complex_smoother >= 0)
      {
        if(f_amg_complex_smoother == 6)
        {
          f_amg_smoother_str = "Schwarz";
        }
        else if(f_amg_complex_smoother == 7)
        {
          f_amg_smoother_str = "Pilut";
        }
        else if(f_amg_complex_smoother == 8)
        {
          f_amg_smoother_str = "ParaSails";
        }
        else if(f_amg_complex_smoother == 9)
        {
          f_amg_smoother_str = "Euclid";
        }
        else
        {
          std::ostringstream err_msg;
          err_msg << "Something wrong setting f COMPLEX smoother string.\n"
            << std::endl;

          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
      }
      else
      {
        std::ostringstream err_msg;
        err_msg << "Something wrong setting f smoother string.\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      f_string += f_amg_smoother_str;
    }
    else if (LPH::F_solver == 8090)
    {
      f_string += "BlkDiagLU";
    }
    else if (LPH::F_solver == 8091)
    {
      f_string += "BlkUpperLU";
    }
    else if (LPH::F_solver == 8092)
    {
      f_string += "BlkLowerLU";
    }
    else if (LPH::F_solver == 8093)
    {
      f_string += "LU";
    }
    else if (LPH::F_solver == 9090)
    {
      f_string += "BlkDiagAMGSim";
    }
    else if (LPH::F_solver == 9091)
    {
      f_string += "BlkUpperAMGSim";
    }
    else if (LPH::F_solver == 9092)
    {
      f_string += "BlkLowerAMGSim";
    }
    else if (LPH::F_solver == 9093)
    {
      f_string += "AMGSim";
    }
    else if (LPH::F_solver == 9190)
    {
      f_string += "BlkDiagAMGStr";
    }
    else if (LPH::F_solver == 9191)
    {
      f_string += "BlkUpperAMGStr";
    }
    else if (LPH::F_solver == 9192)
    {
      f_string += "BlkLowerAMGStr";
    }
    else if (LPH::F_solver == 9193)
    {
      f_string += "AMGStr";
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Something went wrong setting up the F string."
        << "LPH::F_solver is " << LPH::F_solver << ".\n"
        << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Now do the pressure.
    if(LPH::P_solver == 0)
    {
      p_string += "e";
    }
    else if(LPH::P_solver == 1)
    {
      p_string += "1v22Strn0.25CLJPGs";
    }
    else if(LPH::P_solver == 13)
    {
      p_string += "1v22Strn0.7CLJPGs";
    }
    else if(LPH::P_solver == 96)
    {
      p_string += "Pcustom";
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Something went wrong setting up the P string."
        << "LPH::P_solver is " << LPH::P_solver << ".\n"
        << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    } //  P_solver

  } // if NS_solver == 1

  prec_string += f_string;
  prec_string += p_string;
  }
  else
  {
    prec_string = "DirectSolve";
  }


  std::string label = SL::prob_str()
                      + NSPP::create_label()
                      + prec_string
                      + SL::ang_deg_str() + SL::noel_str();
  
  std::replace(label.begin(), label.end(),
                '.','-');
  return label;
}

//===start_of_main======================================================
/// Driver code
//======================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);
#endif
  // Alias the namespace for convenience.
  namespace GV = Global_Variables;


  // Alias the namespace for convenience.
  namespace NSPP = NavierStokesProblemParameters;
  namespace LPH = LagrangianPreconditionerHelpers;
  namespace SL = SquareLagrange;

  // Problem dimension.
  const unsigned dim = 2;

  // Set up doc info - used to store information on solver and 
  // iteration time.
  DocLinearSolverInfo doc_linear_solver_info;
  // Again, pass this to the NSPP and LPH
  NSPP::Doc_linear_solver_info_pt = &doc_linear_solver_info;
  LPH::Doc_linear_solver_info_pt = &doc_linear_solver_info;

  // Set the Label_pt
  LPH::Label_str_pt = &NSPP::Label_str;
  LPH::Vis_pt = &NSPP::Vis;
  SL::Prob_id_pt = &NSPP::Prob_id;

  // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  NSPP::setup_commandline_flags();
  LPH::setup_commandline_flags();
  SL::setup_commandline_flags();

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////

  // dim = 2
  NSPP::generic_problem_setup(dim);
  LPH::generic_setup();
  SL::generic_setup();

  // Solve with Taylor-Hood element, set up problem
  TiltedCavityProblem< QTaylorHoodElement<dim> > problem;

  if(NavierStokesProblemParameters::Distribute_problem)
  {
    problem.distribute();
  }

  //////////////////////////////////////////////////////////////////////////

  // If the Reynolds number is not set, I assume that all the below are set:
  // --rey_start
  // --rey_incre
  // --rey_end
  //
  // This is meticulously checked above.
  // NOTE: This is still not done. I am doing the others first.
  if(!CommandLineArgs::command_line_flag_has_been_set("--rey"))
  {
    unsigned rey_increment = 0;

    for (NSPP::Rey = NSPP::Rey_start; 
        NSPP::Rey <= NSPP::Rey_end; NSPP::Rey += NSPP::Rey_incre)
    {
      // Setup the label. Used for doc solution and preconditioner.
      NSPP::Label_str = create_label();

      time_t rawtime;
      time(&rawtime);

      oomph_info << "RAYDOING: "
        << NSPP::Label_str
        << " on " << ctime(&rawtime) << std::endl;

      // Solve the problem
      problem.newton_solve();

      //Output solution
      if(NSPP::Doc_soln)
      {problem.doc_solution();}

      //////////////////////////////////////////////////////////////////////
      ////////// Outputting results ////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////

      if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
      {
      // Get the global oomph-lib communicator 
      const OomphCommunicator* const comm_pt 
        = MPI_Helpers::communicator_pt();

      // My rank and number of processors. 
      // This is used later for putting the data.
      const unsigned my_rank = comm_pt->my_rank();
      const unsigned nproc = comm_pt->nproc();

      // Variable to indicate if we want to output to a file or not.
      bool output_to_file = false;

      // The output file.
      std::ofstream outfile;

      // If we want to output to a file, we create the outfile.
      if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir"))
      {
        output_to_file = true;
        std::ostringstream filename_stream;
        filename_stream << NSPP::Itstime_dir_str<<"/"
          << NSPP::Label_str
          <<"NP"<<nproc<<"R"<<my_rank;
        outfile.open(filename_stream.str().c_str());
      }


      // Stringstream to hold the results. We do not output the results
      // (timing/iteration counts) as we get it since 
      // it will interlace with the
      // other processors and becomes hard to read.
      std::ostringstream results_stream;

      // Get the 3D vector which holds the iteration counts 
      // and timing results.
      Vector<Vector<Vector<double> > > iters_times
        = NSPP::Doc_linear_solver_info_pt->iterations_and_times();

      ResultsFormat::format_rayits(rey_increment,
                                   &iters_times,&results_stream);
      ResultsFormat::format_prectime(rey_increment,
                                     &iters_times,&results_stream);
      ResultsFormat::format_solvertime(rey_increment,
                                       &iters_times,&results_stream);

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

      rey_increment++;

    }
  }
  else
  {
    // Setup the label. Used for doc solution and preconditioner.
    // NSPP::Label_str 
    // = NSPP::create_label() + LPH::create_label()+SL::create_label();
    NSPP::Label_str = create_label();

    time_t rawtime;
    time(&rawtime);

    oomph_info << "RAYDOING: "
      << NSPP::Label_str
      << " on " << ctime(&rawtime) << std::endl;

    // Solve the problem
    problem.newton_solve();

    //Output solution
    if(NSPP::Doc_soln)
    {problem.doc_solution();}


    ////////////////////////////////////////////////////////////////////////
    //////////// Outputting results ////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
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
    if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir"))
    {
      output_to_file = true;
      std::ostringstream filename_stream;
      filename_stream << NSPP::Itstime_dir_str<<"/"
        << NSPP::Label_str
        <<"NP"<<nproc<<"R"<<my_rank;
      outfile.open(filename_stream.str().c_str());
    }

    // Stringstream to hold the results. We do not output the results
    // (timing/iteration counts) as we get it since 
    // it will interlace with the
    // other processors and becomes hard to read.
    std::ostringstream results_stream;

    // Get the 3D vector which holds the iteration counts 
    // and timing results.
    Vector<Vector<Vector<double> > > iters_times
      = NSPP::Doc_linear_solver_info_pt->iterations_and_times();

    // Since this is a steady state problem, there is only
    // one "time step", thus it is essentially a 2D vector 
    // (the outer-most vector is of size 1).

    // Loop over the time steps and output the iterations, 
    // prec setup time and linear solver time.
    unsigned ntimestep = iters_times.size();
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat::format_rayits(intimestep,&iters_times,&results_stream);
    }

    // Now doing the preconditioner setup time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat::format_prectime(intimestep,
          &iters_times,&results_stream);
    }

    // Now doing the linear solver time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat::format_solvertime(intimestep,
          &iters_times,&results_stream);
    }

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
  } // else do not loop reynolds

  // Print out the iteration counts
  const unsigned num_newton_steps = GV::Iterations.size();
  oomph_info << "RRRNS:\t"<<num_newton_steps << std::endl; 
  for (unsigned stepi = 0; stepi < num_newton_steps; stepi++) 
  {
    oomph_info << "RRRITS:\t" << GV::Iterations[stepi] << std::endl;
  }
  oomph_info << "RRR:\t" << std::endl; 
  

#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main
