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

// My own header
#include "./../rayheader.h"
using namespace std;

using namespace oomph;

// Rx = [1 0   0
//       0 cx -sx
//       0 sx  cx];
//
// Ry = [cy 0 sy
//       0  1 0
//      -sy 0 cy];
//
// Rz = [cz -sz 0
//       sz  cz  0
//       0   0  1];
//
// R = Rz*Ry*Rx
//
// Rx = x_new
void rotate_forward(const double& x, const double& y, const double z,
    const double& phix, const double& phiy, const double& phiz,
    Vector<double>& x_new)
{
  x_new.resize(3,0);
  x_new[0] = (cos(phiy)*cos(phiz))*x 
    + (cos(phiz)*sin(phix)*sin(phiy) - cos(phix)*sin(phiz))*y 
    + (sin(phix)*sin(phiz) + cos(phix)*cos(phiz)*sin(phiy))*z;

  x_new[1] = (cos(phiy)*sin(phiz))*x 
    + (cos(phix)*cos(phiz) + sin(phix)*sin(phiy)*sin(phiz))*y 
    + (cos(phix)*sin(phiy)*sin(phiz) - cos(phiz)*sin(phix))*z;

  x_new[2] = (-sin(phiy))*x 
    + (cos(phiy)*sin(phix))*y
    + ( cos(phix)*cos(phiy))*z;
}
// Rx = [1 0   0
//       0 cx -sx
//       0 sx  cx];
//
// Ry = [cy 0 sy
//       0  1 0
//      -sy 0 cy];
//
// Rz = [cz -sz 0
//       sz  cz  0
//       0   0  1];
//
// R = Rx*Ry*Rz (note the ordering)
//
// Rx = x_new
void rotate_backward(const double& x, const double& y, const double z,
    const double& phix, const double& phiy, const double& phiz,
    Vector<double>& x_new)
{
  x_new.resize(3,0);

  x_new[0] = (cos(phiy)*cos(phiz))*x 
    -(cos(phiy)*sin(phiz))*y 
    + (sin(phiy))*z;
  x_new[1] = (cos(phix)*sin(phiz) + cos(phiz)*sin(phix)*sin(phiy))*x 
    + (cos(phix)*cos(phiz) - sin(phix)*sin(phiy)*sin(phiz))*y 
    -(cos(phiy)*sin(phix))*z;
  x_new[2] = (sin(phix)*sin(phiz) - cos(phix)*cos(phiz)*sin(phiy))*x 
    + (cos(phiz)*sin(phix) + cos(phix)*sin(phiy)*sin(phiz))*y 
    + (cos(phix)*cos(phiy))*z;
}

namespace oomph
{
//========================================================================
/// \short A Sloping Mesh  class.
///
/// derived from RectangularQuadMesh:
/// the same mesh rotated with an angle phi
//========================================================================
 template<class ELEMENT>
 class SlopingCubicMesh : public SimpleCubicMesh<ELEMENT>
 {
 public:

  /// Constructor.
  SlopingCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                   const double& lx,  const double& ly, const double& lz,
                   const double& phix, const double& phiy, const double& phiz,
                   TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper) :
   SimpleCubicMesh<ELEMENT>(nx,ny,nz,lx,ly,lz,time_stepper_pt)
   {
    // Find out how many nodes there are
    unsigned n_node=this->nnode();

    // Loop over all nodes
    for (unsigned n=0;n<n_node;n++)
     {
      // Pointer to node:
      Node* nod_pt=this->node_pt(n);

      // Get the x/y coordinates
      const double x=nod_pt->x(0);
      const double y=nod_pt->x(1);
      const double z=nod_pt->x(2);

      Vector<double> x_new;
      rotate_forward(x,y,z,phix,phiy,phiz,x_new);

      // Set new nodal coordinates
      nod_pt->x(0)=x_new[0];
      nod_pt->x(1)=x_new[1];
      nod_pt->x(2)=x_new[2];
     }
   }
 };
} // end of namespace oomph

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

#ifdef OOMPH_HAS_HYPRE

//=============================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//=============================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
 Preconditioner* set_hypre_preconditioner()
 {
  return new HyprePreconditioner;
 }
}

#endif

 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Variables
{
 double Time_start = 0.0;
 double Time_end = 1.0;
 double Delta_t = 0.0;

 double Ang_deg = 45;

 double Ang = Ang_deg * (MathematicalConstants::Pi / 180.0);

 /// Reynolds number
 double Re=50.0;

 /// Storage for number of iterations during Newton steps 
 Vector<unsigned> Iterations;

 /// Storage for linear solver times during Newton steps 
 Vector<double> Linear_solver_time;

// /// Traction at the outflow boundary
// void prescribed_traction(const double& t,
//                          const Vector<double>& x,
//                          const Vector<double> &n,
//                          Vector<double>& traction)
// {
//  traction.resize(3);
//  traction[0]=1.0;
//  traction[1]=0.0;
//  traction[2]=0.0;
// } 

// /// Traction at the outflow boundary
// void get_prescribed_inflow(const double& t,
//                            const Vector<double>& x,
//                            Vector<double>& presc_inflow)
// {
//   // Get the x and y coordinates.
//   double y = x[1];
//   double z = x[2];
//
//   // For the velocity profile in the x direction.
//   // 1) First form the parabolic profile
//   double ux = 0.0;
//   if((y > 0.5)&&(z > 0.5))
//   {
////     // 2) Now make it move in time
////     const double trig_scaling = 0.025;
////     const double ux_scaling = 1.0 
////                               - cos(trig_scaling
////                                     *MathematicalConstants::Pi
////                                     *t);
////
////     ux = (y-0.5)*(1.0-y)*(z-0.5)*(1.0-z) * ux_scaling;
//     ux = (y-0.5)*(1.0-y)*(z-0.5)*(1.0-z);
//   }
//
//
//  presc_inflow.resize(3);
//  presc_inflow[0]=ux;
//  presc_inflow[1]=0.0;
//  presc_inflow[2]=0.0;
// } 

 void get_prescribed_inflow(const double& t,
                            const double& y,
                            const double& z,
                            double& ux)
 {
   // For the velocity profile in the x direction.
   // 1) First form the parabolic profile
   ux = 0.0;
   if((y > 0.5)&&(z > 0.5))
   {
     // 2) Now make it move in time
//     const double trig_scaling = 0.025;
//     const double ux_scaling = (1.0 
//                               - cos(trig_scaling
//                                     *MathematicalConstants::Pi
//                                     *t)) * 2.0;

     const double ux_scaling = t / Time_end;
     ux = (y-0.5)*(1.0-y)*(z-0.5)*(1.0-z) * ux_scaling;
//     ux = (y-0.5)*(1.0-y)*(z-0.5)*(1.0-z);
   }
 } 


 void get_prescribed_inflow_full(const double& t,
                            const double& y,
                            const double& z,
                            double& ux)
 {
   // For the velocity profile in the x direction.
   // 1) First form the parabolic profile
   ux = 0.0;
   {
     // 2) Now make it move in time
//     const double trig_scaling = 0.025;
//     const double ux_scaling = (1.0 
//                               - cos(trig_scaling
//                                     *MathematicalConstants::Pi
//                                     *t)) * 2.0;

     const double ux_scaling = t / Time_end;
     ux = (y-0.0)*(1.0-y)*(z-0.0)*(1.0-z) * ux_scaling;
//     ux = (y-0.5)*(1.0-y)*(z-0.5)*(1.0-z);
   }
 } 

} // end_of_namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Test problem for Fp/PCD preconditioner
//====================================================================
template<class ELEMENT>
class CubeProblem : public Problem
{

public:


 /// Constructor
 CubeProblem(const unsigned& n_element);

 
 /// Destructor: Cleanup
 ~CubeProblem()
  {
   delete Solver_pt;
   delete Prec_pt;
   delete P_matrix_preconditioner_pt;
   delete F_matrix_preconditioner_pt;
  }

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<NavierStokesEquations<3>*>(Bulk_mesh_pt->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure


//////////////////////////////

 void actions_before_implicit_timestep()
  {
    const double ang = Global_Variables::Ang;
    const double minus_ang = -ang;

   {
    // Inflow in upper half of inflow boundary
    const unsigned ibound=Inflow_boundary; 
    const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      const double x=nod_pt->x(0);
      const double y=nod_pt->x(1);
      const double z=nod_pt->x(2);

      const double time=time_pt()->time();

      // Rotate x y and z back
      Vector<double>x_new;
      rotate_backward(x,y,z,minus_ang,minus_ang,minus_ang,x_new);

      double ux = 0.0;

      Global_Variables::get_prescribed_inflow(time,x_new[1],x_new[2],ux);

      // Now rotate the velocity profile
      Vector<double>u_new;
      rotate_forward(ux,0,0,ang,ang,ang,u_new);

      nod_pt->set_value(0,u_new[0]);
      nod_pt->set_value(1,u_new[1]);
      nod_pt->set_value(2,u_new[2]);
     }
   }
  } // end of actions_before_implicit_timestep



////////////////////////////////


 void actions_before_distribute()
 {
   
   GenericProblemSetup::delete_flux_elements(Surface_mesh_pt);
   rebuild_global_mesh();
 }
 
 void actions_after_distribute()
 {
     create_parall_outflow_lagrange_elements(Outflow_boundary,
                                             Bulk_mesh_pt,
                                             Surface_mesh_pt);
     rebuild_global_mesh();
 }


 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   
   // Now set the first pressure dof in the first element to 0.0
//   if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt
 

 /// Actions after Newton step record number of iterations
 void actions_after_newton_step() 
  {                               
//   Global_Variables::Iterations.push_back(
//    dynamic_cast<IterativeLinearSolver*>
//    (this->linear_solver_pt())->iterations());
//   
//   Global_Variables::Linear_solver_time.push_back(
//    linear_solver_pt()->linear_solver_solution_time());
  }  

 /// Update the after solve (empty)
 void actions_after_newton_solve()
 {


 }

 /// Update the problem specs before solve. 
 void actions_before_newton_solve()
 {

//   {
//    // Inflow in upper half of inflow boundary
//    const unsigned ibound=Inflow_boundary; 
//    const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
//    for (unsigned inod=0;inod<num_nod;inod++)
//     {
//      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
//      const double y=nod_pt->x(1);
//      const double z=nod_pt->x(2);
//      const double time=1.0;
//
//      double ux = 0.0;
//
//      Global_Variables::get_prescribed_inflow_full(time,y,z,ux);
//
//      nod_pt->set_value(0,ux);
//      nod_pt->set_value(1,0.0);
//      nod_pt->set_value(2,0.0);
//     }
//   }


  // Initialise counter for iterations
//  Global_Variables::Iterations.clear();
//  Global_Variables::Linear_solver_time.clear();
  
 } // end_of_actions_before_newton_solve

 /// Run an unsteady simulation
 void unsteady_run(DocInfo& doc_info); 
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Create traction elements on outflow boundary
// template<class ELEMENT>
// void create_traction_elements();

// void create_inflow_traction_elements(const unsigned &b,
//                                      Mesh* const &bulk_mesh_pt,
//                                      Mesh* const &surface_mesh_pt);

 void create_parall_outflow_lagrange_elements(const unsigned &b,
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);
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
 
 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_pt;

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
CubeProblem<ELEMENT>::CubeProblem(const unsigned& n_el)
{ 
  // Alias the namespace for convenience
  namespace NSPP = NavierStokesProblemParameters;
  namespace LPH = LagrangianPreconditionerHelpers;
  namespace CL = CubeLagrange;

  // Assign boundary IDs defined in rayheader.h
  Left_boundary = CL::Left_boundary;
  Right_boundary = CL::Right_boundary;
  Front_boundary = CL::Front_boundary;
  Back_boundary = CL::Back_boundary;
  Bottom_boundary = CL::Bottom_boundary;
  Top_boundary = CL::Top_boundary; 

  Inflow_boundary=Left_boundary;
  Outflow_boundary=Right_boundary;

  add_time_stepper_pt(new BDF<2>);
  // Setup mesh

  // # of elements in x-direction
  unsigned n_x=n_el;

  // # of elements in y-direction
  unsigned n_y=n_el;

  // # of elements in z-direction
  unsigned n_z=n_el;

  // Domain length in x-direction
  double l_x=1.0;

  // Domain length in y-direction
  double l_y=1.0;

  // Domain length in y-direction
  double l_z=1.0;


//  Bulk_mesh_pt = 
//    new SimpleCubicMesh<ELEMENT>(n_x,n_y,n_z,l_x,l_y,l_z,time_stepper_pt());

  const double ang = Global_Variables::Ang;

  Bulk_mesh_pt = 
    new SlopingCubicMesh<ELEMENT>(n_x,n_y,n_z,
                                  l_x,l_y,l_z,
                                  ang,ang,ang,
                                  time_stepper_pt());

//  Bulk_mesh_pt = 
//    new SimpleCubicMesh<ELEMENT>(n_x,n_y,n_z,l_x,l_y,l_z);


  // Create "surface mesh" that will contain only the prescribed-traction 
  // elements.
  Surface_mesh_pt = new Mesh;


//   create_inflow_traction_elements(Inflow_boundary,
//                                   Bulk_mesh_pt,Surface_mesh_pt);

  create_parall_outflow_lagrange_elements(Outflow_boundary,
                                          Bulk_mesh_pt,
                                          Surface_mesh_pt);

  // Add the two sub meshes to the problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);

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
      // Loop over values (u, v and w velocities)
      for (unsigned i=0;i<3;i++)
      {
        Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries!


  {
    const double minus_ang = -Global_Variables::Ang;
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
//          if (!(nod_pt->is_on_boundary(0)) ||
//              !(nod_pt->is_on_boundary(1)) )
          {
            // We need to rotate the coordinates back
            Vector<double> x_old;
            const double x = nod_pt->x(0);
            const double y = nod_pt->x(1);
            const double z = nod_pt->x(2);
            
            rotate_backward(x,y,z,minus_ang,minus_ang,minus_ang,x_old);

            if ((x_old[1]<0.5)&&x_old[2]<0.5)
            {
              nod_pt->unpin(0);
              nod_pt->unpin(1);
              nod_pt->unpin(2);
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
    el_pt->re_pt() = &Global_Variables::Re;
    el_pt->re_st_pt() = &Global_Variables::Re;
  } // end loop over elements


  // Now set the first pressure value in element 0 to 0.0
  // if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);

  // Setup equation numbering scheme
  oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 



  LagrangeEnforcedflowPreconditioner* lgr_prec_pt
      = new LagrangeEnforcedflowPreconditioner;

  Vector<Mesh*> tmp_mesh_pt(2,0);
  tmp_mesh_pt[0] = Bulk_mesh_pt;
  tmp_mesh_pt[1] = Surface_mesh_pt;

  lgr_prec_pt->set_meshes(tmp_mesh_pt);



  // Build preconditoner
  NavierStokesSchurComplementPreconditioner* ns_prec_pt = 
    new NavierStokesSchurComplementPreconditioner(this);
//  Prec_pt=prec_pt;



  // By default, the Schur Complement Preconditioner uses SuperLU as
  // an exact preconditioner (i.e. a solver) for the
  // momentum and Schur complement blocks. 
  // Can overwrite this by passing pointers to 
  // other preconditioners that perform the (approximate)
  // solves of these blocks.


  // Create internal preconditioners used on Schur block
  P_matrix_preconditioner_pt=0;
  // if (use_hypre_for_pressure)
  {
#ifdef OOMPH_HAS_HYPRE

    oomph_info << "Using HYPRE for pressure block" << std::endl; 

    // Create preconditioner
    P_matrix_preconditioner_pt = new HyprePreconditioner;

    // Set parameters for use as preconditioner on Poisson-type problem
    Hypre_default_settings::set_defaults_for_3D_poisson_problem(
        static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt));

    // Use Hypre for the Schur complement block
    ns_prec_pt->set_p_preconditioner(P_matrix_preconditioner_pt);


    // Shut up!
    //   static_cast<HyprePreconditioner*>(P_matrix_preconditioner_pt)->
    //    disable_doc_time();

#endif
  }

  // Create block-diagonal preconditioner used on momentum block
  F_matrix_preconditioner_pt=0;   
  // if (use_block_diagonal_for_momentum)
  {

#ifdef OOMPH_HAS_HYPRE
    F_matrix_preconditioner_pt = new HyprePreconditioner;

    Hypre_default_settings::set_defaults_for_navier_stokes_momentum_block(
        static_cast<HyprePreconditioner*>(F_matrix_preconditioner_pt));
#endif       

    // Use Hypre for momentum block 
    ns_prec_pt->set_f_preconditioner(F_matrix_preconditioner_pt);
  }

  ns_prec_pt->use_lsc();

  // Set Navier Stokes mesh
  ns_prec_pt->set_navier_stokes_mesh(Bulk_mesh_pt);

  lgr_prec_pt->set_navier_stokes_lsc_preconditioner(ns_prec_pt);

  Prec_pt = lgr_prec_pt;




#ifdef OOMPH_HAS_TRILINOS

  // Build iterative linear solver
  oomph_info << "Using Trilinos GMRES\n"; 
  TrilinosAztecOOSolver* iterative_linear_solver_pt = 
    new TrilinosAztecOOSolver;

  Solver_pt=iterative_linear_solver_pt;

#else

  // Build solve and preconditioner
  Solver_pt = new GMRES<CRDoubleMatrix>;
  dynamic_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();

#endif

  // Set solver and preconditioner
  Solver_pt->preconditioner_pt() = Prec_pt;

  linear_solver_pt() = Solver_pt;

} // end_of_constructor

  template<class ELEMENT>
void CubeProblem <ELEMENT>::unsteady_run(DocInfo& doc_info)
{

  //Set value of dt
  double dt = Global_Variables::Delta_t;

  // Initialise all history values for an impulsive start
  assign_initial_values_impulsive(dt);
  cout << "IC = impulsive start" << std::endl;

  //Now do many timesteps
  unsigned ntsteps = Global_Variables::Time_end / dt;
  std::cout << "NTIMESTEP IS: " << ntsteps << std::endl; 


  // Doc initial condition
  doc_solution(doc_info);

  // increment counter
  doc_info.number()++;

  //Loop over the timesteps
  for(unsigned t=1;t<=ntsteps;t++)
  {
    cout << "TIMESTEP " << t << std::endl;

    //Take one fixed timestep
    unsteady_newton_solve(dt);

    //Output the time
    cout << "Time is now " << time_pt()->time() << std::endl;

    // Doc solution
    doc_solution(doc_info);

    // increment counter
    doc_info.number()++;
  }

} // end of unsteady run

//============start_of_fluid_traction_elements==============================
/// Create fluid traction elements 
//=======================================================================
template<class ELEMENT>
void CubeProblem<ELEMENT>::create_parall_outflow_lagrange_elements
(const unsigned &b, Mesh* const &bulk_mesh_pt, Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
 // Loop over the bulk elements adjacent to boundary b
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    Bulk_mesh_pt->boundary_element_pt(b,e));
     
   //What is the index of the face of the element e along boundary b
   int face_index = Bulk_mesh_pt->face_index_at_boundary(b,e);

   // Set the pointer to the prescribed traction function
   {
    // Build the corresponding impose_impenetrability_element
    ImposeParallelOutflowElement<ELEMENT>* flux_element_pt = new
     ImposeParallelOutflowElement<ELEMENT>(bulk_elem_pt,
                                           face_index);

//    flux_element_pt->set_tangent_direction(&Tangent_direction);
    surface_mesh_pt->add_element_pt(flux_element_pt);

    {
      const double minus_ang = -Global_Variables::Ang;
      // Loop over the nodes
      unsigned nnod=flux_element_pt->nnode();
      for (unsigned inod=0;inod<nnod;inod++)
      {
        Node* nod_pt = flux_element_pt->node_pt(inod);

        ///////// THIS IS FOR FULL FLOW
        //      if (  (nod_pt->is_on_boundary(1))||(nod_pt->is_on_boundary(5))
        //            ||(nod_pt->is_on_boundary(3))||(nod_pt->is_on_boundary(0)))
        //       {
        //        // How many nodal values were used by the "bulk" element
        //        // that originally created this node?
        //        unsigned n_bulk_value=flux_element_pt->nbulk_value(inod);
        //
        //        // The remaining ones are Lagrange multipliers and we pin them.
        //        unsigned nval=nod_pt->nvalue();
        //        for (unsigned j=n_bulk_value;j<nval;j++)
        //         {
        //          nod_pt->pin(j);
        //         }
        //       }

        // First, pin all the nodes on two boundaries
        std::set<unsigned>* bnd_pt=0;
        nod_pt->get_boundaries_pt(bnd_pt);
        if (bnd_pt!=0)
        {
          if (bnd_pt->size()>=2)
          {
            // How many nodal values were used by the "bulk" element
            // that originally created this node?
            unsigned n_bulk_value=flux_element_pt->nbulk_value(inod);

            // The remaining ones are Lagrange multipliers and we pin them.
            unsigned nval=nod_pt->nvalue();
            for (unsigned j=n_bulk_value;j<nval;j++)
            {
              nod_pt->pin(j);
            }

          }
        } // multiple boundaries

        // Now do the rest
        const double x = nod_pt->x(0);
        const double y = nod_pt->x(1);
        const double z = nod_pt->x(2);

        Vector<double>x_old;
        rotate_backward(x,y,z,minus_ang,minus_ang,minus_ang,x_old);

        if ((x_old[1] >= 0.5) || (x_old[2] >= 0.5))
        {
          // How many nodal values were used by the "bulk" element
          // that originally created this node?
          unsigned n_bulk_value=flux_element_pt->nbulk_value(inod);

          // The remaining ones are Lagrange multipliers and we pin them.
          unsigned nval=nod_pt->nvalue();
          for (unsigned j=n_bulk_value;j<nval;j++)
          {
            nod_pt->pin(j);
          }
        }

      }
    } // Encapsulation

   } // Encapsulation

  } // Loop over all elements on boundary b

} // end of create_parall_outflow_lagrange_elements


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


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void CubeProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 
 ofstream some_file;
 char filename[100];

 // Number of plot points
 unsigned npts;
 npts=5; 

 // Output solution 
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 Bulk_mesh_pt->output(some_file,npts);
 some_file.close();

} // end_of_doc_solution





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


  //Label for output
  DocInfo doc_info;

  //Set output directory
  doc_info.set_directory("RESLT");



  Global_Variables::Delta_t = atof(argv[1]);

  NavierStokesEquations<3>::Gamma[0]=1.0;
  NavierStokesEquations<3>::Gamma[1]=1.0;
  NavierStokesEquations<3>::Gamma[2]=1.0;

  unsigned nel_1d = atoi(argv[2]);

  // Build the problem 
  CubeProblem <QTaylorHoodElement<3> >problem(nel_1d);

  double re = 200.0;

  // Set Reynolds
  Global_Variables::Re=re;

  // Solve the problem 
//              problem.newton_solve();

  problem.distribute();
  problem.unsteady_run(doc_info);

  problem.doc_solution(doc_info);
  doc_info.number()++;


#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif


} // end_of_main


