//LIC//====================================================================
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

#include <sstream>
#include <iomanip>
#include <ios>

// Oomph-lib includes
#include "generic.h"
#include "navier_stokes.h"

// Ray's own header!!
#include "./../rayheader.h"

// The 2D mesh
#include "meshes/simple_rectangular_quadmesh.h"
#include "meshes/rectangular_quadmesh.h"
#include "meshes/backward_step_mesh.h"

//using namespace std;
using namespace oomph;

namespace oomph
{
//========================================================================
/// \short A Sloping Mesh  class.
///
/// derived from RectangularQuadMesh:
/// the same mesh rotated with an angle phi
//========================================================================
 template<class ELEMENT>
 class SlopingQuadMesh : public BackwardStepQuadMesh<ELEMENT>
 {
 public:

  /// Constructor.
  SlopingQuadMesh(const unsigned& nx, const unsigned& ny,
                  const unsigned& nx_cut_out, const unsigned& ny_cut_out,
                  const double& lx,  const double& ly, 
                  const double& phi ) :
   RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly),
   BackwardStepQuadMesh<ELEMENT>(nx,ny,nx_cut_out,ny_cut_out,lx,ly)
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
class BackwardStepProblem : public Problem
{
public:

 /// \short Constructor: Pass number of elements in x and y directions and
 /// lengths
 BackwardStepProblem();

 /// Update before solve is empty
 void actions_before_newton_solve()
 {
   // Initialise counters for each newton solve.
   Doc_linear_solver_info_pt->setup_new_time_step();
 }

 /// \short Update after solve is empty
 void actions_after_newton_solve()
 {
 }

 void actions_after_newton_step()
 {
   unsigned iters = 0;
   double preconditioner_setup_time = 0.0;
   double solver_time = 0.0;

   // Get the iteration counts.
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
     preconditioner_setup_time 
       = iterative_solver_pt->preconditioner_pt()->setup_time();
   }
#else
   iters = static_cast<IterativeLinearSolver*>
             (this->linear_solver_pt())->iterations();
   preconditioner_setup_time = static_cast<IterativeLinearSolver*>
             (this->linear_solver_pt())->preconditioner_pt()->setup_time();

#endif

   // Get the preconditioner setup time.
   


   // Set the solver time.
   if(StepLagrange::Using_trilinos_solver)
   {
     TrilinosAztecOOSolver* trilinos_solver_pt 
       = dynamic_cast<TrilinosAztecOOSolver*>(this->linear_solver_pt());
     solver_time = trilinos_solver_pt->linear_solver_solution_time();
   }
   else
   {
     solver_time = linear_solver_pt()->linear_solver_solution_time();
   }

   Doc_linear_solver_info_pt->add_iteration_and_time
     (iters,preconditioner_setup_time,solver_time);
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

private:
 void set_prec_and_solver();

 /// Pointer to the "bulk" mesh
 SlopingQuadMesh<ELEMENT>* Bulk_mesh_pt;

 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_T_pt;
 Mesh* Surface_mesh_P_pt;

 // Preconditioner
 Preconditioner* Prec_pt;
 // Solver
 IterativeLinearSolver* Solver_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;
};



//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
template<class ELEMENT> // rrrback - changed here.
BackwardStepProblem<ELEMENT>::BackwardStepProblem()
{
 // Alias the namespace for convenience
 namespace SL = StepLagrange;
 
 Doc_linear_solver_info_pt = SL::Doc_linear_solver_info_pt;

 // Assign the boundaries:
 const unsigned if_b=3;
 //unsigned tf_b=1;
 const unsigned po_b=5;

 /// Setup the mesh
 
 // Domain length in x-direction
 const double lx=6.0;

 // Domain length in y-direction
 const double ly=2.0;

 // # of elements in x-direction
 const unsigned nx=SL::Noel * lx;

 // # of elements in y-direction
 const unsigned ny=SL::Noel * ly;

 // (number of elements to keep).
 const unsigned nx_cut_out = 5 * SL::Noel;
 const unsigned ny_cut_out = 1 * SL::Noel;

 Bulk_mesh_pt =
  new SlopingQuadMesh<ELEMENT>(nx,ny,nx_cut_out,ny_cut_out,lx,ly,SL::Ang);

 // Create a "surface mesh" that will contain only
 // ImposeParallelOutflowElements in boundary 1
 // The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Surface_mesh_P_pt = new Mesh;
 //Surface_mesh_T_pt = new Mesh;

 // Create ImposeParallelOutflowElement from all elements that are
 // adjacent to the Neumann boundary.
 create_parall_outflow_lagrange_elements(po_b,
                                         Bulk_mesh_pt,Surface_mesh_P_pt);
 //create_impenetrable_lagrange_elements(po_b,
 //                                        Bulk_mesh_pt,Surface_mesh_P_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_P_pt);
 //add_sub_mesh(Surface_mesh_T_pt);
 // Combine all submeshes into a single Mesh
 build_global_mesh();

 unsigned num_bound=Bulk_mesh_pt->nboundary();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 for(unsigned ibound=0;ibound<num_bound;ibound++)
 {
   //if((ibound != po_b)&&(ibound != tf_b))
   if(ibound != po_b)
   {
     unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
     {
       // Get node
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

       nod_pt->pin(0);
       nod_pt->pin(1);

     }
   }
 }

 unsigned num_nod= mesh_pt()->nboundary_node(if_b);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(if_b,inod);
   double x=nod_pt->x(0);
   double y=nod_pt->x(1);

   // Tilt it back
   double ytiltedback = x*sin(-SL::Ang)
                        +y*cos(-SL::Ang);
   double u=0.0;
   u=-4.0*(ytiltedback-1.0)*(2.0-ytiltedback);

   // Now apply the rotation to u, using rotation matrices.
   // with x = u and y = 0, i.e. R*[u;0] since we have the
   // velocity in the x direction only. There is no velocity
   // in the y direction.
   double ux=u*cos(SL::Ang);
   double uy=u*sin(SL::Ang);

   nod_pt->set_value(0,ux);
   nod_pt->set_value(1,uy);
 }

 //Complete the problem setup to make the elements fully functional

 //Loop over the elements
 unsigned n_el = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_el;e++)
  {
   //Cast to a fluid element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &SL::Rey;

  } // for(unsigned e=0;e<n_el;e++)

 //Assgn equation numbers
 std::cout << "\n equation numbers : "<< assign_eqn_numbers() << std::endl;

// set_prec_and_solver();

 ////// Build the preconditioner
// LagrangeEnforcedflowPreconditioner* prec_pt
//   = new LagrangeEnforcedflowPreconditioner;
 
// Prec_pt = prec_pt;

// Vector<Mesh*> mesh_pt;
// mesh_pt.resize(2);
// mesh_pt[0] = Bulk_mesh_pt;
// mesh_pt[1] = Surface_mesh_P_pt;
 //meshes_pt[2] = Surface_mesh_T_pt;
// prec_pt->set_meshes(mesh_pt);
 

// if(!SL::Use_axnorm)
// {
//   prec_pt->scaling_sigma() = SL::Scaling_sigma;
// }

 //////////////////////////////////////////////////////////////////////////////
 // Setting up the solver an preconditioners.

// // W solver. Use SuperLU
// if(SL::W_solver == 0)
// {
// }
// else
// {
//   std::cout << "Other W solvers not complemented yet. Using default SuperLU"
//             << std::endl;
// }

 //////////////////////////////////////////////////////////////////////////////
 // NS preconditioner
// ConstrainedNavierStokesSchurComplementPreconditioner* ns_preconditioner_pt =
// new ConstrainedNavierStokesSchurComplementPreconditioner;


// // The preconditioner for the fluid block:
// if(SL::NS_solver == 0) // Exact solve.
// {}
// else if(SL::NS_solver == 1) // LSC
// {
//   NavierStokesSchurComplementPreconditioner* ns_preconditioner_pt =
//     new NavierStokesSchurComplementPreconditioner(this);
//
//   prec_pt->set_navier_stokes_lsc_preconditioner(ns_preconditioner_pt);
//   ns_preconditioner_pt->set_navier_stokes_mesh(Bulk_mesh_pt);
//
//   // F block solve
//   // Preconditioner for the F block:
//   Preconditioner* f_preconditioner_pt = 0;
//   // SL::F_solver == 0 is default, so do nothing.
//   if(SL::F_solver == 11)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_2D_poison_problem();
//#endif
//   }
//   else if(SL::F_solver == 12)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_navier_stokes_momentum_block();
//#endif
//   }
//   else if(SL::F_solver == 13)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_CLJPGSStrn075();
//#endif
//   }
//   else if(SL::F_solver == 14)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_RSGSStrn075();
//#endif
//   }
//   else if(SL::F_solver == 15)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_CLJPPilutStrn075();
//#endif
//   }
//   else if(SL::F_solver == 16)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_RSPilutStrn075();
//#endif
//   }
//   else if(SL::F_solver == 17)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_augmented_momentum_block();
//#endif
//   }
//   else if(SL::F_solver == 81)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_CLJPGSStrn0668();
//#endif
//   }
//   else if(SL::F_solver == 82)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_CLJPJStrn0668();
//#endif
//   }
//   else if(SL::F_solver == 83)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_CLJPPilutStrn0668();
//#endif
//   }
//   else if(SL::F_solver == 84)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_RSGSStrn0668();
//#endif
//   }
//   else if(SL::F_solver == 85)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_RSJStrn0668();
//#endif
//   }
//   else if(SL::F_solver == 86)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // LSC takes type "Preconditioner".
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_for_RSPilutStrn0668();
//#endif
//   }
//   else if(SL::F_solver == 2)
//   {
////     f_preconditioner_pt = new RayBlockDiagonalPreconditioner<CRDoubleMatrix>;
//     f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
//   }
//   else if(SL::F_solver == 3)
//   {
//     f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
//#ifdef OOMPH_HAS_HYPRE
//     dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
//       (f_preconditioner_pt)->set_subsidiary_preconditioner_function
//       (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_for_2D_poison_problem);
//#endif
//   }
//   else if (SL::F_solver == 69)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     // AMG coarsening: Ruge-Stuben
//     RayGlobalAMGParam::amg_coarsening = 1;
//     
//     // AMG smoother: Gauss-Seidel
//     RayGlobalAMGParam::amg_smoother=0;
//     
//     // There is no damping with GS, otherwise we set the parameter:
//     // RayGlobalAMGParam::amg_damping
//
//     // Different amg strength for simple/stress divergence for viscuous term.
//     if(SL::Vis == 0)
//     {
//       // Simple form
//       RayGlobalAMGParam::amg_strength = 0.25;
//     }
//     else
//     {
//       // Stress divergence form
//       RayGlobalAMGParam::amg_strength = 0.668;
//     }
//     
//     // Setup the preconditioner.
//     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
//       set_hypre_using_2D_poisson_base();
//#endif
//   }
//
//   // Set the preconditioner in the LSC preconditioner.
//   ns_preconditioner_pt->set_f_preconditioner(f_preconditioner_pt);
//   
//   // P block solve
//   //SL::P_solver == 0 is default, so do nothing.
//   if(SL::P_solver == 1)
//   {
//#ifdef OOMPH_HAS_HYPRE
//     Preconditioner* p_preconditioner_pt = new HyprePreconditioner;
//
//     HyprePreconditioner* hypre_preconditioner_pt =
//       static_cast<HyprePreconditioner*>(p_preconditioner_pt);
//
//     Hypre_default_settings::
//     set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
//#endif
//   }
// } // if for using LSC as NS prec.
// else
// {
//   pause("There is no solver for NS.");
// }
//
// // Set the doc info for book keeping purposes.
// prec_pt->set_doc_linear_solver_info_pt(SL::Doc_linear_solver_info_pt);
//
// if(SL::Use_diagonal_w_block)
// {
//   prec_pt->use_diagonal_w_block();
// }
// else
// {
//   prec_pt->use_block_diagonal_w_block();
// }

// if(SL::Doc_prec)
// {
//   prec_pt->enable_doc_prec();
// }

 // Set the label, use to output information from the preconditioner, such
 // as the block matrices and the rhs vector
// prec_pt->set_label_pt(&SL::Label);

// // Build solve and preconditioner
//#ifdef OOMPH_HAS_TRILINOS
// TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
// trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
// Solver_pt = trilinos_solver_pt;
// SL::Using_trilinos_solver = true;
//#else
// Solver_pt = new GMRES<CRDoubleMatrix>;
// // We use RHS preconditioning. Note that by default,
// // left hand preconditioning is used.
// static_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();
// SL::Using_trilinos_solver = false;
//#endif
//
// Solver_pt->tolerance() = 1.0e-6;
// this->newton_solver_tolerance() = 1.0e-6;
//
// // Set solver and preconditioner
// Solver_pt->preconditioner_pt() = Prec_pt;
// linear_solver_pt() = Solver_pt;
}

//==start_of_doc_solution=================================================
/// Set up the solver and preconditioners.
//========================================================================
template<class ELEMENT>
void BackwardStepProblem<ELEMENT>::set_prec_and_solver()
{
  // Alias the namespace for convenience
  namespace SL = StepLagrange;

 ////// Build the preconditioner
 LagrangeEnforcedflowPreconditioner* prec_pt
   = new LagrangeEnforcedflowPreconditioner;

 Prec_pt = prec_pt;

 Vector<Mesh*> mesh_pt;
 mesh_pt.resize(2);
 mesh_pt[0] = Bulk_mesh_pt;
 mesh_pt[1] = Surface_mesh_P_pt;

 prec_pt->set_meshes(mesh_pt);

  if(!SL::Use_axnorm)
  {
    prec_pt->scaling_sigma() = SL::Scaling_sigma;
  }

   // W solver. Use SuperLU

   // W solver. Use SuperLU
 if(SL::W_solver == 0)
 {
 }
 else
 {
   std::cout << "Other W solvers not complemented yet. Using default SuperLU"
             << std::endl;
 }

 // The preconditioner for the fluid block:
 if(SL::NS_solver == 0) // Exact solve.
 {}
 else if(SL::NS_solver == 1) // LSC
 {
   NavierStokesSchurComplementPreconditioner* ns_preconditioner_pt =
     new NavierStokesSchurComplementPreconditioner(this);

   prec_pt->set_navier_stokes_lsc_preconditioner(ns_preconditioner_pt);
   ns_preconditioner_pt->set_navier_stokes_mesh(Bulk_mesh_pt);

   // F block solve
   // Preconditioner for the F block:
   Preconditioner* f_preconditioner_pt = 0;
   // SL::F_solver == 0 is default, so do nothing.
   if(SL::F_solver == 11)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_2D_poison_problem();
#endif
   }
   else if(SL::F_solver == 12)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_navier_stokes_momentum_block();
#endif
   }
   else if(SL::F_solver == 13)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPGSStrn075();
#endif
   }
   else if(SL::F_solver == 14)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSGSStrn075();
#endif
   }
   else if(SL::F_solver == 15)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPPilutStrn075();
#endif
   }
   else if(SL::F_solver == 16)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSPilutStrn075();
#endif
   }
   else if(SL::F_solver == 17)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_augmented_momentum_block();
#endif
   }
   else if(SL::F_solver == 81)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPGSStrn0668();
#endif
   }
   else if(SL::F_solver == 82)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPJStrn0668();
#endif
   }
   else if(SL::F_solver == 83)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPPilutStrn0668();
#endif
   }
   else if(SL::F_solver == 84)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSGSStrn0668();
#endif
   }
   else if(SL::F_solver == 85)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSJStrn0668();
#endif
   }
   else if(SL::F_solver == 86)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSPilutStrn0668();
#endif
   }
   else if(SL::F_solver == 2)
   {
//     f_preconditioner_pt = new RayBlockDiagonalPreconditioner<CRDoubleMatrix>;
     f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
   }
   else if(SL::F_solver == 3)
   {
     f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
#ifdef OOMPH_HAS_HYPRE
     dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
       (f_preconditioner_pt)->set_subsidiary_preconditioner_function
       (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_for_2D_poison_problem);
#endif
   }
   else if (SL::F_solver == 69)
   {
#ifdef OOMPH_HAS_HYPRE
     // AMG coarsening: Ruge-Stuben
     RayGlobalAMGParam::amg_coarsening = 1;
     
     // AMG smoother: Gauss-Seidel
     RayGlobalAMGParam::amg_smoother=0;
     
     // There is no damping with GS, otherwise we set the parameter:
     // RayGlobalAMGParam::amg_damping

     // Different amg strength for simple/stress divergence for viscuous term.
     if(SL::Vis == 0)
     {
       // Simple form
       RayGlobalAMGParam::amg_strength = 0.25;
     }
     else
     {
       // Stress divergence form
       RayGlobalAMGParam::amg_strength = 0.668;
     }
     
     // Setup the preconditioner.
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_using_2D_poisson_base();
#endif
   }
   else if (SL::F_solver == 96)
   {
#ifdef OOMPH_HAS_HYPRE
     // AMG coarsening:
     // Set: RayGlobalAMGParam::amg_coarsening = 
     // 0 - CLJP
     // 1 - RS
     
     // AMG smoother:
     // Set: RayGlobalAMGParam::amg_smoother = 
     // 0 - Jacobi (Need to set damping as well)
     // 1 - Gauss-Seidel
     // 2 - Pilut
     
     // There is no damping with GS, otherwise we set the parameter:
     // RayGlobalAMGParam::amg_damping

     RayGlobalAMGParam::amg_strength = SL::f_amg_strength;
     RayGlobalAMGParam::amg_damping = SL::f_amg_damping;
     RayGlobalAMGParam::amg_coarsening = SL::f_amg_coarsening;
     RayGlobalAMGParam::amg_smoother = SL::f_amg_smoother;
     RayGlobalAMGParam::amg_iterations = SL::f_amg_iterations;
     RayGlobalAMGParam::amg_smoother_iterations = SL::f_amg_smoother_iterations;
     RayGlobalAMGParam::print_hypre = SL::Print_hypre;


     // Different amg strength for simple/stress divergence for viscuous term.
     if(RayGlobalAMGParam::amg_strength < 0.0)
     {
       if(SL::Vis == 0)
       {
         // Simple form
         RayGlobalAMGParam::amg_strength = 0.25;
       }
       else
       {
         // Stress divergence form
         RayGlobalAMGParam::amg_strength = 0.668;
       }
     }
     
     // Setup the preconditioner.
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_ray();
#endif
   }

   // Set the preconditioner in the LSC preconditioner.
   ns_preconditioner_pt->set_f_preconditioner(f_preconditioner_pt);
   
   // P block solve
   //SL::P_solver == 0 is default, so do nothing.
   if(SL::P_solver == 1)
   {
#ifdef OOMPH_HAS_HYPRE
     Preconditioner* p_preconditioner_pt = new HyprePreconditioner;

     HyprePreconditioner* hypre_preconditioner_pt =
       static_cast<HyprePreconditioner*>(p_preconditioner_pt);

     Hypre_default_settings::
     set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif
   }
   else if(SL::P_solver == 96)
   {
#ifdef OOMPH_HAS_HYPRE
     Preconditioner* p_preconditioner_pt = 0;
     
//* 
     RayGlobalAMGParam::amg_iterations = SL::p_amg_iterations;
     RayGlobalAMGParam::amg_smoother_iterations = SL::p_amg_smoother_iterations;
     RayGlobalAMGParam::amg_smoother = SL::p_amg_smoother;
     RayGlobalAMGParam::amg_strength = SL::p_amg_strength;
     //RayGlobalAMGParam::amg_damping = SL::p_amg_damping;
     RayGlobalAMGParam::amg_coarsening = SL::p_amg_coarsening;
     RayGlobalAMGParam::print_hypre = SL::Print_hypre;

//     std::cout << "p_amg_iterations:" << SL::p_amg_iterations << std::endl; 
//     std::cout << "p_amg_smoother_iterations" << SL::p_amg_smoother_iterations << std::endl; 
//     std::cout << "p_amg_strength" << SL::p_amg_strength << std::endl;
//     std::cout << "p_amg_coarsening" << SL::p_amg_coarsening << std::endl; 
// */

     p_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::set_hypre_ray();

     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif
   }
   else if(SL::P_solver == 2)
   {
#ifdef OOMPH_HAS_HYPRE
     Preconditioner* p_preconditioner_pt = new HyprePreconditioner;

     HyprePreconditioner* hypre_preconditioner_pt =
       static_cast<HyprePreconditioner*>(p_preconditioner_pt);

     hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

     // Setup v-cycles
     hypre_preconditioner_pt->set_amg_iterations(2);
     hypre_preconditioner_pt->amg_smoother_iterations() = 2;
     
     // Setup smoother
     // simple: 0 - DJ, 1 - GS
     // compelx: Pilut - 7
     hypre_preconditioner_pt->amg_using_simple_smoothing();
     hypre_preconditioner_pt->amg_simple_smoother() = 0;
     // only applicable for DJ
     hypre_preconditioner_pt->amg_damping() = 0.8;

     // Setup coarsening
     // 0 - CLJP
     // 1 - RS
     hypre_preconditioner_pt->amg_coarsening() = 1;

     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif

   }
 } // if for using LSC as NS prec.
 else
 {
   pause("There is no solver for NS.");
 }


 // Set the doc info for book keeping purposes.
 prec_pt->set_doc_linear_solver_info_pt(SL::Doc_linear_solver_info_pt);

 if(SL::Use_block_diagonal_w)
 {
   prec_pt->use_block_diagonal_w_block();
 }
 else
 {
   prec_pt->use_diagonal_w_block();
 }

 if(SL::Doc_prec)
 {
   prec_pt->enable_doc_prec();
 }

 // Set the label, use to output information from the preconditioner, such
 // as the block matrices and the rhs vector
 prec_pt->set_label_pt(&SL::Label);
 prec_pt->set_doc_prec_directory_pt(&SL::Doc_prec_dir);

 // Build solve and preconditioner
#ifdef OOMPH_HAS_TRILINOS
 TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
 trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
 Solver_pt = trilinos_solver_pt;
 SL::Using_trilinos_solver = true;
#else
 Solver_pt = new GMRES<CRDoubleMatrix>;
 // We use RHS preconditioning. Note that by default,
 // left hand preconditioning is used.
 static_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();
 SL::Using_trilinos_solver = false;
#endif

 Solver_pt->tolerance() = 1.0e-6;
 this->newton_solver_tolerance() = 1.0e-6;

 // Set solver and preconditioner
 Solver_pt->preconditioner_pt() = Prec_pt;
 linear_solver_pt() = Solver_pt;




}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void BackwardStepProblem<ELEMENT>::doc_solution()
{

  namespace SL = StepLagrange;
  
  std::ofstream some_file;
  std::stringstream filename;
  filename << SL::Soln_dir<<"/"<<SL::Label<<".dat";

  // Number of plot points
  unsigned npts=5;

  // Output solution
  some_file.open(filename.str().c_str());
  Bulk_mesh_pt->output(some_file,npts);
  some_file.close();
}



//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void BackwardStepProblem<ELEMENT>::
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
                                          face_index);


   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 0 or 4?
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(4)))
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
void BackwardStepProblem<ELEMENT>::
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

//===start_of_main======================================================
/// Driver code
//======================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);
#endif
  // Get the global oomph-lib communicator 
  //  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

  // my rank and number of processors. This is used later for putting the data.
  //  unsigned my_rank = comm_pt->my_rank();
  //unsigned nproc = comm_pt->nproc();


  // Alias the namespace for convenience.
  namespace SL = StepLagrange;

  // Set up doc info
  DocLinearSolverInfo doc_linear_solver_info;

  SL::Doc_linear_solver_info_pt = &doc_linear_solver_info;

  //SL::Soln_dir = "RESLT";

  // Store commandline arguments
  CommandLineArgs::setup(argc,argv);

  CommandLineArgs::specify_command_line_flag("--dist_prob");

  // Flag to output the solution.
  CommandLineArgs::specify_command_line_flag("--doc_soln", &SL::Soln_dir);
  // Flag to output the preconditioner, used for debugging.
  CommandLineArgs::specify_command_line_flag("--doc_prec", &SL::Doc_prec_dir);

  // A problem ID, there are eight different types of problems.
  // Check the header file.
  CommandLineArgs::specify_command_line_flag("--prob_id",&SL::Prob_id);

  CommandLineArgs::specify_command_line_flag("--w_solver", &SL::W_solver);
  CommandLineArgs::specify_command_line_flag("--ns_solver", &SL::NS_solver);
  CommandLineArgs::specify_command_line_flag("--p_solver", &SL::P_solver);
  CommandLineArgs::specify_command_line_flag("--f_solver", &SL::F_solver);
  CommandLineArgs::specify_command_line_flag("--visc", &SL::Vis);
  CommandLineArgs::specify_command_line_flag("--ang", &SL::Ang_deg);
  CommandLineArgs::specify_command_line_flag("--rey", &SL::Rey);
  CommandLineArgs::specify_command_line_flag("--rey_start", &SL::Rey_start);
  CommandLineArgs::specify_command_line_flag("--rey_incre", &SL::Rey_incre);
  CommandLineArgs::specify_command_line_flag("--rey_end", &SL::Rey_end);
  CommandLineArgs::specify_command_line_flag("--noel", &SL::Noel);
  CommandLineArgs::specify_command_line_flag("--sigma",
      &SL::Scaling_sigma);
  CommandLineArgs::specify_command_line_flag("--bdw");

  // Iteration count and times directory.
  CommandLineArgs::specify_command_line_flag("--itstimedir", &SL::Itstime_dir);

  // NS_F block AMG parameters
  CommandLineArgs::specify_command_line_flag("--f_amg_str", 
      &SL::f_amg_strength);
  CommandLineArgs::specify_command_line_flag("--f_amg_damp", 
      &SL::f_amg_damping);
  CommandLineArgs::specify_command_line_flag("--f_amg_coarse", 
      &SL::f_amg_coarsening);
  CommandLineArgs::specify_command_line_flag("--f_amg_smoo", 
      &SL::f_amg_smoother);
  CommandLineArgs::specify_command_line_flag("--f_amg_iter", 
      &SL::f_amg_iterations);
  CommandLineArgs::specify_command_line_flag("--f_amg_smiter", 
      &SL::f_amg_smoother_iterations);

  // NS_P block AMG parameters
  CommandLineArgs::specify_command_line_flag("--p_amg_str", 
      &SL::p_amg_strength);
  CommandLineArgs::specify_command_line_flag("--p_amg_damp", 
      &SL::p_amg_damping);
  CommandLineArgs::specify_command_line_flag("--p_amg_coarse", 
      &SL::p_amg_coarsening);
  CommandLineArgs::specify_command_line_flag("--p_amg_smoo", 
      &SL::p_amg_smoother);
  CommandLineArgs::specify_command_line_flag("--p_amg_iter", 
      &SL::p_amg_iterations);
  CommandLineArgs::specify_command_line_flag("--p_amg_smiter", 
      &SL::p_amg_smoother_iterations);


  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////

  // Do we have to distribute the problem?
  if(CommandLineArgs::command_line_flag_has_been_set("--dist_prob"))
  {
    SL::Distribute_problem = true;
  }
  else
  {
    SL::Distribute_problem = false;
  }

  // Document the solution? Default is false.
  if(CommandLineArgs::command_line_flag_has_been_set("--doc_soln"))
  {
    // The argument immediately after --doc_soln is put into SL::Soln_dir.
    // If this begins with "--", then no solution directory has been provided.
    std::size_t found = SL::Soln_dir.find("--");

    // Check if they have set the solution directory.
    if(found != std::string::npos)
    {
      std::ostringstream err_msg;
      err_msg << "Please provide the doc_soln directory "
        << "after the argument --doc_soln.\n" 
        << "This must not start with \"--\"." << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      SL::Doc_soln = true;
    }
  }


  // Document the preconditioner? Default is false.
  if(CommandLineArgs::command_line_flag_has_been_set("--doc_prec"))
  {
    // The argument immediately after --doc_prec is put into SL::Doc_prec_dir.
    // If this begins with "--", then no prec directory has been provided.
    std::size_t found = SL::Doc_prec_dir.find("--");

    // Check if they have set the doc_prec directory.
    if(found != std::string::npos)
    {
      std::ostringstream err_msg;
      err_msg << "Please provide the doc_prec directory "
        << "after the argument --doc_prec.\n" 
        << "This must not start with \"--\"." << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      SL::Doc_prec = true;
    }
  }


    // Set a problem id to identify the problem.
  // This is used for book keeping purposes.
  if(CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
  {
    // The argument immediately after --prob_id is put into SL::Prob_id.
    // If this begins with "--", then no problem id has been provided.

    // Maybe I should check if SL::Prob_id is a number or a string...

    // We only accept problem IDs as defined below.
    // Creating a set of acceptable IDs
    int prob_id_array[]= {10,11};

    const bool inset = check_if_in_set<int>(prob_id_array,2,SL::Prob_id);

    // Check if they have provided an acceptable ID.
    // If a new element has been inserted, it means the user has provided an
    // ID not in the set.
    if(inset == false)
    {
      std::ostringstream err_msg;
      err_msg << "Please provide a problem id to identify the problem after "
        << "after the argument --prob_id.\n" 
        << "Acceptable IDs are:\n"
        << "10 = (StpTmp) Step, custom stuff...\n"
        << "11 = (StpPo) Step, Parallel outflow (para inflow)\n"
        << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  }


  // Set the viscuous term.
  // Default: 0, Sim
  if(CommandLineArgs::command_line_flag_has_been_set("--visc"))
  {
    if (SL::Vis == 0)
    {
      NavierStokesEquations<2>::Gamma[0]=0.0;
      NavierStokesEquations<2>::Gamma[1]=0.0;

    }
    else if (SL::Vis == 1)
    {
      NavierStokesEquations<2>::Gamma[0]=1.0;
      NavierStokesEquations<2>::Gamma[1]=1.0;
    } // else - setting viscuous term.
    else
    {
      std::ostringstream err_msg;
      err_msg << "Do not recognise viscuous term: " << SL::Vis << ".\n"
        << "Vis = 0 for simple form\n"
        << "Vis = 1 for stress divergence form\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  }

  // Now we need to convert Ang into radians.
  SL::Ang = SL::Ang_deg * (MathematicalConstants::Pi / 180.0);



  // Now we need to convert Ang into radians.
  SL::Ang = SL::Ang_deg * (MathematicalConstants::Pi / 180.0);

  // Check if the Reynolds numbers have been set.
  if(  CommandLineArgs::command_line_flag_has_been_set("--rey_start")
      &&CommandLineArgs::command_line_flag_has_been_set("--rey_incre")
      &&CommandLineArgs::command_line_flag_has_been_set("--rey_end")
      &&CommandLineArgs::command_line_flag_has_been_set("--rey"))
  {
    std::ostringstream err_msg;
    err_msg << "You have set all --rey* argument, please choose carefully!\n"
      << std::endl;
    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }
  else if(  CommandLineArgs::command_line_flag_has_been_set("--rey_start")
      &&CommandLineArgs::command_line_flag_has_been_set("--rey_incre")
      &&CommandLineArgs::command_line_flag_has_been_set("--rey_end"))
  {
    std::cout << "Looping Reynolds: \n"
      << "Rey_start = " << SL::Rey_start << std::endl; 
    std::cout << "Rey_incre = " << SL::Rey_incre << std::endl; 
    std::cout << "Rey_end = " << SL::Rey_end << std::endl; 
  }
  else if(!CommandLineArgs::command_line_flag_has_been_set("--rey"))
  {
    std::ostringstream err_msg;
    err_msg << "No Reynolds numbers have been set.\n"
      << "For a single Reynolds number, use --rey.\n"
      << "For looping through Reynolds numbers, use:\n"
      << "--rey_start --rey_incre --rey_end.\n"
      << std::endl;
    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }







  // Set a string to identify the problem. This is unique to each problem,
  // so we hard code this. 2DStrPo = 2 dimension, straight parallel outflow.
  // straight describes the velocity flow field. Po = Parallel outflow
  // describes the boundary type.
  // SL::Prob_str = "StepPo";

  // Set the strings to identify the preconditioning,
  // This is used purely for book keeping purposes.

  // // Default: W_solver = 0, W_str = We
  // if(CommandLineArgs::command_line_flag_has_been_set("--w_solver"))
  // {
  //  switch(SL::W_solver)
  //  {
  //    case 0:
  //      SL::W_str = "We";
  //      break;
  //    case 1:
  //    {
  //      pause("No other W block solver coded."); 
  //      SL::W_str = "Wa";
  //    }
  //      break;
  //    default:
  //      std::cout << "Do not recognise W: " << SL::W_solver << "\n"
  //                << "Exact preconditioning = 0\n"
  //                << "AMG = 1\n"
  //                << "Using default: Exact (W_solver = 0)"<< std::endl;
  //  }  // switch
  // } // if

  // // Default: NS_solver = 1, NS_str = Nl
  // if(CommandLineArgs::command_line_flag_has_been_set("--ns_solver"))
  // {
  //  switch(SL::NS_solver)
  //  {
  //    case 0:
  //      SL::NS_str = "Ne";
  //      SL::P_str = "";
  //      SL::F_str = "";
  //      break;
  //    case 1:
  //      SL::NS_str = "Nl";
  //      break;
  //    default:
  //      std::cout << "Do not recognise NS: " << SL::NS_solver << "\n"
  //                << "Exact solve = 0\n"
  //                << "LSC = 1\n"
  //                << "Using default: LSC for NS block (NS_solver = 1)"<<std::endl;
  //  }  // switch
  // } // if

  // // Default: This can only be set if NS_solver != 0 i.e. we are using LSC for 
  // // the NS block
  // // Default: P_solver = 0, P_str = Pe
  // if(CommandLineArgs::command_line_flag_has_been_set("--p_solver"))
  // {
  //  if(SL::NS_solver == 0)
  //  {
  //    pause("NS solve is exact. There cannot be a P solver.");
  //  }
  //
  //  switch(SL::P_solver)
  //  {
  //    case 0:
  //      SL::P_str = "Pe";
  //      break;
  //    case 1:
  //      SL::P_str = "Pa";
  //      break;
  //    default:
  //      std::cout << "Do not recognise P: " << SL::P_solver << "\n"
  //                << "Exact preconditioning = 0\n"
  //                << "AMG = 1\n"
  //                << "Using default: Exact P solve (P_solver = 0)"<< std::endl;
  //  }  // switch
  // } // if

  // // This can only be set if we're using LSC for the NS block.
  // // Default: 0, Fe
  // if(CommandLineArgs::command_line_flag_has_been_set("--f_solver"))
  // {
  //  if(SL::NS_solver == 0)
  //  {
  //    pause("NS solve is exact. There cannot be an F solver.");
  //  }
  //
  //  switch(SL::F_solver)
  //  {
  //    case 0:
  //      SL::F_str = "Fe";
  //      break;
  //    case 69:
  //      SL::F_str = "Fa";
  //      break;
  //    case 11:
  //      SL::F_str = "Fh2dp";
  //      break;
  //    case 12:
  //      SL::F_str = "Fhns";
  //      break;
  //    case 13:
  //      SL::F_str = "CLJPGSStrn075";
  //      break;
  //    case 14:
  //      SL::F_str = "FRSGSStrn075";
  //      break;
  //    case 15:
  //      SL::F_str = "FCLJPPilutStrn075";
  //      break;
  //    case 16:
  //      SL::F_str = "FRSPilutStrn075";
  //      break;
  //    case 17:
  //      SL::F_str = "Fray"; // I have no short hand for this...
  //      break;
  //    case 81:
  //      SL::F_str = "CLJPGSStrn0668";
  //      break;
  //    case 82:
  //      SL::F_str = "CLJPJStrn0668";
  //      break;
  //    case 83:
  //      SL::F_str = "CLJPPilutStrn0668";
  //      break;
  //    case 84:
  //      SL::F_str = "RSGSStrn0668";
  //      break;
  //    case 85:
  //      SL::F_str = "RSJStrn0668";
  //      break;
  //    case 86:
  //      SL::F_str = "RSPilutStrn0668";
  //      break;
  //    case 2:
  //      SL::F_str = "Fde";
  //      break;
  //    case 3:
  //      SL::F_str = "Fda";
  //      break;
  //    default:
  //      std::cout << "Do not recognise F: " << SL::F_solver << "\n"
  //                << "Exact preconditioning = 0\n"
  //                << "AMG = xxx Look in the code...\n"
  //                << "Using default: Exact F solve (F_solver = 0)"<< std::endl;
  //  }  // switch
  // } // if

  // // Set the viscuous term.
  // // Default: 0, Sim
  // if(CommandLineArgs::command_line_flag_has_been_set("--visc"))
  // {
  //   if (SL::Vis == 0)
  //   {
  //     SL::Vis_str = "Sim";
  //     NavierStokesEquations<2>::Gamma[0]=0.0;
  //     NavierStokesEquations<2>::Gamma[1]=0.0;
  //
  //   }
  //   else if (SL::Vis == 1)
  //   {
  //     SL::Vis_str = "Str";
  //     NavierStokesEquations<2>::Gamma[0]=1.0;
  //     NavierStokesEquations<2>::Gamma[1]=1.0;
  //   } // else - setting viscuous term.
  //   else
  //   {
  //     std::cout << "There is no such Viscous term, using 0 = simple." 
  //               << std::endl; 
  //   }
  // }

  // // Set Ang_str
  // // Default: A30
  // if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
  // {
  //   std::ostringstream strs;
  //   strs << "A" << SL::Ang;
  //   SL::Ang_str = strs.str();
  // }

  // Now we need to convert Ang into radians.
  //SL::Ang = SL::Ang * (MathematicalConstants::Pi / 180.0);

  // // Set Noel_str, used for book keeping.
  // if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
  // {
  //   std::ostringstream strs;
  //   strs << "N" << SL::Noel;
  //   SL::Noel_str = strs.str();
  // }

  // // Set Use_axnorm, if sigma has not been set, norm os momentum block is used.
  // if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
  // {
  //   SL::Use_axnorm = false;
  //
  //   std::ostringstream strs;
  //   strs << "S" << SL::Scaling_sigma;
  //   SL::Sigma_str = strs.str();
  // }

  // // use the diagonal or block diagonal approximation for W block.
  // if(CommandLineArgs::command_line_flag_has_been_set("--bdw"))
  // {
  //   SL::Use_diagonal_w_block = false;
  //   SL::W_approx_str = "bdw";
  // }
  // else
  // {
  //   SL::Use_diagonal_w_block = true;
  //   SL::W_approx_str = "";
  // }

  // // Solve with Taylor-Hood element, set up problem
  // // Set Rey_str, used for book keeping.
  // if(CommandLineArgs::command_line_flag_has_been_set("--rey"))
  // {
  //   if(SL::Rey < 0)
  //   {
  //     SL::Loop_reynolds = true;
  //   }
  //   else
  //   {
  //     std::ostringstream strs;
  //     strs << "R" << SL::Rey;
  //     SL::Rey_str = strs.str();
  //   }
  // }

  BackwardStepProblem< QTaylorHoodElement<2> > problem;


  ///////////////////////////////////////////////////////////////////////////////
















  // if(SL::Loop_reynolds)
  // {
  //   double Rey_start = 0.0;
  //   double Rey_end = 200.0;
  //   unsigned rey_increment = 0; // used for output of iters/times
  //   for (SL::Rey = Rey_start; SL::Rey <= Rey_end; SL::Rey += 50.0) 
  //   {
  //     std::ostringstream strs;
  //     strs << "R" << SL::Rey;
  //     SL::Rey_str = strs.str();
  //
  //     // Setup the label. Used for doc solution and preconditioner.
  //     SL::Label = SL::Prob_str
  //       + SL::W_str + SL::NS_str + SL::F_str + SL::P_str
  //       + SL::Vis_str + SL::Ang_str + SL::Rey_str
  //       + SL::Noel_str + SL::W_approx_str + SL::Sigma_str;
  //
  //     time_t rawtime;
  //     time(&rawtime);
  //
  //     std::cout << "RAYDOING: "
  //       << SL::Label
  //       << " on " << ctime(&rawtime) << std::endl;
  //
  //
  //     // Solve the problem
  //     problem.newton_solve();
  //     
  //     //Output solution
  //     if(SL::Doc_soln){problem.doc_solution();}
  //     
  //     if(my_rank == 0)
  //     {
  //     // Output the iteration counts and times if my_rank = 0
  //     // Create the File...
  //     std::ostringstream filename_stream;
  //     filename_stream << "runs"<<SL::Prob_str<<"/RAYOUT"<<SL::Label;
  //     
  //     std::ofstream outfile;
  //     outfile.open(filename_stream.str().c_str());
  //
  //     // We now output the iteration and time.
  //     Vector<Vector<Vector<double> > > iters_times
  //       = SL::Doc_linear_solver_info_pt->iterations_and_times();
  //
  //     // Below outputs the iteration counts and time.
  //     // Output the number of iterations
  //     // Since this is a steady state problem, there is only
  //     // one "time step".
  //     //*
  //     // Loop over the time steps and output the iterations, prec setup time and
  //     // linear solver time.
  //     //unsigned ntimestep = iters_times.size();
  //     //for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
  //     {
  //       // New timestep:
  //       outfile << "RAYITS:\t" << rey_increment << "\t";
  //     
  //       // Loop through the Newtom Steps
  //       unsigned nnewtonstep = iters_times[rey_increment].size();
  //       unsigned sum_of_newtonstep_iters = 0;
  //       for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
  //           innewtonstep++)
  //       {
  //         sum_of_newtonstep_iters += iters_times[rey_increment][innewtonstep][0];
  //         outfile << iters_times[rey_increment][innewtonstep][0] << " ";
  //       }
  //       double average_its = ((double)sum_of_newtonstep_iters)
  //         / ((double)nnewtonstep);
  //  
  //       // Print to one decimal place if the average is not an exact
  //       // integer. Otherwise we print normally.
  //       std::streamsize cout_precision = outfile.precision();
  //       ((unsigned(average_its*10))%10)?
  //         outfile << "\t"<< std::fixed << std::setprecision(1)
  //         << average_its << "(" << nnewtonstep << ")" << std::endl:
  //         outfile << "\t"<< average_its << "(" << nnewtonstep << ")" << std::endl;
  //       outfile << std::setprecision(cout_precision);
  //     }
  //
  //     // Now doing the preconditioner setup time.
  //     //for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
  //     {
  //       // New timestep:
  //       outfile << "RAYPRECSETUP:\t" << rey_increment << "\t";
  //       // Loop through the Newtom Steps
  //       unsigned nnewtonstep = iters_times[rey_increment].size();
  //       double sum_of_newtonstep_times = 0;
  //       for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
  //           innewtonstep++)
  //       {
  //         sum_of_newtonstep_times += iters_times[rey_increment][innewtonstep][1];
  //         outfile << iters_times[rey_increment][innewtonstep][1] << " ";
  //       }
  //       double average_time = ((double)sum_of_newtonstep_times)
  //         / ((double)nnewtonstep);
  //  
  //       // Print to one decimal place if the average is not an exact
  //       // integer. Otherwise we print normally.
  //       outfile << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
  //     }
  //
  //     // Now doing the linear solver time.
  //     //for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
  //     {
  //       // New timestep:
  //       outfile << "RAYLINSOLVER:\t" << rey_increment << "\t";
  //       // Loop through the Newtom Steps
  //       unsigned nnewtonstep = iters_times[rey_increment].size();
  //       double sum_of_newtonstep_times = 0;
  //       for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
  //           innewtonstep++)
  //       {
  //         sum_of_newtonstep_times += iters_times[rey_increment][innewtonstep][2];
  //         outfile << iters_times[rey_increment][innewtonstep][2] << " ";
  //       }
  //       double average_time = ((double)sum_of_newtonstep_times)
  //         / ((double)nnewtonstep);
  //
  //       // Print to one decimal place if the average is not an exact
  //       // integer. Otherwise we print normally.
  //       outfile << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
  //     }
  //     } // if rank == 0 output stuff
  //     rey_increment++;
  //   }
  // }
  // else
  // {
  //   // Setup the label. Used for doc solution and preconditioner.
  //   SL::Label = SL::Prob_str
  //               + SL::W_str + SL::NS_str + SL::F_str + SL::P_str
  //               + SL::Vis_str + SL::Ang_str + SL::Rey_str
  //               + SL::Noel_str + SL::W_approx_str + SL::Sigma_str;
  //
  //   time_t rawtime;
  //   time(&rawtime);
  //
  //   std::cout << "RAYDOING: "
  //     << SL::Label
  //     << " on " << ctime(&rawtime) << std::endl;
  //
  //   // Solve the problem
  //   problem.newton_solve();
  //
  //   //Output solution
  //   if(SL::Doc_soln){problem.doc_solution();}
  //
  //   // We now output the iteration and time.
  //   Vector<Vector<Vector<double> > > iters_times
  //     = SL::Doc_linear_solver_info_pt->iterations_and_times();
  //
  //   // Below outputs the iteration counts and time.
  //   // Output the number of iterations
  //   // Since this is a steady state problem, there is only
  //   // one "time step".
  //   //*
  //   // Loop over the time steps and output the iterations, prec setup time and
  //   // linear solver time.
  //   unsigned ntimestep = iters_times.size();
  //   for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
  //   {
  //     // New timestep:
  //     std::cout << "RAYITS:\t" << intimestep << "\t";
  //     
  //     // Loop through the Newtom Steps
  //     unsigned nnewtonstep = iters_times[intimestep].size();
  //     unsigned sum_of_newtonstep_iters = 0;
  //     for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
  //         innewtonstep++)
  //     {
  //       sum_of_newtonstep_iters += iters_times[intimestep][innewtonstep][0];
  //       std::cout << iters_times[intimestep][innewtonstep][0] << " ";
  //     }
  //     double average_its = ((double)sum_of_newtonstep_iters)
  //       / ((double)nnewtonstep);
  //
  //     // Print to one decimal place if the average is not an exact
  //     // integer. Otherwise we print normally.
  //     std::streamsize cout_precision = std::cout.precision();
  //     ((unsigned(average_its*10))%10)?
  //       std::cout << "\t"<< std::fixed << std::setprecision(1)
  //       << average_its << "(" << nnewtonstep << ")" << std::endl:
  //       std::cout << "\t"<< average_its << "(" << nnewtonstep << ")" << std::endl;
  //     std::cout << std::setprecision(cout_precision);
  //   }
  //
  //   // Now doing the preconditioner setup time.
  //   for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
  //   {
  //     // New timestep:
  //     std::cout << "RAYPRECSETUP:\t" << intimestep << "\t";
  //     // Loop through the Newtom Steps
  //     unsigned nnewtonstep = iters_times[intimestep].size();
  //     double sum_of_newtonstep_times = 0;
  //     for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
  //         innewtonstep++)
  //     {
  //       sum_of_newtonstep_times += iters_times[intimestep][innewtonstep][1];
  //       std::cout << iters_times[intimestep][innewtonstep][1] << " ";
  //     }
  //     double average_time = ((double)sum_of_newtonstep_times)
  //       / ((double)nnewtonstep);
  //
  //     // Print to one decimal place if the average is not an exact
  //     // integer. Otherwise we print normally.
  //     std::cout << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
  //   }
  //
  //   // Now doing the linear solver time.
  //   for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
  //   {
  //     // New timestep:
  //     std::cout << "RAYLINSOLVER:\t" << intimestep << "\t";
  //     // Loop through the Newtom Steps
  //     unsigned nnewtonstep = iters_times[intimestep].size();
  //     double sum_of_newtonstep_times = 0;
  //     for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
  //         innewtonstep++)
  //     {
  //       sum_of_newtonstep_times += iters_times[intimestep][innewtonstep][2];
  //       std::cout << iters_times[intimestep][innewtonstep][2] << " ";
  //     }
  //     double average_time = ((double)sum_of_newtonstep_times)
  //       / ((double)nnewtonstep);
  //
  //     // Print to one decimal place if the average is not an exact
  //     // integer. Otherwise we print normally.
  //     std::cout << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
  //   }
  //
  // }


#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main