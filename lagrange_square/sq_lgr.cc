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

// The 2D mesh
#include "meshes/simple_rectangular_quadmesh.h"
#include "meshes/rectangular_quadmesh.h"

// My own header
#include "./../rayheader.h"

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
   if(SquareLagrange::Using_trilinos_solver)
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

 void delete_flux_elements(Mesh* const &surface_mesh_pt);

 void actions_before_distribute()
 {
    delete_flux_elements(Surface_mesh_P_pt);
    rebuild_global_mesh();
 }

 void actions_after_distribute()
 {

   create_parall_outflow_lagrange_elements(1,Bulk_mesh_pt,Surface_mesh_P_pt);
   rebuild_global_mesh();
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

 void set_bc_for_SqTmp();
 void set_mesh_bc_for_SqPo();
 void set_bc_for_SqTf();
 void set_bc_for_SqTfPo();
 void set_bc_for_AwTmp();
 void set_bc_for_AwPo();
 void set_bc_for_AwTf();
 void set_bc_for_AwTfPo();

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
  namespace SL = SquareLagrange;

  Bottom_bound = 0;
  Right_bound = 1;
  Top_bound = 2;
  Left_bound = 3;
  
  Doc_linear_solver_info_pt = SL::Doc_linear_solver_info_pt;
  
  // First we set the mesh.
  if((SL::Prob_id >= 10) && (SL::Prob_id < 19) )
  {
    // This is the tilted cavity mesh.
    /// Setup the mesh
    // # of elements in x-direction
    unsigned nx=SL::Noel;
    
    // # of elements in y-direction
    unsigned ny=SL::Noel;
    
    // Domain length in x-direction
    double lx=SL::Lx;

    // Domain length in y-direction
    double ly=SL::Ly;
    
    Bulk_mesh_pt =
      new SlopingQuadMesh<ELEMENT>(nx,ny,lx,ly,SL::Ang);
      
  }
  else
  {
    pause("Not done yet."); 
    
  }

  set_mesh_bc_for_SqPo();



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
   el_pt->re_pt() = &SL::Rey;

  } // for(unsigned e=0;e<n_el;e++)

 //Assign equation numbers
 std::cout << "\n equation numbers : "<< assign_eqn_numbers() << std::endl;

 set_prec_and_solver();
 

 //////////////////////////////////////////////////////////////////////////////
 // Setting up the solver an preconditioners.


 //////////////////////////////////////////////////////////////////////////////
 // NS preconditioner
// ConstrainedNavierStokesSchurComplementPreconditioner* ns_preconditioner_pt =
// new ConstrainedNavierStokesSchurComplementPreconditioner;


}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::doc_solution()
{

  namespace SL = SquareLagrange;
  
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

//============RAYRAY===========
/// RAYRAY
//=======================================================================
template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::set_prec_and_solver()
{
  // Alias the namespace for convenience
  namespace SL = SquareLagrange;

  ////// Build the preconditioner
  LagrangeEnforcedflowPreconditioner* prec_pt
    = new LagrangeEnforcedflowPreconditioner;

  Prec_pt = prec_pt;

  Vector<Mesh*> mesh_pt;
  mesh_pt.resize(2);
  mesh_pt[0] = Bulk_mesh_pt;
  mesh_pt[1] = Surface_mesh_P_pt;
  //meshes_pt[2] = Surface_mesh_T_pt;
  
  prec_pt->set_meshes(mesh_pt);

  if(!SL::Use_axnorm)
  {
    prec_pt->scaling_sigma() = SL::Scaling_sigma;
  }
  
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
  unsigned po_b = 1; // parallel outflow

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
  unsigned num_bound = mesh_pt()->nboundary();

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
 unsigned current_bound = -1;
 
 // The number of nodes on a boundary.
 unsigned num_nod = -1;

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


template<class ELEMENT>
void TiltedCavityProblem<ELEMENT>::
delete_flux_elements(Mesh* const &surface_mesh_pt)
{
  // How many surface elements are there in the mesh?
  unsigned n_element = surface_mesh_pt->nelement();

  // Loop over the surface elements
  for(unsigned e=0;e<n_element;e++)
  {
    // Kill surface elements
    delete surface_mesh_pt->element_pt(e);
  }

  // Wipe the mesh
  surface_mesh_pt->flush_element_and_node_storage();
} //  end of delete_flux_elements


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

  // Alias the namespace for convenience.
  namespace NSPP = NavierStokesProblemParameters;
  namespace LPH = LagrangianPreconditionerHelpers;
  namespace SL = SquareLagrange;

  const unsigned dim = 2;

  LPH::PrecParam prec_param;

  // Set up doc info
  DocLinearSolverInfo doc_linear_solver_info;
  NSPP::Doc_linear_solver_info_pt = &doc_linear_solver_info;


  // Store commandline arguments
  CommandLineArgs::setup(argc,argv);

  NSPP::setup_commandline_flags();

  LPH::setup_commandline_flags(&prec_param);

  SL::setup_commandline_flags();

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();
  pause("Done new cl flags"); 
  

  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////

  // dim = 2
  NSPP::generic_problem_setup(dim);

  LPH::generic_setup(&prec_param);

  SL::generic_setup(&prec_param);

  pause("Paussso"); 
  
  // Solve with Taylor-Hood element, set up problem
  TiltedCavityProblem< QTaylorHoodElement<dim> > problem;


  //////////////////////////////////////////////////////////////////////////////

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

    for (SL::Rey = SL::Rey_start; 
        SL::Rey <= SL::Rey_end; SL::Rey += SL::Rey_incre)
    {
      std::ostringstream strs;
      strs << "R" << SL::Rey;
      //     SL::Rey_str = strs.str(); RAY RAY FIX THIS

      // Setup the label. Used for doc solution and preconditioner.
      SL::Label = SL::create_label(&prec_param);

      time_t rawtime;
      time(&rawtime);

      std::cout << "RAYDOING: "
        << SL::Label
        << " on " << ctime(&rawtime) << std::endl;


      // Solve the problem
      problem.newton_solve();

      //Output solution
      if(SL::Doc_soln){problem.doc_solution();}

      //////////////////////////////////////////////////////////////////////////
      ////////////// Outputting results ////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////

      // Get the global oomph-lib communicator 
      const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

      // My rank and number of processors. 
      // This is used later for putting the data.
      unsigned my_rank = comm_pt->my_rank();

      // Output the iteration counts and times if my_rank == 0
      if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir") 
          && (my_rank == 0))
      {

        // Create the File...
        std::ostringstream filename_stream;
        filename_stream << SL::Itstime_dir<<"/"<<SL::Label;
        std::ofstream outfile;
        outfile.open(filename_stream.str().c_str());

        // We now output the iteration and time.
        Vector<Vector<Vector<double> > > iters_times
          = SL::Doc_linear_solver_info_pt->iterations_and_times();

        // Below outputs the iteration counts and time.
        // Output the number of iterations
        // Since this is a steady state problem, there is only
        // one "time step".
        //*
        // Loop over the time steps and output the iterations, prec setup time and
        // linear solver time.
        //unsigned ntimestep = iters_times.size();
        //for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
        {
          // New timestep:
          outfile << "RAYITS:\t" << rey_increment << "\t";
          std::cout << "RAYITS:\t" << rey_increment << "\t";

          // Loop through the Newtom Steps
          unsigned nnewtonstep = iters_times[rey_increment].size();
          unsigned sum_of_newtonstep_iters = 0;
          for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
              innewtonstep++)
          {
            sum_of_newtonstep_iters += iters_times[rey_increment][innewtonstep][0];
            outfile << iters_times[rey_increment][innewtonstep][0] << " ";
            std::cout << iters_times[rey_increment][innewtonstep][0] << " ";
          }
          double average_its = ((double)sum_of_newtonstep_iters)
            / ((double)nnewtonstep);

          // Print to one decimal place if the average is not an exact
          // integer. Otherwise we print normally.
          std::streamsize cout_precision = std::cout.precision();
          ((unsigned(average_its*10))%10)?
            outfile << "\t"<< std::fixed << std::setprecision(1)
            << average_its << "(" << nnewtonstep << ")" << std::endl:
            outfile << "\t"<< average_its << "(" << nnewtonstep << ")" << std::endl;
          outfile << std::setprecision(cout_precision);

          ((unsigned(average_its*10))%10)?
            std::cout << "\t"<< std::fixed << std::setprecision(1)
            << average_its << "(" << nnewtonstep << ")" << std::endl:
            std::cout << "\t"<< average_its << "(" << nnewtonstep << ")" << std::endl;
          std::cout << std::setprecision(cout_precision);
        }

        // Now doing the preconditioner setup time.
        //for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
        {
          // New timestep:
          outfile << "RAYPRECSETUP:\t" << rey_increment << "\t";
          std::cout << "RAYPRECSETUP:\t" << rey_increment << "\t";
          // Loop through the Newton Steps
          unsigned nnewtonstep = iters_times[rey_increment].size();
          double sum_of_newtonstep_times = 0;
          for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
              innewtonstep++)
          {
            sum_of_newtonstep_times += iters_times[rey_increment][innewtonstep][1];
            outfile << iters_times[rey_increment][innewtonstep][1] << " ";
            std::cout << iters_times[rey_increment][innewtonstep][1] << " ";
          }
          double average_time = ((double)sum_of_newtonstep_times)
            / ((double)nnewtonstep);

          // Print to one decimal place if the average is not an exact
          // integer. Otherwise we print normally.
          outfile << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
          std::cout << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
        }

        // Now doing the linear solver time.
        //for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
        {
          // New timestep:
          outfile << "RAYLINSOLVER:\t" << rey_increment << "\t";
          std::cout << "RAYLINSOLVER:\t" << rey_increment << "\t";
          // Loop through the Newtom Steps
          unsigned nnewtonstep = iters_times[rey_increment].size();
          double sum_of_newtonstep_times = 0;
          for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
              innewtonstep++)
          {
            sum_of_newtonstep_times += iters_times[rey_increment][innewtonstep][2];
            outfile << iters_times[rey_increment][innewtonstep][2] << " ";
            std::cout << iters_times[rey_increment][innewtonstep][2] << " ";
          }
          double average_time = ((double)sum_of_newtonstep_times)
            / ((double)nnewtonstep);

          // Print to one decimal place if the average is not an exact
          // integer. Otherwise we print normally.
          outfile << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
          std::cout << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
        }
        outfile.close();
      } // if my_rank == 0

      rey_increment++;

    }
  }
  else
  {
    // Setup the label. Used for doc solution and preconditioner.
    SL::Label = SL::create_label(&prec_param);

    time_t rawtime;
    time(&rawtime);

    std::cout << "RAYDOING: "
      << SL::Label
      << " on " << ctime(&rawtime) << std::endl;

    problem.distribute();

    // Solve the problem
    problem.newton_solve();

    //Output solution
    if(SL::Doc_soln)
    {problem.doc_solution();}


    //////////////////////////////////////////////////////////////////////////
    ////////////// Outputting results ////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

    // Get the global oomph-lib communicator 
    const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

    // my rank and number of processors. 
    // This is used later for putting the data.
    unsigned my_rank = comm_pt->my_rank();
    unsigned nproc = comm_pt->nproc();

    // Variable to indicate if we want to output to a file or not.
    bool output_to_file = false;

    // The output file.
    std::ofstream outfile;

    // If we want to output to a file, we create the outfile.
    if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir"))
    {
      output_to_file = true;
      std::ostringstream filename_stream;
      filename_stream << SL::Itstime_dir<<"/"
        << SL::Label
        <<"NP"<<nproc<<"R"<<my_rank;
      outfile.open(filename_stream.str().c_str());
    }

    // Stringstream to hold the results. We do not output the results
    // (timing/iteration counts) as we get it since it will interlace with the
    // other processors and becomes hard to read.
    std::ostringstream results_stream;

    // Get the 3D vector which holds the iteration counts and timing results.
    Vector<Vector<Vector<double> > > iters_times
      = SL::Doc_linear_solver_info_pt->iterations_and_times();

    // Since this is a steady state problem, there is only
    // one "time step", thus it is essentially a 2D vector 
    // (the outer-most vector is of size 1).

    // Loop over the time steps and output the iterations, prec setup time and
    // linear solver time.
    unsigned ntimestep = iters_times.size();
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      // New timestep:
      results_stream << "RAYITS:\t" << intimestep << "\t";

      // Loop through the Newton Steps
      unsigned nnewtonstep = iters_times[intimestep].size();
      unsigned sum_of_newtonstep_iters = 0;
      for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
          innewtonstep++)
      {
        sum_of_newtonstep_iters += iters_times[intimestep][innewtonstep][0];
        results_stream << iters_times[intimestep][innewtonstep][0] << " ";
      }
      double average_its = ((double)sum_of_newtonstep_iters)
        / ((double)nnewtonstep);

      // Print to one decimal place if the average is not an exact
      // integer. Otherwise we print normally.
      ((unsigned(average_its*10))%10)?
        results_stream << "\t"<< std::fixed << std::setprecision(1)
        << average_its << "(" << nnewtonstep << ")" << "\n":
        results_stream << "\t"<< average_its << "(" << nnewtonstep << ")" << "\n";
      results_stream << std::setprecision(std::cout.precision());
    }

    // Now doing the preconditioner setup time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      // New timestep:
      results_stream << "RAYPRECSETUP:\t" << intimestep << "\t";
      // Loop through the Newtom Steps
      unsigned nnewtonstep = iters_times[intimestep].size();
      double sum_of_newtonstep_times = 0;
      for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
          innewtonstep++)
      {
        sum_of_newtonstep_times += iters_times[intimestep][innewtonstep][1];
        results_stream << iters_times[intimestep][innewtonstep][1] << " ";
      }
      double average_time = ((double)sum_of_newtonstep_times)
        / ((double)nnewtonstep);

      // Print to one decimal place if the average is not an exact
      // integer. Otherwise we print normally.
      results_stream << "\t"<< average_time << "(" << nnewtonstep << ")" << "\n";
    }

    // Now doing the linear solver time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      // New timestep:
      results_stream << "RAYLINSOLVER:\t" << intimestep << "\t";
      // Loop through the Newtom Steps
      unsigned nnewtonstep = iters_times[intimestep].size();
      double sum_of_newtonstep_times = 0;
      for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
          innewtonstep++)
      {
        sum_of_newtonstep_times += iters_times[intimestep][innewtonstep][2];
        results_stream << iters_times[intimestep][innewtonstep][2] << " ";
      }
      double average_time = ((double)sum_of_newtonstep_times)
        / ((double)nnewtonstep);

      // Print to one decimal place if the average is not an exact
      // integer. Otherwise we print normally.
      results_stream << "\t"<< average_time << "(" << nnewtonstep << ")" << "\n";
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
  } // else do not loop reynolds


#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main
