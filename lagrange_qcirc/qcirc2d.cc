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
// Driver for adaptive 2D quarter circle driven cavity. Solved with black
// box adaptation, using Taylor Hood and Crouzeix Raviart elements.

// Generic oomph-lib header
#include "generic.h"

// Navier Stokes headers
#include "navier_stokes.h"

// The mesh
#include "meshes/quarter_circle_sector_mesh.h"

// My own header
#include "./../rayheader.h"

using namespace std;

using namespace oomph;


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

//==start_of_problem_class============================================
/// Driven cavity problem in quarter circle domain, templated
/// by element type. 
//====================================================================
template<class ELEMENT>
class QuarterCircleProblem : public Problem
{

public:

 /// Constructor
 QuarterCircleProblem(
  NavierStokesEquations<2>::NavierStokesBodyForceFctPt body_force_fct_pt);

 /// Destructor: Empty
 ~QuarterCircleProblem() {}

 /// Update the after solve (empty)
 void actions_after_newton_solve() {}

 /// \short Update the problem specs before solve. 
 /// (Re-)set velocity boundary conditions just to be on the safe side...
 void actions_before_newton_solve()
  { 
  // Setup tangential flow along boundary 1:
  unsigned ibound=Curved_bound; 
  unsigned num_nod= mesh_pt()->nboundary_node(ibound);
  for (unsigned inod=0;inod<num_nod;inod++)
   {
    // get coordinates
    double x=mesh_pt()->boundary_node_pt(ibound,inod)->x(0);
    double y=mesh_pt()->boundary_node_pt(ibound,inod)->x(1);
    // find Lagrangian coordinate (the angle)
    double zeta=0.0;
    if (x!=0.0)
     {
      zeta=atan(y/x);
     }
    // Tangential flow u0
    unsigned i=0;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,-sin(zeta));
    // Tangential flow u1
    i=1;
    mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,cos(zeta));
   }

  // Overwrite with no flow along all boundaries
  unsigned num_bound = mesh_pt()->nboundary();
  for(unsigned ibound=0;ibound<num_bound;ibound++)
   {
    if (ibound!=Curved_bound)
     {
      unsigned num_nod= mesh_pt()->nboundary_node(ibound);
      for (unsigned inod=0;inod<num_nod;inod++)
       {
        for (unsigned i=0;i<2;i++)
         {
          mesh_pt()->boundary_node_pt(ibound,inod)->set_value(i,0.0);
         }
       }
     }
   }

  } // end_of_actions_before_newton_solve


 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(mesh_pt()->element_pt());

   // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(mesh_pt()->element_pt());
   
   // Now pin the first pressure dof in the first element and set it to 0.0
   fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt
 
 /// Doc the solution
 void doc_solution();
 
private:

 /// Pointer to body force function
 NavierStokesEquations<2>::NavierStokesBodyForceFctPt Body_force_fct_pt;

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure


 Mesh* Bulk_mesh_pt;
 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_T_pt;
 Mesh* Surface_mesh_P_pt;

 // Preconditioner
 Preconditioner* Prec_pt;
 // Solver
 IterativeLinearSolver* Solver_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;

 unsigned Curved_bound;
 unsigned Left_bound;
 unsigned Bottom_bound;



}; // end_of_problem_class



//==start_of_constructor==================================================
/// Constructor for driven cavity problem in quarter circle domain
//========================================================================
template<class ELEMENT>
QuarterCircleProblem<ELEMENT>::QuarterCircleProblem(
 NavierStokesEquations<2>::NavierStokesBodyForceFctPt body_force_fct_pt) :
 Body_force_fct_pt(body_force_fct_pt)
{ 
  // Alias the namespace for convenience
  namespace NSPP = NavierStokesProblemParameters;
  namespace LPH = LagrangianPreconditionerHelpers;
  namespace QCL = QuarterCircleLagrange;

  // Set boundary IDs
  Left_bound = 2;
  Curved_bound = 1;
  Bottom_bound = 0;


 // Build geometric object that parametrises the curved boundary
 // of the domain

 // Half axes for ellipse
 double a_ellipse=1.0;
 double b_ellipse=1.0;

 // Setup elliptical ring 
 GeomObject* Wall_pt=new Ellipse(a_ellipse,b_ellipse); 

 // End points for wall
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);

 //Now create the mesh
 double fract_mid=0.5;
// Problem::mesh_pt() = new 
//  RefineableQuarterCircleSectorMesh<ELEMENT>(
//   Wall_pt,xi_lo,fract_mid,xi_hi);

 Bulk_mesh_pt = new 
  RefineableQuarterCircleSectorMesh<ELEMENT>(
   Wall_pt,xi_lo,fract_mid,xi_hi);

 add_sub_mesh(Bulk_mesh_pt);

 build_global_mesh();

 // Set error estimator
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<RefineableQuarterCircleSectorMesh<ELEMENT>*>(
  Bulk_mesh_pt)->spatial_error_estimator_pt()=error_estimator_pt;
 
 if(NSPP::Prob_id == QCL::PID_QC_VA_DRIVEN)
 {
 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here: All boundaries are Dirichlet boundaries.
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       mesh_pt()->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries
 }
 else
 {
   std::ostringstream err_msg;
   err_msg << "No such boundary conditions set for Problem id:\n"
           << "NSPP::Prob_id: " << NSPP::Prob_id << "\n"
           << std::endl;
   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION);
 }

 //Find number of elements in mesh
 unsigned n_element = mesh_pt()->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &NSPP::Rey;
   //Set the Re/Fr
   //NOTE: We take Fr as 1, therefore re_invFr  = Rey
   el_pt->re_invfr_pt() = &NSPP::Rey;
   //Set Gravity vector
   el_pt->g_pt() = &QCL::Gravity;
   //set body force function
   el_pt->body_force_fct_pt() = Body_force_fct_pt;

  } // end loop over elements
 
 // Initial refinement level
 for (unsigned ref_i = 0; ref_i < QCL::Noref; ref_i++) 
 {
 refine_uniformly();
 }

 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(mesh_pt()->element_pt());
 
 // Now pin the first pressure dof in the first element and set it to 0.0
 fix_pressure(0,0,0.0);
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 



 Vector<Mesh*> mesh_pt;
 mesh_pt.resize(1,0);
 mesh_pt[0] = Bulk_mesh_pt;

 LPH::Mesh_pt = mesh_pt;
 LPH::Problem_pt = this;
 Prec_pt = LPH::get_preconditioner();


 const double solver_tol = 1.0e-8;
 const double newton_tol = 1.0e-8;

 GenericProblemSetup::setup_solver(NSPP::Max_solver_iteration,
                                   solver_tol,newton_tol,
                                   NSPP::Solver_type,this,Prec_pt);


// // Create the preconditioner
// {
//   NavierStokesSchurComplementPreconditioner* ns_preconditioner_pt =
//      new NavierStokesSchurComplementPreconditioner(this);
//
//   Prec_pt = ns_preconditioner_pt;
//
//   ns_preconditioner_pt->set_navier_stokes_mesh(Bulk_mesh_pt);
//
//      TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
//      trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
//      Solver_pt = trilinos_solver_pt;
// }
//
// Solver_pt->preconditioner_pt() = Prec_pt;
// this->linear_solver_pt() = Solver_pt;


} // end_of_constructor



//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
  template<class ELEMENT>
void QuarterCircleProblem<ELEMENT>::doc_solution()
{ 

  namespace NSPP = NavierStokesProblemParameters;

  std::ofstream some_file;
  std::stringstream filename;
  filename << NSPP::Soln_dir_str<<"/"<<NSPP::Label_str<<".dat";

  // Number of plot points
  const unsigned npts=5;

  some_file.open(filename.str().c_str());
  Bulk_mesh_pt->output(some_file,npts);
  some_file.close();

} // end_of_doc_solution




//==start_of_main======================================================
/// Driver for QuarterCircleProblem test problem 
//=====================================================================
int main(int argc, char* argv[])
{
#ifdef OOMPH_HAS_MPI
  // Initialise MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Alias the namespace for convenience.
  namespace NSPP = NavierStokesProblemParameters;
  namespace LPH = LagrangianPreconditionerHelpers;
  namespace QCL = QuarterCircleLagrange;

  // Problem dimension.
  const unsigned dim = 2;

  // Set up doc info - used to store information on solver and iteration time.
  DocLinearSolverInfo doc_linear_solver_info;
  // Again, pass this to the NSPP and LPH
  NSPP::Doc_linear_solver_info_pt = &doc_linear_solver_info;
  LPH::Doc_linear_solver_info_pt = &doc_linear_solver_info;

  // Set the Label_pt
  LPH::Label_str_pt = &NSPP::Label_str;
  LPH::Vis_pt = &NSPP::Vis;
  QCL::Prob_id_pt = &NSPP::Prob_id;

  // Store commandline arguments
  CommandLineArgs::setup(argc,argv);

  NSPP::setup_commandline_flags();
  LPH::setup_commandline_flags();
  QCL::setup_commandline_flags();

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////

  // dim = 2
  NSPP::generic_problem_setup(dim);
  LPH::generic_setup();
  QCL::generic_setup();


  // Set up zero-Gravity vector
  QCL::Gravity[0] = 0.0;
  QCL::Gravity[1] = 0.0;

  // Set up Gamma vector for simplified form
//  NavierStokesEquations<2>::Gamma[0]=0;
//  NavierStokesEquations<2>::Gamma[1]=0;

  // Build problem with body force function and simplified form,
  // using body force function
  QuarterCircleProblem<RefineableQTaylorHoodElement<2> >
    problem(&QCL::zero_body_force);

  NSPP::Label_str = "soln0";

  // Solve the problem with automatic adaptation
  problem.newton_solve();

  if(NSPP::Doc_soln)
  {
  // Output solution
    problem.doc_solution();
  }

#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main


