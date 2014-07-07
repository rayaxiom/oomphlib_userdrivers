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

   namespace NSPP = NavierStokesProblemParameters;
   namespace QCL = QuarterCircleLagrange;

   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
   {
    // Initialise counters for each newton solve.
    Doc_linear_solver_info_pt->setup_new_time_step();
   }

   unsigned current_bound;
   unsigned num_nod;

   if(NSPP::Prob_id == QCL::PID_QC_TF_ALLDIRI)
   {
     // Outflow: bottom boundary
     current_bound = 0;
     num_nod = Bulk_mesh_pt->nboundary_node(current_bound);
     for (unsigned inod = 0; inod < num_nod; inod++) 
     {
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(current_bound,inod);

       nod_pt->pin(0);
       nod_pt->pin(1);

       const double x0 = nod_pt->x(0);

       const double u0 = 0.0;
       const double u1 = -x0*(2.0-x0);

       nod_pt->set_value(0,u0);
       nod_pt->set_value(1,u1);
       //      nod_pt->set_value(1,u1);
     }

     // Inflow, left boundary, boundary 2
     current_bound = 2;
     num_nod = Bulk_mesh_pt->nboundary_node(current_bound);
     for (unsigned inod = 0; inod < num_nod; inod++) 
     {
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(current_bound,inod);

       nod_pt->pin(0);
       nod_pt->pin(1);

       const double x1 = nod_pt->x(1);

       const double u0 = x1 * (2.0 - x1);
       const double u1 = 0.0;

       nod_pt->set_value(0,u0);
       nod_pt->set_value(1,u1);

     }
   }
   else if(NSPP::Prob_id == QCL::PID_QC_TF_BOTNEU)
   {
     // Outflow: bottom boundary
     current_bound = 0;
     num_nod = Bulk_mesh_pt->nboundary_node(current_bound);
     for (unsigned inod = 0; inod < num_nod; inod++) 
     {
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(current_bound,inod);

       nod_pt->pin(0);
       nod_pt->unpin(1);

       nod_pt->set_value(0,0.0);
     }

     // Inflow, left boundary, boundary 2
     current_bound = 2;
     num_nod = Bulk_mesh_pt->nboundary_node(current_bound);
     for (unsigned inod = 0; inod < num_nod; inod++) 
     {
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(current_bound,inod);

       nod_pt->pin(0);
       nod_pt->pin(1);

       const double x1 = nod_pt->x(1);

       const double u0 = x1 * (2.0 - x1);
       const double u1 = 0.0;

       nod_pt->set_value(0,u0);
       nod_pt->set_value(1,u1);

     }
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

 } // end_of_actions_before_newton_solve

 void actions_after_newton_step()
 {
   namespace NSPP = NavierStokesProblemParameters;
   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
   {
     NSPP::doc_iter_times(this,Doc_linear_solver_info_pt);
   }
 }

 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {

    create_impenetrable_lagrange_elements(Curved_bound,
                                          Bulk_mesh_pt,
                                          Surface_mesh_T_pt);

    rebuild_global_mesh();

   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(Bulk_mesh_pt->element_pt());

   // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
   
   // Now pin the first pressure dof in the first element and set it to 0.0
   //fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt

 void actions_before_adapt()
 {
   // Kill the flux elements and wipe the surface mesh
   delete_impenetrable_lagrange_elements(Surface_mesh_T_pt);
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh(); 
 }

 void create_impenetrable_lagrange_elements(const unsigned &b,
                                            Mesh* const &bulk_mesh_pt,
                                            Mesh* const &surface_mesh_pt); 

 void delete_impenetrable_lagrange_elements(Mesh* const &surface_mesh_pt);
 
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
   dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure


 Mesh* Bulk_mesh_pt;
// RefineableQuarterCircleSectorMesh<ELEMENT>* Bulk_mesh_pt;
 /// Pointer to the "surface" mesh
 Mesh* Surface_mesh_T_pt;

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

  Doc_linear_solver_info_pt = NSPP::Doc_linear_solver_info_pt;

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


 Surface_mesh_T_pt = new Mesh;

 create_impenetrable_lagrange_elements(Curved_bound,
                                       Bulk_mesh_pt,
                                       Surface_mesh_T_pt);


 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_T_pt);
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
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  } // end loop over boundaries
 }
 else if(NSPP::Prob_id == QCL::PID_QC_TF_ALLDIRI)
 {
   unsigned current_bound;
   unsigned num_nod;

   // Pin both the x and y velocity of the bottom and left boundary.
   {
     current_bound = Bottom_bound;
     num_nod = Bulk_mesh_pt->nboundary_node(current_bound);
     for (unsigned inod = 0; inod < num_nod; inod++) 
     {
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(current_bound,inod);

       nod_pt->pin(0);
       nod_pt->pin(1);
     }
   }

   {
     current_bound = Left_bound;
     num_nod = Bulk_mesh_pt->nboundary_node(current_bound);
     for (unsigned inod = 0; inod < num_nod; inod++) 
     {
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(current_bound,inod);

       // Pin both velocity components
       nod_pt->pin(0);
       nod_pt->pin(1);

     }
   }


 }
 else if(NSPP::Prob_id == QCL::PID_QC_TF_BOTNEU)
   // The bottom boundary has ux pinned, uy free
 {
   unsigned current_bound;
   unsigned num_nod;

   // Pin both the x and y velocity of the bottom and left boundary.
   {
     current_bound = Bottom_bound;
     num_nod = Bulk_mesh_pt->nboundary_node(current_bound);
     for (unsigned inod = 0; inod < num_nod; inod++) 
     {
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(current_bound,inod);

       nod_pt->pin(0);
       nod_pt->unpin(1);
     }
   }

   {
     current_bound = Left_bound;
     num_nod = Bulk_mesh_pt->nboundary_node(current_bound);
     for (unsigned inod = 0; inod < num_nod; inod++) 
     {
       Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(current_bound,inod);

       // Pin both velocity components
       nod_pt->pin(0);
       nod_pt->pin(1);

     }
   }
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
 unsigned n_element = Bulk_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

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
  pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
 
 // Now pin the first pressure dof in the first element and set it to 0.0
 //fix_pressure(0,0,0.0);
 
 // Setup equation numbering scheme
 cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 



 // Only do this bit if we do NOT have a direct solver.
 if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
 {
   // Create the vector of mesh pointers!
   Vector<Mesh*> mesh_pt;
   if((NSPP::Prob_id == QCL::PID_QC_TF_ALLDIRI) ||
      (NSPP::Prob_id == QCL::PID_QC_TF_BOTNEU))
   {
     mesh_pt.resize(2,0);
     mesh_pt[0] = Bulk_mesh_pt;
     mesh_pt[1] = Surface_mesh_T_pt;
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

/////////////////////////////////////////////////////////////////////////

// Vector<Mesh*> mesh_pt;
// mesh_pt.resize(1,0);
// mesh_pt[0] = Bulk_mesh_pt;
//
// LPH::Mesh_pt = mesh_pt;
// LPH::Problem_pt = this;
// Prec_pt = LPH::get_preconditioner();
//
//
// const double solver_tol = 1.0e-6;
// const double newton_tol = 1.0e-6;
//
// GenericProblemSetup::setup_solver(NSPP::Max_solver_iteration,
//                                   solver_tol,newton_tol,
//                                   NSPP::Solver_type,this,Prec_pt);
//
//
/////////////////////////////////////////////////////////////////////////


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

//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void QuarterCircleProblem<ELEMENT>::
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

   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 0 or 2?
     if ((nod_pt->is_on_boundary(Left_bound))||
         (nod_pt->is_on_boundary(Bottom_bound)))
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

//============start_of_delete_flux_elements==============================
/// Delete Poisson Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void QuarterCircleProblem<ELEMENT>::
delete_impenetrable_lagrange_elements(Mesh* const &surface_mesh_pt)
{
 // How many surface elements are in the surface mesh
 unsigned n_element = surface_mesh_pt->nelement();
 
 std::cout << "No. elements on surface mesh: " << n_element << std::endl; 
 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete surface_mesh_pt->element_pt(e);
  }

 // Wipe the mesh
 surface_mesh_pt->flush_element_and_node_storage();
} // end of delete_flux_elements




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
  namespace QCL = QuarterCircleLagrange;

  std::string label = QCL::prob_str()
                      + NSPP::create_label() 
                      + LPH::create_label() 
                      + QCL::noref_str();
  return label;
}


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


  // Build problem with body force function and simplified form,
  // using body force function
  QuarterCircleProblem<RefineableQTaylorHoodElement<2> >
    problem(&QCL::zero_body_force);

  if(NavierStokesProblemParameters::Distribute_problem)
  {
    problem.distribute();
  }

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

      std::cout << "RAYDOING: "
        << NSPP::Label_str
        << " on " << ctime(&rawtime) << std::endl;

      // Solve the problem
      problem.newton_solve();

      //Output solution
      if(NSPP::Doc_soln)
      {problem.doc_solution();}

      //////////////////////////////////////////////////////////////////////////
      ////////////// Outputting results ////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////

      if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
      {
      // Get the global oomph-lib communicator 
      const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

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
      // (timing/iteration counts) as we get it since it will interlace with the
      // other processors and becomes hard to read.
      std::ostringstream results_stream;

      // Get the 3D vector which holds the iteration counts and timing results.
      Vector<Vector<Vector<double> > > iters_times
        = NSPP::Doc_linear_solver_info_pt->iterations_and_times();

      ResultsFormat::format_rayits(rey_increment,&iters_times,&results_stream);
      ResultsFormat::format_prectime(rey_increment,&iters_times,&results_stream);
      ResultsFormat::format_solvertime(rey_increment,&iters_times,&results_stream);

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
    //    NSPP::Label_str = NSPP::create_label() + LPH::create_label()+SL::create_label();
    NSPP::Label_str = create_label();

    time_t rawtime;
    time(&rawtime);

    std::cout << "RAYDOING: "
      << NSPP::Label_str
      << " on " << ctime(&rawtime) << std::endl;

    // Solve the problem
    problem.newton_solve();

    //Output solution
    if(NSPP::Doc_soln)
    {problem.doc_solution();}


    //////////////////////////////////////////////////////////////////////////
    ////////////// Outputting results ////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////

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
    // (timing/iteration counts) as we get it since it will interlace with the
    // other processors and becomes hard to read.
    std::ostringstream results_stream;

    // Get the 3D vector which holds the iteration counts and timing results.
    Vector<Vector<Vector<double> > > iters_times
      = NSPP::Doc_linear_solver_info_pt->iterations_and_times();

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

    // Now doing the preconditioner setup time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat::format_prectime(intimestep,&iters_times,&results_stream);
    }

    // Now doing the linear solver time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat::format_solvertime(intimestep,&iters_times,&results_stream);
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

////////////////////////////////////////////////////////////////////////
  // Solve the problem with automatic adaptation
  problem.newton_solve();

  if(NSPP::Doc_soln)
  {
  // Output solution
    problem.doc_solution();
  }
/////////////////////////////////////////////////////////////////////////
#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main


