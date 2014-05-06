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
   if(NavierStokesProblemParameters::Using_trilinos_solver)
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
   create_parall_outflow_lagrange_elements(5,Bulk_mesh_pt,Surface_mesh_P_pt);
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
 namespace NSPP = NavierStokesProblemParameters;
 namespace LPH = LagrangianPreconditionerHelpers;
 namespace SL = StepLagrange;
 
 Doc_linear_solver_info_pt = NSPP::Doc_linear_solver_info_pt;

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
   el_pt->re_pt() = &NSPP::Rey;

  } // for(unsigned e=0;e<n_el;e++)

 //Assgn equation numbers
 std::cout << "\n equation numbers : "<< assign_eqn_numbers() << std::endl;

 Vector<Mesh*> mesh_pt(2,0);
 mesh_pt[0] = Bulk_mesh_pt;
 mesh_pt[1] = Surface_mesh_P_pt;

 LPH::Mesh_pt = mesh_pt;
 LPH::Problem_pt = this;

 Prec_pt = LPH::get_preconditioner();


 // Build solve and preconditioner
#ifdef OOMPH_HAS_TRILINOS
 TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
 trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
 Solver_pt = trilinos_solver_pt;
 NSPP::Using_trilinos_solver = true;
#else
 Solver_pt = new GMRES<CRDoubleMatrix>;
 // We use RHS preconditioning. Note that by default,
 // left hand preconditioning is used.
 static_cast<GMRES<CRDoubleMatrix>*>(Solver_pt)->set_preconditioner_RHS();
 NSPP::Using_trilinos_solver = false;
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

template<class ELEMENT>
void BackwardStepProblem<ELEMENT>::
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
  namespace SL = StepLagrange;

  std::string label = SL::Prob_str 
                      + NSPP::create_label() 
                      + LPH::create_label() 
                      + SL::Ang_deg_str + SL::Noel_str;
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
  namespace NSPP = NavierStokesProblemParameters;
  namespace LPH = LagrangianPreconditionerHelpers;
  namespace SL = StepLagrange;

  // Get the global oomph-lib communicator 
  //  const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();

  // my rank and number of processors. This is used later for putting the data.
  //  unsigned my_rank = comm_pt->my_rank();
  //unsigned nproc = comm_pt->nproc();

  const unsigned dim = 2;


  // Set up doc info
  DocLinearSolverInfo doc_linear_solver_info;
  NSPP::Doc_linear_solver_info_pt = &doc_linear_solver_info;
  LPH::Doc_linear_solver_info_pt = &doc_linear_solver_info;

  // Set the Label_pt
  LPH::Label_str_pt = &NSPP::Label_str;
  LPH::Vis_pt = &NSPP::Vis;
  SL::Prob_id_pt = &NSPP::Prob_id;

  // Store commandline arguments
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

  BackwardStepProblem< QTaylorHoodElement<dim> > problem;


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
