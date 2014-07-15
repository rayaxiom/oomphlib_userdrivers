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
   namespace NSPP = NavierStokesProblemParameters;

   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
   {
     NSPP::doc_iter_times(this,Doc_linear_solver_info_pt);
   }
 }

 void actions_before_distribute()
 {
   namespace NSPP = NavierStokesProblemParameters;
   namespace SL = StepLagrange;

   if(NSPP::Distribute_problem)
   {
   if(NSPP::Prob_id == SL::PID_ST_PO)
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
   namespace SL = StepLagrange;

   if(NSPP::Distribute_problem)
   {
   if(NSPP::Prob_id == SL::PID_ST_PO)
   {
     create_parall_outflow_lagrange_elements(Left_b,Bulk_mesh_pt,Surface_mesh_P_pt);
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

private:

 void set_mesh_bc_for_StPo();
 void set_mesh_bc_for_StVa();

 // Boundaries, these are set in the problem initialiser list.
 const unsigned Left_b;
 const unsigned Rightmost_b;
 const unsigned Right_cutout_b;
 const unsigned Bottommost_b;
 const unsigned Bottom_cutout_b;
 const unsigned Top_b;

 /// Pointer to the "bulk" mesh
// SlopingQuadMesh<ELEMENT>* Bulk_mesh_pt;
 Mesh* Bulk_mesh_pt;

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
/// Problem constructor.
///
/// Boundary IDs:
///
///        4 Top_b
///      --------------------------
///Left_b|                        | 3 Rightmost_b
///  5   |                 _______|
///      |                |    2
///      |                | 1
///      ------------------
///        0 Bottommost_b
///
///====================================================================
template<class ELEMENT>
BackwardStepProblem<ELEMENT>::BackwardStepProblem() : 
  Left_b(5), Rightmost_b(3), Right_cutout_b(1),
  Bottommost_b(0), Bottom_cutout_b(2), Top_b(4)
{
 // Alias the namespace for convenience
 namespace NSPP = NavierStokesProblemParameters;
 namespace LPH = LagrangianPreconditionerHelpers;
 namespace SL = StepLagrange;

 Doc_linear_solver_info_pt = NSPP::Doc_linear_solver_info_pt;

 const unsigned if_b = Rightmost_b;
 //unsigned tf_b=1;
 const unsigned po_b = Left_b;

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

 if(NSPP::Prob_id == SL::PID_ST_PO)
 {
 Bulk_mesh_pt =
  new SlopingQuadMesh<ELEMENT>(nx,ny,nx_cut_out,ny_cut_out,lx,ly,SL::Ang);
 }
 else if(NSPP::Prob_id == SL::PID_ST_VA)
 {
   Bulk_mesh_pt = 
     new BackwardStepQuadMesh<ELEMENT>(nx,ny,nx_cut_out,ny_cut_out,lx,ly);
 }
 else
 {
       std::ostringstream err_msg;
   err_msg << "No mesh set up for problem id: "
           << NSPP::Prob_id << ".\n"
           << std::endl;

   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION); 
 }


  if(NSPP::Prob_id == SL::PID_ST_PO)
  {
    set_mesh_bc_for_StPo();
  }
  else if(NSPP::Prob_id == SL::PID_ST_VA)
  {
    set_mesh_bc_for_StVa();
  }
  else
  {

       std::ostringstream err_msg;
   err_msg << "No boundary conditions set up for problem id: "
           << NSPP::Prob_id << ".\n"
           << std::endl;

   throw OomphLibError(err_msg.str(),
       OOMPH_CURRENT_FUNCTION,
       OOMPH_EXCEPTION_LOCATION); 

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

 //Assign equation numbers
 std::cout << "\n equation numbers : "<< assign_eqn_numbers() << std::endl;

 if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
 {
 Vector<Mesh*> mesh_pt;
 if(NSPP::Prob_id == SL::PID_ST_PO)
 {
   mesh_pt.resize(2,0);
   mesh_pt[0] = Bulk_mesh_pt;
   mesh_pt[1] = Surface_mesh_P_pt;
 }
 else if(NSPP::Prob_id == SL::PID_ST_VA)
 {
   mesh_pt.resize(1,0);
   mesh_pt[0] = Bulk_mesh_pt;
 }

 LPH::Mesh_pt = mesh_pt;
 LPH::Problem_pt = this;
 Prec_pt = LPH::get_preconditioner();
 }
 const double solver_tol = 1.0e-6;
 const double newton_tol = 1.0e-6;
 GenericProblemSetup::setup_solver(NSPP::Max_solver_iteration,
                                   solver_tol, newton_tol,
                                   NSPP::Solver_type,this,Prec_pt);

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
void BackwardStepProblem<ELEMENT>::set_mesh_bc_for_StPo()
{
  namespace SL = StepLagrange;

 // Create a "surface mesh" that will contain only
 // ImposeParallelOutflowElements in boundary 1
 // The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
 Surface_mesh_P_pt = new Mesh;
 //Surface_mesh_T_pt = new Mesh;
 //
 // Create ImposeParallelOutflowElement from all elements that are
 // adjacent to the Neumann boundary.
 create_parall_outflow_lagrange_elements(Left_b,
                                         Bulk_mesh_pt,Surface_mesh_P_pt);
 //create_impenetrable_lagrange_elements(po_b,
 //                                        Bulk_mesh_pt,Surface_mesh_P_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh_P_pt);
 //add_sub_mesh(Surface_mesh_T_pt);
 // Combine all submeshes into a single Mesh
 build_global_mesh();



 const unsigned num_bound=Bulk_mesh_pt->nboundary();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 for(unsigned ibound=0;ibound<num_bound;ibound++)
 {
   //if((ibound != po_b)&&(ibound != tf_b))
   if(ibound != Left_b)
   {
     unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
     {
       // Get node
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

       nod_pt->pin(0);
       nod_pt->pin(1);

       nod_pt->set_value(0,0);
       nod_pt->set_value(1,0);

     }
   }
 }

 // Which boundary are we dealing with?
 unsigned current_bound;

 // The number of nodes on the current boundary.
 unsigned num_nod;

 // PARABOLIC INFLOW!
 // Inflow is at the Rightmost_b.
 current_bound = Rightmost_b;
 num_nod= mesh_pt()->nboundary_node(current_bound);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);

   // Pin both velocity components (just to be safe!)
   nod_pt->pin(0);
   nod_pt->pin(1);

   // Get the current Cartesian coordinates
   const double x0 = nod_pt->x(0);
   const double x1 = nod_pt->x(1);

   // When applying the parabolic inflow on a zero degrees tilt, it only
   // depends on the y coordinate. We tilt the y coordinate back, create
   // the velocity profile, them tilt this according to the angle.
   const double x1_old = x0*sin(-SL::Ang) + x1*cos(-SL::Ang);

   // Now calculate the parabolic inflow at this point.
   const double u0_old = -4.0*(x1_old - 1.0) * (2.0 - x1_old);

   // Now apply the rotation to u, using rotation matrices.
   // with x = u and y = 0, i.e. R*[u;0] since we have the
   // velocity in the x direction only. There is no velocity
   // in the y direction.
   const double u0 = u0_old*cos(SL::Ang);
   const double u1 = u0_old*sin(SL::Ang);

   nod_pt->set_value(0,u0);
   nod_pt->set_value(1,u1);
 } // Setting the inflow
} // set_mesh_bc_for_StPo

template<class ELEMENT>
void BackwardStepProblem<ELEMENT>::set_mesh_bc_for_StVa()
{
  namespace SL = StepLagrange;

 // Create a "surface mesh" that will contain only
 // ImposeParallelOutflowElements in boundary 1
 // The constructor just creates the mesh without
 // giving it any elements, nodes, etc.
// Surface_mesh_P_pt = new Mesh;
 //Surface_mesh_T_pt = new Mesh;
 //
 // Create ImposeParallelOutflowElement from all elements that are
 // adjacent to the Neumann boundary.
// create_parall_outflow_lagrange_elements(Left_b,
//                                         Bulk_mesh_pt,Surface_mesh_P_pt);
 //create_impenetrable_lagrange_elements(po_b,
 //                                        Bulk_mesh_pt,Surface_mesh_P_pt);

 // Add the two sub meshes to the problem
 add_sub_mesh(Bulk_mesh_pt);
// add_sub_mesh(Surface_mesh_P_pt);
 //add_sub_mesh(Surface_mesh_T_pt);
 // Combine all submeshes into a single Mesh
 build_global_mesh();

 // For the outflow boundary, we pin the y velocity to 0 and leave x.
 {
   const unsigned num_nod = mesh_pt()->nboundary_node(Left_b);
   for (unsigned inod = 0; inod < num_nod; inod++) 
   {
     Node* nod_pt = mesh_pt()->boundary_node_pt(Left_b,inod);
     
     // Unpin x
     nod_pt->unpin(0);

     // pin y
     nod_pt->pin(1);

     // Set y velocity to zero.
     nod_pt->set_value(1,0);
   }
 }

 const unsigned num_bound=Bulk_mesh_pt->nboundary();

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 for(unsigned ibound=0;ibound<num_bound;ibound++)
 {
   //if((ibound != po_b)&&(ibound != tf_b))
   if(ibound != Left_b)
   {
     unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
     {
       // Get node
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

       nod_pt->pin(0);
       nod_pt->pin(1);

       nod_pt->set_value(0,0);
       nod_pt->set_value(1,0);

     }
   }
 }

 // Which boundary are we dealing with?
 unsigned current_bound;

 // The number of nodes on the current boundary.
 unsigned num_nod;

 // PARABOLIC INFLOW!
 // Inflow is at the Rightmost_b.
 current_bound = Rightmost_b;
 num_nod= mesh_pt()->nboundary_node(current_bound);
 for(unsigned inod=0;inod<num_nod;inod++)
 {
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);

   // Pin both velocity components (just to be safe!)
   nod_pt->pin(0);
   nod_pt->pin(1);

   // Get the y Cartesian coordinates
   const double x1 = nod_pt->x(1);

   // Now calculate the parabolic inflow at this point.
   const double u0 = -4.0*(x1 - 1.0) * (2.0 - x1);

   const double u1 = 0.0;

   nod_pt->set_value(0,u0);
   nod_pt->set_value(1,u1);
 } // Setting the inflow
} // set_mesh_bc_for_StVa

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

  std::string label = SL::prob_str()
                      + NSPP::create_label() 
                      + LPH::create_label() 
                      + SL::ang_deg_str() + SL::noel_str();
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

  // Problem dimension
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

  if(NavierStokesProblemParameters::Distribute_problem)
  {
    problem.distribute();
  }

  ///////////////////////////////////////////////////////////////////////////////



  if(!CommandLineArgs::command_line_flag_has_been_set("--rey"))
  {
    unsigned rey_increment = 0; // used for output of iters/times

    for (NSPP::Rey = NSPP::Rey_start; 
        NSPP::Rey <= NSPP::Rey_end; NSPP::Rey += NSPP::Rey_incre) 
    {
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


#ifdef OOMPH_HAS_MPI
  // finalize MPI
  MPI_Helpers::finalize();
#endif
  return(EXIT_SUCCESS);
} // end_of_main
