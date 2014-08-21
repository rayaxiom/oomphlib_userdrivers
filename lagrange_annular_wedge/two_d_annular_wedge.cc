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
 class PartialAnnulusMesh : public RectangularQuadMesh<ELEMENT>
 {
 public:

  /// Constructor.
  PartialAnnulusMesh(const unsigned& nx, const unsigned& ny,
                  const double& lx,  const double& ly) :
   RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly)
   {
    // Alias the namespace for convenience.
    namespace AWL = AnnularWedgeLagrange;

    // For converting degrees to radians.
    const double deg_to_rad_ratio = MathematicalConstants::Pi / 180.0;

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
      
      // Map x to r
      double r = AWL::R_lo + (AWL::R_hi - AWL::R_lo) * x;

      // Map y to the angle phi
      double phi = (AWL::Phi_lo + (AWL::Phi_hi - AWL::Phi_lo)*y) 
                   * deg_to_rad_ratio;

      // Set new nodal coordinates
      nod_pt->x(0)=r * cos(phi);
      nod_pt->x(1)=r * sin(phi);
     }
   }
 };

} // end of namespace oomph


//===start_of_problem_class=============================================
//======================================================================

template<class ELEMENT>
class PartialAnnulusProblem : public Problem
{
public:

 /// \short Constructor: Pass number of elements in x and y directions and
 /// lengths
 PartialAnnulusProblem();

 /// Update before solve is empty
 void actions_before_newton_solve()
 {
   namespace NSPP = NavierStokesProblemParameters;

   // Initialise counters for each newton solve.
   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
   {
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
   namespace AWL = AnnularWedgeLagrange;

   if(NavierStokesProblemParameters::Distribute_problem)
   {
   if(NSPP::Prob_id == AWL::PID_AW_PO)
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
   namespace AWL = AnnularWedgeLagrange;

     if(NavierStokesProblemParameters::Distribute_problem)
   {
   if(NSPP::Prob_id == AWL::PID_AW_PO)
     {
       create_parall_outflow_lagrange_elements(
           1,Bulk_mesh_pt,Surface_mesh_P_pt);
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

 ///Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to full element type and fix the pressure at that element
   dynamic_cast<NavierStokesEquations<2>*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end of fix_pressure

 DocLinearSolverInfo* doc_linear_solver_info_pt()
 {
   return Doc_linear_solver_info_pt;
 }


private:

 double get_radial_v(const double& radial_distance);
 std::pair<double,double> get_radial_v(const Node* nod_pt);

 void set_mesh_bc_for_AwPo();
 void set_mesh_bc_for_AwPo_left_partial();
 void set_mesh_bc_for_AwPo_left_dirichlet();
 void set_mesh_bc_for_AwPo_left_neumann();
 void set_mesh_bc_for_AwPo_bottom_partial();
 void set_mesh_bc_for_AwPo_bottom_dirichlet();
 void set_mesh_bc_for_AwPo_bottom_neumann();
 void set_mesh_bc_for_AwPo_inflow();

 void set_mesh_bc_for_AwTf();
 void set_mesh_bc_for_AwTf_dirichlet();
 void set_mesh_bc_for_AwTmp();

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
PartialAnnulusProblem<ELEMENT>::PartialAnnulusProblem()
{
  // Alias the namespace for convenience
  namespace NSPP = NavierStokesProblemParameters;
  namespace LPH = LagrangianPreconditionerHelpers;
  namespace AWL = AnnularWedgeLagrange;

  Bottom_bound = 0;
  Right_bound = 1;
  Top_bound = 2;
  Left_bound = 3;

  Doc_linear_solver_info_pt = NSPP::Doc_linear_solver_info_pt;

  /// Setup the mesh
  // # of elements in x-direction
  unsigned nx=AWL::Noel;

  // # of elements in y-direction
  unsigned ny=AWL::Noel;

  // Domain length in x-direction
  double lx=AWL::Lx;

  // Domain length in y-direction
  double ly=AWL::Ly;

  Bulk_mesh_pt =
    new PartialAnnulusMesh<ELEMENT>(nx,ny,lx,ly);

  // Create the surface meshes
  if(NSPP::Prob_id == AWL::PID_AW_PO)
  {
    set_mesh_bc_for_AwPo();
  }
  else if(NSPP::Prob_id == AWL::PID_AW_TF)
  {
    set_mesh_bc_for_AwTf();
  }


  // Add the two meshes to the problem.
  add_sub_mesh(Bulk_mesh_pt);

  if(NSPP::Prob_id == AWL::PID_AW_PO)
  {
    add_sub_mesh(Surface_mesh_P_pt);
  }
  else if(NSPP::Prob_id == AWL::PID_AW_TF)
  {
    add_sub_mesh(Surface_mesh_T_pt);
  }

  // combine all sub-meshes into a single mesh.
  build_global_mesh();

  // Do the rest of the boundary conditions
  if(NSPP::Prob_id == AWL::PID_AW_PO)
  {
    if(AWL::BC_setting == 0)
    {
      set_mesh_bc_for_AwPo_bottom_partial();
      set_mesh_bc_for_AwPo_left_partial();
    }
    else if(AWL::BC_setting == 1)
    {
      set_mesh_bc_for_AwPo_bottom_partial();
      set_mesh_bc_for_AwPo_left_dirichlet();
    }
    else if(AWL::BC_setting == 2)
    {
      set_mesh_bc_for_AwPo_bottom_dirichlet();
      set_mesh_bc_for_AwPo_left_partial();
    }
    else if(AWL::BC_setting == 3)
    {
      set_mesh_bc_for_AwPo_bottom_dirichlet();
      set_mesh_bc_for_AwPo_left_dirichlet();
    }
    else if(AWL::BC_setting == 4)
    {
      set_mesh_bc_for_AwPo_bottom_neumann();
      set_mesh_bc_for_AwPo_left_neumann();
    }

    set_mesh_bc_for_AwPo_inflow();
  }
  else if(NSPP::Prob_id == AWL::PID_AW_TF)
  {
    set_mesh_bc_for_AwTf_dirichlet();
  }
  else if(NSPP::Prob_id == AWL::PID_AW_TMP)
  {
    set_mesh_bc_for_AwTmp();
  }

  //fix_pressure(0,0,0.0);

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

  if(NSPP::Prob_id == AWL::PID_AW_PO)
  {
  mesh_pt.resize(2,0);
  mesh_pt[0] = Bulk_mesh_pt;
  mesh_pt[1] = Surface_mesh_P_pt;
  }
  else if(NSPP::Prob_id == AWL::PID_AW_TF)
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
}

//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::doc_solution()
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
double PartialAnnulusProblem<ELEMENT>::get_radial_v(const double& radial_distance)
{
  namespace AWL = AnnularWedgeLagrange;

  return AWL::V_loR_lo * (1.0/radial_distance);
} // set_mesh_bc_for_AwPo

template<class ELEMENT>
std::pair<double, double> PartialAnnulusProblem<ELEMENT>::get_radial_v(const Node* nod_pt)
{
  const double x0 = nod_pt->x(0);
  const double x1 = nod_pt->x(1);
  const double radial_dist = sqrt(x0*x0 + x1*x1);

  const double phi = atan2(x1,x0);

  const double radial_v = get_radial_v(radial_dist);

  std::cout << "x = " << x0 << ", y = " << x1 << ", radial v: " << radial_v << std::endl; 
  

  std::pair <double,double> u_pair(radial_v*cos(phi),radial_v*sin(phi));

  return u_pair;
} // set_mesh_bc_for_AwPo


template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwPo()
{

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

} // set_mesh_bc_for_AwPo

template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwTf()
{
  namespace AWL = AnnularWedgeLagrange;

  // Assign the boundaries:
  //             2 inflow
  //           ----------
  //           |        |
  // 3 non slip|        |1 TF
  //           |        |
  //           ----------
  //             0 outflow
//  unsigned if_b = 3; // inflow
  const unsigned tf_b = 1; // parallel outflow

  // Create a "surface mesh" that will contain only
  // ImposeParallelOutflowElements in boundary 1
  // The constructor just creates the mesh without
  // giving it any elements, nodes, etc.
  Surface_mesh_T_pt = new Mesh;

  // Create ImposeParallelOutflowElement from all elements that are
  // adjacent to the Neumann boundary.
//  create_parall_outflow_lagrange_elements(po_b,
//                                          Bulk_mesh_pt,Surface_mesh_P_pt);
  create_impenetrable_lagrange_elements(tf_b,
                                        Bulk_mesh_pt,Surface_mesh_T_pt);
} // set_mesh_bc_for_AwTf

template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwTf_dirichlet()
{
  namespace AWL = AnnularWedgeLagrange;

  // Which boundary are we dealing with?
  unsigned current_bound = 1337;

  // The number of the nodes on a boundary.
  unsigned num_nod = 1337;

  // Inflow
  current_bound = 2;
  num_nod = mesh_pt()->nboundary_node(current_bound);
  for (unsigned inod = 0; inod < num_nod; inod++) 
  {
    Node* nod_pt = mesh_pt()->boundary_node_pt(current_bound,inod);

    nod_pt->pin(0);
    nod_pt->pin(1);

    // Get the coordinates
//    const double x0 = nod_pt->x(0);
    const double x1 = nod_pt->x(1);

    // Calculate the velocity
    const double u0 = (x1 - AWL::R_lo) * (2*AWL::R_hi - x1);

    // Set the velocity, push to the left!
    nod_pt->set_value(0,u0);
    nod_pt->set_value(1,0.0);
  }

  // left / non-slip boundary
  current_bound = 3;
  num_nod = mesh_pt()->nboundary_node(current_bound);
  for (unsigned inod = 0; inod < num_nod; inod++) 
  {
    Node* nod_pt = mesh_pt()->boundary_node_pt(current_bound,inod);

    nod_pt->pin(0);
    nod_pt->pin(1);

    nod_pt->set_value(0,0.0);
    nod_pt->set_value(1,0.0);
  }


  // Bottom/outflow
  current_bound = 0;
  num_nod = mesh_pt()->nboundary_node(current_bound);
  for (unsigned inod = 0; inod < num_nod; inod++) 
  {
    Node* nod_pt = mesh_pt()->boundary_node_pt(current_bound,inod);

    nod_pt->pin(0);
    nod_pt->pin(1);

    // Get the coordinates
    const double x0 = nod_pt->x(0);
//    const double x1 = nod_pt->x(1);

    // Calculate the velocity
    const double u1 = (-1) * (x0 - AWL::R_lo) * (2*AWL::R_hi - x0);

    // Set the velocity, push to the left!
    nod_pt->set_value(0,0.0);
    nod_pt->set_value(1,u1);
  }
} // set_mesh_bc_for_AwTf_dirichlet()

template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwTmp()
{
  namespace AWL = AnnularWedgeLagrange;

  // Which boundary are we dealing with?
  unsigned current_bound = 1337;

  // The number of the nodes on a boundary.
  unsigned num_nod = 1337;

  // Inflow
  current_bound = 2;
  num_nod = mesh_pt()->nboundary_node(current_bound);
  for (unsigned inod = 0; inod < num_nod; inod++) 
  {
    Node* nod_pt = mesh_pt()->boundary_node_pt(current_bound,inod);

    nod_pt->unpin(0);
    nod_pt->pin(1);

    // Get the coordinates
//    const double x0 = nod_pt->x(0);
//    const double x1 = nod_pt->x(1);

    // Calculate the velocity
    // R_lo = 1.0, R_hi = 3.
    // So x1 = 1, (1 - 1) * (2*3 - 1) = 0
    //    x1 = 3, (3 - 1) * (2*3 - 3) = 2*3 = 6
    // So we divide by 6 to normalise it.
    //const double u0 = ((x1 - AWL::R_lo) * (2*AWL::R_hi - x1)) / 6.0;

    // Set the velocity, push to the left!
//    nod_pt->set_value(0,u0);
    nod_pt->set_value(1,0.0);
  }

  // Bottom/outflow
  current_bound = 0;
  num_nod = mesh_pt()->nboundary_node(current_bound);
  for (unsigned inod = 0; inod < num_nod; inod++) 
  {
    Node* nod_pt = mesh_pt()->boundary_node_pt(current_bound,inod);

    nod_pt->pin(0);
    nod_pt->unpin(1);

    // Get the coordinates
//    const double x0 = nod_pt->x(0);
//    const double x1 = nod_pt->x(1);

    // Calculate the velocity
//    const double u1 = (-1) * (x0 - AWL::R_lo) * (AWL::R_hi - x0);

    // Set the velocity, push to the left!
    nod_pt->set_value(0,0.0);
//    nod_pt->set_value(1,u1);
  }

  // left / non-slip boundary
  current_bound = 3;
  num_nod = mesh_pt()->nboundary_node(current_bound);
  for (unsigned inod = 0; inod < num_nod; inod++) 
  {
    Node* nod_pt = mesh_pt()->boundary_node_pt(current_bound,inod);

    nod_pt->pin(0);
    nod_pt->pin(1);

    const double x0 = nod_pt->x(0);
    const double x1 = nod_pt->x(1);

    const double x0_x1_total = x0 + x1;

    const double u0 = (x0 / x0_x1_total);
    const double u1 = - (x1 / x0_x1_total);

//    nod_pt->set_value(0,0.0);
//    nod_pt->set_value(1,0.0);

    nod_pt->set_value(0,u0);
    nod_pt->set_value(1,u1);

  }

  // right / non-slip boundary
  current_bound = 1;
  num_nod = mesh_pt()->nboundary_node(current_bound);
  for (unsigned inod = 0; inod < num_nod; inod++) 
  {
    Node* nod_pt = mesh_pt()->boundary_node_pt(current_bound,inod);

    nod_pt->pin(0);
    nod_pt->pin(1);

    const double x0 = nod_pt->x(0);
    const double x1 = nod_pt->x(1);

    // Since -x/y is the gradient of the tangent to the circumference, this
    // gives us the ratio of ux and uy components!
    // We want to impose the value 1, so....
    const double x0_x1_total = x0 + x1;

    const double u0 = (x0 / x0_x1_total);
    const double u1 = - (x1 / x0_x1_total);

//    std::cout << "x: " << x0 << ", y: " << x1 << ", u: " << u0 << ", v:" << u1 << std::endl;
    

    nod_pt->set_value(0,u0);
    nod_pt->set_value(1,u1);
  }


} // set_mesh_bc_for_AwTmp()

template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwPo_bottom_partial()
{

 // Which boundary are we dealing with?
 unsigned current_bound = 1337;
 
 // The number of nodes on a boundary.
 unsigned num_nod = 1337;

 // bottom boundary: pin y, leave x free.
 current_bound = 0;
 num_nod = mesh_pt()->nboundary_node(current_bound);
 for (unsigned inod = 0; inod < num_nod; inod++) 
 {
   Node* nod_pt = mesh_pt()->boundary_node_pt(current_bound,inod);
   
   nod_pt->unpin(0);
   nod_pt->pin(1);

   std::pair<double,double>u_pair = get_radial_v(nod_pt);

   nod_pt->set_value(1,0.0);
 }

} // set_mesh_bc_for_AwPo

template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwPo_bottom_dirichlet()
{

 // Which boundary are we dealing with?
 unsigned current_bound = 1337;
 
 // The number of nodes on a boundary.
 unsigned num_nod = 1337;

 // bottom boundary: pin y, leave x free.
 current_bound = 0;
 num_nod = mesh_pt()->nboundary_node(current_bound);
 for (unsigned inod = 0; inod < num_nod; inod++) 
 {
   Node* nod_pt = mesh_pt()->boundary_node_pt(current_bound,inod);
   
   nod_pt->pin(0);
   nod_pt->pin(1);

   const double x0 = nod_pt->x(0);

   const double u = 1.0/x0;

   nod_pt->set_value(0,u);
   nod_pt->set_value(1,0);
 }
} // set_mesh_bc_for_AwPo


template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwPo_bottom_neumann()
{

 // Which boundary are we dealing with?
 unsigned current_bound = 1337;
 
 // The number of nodes on a boundary.
 unsigned num_nod = 1337;

 // bottom boundary: pin y, leave x free.
 current_bound = 0;
 num_nod = mesh_pt()->nboundary_node(current_bound);
 for (unsigned inod = 0; inod < num_nod; inod++) 
 {
   Node* nod_pt = mesh_pt()->boundary_node_pt(current_bound,inod);
   
   nod_pt->unpin(0);
   nod_pt->unpin(1);
 }
} // set_mesh_bc_for_AwPo


template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwPo_left_partial()
{
 // Which boundary are we dealing with?
 unsigned current_bound = 1337;
 
 // The number of nodes on a boundary.
 unsigned num_nod = 1337;

 // left boundary: pin x, leave y free.
 current_bound = 2;
 num_nod=Bulk_mesh_pt->nboundary_node(current_bound);
 for (unsigned inod=0;inod<num_nod;inod++)
 {
   // Get node
   //Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);

   // Get node
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);
   
   nod_pt->pin(0);
   nod_pt->unpin(1);

//   double y=nod_pt->x(1);

//   double u = 1.0/y;

   std::pair<double,double> u_pair = get_radial_v(nod_pt);

   nod_pt->set_value(0,0.0);
//   nod_pt->set_value(1,u);
 }

} // set_mesh_bc_for_AwPo

template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwPo_left_dirichlet()
{
 // Which boundary are we dealing with?
 unsigned current_bound = 1337;
 
 // The number of nodes on a boundary.
 unsigned num_nod = 1337;

 // left boundary: pin x, leave y free.
 current_bound = 2;
 num_nod=Bulk_mesh_pt->nboundary_node(current_bound);
 for (unsigned inod=0;inod<num_nod;inod++)
 {
   // Get node
   //Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);

   // Get node
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);
   
   nod_pt->pin(0);
   nod_pt->pin(1);

   double y=nod_pt->x(1);

   double u = 1.0/y;

   nod_pt->set_value(0,0);
   nod_pt->set_value(1,u);
 }

} // set_mesh_bc_for_AwPo

template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwPo_left_neumann()
{
 // Which boundary are we dealing with?
 unsigned current_bound = 1337;
 
 // The number of nodes on a boundary.
 unsigned num_nod = 1337;

 // left boundary: pin x, leave y free.
 current_bound = 2;
 num_nod=Bulk_mesh_pt->nboundary_node(current_bound);
 for (unsigned inod=0;inod<num_nod;inod++)
 {
   // Get node
   //Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);

   // Get node
   Node* nod_pt=mesh_pt()->boundary_node_pt(current_bound,inod);
   
   nod_pt->unpin(0);
   nod_pt->unpin(1);
 }

} // set_mesh_bc_for_AwPo

template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::set_mesh_bc_for_AwPo_inflow()
{

  pause("About to do inflow"); 
  
 // Which boundary are we dealing with?
 unsigned current_bound = 1337;
 
 // The number of nodes on a boundary.
 unsigned num_nod = 1337;

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
//   const double x0=nod_pt->x(0);
//   const double x1=nod_pt->x(1);
//   const double phi = atan2(x1,x0);

//   const double u = 1.0;

   std::pair<double,double> u_pair = get_radial_v(nod_pt);





   nod_pt->set_value(0,u_pair.first);
   nod_pt->set_value(1,u_pair.second);
 }

 pause("done inflow"); 
 
} // set_mesh_bc_for_AwPo


//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
template<class ELEMENT>
void PartialAnnulusProblem<ELEMENT>::
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
void PartialAnnulusProblem<ELEMENT>::
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
  namespace AWL = AnnularWedgeLagrange;

  std::string label = AWL::prob_str()
                      + NSPP::create_label() 
                      + LPH::create_label() 
                      + AWL::noel_str();
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
  namespace AWL = AnnularWedgeLagrange;

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
  AWL::Prob_id_pt = &NSPP::Prob_id;

  // Store commandline arguments
  CommandLineArgs::setup(argc,argv);

  NSPP::setup_commandline_flags();
  LPH::setup_commandline_flags();
  AWL::setup_commandline_flags();

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////

  // dim = 2
  NSPP::generic_problem_setup(dim);
  LPH::generic_setup();
  AWL::generic_setup();

  // Solve with Taylor-Hood element, set up problem
  PartialAnnulusProblem< QTaylorHoodElement<dim> > problem;

  if(NavierStokesProblemParameters::Distribute_problem)
  {
    problem.distribute();
  }

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
