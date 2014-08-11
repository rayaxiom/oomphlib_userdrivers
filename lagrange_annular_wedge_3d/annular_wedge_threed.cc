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

// Get the mesh
#include "meshes/tetgen_mesh.h"
#include "meshes/brick_from_tet_mesh.h"

// My own header
#include "./../rayheader.h"
using namespace std;

using namespace oomph;

// Alias the namespace for convenience.
namespace NSPP = NavierStokesProblemParameters;
namespace LPH = LagrangianPreconditionerHelpers;
namespace AW3DL = AnnularWedge3DLagrange;


void convert_to_melon(const double&x, const double& y, const double& z,
    Vector<double>& x_new)
{
  const double r = 1.0 + 2.0*x;
  const double phi = 90.0*y * (MathematicalConstants::Pi / 180.0);
  const double new_z = 2.0*z;

  x_new.resize(3,0);

  x_new[0] = r * cos(phi);
  x_new[1] = r * sin(phi);
  x_new[2] = new_z;
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
 class AnnularWedgeCubicMesh : public SimpleCubicMesh<ELEMENT>
 {
 public:

  /// Constructor.
  AnnularWedgeCubicMesh(const unsigned& nx, const unsigned& ny, const unsigned& nz,
                   const double& lx,  const double& ly, const double& lz,
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
      convert_to_melon(x,y,z,x_new);

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


//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Variables
{
 double Ang_deg = 30;

 double Ang = Ang_deg * (MathematicalConstants::Pi / 180.0);

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




// void get_prescribed_inflow_full(const double& t,
//                            const double& y,
//                            const double& z,
//                            double& ux)
// {
//   // For the velocity profile in the x direction.
//   // 1) First form the parabolic profile
//   ux = 0.0;
//   {
//     // 2) Now make it move in time
////     const double trig_scaling = 0.025;
////     const double ux_scaling = (1.0 
////                               - cos(trig_scaling
////                                     *MathematicalConstants::Pi
////                                     *t)) * 2.0;
//
//     const double ux_scaling = t / Time_end;
//     ux = (y-0.0)*(1.0-y)*(z-0.0)*(1.0-z) * ux_scaling;
////     ux = (y-0.5)*(1.0-y)*(z-0.5)*(1.0-z);
//   }
// } 

} // end_of_namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
namespace RAYRAY
{

  double get_radial_v(const double& radial_distance, const double& t)
  {
    const double scaling = -cos(MathematicalConstants::Pi*t)/2.0 + 0.5;
    //  namespace AWL = AnnularWedgeLagrange;

    //  return AWL::V_loR_lo * (1.0/radial_distance);
    return (1.0 / radial_distance) * scaling;
  } // set_mesh_bc_for_AwPo


std::pair<double, double> get_radial_v(const Node* nod_pt, const double& t)
{
  const double x0 = nod_pt->x(0);
  const double x1 = nod_pt->x(1);
  const double radial_dist = sqrt(x0*x0 + x1*x1);

  const double phi = atan2(x1,x0);

  const double radial_v = get_radial_v(radial_dist,t);

  std::pair <double,double> u_pair(radial_v*cos(phi),radial_v*sin(phi));

  return u_pair;
} // set_mesh_bc_for_AwPo



}

//==start_of_problem_class============================================
/// Test problem for Fp/PCD preconditioner
//====================================================================
template<class ELEMENT>
class CubeProblem : public Problem
{

public:


 /// Constructor
 CubeProblem();

 
 /// Destructor: Cleanup
 virtual ~CubeProblem()
  {
    GenericProblemSetup::clean_up_solver_memory();
    LPH::clean_up_memory();

    delete Surface_mesh_pt;
    delete Bulk_mesh_pt;
//   delete Prec_pt;
//   delete P_matrix_preconditioner_pt;
//   delete F_matrix_preconditioner_pt;
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
    Doc_linear_solver_info_pt->clear_current_time_step();

   {
    // Inflow in upper half of inflow boundary
    const unsigned ibound=Inflow_boundary; 
    const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

      const double time=time_pt()->time();


      std::pair<double,double>u_xy_pair = RAYRAY::get_radial_v(nod_pt,time);

      nod_pt->set_value(0,u_xy_pair.first);
      nod_pt->set_value(1,u_xy_pair.second);
      nod_pt->set_value(2,0.0);

     }
   }

//   {
//    // Inflow in upper half of inflow boundary
//    const unsigned ibound=YZ_boundary; 
//    const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
//    for (unsigned inod=0;inod<num_nod;inod++)
//     {
//      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
//
//      const double time=time_pt()->time();
//
//      std::pair<double,double>u_xy_pair = RAYRAY::get_radial_v(nod_pt,time);
//
//      nod_pt->set_value(0,0.0);
//      nod_pt->set_value(1,u_xy_pair.second);
//      nod_pt->set_value(2,0.0);
//     }
//   }
//   {
//    // Inflow in upper half of inflow boundary
//    const unsigned ibound=XZ_boundary; 
//    const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
//    for (unsigned inod=0;inod<num_nod;inod++)
//     {
//      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
//
//      const double time=time_pt()->time();
//
//      std::pair<double,double>u_xy_pair = RAYRAY::get_radial_v(nod_pt,time);
//
//      nod_pt->set_value(0,u_xy_pair.first);
//      nod_pt->set_value(1,0.0);
//      nod_pt->set_value(2,0.0);
//     }
//   }

  } // end of actions_before_implicit_timestep



////////////////////////////////


 void actions_before_distribute()
 {
   if(NSPP::Distribute_problem)
   {
     GenericProblemSetup::delete_flux_elements(Surface_mesh_pt);
     rebuild_global_mesh();
   }
 }

 void actions_after_distribute()
 {
   if(NSPP::Distribute_problem)
   {
     create_parall_outflow_lagrange_elements(Outflow_boundary,
         Bulk_mesh_pt,
         Surface_mesh_pt);
     rebuild_global_mesh();
   }
 }


 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   
   // Now set the first pressure dof in the first element to 0.0
//   if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt
 

 /// Update the after solve (empty)
 void actions_after_newton_solve()
 {
 }

 /// Before a newton solve, we set up a new "time step" in the
 /// Doc_linear_solver_info_pt.
 void actions_before_newton_solve()
 {
 } // end_of_actions_before_newton_solve


 // After a Newton step, we push in the iteration counts and timing
 // results.
 void actions_after_newton_step()
 {
   if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
   {
     NSPP::doc_iter_times(this,Doc_linear_solver_info_pt);
   }
 }

 /// Global error norm for adaptive time-stepping
 double global_temporal_error_norm();

 void create_parall_outflow_lagrange_elements(const unsigned &b,
                                              Mesh* const &bulk_mesh_pt,
                                              Mesh* const &surface_mesh_pt);
 /// Pointer to the "bulk" mesh
 Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

//private:

 /// Solver
 IterativeLinearSolver* Solver_pt;

 /// Solver
 Preconditioner* Prec_pt;

 /// Inexact solver for P block
 Preconditioner* P_matrix_preconditioner_pt;

 /// Inexact solver for F block
 Preconditioner* F_matrix_preconditioner_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;



 unsigned Small_curve_boundary;
 unsigned Big_curve_boundary;
 unsigned Z_max_boundary;
 unsigned Z_min_boundary;
 unsigned XZ_boundary;
 unsigned YZ_boundary;

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
CubeProblem<ELEMENT>::CubeProblem()
{ 
  // Assign boundary IDs defined in rayheader.h
  Small_curve_boundary = AW3DL::Small_curve_boundary;
  Big_curve_boundary = AW3DL::Big_curve_boundary;
  Z_max_boundary = AW3DL::Z_max_boundary;
  Z_min_boundary = AW3DL::Z_min_boundary;
  XZ_boundary = AW3DL::XZ_boundary;
  YZ_boundary = AW3DL::YZ_boundary; 

  // Identify the inflow and outflow boundaries.
  Inflow_boundary=Small_curve_boundary;
  Outflow_boundary=Big_curve_boundary;

  Doc_linear_solver_info_pt = NSPP::Doc_linear_solver_info_pt;
  if(NSPP::Delta_t < 0.0)
  {
    add_time_stepper_pt(new BDF<2>(true));
  }
  else
  {
    add_time_stepper_pt(new BDF<2>);
  }
  // Setup mesh


  Bulk_mesh_pt = 
    new AnnularWedgeCubicMesh<ELEMENT>(AW3DL::Noel,AW3DL::Noel,AW3DL::Noel,
        AW3DL::Lx, AW3DL::Ly, AW3DL::Lz,
        time_stepper_pt());


  // Create "surface mesh" that will contain only the prescribed-traction 
  // elements.
  Surface_mesh_pt = new Mesh;

  create_parall_outflow_lagrange_elements(Outflow_boundary,
      Bulk_mesh_pt,
      Surface_mesh_pt);

  // Add the two sub meshes to the problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);

  // Combine all submeshes into a single Mesh
  build_global_mesh();

  // Set up the boundary conditions. The order in which we set up the
  // boundaries is important. Since we do not want to unpin nodes which
  // previous boundary conditions has pinned. Pinned nodes always wins
  // over unpinned nodes.
  //
  // So we do the unpinned nodes first. Then work our way to the pinned nodes.

  // Outflow boundary. We unpin all the velocity components which is not
  // on another boundary. I.e., in 2D as want to just pin the "inner" nodes:
  //
  // -----------pinned!
  // |         |unpinned
  // |         |unpinned
  // |         |unpinned 
  // |         |unpinned
  // -----------pinned!
  //
  // Extend the above to 3D surface!
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
          nod_pt->unpin(0);
          nod_pt->unpin(1);
          nod_pt->unpin(2);
        }
      }
    }
  }


  // Now we do the top and bottom boundaries of the SLICE.
  // This corresponds to the front and back of the unit cube.
  // With these boundaries, only the uz velocity component is pinned to 0.
  // we let the x and y free to do what they do.
  {
    unsigned ibound=Z_max_boundary;
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      nod_pt->unpin(0);
      nod_pt->unpin(1);
      nod_pt->pin(2);

      // Set the values for safe measure.
      nod_pt->set_value(2,0.0);
    }
  }
  {
    unsigned ibound=Z_min_boundary;
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      nod_pt->unpin(0);
      nod_pt->unpin(1);
      nod_pt->pin(2);

      // Set the values for safe measure.
      nod_pt->set_value(2,0.0);
    }
  }


  // We now do the two sides where the knife cuts.
  {
    // The top boundary is parallel to the Y axis. As such, only uy is free.
    // The other velocity components should be 0.
    unsigned ibound=YZ_boundary;
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      nod_pt->pin(0);
      nod_pt->unpin(1);
      nod_pt->pin(2);

      // Set the values for safe measure.
      nod_pt->set_value(0,0.0);
      nod_pt->set_value(2,0.0);
    }
  }
  {
    // The bottom boundary of the unit cube. This is parallel to the x axis
    // in the annular wedge. As such, only ux is free, the other components
    // are pinned to zero.
    unsigned ibound=XZ_boundary;
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      nod_pt->unpin(0);
      nod_pt->pin(1);
      nod_pt->pin(2);

      // Set the values for safe measure.
      nod_pt->set_value(1,0.0);
      nod_pt->set_value(2,0.0);
    }
  }

  // Now do the inflow velocity. All of this is pinned.
  {
    // Now do the sides
    unsigned ibound=Inflow_boundary;
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
    {
      Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      nod_pt->pin(0);
      nod_pt->pin(1);
      nod_pt->pin(2);

      // Set the values for safe measure.
      nod_pt->set_value(0,0.0);
      nod_pt->set_value(1,0.0);
      nod_pt->set_value(2,0.0);
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
    el_pt->re_pt() = &NSPP::Rey;
    el_pt->re_st_pt() = &NSPP::Rey;
  } // end loop over elements


  // Setup equation numbering scheme
  oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

  if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
  {
    Vector<Mesh*> mesh_pt;
    if(NSPP::Prob_id == AW3DL::PID_AW3D_PO)
    {
      mesh_pt.resize(2,0);
      mesh_pt[0] = Bulk_mesh_pt;
      mesh_pt[1] = Surface_mesh_pt;
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
} // end_of_constructor


  template<class ELEMENT>
double CubeProblem<ELEMENT>::global_temporal_error_norm()
{
  return GenericProblemSetup::global_temporal_error_norm(this,
      3,Bulk_mesh_pt);
} // end of global_temporal_error_norm


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


      }
    } // Encapsulation

   } // Encapsulation

  } // Loop over all elements on boundary b

} // end of create_parall_outflow_lagrange_elements



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

  std::string label = AW3DL::prob_str()
    + NSPP::create_label() 
    + LPH::create_label() 
    + AW3DL::noel_str();
  return label;
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

  // Problem dimension.
  const unsigned dim = 3;

  // Set up doc info - used to store information on solver and iteration time.
  DocLinearSolverInfo doc_linear_solver_info;
  // Again, pass this to the NSPP and LPH
  NSPP::Doc_linear_solver_info_pt = &doc_linear_solver_info;
  LPH::Doc_linear_solver_info_pt = &doc_linear_solver_info;

  // Set the Label_pt
  LPH::Label_str_pt = &NSPP::Label_str;
  LPH::Vis_pt = &NSPP::Vis;
  AW3DL::Prob_id_pt = &NSPP::Prob_id;

  // Store commandline arguments
  CommandLineArgs::setup(argc,argv);

  NSPP::setup_commandline_flags();
  LPH::setup_commandline_flags();
  AW3DL::setup_commandline_flags(); 

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////

  // dim = 3
  NSPP::generic_problem_setup(dim);
  LPH::generic_setup();
  AW3DL::generic_setup(); 

  //////////////////////////////////////  

  {
  // Build the problem 
  CubeProblem <QTaylorHoodElement<3> >problem;


  // Solve the problem 
  //              problem.newton_solve();

  if(NSPP::Distribute_problem)
  {
    problem.distribute();
  }


  NSPP::Label_str = create_label();

  time_t rawtime;
  time(&rawtime);

  std::cout << "RAYDOING: "
    << NSPP::Label_str
    << " on " << ctime(&rawtime) << std::endl;

  GenericProblemSetup::doc_solution(problem.bulk_mesh_pt(),0);

  GenericProblemSetup::unsteady_run(&problem,
                                    &doc_linear_solver_info,
                                    problem.bulk_mesh_pt());
  }

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
    
    ResultsFormat::format_rayavgits(&iters_times,&results_stream);
    ResultsFormat::format_rayavavgits(&iters_times,&results_stream);
    
    // Now doing the preconditioner setup time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat::format_prectime(intimestep,&iters_times,&results_stream);
    }

    ResultsFormat::format_avgprectime(&iters_times,&results_stream);
    ResultsFormat::format_avavgprectime(&iters_times,&results_stream);

    // Now doing the linear solver time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat::format_solvertime(intimestep,&iters_times,&results_stream);
    }
    
    ResultsFormat::format_avgsolvertime(&iters_times,&results_stream);
    ResultsFormat::format_avavgsolvertime(&iters_times,&results_stream);
    
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

#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

} // end_of_main


