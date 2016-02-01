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
//#include "./../rayheader.h"

#include "./../ray_preconditioner_creation.h"
#include "./../ray_navier_stokes_parameters.h"
#include "./../ray_general_problem_parameters.h"



using namespace std;

using namespace oomph;


namespace GenProbHelpers = GeneralProblemHelpers;
namespace PrecHelpers = PreconditionerHelpers;
namespace NSHelpers = NavierStokesHelpers;


namespace ProblemHelpers
{
  unsigned Noel = 0; //set via commandline
  int Prob_id = -1; // set via commandline
  double Ang_deg = 0.0; // set via commandline

  // set via ang_deg
  double Ang_rad = 0.0;

  const double Length = 1.0;

  bool Vanilla = false;



  inline void specify_command_line_flags()
  {
    CommandLineArgs::specify_command_line_flag("--prob_id", &Prob_id);
    CommandLineArgs::specify_command_line_flag("--ang", &Ang_deg);
    CommandLineArgs::specify_command_line_flag("--noel",&Noel);
    CommandLineArgs::specify_command_line_flag("--vanilla");
  }

  inline void setup_command_line_flags()
  {
    // Set the Radian.
    if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
    {
      Ang_rad  = Ang_deg * (MathematicalConstants::Pi / 180.0);
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please set --ang" << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Check the problem id
    if(!CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
    {
      std::ostringstream err_msg;
      err_msg << "Please set --prob_id" << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION); 
    }

    if(CommandLineArgs::command_line_flag_has_been_set("--vanilla"))
    {
      Vanilla = true;
    }
    else
    {
      Vanilla = false;
    }
  }


  inline std::string prob_str()
  {
    std::string prob_str = "";

    if(Prob_id == 0)
    {
      prob_str = "CuPo";
    }
    else
    {
      std::ostringstream error_message;
      error_message << "No other problem ids done." << std::endl;

      throw OomphLibError(error_message.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    return prob_str;
  }

  inline std::string ang_deg_str()
  {
    std::ostringstream strs;
    strs << "A" << Ang_deg;
    return strs.str();
  }

  inline std::string noel_str()
  {
    std::ostringstream strs;
    strs << "N" << Noel;
    return strs.str();
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


  inline double get_prescribed_inflow_for_quarter(const double& y, 
      const double& z)
  {
#ifdef PARANOID
    // Quick check that we are in the correct range of
    // y and z coordinates.
    if( !((y > 0.5)&&(z > 0.5)) )
    {
      std::ostringstream error_message;
      error_message << "Prescribed inflow: incorrect range of y, z.\n"
        << "y = " << y << ", z = " << z << std::endl;

      throw OomphLibError(error_message.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
#endif
    return (y-0.5)*(1.0-y)*(z-0.5)*(1.0-z);
  }

  inline double get_prescribed_inflow_for_quarter(const double& t,
      const double& y,
      const double& z)
  {
    // For the velocity profile in the x direction.
    // 1) First form the parabolic profile

    // Note: +0.51, so at time = 0, there is still come velocity
    const double ux_scaling = -cos(MathematicalConstants::Pi*t)/2.0 + 0.51;
    return get_prescribed_inflow_for_quarter(y,z) * ux_scaling;
  } 


void get_prescribed_inflow(const double& t,
                           const double& y,
                           const double& z,
                           double& ux)
{    
  const double time_end = 1.0;
     
  // For the velocity profile in the x direction.
  // 1) First form the parabolic profile
  ux = 0.0;
  if((y > 0.5)&&(z > 0.5))
  {  
    const double ux_scaling = t / time_end;
    ux = (y-0.5)*(1.0-y)*(z-0.5)*(1.0-z) * ux_scaling;
  }  
} 


}

namespace ProbHelpers = ProblemHelpers;




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
      ProbHelpers::rotate_forward(x,y,z,phix,phiy,phiz,x_new);

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
//    GenericProblemSetup::clean_up_solver_memory();
//    LPH::clean_up_memory();

//    delete Surface_mesh_pt;
//    delete Bulk_mesh_pt;
//   delete Prec_pt;
//   delete P_matrix_preconditioner_pt;
//   delete F_matrix_preconditioner_pt;
  }

// ///Fix pressure in element e at pressure dof pdof and set to pvalue
// void fix_pressure(const unsigned &e, const unsigned &pdof, 
//                   const double &pvalue)
//  {
//   //Cast to full element type and fix the pressure at that element
//   dynamic_cast<NavierStokesEquations<3>*>(Bulk_mesh_pt->element_pt(e))->
//                          fix_pressure(pdof,pvalue);
//  } // end of fix_pressure


//////////////////////////////

 void actions_before_implicit_timestep()
 {
   if(GenProbHelpers::Solver_type !=GenProbHelpers::Solver_type_DIRECT_SOLVE)
   {
   if(GenProbHelpers::Time_type != GenProbHelpers::Time_type_STEADY)
   {
     // NOTE: before an implicit time step, we clear the previous times.
     // (This has no effect if no times were added i.e.
     // we are starting the first implicit time step)
     //
     // The logic for this is that we want to keep only the time step we
     // use, not the rejected ones.
     //
     // In the unsteady solve loop, it is there at we add new storage for
     // a new time step.
     Doc_linear_solver_info_pt->clear_current_time_step();
   }
   }

     if(!ProbHelpers::Vanilla)
     {
     // Inflow in upper half of inflow boundary
     const unsigned ibound=Inflow_boundary; 
     const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
     {
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

       // Only set if it is on a single boundary

//       std::set<unsigned>* bnd_pt = 0;

//       nod_pt->get_boundaries_pt(bnd_pt);

//       if(bnd_pt != 0)
       {
//         if(bnd_pt->size() < 2)
         {
           const double x=nod_pt->x(0);
           const double y=nod_pt->x(1);
           const double z=nod_pt->x(2);


           // Locally cache the angle for convenience.
           double ang_rad = ProbHelpers::Ang_rad;

           // Rotate x y and z back
           Vector<double>x_new;
           ProbHelpers::rotate_backward(
               x,y,z,-ang_rad,-ang_rad,-ang_rad,x_new);

           if((x_new[1] > 0.5) && (x_new[2] > 0.5))
           {
             const double time=time_pt()->time();
             double ux = ProbHelpers::get_prescribed_inflow_for_quarter
               (time,x_new[1],x_new[2]);

             // Now rotate the velocity profile
             Vector<double>u_new;
             ProbHelpers::rotate_forward(
                 ux,0,0,ang_rad,ang_rad,ang_rad,u_new);

             nod_pt->set_value(0,u_new[0]);
             nod_pt->set_value(1,u_new[1]);
             nod_pt->set_value(2,u_new[2]);
           }
           else
           {
             nod_pt->set_value(0,0.0);
             nod_pt->set_value(1,0.0);
             nod_pt->set_value(2,0.0);
           }
         }
       }
     }
     }// if !Vanilla
     else
     {
     // Inflow in upper half of inflow boundary
     const unsigned ibound=Inflow_boundary; 
     const unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
     {
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);

       // Only set if it is on a single boundary

//       std::set<unsigned>* bnd_pt = 0;

//       nod_pt->get_boundaries_pt(bnd_pt);

//       if(bnd_pt != 0)
       {
//         if(bnd_pt->size() < 2)
         {
           const double x=nod_pt->x(0);
           const double y=nod_pt->x(1);
           const double z=nod_pt->x(2);


           // Locally cache the angle for convenience.
//           double ang_rad = ProbHelpers::Ang_rad;

           // Rotate x y and z back
//           Vector<double>x_new;
//           ProbHelpers::rotate_backward(
//               x,y,z,-ang_rad,-ang_rad,-ang_rad,x_new);

           if((y > 0.5) && (z > 0.5))
           {
             const double time=time_pt()->time();
             double ux = ProbHelpers::get_prescribed_inflow_for_quarter
               (time,y,z);

             // Now rotate the velocity profile
//             Vector<double>u_new;
//             ProbHelpers::rotate_forward(
//                 ux,0,0,ang_rad,ang_rad,ang_rad,u_new);

             nod_pt->set_value(0,ux);
           }
           else
           {
             nod_pt->set_value(0,0.0);
           }
      nod_pt->set_value(1,0.0);
      nod_pt->set_value(2,0.0);
         }
       }
     }

     } // else Vanilla
 } // end of actions_before_implicit_timestep



////////////////////////////////


 void actions_before_distribute()
 {
   if(!ProbHelpers::Vanilla)
   {
   GenProbHelpers::delete_flux_elements(Surface_mesh_pt);
   rebuild_global_mesh();
   }
 }

 void actions_after_distribute()
 {
   if(!ProbHelpers::Vanilla)
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
 void actions_after_newton_step()
 {
   if(GenProbHelpers::Solver_type 
       != GenProbHelpers::Solver_type_DIRECT_SOLVE)
   {
     GenProbHelpers::doc_iter_times(this,Doc_linear_solver_info_pt);
   }
 }
 
 /// Global error norm for adaptive time-stepping
 double global_temporal_error_norm();

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

//private:

 /// Solver
 IterativeLinearSolver* Solver_pt;

 /// Solver
 Preconditioner* Prec_pt;

 Preconditioner* NS_matrix_preconditioner_pt;

 /// Inexact solver for P block
 Preconditioner* P_matrix_preconditioner_pt;

 /// Inexact solver for F block
 Preconditioner* F_matrix_preconditioner_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;



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
CubeProblem<ELEMENT>::CubeProblem()
{ 
  Left_boundary = 4;
  Right_boundary = 2;
  Front_boundary = 5;
  Back_boundary = 0;
  Bottom_boundary = 1;
  Top_boundary = 3;

  Inflow_boundary=Left_boundary;
  Outflow_boundary=Right_boundary;
  Doc_linear_solver_info_pt = GenProbHelpers::Doc_linear_solver_info_pt;


  if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_STEADY)
  {
    oomph_info 
      << "RAYINFO: Doing steady state, no time stepper added." << std::endl;
  }
  else if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_ADAPT)
  {
    oomph_info 
      << "RAYINFO: Adding adaptive time stepper" << std::endl;
    add_time_stepper_pt(new BDF<2>(true));
  }
  else if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_FIXED)
  {
    oomph_info << "RAYINFO: Adding non-adaptive time stepper" << std::endl;
    add_time_stepper_pt(new BDF<2>);
  }
  else
  {
    std::ostringstream err_msg;
    err_msg << "Time stepper for Time_type: "
      << GenProbHelpers::Time_type << std::endl;

    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }


//  Bulk_mesh_pt = 
//    new SimpleCubicMesh<ELEMENT>(n_x,n_y,n_z,l_x,l_y,l_z,time_stepper_pt());

  const unsigned noel = ProbHelpers::Noel;
  const double length = ProbHelpers::Length;
  const double ang_rad = ProbHelpers::Ang_rad;

  if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_STEADY)
  {
    Bulk_mesh_pt = 
      new SlopingCubicMesh<ELEMENT>(noel,noel,noel,
          length, length, length,
          ang_rad, ang_rad, ang_rad);
  }
  else if((GenProbHelpers::Time_type == GenProbHelpers::Time_type_ADAPT)
      ||GenProbHelpers::Time_type == GenProbHelpers::Time_type_FIXED)
  {
    if(ProbHelpers::Vanilla)
    {
      Bulk_mesh_pt = new SimpleCubicMesh<ELEMENT>(noel,noel,noel,
          length,length,length,time_stepper_pt());
    }
    else
    {
    Bulk_mesh_pt = 
      new SlopingCubicMesh<ELEMENT>(noel,noel,noel, 
          length, length, length,
          ang_rad,ang_rad,ang_rad,
          time_stepper_pt());
    }
  }
  else
  {
    std::ostringstream err_msg;
    err_msg << "No mesh set up for Time_type: "
      << GenProbHelpers::Time_type << std::endl;

    throw OomphLibError(err_msg.str(),
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
  }




//  {
//    const unsigned n_bound = Bulk_mesh_pt->nboundary();
//    oomph_info << "n_bound = " << n_bound << std::endl;
//
//    for(unsigned ibound=0;ibound<n_bound;ibound++)
//    {
//      oomph_info << "Boundary no: " << ibound << std::endl; 
//
//      unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
//      for (unsigned inod=0;inod<num_nod;inod++)
//      {
//        // Loop over values (u, v and w velocities)
//        Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);
//
//        const double x = nod_pt->x(0);
//        const double y = nod_pt->x(1);
//        const double z = nod_pt->x(2);
//
//        oomph_info << x << " " << y << " " << z << std::endl;
//
//      }
//    } // end loop over boundaries!
//  }
  

//  Bulk_mesh_pt = 
//    new SimpleCubicMesh<ELEMENT>(n_x,n_y,n_z,l_x,l_y,l_z);


  if(!ProbHelpers::Vanilla)
  {
  // Create "surface mesh" that will contain only the prescribed-traction 
  // elements.
  Surface_mesh_pt = new Mesh;


//   create_inflow_traction_elements(Inflow_boundary,
//                                   Bulk_mesh_pt,Surface_mesh_pt);

  create_parall_outflow_lagrange_elements(Outflow_boundary,
                                          Bulk_mesh_pt,
                                          Surface_mesh_pt);
  }

  if(ProbHelpers::Vanilla)
  {
  add_sub_mesh(Bulk_mesh_pt);
  }
  else
  {
  // Add the two sub meshes to the problem
  add_sub_mesh(Bulk_mesh_pt);
  add_sub_mesh(Surface_mesh_pt);
  }

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


  if(!ProbHelpers::Vanilla)
  {
    unsigned ibound=Outflow_boundary;
    unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
    const double ang_rad = ProbHelpers::Ang_rad;
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
            
            ProbHelpers::rotate_backward(
                x,y,z,
                -ang_rad,-ang_rad,-ang_rad,
                x_old);

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
  else
  {
    unsigned ibound = Outflow_boundary;
    unsigned num_nod = Bulk_mesh_pt->nboundary_node(ibound);
    for(unsigned inod = 0; inod < num_nod; inod++)
    {
      Node* nod_pt = Bulk_mesh_pt->boundary_node_pt(ibound,inod);
      // Only free if node is ONLY on a single boundary
      std::set<unsigned>* bnd_pt=0;
      nod_pt->get_boundaries_pt(bnd_pt);
      if(bnd_pt != 0)
      {
        if(!(nod_pt->is_on_boundary(0)))
        {
          if ((nod_pt->x(1)<0.5)&&nod_pt->x(2)<0.5) nod_pt->unpin(0); 
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
    el_pt->re_pt() = &NSHelpers::Rey;
    el_pt->re_st_pt() = &NSHelpers::Rey;
//    el_pt->re_invfr_pt() = &NSHelpers::Rey;
  } // end loop over elements


  // Now set the first pressure value in element 0 to 0.0
  // if (Problem_id==Global_Variables::Driven_cavity) fix_pressure(0,0,0.0);

  // Setup equation numbering scheme
  oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 


  F_matrix_preconditioner_pt = 0;
  P_matrix_preconditioner_pt = 0;
  NS_matrix_preconditioner_pt = 0;
  if(GenProbHelpers::Solver_type != GenProbHelpers::Solver_type_DIRECT_SOLVE)
  {
  if(!ProbHelpers::Vanilla)
  {
  F_matrix_preconditioner_pt
    = PrecHelpers::create_f_p_amg_preconditioner(PrecHelpers::F_amg_param,0);
  P_matrix_preconditioner_pt
    = PrecHelpers::create_f_p_amg_preconditioner(PrecHelpers::P_amg_param,1);

  NS_matrix_preconditioner_pt = PrecHelpers::create_lsc_preconditioner(
      this,
      Bulk_mesh_pt,
      F_matrix_preconditioner_pt,
      P_matrix_preconditioner_pt);

  Vector<Mesh*> mesh_pt(2,0);
  mesh_pt[0] = Bulk_mesh_pt;
  mesh_pt[1] = Surface_mesh_pt;
  Prec_pt = PrecHelpers::create_lgr_precondiitoner(
      mesh_pt,
      PrecHelpers::W_solver,
      NS_matrix_preconditioner_pt);

  }
  else
  {
  F_matrix_preconditioner_pt
    = PrecHelpers::create_f_p_amg_preconditioner(PrecHelpers::F_amg_param,0);
  P_matrix_preconditioner_pt
    = PrecHelpers::create_f_p_amg_preconditioner(PrecHelpers::P_amg_param,1);

  Prec_pt = PrecHelpers::create_lsc_preconditioner(
      this,
      Bulk_mesh_pt,
      F_matrix_preconditioner_pt,
      P_matrix_preconditioner_pt);


  }
  }
  const double solver_tol = 1.0e-6;
  const double newton_tol = 1.0e-6;

  Solver_pt = GenProbHelpers::setup_solver(
      GenProbHelpers::Max_solver_iteration,
      solver_tol,newton_tol,
      this,Prec_pt);




//  if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
//  {
//    Vector<Mesh*> mesh_pt;
//    if(NSPP::Prob_id == CL::PID_CU_PO_QUARTER)
//    {
//      mesh_pt.resize(2,0);
//      mesh_pt[0] = Bulk_mesh_pt;
//      mesh_pt[1] = Surface_mesh_pt;
//    }
//
//    LPH::Mesh_pt = mesh_pt;
//    LPH::Problem_pt = this;
//    Prec_pt = LPH::get_preconditioner();
//  }
//
// const double solver_tol = 1.0e-6;
// const double newton_tol = 1.0e-6;
// GenericProblemSetup::setup_solver(NSPP::Max_solver_iteration,
//                                   solver_tol,newton_tol,
//                                   NSPP::Solver_type,this,Prec_pt);
} // end_of_constructor

  template<class ELEMENT>
double CubeProblem<ELEMENT>::global_temporal_error_norm()
{
  return NSHelpers::global_temporal_error_norm(
      this,
      NSHelpers::Dim,
      Bulk_mesh_pt);
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

        // Now do the rest
        const double x = nod_pt->x(0);
        const double y = nod_pt->x(1);
        const double z = nod_pt->x(2);

        const double ang_rad = ProbHelpers::Ang_rad;

        Vector<double>x_old;
        ProbHelpers::rotate_backward(
            x,y,z,
            -ang_rad,-ang_rad,-ang_rad,
            x_old);

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
  
//  std::string label = CL::prob_str()
//                      + NSPP::create_label() 
//                      + LPH::create_label() 
//                      + CL::ang_deg_str() + CL::noel_str();
std::string label = "";
  if(!ProbHelpers::Vanilla)
  {
  label = ProbHelpers::prob_str()
     + NSHelpers::create_label()
     + PrecHelpers::Lgr_prec_str
     + ProbHelpers::ang_deg_str()
     + ProbHelpers::noel_str();
  }
  else
  {
    label = "testingvan";
  }

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

  // Set up doc info - used to store information on solver and iteration time.
  DocLinearSolverInfo doc_linear_solver_info;

  // The TWO things:
  GenProbHelpers::Doc_linear_solver_info_pt = &doc_linear_solver_info;
  NSHelpers::Dim = 3;

  // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  // From GeneralProblemHelpers, we need to set:
  // --time_type, 
  //   0 = steady state, 
  //   1 = adaptive, 
  //   2 = fixed.
  // --solver_type, 
  //   0 = exact solver, 
  //   1 = OOMPH-LIB's GMRES, 
  //   2 = trilinos GMRES
  // -- dist_prob (no arguments)
  //
  // --max_solver_iter - set an integer
  //
  // --dt - only set if doing fixed time stepping.
  // --time_start - always set if not steady state
  // --time_end - always set if not steady state
  //
  // -- doc_soln [soln_dir]
  //
  // --itstimedir [results_dir]
  //
  GenProbHelpers::specify_command_line_flags();

  // From NavierStokesHelpers, set:
  // --visc 0 or 1
  // --rey a double
  // --rey_start - only set if looping through reynolds numbers
  // --rey_incre - same as above
  // --rey_end - same as above
  NSHelpers::specify_command_line_flags();

  // From PreconditionerHelpers, set:
  // --f_solver 0 (exact) or 1 (amg)
  // --p_solver 0 (exact) or 1 (amg)
  //
  // if f_solver or p_solver is 1, we NEED to additionally set (for f_solver)
  // --f_amg_iter (usually 1)
  //
  // --f_amg_smiter (usually 2)
  //
  // --f_amg_sim_smoo
  //   0 = Jacobi, IMPORTANT: set --f_amg_damp (to something like 1).
  //   1 = Gauss-Seidel, sequential
  //       (very slow in parallel!)
  //   2 = Gauss-Seidel, interior points in parallel, boundary sequential
  //       (slow in parallel!)
  //   3 = hybrid Gauss-Seidel or SOR, forward solve
  //   4 = hybrid Gauss-Seidel or SOR, backward solve
  //   6 = hybrid symmetric Gauss-Seidel or SSOR 
  //
  // OR set this (NOT BOTH):
  // --f_amg_com_smoo
  //    
  //   6 = Schwarz
  //   7 = Pilut
  //   8 = ParaSails
  //   9 = Euclid
  //
  // --f_amg_str - strength of dependence
  //
  // --f_amg_coarse:
  //    0 = CLJP (parallel coarsening using independent sets)
  //    1 = classical RS with no boundary treatment (not recommended
  //        in parallel)
  //    3 = modified RS with 3rd pass to add C points on the boundaries
  //    6 = Falgout (uses 1 then CLJP using interior coarse points as
  //        first independent set) THIS IS DEFAULT ON DOCUMENTATION
  //    8 = PMIS (parallel coarsening using independent sets - lower
  //        complexities than 0, maybe also slower convergence)
  //    10= HMIS (one pass RS on each processor then PMIS on interior
  //        coarse points as first independent set)
  //    11= One pass RS on each processor (not recommended)
  //
  // --print_f_hypre - to print the hypre parameters to confirm.
  //    The information is extracted from the preconditioner after it has
  //    been created, so it is good to always do this.
  //
  // REPEAT FOR p solver if it is also using AMG
  // THERE ARE 6 things to set!
  //
  PrecHelpers::specify_command_line_flags();



  // --prob_id
  // --ang
  // --noel
  ProbHelpers::specify_command_line_flags();

  // --noel number of elements in 1D
  // --prob_id - currently, the only prob id is 0
//  problem_specific_setup_commandline_flags();

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  GenProbHelpers::setup_command_line_flags();
  NSHelpers::setup_command_line_flags();
  PrecHelpers::setup_command_line_flags();
  ProbHelpers::setup_command_line_flags();


  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////
  std::string label = "";

  // Build the problem
  CubeProblem <QTaylorHoodElement<3> >problem;

  // Solve the problem 
  //              problem.newton_solve();

  if(GenProbHelpers::Distribute_problem)
  {
    const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();
    const unsigned nproc = comm_pt->nproc();

    if(nproc == 1)
    {
      oomph_info << "RAYINFO: only 1 core, "
                 << "not distributing problem." << std::endl;
    }
    else
    {
      oomph_info << "RAYINFO: I am distributing the problem" << std::endl;

      problem.distribute();
    }
  }

  GenProbHelpers::Distribute_problem = problem.distributed();
  oomph_info << "Problem.distributed() is " 
             << GenProbHelpers::Distribute_problem << std::endl; 

  label = create_label();

  time_t rawtime;
  time(&rawtime);

  std::cout << "RAYDOING: "
    << label
    << " on " << ctime(&rawtime) << std::endl;

  // There are two types of solves, one for steady state, another for
  // time stepping.
  if(GenProbHelpers::Time_type == GenProbHelpers::Time_type_STEADY)
  {
    problem.newton_solve();

    if(GenProbHelpers::Doc_soln_flag)
    {
      GenProbHelpers::doc_solution(problem.bulk_mesh_pt(),
          GenProbHelpers::Soln_dir_str,
          label);
    }
  }
  else
  {
    GenProbHelpers::unsteady_run(&problem,
        problem.bulk_mesh_pt(),
        &doc_linear_solver_info,
        label,
        GenProbHelpers::Soln_dir_str);
  }


  //////////////////////////////////////////////////////////////////////////
  ////////////// Outputting results ////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////

  if(GenProbHelpers::Solver_type != 
     GenProbHelpers::Solver_type_DIRECT_SOLVE)
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
    if(GenProbHelpers::Doc_time_flag)
    {
      output_to_file = true;
      std::ostringstream filename_stream;
      filename_stream << GenProbHelpers::Itstime_dir_str<<"/"
        << label
        <<"NP"<<nproc<<"R"<<my_rank;
      outfile.open(filename_stream.str().c_str());
    }

    // Stringstream to hold the results. We do not output the results
    // (timing/iteration counts) as we get it since it will interlace with the
    // other processors and becomes hard to read.
    std::ostringstream results_stream;

    // Get the 3D vector which holds the iteration counts and timing results.
    Vector<Vector<Vector<double> > > iters_times
      = GenProbHelpers::Doc_linear_solver_info_pt->iterations_and_times();

    // Since this is a steady state problem, there is only
    // one "time step", thus it is essentially a 2D vector 
    // (the outer-most vector is of size 1).

    // Loop over the time steps and output the iterations, prec setup time and
    // linear solver time.
    unsigned ntimestep = iters_times.size();
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat2::format_rayits(intimestep,&iters_times,&results_stream);
    }
    
    ResultsFormat2::format_rayavgits(&iters_times,&results_stream);
    ResultsFormat2::format_rayavavgits(&iters_times,&results_stream);

    // Now doing the preconditioner setup time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat2::format_prectime(intimestep,&iters_times,&results_stream);
    }

    ResultsFormat2::format_avgprectime(&iters_times,&results_stream);
    ResultsFormat2::format_avavgprectime(&iters_times,&results_stream);
    // Now doing the linear solver time.
    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      ResultsFormat2::format_solvertime(intimestep,&iters_times,&results_stream);
    }
    
    ResultsFormat2::format_avgsolvertime(&iters_times,&results_stream);
    ResultsFormat2::format_avavgsolvertime(&iters_times,&results_stream); 

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


