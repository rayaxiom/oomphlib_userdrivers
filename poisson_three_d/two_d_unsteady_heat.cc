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
//Driver for 2D unsteady heat problem 


#include <fenv.h>

//Generic includes
#include "generic.h"
#include "navier_stokes.h"

// The Poisson equations
#include "poisson.h"


#include "meshes/simple_cubic_mesh.h"

#include "./../ray_general_problem_parameters.h"
#include "./../ray_preconditioner_creation.h"


// The unsteady heat equations
#include "unsteady_heat.h"

// Mesh
#include "meshes/rectangular_quadmesh.h"

using namespace std;

using namespace oomph;

using namespace MathematicalConstants;

namespace GenProbHelpers = GeneralProblemHelpers;
namespace PrecHelpers = PreconditionerHelpers;





// Namespace extension
namespace oomph
{

 namespace RayPreconditionerCreationFunctions
 {


  // AMG parameters:
  int AMG_iterations = -1;
  int AMG_smoother_iterations = -1;
  int AMG_simple_smoother = -1;
  int AMG_complex_smoother = -1;
  double AMG_damping = -1.0;
  double AMG_strength = -1.0;
  int AMG_coarsening = -1.0;
  unsigned AMG_print_level = 0;

  /// \short Helper function to create a SuperLu preconditioner (for use as
  /// the default subsididary preconditioner creator in
  /// GeneralPurposeBlockPreconditioners).
  inline Preconditioner* create_hypre_preconditioner()
  {
   Preconditioner* prec_pt = new HyprePreconditioner;

  // Pointless cast because I want to.
  HyprePreconditioner* hypre_preconditioner_pt = 
      checked_static_cast<HyprePreconditioner*>(prec_pt);

     // Set the hypre_method to BoomerAMG. This is hard coded.
    hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

    hypre_preconditioner_pt->set_amg_iterations(AMG_iterations);
    
    hypre_preconditioner_pt->amg_smoother_iterations() 
      = AMG_smoother_iterations;

        if(AMG_simple_smoother >= 0)
    {
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother()
        = AMG_simple_smoother;
    }
    else if(AMG_complex_smoother >=0)
    {
      hypre_preconditioner_pt->amg_using_complex_smoothing();
      hypre_preconditioner_pt->amg_complex_smoother()
        = AMG_complex_smoother;
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "You have not supplied a valid smoother.\n";
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

     // Set the damping parameter.
    hypre_preconditioner_pt->amg_damping() = AMG_damping;

    // Now set the AMG strength parameter.
    hypre_preconditioner_pt->amg_strength() = AMG_strength;

    // AMG coarsening strategy.
    hypre_preconditioner_pt->amg_coarsening() =AMG_coarsening;

    hypre_preconditioner_pt->amg_print_level() = AMG_print_level;

    PrecHelpers::print_hypre_parameters(
            hypre_preconditioner_pt);
    return prec_pt;
  }

 }


  //==start_of_mylinearelasticityelement===============================
/// Wrapper to make quadratic linear elasticity element block
/// preconditionable 
//===================================================================
template<unsigned DIM, unsigned NNODE_1D>
class MyUnsteadyHeatElement : public virtual QUnsteadyHeatElement<DIM,NNODE_1D>
{
 
public: 

 /// \short The number of "DOF types" that degrees of freedom in this element
 /// are sub-divided into: heat, one dimension only.
 unsigned ndof_types() const
  {
   return 1;
  }
 
/// Create a list of pairs for all unknowns in this element,
/// so the first entry in each pair contains the global equation
/// number of the unknown, while the second one contains the number
/// of the "DOF type" that this unknown is associated with.
/// (Function can obviously only be called if the equation numbering
/// scheme has been set up.)
/// 
/// The dof type enumeration (in 3D) is as follows:
/// u = 0
/// 
 void get_dof_numbers_for_unknowns(
  std::list<std::pair<unsigned long,unsigned> >& dof_lookup_list) const
  {
   // number of nodes
   unsigned n_node = this->nnode();
   
   // temporary pair (used to store dof lookup prior to being added to list)
   std::pair<unsigned,unsigned> dof_lookup;
   
   // loop over the nodes
   for (unsigned j=0;j<n_node;j++)
    {
     //loop over displacement components
     for (unsigned i=0;i<1;i++)
      {
       // determine local eqn number
       int local_eqn_number = this->nodal_local_eqn(j,i);
       
       // ignore pinned values - far away degrees of freedom resulting 
       // from hanging nodes can be ignored since these are be dealt
       // with by the element containing their master nodes
       if (local_eqn_number >= 0)
        {
         // store dof lookup in temporary pair: Global equation number
         // is the first entry in pair
         dof_lookup.first = this->eqn_number(local_eqn_number);
         
         // set dof numbers: Dof number is the second entry in pair
         dof_lookup.second = i;
         
         // add to list
         dof_lookup_list.push_front(dof_lookup);
        }
      }
    }
  }

};


//=======================================================================
/// Face geometry for element is the same as that for the underlying
/// wrapped element
//=======================================================================
template<unsigned DIM,unsigned NNODE_1D>
class FaceGeometry<MyUnsteadyHeatElement<DIM,NNODE_1D> >
 : public virtual QElement<DIM-1,NNODE_1D> 
 {
 public:
  FaceGeometry() : QElement<DIM-1,NNODE_1D>() {}
 };




}






/////////////////////////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////// 


//============ start_of_namespace=====================================
/// Namespace for const source term in Poisson equation
//====================================================================
namespace ConstSourceForPoisson
{ 
 
 /// Strength of source function: default value -1.0
 const double Strength=-1.00;

/// Const source function
 void source_function(const double& time, const Vector<double>& x, double& source)
 {
  source = Strength;
 }

} // end of namespace


namespace Global_Parameters
{
  unsigned Noel = 0;

  const unsigned Length = 1;

  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;
}






//======start_of_ExactSolnForUnsteadyHeat=====================
/// Namespace for unforced exact solution for UnsteadyHeat equation 
//====================================================================
//namespace ExactSolnForUnsteadyHeat
//{
//
// /// Decay factor
// double K=10;
//
// /// Angle of bump
// double Phi=1.0; 
// 
// /// Exact solution as a Vector
// void get_exact_u(const double& time, const Vector<double>& x, 
//                  Vector<double>& u)
// {
//  double zeta=cos(Phi)*x[0]+sin(Phi)*x[1];
//  u[0]=exp(-K*time)*sin(zeta*sqrt(K));
// }
// 
// /// Exact solution as a scalar
// void get_exact_u(const double& time, const Vector<double>& x, double& u)
// {
//  double zeta=cos(Phi)*x[0]+sin(Phi)*x[1];
//  u=exp(-K*time)*sin(zeta*sqrt(K));
// }
// 
// /// Source function to make it an exact solution 
// void get_source(const double& time, const Vector<double>& x, double& source)
// {
//  source = 0.0;
// }
//
//} // end of ExactSolnForUnsteadyHeat

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//=====start_of_problem_class=========================================
/// UnsteadyHeat problem 
//====================================================================
template<class ELEMENT>
class UnsteadyHeatProblem : public Problem
{

public:

 /// Constructor
 UnsteadyHeatProblem();

 /// Destructor (empty)
 ~UnsteadyHeatProblem(){}

 /// Update the problem specs after solve (empty)
 void actions_after_newton_solve() {}
 // After a Newton step, we push in the iteration counts and timing
 // results.
 void actions_after_newton_step()
 {
     GenProbHelpers::doc_iter_times(this,Doc_linear_solver_info_pt);
 }

 /// \short Update the problem specs before solve (empty)
 void actions_before_newton_solve() 
 {
   GenProbHelpers::doc_iter_times(this,Global_Parameters::Doc_linear_solver_info_pt);
 }

 /// Update the problem specs after solve (empty)
 void actions_after_implicit_timestep() {}

 /// \short Update the problem specs before next timestep: 
 /// Set Dirchlet boundary conditions from exact solution.
 void actions_before_implicit_timestep();

 /// \short Set initial condition (incl previous timesteps) according
 /// to specified function. 
// void set_initial_condition();

 /// Doc the solution
 void doc_solution(DocInfo& doc_info, ofstream& trace_file);


 /// Pointer to the "bulk" mesh
 Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

 
private:

 /// Pointer to source function
 UnsteadyHeatEquations<2>::UnsteadyHeatSourceFctPt Source_fct_pt;
 
 /// Pointer to control node at which the solution is documented
 Node* Control_node_pt;


 /// Solver
 IterativeLinearSolver* Solver_pt;

 /// Solver
 Preconditioner* Prec_pt;

 DocLinearSolverInfo* Doc_linear_solver_info_pt;

 unsigned Left_boundary;
 unsigned Right_boundary;
 unsigned Front_boundary;
 unsigned Back_boundary;
 unsigned Bottom_boundary;
 unsigned Top_boundary;


 /// Pointer to the "bulk" mesh
 Mesh* Bulk_mesh_pt;


}; // end of problem class

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
//========start_of_constructor============================================
/// Constructor for UnsteadyHeat problem in square domain
//========================================================================
template<class ELEMENT>
UnsteadyHeatProblem<ELEMENT>::UnsteadyHeatProblem()
{ 
  Doc_linear_solver_info_pt = GenProbHelpers::Doc_linear_solver_info_pt;
 // Allocate the timestepper -- this constructs the Problem's 
 // time object with a sufficient amount of storage to store the
 // previous timsteps. 
 add_time_stepper_pt(new BDF<2>);

  // Assign boundary IDs defined in rayheader.h
  Left_boundary = 4;
  Right_boundary = 2;
  Front_boundary = 5;
  Back_boundary = 0;
  Bottom_boundary = 1;
  Top_boundary = 3;

  const unsigned length = Global_Parameters::Length;
  const unsigned noel = Global_Parameters::Noel;


  Bulk_mesh_pt = 
      new SimpleCubicMesh<ELEMENT>(noel,noel,noel,length,length,length,
                                   time_stepper_pt());

  add_sub_mesh(Bulk_mesh_pt);

  // Combine all submeshes into a single Mesh
  build_global_mesh();


 // Setup parameters for exact solution
 // -----------------------------------

 // Decay parameter
// ExactSolnForUnsteadyHeat::K=5.0;

 // Setup mesh
 //-----------

 // Number of elements in x and y directions
// unsigned nx=5;
// unsigned ny=5;

 // Lengths in x and y directions
// double lx=1.0;
// double ly=1.0;

 // Build mesh
// mesh_pt() = new RectangularQuadMesh<ELEMENT>(nx,ny,lx,ly,time_stepper_pt());

 // Choose a control node at which the solution is documented
 //----------------------------------------------------------
 // Total number of elements
// unsigned n_el=mesh_pt()->nelement();

 // Choose an element in the middle
// unsigned control_el=unsigned(n_el/2);

 // Choose its first node as the control node
// Control_node_pt=mesh_pt()->finite_element_pt(control_el)->node_pt(0);

// cout << "Recording trace of the solution at: " 
//      << Control_node_pt->x(0) << " "
//      << Control_node_pt->x(1) << std::endl;


 // Set the boundary conditions for this problem: 
 // ---------------------------------------------
 // All nodes are free by default -- just pin the ones that have 
 // Dirichlet conditions here. 
 unsigned n_bound = mesh_pt()->nboundary();
 for(unsigned b=0;b<n_bound;b++)
  {
   unsigned n_node = mesh_pt()->nboundary_node(b);
   for (unsigned n=0;n<n_node;n++)
    {
     Bulk_mesh_pt->boundary_node_pt(b,n)->pin(0); 

     // Assign the homogenous boundary condition for the one and only
     // nodal value
     Bulk_mesh_pt->boundary_node_pt(b,n)->set_value(0,0.0);
    }
  } // end of set boundary conditions


 // Complete the build of all elements so they are fully functional
 //----------------------------------------------------------------

 // Find number of elements in mesh
 unsigned n_element = Bulk_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));

   //Set the source function pointer
   el_pt->source_fct_pt() = &ConstSourceForPoisson::source_function;
  }

 // Do equation numbering
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 


 if(CommandLineArgs::command_line_flag_has_been_set("--use_bpf"))
 {
   Prec_pt = new ExactBlockPreconditioner<CRDoubleMatrix>;
   ExactBlockPreconditioner<CRDoubleMatrix>* exact_block_prec_pt = 
     checked_static_cast<ExactBlockPreconditioner<CRDoubleMatrix>* >(Prec_pt);

//   exact_block_prec_pt->set_nmesh(1);
   exact_block_prec_pt->push_back_mesh(Bulk_mesh_pt);
   exact_block_prec_pt
     ->set_subsidiary_preconditioner_function(
         &RayPreconditionerCreationFunctions::create_hypre_preconditioner);
 }
 else
 {
   Prec_pt 
     = RayPreconditionerCreationFunctions::create_hypre_preconditioner();
 }


    // Now set up the solver.
  TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
//  trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
  trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::CG;
  Solver_pt = trilinos_solver_pt;

  Solver_pt->tolerance() = 1e-8;
  Solver_pt->max_iter() = 100;
  Solver_pt->preconditioner_pt() = Prec_pt;

  this->linear_solver_pt() = Solver_pt;

  this->newton_solver_tolerance() = 1e-8;




} // end of constructor



//=========start of actions_before_implicit_timestep===============================
/// \short Actions before timestep: update the domain, then reset the 
/// boundary conditions for the current time.
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::actions_before_implicit_timestep()
{
 // Get current time
 double time=time_pt()->time();


 // Recall that we have this:
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
 
 // Let's heat up the bottom only
 unsigned bot_bound = 1;
 //Loop over the boundaries
// unsigned num_bound = mesh_pt()->nboundary();
// for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Loop over the nodes on boundary
//   unsigned num_nod=mesh_pt()->nboundary_node(ibound);
   unsigned num_nod=mesh_pt()->nboundary_node(bot_bound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     Node* nod_pt=mesh_pt()->boundary_node_pt(bot_bound,inod);
     double u;
//     Vector<double> x(2);
//     x[0]=nod_pt->x(0);
//     x[1]=nod_pt->x(1);
     unsigned x_cor = nod_pt->x(0);
     unsigned z_cor = nod_pt->x(2);
     u = (1.0 - x_cor)*(1-z_cor)*time;


     // Get current values of the boundary conditions from the
     // exact solution
//     ExactSolnForUnsteadyHeat::get_exact_u(time,x,u);
     nod_pt->set_value(0,u);
    }
  }
} // end of actions_before_implicit_timestep



////======================start_of_set_initial_condition====================
///// \short Set initial condition: Assign previous and current values
///// from exact solution.
////========================================================================
//template<class ELEMENT>
//void UnsteadyHeatProblem<ELEMENT>::set_initial_condition()
//{ 
// // Backup time in global Time object
// double backed_up_time=time_pt()->time();
//         
// // Past history needs to be established for t=time0-deltat, ...
// // Then provide current values (at t=time0) which will also form
// // the initial guess for the first solve at t=time0+deltat
// 
// // Vector of exact solution value
// Vector<double> soln(1);
// Vector<double> x(2);
//
// //Find number of nodes in mesh
// unsigned num_nod = mesh_pt()->nnode();
//
// // Set continuous times at previous timesteps:
// // How many previous timesteps does the timestepper use?
// int nprev_steps=time_stepper_pt()->nprev_values();
// Vector<double> prev_time(nprev_steps+1);
// for (int t=nprev_steps;t>=0;t--)
//  {
//   prev_time[t]=time_pt()->time(unsigned(t));
//  } 
//
// // Loop over current & previous timesteps
// for (int t=nprev_steps;t>=0;t--)
//  {
//   // Continuous time
//   double time=prev_time[t];
//   cout << "setting IC at time =" << time << std::endl;
//   
//   // Loop over the nodes to set initial guess everywhere
//   for (unsigned n=0;n<num_nod;n++)
//    {
//     // Get nodal coordinates
//     x[0]=mesh_pt()->node_pt(n)->x(0);
//     x[1]=mesh_pt()->node_pt(n)->x(1);
//
//     // Get exact solution at previous time
//     ExactSolnForUnsteadyHeat::get_exact_u(time,x,soln);
//     
//     // Assign solution
//     mesh_pt()->node_pt(n)->set_value(t,0,soln[0]);
//     
//     // Loop over coordinate directions: Mesh doesn't move, so 
//     // previous position = present position
//     for (unsigned i=0;i<2;i++)
//      {
//       mesh_pt()->node_pt(n)->x(t,i)=x[i];
//      }
//    } 
//  }
//
// // Reset backed up time for global timestepper
// time_pt()->time()=backed_up_time;
//
//} // end of set_initial_condition



//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnsteadyHeatProblem<ELEMENT>::
doc_solution(DocInfo& doc_info,ofstream& trace_file)
{ 
// ofstream some_file;
// char filename[100];
//
// // Number of plot points
// unsigned npts;
// npts=5;
//
//
// cout << std::endl;
// cout << "=================================================" << std::endl;
// cout << "Docing solution for t=" << time_pt()->time() << std::endl;
// cout << "=================================================" << std::endl;


 // Output solution 
 //-----------------
// sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
//         doc_info.number());
// some_file.open(filename);
// mesh_pt()->output(some_file,npts);

 /// Write file as a tecplot text object
// some_file << "TEXT X=2.5,Y=93.6,F=HELV,HU=POINT,C=BLUE,H=26,T=\"time = " 
//           << time_pt()->time() << "\"";
// // ...and draw a horizontal line whose length is proportional
// // to the elapsed time
// some_file << "GEOMETRY X=2.5,Y=98,T=LINE,C=BLUE,LT=0.4" << std::endl;
// some_file << "1" << std::endl;
// some_file << "2" << std::endl;
// some_file << " 0 0" << std::endl;
// some_file << time_pt()->time()*20.0 << " 0" << std::endl;
// some_file.close();


 // Output exact solution 
 //----------------------
// sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
//         doc_info.number());
// some_file.open(filename);
// mesh_pt()->output_fct(some_file,npts,time_pt()->time(),
//                       ExactSolnForUnsteadyHeat::get_exact_u); 
// some_file.close();

 // Doc error
 //----------
// double error,norm;
// sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
//         doc_info.number());
// some_file.open(filename);
// mesh_pt()->compute_error(some_file,
//                          ExactSolnForUnsteadyHeat::get_exact_u,
//                          time_pt()->time(),
//                          error,norm); 
// some_file.close();

 // Doc solution and error
 //-----------------------
// cout << "error: " << error << std::endl; 
// cout << "norm : " << norm << std::endl << std::endl;

 // Get exact solution at control node
// Vector<double> x_ctrl(2);
// x_ctrl[0]=Control_node_pt->x(0);
// x_ctrl[1]=Control_node_pt->x(1);
// double u_exact;
// ExactSolnForUnsteadyHeat::get_exact_u(time_pt()->time(),x_ctrl,u_exact);
// trace_file << time_pt()->time() << " " 
//            << Control_node_pt->value(0) << " " 
//            << u_exact << " " 
//            << error   << " " 
//            << norm    << std::endl;

} // end of doc_solution



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//=======start_of_main====================================================
/// \short Driver code for unsteady heat equation
//========================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc,argv);
#endif


  // Set up doc info - used to store information on solver and iteration time.
  DocLinearSolverInfo doc_linear_solver_info;

  Global_Parameters::Doc_linear_solver_info_pt = &doc_linear_solver_info;
GenProbHelpers::Doc_linear_solver_info_pt = &doc_linear_solver_info;
  // Store commandline arguments
  CommandLineArgs::setup(argc,argv);

  CommandLineArgs::specify_command_line_flag("--noel", 
    &Global_Parameters::Noel);
  CommandLineArgs::specify_command_line_flag("--amg_iter", 
    &RayPreconditionerCreationFunctions::AMG_iterations);
  CommandLineArgs::specify_command_line_flag("--amg_smiter", 
    &RayPreconditionerCreationFunctions::AMG_smoother_iterations);
  CommandLineArgs::specify_command_line_flag("--amg_sim_smoo", 
    &RayPreconditionerCreationFunctions::AMG_simple_smoother);
  CommandLineArgs::specify_command_line_flag("--amg_com_smoo", 
    &RayPreconditionerCreationFunctions::AMG_complex_smoother);
  CommandLineArgs::specify_command_line_flag("--amg_damp", 
      &RayPreconditionerCreationFunctions::AMG_damping);
  CommandLineArgs::specify_command_line_flag("--amg_strn", 
      &RayPreconditionerCreationFunctions::AMG_strength);
  CommandLineArgs::specify_command_line_flag("--amg_coarse", 
      &RayPreconditionerCreationFunctions::AMG_coarsening);
  CommandLineArgs::specify_command_line_flag("--amg_print_level", 
      &RayPreconditionerCreationFunctions::AMG_print_level);
  CommandLineArgs::specify_command_line_flag("--use_bpf");




  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  GenProbHelpers::Solver_type = 2;
  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////



 // Build problem
 UnsteadyHeatProblem<MyUnsteadyHeatElement<3,3> >
  problem;
 
 // Setup labels for output
// DocInfo doc_info;

 // Output directory
// doc_info.set_directory("RESLT");
 
 // Output number
// doc_info.number()=0;
 
 // Open a trace file
// ofstream trace_file;
// char filename[100];   
// sprintf(filename,"%s/trace.dat",doc_info.directory().c_str());
// trace_file.open(filename);
// trace_file << "VARIABLES=\"time\",\"u<SUB>FE</SUB>\","
//            << "\"u<SUB>exact</SUB>\",\"norm of error\",\"norm of solution\""
//            << std::endl;

 // Choose simulation interval and timestep
 double t_max=0.5;
 double dt=0.01;

 // Initialise timestep -- also sets the weights for all timesteppers
 // in the problem.
// problem.initialise_dt(dt);
 
 // Set IC
// problem.set_initial_condition();
 problem.assign_initial_values_impulsive(dt);

 //Output initial condition
// problem.doc_solution(doc_info,trace_file);
 
 //Increment counter for solutions 
// doc_info.number()++;

 // Find number of steps
 unsigned nstep = unsigned(t_max/dt);

 // Timestepping loop
 for (unsigned istep=0;istep<nstep;istep++)
  {
   oomph_info << "Timestep " << istep << std::endl;

   doc_linear_solver_info.setup_new_time_step();

   
   // Take timestep
   problem.unsteady_newton_solve(dt);


   oomph_info << "Time is now " << dt << std::endl;
   
   //Output solution
//   problem.doc_solution(doc_info,trace_file);
   
   //Increment counter for solutions 
//   doc_info.number()++;
  }


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

 }
 
 // Close trace file
// trace_file.close();
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

}; // end of main
