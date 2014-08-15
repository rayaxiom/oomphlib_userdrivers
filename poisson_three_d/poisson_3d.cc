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

// The Poisson equations
#include "poisson.h"

#include "meshes/simple_cubic_mesh.h"

#include "./../rayheader.h"


namespace NSPP = NavierStokesProblemParameters;

using namespace std;

using namespace oomph;


//============ start_of_namespace=====================================
/// Namespace for const source term in Poisson equation
//====================================================================
namespace ConstSourceForPoisson
{ 
 
 /// Strength of source function: default value -1.0
 const double Strength=-1.00;

/// Const source function
 void source_function(const Vector<double>& x, double& source)
 {
  source = Strength;
 }

} // end of namespace


namespace GlobalParam
{
  unsigned Noel = 0;

  const unsigned Length = 1;

  // AMG parameters:
  int AMG_iterations = -1;
  int AMG_smoother_iterations = -1;
  int AMG_simple_smoother = -1;
  int AMG_complex_smoother = -1;
  double AMG_damping = -1.0;
  double AMG_strength = -1.0;
  int AMG_coarsening = -1.0;

  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;
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
    delete Solver_pt;
    delete Prec_pt;
    delete Bulk_mesh_pt;
  }


//////////////////////////////

 void actions_before_implicit_timestep()
  {
  } // end of actions_before_implicit_timestep



////////////////////////////////


 void actions_before_distribute()
 {
 }

 void actions_after_distribute()
 {
 }


 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
  } // end_of_actions_after_adapt
 

 /// Update the after solve (empty)
 void actions_after_newton_solve()
 {
 }

 /// Update the problem specs before solve. 
 void actions_before_newton_solve()
 {
   GlobalParam::Doc_linear_solver_info_pt->setup_new_time_step();

 } // end_of_actions_before_newton_solve

 void actions_after_newton_step()
 {
   NSPP::doc_iter_times(this,GlobalParam::Doc_linear_solver_info_pt);
 }

 
 /// Pointer to the "bulk" mesh
 Mesh*& bulk_mesh_pt() {return Bulk_mesh_pt;}

private:

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
  Left_boundary = 4;
  Right_boundary = 2;
  Front_boundary = 5;
  Back_boundary = 0;
  Bottom_boundary = 1;
  Top_boundary = 3;


  const unsigned length = GlobalParam::Length;
  const unsigned noel = GlobalParam::Noel;

  Bulk_mesh_pt = 
      new SimpleCubicMesh<ELEMENT>(noel,noel,noel,length,length,length);

  
  add_sub_mesh(Bulk_mesh_pt);

  // Combine all submeshes into a single Mesh
  build_global_mesh();


 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here. Since the boundary values are never changed, we set
 // them here rather than in actions_before_solve(). 
 unsigned n_bound = Bulk_mesh_pt->nboundary();
 for(unsigned i=0;i<n_bound;i++)
  {
   unsigned n_node = mesh_pt()->nboundary_node(i);
   for (unsigned n=0;n<n_node;n++)
    {
     // Pin the single scalar value at this node
     Bulk_mesh_pt->boundary_node_pt(i,n)->pin(0); 

     // Assign the homogenous boundary condition for the one and only
     // nodal value
     Bulk_mesh_pt->boundary_node_pt(i,n)->set_value(0,0.0);
    }
  }

 // Loop over elements and set pointers to source function
 unsigned n_element = Bulk_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from FiniteElement to the present element
   ELEMENT *el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the source function pointer
   el_pt->source_fct_pt() = &ConstSourceForPoisson::source_function;
  }

 // Setup the equation numbering scheme
 std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 


 // Set the linear solver and preconditioner.
 
  
  // Create a new hypre preconditioner
  Prec_pt = new HyprePreconditioner;

  // Pointless cast because I want to.
  HyprePreconditioner* hypre_preconditioner_pt = 
      checked_static_cast<HyprePreconditioner*>(Prec_pt);

  // Set up the hypre preconditioner.

    // Set the hypre_method to BoomerAMG. This is hard coded.
    hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

    hypre_preconditioner_pt->set_amg_iterations(GlobalParam::AMG_iterations);

    hypre_preconditioner_pt->amg_smoother_iterations() 
      = GlobalParam::AMG_smoother_iterations;

    if(GlobalParam::AMG_simple_smoother >= 0)
    {
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother()
        = GlobalParam::AMG_simple_smoother;
    }
    else if(GlobalParam::AMG_complex_smoother >=0)
    {
      hypre_preconditioner_pt->amg_using_complex_smoothing();
      hypre_preconditioner_pt->amg_complex_smoother()
        = GlobalParam::AMG_complex_smoother;
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
    hypre_preconditioner_pt->amg_damping() = GlobalParam::AMG_damping;

    // Now set the AMG strength parameter.
    hypre_preconditioner_pt->amg_strength() = GlobalParam::AMG_strength;

    // AMG coarsening strategy.
    hypre_preconditioner_pt->amg_coarsening() = GlobalParam::AMG_coarsening;

    Hypre_Subsidiary_Preconditioner_Helper::print_hypre_settings(
            hypre_preconditioner_pt);



    // Now set up the solver.
  TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
  trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
  Solver_pt = trilinos_solver_pt;

  Solver_pt->tolerance() = 1e-8;
  Solver_pt->max_iter() = 100;
  Solver_pt->preconditioner_pt() = Prec_pt;

  this->linear_solver_pt() = Solver_pt;

  this->newton_solver_tolerance() = 1e-8;

} // end_of_constructor


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

  GlobalParam::Doc_linear_solver_info_pt = &doc_linear_solver_info;

  // Store commandline arguments
  CommandLineArgs::setup(argc,argv);

  CommandLineArgs::specify_command_line_flag("--noel", 
                                             &GlobalParam::Noel);
  CommandLineArgs::specify_command_line_flag("--amg_iter", 
                                             &GlobalParam::AMG_iterations);
  CommandLineArgs::specify_command_line_flag("--amg_smiter", 
                                             &GlobalParam::AMG_smoother_iterations);
  CommandLineArgs::specify_command_line_flag("--amg_sim_smoo", 
      &GlobalParam::AMG_simple_smoother);
  CommandLineArgs::specify_command_line_flag("--amg_com_smoo", 
      &GlobalParam::AMG_complex_smoother);
  CommandLineArgs::specify_command_line_flag("--amg_damp", 
      &GlobalParam::AMG_damping);
  CommandLineArgs::specify_command_line_flag("--amg_strn", 
      &GlobalParam::AMG_strength);
  CommandLineArgs::specify_command_line_flag("--amg_coarse", 
      &GlobalParam::AMG_coarsening);

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  ////////////////////////////////////////////////////
  // Now set up the flags/parameters for the problem//
  ////////////////////////////////////////////////////


  CubeProblem<QPoissonElement<3,2> > problem;

  problem.distribute();
  problem.newton_solve();



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

//    // If we want to output to a file, we create the outfile.
//    if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir"))
//    {
//      output_to_file = true;
//      std::ostringstream filename_stream;
//      filename_stream << NSPP::Itstime_dir_str<<"/"
//        << NSPP::Label_str
//        <<"NP"<<nproc<<"R"<<my_rank;
//      outfile.open(filename_stream.str().c_str());
//    }

    // Stringstream to hold the results. We do not output the results
    // (timing/iteration counts) as we get it since it will interlace with the
    // other processors and becomes hard to read.
    std::ostringstream results_stream;

    // Get the 3D vector which holds the iteration counts and timing results.
    Vector<Vector<Vector<double> > > iters_times
      = GlobalParam::Doc_linear_solver_info_pt->iterations_and_times();

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


//  //////////////////////////////////////////////////////////////////////////
//  ////////////// Outputting results ////////////////////////////////////////
//  //////////////////////////////////////////////////////////////////////////
//
//  if(NSPP::Solver_type != NSPP::Solver_type_DIRECT_SOLVE)
//  {
//    // Get the global oomph-lib communicator 
//    const OomphCommunicator* const comm_pt = MPI_Helpers::communicator_pt();
//
//    // my rank and number of processors. 
//    // This is used later for putting the data.
//    const unsigned my_rank = comm_pt->my_rank();
//    const unsigned nproc = comm_pt->nproc();
//
//    // Variable to indicate if we want to output to a file or not.
//    bool output_to_file = false;
//
//    // The output file.
//    std::ofstream outfile;
//
//    // If we want to output to a file, we create the outfile.
//    if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir"))
//    {
//      output_to_file = true;
//      std::ostringstream filename_stream;
//      filename_stream << NSPP::Itstime_dir_str<<"/"
//        << NSPP::Label_str
//        <<"NP"<<nproc<<"R"<<my_rank;
//      outfile.open(filename_stream.str().c_str());
//    }
//
//    // Stringstream to hold the results. We do not output the results
//    // (timing/iteration counts) as we get it since it will interlace with the
//    // other processors and becomes hard to read.
//    std::ostringstream results_stream;
//
//    // Get the 3D vector which holds the iteration counts and timing results.
//    Vector<Vector<Vector<double> > > iters_times
//      = NSPP::Doc_linear_solver_info_pt->iterations_and_times();
//
//    // Since this is a steady state problem, there is only
//    // one "time step", thus it is essentially a 2D vector 
//    // (the outer-most vector is of size 1).
//
//    // Loop over the time steps and output the iterations, prec setup time and
//    // linear solver time.
//    unsigned ntimestep = iters_times.size();
//    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
//    {
//      ResultsFormat::format_rayits(intimestep,&iters_times,&results_stream);
//    }
//    
//    ResultsFormat::format_rayavgits(&iters_times,&results_stream);
//    
//    // Now doing the preconditioner setup time.
//    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
//    {
//      ResultsFormat::format_prectime(intimestep,&iters_times,&results_stream);
//    }
//
//    ResultsFormat::format_avgprectime(&iters_times,&results_stream);
//
//    // Now doing the linear solver time.
//    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
//    {
//      ResultsFormat::format_solvertime(intimestep,&iters_times,&results_stream);
//    }
//    
//    ResultsFormat::format_avgsolvertime(&iters_times,&results_stream);
//    
//    // Print the result to oomph_info one processor at a time...
//    // This still doesn't seem to always work, since there are other calls
//    // to oomph_info before this one...
//    for (unsigned proc_i = 0; proc_i < nproc; proc_i++) 
//    {
//      if(proc_i == my_rank)
//      {
//        oomph_info << "\n" 
//          << "========================================================\n"
//          << results_stream.str()
//          << "========================================================\n"
//          << "\n" << std::endl;
//      }
//      MPI_Barrier(MPI_COMM_WORLD);
//    }
//
//    if(output_to_file)
//    {
//      outfile << "\n" << results_stream.str();
//      outfile.close();
//    }
//  }



#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif

} // end_of_main


