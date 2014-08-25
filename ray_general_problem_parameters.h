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


#ifndef RAY_GENERAL_PROBLEM_HEADER
#define RAY_GENERAL_PROBLEM_HEADER

// Oomph-lib includes
#include "generic.h"

using namespace oomph;

namespace GeneralProblemHelpers
{

  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;

  const static int Time_type_STEADY = 0;
  const static int Time_type_ADAPT = 1;
  const static int Time_type_FIXED = 2;

  int Time_type = -1; // SET FROM CL

  const static int Solver_type_DIRECT_SOLVE = 0;
  const static int Solver_type_OOMPHLIB_GMRES = 1;
  const static int Solver_type_TRILINOS_GMRES = 2;

  int Solver_type = -1; // SET FROM CL


  // Usually from NSPP:
  bool Distribute_problem = true;


  int Max_solver_iteration = -1;

  double Delta_t = -1.0;
  double Time_start = -1.0;
  double Time_end = -1.0;

  bool Doc_soln_flag = false;
  std::string Soln_dir_str = "";

  bool Doc_time_flag = false;
  std::string Itstime_dir_str = "";

  inline void setup_commandline_flags()
  {
    // This is used within the driver code and within here.
    // In the driver code, we may want to set different boundary conditions
    // or set the boundary conditions in different places 
    // (i.e. actions_before_newton_solver for steady state, and 
    //       actions_before_implicit_time_step for time stepping)
    CommandLineArgs::specify_command_line_flag("--time_type",
        &Time_type);

    // Used to set the solver in setup_solver() below, and to detect if we
    // need to doc the times and iteration counts or not.
    CommandLineArgs::specify_command_line_flag("--solver_type",
        &Solver_type);


    CommandLineArgs::specify_command_line_flag("--dist_prob");


    CommandLineArgs::specify_command_line_flag("--max_solver_iter", 
        &Max_solver_iteration);

    CommandLineArgs::specify_command_line_flag("--dt", &Delta_t);
    CommandLineArgs::specify_command_line_flag("--time_start", &Time_start);
    CommandLineArgs::specify_command_line_flag("--time_end", &Time_end);


    // Flag to output the solution.
    CommandLineArgs::specify_command_line_flag("--doc_soln", 
        &Soln_dir_str);

    // Iteration count and times directory.
    CommandLineArgs::specify_command_line_flag("--itstimedir", 
        &Itstime_dir_str);
  }


  inline void generic_setup()
  {
    if(!CommandLineArgs::command_line_flag_has_been_set("--time_type"))
    {
      std::ostringstream err_msg;
      err_msg << "Please set --time_type:\n"
        << "0 - steady state\n"
        << "1 - adaptive\n"
        << "2 - fixed time step" << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    if(!CommandLineArgs::command_line_flag_has_been_set("--solver_type"))
    {
      std::ostringstream err_msg;
      err_msg << "Please set --solver_type:\n"
        << "0 - direct solve\n"
        << "1 - OOMPH-LIB GMRES\n"
        << "2 - Trilinos GMRES" << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    if(CommandLineArgs::command_line_flag_has_been_set("--dist_prob"))
    {
      Distribute_problem = true;
    }
    else
    {
      Distribute_problem = false;
    }

    /////////////////

    if(!CommandLineArgs::command_line_flag_has_been_set("--max_solver_iter"))
    {
      std::ostringstream err_msg;
      err_msg << "Please set --max_solver_iter" << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Not check if the timing parameters has been set properly.
    if(Time_type == Time_type_ADAPT)
    {
      // start and end time must be set.
      if(!CommandLineArgs::command_line_flag_has_been_set("--time_start"))
      {
        std::ostringstream err_msg;
        err_msg << "Doing adaptive time stepping but\n"
          << "--time_start is not specified." << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      if(!CommandLineArgs::command_line_flag_has_been_set("--time_end"))
      {
        std::ostringstream err_msg;
        err_msg << "Doing adaptive time stepping but\n"
          << "--time_end is not specified." << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
    }
    else if (Time_type == Time_type_FIXED)
    {
      // All three time parameters must be specified.
      if(!CommandLineArgs::command_line_flag_has_been_set("--time_start"))
      {
        std::ostringstream err_msg;
        err_msg << "Doing fixed time stepping but\n"
          << "--time_start is not specified." << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      if(!CommandLineArgs::command_line_flag_has_been_set("--time_end"))
      {
        std::ostringstream err_msg;
        err_msg << "Doing fixed time stepping but\n"
          << "--time_end is not specified." << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }


      if(!CommandLineArgs::command_line_flag_has_been_set("--dt"))
      {
        std::ostringstream err_msg;
        err_msg << "Doing fixed time stepping but\n"
          << "--dt is not specified." << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

    }

    // Document the solution? Default is false.
    Doc_soln_flag = false;
    if(CommandLineArgs::command_line_flag_has_been_set("--doc_soln"))
    {
      // The argument immediately after --doc_soln is put into NSPP::Soln_dir_str.
      // If this begins with "--", then no solution directory has been provided.
      std::size_t found = Soln_dir_str.find("--");

      // Check if they have set the solution directory.
      if(found != std::string::npos)
      {
        std::ostringstream err_msg;
        err_msg << "Please provide the doc_soln directory "
          << "after the argument --doc_soln.\n" 
          << "This must not start with \"--\"." << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        Doc_soln_flag = true;
      }
    }


    Doc_time_flag = false;
    // Store the iteration and timing results in a file?
    // The its and time are always outputted in cout, but maybe we would like
    // to output it to a file.
    if(CommandLineArgs::command_line_flag_has_been_set("--itstimedir"))
    {
      // The argument immediately after --itstimedir is put into 
      // NSPP::Itstime_dir_str.
      // If this begins with "--", then no solution directory has been provided.
      std::size_t found = Itstime_dir_str.find("--");

      // Check if they have set the solution directory.
      if(found != std::string::npos)
      {
        std::ostringstream err_msg;
        err_msg << "Please provide the itstimedir directory "
          << "after the argument --itstimedir.\n" 
          << "This must not start with \"--\"." << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        Doc_time_flag = true;
      }
    }


    if(Doc_linear_solver_info_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "Please set GPH::Doc_linear_solver_info_pt from the\n"
        << "main function.\n"  << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  }


  // Pushes the iteration count and timing results in to the 
  // DocLinearSolverInfo object.
  inline void doc_iter_times(Problem* problem_pt,
      DocLinearSolverInfo* doc_linear_solver_info_pt)
  {
    unsigned iters = 0;
    double preconditioner_setup_time = 0.0;
    double solver_time = 0.0;

    // Get the iteration counts and preconditioner setup time
#ifdef PARANOID
    IterativeLinearSolver* iterative_solver_pt
      = dynamic_cast<IterativeLinearSolver*>
      (problem_pt->linear_solver_pt());
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
      (problem_pt->linear_solver_pt())->iterations();
    preconditioner_setup_time = static_cast<IterativeLinearSolver*>
      (problem_pt->linear_solver_pt())->preconditioner_pt()->setup_time();
#endif

    // Set the solver time.
    if(Solver_type == Solver_type_TRILINOS_GMRES)
    {
      TrilinosAztecOOSolver* trilinos_solver_pt 
        = dynamic_cast<TrilinosAztecOOSolver*>(problem_pt->linear_solver_pt());
      solver_time = trilinos_solver_pt->linear_solver_solution_time();
    }
    else
    {
      solver_time 
        = problem_pt->linear_solver_pt()->linear_solver_solution_time();
    }

    doc_linear_solver_info_pt->add_iteration_and_time
      (iters,preconditioner_setup_time,solver_time);
  }

  // A more generalised doc_solution function. Note: This is not 
  // parallelised yet... (maybe?)
  inline void doc_solution(Mesh* bulk_mesh_pt,
      const std::string& soln_dir,
      const std::string& label,
      const int& nt = -1)
  {
    std::ofstream some_file;
    std::stringstream filename;
    if(nt < 0)
    {
      filename << soln_dir<<"/"<< label <<".dat";
    }
    else
    {
      filename << soln_dir<<"/"<< label <<"t"<<nt<<".dat";
    }

    // Number of plot points
    const unsigned npts=5;

    // Output solution
    some_file.open(filename.str().c_str());
    bulk_mesh_pt->output(some_file,npts);
    some_file.close();
  }


  // setup the linear solver. NOTE: The solver is created by this function.
  inline IterativeLinearSolver* setup_solver(const int& max_solver_iter,
      const double& solver_tol, const double& newton_tol,
      Problem* problem_pt, Preconditioner* prec_pt)
  { 
    IterativeLinearSolver* iterative_linear_solver_pt = 0;
#ifdef OOMPH_HAS_TRILINOS
    if(Solver_type == Solver_type_TRILINOS_GMRES)
    {
      TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
      trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
      iterative_linear_solver_pt = trilinos_solver_pt;
    }
    else if(Solver_type == Solver_type_OOMPHLIB_GMRES)
    {
      iterative_linear_solver_pt = new GMRES<CRDoubleMatrix>;
      // We use RHS preconditioning. Note that by default,
      // left hand preconditioning is used.
      static_cast<GMRES<CRDoubleMatrix>*>(iterative_linear_solver_pt)
        ->set_preconditioner_RHS();
    }
#else
    if(Solver_type == Solver_type_TRILINOS_GMRES)
    {
      std::ostringstream err_msg;
      err_msg << "You have set --solver_type 2\n"
        << "But OOMPH does not have trilinos" << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else if(Solver_type == Solver_type_OOMPHLIB_GMRES)
    {
      iterative_linear_solver_pt = new GMRES<CRDoubleMatrix>;
      // We use RHS preconditioning. Note that by default,
      // left hand preconditioning is used.
      static_cast<GMRES<CRDoubleMatrix>*>(iterative_linear_solver_pt)
        ->set_preconditioner_RHS();
    }
#endif

    if(max_solver_iter < 0)
    {
      std::ostringstream err_msg;
      err_msg << "Max solver iteration is " << max_solver_iter << ".\n"
        << "Something has gone wrong. Have you set the flag\n"
        << "--max_solver_iter ?" << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Now set everything!
    if(iterative_linear_solver_pt != 0)
    {
      iterative_linear_solver_pt->tolerance() = solver_tol;
      iterative_linear_solver_pt->max_iter() = max_solver_iter;
      iterative_linear_solver_pt->preconditioner_pt() = prec_pt;
      problem_pt->linear_solver_pt() = iterative_linear_solver_pt;
    }

    problem_pt->newton_solver_tolerance() = newton_tol;

    return iterative_linear_solver_pt;
  }


  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  //Delete flux element function
  inline void delete_flux_elements(Mesh* const &surface_mesh_pt)
  {
    // How many surface elements are there in the mesh?
    const unsigned n_element = surface_mesh_pt->nelement();

    // Loop over the surface elements
    for(unsigned e=0;e<n_element;e++)
    {
      // Kill surface elements
      delete surface_mesh_pt->element_pt(e);
    }

    // Wipe the mesh
    surface_mesh_pt->flush_element_and_node_storage();
  }

  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////
  // Used for time stepping.
  inline void unsteady_run(Problem* problem_pt,
      Mesh* mesh_pt,
      DocLinearSolverInfo* doc_linear_solver_info_pt,
      const std::string& label_str,
      const std::string& soln_dir_str = std::string())
  {
    double dt = 0.0;
    if(Time_type == Time_type_ADAPT)
    {
      dt = 1e-2;
    }
    else
    {
      dt = Delta_t;
    }

    // Initialise all history values for an impulsive start
    problem_pt->assign_initial_values_impulsive(dt);
    oomph_info << "RAYINFO: IC = Impulsive start" << std::endl;

    // Now do many time steps
    if (Time_type == Time_type_FIXED)
    {
      const unsigned nsteps = unsigned(std::ceil((Time_end 
              - Time_start) / dt));

      oomph_info << "Taking constant time steps of: " << dt << std::endl;
      oomph_info << "NTIMESTEP is: " << nsteps << std::endl;
    }

    unsigned current_time_step = 0;

    if(Doc_soln_flag)
    {
      doc_solution(mesh_pt,soln_dir_str,label_str,current_time_step);
    }

    const double time_tol = 1e-4;

    while(problem_pt->time_pt()->time() < Time_end)
    {
      oomph_info << "TIMESTEP: " << current_time_step << std::endl;

      // Setup storage for a new time step
      if(Solver_type != Solver_type_DIRECT_SOLVE)
      {
        // Initialise counters for each newton solve.
        doc_linear_solver_info_pt->setup_new_time_step();
      }

      if (Time_type == Time_type_ADAPT) 
      {
        oomph_info << "DELTA_T: " << dt << std::endl;

        // Calculate the next time step.
        dt = problem_pt->adaptive_unsteady_newton_solve(dt,time_tol);
      }
      else
      {
        // Take one fixed time step
        problem_pt->unsteady_newton_solve(dt);
      }

      oomph_info << "Time is now: " 
        << problem_pt->time_pt()->time() << std::endl;

      if(Doc_soln_flag)
      {
        doc_solution(mesh_pt,soln_dir_str,label_str,current_time_step);
      }
      current_time_step++;
    }
  } // EoF unsteady_run



}


//=============================================================================
/// Namespace to format results.
//=============================================================================
namespace ResultsFormat
{
  inline void format_rayits(
      const unsigned& intimestep,      
      const Vector<Vector<Vector<double> > >* iters_times_pt,
      std::ostringstream* results_stream_pt)
  {
    // New timestep:
    (*results_stream_pt) << "RAYITS:\t" << intimestep << "\t";

    // Loop through the Newton Steps
    unsigned nnewtonstep = (*iters_times_pt)[intimestep].size();
    unsigned sum_of_newtonstep_iters = 0;
    for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
        innewtonstep++)
    {
      sum_of_newtonstep_iters += (*iters_times_pt)[intimestep][innewtonstep][0];
      (*results_stream_pt) << (*iters_times_pt)[intimestep][innewtonstep][0] << " ";
    }
    double average_its = ((double)sum_of_newtonstep_iters)
      / ((double)nnewtonstep);

    // Print to one decimal place if the average is not an exact
    // integer. Otherwise we print normally.
    //      ((unsigned(average_its*10))%10)?
    //        (*results_stream_pt) << "\t"<< std::fixed << std::setprecision(1)
    //        << average_its << "(" << nnewtonstep << ")" << "\n":
    //        (*results_stream_pt) << "\t"<< average_its << "(" << nnewtonstep << ")" << "\n";

    std::streamsize tmp_precision = results_stream_pt->precision();

    (*results_stream_pt) << "\t" << std::fixed << std::setprecision(1)
      << average_its << "(" << nnewtonstep << ")" << "\n";


    (*results_stream_pt) << std::setprecision(tmp_precision);
  }

  inline void format_rayavavgits(
      const Vector<Vector<Vector<double> > >* iters_times_pt,
      std::ostringstream* results_stream_pt)
  {
    // Loop through all the time steps
    const unsigned ntimestep = iters_times_pt->size();

    unsigned total_nnewton_step = 0;

    unsigned total_its = 0;
    unsigned n_total_its = 0;

    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      // Loop through the Newton Steps
      unsigned nnewtonstep = (*iters_times_pt)[intimestep].size();
      total_nnewton_step += nnewtonstep;

      for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
          innewtonstep++)
      {
        total_its += (*iters_times_pt)[intimestep][innewtonstep][0];
        n_total_its++;
      }
    }

    double average_its = ((double)total_its)
      / ((double)n_total_its);

    double average_n_newton_step = ((double)total_nnewton_step)
      / ((double)ntimestep);

    (*results_stream_pt) << "RAYAVGAVGITS:\t";
    // Print to one decimal place if the average is not an exact
    // integer. Otherwise we print normally.
    //    ((unsigned(average_its*10))%10)?
    //      (*results_stream_pt) << "\t"<< std::fixed << std::setprecision(1)
    //      << average_its << "(" << n_total_its << ")" << "\n":
    //     (*results_stream_pt) << "\t"<< average_its << "(" << n_total_its << ")" << "\n";
    std::streamsize tmp_precision = results_stream_pt->precision();

    (*results_stream_pt) << "\t" << std::fixed << std::setprecision(1)
      << average_its << "(" << average_n_newton_step << ")"
      << "(" << ntimestep << ")\n";

    // reset the precision
    (*results_stream_pt) << std::setprecision(tmp_precision);
  }


  inline void format_rayavgits(
      const Vector<Vector<Vector<double> > >* iters_times_pt,
      std::ostringstream* results_stream_pt)
  {
    // Loop through all the time steps
    const unsigned ntimestep = iters_times_pt->size();
    unsigned total_its = 0;
    unsigned n_total_its = 0;

    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {

      // Loop through the Newton Steps
      unsigned nnewtonstep = (*iters_times_pt)[intimestep].size();
      for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
          innewtonstep++)
      {
        total_its += (*iters_times_pt)[intimestep][innewtonstep][0];
        n_total_its++;
      }
    }

    double average_its = ((double)total_its)
      / ((double)n_total_its);

    (*results_stream_pt) << "RAYAVGITS:\t";
    // Print to one decimal place if the average is not an exact
    // integer. Otherwise we print normally.
    std::streamsize tmp_precision = results_stream_pt->precision();
    (*results_stream_pt) << "\t" << std::fixed << std::setprecision(1)
      << average_its << "(" << n_total_its << ")" << "\n";

    // reset the precision
    (*results_stream_pt) << std::setprecision(tmp_precision);
  }

  inline void format_prectime(
      const unsigned& intimestep,      
      const Vector<Vector<Vector<double> > >* iters_times_pt,
      std::ostringstream* results_stream_pt)
  {
    // New timestep:
    (*results_stream_pt) << "RAYPRECSETUP:\t" << intimestep << "\t";
    // Loop through the Newtom Steps
    unsigned nnewtonstep = (*iters_times_pt)[intimestep].size();
    double sum_of_newtonstep_times = 0;
    for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
        innewtonstep++)
    {
      sum_of_newtonstep_times += (*iters_times_pt)[intimestep][innewtonstep][1];
      (*results_stream_pt) << (*iters_times_pt)[intimestep][innewtonstep][1] << " ";
    }
    double average_time = ((double)sum_of_newtonstep_times)
      / ((double)nnewtonstep);

    // Print to one decimal place if the average is not an exact
    // integer. Otherwise we print normally.
    (*results_stream_pt) << "\t"<< average_time << "(" << nnewtonstep << ")" << "\n";
  }



  inline void format_avavgprectime(
      const Vector<Vector<Vector<double> > >* iters_times_pt,
      std::ostringstream* results_stream_pt)
  {

    // Loop through all the time steps
    const unsigned ntimestep = iters_times_pt->size();

    unsigned total_nnewton_step = 0;

    double total_time = 0.0;
    unsigned n_total_time = 0;

    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      // Loop through the Newton Steps
      unsigned nnewtonstep = (*iters_times_pt)[intimestep].size();
      total_nnewton_step += nnewtonstep;

      for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
          innewtonstep++)
      {
        total_time += (*iters_times_pt)[intimestep][innewtonstep][1];
        n_total_time++;
      }
    }

    (*results_stream_pt) << "RAYAVGAVGPRECSETUP:\t";

    double average_time = ((double)total_time)
      / ((double)n_total_time);

    double average_n_newton_step = ((double)total_nnewton_step)
      / ((double)ntimestep);


    // Print to one decimal place if the average is not an exact
    // integer. Otherwise we print normally.
    (*results_stream_pt) << "\t"
      << average_time 
      << "(" <<  average_n_newton_step << ")" 
      << "(" << ntimestep << ")" << "\n";
  }

  inline void format_avgprectime(
      const Vector<Vector<Vector<double> > >* iters_times_pt,
      std::ostringstream* results_stream_pt)
  {

    // Loop through all the time steps
    const unsigned ntimestep = iters_times_pt->size();
    double total_time = 0.0;
    unsigned n_total_time = 0;

    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    {
      // New timestep:
      //      (*results_stream_pt) << "RAYAVGPRECSETUP:\t" << intimestep << "\t";
      // Loop through the Newtom Steps
      unsigned nnewtonstep = (*iters_times_pt)[intimestep].size();
      for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
          innewtonstep++)
      {
        total_time += (*iters_times_pt)[intimestep][innewtonstep][1];
        n_total_time++;
      }
    }
    (*results_stream_pt) << "RAYAVGPRECSETUP:\t";
    double average_time = ((double)total_time)
      / ((double)n_total_time);

    // Print to one decimal place if the average is not an exact
    // integer. Otherwise we print normally.
    (*results_stream_pt) << "\t"<< average_time << "(" << n_total_time << ")" << "\n";
  }


  inline void format_solvertime(
      const unsigned& intimestep,      
      const Vector<Vector<Vector<double> > >* iters_times_pt,
      std::ostringstream* results_stream_pt)
  {
    // New timestep:
    (*results_stream_pt) << "RAYLINSOLVER:\t" << intimestep << "\t";
    // Loop through the Newtom Steps
    unsigned nnewtonstep = (*iters_times_pt)[intimestep].size();
    double sum_of_newtonstep_times = 0;
    for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
        innewtonstep++)
    {
      sum_of_newtonstep_times += (*iters_times_pt)[intimestep][innewtonstep][2];
      (*results_stream_pt) << (*iters_times_pt)[intimestep][innewtonstep][2] << " ";
    }
    double average_time = ((double)sum_of_newtonstep_times)
      / ((double)nnewtonstep);

    // Print to one decimal place if the average is not an exact
    // integer. Otherwise we print normally.
    (*results_stream_pt) << "\t"<< average_time << "(" << nnewtonstep << ")" << "\n";
  }



  inline void format_avavgsolvertime(
      const Vector<Vector<Vector<double> > >* iters_times_pt,
      std::ostringstream* results_stream_pt)
  {
    const unsigned ntimestep = iters_times_pt->size();

    unsigned total_nnewton_step = 0;

    double total_time = 0.0;
    unsigned n_total_time = 0;

    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    { 
      // Loop through the Newtom Steps
      unsigned nnewtonstep = (*iters_times_pt)[intimestep].size();

      total_nnewton_step += nnewtonstep;

      for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
          innewtonstep++)
      {
        total_time += (*iters_times_pt)[intimestep][innewtonstep][2];
        n_total_time++;
      }
    }
    (*results_stream_pt) << "RAYAVGAVGLINSOLVER:\t";


    double average_time = ((double)total_time)
      / ((double)n_total_time);
    double average_n_newton_step = ((double)total_nnewton_step)
      / ((double)ntimestep);


    // Print to one decimal place if the average is not an exact
    // integer. Otherwise we print normally.
    (*results_stream_pt) << "\t"
      << average_time 
      << "(" << average_n_newton_step << ")"
      << "(" << ntimestep << ")" << "\n";
  }


  inline void format_avgsolvertime(
      const Vector<Vector<Vector<double> > >* iters_times_pt,
      std::ostringstream* results_stream_pt)
  {
    // New timestep:
    //      (*results_stream_pt) << "RAYLINSOLVER:\t" << intimestep << "\t";
    // Loop through all the time steps
    const unsigned ntimestep = iters_times_pt->size();
    double total_time = 0.0;
    unsigned n_total_time = 0;

    for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
    { 

      // Loop through the Newtom Steps
      unsigned nnewtonstep = (*iters_times_pt)[intimestep].size();
      for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
          innewtonstep++)
      {
        total_time += (*iters_times_pt)[intimestep][innewtonstep][2];
        n_total_time++;
      }
    }
    (*results_stream_pt) << "RAYAVGLINSOLVER:\t";
    double average_time = ((double)total_time)
      / ((double)n_total_time);

    // Print to one decimal place if the average is not an exact
    // integer. Otherwise we print normally.
    (*results_stream_pt) << "\t"<< average_time << "(" << n_total_time << ")" << "\n";
  }
} // namespace ResultsFormat



#endif



