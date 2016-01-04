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


#ifndef RAY_PRECONDITIONER_HEADER
#define RAY_PRECONDITIONER_HEADER

// Oomph-lib includes
#include "generic.h"

using namespace oomph;


#ifdef OOMPH_HAS_HYPRE

//=============================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//=============================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
 Preconditioner* set_hypre_JhalfStrnSimOneVTwoTwoRS()
 {
  // Create a new HyprePreconditioner
  HyprePreconditioner* hypre_preconditioner_pt = new HyprePreconditioner;

  // Set the smoothers
  //   0 = Jacobi 
  //   1 = Gauss-Seidel, sequential
  //       (very slow in parallel!)
  //   2 = Gauss-Seidel, interior points in parallel, boundary sequential
  //       (slow in parallel!)
  //   3 = hybrid Gauss-Seidel or SOR, forward solve
  //   4 = hybrid Gauss-Seidel or SOR, backward solve
  //   6 = hybrid symmetric Gauss-Seidel or SSOR
  //
  //   As per Richard p91, we use Jacobi with a damping param if 0.5
  hypre_preconditioner_pt->amg_simple_smoother() = 0;
  // Set smoother damping
  hypre_preconditioner_pt->amg_damping() = 0.5;


  // Set the strength to 0.25
  hypre_preconditioner_pt->amg_strength() = 0.25;

  // Now set the coarsening
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
  hypre_preconditioner_pt->amg_coarsening() = 1;

  // Now create 1v22 cycles
  
  // Set the number of amg smoother iterations. This is usually set to 2
  hypre_preconditioner_pt->amg_smoother_iterations()=2;

  // Set the amg_iterations, this is usually set to 1, for V(2,2)
  // there is only one of them... I think that's what I means, 
  // the two is the smoother... pre and post smoothing.
  hypre_preconditioner_pt
      ->set_amg_iterations(1);

  return hypre_preconditioner_pt;
 }

 Preconditioner* set_hypre_JhalfStrnStrOneVTwoTwoRS()
 {
  // Create a new HyprePreconditioner
  HyprePreconditioner* hypre_preconditioner_pt = new HyprePreconditioner;

  // Set the smoothers
  //   0 = Jacobi 
  //   1 = Gauss-Seidel, sequential
  //       (very slow in parallel!)
  //   2 = Gauss-Seidel, interior points in parallel, boundary sequential
  //       (slow in parallel!)
  //   3 = hybrid Gauss-Seidel or SOR, forward solve
  //   4 = hybrid Gauss-Seidel or SOR, backward solve
  //   6 = hybrid symmetric Gauss-Seidel or SSOR
  //
  //   As per Richard p91, we use Jacobi with a damping param if 0.5
  hypre_preconditioner_pt->amg_simple_smoother() = 0;
  // Set smoother damping
  hypre_preconditioner_pt->amg_damping() = 0.5;


  // Set the strength to 0.25
  hypre_preconditioner_pt->amg_strength() = 0.668;

  // Now set the coarsening
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
  hypre_preconditioner_pt->amg_coarsening() = 1;

  // Now create 1v22 cycles
  
  // Set the number of amg smoother iterations. This is usually set to 2
  hypre_preconditioner_pt->amg_smoother_iterations()=2;

  // Set the amg_iterations, this is usually set to 1, for V(2,2)
  // there is only one of them... I think that's what I means, 
  // the two is the smoother... pre and post smoothing.
  hypre_preconditioner_pt
      ->set_amg_iterations(1);

  return hypre_preconditioner_pt;
 }


}

#endif



namespace PreconditionerHelpers
{


  const static int W_solver_exact = 0;
  const static int W_solver_cg = 1;
  int W_solver = -1;

  const static int NS_solver_exact = 0;
  const static int NS_solver_lsc = 1;
  int NS_solver = -1;

  const static int F_solver_exact = 0;
  const static int F_solver_amg = 1;
  const static int F_solver_blkdiag_amg = 9090;
  const static int F_solver_upper_amg = 9091;
  const static int F_solver_lower_amg = 9092;
  const static int F_solver_full_amg = 9093;
  int F_solver = -1;

  const static int P_solver_exact = 0;
  const static int P_solver_amg = 1;
  int P_solver = -1;

  double Scaling_sigma = 0.0;
  double Scaling_sigma_multiplier = 0.0;





  // Variations:
  // W solver: We (exact), Wcj (conjugate gradient)
  // NS solver: Ne (exact), Nl (lsc)
  // So we can have WeNl for example.
  std::string Lgr_prec_str="";

  // Variations:
  // F_solver: Fe (exact), Fa (amg),
  // P_solver: Pe (exact), Pa (amg)
  std::string NS_prec_str="";

  // These are non-null if and only if AMG as been set. Then it takes the
  // form:
  // [cycle][strn][coarsener][smoother]
  std::string NS_f_prec_str="";
  std::string NS_p_prec_str="";

  // These are the parameters which will be used when the 
  // create_hypre_preconditioner function is used.
  int G_AMG_iterations = -1;
  int G_AMG_smoother_iterations = -1;
  int G_AMG_simple_smoother = -1;
  int G_AMG_complex_smoother = -1;
  double G_AMG_damping = -1.0;
  double G_AMG_strength = -1.0;
  int G_AMG_coarsening = -1.0;

  // Not set as an AMG parameter, but used for book keeping.
  bool G_AMG_print_parameters = false;


  // These are not set as AMG parameters, these are additional stuff used
  // by me for book keeping.
  bool G_AMG_parameters_has_been_set = false;

  // A struct to store custom parameters.
  struct AMGParameters
  {
    int Iterations;
    int Smoother_iterations;
    int Simple_smoother;
    int Complex_smoother;
    double Damping;
    double Strength;
    int Coarsening;
    bool Print_parameters;


    AMGParameters()
    {
      // AMG parameters:
      Iterations = -1;
      Smoother_iterations = -1;
      Simple_smoother = -1;
      Complex_smoother = -1;
      Damping = -1.0;
      Strength = -1.0;
      Coarsening = -1.0;
      Print_parameters = false;
    }
  };

  AMGParameters F_amg_param;
  AMGParameters P_amg_param;

  // Reset all of the 
  inline void reset_amg_parameters()
  {
    G_AMG_iterations = -1;
    G_AMG_smoother_iterations = -1;
    G_AMG_simple_smoother = -1;
    G_AMG_complex_smoother = -1;
    G_AMG_damping = -1.0;
    G_AMG_strength = -1.0;
    G_AMG_coarsening = -1.0;
    G_AMG_print_parameters = false;
    G_AMG_parameters_has_been_set = false;
  }

  // reset_amg_parameters() must be called before this can be called.
  inline void set_amg_parameters(AMGParameters& amg_param)
  {
    if(G_AMG_parameters_has_been_set)
    {
      std::ostringstream err_msg;
      err_msg << "AMG parameters has not been reset. Please reset them\n";
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      G_AMG_iterations = amg_param.Iterations;
      G_AMG_smoother_iterations = amg_param.Smoother_iterations;
      G_AMG_simple_smoother = amg_param.Simple_smoother;
      G_AMG_complex_smoother = amg_param.Complex_smoother;
      G_AMG_damping = amg_param.Damping;
      G_AMG_strength = amg_param.Strength;
      G_AMG_coarsening = amg_param.Coarsening;
      G_AMG_print_parameters = amg_param.Print_parameters;
      G_AMG_parameters_has_been_set = true;
    }
  }


  inline void print_hypre_parameters(HyprePreconditioner* h_prec_pt)
  {
    // Print AMG iterations:
    oomph_info << "HyprePreconditioner settings are: " << std::endl;
    oomph_info << "Max_iter: " << h_prec_pt->amg_iterations() << std::endl;
    oomph_info << "smoother iter: " << h_prec_pt->amg_smoother_iterations() 
      << std::endl;
    oomph_info << "Hypre_method: " << h_prec_pt->hypre_method() << std::endl;
    oomph_info << "internal_preconditioner: " 
      << h_prec_pt->internal_preconditioner() << std::endl;
    oomph_info << "AMG_using_simple_smoothing: " 
      << h_prec_pt->amg_using_simple_smoothing_flag() << std::endl; 
    oomph_info << "AMG_simple_smoother: " 
      << h_prec_pt->amg_simple_smoother() << std::endl; 
    oomph_info << "AMG_complex_smoother: " 
      << h_prec_pt->amg_complex_smoother() << std::endl; 
    oomph_info << "AMG_coarsening: " 
      << h_prec_pt->amg_coarsening() << std::endl;
    oomph_info << "AMG_max_levels: " 
      << h_prec_pt->amg_max_levels() << std::endl;
    oomph_info << "AMG_damping: " 
      << h_prec_pt->amg_damping() << std::endl;
    oomph_info << "AMG_strength: " 
      << h_prec_pt->amg_strength() << std::endl;;
    oomph_info << "AMG_max_row_sum: " 
      << h_prec_pt->amg_max_row_sum() << std::endl;
    oomph_info << "AMG_truncation: " 
      << h_prec_pt->amg_truncation() << std::endl;
    oomph_info << "\n" << std::endl;

  }


  ///////////////////////////////////////////////////////////////////////////
  // Function to create a hypre preconditioner based on the global hypre
  // parameters, these begin with G_AMG_
  // 
  // After creating the preconditioner, the global parameters are reset
  // to the default, and would have to be set again.
  //
  // Forcing the user to call set_amg_parameters just before he calls 
  // create_hypre_preconditioner, we ensure that the user knows which
  // preconditioner he is 
  //
  ///////////////////////////////////////////////////////////////////////////
  inline HyprePreconditioner* create_hypre_preconditioner(AMGParameters &amg_param)
  {
    HyprePreconditioner* hypre_preconditioner_pt = 
      new HyprePreconditioner;

    // Set the hypre_method to BoomerAMG. This is hard coded.
    hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

    // Set the parameters above.
    hypre_preconditioner_pt->set_amg_iterations(amg_param.Iterations);

    hypre_preconditioner_pt->amg_smoother_iterations() 
      = amg_param.Smoother_iterations;

    // Check if both simple and complex smoother are greater than 0, only 
    // one should be greater than 0.
    const int simple_smoother = amg_param.Simple_smoother;
    const int complex_smoother = amg_param.Complex_smoother;

    if((simple_smoother >= 0) 
        && (complex_smoother >= 0))
    {
      // Both are >= 0, only 1 should be.
      std::ostringstream err_msg;
      err_msg << "Both simple and complex smoother is set.\n"
        << "simple_smoother: " << simple_smoother << "\n"
        << "complex_smoother: " << complex_smoother 
        << "\n";
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else if(simple_smoother >= 0)
    {
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother()
        = simple_smoother;
    }
    else if(complex_smoother >= 0)
    {
      hypre_preconditioner_pt->amg_using_complex_smoothing();
      hypre_preconditioner_pt->amg_complex_smoother()
        = complex_smoother;
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Both smoother parameters are less than 0.\n"
        << "Please set at most one of them." << "\n"
        << "simple_smoother: " << simple_smoother << "\n"
        << "complex_smoother: " << complex_smoother <<"\n";
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Set the damping parameter
    hypre_preconditioner_pt->amg_damping() = amg_param.Damping;

    // Set the strength of dependence.
    hypre_preconditioner_pt->amg_strength() = amg_param.Strength;

    // Set the coarsening strategy
    hypre_preconditioner_pt->amg_coarsening() = amg_param.Coarsening;

    if(amg_param.Print_parameters)
    {
      print_hypre_parameters(hypre_preconditioner_pt);
    }

    return hypre_preconditioner_pt;
  } // create_hypre_preconditioner




  inline void specify_command_line_flags()
  {
    // Flag to output the preconditioner, used for debugging.
    // string
    //    CommandLineArgs::specify_command_line_flag(
    //        "--doc_prec",&Doc_prec_dir_str);

    // No parameter after.
    //    CommandLineArgs::specify_command_line_flag(
    //        "--lsc_only");

    // double
    CommandLineArgs::specify_command_line_flag(
        "--sigma",&Scaling_sigma);

    // double
    CommandLineArgs::specify_command_line_flag(
        "--sigma_mult",&Scaling_sigma_multiplier);


    // int
    CommandLineArgs::specify_command_line_flag(
        "--w_solver",&W_solver);

    // Nothing set
    CommandLineArgs::specify_command_line_flag("--bdw");

    // int
    CommandLineArgs::specify_command_line_flag(
        "--ns_solver",&NS_solver);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--f_solver",&F_solver);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--p_solver",&P_solver);


    // NS_F block AMG parameters
    // int
    CommandLineArgs::specify_command_line_flag(
        "--f_amg_iter",&F_amg_param.Iterations);

    //int
    CommandLineArgs::specify_command_line_flag("--f_amg_smiter",
        &F_amg_param.Smoother_iterations);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--f_amg_sim_smoo",&F_amg_param.Simple_smoother);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--f_amg_com_smoo",&F_amg_param.Complex_smoother);

    // double
    CommandLineArgs::specify_command_line_flag(
        "--f_amg_damp",&F_amg_param.Damping);

    // double
    CommandLineArgs::specify_command_line_flag(
        "--f_amg_str",&F_amg_param.Strength);

    // int
    CommandLineArgs::specify_command_line_flag("--f_amg_coarse",
        &F_amg_param.Coarsening);

    CommandLineArgs::specify_command_line_flag("--print_f_hypre");


    // NS_F block AMG parameters
    // int
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_iter",&P_amg_param.Iterations);

    //int
    CommandLineArgs::specify_command_line_flag("--p_amg_smiter",
        &P_amg_param.Smoother_iterations);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_sim_smoo",&P_amg_param.Simple_smoother);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_com_smoo",&P_amg_param.Complex_smoother);

    // double
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_damp",&P_amg_param.Damping);

    // double
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_str",&P_amg_param.Strength);

    // int
    CommandLineArgs::specify_command_line_flag("--p_amg_coarse",
        &P_amg_param.Coarsening);

    CommandLineArgs::specify_command_line_flag("--print_p_hypre");


    // For general print hypres.
    CommandLineArgs::specify_command_line_flag("--print_hypre");
  }


  inline void setup_command_line_flags()
  {
    if(CommandLineArgs::command_line_flag_has_been_set("--print_f_hypre"))
    {
      F_amg_param.Print_parameters = true;
    }
    if(CommandLineArgs::command_line_flag_has_been_set("--print_p_hypre"))
    {
      P_amg_param.Print_parameters = true;
    }
  }

  inline std::string 
    get_hypre_preconditioner_string(HyprePreconditioner* h_prec_pt)
    {

      std::string hypre_string = "";
      {

        // Determine the cycle
        std::stringstream cycle_stream;
        cycle_stream << h_prec_pt->amg_iterations() << "v" 
          << h_prec_pt->amg_smoother_iterations();

        hypre_string += cycle_stream.str();

        // Determine strength of dependence
        std::stringstream strn_stream;
        strn_stream << "Strn" << h_prec_pt->amg_strength();

        hypre_string += strn_stream.str();

        // Determine the coarsening strategy
        const int amg_coarsening = h_prec_pt->amg_coarsening();
        std::string amg_coarsening_str = "";
        switch (amg_coarsening)
        {
          case 0:
            amg_coarsening_str = "CLPJ";
            break;
          case 1:
            amg_coarsening_str = "cRS";
            break;
          case 3:
            amg_coarsening_str = "mRS";
            break;
          case 6:
            amg_coarsening_str = "Falgout";
            break;
          case 8:
            amg_coarsening_str = "PMIS";
            break;
          case 10:
            amg_coarsening_str = "HMIS";
            break;
          case 11:
            amg_coarsening_str = "oRS";
            break;
          default:            
            {
              std::ostringstream err_msg;
              err_msg << "Something wrong setting coarsening string."
                << "amg_coarsening is " << amg_coarsening << ".\n"
                << std::endl;

              throw OomphLibError(err_msg.str(),
                  OOMPH_CURRENT_FUNCTION,
                  OOMPH_EXCEPTION_LOCATION);
            }
            break;
        }

        hypre_string += amg_coarsening_str;

        // Now do the smoother string.
        const bool amg_using_simple_smoothing_flag 
          = h_prec_pt->amg_using_simple_smoothing_flag();

        const int amg_simple_smoother = h_prec_pt->amg_simple_smoother();
        const int amg_complex_smoother = h_prec_pt->amg_complex_smoother();
        std::string amg_smoother_str = "";

        if(amg_using_simple_smoothing_flag)
        {
          if(amg_simple_smoother == 0)
          {
            amg_smoother_str = "Jac";

            // If it is Jacobi, we have the damping!
            std::ostringstream dampstream;
            dampstream << h_prec_pt->amg_damping();
            amg_smoother_str += dampstream.str();
          }
          else if(amg_simple_smoother == 1)
          {
            amg_smoother_str = "Gs";
          }
          else if(amg_simple_smoother == 2)
          {
            amg_smoother_str = "Gspinter";
          }
          else if(amg_simple_smoother == 3)
          {
            amg_smoother_str = "SORfs";
            std::ostringstream dampstream;
            dampstream << h_prec_pt->amg_damping();
            amg_smoother_str += dampstream.str();
          }
          else if(amg_simple_smoother == 4)
          {
            amg_smoother_str = "SORbs";
            std::ostringstream dampstream;
            dampstream << h_prec_pt->amg_damping();
            amg_smoother_str += dampstream.str();
          }
          else if(amg_simple_smoother == 6)
          {
            amg_smoother_str = "SSOR";
            std::ostringstream dampstream;
            dampstream << h_prec_pt->amg_damping();
            amg_smoother_str += dampstream.str();
          }
          else
          {
            std::ostringstream err_msg;
            err_msg << "Something wrong setting smoother string."
              << "amg_simple_smoother is " << amg_simple_smoother << ".\n"
              << std::endl;

            throw OomphLibError(err_msg.str(),
                OOMPH_CURRENT_FUNCTION,
                OOMPH_EXCEPTION_LOCATION);
          }
        }
        else
        {
          if(amg_complex_smoother == 6)
          {
            amg_smoother_str = "Schwarz";
          }
          else if(amg_complex_smoother == 7)
          {
            amg_smoother_str = "Pilut";
          }
          else if(amg_complex_smoother == 8)
          {
            amg_smoother_str = "ParaSails";
          }
          else if(amg_complex_smoother == 9)
          {
            amg_smoother_str = "Euclid";
          }
          else
          {
            std::ostringstream err_msg;
            err_msg << "Something wrong setting COMPLEX smoother string.\n"
              << std::endl;

            throw OomphLibError(err_msg.str(),
                OOMPH_CURRENT_FUNCTION,
                OOMPH_EXCEPTION_LOCATION);
          }
        }

        hypre_string += amg_smoother_str;
      }
  std::replace(hypre_string.begin(), hypre_string.end(),
                '.','-');
      return hypre_string;
    }

  inline Preconditioner* 
    create_f_p_amg_preconditioner(AMGParameters& amg_param,
        const int& f_or_p)
    {
      // f_or_p = 0 means do f, 1 means do p
      if(f_or_p == 0)
      {
        if((F_solver == F_solver_exact) || (NS_solver == NS_solver_exact))
        {
          NS_f_prec_str = "";
          return 0;
        }
        else if(F_solver == F_solver_amg)
        {
          if(amg_param.Print_parameters)
          {
            oomph_info << "RAYINFO: F amg parameters are: " << std::endl; 
          }
          
          HyprePreconditioner* h_prec_pt 
            = create_hypre_preconditioner(amg_param);
          NS_f_prec_str = get_hypre_preconditioner_string(h_prec_pt);
          return h_prec_pt;
        }
        else if(F_solver == F_solver_blkdiag_amg)
        {
#ifdef OOMPH_HAS_HYPRE
      // Create a block diagonal preconditioner.
      Preconditioner* f_preconditioner_pt =
        new BlockDiagonalPreconditioner<CRDoubleMatrix>;

      // Now, since f_precondtioner_pt is a Preconditioner*, it needs to be
      // caste to a block diagonal one if we want to call functions from
      // that class.
      dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
      (f_preconditioner_pt)->set_subsidiary_preconditioner_function
      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_JhalfStrnSimOneVTwoTwoRS);


      // All done. The below is prints the Hypre settings.
      //if(amg_param.Print_parameters)
      {
        oomph_info << "RAYINFO: F amg parameters are: " << std::endl; 
        oomph_info << "RAYINFO: Doing DIAGONAL " << std::endl; 
      

      // Check the Hypre values used, we encapsulate this so we can easily 
      // take it out later.
  
        // Create a new preconditioner with the above function we set.
        Preconditioner* check_prec_pt = 
          Hypre_Subsidiary_Preconditioner_Helper::
          set_hypre_JhalfStrnSimOneVTwoTwoRS();
        // Convert the above to a hypre preconditioner.
       HyprePreconditioner* h_prec_pt 
         = checked_dynamic_cast<HyprePreconditioner*>(check_prec_pt);

        // Now print it out to see the settings!
        print_hypre_parameters(h_prec_pt);
        NS_f_prec_str = get_hypre_preconditioner_string(h_prec_pt);
      }

      return f_preconditioner_pt;
#endif
 
        }
        else if(F_solver == F_solver_upper_amg)
        {
#ifdef OOMPH_HAS_HYPRE

      // Create a block triangular preconditioner.
      Preconditioner* f_preconditioner_pt =
        new BlockTriangularPreconditioner<CRDoubleMatrix>;

      // Use upper triangular preconditioner.
      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
      (f_preconditioner_pt)->upper_triangular();

      // Set the Hypre preconditioner.
      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
      (f_preconditioner_pt)->set_subsidiary_preconditioner_function
      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_JhalfStrnSimOneVTwoTwoRS);

      // All done. The below is prints the Hypre settings.
      //if(amg_param.Print_parameters)
      {
        oomph_info << "RAYINFO: F amg parameters are: " << std::endl; 
      

      // Check the Hypre values used, we encapsulate this so we can easily 
      // take it out later.
  
        // Create a new preconditioner with the above function we set.
        Preconditioner* check_prec_pt = 
          Hypre_Subsidiary_Preconditioner_Helper::
          set_hypre_JhalfStrnSimOneVTwoTwoRS();
        // Convert the above to a hypre preconditioner.
       HyprePreconditioner* h_prec_pt 
         = checked_dynamic_cast<HyprePreconditioner*>(check_prec_pt);

        // Now print it out to see the settings!
        print_hypre_parameters(h_prec_pt);
        oomph_info << "RAYINFO: Doing UPPER triangular AMG" << std::endl;
        NS_f_prec_str = get_hypre_preconditioner_string(h_prec_pt);
      }

      return f_preconditioner_pt;

#endif
 
        }
        else if(F_solver == F_solver_lower_amg)
        {
#ifdef OOMPH_HAS_HYPRE

      // Create a block triangular preconditioner.
      Preconditioner* f_preconditioner_pt =
        new BlockTriangularPreconditioner<CRDoubleMatrix>;

      // Use upper triangular preconditioner.
      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
      (f_preconditioner_pt)->lower_triangular();

      // Set the Hypre preconditioner.
      dynamic_cast<BlockTriangularPreconditioner<CRDoubleMatrix>* >
      (f_preconditioner_pt)->set_subsidiary_preconditioner_function
      (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_JhalfStrnSimOneVTwoTwoRS);

      // All done. The below is prints the Hypre settings.
      //if(amg_param.Print_parameters)
      {
        oomph_info << "RAYINFO: F amg parameters are: " << std::endl; 
      

      // Check the Hypre values used, we encapsulate this so we can easily 
      // take it out later.
  
        // Create a new preconditioner with the above function we set.
        Preconditioner* check_prec_pt = 
          Hypre_Subsidiary_Preconditioner_Helper::
          set_hypre_JhalfStrnSimOneVTwoTwoRS();
        // Convert the above to a hypre preconditioner.
       HyprePreconditioner* h_prec_pt 
         = checked_dynamic_cast<HyprePreconditioner*>(check_prec_pt);

        // Now print it out to see the settings!
        oomph_info << "RAYINFO: Doing LOWER triangular AMG" << std::endl;
        print_hypre_parameters(h_prec_pt);
        NS_f_prec_str = get_hypre_preconditioner_string(h_prec_pt);
      }

      return f_preconditioner_pt;

#endif
 
        }
        else if(F_solver == F_solver_full_amg)
        {
      Preconditioner* f_preconditioner_pt 
        = Hypre_Subsidiary_Preconditioner_Helper::
          set_hypre_JhalfStrnSimOneVTwoTwoRS();
        // Convert the above to a hypre preconditioner.
       HyprePreconditioner* h_prec_pt 
         = checked_dynamic_cast<HyprePreconditioner*>(
             f_preconditioner_pt);

        // Now print it out to see the settings!
        oomph_info << "RAYINFO: F amg parameters are: " << std::endl; 
        oomph_info << "RAYINFO: Doing FULL AMG" << std::endl;
        print_hypre_parameters(h_prec_pt);
        NS_f_prec_str = get_hypre_preconditioner_string(h_prec_pt);


      return f_preconditioner_pt;
        }

      }
      else
      {
        if((P_solver == P_solver_exact) || (NS_solver == NS_solver_exact))
        {
          NS_p_prec_str = "";
          return 0;
        }
        else
        {
          if(amg_param.Print_parameters)
          {
            oomph_info << "RAYINFO: P amg parameters are: " << std::endl;
          }
          HyprePreconditioner* h_prec_pt
            = create_hypre_preconditioner(amg_param);
          NS_p_prec_str = get_hypre_preconditioner_string(h_prec_pt);

          return h_prec_pt;
        }
      }

    }



  inline NavierStokesSchurComplementPreconditioner*
    create_lsc_preconditioner(Problem* problem_pt,
        Mesh* mesh_pt,
        Preconditioner* f_prec_pt,
        Preconditioner* p_prec_pt)
    {
      NS_prec_str = "";

      // Create the NS LSC preconditioner.
      NavierStokesSchurComplementPreconditioner* ns_prec_pt = 0;

      if(NS_solver != NS_solver_exact)
      {
        ns_prec_pt = new NavierStokesSchurComplementPreconditioner(problem_pt);

        // Give LSC the bulk mesh (Navier-Stokes mesh).
        ns_prec_pt->set_navier_stokes_mesh(mesh_pt);

        if(f_prec_pt != 0)
        {
          ns_prec_pt->set_f_preconditioner(f_prec_pt);
          NS_prec_str+="Fa" + NS_f_prec_str;
        }
        else
        {
          NS_prec_str+="Fe" + NS_f_prec_str;
        }

        if(p_prec_pt != 0)
        {
          ns_prec_pt->set_p_preconditioner(p_prec_pt);
          NS_prec_str+="Pa" + NS_p_prec_str;
        }
        else
        {
          NS_prec_str+="Pe" + NS_p_prec_str;
        }
      }

      return ns_prec_pt;
    }


  inline LagrangeEnforcedflowPreconditioner* 
    create_lgr_precondiitoner(Vector<Mesh*>& mesh_pt,
        const int w_solver, 
        Preconditioner* ns_prec_pt)
    {
      Lgr_prec_str = "";

      LagrangeEnforcedflowPreconditioner* lgr_prec_pt
        = new LagrangeEnforcedflowPreconditioner;

      // First set the meshes
      lgr_prec_pt->set_meshes(mesh_pt);

      if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
      {
        lgr_prec_pt->scaling_sigma() = Scaling_sigma;
        std::ostringstream strs;
        strs << "Sig" << Scaling_sigma;
        Lgr_prec_str += strs.str();
      }

      if(CommandLineArgs::command_line_flag_has_been_set("--sigma_mult"))
      {
        lgr_prec_pt->scaling_sigma_multiplier() = Scaling_sigma_multiplier;

        std::ostringstream strs;
        strs << "Sigmult" << Scaling_sigma_multiplier;
        Lgr_prec_str += strs.str();
      }

      if(CommandLineArgs::command_line_flag_has_been_set("--bdw"))
      {
        lgr_prec_pt->use_block_diagonal_w_block();
      }

      // Now do the W solver
      if(w_solver == -1)
      {
        std::ostringstream err_msg;
        err_msg << "There W_solver has not been set.\n"
          << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      else if(w_solver == 0)
      {
        Lgr_prec_str += "We";
        // Using SuperLU, this is the default, do nothing.
      }
      else if (w_solver ==1)
      {
        lgr_prec_pt->set_lagrange_multiplier_subsidiary_preconditioner
          (Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper
           ::get_lagrange_multiplier_preconditioner);
        Lgr_prec_str+="Wcg";
      }
      else
      {
        std::ostringstream err_msg;
        err_msg << "There is no other W solver set.\n"
          << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }


      // Now set the NS solver
      if(ns_prec_pt != 0)
      {
        lgr_prec_pt->set_navier_stokes_lsc_preconditioner(ns_prec_pt);
        Lgr_prec_str+="Nl" + NS_prec_str;
      }
      else
      {
        Lgr_prec_str+="Ne";
      }
      return lgr_prec_pt;
    }


}





#endif



