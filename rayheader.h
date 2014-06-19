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
//
// This is a header for shared data I keep.
//  - Ray
//
#ifndef RAY_LAGRANGE_HEADER
#define RAY_LAGRANGE_HEADER

// Oomph-lib includes
#include "generic.h"

using namespace oomph;

namespace RayGen
{
  //=============================================================================
  /// Check if a number is in a set.
  /// Input: myarray = The set of numbers
  ///        nval = The number of values in the set
  ///        number = The number to check
  ///
  /// Example: Checking if the number 91 is in the set 
  ///          10, 11, 12, 13, 20, 21, 22, 23
  ///   int prob_id_array[]= {10,11,12,13,
  ///                         20,21,22,23};
  ///
  ///   bool inset = check_if_in_set<int>(prob_id_array,8,91);
  ///
  //=============================================================================
  template <typename T>
    bool check_if_in_set(const T myarray[], const unsigned& nval, const T number)
    {
      // Convert the array of numbers into a set.
      std::set<T> myset(myarray,myarray+nval);

      // Create a pair of iterator and boolean. Must declare std::set<T>::iterator
      // as a typename because the compiler gets confused. It is not templated at 
      // compile time... so the compiler thinks it's a static field and will
      // complain.
      std::pair<typename std::set<T>::iterator,bool> ret;

      // The second in the pair is false if not added. 
      // This means it is in the set, so we negate this to make it true.
      ret = myset.insert(number);

      return (!ret.second);
    } // check_if_in_set(...)

  bool check_amg_coarsener(const int &amg_coarsener)
  {
    static const unsigned n_coarseners = 7;
    static const int coarsener_array[] = {0,1,3,6,8,10,11};
    return check_if_in_set<int>(coarsener_array,n_coarseners,amg_coarsener);
  }
  
  bool check_amg_sim_smoother(const int &amg_sim_smoother)
  {
    static const unsigned n_sim_smoothers = 6;
    static const int sim_smoother_array[] = {0,1,2,3,4,6};
    return check_if_in_set<int>(sim_smoother_array,n_sim_smoothers,
                                amg_sim_smoother);
  }

  bool check_amg_com_smoother(const int &amg_com_smoother)
  {
    static const unsigned n_com_smoothers = 4;
    static const int com_smoother_array[] = {6,7,8,9};
    return check_if_in_set<int>(com_smoother_array,n_com_smoothers,
                                amg_com_smoother);
  }

  bool check_amg_sim_smoother_damping_required(const int &amg_sim_smoother)
  {
    static const unsigned n_sim_smoothers = 4;
    static const int sim_smoother_array[] = {0,3,4,6};
    return check_if_in_set<int>(sim_smoother_array,n_sim_smoothers,
                                amg_sim_smoother);
  }
}


#ifdef OOMPH_HAS_HYPRE
//=============================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//=============================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
  void check_amg_parameters(
      const int& amg_iterations, const int& amg_smoother_iterations,
      const int& amg_simple_smoother, const int& amg_complex_smoother,
      const double& amg_damping, const double& amg_strength,
      const int& amg_coarsening)
  {
    bool coarsening_ok = RayGen::check_amg_coarsener(amg_coarsening);
    if((amg_coarsening < 0) || (!coarsening_ok))
    {
      std::ostringstream err_msg;
      err_msg 
        << "Please set a coarsening strategy with --x_amg_coarse.\n"
        << "Current coarsening strategy (amg_coarsening) is: " 
        << amg_coarsening << "\n\n"

        << "You have either not set it or it is not valid.\n\n"

        << "Valid IDs are:\n"
        << "0 - CLJP\n"
        << "1 - Classical RS\n"
        << "3 - modified RS\n"
        << "6 - Falgout (default in documentation)\n"
        << "8 - PMIS\n"
        << "10 - HMIS\n"
        << "11 - One pass on RS coarsening on each processor, "
        << "not recommended.\n"
        << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // First check if the user has selected both simple and complex
    // smoothing, this is not allowed, either one or the other.
    if((amg_simple_smoother >= 0) && 
        (amg_complex_smoother >= 0))
    {
      std::ostringstream err_msg;
      err_msg 
        << "Both simple and complex smoothers are set.\n"
        << "Please choose one or the other.\n\n" 
        << "amg_simple_smoother is " << amg_simple_smoother << "\n"
        << "amg_complex_smoother is " << amg_complex_smoother << "\n\n"

        << "Simple smoother IDs, set with --x_amg_sim_smoo\n"
        << "0 - Jacobi (Need to set damping as well)\n"
        << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
        << "2 - GS - interior parallel, serial on boundary\n"
        << "3 - hybrid GS or SOR, forward solve\n"
        << "4 - hybrid GS or SOR, backwards solve\n"
        << "6 - hybrid symmetric GS or SSOR.\n\n"

        << "Complex smoother IDs, set with --x_amg_com_smoo\n"
        << "6 - Schwarz (default in documentation)\n"
        << "7 - Pilut\n"
        << "8 - ParaSails\n"
        << "9 - Euclid\n" 
        << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else if((amg_simple_smoother < 0) && (amg_complex_smoother < 0))
    {
      std::ostringstream err_msg;
      err_msg 
        << "Please select a smoother for the f block.\n"
        << "Use --x_amg_sim_smoo or --x_amg_com_smoo flag.\n\n"

        << "Simple smoother IDs, set with --x_amg_sim_smoo\n"
        << "0 - Jacobi (Need to set damping as well)\n"
        << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
        << "2 - GS - interior parallel, serial on boundary\n"
        << "3 - hybrid GS or SOR, forward solve (damping needed?)\n"
        << "4 - hybrid GS or SOR, backwards solve (damping needed?)\n"
        << "6 - hybrid symmetric GS or SSOR (damping needed?).\n\n"

        << "Complex smoother IDs, set with --x_amg_com_smoo\n"
        << "6 - Schwarz (default in documentation)\n"
        << "7 - Pilut\n"
        << "8 - ParaSails\n"
        << "9 - Euclid\n" 
        << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }


    // Only one of the two smoothers have been set. We see if these are
    // valid smoothing IDs.
    if(amg_simple_smoother >= 0)
    {
      // check if simple smoother is okay.
      bool sim_smoother_ok 
        = RayGen::check_amg_sim_smoother(amg_simple_smoother);

      if(!sim_smoother_ok)
      {
        std::ostringstream err_msg;
        err_msg 
          << "Please provide a valid simple smoother \n"
          << "using --x_amg_sim_smoo. You have provided: " 
          << amg_simple_smoother <<"\n\n"
          << "Valid IDs are:\n"
          << "0 - Jacobi (Need to set damping as well)\n"
          << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
          << "2 - GS - interior parallel, serial on boundary\n"
          << "3 - hybrid GS or SOR, forward solve\n"
          << "4 - hybrid GS or SOR, backwards solve\n"
          << "6 - hybrid symmetric GS or SSOR.\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
    }
    else if(amg_complex_smoother >=0)
    {
      // check if complex smoother is valid.
      bool com_smoother_ok
        = RayGen::check_amg_com_smoother(amg_complex_smoother);

      if(!com_smoother_ok)
      {
        std::ostringstream err_msg;
        err_msg 
          << "Please provide a valid complex smoother \n"
          << "using --x_amg_com_smoo. You have provided: " 
          << amg_complex_smoother <<"\n\n"
          << "Valid IDs are:\n"
          << "6 - Schwarz (default in documentation)\n"
          << "7 - Pilut\n"
          << "8 - ParaSails\n"
          << "9 - Euclid\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
    }
    else
    {
      // Should never get here, something has gone wrong.
      std::ostringstream err_msg;
      err_msg 
        << "Something when wrong with your selection of smoother.\n"
        << "Choose to use either a simple smoother or complex smoother\n"
        << "with the flags --x_amg_sim_smoo and --x_amg_com_smoo.\n\n"

        << "Current amg_simple_smoother is " 
        << amg_simple_smoother << "\n"
        << "Current amg_complex_smoother is " 
        << amg_complex_smoother << "\n\n"

        << "Simple smoother IDs, set with --x_amg_sim_smoo\n"
        << "0 - Jacobi (Need to set damping as well)\n"
        << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
        << "2 - GS - interior parallel, serial on boundary\n"
        << "3 - hybrid GS or SOR, forward solve\n"
        << "4 - hybrid GS or SOR, backwards solve\n"
        << "6 - hybrid symmetric GS or SSOR.\n\n"

        << "Complex smoother IDs, set with --x_amg_com_smoo\n"
        << "6 - Schwarz (default in documentation)\n"
        << "7 - Pilut\n"
        << "8 - ParaSails\n"
        << "9 - Euclid\n" 
        << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Damping is required for Jacobi or hybrid SOR smoothers.
    // First check if amg_simple_smoother is one of those.
    bool damping_required
      = RayGen::check_amg_sim_smoother_damping_required(amg_simple_smoother);

    if(damping_required && (amg_damping < 0))
    {
      std::ostringstream err_msg;
      err_msg 
        << "Have selected simple smoother: " << amg_simple_smoother << "\n"
        << "Damping parameter is required for this smoother.\n"
        << "Please set it via --x_amg_damp.\n"
        << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    // Note that we haven't checked the reverse, i.e. if a smoother
    // which does not required damping is set, and the damping parameter
    // is set, we allow this... if implemented properly, the hypre
    // code should just ignore it. It will also be interesting to find
    // smoothers which the damping parameter affects the result.

    // Now check if the strength parameter is set. 
    if(amg_strength < 0)
    {
      std::ostringstream err_msg;
      err_msg 
        << "You have not set the strength parameter via --x_amg_str\n"
        << std::endl;

      throw OomphLibWarning(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Now check the f_amg_iterations
    if(amg_iterations < 0)
    {
      std::ostringstream err_msg;
      err_msg 
        << "Have not set f_amg_iterations via --x_amg_iter\n"
        << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    if(amg_smoother_iterations < 0)
    {
      std::ostringstream err_msg;
      err_msg 
        << "Have not set f_amg_smoother_iterations via --x_amg_smiter\n"
        << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  }

  //===========================================================================
  /// Given a preconditioner, it will print out the hypre settings.
  //===========================================================================
  void print_hypre_settings(Preconditioner* preconditioner_pt)
  {
    HyprePreconditioner* h_prec_pt 
      = checked_dynamic_cast<HyprePreconditioner*>(preconditioner_pt);

    // Print AMG iterations:
    std::cout << "HyprePreconditioner settings are: " << std::endl;
    std::cout << "Max_iter: " << h_prec_pt->amg_iterations() << std::endl;
    std::cout << "smoother iter: " << h_prec_pt->amg_smoother_iterations() 
      << std::endl;
    std::cout << "Hypre_method: " << h_prec_pt->hypre_method() << std::endl;
    std::cout << "internal_preconditioner: " 
      << h_prec_pt->internal_preconditioner() << std::endl;
    std::cout << "AMG_using_simple_smoothing: " 
      << h_prec_pt->amg_using_simple_smoothing_flag() << std::endl; 
    std::cout << "AMG_simple_smoother: " 
      << h_prec_pt->amg_simple_smoother() << std::endl; 
    std::cout << "AMG_complex_smoother: " 
      << h_prec_pt->amg_complex_smoother() << std::endl; 
    std::cout << "AMG_coarsening: " 
      << h_prec_pt->amg_coarsening() << std::endl;
    std::cout << "AMG_max_levels: " 
      << h_prec_pt->amg_max_levels() << std::endl;
    std::cout << "AMG_damping: " 
      << h_prec_pt->amg_damping() << std::endl;
    std::cout << "AMG_strength: " 
      << h_prec_pt->amg_strength() << std::endl;;
    std::cout << "AMG_max_row_sum: " 
      << h_prec_pt->amg_max_row_sum() << std::endl;
    std::cout << "AMG_truncation: " 
      << h_prec_pt->amg_truncation() << std::endl;
    std::cout << "\n" << std::endl;
  }

  // This function assumes that the global amg parameters are set
  // correctly and will do minimal error checking.
  std::string get_amg_label(Preconditioner* preconditioner_pt)
  {
    // Cast to Hypre preconditioner.
    HyprePreconditioner* h_prec_pt 
      = checked_dynamic_cast<HyprePreconditioner*>(preconditioner_pt);

    // Start extracting the parameters.
    const int amg_iterations = h_prec_pt->amg_iterations();
    const int amg_smoother_iterations = h_prec_pt->amg_smoother_iterations();

    // One of these have to remain -1, for logic tests below.
    int amg_simple_smoother = -1;
    int amg_complex_smoother = -1;

    if(h_prec_pt->amg_using_simple_smoothing_flag())
    {
      amg_simple_smoother = h_prec_pt->amg_simple_smoother();
    }
    else
    {
      amg_complex_smoother = h_prec_pt->amg_complex_smoother();
    }

    // Continue getting the rest of the parameters.
    const double amg_damping = h_prec_pt->amg_damping();
    const double amg_strength = h_prec_pt->amg_strength();
    const int amg_coarsening = h_prec_pt->amg_coarsening();

    // String to glue together at the end.
    std::string strength_str = "";
    std::string coarsening_str = "";
    std::string smoother_str = "";
    std::string cycle_str = "";

    // Set the cycle.
    std::stringstream cycle_ss;
    cycle_ss << "It" << amg_iterations
      << "Sit" << amg_smoother_iterations;
    cycle_str = cycle_ss.str();

    // Set the smoothers
    //   0 = Jacobi 
    //   1 = Gauss-Seidel, sequential
    //       (very slow in parallel!)
    //   2 = Gauss-Seidel, interior points in parallel, boundary sequential
    //       (slow in parallel!)
    //   3 = hybrid Gauss-Seidel or SOR, forward solve
    //   4 = hybrid Gauss-Seidel or SOR, backward solve
    //   6 = hybrid symmetric Gauss-Seidel or SSOR 
    std::stringstream smoother_ss;

    switch (amg_simple_smoother)
    {
      case 0:
        smoother_ss << "Jac";
        break;
      case 1:
        smoother_ss << "GSs";
        break;
      case 2:
        smoother_ss << "GSp";
        break;
      case 3:
        smoother_ss << "SORf";
        break;
      case 4:
        smoother_ss << "SORb";
        break;
      case 6:
        smoother_ss << "SSOR";
        break;
    }

    // Set the damping if it's set (greater than 0)
    if(amg_damping >= 0)
    {
      smoother_ss << "D" << amg_damping;
    }

    // Set complex smoother string
    //   6 = Schwarz
    //   7 = Pilut
    //   8 = ParaSails
    //   9 = Euclid
    switch (amg_complex_smoother)
    {
      case 6:
        {
          smoother_ss << "Schwarz";
        }
        break;
      case 7:
        {
          smoother_ss << "Pilut";
        }
        break;
      case 8:
        {
          smoother_ss << "ParaSails";
        }
        break;
      case 9:
        {
          smoother_ss << "Euclid";
        }
        break;
    }

    smoother_str = smoother_ss.str();

    // Now set the coarsening string
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
    std::stringstream coarsening_ss;
    switch (amg_coarsening)
    {
      case 0:
        coarsening_ss << "CLJP";
        break;
      case 1:
        coarsening_ss << "RSc";
        break;
      case 3:
        coarsening_ss << "RSm";
        break;
      case 6:
        coarsening_ss << "Falgout";
        break;
      case 8:
        coarsening_ss << "PMIS";
        break;
      case 10:
        coarsening_ss << "HMIS";
        break;
      case 11:
        coarsening_ss << "RS1p";
        break;
    }
    coarsening_str = coarsening_ss.str();

    // Create the string for strength
    std::stringstream strength_ss;
    strength_ss << "Strn" << amg_strength;
    strength_str = strength_ss.str();

    // Concatenate everything together
    std::string final_str = cycle_str + "_" + 
      smoother_str + "_" +
      coarsening_str + "_" +
      strength_str;

    return final_str;
  }

  //===========================================================================
  /// Quick reminder of the AMG parameters:
  ///
  /// \short Default for AMG strength (0.25 recommended for 2D problems;
  /// larger (0.5-0.75, say) for 3D
  /// extern double AMG_strength;
  ///
  /// \short Default AMG coarsening strategy. Coarsening types include:
  ///    0 = CLJP (parallel coarsening using independent sets)
  ///    1 = classical RS with no boundary treatment (not recommended
  ///        in parallel)
  ///    3 = modified RS with 3rd pass to add C points on the boundaries
  ///    6 = Falgout (uses 1 then CLJP using interior coarse points as
  ///        first independent set) THIS IS DEFAULT ON DOCUMENTATION
  ///    8 = PMIS (parallel coarsening using independent sets - lower
  ///        complexities than 0, maybe also slower convergence)
  ///    10= HMIS (one pass RS on each processor then PMIS on interior
  ///        coarse points as first independent set)
  ///    11= One pass RS on each processor (not recommended)
  ///    extern unsigned AMG_coarsening;/
  ///
  ///
  /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
  /// include:
  ///   0 = Jacobi 
  ///   1 = Gauss-Seidel, sequential
  ///       (very slow in parallel!)
  ///   2 = Gauss-Seidel, interior points in parallel, boundary sequential
  ///       (slow in parallel!)
  ///   3 = hybrid Gauss-Seidel or SOR, forward solve
  ///   4 = hybrid Gauss-Seidel or SOR, backward solve
  ///   6 = hybrid symmetric Gauss-Seidel or SSOR
  /// To use these methods set AMG_using_simple_smoothing to true
  /// unsigned AMG_simple_smoother;
  /// 
  /// \short Complex smoothing methods used in BoomerAMG. Relaxation types
  /// are:
  ///   6 = Schwarz
  ///   7 = Pilut
  ///   8 = ParaSails
  ///   9 = Euclid
  /// To use these methods set AMG_using_simple_smoothing to false
  /// unsigned AMG_complex_smoother; 
  //===========================================================================

  //===========================================================================
  /// Preconditioner* set_hypre_for_2D_poison_problem()
  ///
  /// Returns a hypre preconditioner with the following settings:
  ///
  ///   // Set iterations to 1
  ///   hypre_preconditioner_pt->set_amg_iterations(1);
  ///
  ///   // Use simple smoother
  ///   hypre_preconditioner_pt->amg_using_simple_smoothing();
  ///
  ///   // Smoother types:
  ///           0=Jacobi
  ///           1 = Gauss-Seidel, sequential
  ///   hypre_preconditioner_pt->amg_simple_smoother() = 1;
  ///
  ///  // AMG preconditioner
  ///  hypre_preconditioner_pt->hypre_method() 
  ///    = HyprePreconditioner::BoomerAMG;
  ///
  ///  // Choose strength parameter for amg
  ///  hypre_preconditioner_pt->amg_strength() = 0.25;
  ///
  ///  // Coarsening type 0 = CLJP
  ///  hypre_preconditioner_pt->amg_coarsening() = 0; 
  //===========================================================================
  Preconditioner* set_hypre_for_2D_poison_problem()
  {
    // Create a new hypre preconditioner
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
 //   HyprePreconditioner* hypre_preconditioner_pt = 
 //     checked_static_cast<HyprePreconditioner*>(another_preconditioner_pt);

 //   // 2D Poisson problem defaults as defined above.
 //   Hypre_default_settings::
 //     set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

 //   // If the strength parameter has been set (>=0), then we set it.
 //   if(RayGlobalAMGParam::amg_strength >= 0)
 //   {
 //     hypre_preconditioner_pt->amg_strength() 
 //       = RayGlobalAMGParam::amg_strength;
 //   }

 //   // reset the amg parameters.
 //   reset_amg_param();

 //   // Print to confirm?
 //   if(RayGlobalAMGParam::print_hypre == true)
 //   {
 //     print_hypre_settings(hypre_preconditioner_pt);
 //   }

    // Return the hypre preconditioner.
    return another_preconditioner_pt;
  } // set_hypre_for_2D_poison_problem()

  //===========================================================================
  /// Preconditioner* set_hypre_for_navier_stokes_momentum_block()
  ///
  /// Returns a hypre preconditioner with the following settings   
  ///
  ///  // Set default settings as for 2D Poisson
  ///  set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
  /// 
  ///  // Change smoother type:
  ///  //           0=Jacobi
  ///  //           1=Gauss-Seidel
  ///  hypre_preconditioner_pt->amg_simple_smoother() = 0;
  ///  
  ///  // Set smoother damping
  ///  hypre_preconditioner_pt->amg_damping() = 0.5;
  ///
  ///  // Change strength parameter for amg
  ///  hypre_preconditioner_pt->amg_strength() = 0.75;:
  //===========================================================================
  Preconditioner* set_hypre_for_navier_stokes_momentum_block()
  {
    // Create a new hypre preconditioner
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      checked_static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    // Defaults for Navier Stokes momentum block as defined above.
//    Hypre_default_settings::
//      set_defaults_for_navier_stokes_momentum_block(hypre_preconditioner_pt);
//
//    // Have the strength parameter been set?
//    if(RayGlobalAMGParam::amg_strength >= 0)
//    {
//      hypre_preconditioner_pt->amg_strength() 
//        = RayGlobalAMGParam::amg_strength;
//    }
//
//    // Has the damping parameter been set?
//    if(RayGlobalAMGParam::amg_damping >= 0)
//    {
//      hypre_preconditioner_pt->amg_damping() = RayGlobalAMGParam::amg_damping;
//    }
//
//    // Reset the amg parameters.
//    reset_amg_param();
//
//    // Print to confirm?
//    if(RayGlobalAMGParam::print_hypre == true)
//    {
//      print_hypre_settings(hypre_preconditioner_pt);
//    }

    // Return the newly created preconditioner.
    return another_preconditioner_pt;
  } // set_hypre_for_navier_stokes_momentum_block()

  //===========================================================================
  /// Preconditioner* get_custom_hypre_preconditioner()
  ///
  /// Returns a hypre preconditioner where all the settings must be set.
  ///
  //
     /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
   /// include:
   ///  0 = Jacobi 
   ///  1 = Gauss-Seidel, sequential
   ///      (very slow in parallel!)
   ///  2 = Gauss-Seidel, interior points in parallel, boundary sequential
   ///      (slow in parallel!)
   ///  3 = hybrid Gauss-Seidel or SOR, forward solve
   ///  4 = hybrid Gauss-Seidel or SOR, backward solve
   ///  6 = hybrid symmetric Gauss-Seidel or SSOR
   /// To use these methods set AMG_using_simple_smoothing to true
  //
    /// \short Complex smoothing methods used in BoomerAMG. Relaxation types
   /// are:
   ///  6 = Schwarz
   ///  7 = Pilut
   ///  8 = ParaSails
   ///  9 = Euclid
   /// To use these methods set AMG_using_simple_smoothing to false
  //
  //===========================================================================
  Preconditioner* get_custom_hypre_preconditioner(
      const int& amg_iterations, const int& amg_smoother_iterations,
      const int& amg_simple_smoother, const int& amg_complex_smoother,
      const double& amg_damping, const double& amg_strength,
      const int& amg_coarsening)
  {
    // Create a new hypre preconditioner
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      checked_static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    // Set the hypre_method to BoomerAMG. This is hard coded.
    hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

    // Check that all the parameters are valid and correct.
    check_amg_parameters(
        amg_iterations, amg_smoother_iterations, 
        amg_simple_smoother, amg_complex_smoother,
        amg_damping, amg_strength,
        amg_coarsening);

    // Set the amg_iterations, this is usually set to 1, for V(2,2)
    // there is only one of them... I think that's what I means, 
    // the two is the smoother... pre and post smoothing.
    hypre_preconditioner_pt
      ->set_amg_iterations(amg_iterations);

    // Set the number of amg smoother iterations. This is usually set to 2
    hypre_preconditioner_pt->amg_smoother_iterations() 
      = amg_smoother_iterations;

    // Now set the smoother type. Remember to use damping for Jacobi.
    if(amg_simple_smoother >= 0)
    {
      // Set simple smoothing flag
      hypre_preconditioner_pt->amg_using_simple_smoothing();

      // Set the simple smoother type
      hypre_preconditioner_pt->amg_simple_smoother() = amg_simple_smoother;
    }
    else if(amg_complex_smoother >= 0)
      // Otherwise we have complex smoothing
    {
      // Set the complex smoothing flag.
      hypre_preconditioner_pt->amg_using_complex_smoothing();

      // Set the complex smoother type.
      hypre_preconditioner_pt->amg_complex_smoother() = amg_complex_smoother;
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
    hypre_preconditioner_pt->amg_damping() = amg_damping;

    // Now set the AMG strength parameter.
    hypre_preconditioner_pt->amg_strength() = amg_strength;

    // AMG coarsening strategy.
    hypre_preconditioner_pt->amg_coarsening() = amg_coarsening;

    std::string amg_label_str = get_amg_label(hypre_preconditioner_pt);

    std::cout << "RAYHYPRE: " << amg_label_str << std::endl;

    return another_preconditioner_pt;
  } // get_custom_hypre_preconditioner

  ///////////////////////////////////////////////
  ///////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  Preconditioner* set_hypre_using_2D_poisson_base()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    string smoother_str = "";
//    string damping_str = "";
//    string coarsening_str = "";
//    string strength_str = "";
//    ///////////////////////////////////////////////////////////////////////////
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    // Setting the smoother:
//    if(RayGlobalAMGParam::amg_smoother == 0)
//    {
//      smoother_str = "GS";
//      // Setting up Gauss-Seidel
//      hypre_preconditioner_pt->amg_using_simple_smoothing();
//      hypre_preconditioner_pt->amg_simple_smoother() = 1;
//    }
//    else if(RayGlobalAMGParam::amg_smoother == 1)
//    {
//      smoother_str = "J";
//      hypre_preconditioner_pt->amg_damping() = RayGlobalAMGParam::amg_damping;
//
//      // Setting up Jacobi with damping.
//      hypre_preconditioner_pt->amg_using_simple_smoothing();
//      hypre_preconditioner_pt->amg_simple_smoother() = 0;
//      if(!(RayGlobalAMGParam::amg_damping < 0))
//      {
//        std::ostringstream strs;
//        strs << "Dmp" << hypre_preconditioner_pt->amg_damping();
//        damping_str = strs.str(); 
//      }
//      else
//      {
//        std::cout << "Please set your damping using --amg_damping" << std::endl; 
//        pause("Please do not continue."); 
//      }
//    }
//    else if(RayGlobalAMGParam::amg_smoother == 2)
//    {
//      smoother_str = "Pilut";
//      hypre_preconditioner_pt->amg_using_complex_smoothing();
//      hypre_preconditioner_pt->amg_complex_smoother() = 7;
//    }
//    else
//    {
//      std::cout << "You supplied smoother: " << RayGlobalAMGParam::amg_smoother << std::endl;
//      std::cout << "No such smoother. 0 is GS, 1 is Jacobi, 2 is Pilut" << std::endl;
//      std::cout << "Please set your smoother using --smoother" << std::endl;
//      pause("Please do not continue.");
//    }
//
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    if(RayGlobalAMGParam::amg_coarsening == 0)
//    {
//      coarsening_str = "CLJP";
//      hypre_preconditioner_pt->amg_coarsening() = 0;
//    }
//    else if(RayGlobalAMGParam::amg_coarsening == 1)
//    {
//      coarsening_str = "RS";
//      hypre_preconditioner_pt->amg_coarsening() = 1;
//    }
//    else
//    {
//      std::cout << "There is no such coarsening: " << RayGlobalAMGParam::amg_coarsening << std::endl;
//      std::cout << "0 - CLJP, 1 - RS, use --amg_coarsening" << std::endl;
//      pause("Do not continue"); 
//
//    }
//
//    if(!(RayGlobalAMGParam::amg_strength < 0))
//    {
//      hypre_preconditioner_pt->amg_strength() = RayGlobalAMGParam::amg_strength;
//      std::ostringstream strs;
//      strs << "Strn" << hypre_preconditioner_pt->amg_strength();
//      strength_str = strs.str(); 
//    }
//    else
//    {
//      std::cout << "Please set the amg_strengh using --amg_strength" << std::endl;
//      pause("Do not continue");
//    }
//
//    std::cout << "RAYHYPRE: " << coarsening_str 
//      << smoother_str
//      << damping_str
//      << strength_str
//      << std::endl;

    return another_preconditioner_pt;
  }
  /////////////////////////////////////////////////////////////////////////////

  Preconditioner* set_hypre_for_augmented_momentum_block()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_navier_stokes_momentum_block(hypre_preconditioner_pt);
//
//    hypre_preconditioner_pt->amg_strength() = 0.668;
//
//    hypre_preconditioner_pt->amg_damping() = 1.0;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_CLJPGSStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    ///////////////////////////////////////////////////////////////////////////
//    // Parameters from Milan.
//
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    //hypre_preconditioner_pt->amg_using_complex_smoothing();
//    //hypre_preconditioner_pt->amg_complex_smoother() = 7;
//    hypre_preconditioner_pt->amg_using_simple_smoothing();
//    hypre_preconditioner_pt->amg_simple_smoother() = 1;
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    hypre_preconditioner_pt->amg_coarsening() = 0;
//    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_RSGSStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    ///////////////////////////////////////////////////////////////////////////
//    // Parameters from Milan.
//
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    //hypre_preconditioner_pt->amg_using_complex_smoothing();
//    //hypre_preconditioner_pt->amg_complex_smoother() = 7;
//    hypre_preconditioner_pt->amg_using_simple_smoothing();
//    hypre_preconditioner_pt->amg_simple_smoother() = 1;
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    hypre_preconditioner_pt->amg_coarsening() = 1;
//    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_CLJPPilutStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    ///////////////////////////////////////////////////////////////////////////
//    // Parameters from Milan.
//
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    hypre_preconditioner_pt->amg_using_complex_smoothing();
//    hypre_preconditioner_pt->amg_complex_smoother() = 7;
//    //hypre_preconditioner_pt->amg_using_simple_smoothing();
//    //hypre_preconditioner_pt->amg_simple_smoother() = 1;
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    hypre_preconditioner_pt->amg_coarsening() = 0;
//    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_RSPilutStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    ///////////////////////////////////////////////////////////////////////////
//    // Parameters from Milan.
//
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    hypre_preconditioner_pt->amg_using_complex_smoothing();
//    hypre_preconditioner_pt->amg_complex_smoother() = 7;
//    //hypre_preconditioner_pt->amg_using_simple_smoothing();
//    //hypre_preconditioner_pt->amg_simple_smoother() = 1;
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    hypre_preconditioner_pt->amg_coarsening() = 1;
//    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  }
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  Preconditioner* set_hypre_for_CLJPGSStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    ///////////////////////////////////////////////////////////////////////////
//    // Parameters from Milan.
//
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    hypre_preconditioner_pt->amg_using_simple_smoothing();
//    hypre_preconditioner_pt->amg_simple_smoother() = 1;
//    //hypre_preconditioner_pt->amg_damping() = 1.0;
//
//    //hypre_preconditioner_pt->amg_using_complex_smoothing();
//    //hypre_preconditioner_pt->amg_complex_smoother() = 7;
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    hypre_preconditioner_pt->amg_coarsening() = 0;
//    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_CLJPJStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    ///////////////////////////////////////////////////////////////////////////
//    // Parameters from Milan.
//
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    hypre_preconditioner_pt->amg_using_simple_smoothing();
//    hypre_preconditioner_pt->amg_simple_smoother() = 0;
//    hypre_preconditioner_pt->amg_damping() = 1.0;
//
//    //hypre_preconditioner_pt->amg_using_complex_smoothing();
//    //hypre_preconditioner_pt->amg_complex_smoother() = 7;
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    hypre_preconditioner_pt->amg_coarsening() = 0;
//    hypre_preconditioner_pt->amg_strength() = 0.668;
//
    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_CLJPPilutStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    ///////////////////////////////////////////////////////////////////////////
//    // Parameters from Milan.
//
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    //hypre_preconditioner_pt->amg_using_simple_smoothing();
//    //hypre_preconditioner_pt->amg_simple_smoother() = 0;
//    //hypre_preconditioner_pt->amg_damping() = 1.0;
//
//    hypre_preconditioner_pt->amg_using_complex_smoothing();
//    hypre_preconditioner_pt->amg_complex_smoother() = 7;
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    hypre_preconditioner_pt->amg_coarsening() = 0;
//    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  //
  Preconditioner* set_hypre_for_RSGSStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    ///////////////////////////////////////////////////////////////////////////
//    // Parameters from Milan.
//
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    hypre_preconditioner_pt->amg_using_simple_smoothing();
//    hypre_preconditioner_pt->amg_simple_smoother() = 1;
//    //hypre_preconditioner_pt->amg_damping() = 1.0;
//
//    //hypre_preconditioner_pt->amg_using_complex_smoothing();
//    //hypre_preconditioner_pt->amg_complex_smoother() = 7;
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    hypre_preconditioner_pt->amg_coarsening() = 1;
//    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_RSJStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    ///////////////////////////////////////////////////////////////////////////
//    // Parameters from Milan.
//
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    hypre_preconditioner_pt->amg_using_simple_smoothing();
//    hypre_preconditioner_pt->amg_simple_smoother() = 0;
//    hypre_preconditioner_pt->amg_damping() = 1.0;
//
//    //hypre_preconditioner_pt->amg_using_complex_smoothing();
//    //hypre_preconditioner_pt->amg_complex_smoother() = 7;
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    hypre_preconditioner_pt->amg_coarsening() = 1;
//    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_RSPilutStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
//    HyprePreconditioner* hypre_preconditioner_pt = 
//      static_cast<HyprePreconditioner*>(another_preconditioner_pt);
//
//    Hypre_default_settings::
//      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);
//
//    ///////////////////////////////////////////////////////////////////////////
//    // Parameters from Milan.
//
//    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
//    /// include:
//    ///  0 = Jacobi 
//    ///  1 = Gauss-Seidel, sequential
//    ///      (very slow in parallel!)
//    /// To use these methods set AMG_using_simple_smoothing to true
//    //hypre_preconditioner_pt->amg_using_simple_smoothing();
//    //hypre_preconditioner_pt->amg_simple_smoother() = 0;
//    //hypre_preconditioner_pt->amg_damping() = 1.0;
//
//    hypre_preconditioner_pt->amg_using_complex_smoothing();
//    hypre_preconditioner_pt->amg_complex_smoother() = 7;
//
//    // AMG coarsening strategy. Coarsening types include:
//    ///  0 = CLJP (parallel coarsening using independent sets)
//    ///  1 = classical RS with no boundary treatment (not recommended
//    ///      in parallel)
//    hypre_preconditioner_pt->amg_coarsening() = 1;
//    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
} // end of namespace Hypre_Subsidiary_Preconditioner_Helper
#endif // End of ifdef OOMPH_HAS_HYPRE

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
namespace AnnularWedgeLagrange
{

  std::map<int,std::string> valid_prob_id_map;

  // Prob id, set by main method
  const int* Prob_id_pt = 0;

  std::string Prob_str = "";
  std::string Ang_deg_str = "";
  std::string Noel_str = "";

  ///////////////////////
  // Domain dimensions.//
  ///////////////////////
  //
  // This is a square domain: x,y \in [0,1]
  //

  // Min and max x value respectively.
  static const double X_min = 0.0;
  static const double X_max = 1.0;

  // Min and max y value respectively.
  static const double Y_min = 0.0;
  static const double Y_max = 1.0;

  // The length in the x and y direction respectively.
  static const double Lx = X_max - X_min;
  static const double Ly = Y_max - Y_min;


  /////////////////////////////////////////////////////////////////////////////

  // CL - set directly from the command line.
  // To set from CL - a CL value is set, this is changed depending on that
  // value.
  //
  // Problem parameter overview:
  //
  // // Solvers:
  //
  //
  // F_ns + L^T inv(W) L | B^T
  // --------------------------
  //                     | W
  //
  // W = 0 (SuperLU)
  // NS_solver = 0 (SuperLU) or 1 (LSC)
  // 
  // If NS_solver = 1, then we have:
  //
  // | F | B^T |   |
  // |----------   |
  // |   |-M_s |   |
  // |-------------|
  // |         | W |
  //
  // F_solver = 0 (SuperLU) or 1 (AMG)
  // P_solver = 0 (SuperLU) or 1 (AMG)
  // 

  // All problems based on the square domain will be in the same file.
  // Each unique problem will have an id.
  // 00 = (SqTmp) Square, custom stuff...
  // 01 = (SqPo) Square, Parallel outflow (para inflow) 
  // 02 = (SqTf) Square, Tangential flow (Semi para inflow)
  // 03 = (SqTfPo) Square, Tangential flow, Parallel outflow (semi para inflow)
  //
  // 10 = (AwTmp) Annulus wedge, custom stuff...
  // 11 = (AwPo) Annulus wedge, Parallel outflow (para inflow)
  // 12 = (AwTf) Annulus wedge, Tangential flow (semi para inflow)
  // 13 = (AwTfPo) Annulus wedge, Tan. flow, Para. outflow (semi para inflow)

  // These are self explanatory:
  double Ang_deg = -1.0; //CL, Angle in degrees
  double Ang = -1.0; //CL, Angle in degrees
  unsigned Noel = 4; //CL, Number of elements in 1D
  double Phi_lo = 0.0; // CL, Low angle of annulus
  double Phi_hi = 90.0; // CL, Upper angle of annulus
  double R_lo = 1.0; // CL, lower radius from the origin
  double R_hi = 3.0; // CL, higher radius from the origin
  double V_lo = 1.0; // CL, higher radius from the origin
  unsigned BC_setting = 0; // CL, higher radius from the origin


  double V_loR_lo = -1.0;

  inline void setup_commandline_flags()
  {
    CommandLineArgs::specify_command_line_flag("--ang", &Ang_deg);

    CommandLineArgs::specify_command_line_flag("--phi_lo", &Phi_lo);
    CommandLineArgs::specify_command_line_flag("--phi_hi", &Phi_hi);
    CommandLineArgs::specify_command_line_flag("--r_lo", &R_lo);
    CommandLineArgs::specify_command_line_flag("--r_hi", &R_hi);

    CommandLineArgs::specify_command_line_flag("--bc", &BC_setting);

    CommandLineArgs::specify_command_line_flag("--noel", &Noel);
  }

  inline void set_prob_str()
  {

    std::map<int,std::string>::iterator prob_id_it;

    // Set a problem id to identify the problem.
    // This is used for book keeping purposes.
    if(CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
    {
      if(Prob_id_pt == 0)
      {
        std::ostringstream err_msg;
        err_msg << "Please set Prob_id_pt from NSPP." << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      // Locally cache the problem id for convenience.
      const int prob_id = *Prob_id_pt;

      // The argument immediately after --prob_id is put into SL::Prob_id.
      // If this begins with "--", then no problem id has been provided.

      // Maybe I should check if SL::Prob_id is a number or a string...

      // We only accept problem IDs as defined below.
      // Creating a set of acceptable IDs
      // 10,11,12,13
      prob_id_it = valid_prob_id_map.find(prob_id);

      // Check if they have provided an acceptable ID.
      if(prob_id_it == valid_prob_id_map.end())
      {
        std::ostringstream err_msg;
        err_msg << "Please provide a problem id to identify the problem after "
          << "after the argument --prob_id.\n" 
          << "Acceptable IDs are:\n"
          << "10 = (AwTmp) Square, custom stuff...\n"
          << "11 = (AwPo) Square, Parallel outflow (para inflow)\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        Prob_str = prob_id_it->second;
      }
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please set Prob_id with: \n"
        << "--prob_id \n"; 
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  } // set_prob_str

  inline void set_ang_str()
  {
    // Check if the problem id has been set. This will be used to see
    // if we must or must not supply the angle.
    if(Prob_id_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "Oh dear, Prob_id_pt is null. Please set this in main().\n"
        << "This should be stored in NSPP::Prob_id, and set by cmd via\n"
        << "--prob_id \n"; 
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    } 

    // This problem requires tilting, thus we set the Ang and Ang_deg_str.
    // But first we check that --ang has been set.
    // Check that Ang has been set.
    if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
    {
      std::ostringstream err_msg;
      err_msg << "Angle has not been set. Set (in degrees) with: \n"
        << "--ang \n"; 
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
      // Now we need to convert Ang_deg into radians.
      Ang = Ang_deg * (MathematicalConstants::Pi / 180.0);

      // Now we set the Ang_deg_str.
      std::ostringstream strs;
      strs << "A" << Ang_deg;
      Ang_deg_str = strs.str();
    }
    else
    {
      Ang = -1.0;
      std::ostringstream strs;
      strs << "A_";
      Ang_deg_str = strs.str();
    }
  } // set_ang_str


    inline void set_noel_str()
  {
    // Set Noel_str, used for book keeping.
    if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
    {
      std::ostringstream strs;
      strs << "N" <<Noel;
      Noel_str = strs.str();
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please supply the number of elements in 1D using --noel.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  }


  inline void generic_setup()
  {
    // Insert the prob id and string pairs.
    valid_prob_id_map.insert(std::pair<int,std::string>(10,"AwTmp"));
    valid_prob_id_map.insert(std::pair<int,std::string>(11,"AwPo"));

    V_loR_lo = V_lo*R_lo;

    set_prob_str();
    set_ang_str();
    set_noel_str();
  }

  inline std::string prob_str()
  {
    set_prob_str();
    return Prob_str;
  }

  inline std::string ang_deg_str()
  {
    set_ang_str();
    return Ang_deg_str;
  }

  inline std::string noel_str()
  {
    set_noel_str();
    return Noel_str;
  }

  inline std::string create_label()
  {
    std::string label = prob_str() + ang_deg_str() + noel_str();
    return label; 
  } // inlined function create_label
} // Namespace AnnularWedgeLagrange

namespace StepLagrange
{

  std::map<int,std::string> valid_prob_id_map;

  // Prob id, set by main method
  const int* Prob_id_pt = 0;

  std::string Prob_str = "";
  std::string Ang_deg_str = "";
  std::string Noel_str = "";


  ///////////////////////
  // Domain dimensions.//
  ///////////////////////
  //
  // This is a square domain: x,y \in [0,1]
  //

  /////////////////////////////////////////////////////////////////////////////

  // CL - set directly from the command line.
  // To set from CL - a CL value is set, this is changed depending on that
  // value.
  //
  // Problem parameter overview:
  //
  // // Solvers:
  //
  //
  // F_ns + L^T inv(W) L | B^T
  // --------------------------
  //                     | W
  //
  // W = 0 (SuperLU)
  // NS_solver = 0 (SuperLU) or 1 (LSC)
  // 
  // If NS_solver = 1, then we have:
  //
  // | F | B^T |   |
  // |----------   |
  // |   |-M_s |   |
  // |-------------|
  // |         | W |
  //
  // F_solver = 0 (SuperLU) or 1 (AMG)
  // P_solver = 0 (SuperLU) or 1 (AMG)
  // 

  // All problems based on the square domain will be in the same file.
  // Each unique problem will have an id.
  // 00 = (SqTmp) Square, custom stuff...
  // 01 = (SqPo) Square, Parallel outflow (para inflow) 
  // 02 = (SqTf) Square, Tangential flow (Semi para inflow)
  // 03 = (SqTfPo) Square, Tangential flow, Parallel outflow (semi para inflow)
  //
  // 10 = (AwTmp) Annulus wedge, custom stuff...
  // 11 = (AwPo) Annulus wedge, Parallel outflow (para inflow)
  // 12 = (AwTf) Annulus wedge, Tangential flow (semi para inflow)
  // 13 = (AwTfPo) Annulus wedge, Tan. flow, Para. outflow (semi para inflow)

  // These are self explanatory:
  double Ang_deg = 30.0; //CL, Angle in degrees
  double Ang = 0.0; //CL, Angle in degrees
  unsigned Noel = 4; //CL, Number of elements in 1D
  // the default is the norm of the momentum block.


  inline void setup_commandline_flags()
  {
    CommandLineArgs::specify_command_line_flag("--ang", &Ang_deg);

    CommandLineArgs::specify_command_line_flag("--noel", &Noel);
  }

  inline void set_prob_str()
  {

    std::map<int,std::string>::iterator prob_id_it;

    // Set a problem id to identify the problem.
    // This is used for book keeping purposes.
    if(CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
    {
      if(Prob_id_pt == 0)
      {
        std::ostringstream err_msg;
        err_msg << "Please set Prob_id_pt from NSPP." << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      // Locally cache the problem id for convenience.
      const int prob_id = *Prob_id_pt;

      // The argument immediately after --prob_id is put into SL::Prob_id.
      // If this begins with "--", then no problem id has been provided.

      // Maybe I should check if SL::Prob_id is a number or a string...

      // We only accept problem IDs as defined below.
      // Creating a set of acceptable IDs
      // 10,11,12,13
      prob_id_it = valid_prob_id_map.find(prob_id);

      // Check if they have provided an acceptable ID.
      if(prob_id_it == valid_prob_id_map.end())
      {
        std::ostringstream err_msg;
        err_msg << "Please provide a problem id to identify the problem after "
          << "after the argument --prob_id.\n" 
          << "Acceptable IDs are:\n"
          << "10 = (StTmp) Square, custom stuff...\n"
          << "11 = (StPo) Square, Parallel outflow (para inflow)\n"
          << "12 = (StTf) Square, Tangential flow (Semi para inflow)\n"
          << "13 = (StTfPo) Square, Tangential flow, Parallel outflow (semi para inflow)\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        Prob_str = prob_id_it->second;
      }
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please set Prob_id with: \n"
        << "--prob_id \n"; 
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  } // set_prob_str

  inline void set_ang_str()
  {
    // Check if the problem id has been set. This will be used to see
    // if we must or must not supply the angle.
    if(Prob_id_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "Oh dear, Prob_id_pt is null. Please set this in main().\n"
        << "This should be stored in NSPP::Prob_id, and set by cmd via\n"
        << "--prob_id \n"; 
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    } 

    // If this is the vanilla problem, we set the angle as -1 and set the
    // string as "A_". This would indicate that no angle is used.
    // Furthermore, we ensure that no --ang is set.
    if(Prob_str.compare("StVa") == 0)
    {
      if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream err_msg;
        err_msg << "prob_id is 88, doing vanilla LSC with no tilt.\n"
          << "But you have set --ang, please do not set this."; 
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      Ang = -1;
      // Now we set the Ang_deg_str.
      std::ostringstream strs;
      strs << "A_";
      Ang_deg_str = strs.str();
    }
    else
    // This problem requires tilting, thus we set the Ang and Ang_deg_str.
    {
      // But first we check that --ang has been set.
      // Check that Ang has been set.
      if(!CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream err_msg;
        err_msg << "Angle has not been set. Set (in degrees) with: \n"
          << "--ang \n"; 
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      // Now we need to convert Ang_deg into radians.
      Ang = Ang_deg * (MathematicalConstants::Pi / 180.0);

      // Now we set the Ang_deg_str.
      std::ostringstream strs;
      strs << "A" << Ang_deg;
      Ang_deg_str = strs.str();
    }
  }


    inline void set_noel_str()
  {
    // Set Noel_str, used for book keeping.
    if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
    {
      std::ostringstream strs;
      strs << "N" <<Noel;
      Noel_str = strs.str();
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please supply the number of elements in 1D using --noel.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  }


  inline void generic_setup()
  {
    // Insert the prob id and string pairs.
    valid_prob_id_map.insert(std::pair<int,std::string>(10,"StTmp"));
    valid_prob_id_map.insert(std::pair<int,std::string>(11,"StPo"));
    valid_prob_id_map.insert(std::pair<int,std::string>(12,"StTf"));
    valid_prob_id_map.insert(std::pair<int,std::string>(13,"StTfPo"));
    valid_prob_id_map.insert(std::pair<int,std::string>(88,"StVa"));

    set_prob_str();
    set_ang_str();
    set_noel_str();
  }

  inline std::string prob_str()
  {
    set_prob_str();
    return Prob_str;
  }

  inline std::string ang_deg_str()
  {
    set_ang_str();
    return Ang_deg_str;
  }

  inline std::string noel_str()
  {
    set_noel_str();
    return Noel_str;
  }

  inline std::string create_label()
  {
    std::string label = prob_str() + ang_deg_str() + noel_str();
    return label; 
  } // inlined function create_label
} // Namespace StepLagrange


///////////////////////////////////////////////////////////////////////////////
///
/// Namepace for all things related to the problem where the domain is a 
/// unit square.
///
///////////////////////////////////////////////////////////////////////////////
namespace SquareLagrange
{
  static int PID_SQ_TMP = 10;
  static int PID_SQ_PO = 11;
  static int PID_SQ_TF = 12;
  static int PID_SQ_TFPO = 13;
  static int PID_SQ_VA = 88;

  std::map<int,std::string> valid_prob_id_map;
  
  // Prob id, set by main method
  const int* Prob_id_pt = 0;

  std::string Prob_str = "";
  std::string Ang_deg_str = "";
  std::string Noel_str = "";


  ///////////////////////
  // Domain dimensions.//
  ///////////////////////
  //
  // This is a square domain: x,y \in [0,1]
  //

  // Min and max x value respectively.
  static const double X_min = 0.0;
  static const double X_max = 1.0;

  // Min and max y value respectively.
  static const double Y_min = 0.0;
  static const double Y_max = 1.0;

  // The length in the x and y direction respectively.
  static const double Lx = X_max - X_min;
  static const double Ly = Y_max - Y_min;

  /////////////////////////////////////////////////////////////////////////////

  // CL - set directly from the command line.
  // To set from CL - a CL value is set, this is changed depending on that
  // value.
  //
  // Problem parameter overview:
  //
  // // Solvers:
  //
  //
  // F_ns + L^T inv(W) L | B^T
  // --------------------------
  //                     | W
  //
  // W = 0 (SuperLU)
  // NS_solver = 0 (SuperLU) or 1 (LSC)
  // 
  // If NS_solver = 1, then we have:
  //
  // | F | B^T |   |
  // |----------   |
  // |   |-M_s |   |
  // |-------------|
  // |         | W |
  //
  // F_solver = 0 (SuperLU) or 1 (AMG)
  // P_solver = 0 (SuperLU) or 1 (AMG)
  // 

  // All problems based on the square domain will be in the same file.
  // Each unique problem will have an id.
  // 00 = (SqTmp) Square, custom stuff...
  // 01 = (SqPo) Square, Parallel outflow (para inflow) 
  // 02 = (SqTf) Square, Tangential flow (Semi para inflow)
  // 03 = (SqTfPo) Square, Tangential flow, Parallel outflow (semi para inflow)
  //
  // 10 = (AwTmp) Annulus wedge, custom stuff...
  // 11 = (AwPo) Annulus wedge, Parallel outflow (para inflow)
  // 12 = (AwTf) Annulus wedge, Tangential flow (semi para inflow)
  // 13 = (AwTfPo) Annulus wedge, Tan. flow, Para. outflow (semi para inflow)

  // These are self explanatory:
  double Ang_deg = 30.0; //CL, Angle in degrees
  double Ang = 0.0; //CL, Angle in degrees
  unsigned Noel = 4; //CL, Number of elements in 1D
  // the default is the norm of the momentum block.

  inline void setup_commandline_flags()
  {
    CommandLineArgs::specify_command_line_flag("--ang", &Ang_deg);

    CommandLineArgs::specify_command_line_flag("--noel", &Noel);
  }

  inline void set_prob_str()
  {
    std::map<int,std::string>::iterator prob_id_it;

    // Set a problem id to identify the problem.
    // This is used for book keeping purposes.
    if(CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
    {
      if(Prob_id_pt == 0)
      {
        std::ostringstream err_msg;
        err_msg << "Please set Prob_id_pt from NSPP." << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      // The argument immediately after --prob_id is put into SL::Prob_id.
      // If this begins with "--", then no problem id has been provided.

      // Maybe I should check if SL::Prob_id is a number or a string...

      // We only accept problem IDs as defined below.
      // Creating a set of acceptable IDs
      const int prob_id = *Prob_id_pt;
      prob_id_it = valid_prob_id_map.find(prob_id);

      // Check if they have provided an acceptable ID.
      // If a new element has been inserted, it means the user has provided an
      // ID not in the set.
      if(prob_id_it == valid_prob_id_map.end())
      {
        std::ostringstream err_msg;
        err_msg << "Please provide a problem id to identify the problem after "
          << "after the argument --prob_id.\n" 
          << "Acceptable IDs are:\n"
          << "10 = (SqTmp) Square, custom stuff...\n"
          << "11 = (SqPo) Square, Parallel outflow (para inflow)\n"
          << "12 = (SqTf) Square, Tangential flow (Semi para inflow)\n"
          << "13 = (SqTfPo) Square, Tangential flow, Parallel outflow (semi para inflow)\n"
          << "\n"
          << "88 = (SqVa) Square, Vanilla LSC\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        Prob_str = prob_id_it->second;
      }
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please set Prob_id with: \n"
        << "--prob_id \n"; 
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  } // set_prob_str

  inline void set_ang_str()
  {
    if(Prob_id_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "Oh dear, Prob_id_pt is null. Please set this in main().\n"
        << "This should be stored in NSPP::Prob_id, and set by cmd via\n"
        << "--prob_id \n"; 
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // If this is the vanilla problem, we set the angle as -1 and set the
    // string as "A_". This would indicate that no angle is used.
    // Furthermore, we ensure that no --ang is set.
    if(Prob_str.compare("SqVa") == 0)
    {
      if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream err_msg;
        err_msg << "prob_id is 88, doing vanilla LSC with no tilt.\n"
          << "But you have set --ang, please do not set this."; 
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      Ang = -1.0;
      // Now we set the Ang_deg_str.
      std::ostringstream strs;
      strs << "A_";
      Ang_deg_str = strs.str();
    }
    else
    // This problem requires tilting, thus we set the Ang and Ang_deg_str.
    {
      // But first we check that --ang has been set.
      // Check that Ang has been set.
      if(!CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream err_msg;
        err_msg << "Angle has not been set. Set (in degrees) with: \n"
          << "--ang \n"; 
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      // Now we need to convert Ang_deg into radians.
      Ang = Ang_deg * (MathematicalConstants::Pi / 180.0);

      // Now we set the Ang_deg_str.
      std::ostringstream strs;
      strs << "A" << Ang_deg;
      Ang_deg_str = strs.str();
    }
  } // set_ang_str

  inline void set_noel_str()
  {
    // Set Noel_str, used for book keeping.
    if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
    {
      std::ostringstream strs;
      strs << "N" <<Noel;
      Noel_str = strs.str();
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please supply the number of elements in 1D using --noel.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
  }

  inline void generic_setup()
  {
    // Insert the prob id and string pairs.
    valid_prob_id_map.insert(std::pair<int,std::string>(PID_SQ_TMP,"SqTmp"));
    valid_prob_id_map.insert(std::pair<int,std::string>(PID_SQ_PO,"SqPo"));
    valid_prob_id_map.insert(std::pair<int,std::string>(PID_SQ_TF,"SqTf"));
    valid_prob_id_map.insert(std::pair<int,std::string>(PID_SQ_TFPO,"SqTfPo"));

    valid_prob_id_map.insert(std::pair<int,std::string>(PID_SQ_VA,"SqVa"));

    set_prob_str();
    set_ang_str();
    set_noel_str();
  }

  inline std::string prob_str()
  {
    set_prob_str();
    return Prob_str;
  }

  inline std::string ang_deg_str()
  {
    set_ang_str();
    return Ang_deg_str;
  }

  inline std::string noel_str()
  {
    set_noel_str();
    return Noel_str;
  }

  inline std::string create_label()
  {
    std::string label = prob_str() + ang_deg_str() + noel_str();
    return label; 
  } // inlined function create_label

} // Namespace SquareLagrange


namespace LagrangianPreconditionerHelpers
{
  //////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////

  // Set my problem constructor.
  Vector<Mesh*> Mesh_pt;

  // Set by main()
  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;

  // Set by main()
  Problem* Problem_pt = 0;

  // Set by main(), from NSPP
  int* Vis_pt = 0;

  std::string* Label_str_pt = 0;

  ///////////////////////////////////////////////////
  // Set directly by CL
  int W_solver = -1;
  int NS_solver = -1;
  int F_solver = -1;
  int P_solver = -1;

  double f_amg_strength = -1.0;
  double f_amg_damping = -1.0;
  int f_amg_coarsening = -1;
  int f_amg_simple_smoother = -1;
  int f_amg_complex_smoother = -1;
  int f_amg_iterations = -1;
  int f_amg_smoother_iterations = -1;

  double p_amg_strength = -1.0;
  double p_amg_damping = -1.0;
  int p_amg_coarsening = -1;
  int p_amg_simple_smoother = -1;
  int p_amg_complex_smoother = -1;
  int p_amg_iterations = -1;
  int p_amg_smoother_iterations = -1;

  std::string Doc_prec_dir_str = "";
  double Scaling_sigma = 0.0;

  // Set indirectly by CL

  // Set to true if --doc_prec is provided with a directory.
  bool Doc_prec = false;

  // Set to true if --print_hypre is set.
  bool Print_hypre = false;

  // Set to false if --sigma value is set.
  bool Use_axnorm = true;

  // Set to true if --bdw is set.
  bool Use_block_diagonal_w = false;

  // Set to true if --lsc_only is set.
  bool Lsc_only = false;


  inline void setup_commandline_flags()
  {
    // Flag to output the preconditioner, used for debugging.
    // string
    CommandLineArgs::specify_command_line_flag(
        "--doc_prec",&Doc_prec_dir_str);

    // No parameter after.
    CommandLineArgs::specify_command_line_flag(
        "--lsc_only");

    // double
    CommandLineArgs::specify_command_line_flag(
        "--sigma",&Scaling_sigma);

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
        "--p_solver",&P_solver);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--f_solver",&F_solver);

    // NS_F block AMG parameters
    // double
    CommandLineArgs::specify_command_line_flag(
        "--f_amg_str",&f_amg_strength);
    // double
    CommandLineArgs::specify_command_line_flag(
        "--f_amg_damp",&f_amg_damping);

    // int
    CommandLineArgs::specify_command_line_flag("--f_amg_coarse",
        &f_amg_coarsening);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--f_amg_sim_smoo",&f_amg_simple_smoother);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--f_amg_com_smoo",&f_amg_complex_smoother);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--f_amg_iter",&f_amg_iterations);

    //int
    CommandLineArgs::specify_command_line_flag("--f_amg_smiter",
        &f_amg_smoother_iterations);

    // NS_P block AMG parameters
    // double
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_str",&p_amg_strength);
    // double
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_damp",&p_amg_damping);
    // int
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_coarse",&p_amg_coarsening);
    // int
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_sim_smoo",&p_amg_simple_smoother);
    // int
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_com_smoo",&p_amg_complex_smoother);

    // int
    CommandLineArgs::specify_command_line_flag(
        "--p_amg_iter",&p_amg_iterations);
    // int
    CommandLineArgs::specify_command_line_flag("--p_amg_smiter",
        &p_amg_smoother_iterations);

    CommandLineArgs::specify_command_line_flag("--print_hypre");
  }

  inline void generic_setup()
  {
    if(CommandLineArgs::command_line_flag_has_been_set("--lsc_only"))
    {
      Lsc_only = true;
    }
    else
    {
      Lsc_only = false;
    }

    // Document the preconditioner? Default is false.
    if(CommandLineArgs::command_line_flag_has_been_set("--doc_prec"))
    {
      // The argument immediately after --doc_prec is put into SL::Doc_prec_dir.
      // If this begins with "--", then no prec directory has been provided.
      std::size_t found = Doc_prec_dir_str.find("--");

      // Check if they have set the doc_prec directory.
      if(found != std::string::npos)
      {
        std::ostringstream err_msg;
        err_msg 
          << "Please provide the doc_prec directory "
          << "after the argument --doc_prec.\n" 
          << "This must not start with \"--\"." << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      else
      {
        Doc_prec = true;
      }
    }
    else
    {
      Doc_prec = false;
    }


    // If we are using SuperLU for the Navier-Stokes block (NS_solver = 0), 
    // the F_solver and P_solver should not be set.
    if((NS_solver == 0) && 
        (CommandLineArgs::command_line_flag_has_been_set("--p_solver") ||
         CommandLineArgs::command_line_flag_has_been_set("--f_solver")))
    {
      std::ostringstream err_msg;
      err_msg << "NS_solver = 0, using SuperLU for the Navier-Stokes block.\n"
        << "But you have either --f_solver or --p_solver.\n"
        << "These should NOT be set, unless you want to use LSC.\n"
        << "In which case, set --ns_solver 1.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
    {
      Use_axnorm = false;
    }
    else
    {
      Use_axnorm = true;
    }

    if(CommandLineArgs::command_line_flag_has_been_set("--bdw"))
    {
      Use_block_diagonal_w = true;
    }
    else
    {
      Use_block_diagonal_w = false;
    }

    if(CommandLineArgs::command_line_flag_has_been_set("--print_hypre"))
    {
      Print_hypre = true;
    }
    else
    {
      Print_hypre = false;
    }
  } // LPH::generic_setup()


  inline Preconditioner* get_lsc_preconditioner()
  {
    if(Lsc_only && (Mesh_pt.size() != 1))
    {
      std::ostringstream err_msg;
      err_msg << "You have chosen to use only the LSC preconditioner\n"
        << "Thus the Vector Mesh_pt must contain exactly one mesh,\n"
        << "the Navier Stokes bulk mesh.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }


    // If ns_solver is 1, this means we want to use LSC.
    // So the F_solver and P_solver must be set.
    if((F_solver == -1) || (P_solver == -1))
    {
      std::ostringstream err_msg;
      err_msg << "Getting LSC preconditioner for NS block.\n"
        << "But --f_solver and --p_solver have not been set.\n"
        << "0 - Exact (SuperLU)\n"
        << "xy - please check the code for more details.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Check that the problem pointer is set (not null).
    // LSC requires a problem pointer.
    if(Problem_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "Please set the Problem_pt variable.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Create the NS LSC preconditioner.
    NavierStokesSchurComplementPreconditioner* ns_preconditioner_pt =
      new NavierStokesSchurComplementPreconditioner(Problem_pt);


    // Give LSC the bulk mesh (Navier-Stokes mesh).
    ns_preconditioner_pt->set_navier_stokes_mesh(Mesh_pt[0]);



    //// Setting the F solver within the NS block
    /////////////////////////////////////////////

    // Preconditioner for the F block:
    Preconditioner* f_preconditioner_pt = 0;

    // f_solver == 0 is default, so do nothing.
    //
    // AMG depends on the Reynolds number so we check that the Reynolds number
    // is set if LSC is switched on - even if we do not want to use AMG...
    // for consistency.
    if(Vis_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "Please set your Vis_pt variable.\n"
        << "E.g. in main method do: LPH::Vis_pt = &NSPP::Vis;"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    if((*Vis_pt) == -1)
    {
      std::ostringstream err_msg;
      err_msg << "Please set your viscuous term in NSPP header.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }


    ///////////////////////////////////////// FFFFFFFFFFFFFFFFFF


    if(F_solver == 11)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_2D_poison_problem();
#endif
    }
    else if(F_solver == 12)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_navier_stokes_momentum_block();
#endif
    }
    else if(F_solver == 13)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_CLJPGSStrn075();
#endif
    }
    else if(F_solver == 14)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_RSGSStrn075();
#endif
    }
    else if(F_solver == 15)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_CLJPPilutStrn075();
#endif
    }
    else if(F_solver == 16)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_RSPilutStrn075();
#endif
    }
    else if(F_solver == 17)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_augmented_momentum_block();
#endif
    }
    else if(F_solver == 81)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_CLJPGSStrn0668();
#endif
    }
    else if(F_solver == 82)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_CLJPJStrn0668();
#endif
    }
    else if(F_solver == 83)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_CLJPPilutStrn0668();
#endif
    }
    else if(F_solver == 84)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_RSGSStrn0668();
#endif
    }
    else if(F_solver == 85)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_RSJStrn0668();
#endif
    }
    else if(F_solver == 86)
    {
#ifdef OOMPH_HAS_HYPRE
      // LSC takes type "Preconditioner".
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        set_hypre_for_RSPilutStrn0668();
#endif
    }
    else if(F_solver == 2)
    {
      //f_preconditioner_pt = new RayBlockDiagonalPreconditioner<CRDoubleMatrix>;
      f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
    }
    else if(F_solver == 3)
    {
      f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
#ifdef OOMPH_HAS_HYPRE
      dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
        (f_preconditioner_pt)->set_subsidiary_preconditioner_function
        (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_for_2D_poison_problem);
#endif
    }
    else if (F_solver == 69)
    {
#ifdef OOMPH_HAS_HYPRE
      // This is what set_defaults_for_2D_poisson_problem() does:
      f_amg_iterations = 1;
      //f_amg_simple_smoother = 1; // GS - commented out since it is reset below
      //f_amg_strength = 0.25; // commented out since it is reset below, depending on the viscous term.
      //f_amg_coarsening = 0; // CLJP - commented out since it is reset below.
      // END OF 2D poisson stuff.

      // Set this to -1 to make sure complex isn't used
      f_amg_complex_smoother = -1;

      // The default smoother iterations is two.
      f_amg_smoother_iterations = 2;

      // Now set my own stuff.
      f_amg_simple_smoother = 1; // GS
      f_amg_coarsening = 1; // RSs.

      // There is no damping with GS, otherwise we set the parameter:
      f_amg_damping = -1.0;

      // Different amg strength for simple/stress divergence for viscous term.
      // This is only set to the below defaults if the strength parameter is
      // not already set.
      if(f_amg_strength < 0)
      {
        const int vis = *Vis_pt;
        if(vis == 0)
        {
          // Simple form
          f_amg_strength = 0.25;
        }
        else if (vis == 1)
        {
          // Stress divergence form
          f_amg_strength = 0.668;
        }
        else
        {
          std::ostringstream err_msg;
          err_msg << "Do not recognise viscuous term: " << vis << std::endl;

          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Setup the preconditioner.
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        get_custom_hypre_preconditioner(
            f_amg_iterations, f_amg_smoother_iterations, 
            f_amg_simple_smoother, f_amg_complex_smoother,
            f_amg_damping, f_amg_strength,
            f_amg_coarsening);

      if(Print_hypre)
      {
        Hypre_Subsidiary_Preconditioner_Helper::print_hypre_settings(
            f_preconditioner_pt);
      }
#endif
    }
    else if (F_solver == 96)
    {
#ifdef OOMPH_HAS_HYPRE
      // In this method we pass everything to the global AMG parameters and
      // let that handle the work.

      // AMG coarsening:
      // Set: RayGlobalAMGParam::amg_coarsening = 
      //
      // 0 - CLJP
      // 1 - Classical RS
      // 3 - modified RS
      // 6 - Falgout - default on documentation
      // 8 - PMIS
      // 10 - HMIS
      // 11 - One pass on RS coarsening on each processor, not recommended.

      // AMG smoother:
      // Set: RayGlobalAMGParam::amg_smoother = 
      // 0 - Jacobi (Need to set damping as well)
      // 1 - Gauss-Seidel, sequential, very slow in parallel
      // 2 - GS - interior parallel, serial on boundary.
      // 3 - hybrid GS or SOR, forward solve
      // 4 - hybrid GS or SOR, backwards solve.
      // 6 - hybrid symmetric GS or SSOR.
      // REMEMBER TO SET SIMPLE SMOOTHING TO TRUE. This should be handled
      // by your hypre preconditioner creation function.
      //
      // Complex smoothing:
      // 6 - Schwarz
      // 7 - Pilut
      // 8 - ParaSails
      // 9 - Euclid
      // REMEMBER TO SET SIMPLE SMOOTHING TO FALSE, this should be handled
      // by the hypre subsidiary helper function.

      // First check that the user has selected a coarsening strategy
      // using --f_amg_coarse.
      // If the user has, then f_amg_coarsening should no longer be -1.
      bool coarsening_ok 
        = RayGen::check_amg_coarsener(f_amg_coarsening);
      if((f_amg_coarsening < 0) || (!coarsening_ok))
      {
        std::ostringstream err_msg;
        err_msg 
          << "Please set a coarsening strategy with --f_amg_coarse.\n"
          << "Current coarsening strategy (f_amg_coarsening) is: " 
          << f_amg_coarsening << "\n\n"

          << "You have either not set it or it is not valid.\n\n"

          << "Valid IDs are:\n"
          << "0 - CLJP\n"
          << "1 - Classical RS\n"
          << "3 - modified RS\n"
          << "6 - Falgout (default in documentation)\n"
          << "8 - PMIS\n"
          << "10 - HMIS\n"
          << "11 - One pass on RS coarsening on each processor, "
          << "not recommended.\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      // First check if the user has selected both simple and complex
      // smoothing, this is not allowed, either one or the other.
      if((f_amg_simple_smoother >= 0) && 
          (f_amg_complex_smoother >= 0))
      {
        std::ostringstream err_msg;
        err_msg 
          << "Both simple and complex smoothers are set.\n"
          << "Please choose one or the other.\n\n" 
          << "f_amg_simple_smoother is " << f_amg_simple_smoother << "\n"
          << "f_amg_complex_smoother is " << f_amg_complex_smoother << "\n\n"

          << "Simple smoother IDs, set with --f_amg_sim_smoo\n"
          << "0 - Jacobi (Need to set damping as well)\n"
          << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
          << "2 - GS - interior parallel, serial on boundary\n"
          << "3 - hybrid GS or SOR, forward solve\n"
          << "4 - hybrid GS or SOR, backwards solve\n"
          << "6 - hybrid symmetric GS or SSOR.\n\n"

          << "Complex smoother IDs, set with --f_amg_com_smoo\n"
          << "6 - Schwarz (default in documentation)\n"
          << "7 - Pilut\n"
          << "8 - ParaSails\n"
          << "9 - Euclid\n" 
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      else if((f_amg_simple_smoother < 0) && (f_amg_complex_smoother < 0))
      {
        std::ostringstream err_msg;
        err_msg 
          << "Please select a smoother for the f block.\n"
          << "Use --f_amg_sim_smoo or --f_amg_com_smoo flag.\n\n"

          << "Simple smoother IDs, set with --f_amg_sim_smoo\n"
          << "0 - Jacobi (Need to set damping as well)\n"
          << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
          << "2 - GS - interior parallel, serial on boundary\n"
          << "3 - hybrid GS or SOR, forward solve\n"
          << "4 - hybrid GS or SOR, backwards solve\n"
          << "6 - hybrid symmetric GS or SSOR.\n\n"

          << "Complex smoother IDs, set with --f_amg_com_smoo\n"
          << "6 - Schwarz (default in documentation)\n"
          << "7 - Pilut\n"
          << "8 - ParaSails\n"
          << "9 - Euclid\n" 
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }


      // Only one of the two smoothers have been set. We see if these are
      // valid smoothing IDs.
      if(f_amg_simple_smoother >= 0)
      {
        // check if simple smoother is okay.
        bool sim_smoother_ok
          = RayGen::check_amg_sim_smoother(f_amg_simple_smoother);

        if(!sim_smoother_ok)
        {
          std::ostringstream err_msg;
          err_msg 
            << "Please provide a valid simple smoother \n"
            << "using --f_amg_sim_smoo. You have provided: " 
            << f_amg_simple_smoother <<"\n\n"
            << "Valid IDs are:\n"
            << "0 - Jacobi (Need to set damping as well)\n"
            << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
            << "2 - GS - interior parallel, serial on boundary\n"
            << "3 - hybrid GS or SOR, forward solve\n"
            << "4 - hybrid GS or SOR, backwards solve\n"
            << "6 - hybrid symmetric GS or SSOR.\n"
            << std::endl;

          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
      }
      else if(f_amg_complex_smoother >=0)
      {
        // check if complex smoother is valid.
        bool com_smoother_ok
          = RayGen::check_amg_com_smoother(f_amg_complex_smoother);
        if(!com_smoother_ok)
        {
          std::ostringstream err_msg;
          err_msg 
            << "Please provide a valid complex smoother \n"
            << "using --f_amg_com_smoo. You have provided: " 
            << f_amg_complex_smoother <<"\n\n"
            << "Valid IDs are:\n"
            << "6 - Schwarz (default in documentation)\n"
            << "7 - Pilut\n"
            << "8 - ParaSails\n"
            << "9 - Euclid\n"
            << std::endl;

          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
      }
      else
      {
        // Should never get here, something has gone wrong.
        std::ostringstream err_msg;
        err_msg 
          << "Something when wrong with your selection of smoother.\n"
          << "Choose to use either a simple smoother or complex smoother\n"
          << "with the flags --f_amg_sim_smoo and --f_amg_com_smoo.\n\n"

          << "Current f_amg_simple_smoother is " 
          << f_amg_simple_smoother << "\n"
          << "Current f_amg_complex_smoother is " 
          << f_amg_complex_smoother << "\n\n"

          << "Simple smoother IDs, set with --f_amg_sim_smoo\n"
          << "0 - Jacobi (Need to set damping as well)\n"
          << "1 - Gauss-Seidel, sequential, very slow in parallel\n"
          << "2 - GS - interior parallel, serial on boundary\n"
          << "3 - hybrid GS or SOR, forward solve\n"
          << "4 - hybrid GS or SOR, backwards solve\n"
          << "6 - hybrid symmetric GS or SSOR.\n\n"

          << "Complex smoother IDs, set with --f_amg_com_smoo\n"
          << "6 - Schwarz (default in documentation)\n"
          << "7 - Pilut\n"
          << "8 - ParaSails\n"
          << "9 - Euclid\n" 
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }


      // Damping is required for Jacobi or hybrid SOR smoothers.
      // First check if amg_simple_smoother is one of those.
      bool damping_required
        = RayGen::check_amg_sim_smoother_damping_required(f_amg_simple_smoother);

      if(damping_required && (f_amg_damping < 0))
      {
        std::ostringstream err_msg;
        err_msg 
          << "Have selected simple smoother: " << f_amg_simple_smoother << "\n"
          << "Damping parameter is required for this smoother.\n"
          << "Please set it via --f_amg_damp.\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
      // Note that we haven't checked the reverse, i.e. if a smoother
      // which does not required damping is set, and the damping parameter
      // is set, we allow this... if implemented properly, the hypre
      // code should just ignore it. It will also be interesting to find
      // smoothers which the damping parameter affects the result.


      // Now check if the strength parameter is set. 
      // But don't throw an error, just a warning.
      if(f_amg_strength < 0)
      {
        std::ostringstream warning_msg;
        warning_msg 
          << "You have not set the strength parameter for --f_amg_str\n"
          << "The default settings will be used.\n"
          << "0.25 for simple form of the viscous term or \n"
          << "0.668 for stress divergence form of the viscous term.\n"
          << std::endl;

        throw OomphLibWarning(warning_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);

        const int vis = (*Vis_pt);
        if(vis == 0)
        {
          // Simple form
          f_amg_strength = 0.25;
        }
        else if(vis == 1)
        {
          // Stress divergence form
          f_amg_strength = 0.668;
        }
        else
        {
          std::ostringstream err_msg;
          err_msg 
            << "Do not recognise viscuous term: " << vis << "\n"
            << "Please point it to NSPP::Vis\n" 
            << std::endl;

          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
      }

      // Now check the f_amg_iterations
      if(f_amg_iterations < 0)
      {
        std::ostringstream err_msg;
        err_msg 
          << "Have not set f_amg_iterations via --f_amg_iter\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      if(f_amg_smoother_iterations < 0)
      {
        std::ostringstream err_msg;
        err_msg 
          << "Have not set f_amg_smoother_iterations via --f_amg_smiter\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }

      // Setup the preconditioner.
      f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        get_custom_hypre_preconditioner(
            f_amg_iterations, f_amg_smoother_iterations, 
            f_amg_simple_smoother, f_amg_complex_smoother,
            f_amg_damping, f_amg_strength,
            f_amg_coarsening);

      if(Print_hypre)
      {
        Hypre_Subsidiary_Preconditioner_Helper::print_hypre_settings(
            f_preconditioner_pt);
      }
#endif
    }

    // Set the preconditioner in the LSC preconditioner.
    ns_preconditioner_pt->set_f_preconditioner(f_preconditioner_pt);

    ///////////////////////////////////////// FFFFFFFFFFFFFFFFFF


    // P block solve.
    ///////////////////////////////////////////

    // Pointer to the preconditioner.
    Preconditioner * p_preconditioner_pt = 0;

    //SL::P_solver == 0 is default, so do nothing.
    if(P_solver == 1) 
    {
#ifdef OOMPH_HAS_HYPRE

      p_preconditioner_pt = new HyprePreconditioner;

      // Cast it to a Hypre preconditioner so we can set AMG settings.
      HyprePreconditioner* hypre_preconditioner_pt =
        static_cast<HyprePreconditioner*>(p_preconditioner_pt);

      Hypre_default_settings::
        set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

      if(Print_hypre)
      {
        Hypre_Subsidiary_Preconditioner_Helper::print_hypre_settings(
            p_preconditioner_pt);
      }

      // Set it as the p preconditioner for LSC
      //     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif
    }
    else if(P_solver == 96)
    {
#ifdef OOMPH_HAS_HYPRE

      p_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
        get_custom_hypre_preconditioner(
            p_amg_iterations, p_amg_smoother_iterations, 
            p_amg_simple_smoother, p_amg_complex_smoother,
            p_amg_damping, p_amg_strength,
            p_amg_coarsening);

      if(Print_hypre)
      {
        Hypre_Subsidiary_Preconditioner_Helper::print_hypre_settings(
            p_preconditioner_pt);
      }
#endif
    }
    else if(P_solver == 2)
    {
#ifdef OOMPH_HAS_HYPRE
      p_preconditioner_pt = new HyprePreconditioner;

      HyprePreconditioner* hypre_preconditioner_pt =
        static_cast<HyprePreconditioner*>(p_preconditioner_pt);

      hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

      // Setup v-cycles
      hypre_preconditioner_pt->set_amg_iterations(2);
      hypre_preconditioner_pt->amg_smoother_iterations() = 2;

      // Setup smoother
      // simple: 0 - DJ, 1 - GS
      // compelx: Pilut - 7
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother() = 0;
      // only applicable for DJ
      hypre_preconditioner_pt->amg_damping() = 0.8;

      // Setup coarsening
      // 0 - CLJP
      // 1 - RS
      hypre_preconditioner_pt->amg_coarsening() = 1;
#endif
    }

    ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);


    return ns_preconditioner_pt;



  } // EoFunc get_lsc_preconditioner(...)

  inline Preconditioner* get_lgr_preconditioner()
  {

    LagrangeEnforcedflowPreconditioner* prec_pt
      = new LagrangeEnforcedflowPreconditioner;

   // Set the mesh
    if(Mesh_pt.size() < 2)
    {
      std::ostringstream err_msg;
      err_msg << "There must be at least two meshes.\n"
        << "Since the mesh size is 1, did you mean to use lsc only?\n"
        << "If so, set --lsc_only and --f_solver --p_solver\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    prec_pt->set_meshes(Mesh_pt);

    // Set W solver.
    if(W_solver == -1)
    {
      std::ostringstream err_msg;
      err_msg << "There W_solver has not been set.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else if(W_solver == 0)
    {
      // Using SuperLU, this is the default, do nothing.
    }
    else if (W_solver ==1)
    {
      prec_pt->set_lagrange_multiplier_subsidiary_preconditioner
        (Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper
         ::get_lagrange_multiplier_preconditioner);
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

    // The preconditioner for the fluid block:
    if(NS_solver == -1)
    {
      std::ostringstream err_msg;
      err_msg << "The NS solver has not been set.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else if(NS_solver == 0) // Exact solve.
    {
      // This is the default, do nothing.
      // But the param's F_solver and P_solver should not have been set,
      // i.e. it should stay as -1.

      if((F_solver != -1) || (P_solver != -1))
      {
        std::ostringstream err_msg;
        err_msg << "Doing exact NS solve. (NS_solver is 0)\n"
          << "but you have set F_solver and P_solver as well."
          << "Please leave these as -1.\n"
          << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
    }
    else if(NS_solver == 1) // LSC
    {


      // Set the NS preconditioner as LSC.
      prec_pt->set_navier_stokes_lsc_preconditioner(get_lsc_preconditioner());


     } // if for using LSC as NS prec.
    else
    {
      pause("There is no solver for NS.");
    }

    if(!Use_axnorm)
    {
      prec_pt->scaling_sigma() = Scaling_sigma;
    }

    // Set the doc info for book keeping purposes. First we check that it is
    // actually set.
    if(Doc_linear_solver_info_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "Please set Doc_linear_solver_info_pt\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    prec_pt->set_doc_linear_solver_info_pt(Doc_linear_solver_info_pt);

    if(Use_block_diagonal_w)
    {
      prec_pt->use_block_diagonal_w_block();
    }
    else
    {
      prec_pt->use_diagonal_w_block();
    }

    if(Doc_prec)
    {
      prec_pt->enable_doc_prec();
    }

    //     Set the label, use to output information from the preconditioner, such
    //     as the block matrices and the rhs vector
    if(Label_str_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "Please set Label_str_pt, this should point\n"
        << "to the Label_str in NSPP namespace." << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    prec_pt->set_label_pt(Label_str_pt);
    prec_pt->set_doc_prec_directory_pt(&Doc_prec_dir_str);
    return prec_pt;
  }

  inline Preconditioner* get_preconditioner()
  {
    if(Lsc_only)
    {
      return get_lsc_preconditioner();
    }
    else
    {
      return get_lgr_preconditioner();
    }
  } // LPH::setup_preconditioner


  inline std::string create_lsc_label()
  {
    std::string f_str = "";
    std::string p_str = "";
    // Now we continue with setting the string for the solvers.
    // Only set the f_str if NS_solver > 0
    if(NS_solver == 1 || Lsc_only)
    {
      switch(F_solver)
      {
        case 0:
          f_str = "Fe";
          break;
        case 69:
          f_str = "Fa";
          break;
        case 96:
          f_str = "Fray";
          break;
        case 11:
          f_str = "Fh2dp";
          break;
        case 12:
          f_str = "Fhns";
          break;
        case 13:
          f_str = "CLJPGSStrn075";
          break;
        case 14:
          f_str = "FRSGSStrn075";
          break;
        case 15:
          f_str = "FCLJPPilutStrn075";
          break;
        case 16:
          f_str = "FRSPilutStrn075";
          break;
        case 17:
          f_str = "Fray_old"; // I have no short hand for this...
          break;
        case 81:
          f_str = "CLJPGSStrn0668";
          break;
        case 82:
          f_str = "CLJPJStrn0668";
          break;
        case 83:
          f_str = "CLJPPilutStrn0668";
          break;
        case 84:
          f_str = "RSGSStrn0668";
          break;
        case 85:
          f_str = "RSJStrn0668";
          break;
        case 86:
          f_str = "RSPilutStrn0668";
          break;
        case 2:
          f_str = "Fde";
          break;
        case 3:
          f_str = "Fda";
          break;
        default:
          {
            std::ostringstream err_msg;
            err_msg << "There is an unrecognised F_solver, recognised F_solver:\n"
              << "Look at rayheader.h\n"
              << std::endl;
            throw OomphLibError(err_msg.str(),
                OOMPH_CURRENT_FUNCTION,
                OOMPH_EXCEPTION_LOCATION);
          }
      }  // switch for f_solver

      switch(P_solver)
      {
        case 0:
          p_str = "Pe";
          break;
        case 1:
          p_str = "Pa";
          break;
        case 96:
          p_str = "Pray";
          break;
        default:
          {
            std::ostringstream err_msg;
            err_msg << "There is an unrecognised P_solver, recognised P_solver:\n"
              << "Look at rayheader.h\n"
              << std::endl;
            throw OomphLibError(err_msg.str(),
                OOMPH_CURRENT_FUNCTION,
                OOMPH_EXCEPTION_LOCATION);
          }
      } // switch for p_solver
    } // if ns_solve > 0


    std::string prec_str = f_str + p_str;
    return prec_str;
  }

  inline std::string create_lgr_label()
  {
    std::string w_str = "";

    std::string ns_str = "";

    std::string sigma_str = "";

    // Set the string for W_solver.
    switch(W_solver)
    {
      case 0:
        w_str = "We";
        break;
      case 1:
        w_str = "Wc";
        break;
      default:
        {
          std::ostringstream err_msg;
          err_msg << "There is an unrecognised W_solver,\n"
            << "recognised W_solver:\n"
            << "0 = (We) SuperLU solve\n"
            << std::endl;
          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
    } // switch

    if(Use_block_diagonal_w)
    {
      w_str += "bd";
    }
    else
    {
      w_str += "d";
    }

    // Set the string for NS_solver
    switch(NS_solver)
    {
      case 0:
        {
          ns_str = "Ne";
        }
        break;
      case 1:
        ns_str = "Nl";
        break;
      default:
        {
          std::ostringstream err_msg;
          err_msg << "There is an unrecognised NS_solver.\n"
            << "Recognised NS_solver:\n"
            << "0 = (Ne) SuperLU\n"
            << "1 = (Nl) LSC preconditioner\n"
            << std::endl;
          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
    } // switch NS_solver

    if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
    {
      std::ostringstream strs;
      strs << "S" << Scaling_sigma;
      sigma_str = strs.str();
    }

    std::string prec_str = w_str + ns_str + create_lsc_label() + sigma_str;
    return prec_str;
  }

  inline std::string create_label()
  {
    std::string prec_str;
    if(Lsc_only)
    {
      prec_str = create_lsc_label();
    }
    else
    {
      prec_str = create_lgr_label();
    }

    return prec_str;
  } // LPH::create_label()

} // end of namespace LagrangianPreconditionerHelpers




//=============================================================================
/// Namespace to hold variables common to all Navier Stokes problems.
/// Contains the
//=============================================================================
namespace NavierStokesProblemParameters
{
  // From Commandline
  int Prob_id = -1;
  bool Distribute_problem = false;
  int Vis = -1;
  double Rey = -1.0;

  double Rey_start = -1.0;
  double Rey_incre = -1.0;
  double Rey_end = -1.0;

  int Max_solver_iteration = -1;


  bool Using_trilinos_solver = false;

  std::string Soln_dir_str = "";
  std::string Itstime_dir_str = "";

  // Inferred from commandline:
  bool Doc_soln = false;
  std::string Label_str = "";

  // From code:
  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;


  inline void setup_commandline_flags()
  {
    CommandLineArgs::specify_command_line_flag("--dist_prob");

    // A problem ID, there are eight different types of problems.
    // Check the header file.
    CommandLineArgs::specify_command_line_flag("--prob_id",&Prob_id);

    // Flag to output the solution.
    CommandLineArgs::specify_command_line_flag("--doc_soln", 
        &Soln_dir_str);
    CommandLineArgs::specify_command_line_flag("--visc", 
        &Vis);
    CommandLineArgs::specify_command_line_flag("--rey", &Rey);
    CommandLineArgs::specify_command_line_flag("--rey_start", &Rey_start);
    CommandLineArgs::specify_command_line_flag("--rey_incre", &Rey_incre);
    CommandLineArgs::specify_command_line_flag("--rey_end", &Rey_end);
    CommandLineArgs::specify_command_line_flag("--max_solver_iter", 
                                               &Max_solver_iteration);
    // Iteration count and times directory.
    CommandLineArgs::specify_command_line_flag("--itstimedir", 
        &Itstime_dir_str);

    CommandLineArgs::specify_command_line_flag("--trilinos_solver");
  }

  inline void generic_problem_setup(const unsigned& dim)
  {
    // Do we have to distribute the problem?
    if(CommandLineArgs::command_line_flag_has_been_set("--dist_prob"))
    {
      Distribute_problem = true;
    }
    else
    {
      Distribute_problem = false;
    }

    if(!CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
    {
      std::ostringstream err_msg;
      err_msg << "Please set --prob_id." << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Document the solution? Default is false.
    Doc_soln = false;
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
        Doc_soln = true;
      }
    }


    // Set the viscuous term.
    // Default: 0, Sim
    if(CommandLineArgs::command_line_flag_has_been_set("--visc"))
    {
      if (Vis == 0)
      {
        if(dim == 2)
        {
          for (unsigned d = 0; d < 2; d++) 
          {
            NavierStokesEquations<2>::Gamma[d] = 0.0;
          }
        }
        else
        {
          for (unsigned d = 0; d < 3; d++) 
          {
            NavierStokesEquations<3>::Gamma[d] = 0.0;
          }
        }
      }
      else if (Vis == 1)
      {
        if(dim == 2)
        {
          for (unsigned d = 0; d < 2; d++) 
          {
            NavierStokesEquations<2>::Gamma[d] = 1.0;
          }
        }
        else
        {
          for (unsigned d = 0; d < 3; d++) 
          {
            NavierStokesEquations<3>::Gamma[d] = 1.0;
          }
        }
      }
      else
      {
        std::ostringstream err_msg;
        err_msg << "Do not recognise viscuous term: " << Vis << ".\n"
          << "Vis = 0 for simple form\n"
          << "Vis = 1 for stress divergence form\n"
          << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please set --visc to either 0 or 1.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Check if the Reynolds numbers have been set.
    if(  CommandLineArgs::command_line_flag_has_been_set("--rey_start")
        &&CommandLineArgs::command_line_flag_has_been_set("--rey_incre")
        &&CommandLineArgs::command_line_flag_has_been_set("--rey_end")
        &&CommandLineArgs::command_line_flag_has_been_set("--rey"))
    {
      std::ostringstream err_msg;
      err_msg << "You have set all --rey* argument, please choose carefully!\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else if(  CommandLineArgs::command_line_flag_has_been_set("--rey_start")
        &&CommandLineArgs::command_line_flag_has_been_set("--rey_incre")
        &&CommandLineArgs::command_line_flag_has_been_set("--rey_end"))
    {
      std::cout << "Looping Reynolds: \n"
        << "Rey_start = " << Rey_start << std::endl; 
      std::cout << "Rey_incre = " << Rey_incre << std::endl; 
      std::cout << "Rey_end = " << Rey_end << std::endl; 
    }
    else if(!CommandLineArgs::command_line_flag_has_been_set("--rey"))
    {
      std::ostringstream err_msg;
      err_msg << "No Reynolds numbers have been set.\n"
        << "For a single Reynolds number, use --rey.\n"
        << "For looping through Reynolds numbers, use:\n"
        << "--rey_start --rey_incre --rey_end.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

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
    }

    if(CommandLineArgs::command_line_flag_has_been_set("--trilinos_solver"))
    {
#ifdef OOMPH_HAS_TRILINOS
      Using_trilinos_solver = true;
#else
      {
        std::ostringstream err_msg;
        err_msg << "The flag --trilinos_solver has been set.\n"
          << "But OOMPH-LIB does not have trilinos!" << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
#endif
    }
    else
    {
      Using_trilinos_solver = false;
    }

    if(!CommandLineArgs::command_line_flag_has_been_set("--max_solver_iter"))
    {
      std::ostringstream err_msg;
      err_msg << "Please set --max_solver_iter." << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

  } // NSPP::generic_problem_setup()

  inline std::string create_label()
  {
    std::string vis_str = "";
    std::string rey_str = "";

    if(Rey >= 0)
    {
      std::ostringstream strs;
      strs << "R" << Rey;
      rey_str = strs.str();
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Something has gone wrong, the Reynolds number is negative\n"
        << "Rey = " << Rey << "\n"
        << "Please set it again using --rey.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }


    // Set the string for viscuous term.
    if (Vis == 0)
    {
      vis_str = "Sim";
    }
    else if (Vis == 1)
    {
      vis_str = "Str";
    } // else - setting viscuous term.
    else
    {
      std::ostringstream err_msg;
      err_msg << "Do not recognise viscuous term: " << Vis << ".\n"
        << "Vis = 0 for simple form\n"
        << "Vis = 1 for stress divergence form\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    std::string label = vis_str + rey_str;

    return label;
  } // create_label()
} // NavierStokesProblemParameters



//=============================================================================
/// Namespace to set up stuff common to all problems.
//=============================================================================
namespace GenericProblemSetup
{
  inline void setup_solver(const int& max_solver_iter,
                           const double& solver_tol, const double& newton_tol,
                           const bool& use_trilinos_solver,
                           Problem* problem_pt, Preconditioner* prec_pt)
  {
    IterativeLinearSolver* solver_pt = 0;

#ifdef OOMPH_HAS_TRILINOS
    if(use_trilinos_solver)
    {
      TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
      trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;
      solver_pt = trilinos_solver_pt;
    }
    else
    {
      solver_pt = new GMRES<CRDoubleMatrix>;
      // We use RHS preconditioning. Note that by default,
      // left hand preconditioning is used.
      static_cast<GMRES<CRDoubleMatrix>*>(solver_pt)->set_preconditioner_RHS();
    }
#else
    if(use_trilinos_solver)
    {
      std::ostringstream err_msg;
      err_msg << "You have set --trilinos_solver\n"
        << "But OOMPH does not have trilinos" << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      solver_pt = new GMRES<CRDoubleMatrix>;
      // We use RHS preconditioning. Note that by default,
      // left hand preconditioning is used.
      static_cast<GMRES<CRDoubleMatrix>*>(solver_pt)->set_preconditioner_RHS();
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

    solver_pt->tolerance() = solver_tol;
    solver_pt->max_iter() = max_solver_iter;
    problem_pt->newton_solver_tolerance() = newton_tol;

    solver_pt->preconditioner_pt() = prec_pt;
    problem_pt->linear_solver_pt() = solver_pt;
  }
} // namespace GenericProblemSetup


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
      ((unsigned(average_its*10))%10)?
        (*results_stream_pt) << "\t"<< std::fixed << std::setprecision(1)
        << average_its << "(" << nnewtonstep << ")" << "\n":
        (*results_stream_pt) << "\t"<< average_its << "(" << nnewtonstep << ")" << "\n";
      (*results_stream_pt) << std::setprecision(std::cout.precision());
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


} // namespace ResultsFormat
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
#endif // End of ifndef RAY_LAGRANGE_HEADER

