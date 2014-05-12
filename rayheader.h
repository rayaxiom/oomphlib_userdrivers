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
bool check_if_in_set(T myarray[], unsigned nval, T number)
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

//=============================================================================
/// Namespace to hold temporary AMG parameters for customised AMG settings
/// to be used later on. Functions which modifies/uses these variables are:
///   reset_amg_param()
///   set_hypre_for_2D_poison_problem() (might modify str)
///   set_hypre_for_navier_stokes_momentum_block() (might modify str and dmp)
///   set_hypre_ray(...) (might modify everything)
///   
///   There are others... might need to clean this up.
//=============================================================================
namespace RayGlobalAMGParam
{
  double amg_strength = -1.0;
  double amg_damping = -1.0;
  int amg_coarsening = -1;
  int amg_smoother = -1;

  int amg_iterations = -1;
  int amg_smoother_iterations = -1;
  bool print_hypre = false;
}

#ifdef OOMPH_HAS_HYPRE
//=============================================================================
/// helper method for the block diagonal F block preconditioner to allow 
/// hypre to be used for as a subsidiary block preconditioner
//=============================================================================
namespace Hypre_Subsidiary_Preconditioner_Helper
{
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

  //===========================================================================
  /// Reset the variables in RayGlobalAMGParam
  //===========================================================================
  void reset_amg_param()
  {
    RayGlobalAMGParam::amg_strength = -1.0;
    RayGlobalAMGParam::amg_damping = -1.0;
    RayGlobalAMGParam::amg_coarsening = -1;
    RayGlobalAMGParam::amg_smoother = -1;
    RayGlobalAMGParam::amg_iterations = -1;
    RayGlobalAMGParam::amg_smoother_iterations = -1;
    RayGlobalAMGParam::print_hypre = false;
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
  ///        first independent set)
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
    HyprePreconditioner* hypre_preconditioner_pt = 
      checked_static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    // 2D Poisson problem defaults as defined above.
    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    // If the strength parameter has been set (>=0), then we set it.
    if(RayGlobalAMGParam::amg_strength >= 0)
    {
      hypre_preconditioner_pt->amg_strength() 
        = RayGlobalAMGParam::amg_strength;
    }

    // reset the amg parameters.
    reset_amg_param();

    // Print to confirm?
    if(RayGlobalAMGParam::print_hypre == true)
    {
      print_hypre_settings(hypre_preconditioner_pt);
    }

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
    HyprePreconditioner* hypre_preconditioner_pt = 
      checked_static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    // Defaults for Navier Stokes momentum block as defined above.
    Hypre_default_settings::
      set_defaults_for_navier_stokes_momentum_block(hypre_preconditioner_pt);

    // Have the strength parameter been set?
    if(RayGlobalAMGParam::amg_strength >= 0)
    {
      hypre_preconditioner_pt->amg_strength() 
        = RayGlobalAMGParam::amg_strength;
    }

    // Has the damping parameter been set?
    if(RayGlobalAMGParam::amg_damping >= 0)
    {
      hypre_preconditioner_pt->amg_damping() = RayGlobalAMGParam::amg_damping;
    }

    // Reset the amg parameters.
    reset_amg_param();

    // Print to confirm?
    if(RayGlobalAMGParam::print_hypre == true)
    {
      print_hypre_settings(hypre_preconditioner_pt);
    }

    // Return the newly created preconditioner.
    return another_preconditioner_pt;
  } // set_hypre_for_navier_stokes_momentum_block()

  //===========================================================================
  /// Preconditioner* set_hypre_ray()
  ///
  /// Returns a hypre preconditioner where all the settings must be set.
  ///
  //===========================================================================
  Preconditioner* set_hypre_ray()
  {
    // Create a new hypre preconditioner
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      checked_static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    // Set the hypre_method to BoomerAMG. This is hard coded.
    hypre_preconditioner_pt->hypre_method() = HyprePreconditioner::BoomerAMG;

    // Set the amg_iterations.
    // This is usually set to 1
    if(RayGlobalAMGParam::amg_iterations == -1)
    {
      std::ostringstream err_msg;
      err_msg << "The number of AMG iterations must be set." << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      hypre_preconditioner_pt
        ->set_amg_iterations(RayGlobalAMGParam::amg_iterations);
    }


    // Set the number of amg smoother iterations.
    // This is usually set to 2
    if(RayGlobalAMGParam::amg_smoother_iterations == -1)
    {
      std::ostringstream err_msg;
      err_msg << "The number of smoother iterations must be set." << std::endl;

      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      hypre_preconditioner_pt
        ->amg_smoother_iterations() 
        = RayGlobalAMGParam::amg_smoother_iterations;
    }


    // Store this information.
    std::stringstream cyclestream;
    cyclestream << "its"<<RayGlobalAMGParam::amg_iterations
      << "smits" << RayGlobalAMGParam::amg_smoother_iterations;

    string cycle_str = cyclestream.str();
    string smoother_str = "";
    string damping_str = "";
    string coarsening_str = "";
    string strength_str = "";

    ///////////////////////////////////////////////////////////////////////////
    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    // Setting the smoother:
    if(RayGlobalAMGParam::amg_smoother == 1)
    {
      smoother_str = "GS";
      // Setting up Gauss-Seidel
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother() 
        = RayGlobalAMGParam::amg_smoother;
    }
    else if(RayGlobalAMGParam::amg_smoother == 0)
    {
      smoother_str = "J";
      hypre_preconditioner_pt->amg_damping() = RayGlobalAMGParam::amg_damping;

      // Setting up Jacobi with damping.
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother() 
        = RayGlobalAMGParam::amg_smoother;
      if(RayGlobalAMGParam::amg_damping >= 0)
      {
        std::ostringstream strs;
        strs << "Dmp" << hypre_preconditioner_pt->amg_damping();
        damping_str = strs.str(); 
      }
      else
      {
        std::ostringstream err_msg;
        err_msg << "Please set the damping parameter for Jacobi" << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
    }
    else if(RayGlobalAMGParam::amg_smoother == 2)
    {
      smoother_str = "Pilut";
      hypre_preconditioner_pt->amg_using_complex_smoothing();
      hypre_preconditioner_pt->amg_complex_smoother() = 7;
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "You have supplied smoother: " 
        << RayGlobalAMGParam::amg_smoother << "\n";

      err_msg << "No such smoother. 0 = Jacobi, 1 = GS, 2 = Pilut\n";

      err_msg << "Please set RayGlobalAMGParam::smoother\n";
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    if(RayGlobalAMGParam::amg_coarsening == 0)
    {
      coarsening_str = "CLJP";
      hypre_preconditioner_pt->amg_coarsening() 
        = RayGlobalAMGParam::amg_coarsening;
    }
    else if(RayGlobalAMGParam::amg_coarsening == 1)
    {
      coarsening_str = "RS";
      hypre_preconditioner_pt->amg_coarsening() 
        = RayGlobalAMGParam::amg_coarsening;
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "There is no such coarsening: " 
        << RayGlobalAMGParam::amg_coarsening << "\n";
      err_msg << "0 - CLJP, 1 - RS\n";

      err_msg << "Please set RayGlobalAMGParam::amg_coarsening\n";
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Now set the AMG strength parameter.
    if(RayGlobalAMGParam::amg_strength >= 0)
    {
      hypre_preconditioner_pt->amg_strength() 
        = RayGlobalAMGParam::amg_strength;
      std::ostringstream strs;
      strs << "Strn" << hypre_preconditioner_pt->amg_strength();
      strength_str = strs.str(); 
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "Please set RayGlobalAMGParam::amg_strength\n"; 
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Reset the global AMG parameters
    reset_amg_param();

    // Print out the settings?
    if(RayGlobalAMGParam::print_hypre)
    {
      print_hypre_settings(another_preconditioner_pt);
    }

    std::cout << "RAYHYPRE: " << cycle_str
      << coarsening_str 
      << smoother_str
      << damping_str
      << strength_str
      << std::endl;

    return another_preconditioner_pt;
  } // set_hypre_ray

  ///////////////////////////////////////////////
  ///////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  Preconditioner* set_hypre_using_2D_poisson_base()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    string smoother_str = "";
    string damping_str = "";
    string coarsening_str = "";
    string strength_str = "";
    ///////////////////////////////////////////////////////////////////////////
    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    // Setting the smoother:
    if(RayGlobalAMGParam::amg_smoother == 0)
    {
      smoother_str = "GS";
      // Setting up Gauss-Seidel
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother() = 1;
    }
    else if(RayGlobalAMGParam::amg_smoother == 1)
    {
      smoother_str = "J";
      hypre_preconditioner_pt->amg_damping() = RayGlobalAMGParam::amg_damping;

      // Setting up Jacobi with damping.
      hypre_preconditioner_pt->amg_using_simple_smoothing();
      hypre_preconditioner_pt->amg_simple_smoother() = 0;
      if(!(RayGlobalAMGParam::amg_damping < 0))
      {
        std::ostringstream strs;
        strs << "Dmp" << hypre_preconditioner_pt->amg_damping();
        damping_str = strs.str(); 
      }
      else
      {
        std::cout << "Please set your damping using --amg_damping" << std::endl; 
        pause("Please do not continue."); 
      }
    }
    else if(RayGlobalAMGParam::amg_smoother == 2)
    {
      smoother_str = "Pilut";
      hypre_preconditioner_pt->amg_using_complex_smoothing();
      hypre_preconditioner_pt->amg_complex_smoother() = 7;
    }
    else
    {
      std::cout << "You supplied smoother: " << RayGlobalAMGParam::amg_smoother << std::endl;
      std::cout << "No such smoother. 0 is GS, 1 is Jacobi, 2 is Pilut" << std::endl;
      std::cout << "Please set your smoother using --smoother" << std::endl;
      pause("Please do not continue.");
    }


    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    if(RayGlobalAMGParam::amg_coarsening == 0)
    {
      coarsening_str = "CLJP";
      hypre_preconditioner_pt->amg_coarsening() = 0;
    }
    else if(RayGlobalAMGParam::amg_coarsening == 1)
    {
      coarsening_str = "RS";
      hypre_preconditioner_pt->amg_coarsening() = 1;
    }
    else
    {
      std::cout << "There is no such coarsening: " << RayGlobalAMGParam::amg_coarsening << std::endl;
      std::cout << "0 - CLJP, 1 - RS, use --amg_coarsening" << std::endl;
      pause("Do not continue"); 

    }

    if(!(RayGlobalAMGParam::amg_strength < 0))
    {
      hypre_preconditioner_pt->amg_strength() = RayGlobalAMGParam::amg_strength;
      std::ostringstream strs;
      strs << "Strn" << hypre_preconditioner_pt->amg_strength();
      strength_str = strs.str(); 
    }
    else
    {
      std::cout << "Please set the amg_strengh using --amg_strength" << std::endl;
      pause("Do not continue");
    }

    std::cout << "RAYHYPRE: " << coarsening_str 
      << smoother_str
      << damping_str
      << strength_str
      << std::endl;

    return another_preconditioner_pt;
  }
  /////////////////////////////////////////////////////////////////////////////

  Preconditioner* set_hypre_for_augmented_momentum_block()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_navier_stokes_momentum_block(hypre_preconditioner_pt);

    hypre_preconditioner_pt->amg_strength() = 0.668;

    hypre_preconditioner_pt->amg_damping() = 1.0;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_CLJPGSStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 1;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 0;
    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_RSGSStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 1;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 1;
    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_CLJPPilutStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_complex_smoothing();
    hypre_preconditioner_pt->amg_complex_smoother() = 7;
    //hypre_preconditioner_pt->amg_using_simple_smoothing();
    //hypre_preconditioner_pt->amg_simple_smoother() = 1;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 0;
    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  } 

  Preconditioner* set_hypre_for_RSPilutStrn075()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_complex_smoothing();
    hypre_preconditioner_pt->amg_complex_smoother() = 7;
    //hypre_preconditioner_pt->amg_using_simple_smoothing();
    //hypre_preconditioner_pt->amg_simple_smoother() = 1;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 1;
    hypre_preconditioner_pt->amg_strength() = 0.75;

    return another_preconditioner_pt;
  }
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  Preconditioner* set_hypre_for_CLJPGSStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 1;
    //hypre_preconditioner_pt->amg_damping() = 1.0;

    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 0;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_CLJPJStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 0;
    hypre_preconditioner_pt->amg_damping() = 1.0;

    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 0;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_CLJPPilutStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    //hypre_preconditioner_pt->amg_using_simple_smoothing();
    //hypre_preconditioner_pt->amg_simple_smoother() = 0;
    //hypre_preconditioner_pt->amg_damping() = 1.0;

    hypre_preconditioner_pt->amg_using_complex_smoothing();
    hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 0;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  //
  Preconditioner* set_hypre_for_RSGSStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 1;
    //hypre_preconditioner_pt->amg_damping() = 1.0;

    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 1;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_RSJStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    hypre_preconditioner_pt->amg_using_simple_smoothing();
    hypre_preconditioner_pt->amg_simple_smoother() = 0;
    hypre_preconditioner_pt->amg_damping() = 1.0;

    //hypre_preconditioner_pt->amg_using_complex_smoothing();
    //hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 1;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
  Preconditioner* set_hypre_for_RSPilutStrn0668()
  {
    Preconditioner* another_preconditioner_pt =  
      new HyprePreconditioner;
    HyprePreconditioner* hypre_preconditioner_pt = 
      static_cast<HyprePreconditioner*>(another_preconditioner_pt);

    Hypre_default_settings::
      set_defaults_for_2D_poisson_problem(hypre_preconditioner_pt);

    ///////////////////////////////////////////////////////////////////////////
    // Parameters from Milan.

    /// \short Simple smoothing methods used in BoomerAMG. Relaxation types
    /// include:
    ///  0 = Jacobi 
    ///  1 = Gauss-Seidel, sequential
    ///      (very slow in parallel!)
    /// To use these methods set AMG_using_simple_smoothing to true
    //hypre_preconditioner_pt->amg_using_simple_smoothing();
    //hypre_preconditioner_pt->amg_simple_smoother() = 0;
    //hypre_preconditioner_pt->amg_damping() = 1.0;

    hypre_preconditioner_pt->amg_using_complex_smoothing();
    hypre_preconditioner_pt->amg_complex_smoother() = 7;

    // AMG coarsening strategy. Coarsening types include:
    ///  0 = CLJP (parallel coarsening using independent sets)
    ///  1 = classical RS with no boundary treatment (not recommended
    ///      in parallel)
    hypre_preconditioner_pt->amg_coarsening() = 1;
    hypre_preconditioner_pt->amg_strength() = 0.668;

    return another_preconditioner_pt;
  }
} // end of namespace Hypre_Subsidiary_Preconditioner_Helper
#endif // End of ifdef OOMPH_HAS_HYPRE

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

namespace StepLagrange
{
  // Prob id, set by main method
  int* Prob_id_pt = 0;

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

  inline void generic_setup()
  {
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
      int prob_id_array[]= {10,11,12,13};

      bool inset = check_if_in_set<int>(prob_id_array,4,(*Prob_id_pt));

      // Check if they have provided an acceptable ID.
      // If a new element has been inserted, it means the user has provided an
      // ID not in the set.
      if(inset == false)
      {
        std::ostringstream err_msg;
        err_msg << "Please provide a problem id to identify the problem after "
          << "after the argument --prob_id.\n" 
          << "Acceptable IDs are:\n"
          << "10 = (SqTmp) Square, custom stuff...\n"
          << "11 = (SqPo) Square, Parallel outflow (para inflow)\n"
          << "12 = (SqTf) Square, Tangential flow (Semi para inflow)\n"
          << "13 = (SqTfPo) Square, Tangential flow, Parallel outflow (semi para inflow)\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
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
    // We know that --prob_id is set and is set up properly. We now set the
    // Prob_str.
    const int prob_id = *Prob_id_pt;
    switch(prob_id)
    {
      case 10:
        Prob_str = "StTmp";
        break;
      case 11:
        Prob_str = "StPo";
        break;
      case 12:
        Prob_str = "StTf";
        break;
      case 13:
        Prob_str = "StTfPo";
        break;
      default:
        {
          std::ostringstream err_msg;
          err_msg << "There is an unrecognised Prob_id, recognised Prob_id:\n"
            << "10 = (StTmp) Step, custom stuff...\n"
            << "11 = (StPo) Step, Parallel outflow (para inflow)\n"
            << "12 = (StTf) Step, Tangential flow (Semi para inflow)\n"
            << "13 = (StTfPo) Step, Tangential flow, Parallel outflow (semi para inflow)\n"
            << std::endl;
          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        } // Default case
    } // switch Prob_id

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
    // Now we need to convert Ang into radians.
    Ang = Ang_deg * (MathematicalConstants::Pi / 180.0);

    // Now we set the Ang_deg_str.
    // this exists for only the Sq problems, not Aw.
    {
      if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream strs;
        strs << "A" << Ang_deg;
        Ang_deg_str = strs.str();
      }
      else
      {
        std::ostringstream err_msg;
        err_msg << "You have selected an Sq problem."
          << "Please supply the tilting angle with --ang.\n"
          << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
    }

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

  inline string create_label()
  {
    std::string label = Prob_str + Ang_deg_str + Noel_str;

    return label; 
  } // inlined function create_label

} // Namespace SquareLagrange



namespace SquareLagrange
{


  // Prob id, set by main method
  int* Prob_id_pt = 0;

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
  double X_min = 0.0;
  double X_max = 1.0;

  // Min and max y value respectively.
  double Y_min = 0.0;
  double Y_max = 1.0;

  // The length in the x and y direction respectively.
  double Lx = X_max - X_min;
  double Ly = Y_max - Y_min;

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
      int prob_id_array[]= {10,11,12,13,
        20,21,22,23};

      bool inset = check_if_in_set<int>(prob_id_array,8,(*Prob_id_pt));

      // Check if they have provided an acceptable ID.
      // If a new element has been inserted, it means the user has provided an
      // ID not in the set.
      if(inset == false)
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
          << "20 = (AwTmp) Annulus wedge, custom stuff...\n"
          << "21 = (AwPo) Annulus wedge, Parallel outflow (para inflow)\n"
          << "22 = (AwTf) Annulus wedge, Tangential flow (semi para inflow)\n"
          << "23 = (AwTfPo) Annulus wedge, Tan. flow, Para. outflow (semi para inflow)\n"
          << std::endl;

        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
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
    // We know that --prob_id is set and is set up properly. We now set the
    // Prob_str.
    const int prob_id = *Prob_id_pt;
    switch(prob_id)
    {
      case 10:
        Prob_str = "SqTmp";
        break;
      case 11:
        Prob_str = "SqPo";
        break;
      case 12:
        Prob_str = "SqTf";
        break;
      case 13:
        Prob_str = "SqTfPo";
        break;
      case 20:
        Prob_str = "AwTmp";
        break;
      case 21:
        Prob_str = "AwPo";
        break;
      case 22:
        Prob_str = "AwTf";
        break;
      case 23:
        Prob_str = "AwTfPo";
        break;
      default:
        {
          std::ostringstream err_msg;
          err_msg << "There is an unrecognised Prob_id, recognised Prob_id:\n"
            << "10 = (SqTmp) Square, custom stuff...\n"
            << "11 = (SqPo) Square, Parallel outflow (para inflow)\n"
            << "12 = (SqTf) Square, Tangential flow (Semi para inflow)\n"
            << "13 = (SqTfPo) Square, Tangential flow, Parallel outflow (semi para inflow)\n"
            << "\n"
            << "20 = (AwTmp) Annulus wedge, custom stuff...\n"
            << "21 = (AwPo) Annulus wedge, Parallel outflow (para inflow)\n"
            << "22 = (AwTf) Annulus wedge, Tangential flow (semi para inflow)\n"
            << "23 = (AwTfPo) Annulus wedge, Tan. flow, Para. outflow (semi para inflow)\n"
            << std::endl;
          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        } // Default case
    } // switch Prob_id
  }

  inline void set_ang_str()
  {
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

    // Now we need to convert Ang into radians.
    Ang = Ang_deg * (MathematicalConstants::Pi / 180.0);
    // Now we set the Ang_deg_str.
    {
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
    set_prob_str();
    set_ang_str();
    set_noel_str();
  }

  inline string prob_str()
  {
    set_prob_str();
    return Prob_str;
  }

  inline string ang_deg_str()
  {
    set_ang_str();
    return Ang_deg_str;
  }

  inline string noel_str()
  {
    set_noel_str();
    return Noel_str;
  }

  inline string create_label()
  {
    std::string label = prob_str() + ang_deg_str() + noel_str();
    return label; 
  } // inlined function create_label

} // Namespace SquareLagrange


namespace LagrangianPreconditionerHelpers
{
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
  int f_amg_smoother = -1;
  int f_amg_iterations = -1;
  int f_amg_smoother_iterations = -1;

  double p_amg_strength = -1.0;
  double p_amg_damping = -1.0;
  int p_amg_coarsening = -1;
  int p_amg_smoother = -1;
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


  inline void setup_commandline_flags()
  {
    // Flag to output the preconditioner, used for debugging.
    // string
    CommandLineArgs::specify_command_line_flag("--doc_prec",&Doc_prec_dir_str);

    // double
    CommandLineArgs::specify_command_line_flag("--sigma",&Scaling_sigma);

    // int
    CommandLineArgs::specify_command_line_flag("--w_solver",&W_solver);

    // Nothing set
    CommandLineArgs::specify_command_line_flag("--bdw");

    // int
    CommandLineArgs::specify_command_line_flag("--ns_solver",&NS_solver);

    // int
    CommandLineArgs::specify_command_line_flag("--p_solver",&P_solver);

    // int
    CommandLineArgs::specify_command_line_flag("--f_solver",&F_solver);

    // NS_F block AMG parameters
    // double
    CommandLineArgs::specify_command_line_flag("--f_amg_str",&f_amg_strength);
    // double
    CommandLineArgs::specify_command_line_flag("--f_amg_damp",&f_amg_damping);

    // int
    CommandLineArgs::specify_command_line_flag("--f_amg_coarse",
        &f_amg_coarsening);

    // int
    CommandLineArgs::specify_command_line_flag("--f_amg_smoo",&f_amg_smoother);

    // int
    CommandLineArgs::specify_command_line_flag("--f_amg_iter",&f_amg_iterations);

    //int
    CommandLineArgs::specify_command_line_flag("--f_amg_smiter",
        &f_amg_smoother_iterations);

    // NS_P block AMG parameters
    // double
    CommandLineArgs::specify_command_line_flag("--p_amg_str",&p_amg_strength);
    // double
    CommandLineArgs::specify_command_line_flag("--p_amg_damp",&p_amg_damping);
    // int
    CommandLineArgs::specify_command_line_flag("--p_amg_coarse",&p_amg_coarsening);
    // int
    CommandLineArgs::specify_command_line_flag("--p_amg_smoo",&p_amg_smoother);
    // int
    CommandLineArgs::specify_command_line_flag("--p_amg_iter",&p_amg_iterations);
    // int
    CommandLineArgs::specify_command_line_flag("--p_amg_smiter",
        &p_amg_smoother_iterations);

    CommandLineArgs::specify_command_line_flag("--print_hypre");
  }

  inline void generic_setup()
  {
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
        err_msg << "Please provide the doc_prec directory "
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

  inline Preconditioner* get_preconditioner()
  {
    LagrangeEnforcedflowPreconditioner* prec_pt
      = new LagrangeEnforcedflowPreconditioner;

    // Set the mesh
    if(Mesh_pt.size() == 0)
    {
      std::ostringstream err_msg;
      err_msg << "There is no Mesh_pt set.\n"
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
      // If ns_solver is 1, this means we want to use LSC.
      // So the F_solver and P_solver must be set.
      if((F_solver == -1) || (P_solver == -1))
      {
        std::ostringstream err_msg;
        err_msg << "Doing LSC NS solve. (NS_solver is 1)\n"
          << "but F_solver and P_solver have not been set.\n"
          << "0 - Exact (SuperLU)"
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

      // Set the NS preconditioner as LSC.
      prec_pt->set_navier_stokes_lsc_preconditioner(ns_preconditioner_pt);

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
        // AMG coarsening: Ruge-Stuben
        RayGlobalAMGParam::amg_coarsening = 1;

        // AMG smoother: Gauss-Seidel
        RayGlobalAMGParam::amg_smoother=0;

        // There is no damping with GS, otherwise we set the parameter:
        // RayGlobalAMGParam::amg_damping

        // Different amg strength for simple/stress divergence for viscuous term.
        const int vis = *Vis_pt;
        if(vis == 0)
        {
          // Simple form
          RayGlobalAMGParam::amg_strength = 0.25;
        }
        else if (vis == 1)
        {
          // Stress divergence form
          RayGlobalAMGParam::amg_strength = 0.668;
        }
        else
        {
          std::ostringstream err_msg;
          err_msg << "Do not recognise viscuous term: " << vis << std::endl;

          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }

        // Setup the preconditioner.
        f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
          set_hypre_using_2D_poisson_base();

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
        // AMG coarsening:
        // Set: RayGlobalAMGParam::amg_coarsening = 
        // 0 - CLJP
        // 1 - RS

        // AMG smoother:
        // Set: RayGlobalAMGParam::amg_smoother = 
        // 0 - Jacobi (Need to set damping as well)
        // 1 - Gauss-Seidel
        // 2 - Pilut

        // There is no damping with GS, otherwise we set the parameter:
        // RayGlobalAMGParam::amg_damping

        RayGlobalAMGParam::amg_strength = f_amg_strength;
        RayGlobalAMGParam::amg_damping = f_amg_damping;
        RayGlobalAMGParam::amg_coarsening = f_amg_coarsening;
        RayGlobalAMGParam::amg_smoother = f_amg_smoother;
        RayGlobalAMGParam::amg_iterations = f_amg_iterations;
        RayGlobalAMGParam::amg_smoother_iterations 
          = f_amg_smoother_iterations;
        RayGlobalAMGParam::print_hypre = Print_hypre;


        // Different amg strength for simple/stress divergence for viscuous term.
        const int vis = (*Vis_pt);
        if(RayGlobalAMGParam::amg_strength < 0.0)
        {
          if(vis == 0)
          {
            // Simple form
            RayGlobalAMGParam::amg_strength = 0.25;
          }
          else if(vis == 1)
          {
            // Stress divergence form
            RayGlobalAMGParam::amg_strength = 0.668;
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
          set_hypre_ray();
#endif
      }

      // Set the preconditioner in the LSC preconditioner.
      ns_preconditioner_pt->set_f_preconditioner(f_preconditioner_pt);

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

        // Set it as the p preconditioner for LSC
        //     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif
      }
      else if(P_solver == 96)
      {
#ifdef OOMPH_HAS_HYPRE
        //* 
        RayGlobalAMGParam::amg_iterations = p_amg_iterations;
        RayGlobalAMGParam::amg_smoother_iterations = p_amg_smoother_iterations;
        RayGlobalAMGParam::amg_smoother = p_amg_smoother;
        RayGlobalAMGParam::amg_strength = p_amg_strength;
        //RayGlobalAMGParam::amg_damping = param_pt->p_amg_damping;
        RayGlobalAMGParam::amg_coarsening = p_amg_coarsening;
        RayGlobalAMGParam::print_hypre = Print_hypre;

        //     std::cout << "p_amg_iterations:" << SL::p_amg_iterations << std::endl; 
        //     std::cout << "p_amg_smoother_iterations" << SL::p_amg_smoother_iterations << std::endl; 
        //     std::cout << "p_amg_strength" << SL::p_amg_strength << std::endl;
        //     std::cout << "p_amg_coarsening" << SL::p_amg_coarsening << std::endl; 
        // */

        p_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::set_hypre_ray();

        //     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
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

        //     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif

      }

      ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
    } // if for using LSC as NS prec.
    else
    {
      pause("There is no solver for NS.");
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
  } // LPH::setup_preconditioner


  inline string create_label()
  {
    std::string w_str = "";
    std::string ns_str = "";
    std::string f_str = "";
    std::string p_str = "";
    std::string sigma_str = "";

    // Set the string for W_solver.
    switch(W_solver)
    {
      case 0:
        w_str = "We";
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
          p_str = "";
          f_str = "";
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

    // Now we continue with setting the string for the solvers.
    // Only set the f_str if NS_solver > 0
    if(NS_solver > 0)
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

    if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
    {
      std::ostringstream strs;
      strs << "S" << Scaling_sigma;
      sigma_str = strs.str();
    }

    std::string prec_str = w_str + ns_str + f_str + p_str + sigma_str;

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

    // Iteration count and times directory.
    CommandLineArgs::specify_command_line_flag("--itstimedir", 
        &Itstime_dir_str);

    CommandLineArgs::specify_command_line_flag("--trilinos_solver");
  }

  inline void generic_problem_setup(const unsigned dim)
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
      if(dim == 2)
      {
        if (Vis == 0)
        {
          NavierStokesEquations<2>::Gamma[0]=0.0;
          NavierStokesEquations<2>::Gamma[1]=0.0;
        }
        else if (Vis == 1)
        {
          NavierStokesEquations<2>::Gamma[0]=1.0;
          NavierStokesEquations<2>::Gamma[1]=1.0;
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
      }
      else
      {
        if (Vis == 0)
        {
          NavierStokesEquations<3>::Gamma[0]=0.0;
          NavierStokesEquations<3>::Gamma[1]=0.0;
          NavierStokesEquations<3>::Gamma[2]=0.0;
        }
        else if (Vis == 1)
        {
          NavierStokesEquations<3>::Gamma[0]=1.0;
          NavierStokesEquations<3>::Gamma[1]=1.0;
          NavierStokesEquations<3>::Gamma[2]=1.0;
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
  } // NSPP::generic_problem_setup()

  inline string create_label()
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
  inline void setup_solver(const double& solver_tol, const double& newton_tol,
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

    solver_pt->tolerance() = solver_tol;
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

