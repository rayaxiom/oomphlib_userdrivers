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

namespace NavierStokesProblemParameters
{
  int Prob_id = -1;
  bool Distribute_problem = false;
  int Vis = -1;
  double Rey = -1.0;

  double Rey_start = -1.0;
  double Rey_incre = -1.0;
  double Rey_end = -1.0;

  bool Doc_soln = false;

  std::string Soln_dir_str = "";
  std::string Itstime_dir_str = "";

  bool Using_trilinos_solver = false;

  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;


//  void setup_commandline_flags();

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
  }

  void generic_problem_setup(const unsigned dim);

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

    // Document the solution? Default is false.
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

//// RAYRAY here. Do generic Reynolds number ssetup!!!
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

  } // generic_problem_setup()

  string create_label();

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
}


namespace LagrangianPreconditionerHelpers
{
  struct PrecParam
  {
    PrecParam()
    {
      W_solver = -1;
      NS_solver = -1;
      F_solver = -1;
      P_solver = -1;

      f_amg_strength = -1.0;
      f_amg_damping = -1.0;
      f_amg_coarsening = -1;
      f_amg_smoother = -1;
      f_amg_iterations = -1;
      f_amg_smoother_iterations = -1;

      p_amg_strength = -1.0;
      p_amg_damping = -1.0;
      p_amg_coarsening = -1;
      p_amg_smoother = -1;
      p_amg_iterations = -1;
      p_amg_smoother_iterations = -1;

      Mesh_pt.resize(0);

      Scaling_sigma = 0;

      Use_axnorm = true;

      Use_block_diagonal_w = false;

      Doc_prec = false;
      Doc_prec_dir_str = ""; // done

      Print_hypre = true;

      Problem_pt = 0;

      Vis = -1;
    }

    Vector<Mesh*> Mesh_pt;

    double Scaling_sigma;

    bool Use_axnorm;
    bool Use_block_diagonal_w;
    bool Doc_prec;
    bool Print_hypre;

    std::string Doc_prec_dir_str;

    DocLinearSolverInfo* Doc_linear_solver_info_pt;

    Problem* Problem_pt;

    int Vis;

    int W_solver;
    int NS_solver;
    int F_solver;
    int P_solver;

    double f_amg_strength;
    double f_amg_damping;
    int f_amg_coarsening;
    int f_amg_smoother;
    int f_amg_iterations;
    int f_amg_smoother_iterations;

    double p_amg_strength;
    double p_amg_damping;
    int p_amg_coarsening;
    int p_amg_smoother;
    int p_amg_iterations;
    int p_amg_smoother_iterations;
  };


  inline void setup_commandline_flags(PrecParam* param_pt)
  {
  // Flag to output the preconditioner, used for debugging.
  CommandLineArgs::specify_command_line_flag("--doc_prec", 
                                             &param_pt->Doc_prec_dir_str);

  CommandLineArgs::specify_command_line_flag("--sigma",
                                             &param_pt->Scaling_sigma);
  
  CommandLineArgs::specify_command_line_flag("--w_solver", 
                                             &param_pt->W_solver);
  CommandLineArgs::specify_command_line_flag("--bdw");

  CommandLineArgs::specify_command_line_flag("--ns_solver", 
                                             &param_pt->NS_solver);
  CommandLineArgs::specify_command_line_flag("--p_solver", 
                                             &param_pt->P_solver);
  CommandLineArgs::specify_command_line_flag("--f_solver", 
                                             &param_pt->F_solver);

  // NS_F block AMG parameters
  CommandLineArgs::specify_command_line_flag("--f_amg_str", 
                                             &param_pt->f_amg_strength);
  CommandLineArgs::specify_command_line_flag("--f_amg_damp", 
                                             &param_pt->f_amg_damping);
  CommandLineArgs::specify_command_line_flag("--f_amg_coarse", 
                                             &param_pt->f_amg_coarsening);
  CommandLineArgs::specify_command_line_flag("--f_amg_smoo", 
                                             &param_pt->f_amg_smoother);
  CommandLineArgs::specify_command_line_flag("--f_amg_iter", 
                                             &param_pt->f_amg_iterations);
  CommandLineArgs::specify_command_line_flag(
      "--f_amg_smiter", &param_pt->f_amg_smoother_iterations);

  // NS_P block AMG parameters
  CommandLineArgs::specify_command_line_flag("--p_amg_str", 
                                             &param_pt->p_amg_strength);
  CommandLineArgs::specify_command_line_flag("--p_amg_damp", 
                                             &param_pt->p_amg_damping);
  CommandLineArgs::specify_command_line_flag("--p_amg_coarse", 
                                             &param_pt->p_amg_coarsening);
  CommandLineArgs::specify_command_line_flag("--p_amg_smoo", 
                                             &param_pt->p_amg_smoother);
  CommandLineArgs::specify_command_line_flag("--p_amg_iter", 
                                             &param_pt->p_amg_iterations);
  CommandLineArgs::specify_command_line_flag(
      "--p_amg_smiter", &param_pt->p_amg_smoother_iterations);
  }

  inline void generic_setup(PrecParam* param_pt)
  {

    if(param_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "No param_pt set.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    // Document the preconditioner? Default is false.
    if(CommandLineArgs::command_line_flag_has_been_set("--doc_prec"))
    {
      // The argument immediately after --doc_prec is put into SL::Doc_prec_dir.
      // If this begins with "--", then no prec directory has been provided.
      const std::string doc_prec_dir_str = param_pt->Doc_prec_dir_str;

      std::size_t found = doc_prec_dir_str.find("--");

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
        param_pt->Doc_prec = true;
      }
    }

  } // generic_setup()

  Preconditioner* setup_preconditioner(const PrecParam*);

  inline Preconditioner* setup_preconditioner(const PrecParam* param_pt = 0)
  {
    LagrangeEnforcedflowPreconditioner* prec_pt
      = new LagrangeEnforcedflowPreconditioner;


    // Set the mesh
    Vector<Mesh*> mesh_pt = param_pt->Mesh_pt;
    if(mesh_pt.size() == 0)
    {
      std::ostringstream err_msg;
      err_msg << "There is no Mesh_pt set.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    prec_pt->set_meshes(mesh_pt);

    // Set W solver.
    const int w_solver = param_pt->W_solver;
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
 const int ns_solver = param_pt->NS_solver;
   const int f_solver = param_pt->F_solver;
   const int p_solver = param_pt->P_solver;
 if(ns_solver == -1)
 {
      std::ostringstream err_msg;
      err_msg << "The NS solver has not been set.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
 }
 else if(ns_solver == 0) // Exact solve.
 {
   // This is the default, do nothing.
   // But the param's F_solver and P_solver should not have been set,
   // i.e. it should stay as -1.

   if((f_solver != -1) || (p_solver != -1))
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
 else if(ns_solver == 1) // LSC
 {
   // If ns_solver is 1, this means we want to use LSC.
   // So the F_solver and P_solver must be set.
   if((f_solver == -1) || (p_solver == -1))
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

   // Check if the problem pointer is set.
   Problem* problem_pt = param_pt->Problem_pt;

   // Check that the problem pointer is set (not null).
   // LSC requires a problem pointer.
   if(problem_pt != 0)
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
     new NavierStokesSchurComplementPreconditioner(problem_pt);

   // Set the NS preconditioner as LSC.
   prec_pt->set_navier_stokes_lsc_preconditioner(ns_preconditioner_pt);

   // Give LSC the bulk mesh (Navier-Stokes mesh).
   ns_preconditioner_pt->set_navier_stokes_mesh(mesh_pt[0]);


   //// Setting the F solver within the NS block
   /////////////////////////////////////////////

   // Preconditioner for the F block:
   Preconditioner* f_preconditioner_pt = 0;

   const int f_solver = param_pt->F_solver;

   // f_solver == 0 is default, so do nothing.
   //
   // AMG depends on the Reynolds number so we check that the Reynolds number
   // is set if LSC is switched on - even if we do not want to use AMG...
   // for consistency.
   const int vis = param_pt->Vis;
   if(vis == -1)
   {
      std::ostringstream err_msg;
      err_msg << "Please set your viscuous term.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
   }

   if(f_solver == 11)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_2D_poison_problem();
#endif
   }
   else if(f_solver == 12)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_navier_stokes_momentum_block();
#endif
   }
   else if(f_solver == 13)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPGSStrn075();
#endif
   }
   else if(f_solver == 14)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSGSStrn075();
#endif
   }
   else if(f_solver == 15)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPPilutStrn075();
#endif
   }
   else if(f_solver == 16)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSPilutStrn075();
#endif
   }
   else if(f_solver == 17)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_augmented_momentum_block();
#endif
   }
   else if(f_solver == 81)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPGSStrn0668();
#endif
   }
   else if(f_solver == 82)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPJStrn0668();
#endif
   }
   else if(f_solver == 83)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_CLJPPilutStrn0668();
#endif
   }
   else if(f_solver == 84)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSGSStrn0668();
#endif
   }
   else if(f_solver == 85)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSJStrn0668();
#endif
   }
   else if(f_solver == 86)
   {
#ifdef OOMPH_HAS_HYPRE
     // LSC takes type "Preconditioner".
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_for_RSPilutStrn0668();
#endif
   }
   else if(f_solver == 2)
   {
//     f_preconditioner_pt = new RayBlockDiagonalPreconditioner<CRDoubleMatrix>;
     f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
   }
   else if(f_solver == 3)
   {
     f_preconditioner_pt = new BlockDiagonalPreconditioner<CRDoubleMatrix>;
#ifdef OOMPH_HAS_HYPRE
     dynamic_cast<BlockDiagonalPreconditioner<CRDoubleMatrix>* >
       (f_preconditioner_pt)->set_subsidiary_preconditioner_function
       (Hypre_Subsidiary_Preconditioner_Helper::set_hypre_for_2D_poison_problem);
#endif
   }
   else if (f_solver == 69)
   {
#ifdef OOMPH_HAS_HYPRE
     // AMG coarsening: Ruge-Stuben
     RayGlobalAMGParam::amg_coarsening = 1;
     
     // AMG smoother: Gauss-Seidel
     RayGlobalAMGParam::amg_smoother=0;
     
     // There is no damping with GS, otherwise we set the parameter:
     // RayGlobalAMGParam::amg_damping

     // Different amg strength for simple/stress divergence for viscuous term.
     if(vis == 0)
     {
       // Simple form
       RayGlobalAMGParam::amg_strength = 0.25;
     }
     else
     {
       // Stress divergence form
       RayGlobalAMGParam::amg_strength = 0.668;
     }
     
     // Setup the preconditioner.
     f_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::
       set_hypre_using_2D_poisson_base();
#endif
   }
   else if (f_solver == 96)
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

     RayGlobalAMGParam::amg_strength = param_pt->f_amg_strength;
     RayGlobalAMGParam::amg_damping = param_pt->f_amg_damping;
     RayGlobalAMGParam::amg_coarsening = param_pt->f_amg_coarsening;
     RayGlobalAMGParam::amg_smoother = param_pt->f_amg_smoother;
     RayGlobalAMGParam::amg_iterations = param_pt->f_amg_iterations;
     RayGlobalAMGParam::amg_smoother_iterations 
       = param_pt->f_amg_smoother_iterations;
     RayGlobalAMGParam::print_hypre = param_pt->Print_hypre;


     // Different amg strength for simple/stress divergence for viscuous term.
     if(RayGlobalAMGParam::amg_strength < 0.0)
     {
       if(vis == 0)
       {
         // Simple form
         RayGlobalAMGParam::amg_strength = 0.25;
       }
       else
       {
         // Stress divergence form
         RayGlobalAMGParam::amg_strength = 0.668;
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

   const int p_solver = param_pt->P_solver;

   //SL::P_solver == 0 is default, so do nothing.
   if(p_solver == 1) 
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
   else if(p_solver == 96)
   {
#ifdef OOMPH_HAS_HYPRE
//* 
     RayGlobalAMGParam::amg_iterations = param_pt->p_amg_iterations;
     RayGlobalAMGParam::amg_smoother_iterations = param_pt->p_amg_smoother_iterations;
     RayGlobalAMGParam::amg_smoother = param_pt->p_amg_smoother;
     RayGlobalAMGParam::amg_strength = param_pt->p_amg_strength;
     //RayGlobalAMGParam::amg_damping = param_pt->p_amg_damping;
     RayGlobalAMGParam::amg_coarsening = param_pt->p_amg_coarsening;
     RayGlobalAMGParam::print_hypre = param_pt->Print_hypre;

//     std::cout << "p_amg_iterations:" << SL::p_amg_iterations << std::endl; 
//     std::cout << "p_amg_smoother_iterations" << SL::p_amg_smoother_iterations << std::endl; 
//     std::cout << "p_amg_strength" << SL::p_amg_strength << std::endl;
//     std::cout << "p_amg_coarsening" << SL::p_amg_coarsening << std::endl; 
// */

     p_preconditioner_pt = Hypre_Subsidiary_Preconditioner_Helper::set_hypre_ray();

//     ns_preconditioner_pt->set_p_preconditioner(p_preconditioner_pt);
#endif
   }
   else if(p_solver == 2)
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
//// else
//// {
////   pause("There is no solver for NS.");
//// }      std::ostringstream err_msg;

    return prec_pt;
  } // setup_preconditioner


  string create_label(const PrecParam* param_pt);

  inline string create_label(const PrecParam* param_pt)
  {
    if(param_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "No param_pt set.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    std::string w_str = "";
    std::string ns_str = "";
    std::string f_str = "";
    std::string p_str = "";
    std::string sigma_str = "";

    const unsigned w_solver = param_pt->W_solver;
    // Set the string for W_solver.
    switch(w_solver)
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

    const bool use_block_diagonal_w = param_pt->Use_block_diagonal_w;
    if(use_block_diagonal_w)
    {
      w_str += "bd";
    }
    else
    {
      w_str += "d";
    }


   


    std::string prec_str = "";

    return prec_str;

  }

} // end of namespace LagrangianPreconditionerHelpers

namespace SquareLagrange
{

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

  // Default configuration
  int W_solver = 0; //CL, 0 = SuperLU, no other W solver coded.
  int NS_solver = 1; //CL, 0 = SuperLU, 1 - LSC
  int F_solver = 0; //CL, 0 - SuperLU, 1 - AMG
  int P_solver = 0; //CL, 0 - SuperLU, 1 - AMG

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
  int Prob_id = -1;

  // These are self explanatory:
  unsigned Vis = 0; //CL, 0 - Simple, 1 - Stress divergence
  double Ang_deg = 30.0; //CL, Angle in degrees
  double Ang = 0.0; //CL, Angle in degrees
  double Rey = 100.0; //CL, Reynolds number
  unsigned Noel = 4; //CL, Number of elements in 1D
  double Scaling_sigma = 0; //CL, If the scaling sigma is not set, then
  // the default is the norm of the momentum block.

  // These are strings for 
  bool Use_axnorm = true; //Set from CL, use norm for sigma?
  bool Use_block_diagonal_w = false; // To set from CL
  bool Doc_prec = false; // To set from CL
  bool Doc_soln = false; // To set from CL
  bool Print_hypre = true;

  std::string Label = ""; // To be set as the label for this problem. Contains
  // all the information for this run.
  std::string Soln_dir = ""; // Where to put the solution.
  std::string Doc_prec_dir = ""; // Where to put the solution.

  std::string Itstime_dir = ""; //Set from CL, directory to output the 
  // iteration counts and timing results.

  // Used to determine if we are using the TrilinosAztecOOSolver solver or not.
  // This cannot be determined by the OOMPH_HAS_TRILINOS ifdef since we may be
  // using OOMPH-LIB's GMRES even if we have Trilinos. This should be set in
  // the problem constuctor as soon as we set the linear_solver_pt() for the
  // problem.
  bool Using_trilinos_solver = false;

  double Rey_start = 0.0;
  double Rey_incre = 50.0;
  double Rey_end = 500.0;

  // Object to store the linear solver iterations and times.
  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;

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

  inline void setup_commandline_flags()
  {
  CommandLineArgs::specify_command_line_flag("--ang", &Ang_deg);

  CommandLineArgs::specify_command_line_flag("--noel", &Noel);
  }

  inline void generic_setup(LagrangianPreconditionerHelpers::PrecParam* param_pt)
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

    if(CommandLineArgs::command_line_flag_has_been_set("--visc"))
    {
      param_pt->Vis = NavierStokesProblemParameters::Vis;
    }
    else
    {
        std::ostringstream err_msg;
        err_msg << "Viscuous term has not been set. Please set: \n"
          << "--visc \n"; 
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
    }


  // Set a problem id to identify the problem.
  // This is used for book keeping purposes.
  if(CommandLineArgs::command_line_flag_has_been_set("--prob_id"))
  {
    // The argument immediately after --prob_id is put into SL::Prob_id.
    // If this begins with "--", then no problem id has been provided.

    // Maybe I should check if SL::Prob_id is a number or a string...

    // We only accept problem IDs as defined below.
    // Creating a set of acceptable IDs
    int prob_id_array[]= {10,11,12,13,
      20,21,22,23};

    namespace NSPP = NavierStokesProblemParameters;
    bool inset = check_if_in_set<int>(prob_id_array,8,NSPP::Prob_id);

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

  }

  string create_label(
      const LagrangianPreconditionerHelpers::PrecParam* prec_param_pt);

  inline string create_label(
      const LagrangianPreconditionerHelpers::PrecParam* prec_param_pt)
  {

    if(prec_param_pt == 0)
    {
      std::ostringstream err_msg;
      err_msg << "No param_pt set.\n"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    std::string prob_str = "";
    std::string prec_str = "";
    std::string w_str = "";
    std::string ns_str = "";
    std::string f_str = "";
    std::string p_str = "";
    std::string vis_str = "";
    std::string ang_str = "";
    std::string rey_str = "";
    std::string noel_str = "";
    std::string w_approx_str = "";
    std::string sigma_str = "";

    switch(Prob_id)
    {
      case 10:
        prob_str = "SqTmp";
        break;
      case 11:
        prob_str = "SqPo";
        break;
      case 12:
        prob_str = "SqTf";
        break;
      case 13:
        prob_str = "SqTfPo";
        break;
      case 20:
        prob_str = "AwTmp";
        break;
      case 21:
        prob_str = "AwPo";
        break;
      case 22:
        prob_str = "AwTf";
        break;
      case 23:
        prob_str = "AwTfPo";
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


    prec_str = LagrangianPreconditionerHelpers::create_label(prec_param_pt);
    pause("pp from SL create label"); 
    

    // Set the string for W_solver.
    switch(W_solver)
    {
      case 0:
        w_str = "We";
        break;
      default:
        {
          std::ostringstream err_msg;
          err_msg << "There is an unrecognised W_solver, recognised W_solver:\n"
            << "0 = (We) SuperLU solve\n"
            << std::endl;
          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
    } // switch

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
          err_msg << "There is an unrecognised NS_solver, recognised NS_solver:\n"
            << "0 = (Ne) SuperLU\n"
            << "1 = (Nl) LSC preconditioner\n"
            << std::endl;
          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
    } // switch NS_solver

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
      }  // switch
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
      } // switch
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

    // Set Ang_str, this exists for only the Sq problems, not Aw.
    std::size_t found = prob_str.find("Sq");
    if(found != std::string::npos)
    {
      if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream strs;
        strs << "A" << Ang_deg;
        ang_str = strs.str();
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
    else
    {
      if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream err_msg;
        err_msg << "You have selected a Aw problem, there is no tilting angle.\n"
          << "Please take off the --ang command line argument.\n"
          << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Set the Reynolds number.
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

    // Set Noel_str, used for book keeping.
    if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
    {
      std::ostringstream strs;
      strs << "N" <<Noel;
      noel_str = strs.str();
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

    // Use the diagonal or block diagonal approximation for W block.
    if(CommandLineArgs::command_line_flag_has_been_set("--bdw"))
    {
      Use_block_diagonal_w = true;
      w_approx_str = "Wbd";
    }
    else
    {
      Use_block_diagonal_w = false;
      w_approx_str = "Wd";
    }

    // Set Use_axnorm, if sigma has not been set, 
    // norm os momentum block is used.
    if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
    {
      Use_axnorm = false;

      std::ostringstream strs;
      strs << "S" << Scaling_sigma;
      sigma_str = strs.str();
    }

    std::string label = prob_str
      + w_str + ns_str + f_str + p_str
      + vis_str + ang_str + rey_str
      + noel_str + w_approx_str + sigma_str;

    return label; 
  } // inlined function create_label

} // Namespace SquareLagrange

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

namespace StepLagrange
{

  bool Distribute_problem = false;

  ///////////////////////
  // Domain dimensions.//
  ///////////////////////
  //
  // This is a square domain: x,y \in [0,1]
  //

  // Min and max x value respectively.
//  double X_min = 0.0;
//  double X_max = 1.0;

  // Min and max y value respectively.
//  double Y_min = 0.0;
//  double Y_max = 1.0;

  // The length in the x and y direction respectively.
//  double Lx = X_max - X_min;
//  double Ly = Y_max - Y_min;

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

  // Default configuration
  unsigned W_solver = 0; //CL, 0 = SuperLU, no other W solver coded.
  unsigned NS_solver = 1; //CL, 0 = SuperLU, 1 - LSC
  unsigned F_solver = 0; //CL, 0 - SuperLU, 1 - AMG
  unsigned P_solver = 0; //CL, 0 - SuperLU, 1 - AMG

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
  int Prob_id = -1;

  // These are self explanatory:
  unsigned Vis = 0; //CL, 0 - Simple, 1 - Stress divergence
  double Ang_deg = 30.0; //CL, Angle in degrees
  double Ang = 0.0; //CL, Angle in degrees
  double Rey = 100.0; //CL, Reynolds number
  unsigned Noel = 4; //CL, Number of elements in 1D
  double Scaling_sigma = 0; //CL, If the scaling sigma is not set, then
  // the default is the norm of the momentum block.

  // These are strings for 
  bool Use_axnorm = true; //Set from CL, use norm for sigma?
  bool Use_block_diagonal_w = false; // To set from CL
  bool Doc_prec = false; // To set from CL
  bool Doc_soln = false; // To set from CL
  bool Print_hypre = true;

  std::string Label = ""; // To be set as the label for this problem. Contains
  // all the information for this run.
  std::string Soln_dir = ""; // Where to put the solution.
  std::string Doc_prec_dir = ""; // Where to put the solution.

  std::string Itstime_dir = ""; //Set from CL, directory to output the 
  // iteration counts and timing results.

  // Used to determine if we are using the TrilinosAztecOOSolver solver or not.
  // This cannot be determined by the OOMPH_HAS_TRILINOS ifdef since we may be
  // using OOMPH-LIB's GMRES even if we have Trilinos. This should be set in
  // the problem constuctor as soon as we set the linear_solver_pt() for the
  // problem.
  bool Using_trilinos_solver = false;

  double Rey_start = 0.0;
  double Rey_incre = 50.0;
  double Rey_end = 500.0;

  // Object to store the linear solver iterations and times.
  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;

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

  string create_label();

  inline string create_label()
  {
    std::string prob_str = "";
    std::string w_str = "";
    std::string ns_str = "";
    std::string f_str = "";
    std::string p_str = "";
    std::string vis_str = "";
    std::string ang_str = "";
    std::string rey_str = "";
    std::string noel_str = "";
    std::string w_approx_str = "";
    std::string sigma_str = "";

    switch(Prob_id)
    {
      case 10:
        prob_str = "SqTmp";
        break;
      case 11:
        prob_str = "SqPo";
        break;
      case 12:
        prob_str = "SqTf";
        break;
      case 13:
        prob_str = "SqTfPo";
        break;
      case 20:
        prob_str = "AwTmp";
        break;
      case 21:
        prob_str = "AwPo";
        break;
      case 22:
        prob_str = "AwTf";
        break;
      case 23:
        prob_str = "AwTfPo";
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


    // Set the string for W_solver.
    switch(W_solver)
    {
      case 0:
        w_str = "We";
        break;
      default:
        {
          std::ostringstream err_msg;
          err_msg << "There is an unrecognised W_solver, recognised W_solver:\n"
            << "0 = (We) SuperLU solve\n"
            << std::endl;
          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
    } // switch

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
          err_msg << "There is an unrecognised NS_solver, recognised NS_solver:\n"
            << "0 = (Ne) SuperLU\n"
            << "1 = (Nl) LSC preconditioner\n"
            << std::endl;
          throw OomphLibError(err_msg.str(),
              OOMPH_CURRENT_FUNCTION,
              OOMPH_EXCEPTION_LOCATION);
        }
    } // switch NS_solver

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
      }  // switch
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
      } // switch
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

    // Set Ang_str, this exists for only the Sq problems, not Aw.
    std::size_t found = prob_str.find("Sq");
    if(found != std::string::npos)
    {
      if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream strs;
        strs << "A" << Ang_deg;
        ang_str = strs.str();
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
    else
    {
      if(CommandLineArgs::command_line_flag_has_been_set("--ang"))
      {
        std::ostringstream err_msg;
        err_msg << "You have selected a Aw problem, there is no tilting angle.\n"
          << "Please take off the --ang command line argument.\n"
          << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Set the Reynolds number.
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

    // Set Noel_str, used for book keeping.
    if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
    {
      std::ostringstream strs;
      strs << "N" <<Noel;
      noel_str = strs.str();
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

    // Use the diagonal or block diagonal approximation for W block.
    if(CommandLineArgs::command_line_flag_has_been_set("--bdw"))
    {
      Use_block_diagonal_w = true;
      w_approx_str = "Wbd";
    }
    else
    {
      Use_block_diagonal_w = false;
      w_approx_str = "Wd";
    }

    // Set Use_axnorm, if sigma has not been set, 
    // norm os momentum block is used.
    if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
    {
      Use_axnorm = false;

      std::ostringstream strs;
      strs << "S" << Scaling_sigma;
      sigma_str = strs.str();
    }

    std::string label = prob_str
      + w_str + ns_str + f_str + p_str
      + vis_str + ang_str + rey_str
      + noel_str + w_approx_str + sigma_str;

    return label; 
  } // inlined function create_label

} // Namespace SquareLagrange

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

namespace CubeLagrange
{
  bool Distribute_problem = false;

  ///////////////////////
  // Domain dimensions.//
  /////////////////////// 
  //
  // This is a unit cube: x,y,z \in [0,1]  
  //

  // Min and max x, y and z values respectively.
  double X_min = 0.0;
  double X_max = 1.0;
  double Y_min = 0.0;
  double Y_max = 1.0;
  double Z_min = 0.0;
  double Z_max = 1.0;
  
  // Length in the x, y and z direction respectively.
  double Lx = X_max - X_min;
  double Ly = Y_max - Y_min;
  double Lz = Z_max - Z_min;
  
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

  // Default configuration
  unsigned W_solver = 0; //CL, 0 = SuperLU, no other W solver coded.
  unsigned NS_solver = 1; //CL, 0 = SuperLU, 1 - LSC
  unsigned F_solver = 0; //CL, 0 - SuperLU, 1 - AMG
  unsigned P_solver = 0; //CL, 0 - SuperLU, 1 - AMG

  // All problems based on the cube domain will be in the same file.
  // Each unique problem will have an id.
  // 00 = (CuTmp) Square, custom stuff...
  // 01 = (CuPo) Square, Parallel outflow (para inflow) 
  // 02 = (CuTf) Square, Tangential flow (Semi para inflow)
  // 03 = (CuTfPo) Square, Tangential flow, Parallel outflow (semi para inflow)

  int Prob_id = -1;

  unsigned Lagrange_multiplier_id = 42;

  // These are self explanatory:
  unsigned Vis = 0; //CL, 0 - Simple, 1 - Stress divergence
  double Ang_deg = 30.0; //CL, Angle in degrees
  double Angx_deg = 30.0; //CL, Angle in degrees
  double Angy_deg = 30.0; //CL, Angle in degrees
  double Angz_deg = 30.0; //CL, Angle in degrees
  double Angx = 0.0; //CL, Angle in degrees
  double Angy = 0.0; //CL, Angle in degrees
  double Angz = 0.0; //CL, Angle in degrees
  double Rey = 100.0; //CL, Reynolds number
  unsigned Noel = 4; //CL, Number of elements in 1D
  unsigned Noelx = 4; //CL, Number of elements in 1D
  unsigned Noely = 4; //CL, Number of elements in 1D
  unsigned Noelz = 4; //CL, Number of elements in 1D
  double Scaling_sigma = 0; //CL, If the scaling sigma is not set, then
                            // the default is the norm of the momentum block. 
  
  // These are strings for 
  bool Use_axnorm = true; //Set from CL, use norm for sigma?
  bool Use_block_diagonal_w = false; // To set from CL
  bool Doc_prec = false; // To set from CL
  bool Doc_soln = false; // To set from CL
  bool Print_hypre = true;

  std::string Label = ""; // To be set as the label for this problem. Contains
                          // all the information for this run.
  std::string Soln_dir = ""; // Where to put the solution.
  std::string Doc_prec_dir = ""; // Where to put the solution.

  std::string Itstime_dir = ""; //Set from CL, directory to output the 
                                // iteration counts and timing results.

  // Used to determine if we are using the TrilinosAztecOOSolver solver or not.
  // This cannot be determined by the OOMPH_HAS_TRILINOS ifdef since we may be
  // using OOMPH-LIB's GMRES even if we have Trilinos. This should be set in
  // the problem constuctor as soon as we set the linear_solver_pt() for the
  // problem.
  bool Using_trilinos_solver = false;

  double Rey_start = 0.0;
  double Rey_incre = 50.0;
  double Rey_end = 500.0;

  // Object to store the linear solver iterations and times.
  DocLinearSolverInfo* Doc_linear_solver_info_pt = 0;

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

  string create_label();

  inline string create_label()
  {
    std::string prob_str = "";
    std::string w_str = "";
    std::string ns_str = "";
    std::string f_str = "";
    std::string p_str = "";
    std::string vis_str = "";
    std::string ang_str = "";
    std::string rey_str = "";
    std::string noel_str = "";
    std::string w_approx_str = "";
    std::string sigma_str = "";

    switch(Prob_id)
    {
      case 10:
        prob_str = "CuTmp";
        break;
      case 11:
        prob_str = "CuPo";
        break;
      case 12:
        prob_str = "SqTf";
        break;
      case 13:
        prob_str = "SqTfPo";
        break;
      case 20:
        prob_str = "AwTmp";
        break;
      case 21:
        prob_str = "AwPo";
        break;
      case 22:
        prob_str = "AwTf";
        break;
      case 23:
        prob_str = "AwTfPo";
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

     
     
    // Set the string for W_solver.
    switch(W_solver)
    {
      case 0:
        w_str = "We";
        break;
      default:
      {
        std::ostringstream err_msg;
        err_msg << "There is an unrecognised W_solver, recognised W_solver:\n"
                << "0 = (We) SuperLU solve\n"
                << std::endl;
        throw OomphLibError(err_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    } // switch

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
        err_msg << "There is an unrecognised NS_solver, recognised NS_solver:\n"
                << "0 = (Ne) SuperLU\n"
                << "1 = (Nl) LSC preconditioner\n"
                << std::endl;
        throw OomphLibError(err_msg.str(),
                            OOMPH_CURRENT_FUNCTION,
                            OOMPH_EXCEPTION_LOCATION);
      }
    } // switch NS_solver


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

    // Now we continue with setting the string for the solvers.
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
    }  // switch

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
    } // switch

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

    // Set ang_str
    {
      std::ostringstream strs;
      strs << "Ax" << Angx_deg << "y" << Angy_deg << "z" << Angz_deg;
      ang_str = strs.str();
    }

    // Set the Reynolds number.
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

    // Set Noel_str, used for book keeping.
    if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
    {
      std::ostringstream strs;
      strs << "N" <<Noel;
      noel_str = strs.str();
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

    // Use the diagonal or block diagonal approximation for W block.
    if(CommandLineArgs::command_line_flag_has_been_set("--bdw"))
    {
      Use_block_diagonal_w = true;
      w_approx_str = "Wbd";
    }
    else
    {
      Use_block_diagonal_w = false;
      w_approx_str = "";
    }

    // Set Use_axnorm, if sigma has not been set, 
    // norm os momentum block is used.
    if(CommandLineArgs::command_line_flag_has_been_set("--sigma"))
    {
      Use_axnorm = false;
      
      std::ostringstream strs;
      strs << "S" << Scaling_sigma;
      sigma_str = strs.str();
    }

    std::string label = prob_str
                        + w_str + ns_str + f_str + p_str
                        + vis_str + ang_str 
                        + rey_str + noel_str + w_approx_str + sigma_str;
     
    return label; 
  } // inlined function create_label


  // Mesh info:
  unsigned Mesh_type = 0; // SL, 0 = TET, 1 = HEX
  double Mesh_area = 0.0;
  std::string Mesh_folder_str = "";

  /////////// Time stepping stuff.


  bool Impulsive_start_flag = true;

  unsigned Soln_num = 0;

  // The are the time step parameters. 
  double T_min = 0; // Max time.
  double T_max = 0.5; // Min time.
  double Dt = -1; // if this is set, Ntsteps must not be set.
  int Ntsteps = -1; // if this is set, Dt must not be set.


 /// Magnitude of fluid pressure on inflow boundary
 double P_in=1.0;
 // double P_in=0.5;

 /// Period
 double Period = 1.0;

 /// Applied traction on fluid at the inflow boundary
 void prescribed_inflow_traction(const double& t,
                                 const Vector<double>& x,
                                 const Vector<double>& n,
                                 Vector<double>& traction)
 {
  traction[0]= -P_in*(1.0 - cos((2.0*MathematicalConstants::Pi*t)/Period));
  traction[1]=0.0;
  traction[2]=0.0;
 }  
  
///// Traction at the outflow boundary
// void prescribed_traction(const double& t,
//                          const Vector<double>& x,
//                          Vector<double>& traction)
// {
//  traction.resize(3);
//  traction[0]=1.0;
//  traction[1]=0.0;
//  traction[2]=0.0;
// } 

} // End of namespace CubeLagrange


#endif // End of ifndef RAY_LAGRANGE_HEADER

