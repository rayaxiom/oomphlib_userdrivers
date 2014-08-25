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


#ifndef RAY_NAVIER_STOKES_HEADER
#define RAY_NAVIER_STOKES_HEADER

// Oomph-lib includes
#include "generic.h"

using namespace oomph;

namespace NavierStokesHelpers
{

  int Vis = -1;
  double Rey = -1.0;
  double Re_invFr = -1.0; // NOT SET FROM CL

  double Rey_start = -1.0;
  double Rey_incre = -1.0;
  double Rey_end = -1.0;

  double Delta_t = -1.0;
  double Time_start = -1.0;
  double Time_end = -1.0;


  int Dim = -1;

  inline void setup_commandline_flags()
  {
    CommandLineArgs::specify_command_line_flag("--visc", 
        &Vis);
    CommandLineArgs::specify_command_line_flag("--rey", &Rey);


    CommandLineArgs::specify_command_line_flag("--rey_start", &Rey_start);
    CommandLineArgs::specify_command_line_flag("--rey_incre", &Rey_incre);
    CommandLineArgs::specify_command_line_flag("--rey_end", &Rey_end);
  }


  inline void generic_setup()
  {

    /////////////////

    // Setting up the viscous term depends on the dimension of the problem.
    // So we first have to check if this has been set.
    if(Dim == -1)
    {
      std::ostringstream err_msg;
      err_msg << "Please set the dimension Dim in the main function\n"
        << "to either 2 or 3"
        << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    if(!CommandLineArgs::command_line_flag_has_been_set("--visc"))
    {
      std::ostringstream err_msg;
      err_msg << "Please set --visc to 0 or 1" << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {
      // Vis is either 0 or 1.
      if(Vis == 0)
      {
        // Dim is guaranteed to be either 2 or 3
        if(Dim == 2)
        {
          oomph_info << "RAYINFO: setting SIMPLE form for 2D" << std::endl; 
          
          NavierStokesEquations<2>::Gamma[0] = 0.0;
          NavierStokesEquations<2>::Gamma[1] = 0.0;
        }
        else
        {
          oomph_info << "RAYINFO: setting SIMPLE form for 3D" << std::endl; 
          NavierStokesEquations<3>::Gamma[0] = 0.0;
          NavierStokesEquations<3>::Gamma[1] = 0.0;
          NavierStokesEquations<3>::Gamma[2] = 0.0;
        }
      }
      else if(Vis == 1)
      {
        if(Dim == 2)
        {
          oomph_info << "RAYINFO: setting STRESS DIVERGENCE form for 2D" << std::endl; 
          NavierStokesEquations<2>::Gamma[0] = 1.0;
          NavierStokesEquations<2>::Gamma[1] = 1.0;
        }
        else
        {
          oomph_info << "RAYINFO: setting STRESS DIVERGENCE form for 3D" << std::endl; 
          NavierStokesEquations<3>::Gamma[0] = 1.0;
          NavierStokesEquations<3>::Gamma[1] = 1.0;
          NavierStokesEquations<3>::Gamma[2] = 1.0;
        }
      }
      else
      {
        std::ostringstream err_msg;
        err_msg << "Vis is not either 0 or 1" << std::endl;
        throw OomphLibError(err_msg.str(),
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
      }
    }

    // Reynolds number is handled flexibly.
  }

  inline double global_temporal_error_norm(Problem* problem_pt,
      const unsigned& dim,
      Mesh* const & bulk_mesh_pt)
  {
#ifdef OOMPH_HAS_MPI
    double global_error = 0.0;

    // Find out how many nodes there are in the problem.
    unsigned n_node = bulk_mesh_pt->nnode();

    // Loop over the nodes and calculate the estimate error in the values
    // for non-haloes
    int count = 0;
    for (unsigned nod_i = 0; nod_i < n_node; nod_i++) 
    {
      Node* nod_pt = bulk_mesh_pt->node_pt(nod_i);
      if(!(nod_pt->is_halo()))
      {
        // Get the error in solution: Difference between the predicted and
        // actual value for nodal values 0, ... dim-1
        double node_error = 0.0;
        for (unsigned nodal_i = 0; nodal_i < dim; nodal_i++) 
        {
          double error = nod_pt->time_stepper_pt()->
            temporal_error_in_value(nod_pt,nodal_i);
          node_error += error*error;
          count++;
        } // for loop over the nodal values 0,..,dim-1
        global_error += node_error;
      } // if the node is not halo
    } // for loop over nodes

    // Accumulate
    int n_node_local = count;
    int n_node_total = 0;

    MPI_Allreduce(&n_node_local,&n_node_total,1,MPI_INT,MPI_SUM,
        problem_pt->communicator_pt()->mpi_comm());

    double global_error_total = 0.0;
    MPI_Allreduce(&global_error,&global_error_total,1,MPI_DOUBLE,MPI_SUM,
        problem_pt->communicator_pt()->mpi_comm());

    // Divide by the number of nodes
    global_error_total /= double(n_node_total);

    // Return square root...
    // RAYRAY remove this output?
    oomph_info << "Total error " << n_node_total << " " 
      <<  sqrt(global_error_total) << std::endl;
    return sqrt(global_error_total);
#else
    double global_error = 0.0;

    // Find out how many nodes there are in the problem
    unsigned n_node = bulk_mesh_pt->nnode();

    // Loop over the nodes and calculate the errors in the values.
    for (unsigned node_i = 0; node_i < n_node; node_i++) 
    {
      Node* nod_pt = bulk_mesh_pt->node_pt(node_i);

      // Get the error in the solution: Difference between predicted and 
      // actual value for nodal values 0,..,dim-1
      double nodal_error = 0.0;
      for (unsigned nodal_i = 0; nodal_i < dim; nodal_i++) 
      {
        double error = nod_pt->time_stepper_pt()->
          temporal_error_in_value(nod_pt,nodal_i);

        // Add the squared value to the nodal error.
        nodal_error += error*error;
      }

      // add the nodal error to the global error.
      global_error += nodal_error;
    }

    global_error /= double(n_node*3);

    // Return square root
    return sqrt(global_error);
#endif
  } // EoF global_temporal_error_norm


  inline std::string create_label()
  {
    std::string label = "";

    std::ostringstream nsp_stream;
    if(Vis == 0)
    {
      nsp_stream << "Sim";
    }
    else if (Vis == 1)
    {
      nsp_stream << "Str";
    }
    else
    {
      std::ostringstream err_msg;
      err_msg << "No such Vis string"   << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

    if(Rey < 0)
    {
      std::ostringstream err_msg;
      err_msg << "Rey is negative."   << std::endl;
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }
    else
    {

      nsp_stream << "R" << Rey;
    }

    return nsp_stream.str();


  }


}





#endif



