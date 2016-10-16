//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision$
//LIC//
//LIC// $LastChangedDate$
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
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
// Driver code for a simple unstructured fluid problem using a mesh
// generated from an input file generated by the 3d mesh generator
// tetgen

#include <fenv.h>
#include <sstream>
#include <iomanip>
#include <ios>


//Generic routines
#include "generic.h"
#include "constitutive.h"
#include "navier_stokes.h"

// Get the mesh
#include "meshes/tetgen_mesh.h"
#include "meshes/brick_from_tet_mesh.h"

using namespace std;
using namespace oomph;


//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////


//=======start_namespace==========================================
/// Global variables
//================================================================
namespace Global_Parameters
{

 /// Fluid pressure on inflow boundary
 double P_in=0.5;

 /// Applied traction on fluid at the inflow boundary
 void prescribed_inflow_traction(const double& t,
                                 const Vector<double>& x,
                                 const Vector<double>& n,
                                 Vector<double>& traction)
 {
  traction[0]=0.0;
  traction[1]=0.0;
  traction[2]=P_in;
 } 


 /// Fluid pressure on outflow boundary
 double P_out=-0.5; 

 /// Applied traction on fluid at the inflow boundary
 void prescribed_outflow_traction(const double& t,
                                  const Vector<double>& x,
                                  const Vector<double>& n,
                                  Vector<double>& traction)
 {
  traction[0]=0.0;
  traction[1]=0.0;
  traction[2]=-P_out;
 } 
 
 // These are the pre-sets (defaults) for problem parameters.

 /// Problem Dimension
 const unsigned Dim = 3;

 /// Switch for mesh: --use_brick
 bool Use_brick = false;

 /// Switch for linear solver type: --use_iterative_lin_solver
 bool Use_iterative_lin_solver = false;

 /// Switch for iterative linear solver: --use_trilinos
 bool Use_trilinos = false;

 /// Prec. for Navier-Stokes block: --use_lsc
// bool Use_lsc = false; we always use lsc as top level prec

 /// Prec. for velocity block: --use_amg_for_f
 bool Use_amg_for_f = false;

 /// Use Boomer AMG for the pressure block?
 bool Use_amg_for_p = false;

 /// Using stress divergence viscous term?
 bool Use_stress_div = true;

 /// Reynolds number
 double Re = 100.0;

 /// Doc number (for doc_solution)
 unsigned Doc_num = 0;

 /// Doc directory (for doc_solution)
 std::string Doc_dir = "RESLT";

 /// Label for doc solution
 std::string Doc_label = "fluid_soln";

 /// Label for tetgen file (incl folder)
 std::string Tetgen_label = "tetgen_original/fsi_bifurcation_fluid";

 /// Number for tetgen file
 unsigned Tetgen_num = 1;

 /// Storage for number of iterations during Newton steps 
 Vector<unsigned> Iterations;

} //end namespace


namespace DriverCodeHelpers
{
inline void specify_command_line_flag_helper()
{
  // Alias the namespace for convenience.
  namespace GP = Global_Parameters;

  // Use brick mesh?
  CommandLineArgs::specify_command_line_flag("--use_brick");

  // Use iterative linear solver?
  CommandLineArgs::specify_command_line_flag("--use_iterative_lin_solver");

  // Use trilinos GMRES?
  CommandLineArgs::specify_command_line_flag("--use_trilinos");


  // Use lsc solver?
//  CommandLineArgs::specify_command_line_flag("--use_lsc");

  // Use AMG for F block?
  CommandLineArgs::specify_command_line_flag("--use_amg_for_f");

  // Use AMG for P block?
  CommandLineArgs::specify_command_line_flag("--use_amg_for_p");

  // Use stress divergence viscous term?
  CommandLineArgs::specify_command_line_flag("--use_stress_div");

  // Reynolds number.
  CommandLineArgs::specify_command_line_flag("--re",&GP::Re);

  // Setting for Doc_info
  CommandLineArgs::specify_command_line_flag("--doc_dir",&GP::Doc_dir);
  CommandLineArgs::specify_command_line_flag("--doc_num",&GP::Doc_num);
  CommandLineArgs::specify_command_line_flag("--doc_label",&GP::Doc_label);
 
  // Label for tetgen file
  CommandLineArgs::specify_command_line_flag("--tetgen_label",
                                             &GP::Tetgen_label);

  // Number for tetgen file
  CommandLineArgs::specify_command_line_flag("--tetgen_num",
                                             &GP::Tetgen_num);
}

inline void setup_command_line_flags(DocInfo& doc_info)
{
  namespace GP = Global_Parameters;

  if(CommandLineArgs::command_line_flag_has_been_set("--use_brick"))
  {
    GP::Use_brick = true;
  }
  else
  {
    GP::Use_brick = false;
  }

  if(CommandLineArgs::command_line_flag_has_been_set(
        "--use_iterative_lin_solver"))
  {
    GP::Use_iterative_lin_solver = true;
  }
  else
  {
    GP::Use_iterative_lin_solver = false;
  }

  // Set the flag for trilinos solver
  if(CommandLineArgs::command_line_flag_has_been_set("--use_trilinos"))
  {
#ifdef OOMPH_HAS_TRILINOS
    GP::Use_trilinos = true;
#else
    std::ostringstream warning_stream;
    warning_stream << "WARNING: \n"
      << "No trilinos installed, using oomphlib's GMRES solver.\n";
    OomphLibWarning(warning_stream.str(),
      OOMPH_CURRENT_FUNCTION,
      OOMPH_EXCEPTION_LOCATION);

    GP::Use_trilinos = false;
#endif
  }
  else
  {
    GP::Use_trilinos = false;
  }

//  if(CommandLineArgs::command_line_flag_has_been_set("--use_lsc"))
//  {
//    GP::Use_lsc = true;
//  }
//  else
//  {
//    GP::Use_lsc = false;
//  }

  // Set the flag for amg for the momentum block
  if(CommandLineArgs::command_line_flag_has_been_set("--use_amg_for_f"))
  {
    GP::Use_amg_for_f = true;
  }
  else
  {
    GP::Use_amg_for_f = false;
  }

  // Set the flag for amg for the pressure block
  if(CommandLineArgs::command_line_flag_has_been_set("--use_amg_for_p"))
  {
    GP::Use_amg_for_p = true;
  }
  else
  {
    GP::Use_amg_for_f = false;
  }

  if(CommandLineArgs::command_line_flag_has_been_set("--use_stress_div"))
  {
    oomph_info << "Using stress divergence viscous term" << std::endl; 
    for (unsigned d = 0; d < GP::Dim; d++)
    {
      NavierStokesEquations<GP::Dim>::Gamma[d] = 1.0;
    }
  }
  else
  {
    oomph_info << "Using simple viscous term" << std::endl; 
    for (unsigned d = 0; d < GP::Dim; d++)
    {
      NavierStokesEquations<GP::Dim>::Gamma[d] = 0.0;
    }
  }
 
  if(!CommandLineArgs::command_line_flag_has_been_set("--re"))
  {
    oomph_info << "--re has not been set. Using default Re=" 
               << GP::Re << std::endl; 
  }
  
  
  if(!CommandLineArgs::command_line_flag_has_been_set("--doc_dir"))
  {
    oomph_info << "--doc_dir has not been set. Using default Doc_dir=" 
               << GP::Doc_dir << std::endl;
  }
  doc_info.set_directory(GP::Doc_dir);

  if(!CommandLineArgs::command_line_flag_has_been_set("--doc_num"))
  {
    oomph_info << "--doc_num has not been set. Using default Doc_num=" 
               << GP::Doc_num << std::endl; 
  }
  doc_info.number()=GP::Doc_num;

  if(!CommandLineArgs::command_line_flag_has_been_set("--doc_label"))
  {
    oomph_info << "--doc_label has not been set. Using default Doc_label=" 
               << GP::Doc_label << std::endl; 
  }
  doc_info.label()=GP::Doc_label;

  
  if(!CommandLineArgs::command_line_flag_has_been_set("--tetgen_label"))
  {
    oomph_info 
      << "--tetgen_num has not been set. Using default Tetgen_label=" 
      << GP::Tetgen_label << std::endl; 
  }

  if(!CommandLineArgs::command_line_flag_has_been_set("--tetgen_num"))
  {
    oomph_info << "--tetgen_num has not been set. Using default Tetgen_num=" 
               << GP::Tetgen_num << std::endl; 
  }
  else
  {
    if((GP::Tetgen_num < 1)||(GP::Tetgen_num > 13))
    {
      std::ostringstream err_msg;
      err_msg << "Tetgen_num=" << GP::Tetgen_num << "\n"
              << "Must be between 1 and 13, use --tetgen_num\n";
      throw OomphLibError(err_msg.str(),
          OOMPH_CURRENT_FUNCTION,
          OOMPH_EXCEPTION_LOCATION);
    }

  }

} // setup_command_line_flags()


} // namespace DriverCodeHelpers




//======start_problem_class===========================================
/// Unstructured fluid problem
//====================================================================
template<class ELEMENT>
class UnstructuredFluidProblem : public Problem
{

public:

 /// Constructor: 
 UnstructuredFluidProblem();

 /// Destructor (empty)
 ~UnstructuredFluidProblem(){}


 /// \short Update before Newton solve.
 void actions_before_newton_solve()
 {
   if(Global_Parameters::Use_iterative_lin_solver)
   {
     Global_Parameters::Iterations.clear();
   }
 }

 /// \short Update after Newton step - document the number of iterations 
 /// required for the iterative solver to converge.
 void actions_after_newton_step()
 {
   // Get the iteration counts if using iterative linear solver
   if(Global_Parameters::Use_iterative_lin_solver)
   {
     unsigned iter = static_cast<IterativeLinearSolver*>
      (this->linear_solver_pt())->iterations();

     Global_Parameters::Iterations.push_back(iter);
   }
 }


 /// Doc the solution
 void doc_solution(DocInfo& doc_info);
 
 /// Return total number of fluid inflow traction boundaries
 unsigned nfluid_inflow_traction_boundary()
  {
   return Inflow_boundary_id.size();
  }

 /// Return total number of fluid outflow traction boundaries
 unsigned nfluid_outflow_traction_boundary()
  {
   return Outflow_boundary_id.size();
  }

 /// Return total number of fluid outflow traction boundaries
 unsigned nfluid_traction_boundary()
  {
   return Inflow_boundary_id.size()+Outflow_boundary_id.size();
  }

 //private:

 /// Create fluid traction elements at inflow
 void create_fluid_traction_elements();

 /// Bulk fluid mesh
 Mesh* Fluid_mesh_pt;

 /// Meshes of fluid traction elements that apply pressure at in/outflow
 Vector<Mesh*> Fluid_traction_mesh_pt;

 /// \short IDs of fluid mesh boundaries along which inflow boundary conditions
 /// are applied
 Vector<unsigned> Inflow_boundary_id;

 /// \short IDs of fluid mesh boundaries along which inflow boundary conditions
 /// are applied
 Vector<unsigned> Outflow_boundary_id;

 /// Preconditioner (master LSC preconditioner)
 Preconditioner* Prec_pt;

 /// Preconditioner for the momentum block
 Preconditioner* F_preconditioner_pt;

 /// Preconditioner for the pressure block
 Preconditioner* P_preconditioner_pt;

 /// Iterative linear solver
 IterativeLinearSolver* Solver_pt;

};



//==========start_constructor=============================================
/// Constructor for unstructured 3D fluid problem
//========================================================================
template<class ELEMENT>
UnstructuredFluidProblem<ELEMENT>::UnstructuredFluidProblem()
{ 
 
 // Set up the tetgen files.
 std::stringstream ss;
 ss<<Global_Parameters::Tetgen_num;

 //Create fluid bulk mesh, sub-dividing "corner" elements
 string node_file_name=Global_Parameters::Tetgen_label
                       +"."+ss.str()+".node";
 string element_file_name=Global_Parameters::Tetgen_label
                          +"."+ss.str()+".ele";
 string face_file_name=Global_Parameters::Tetgen_label
                       +"."+ss.str()+".face";
 bool split_corner_elements=true;
 if(Global_Parameters::Use_brick)
 {
   Fluid_mesh_pt = new BrickFromTetMesh<ELEMENT>(node_file_name,
                                                 element_file_name,
                                                 face_file_name,
                                                 split_corner_elements);
 }
 else
 {
   Fluid_mesh_pt =  new TetgenMesh<ELEMENT>(node_file_name,
                                            element_file_name,
                                            face_file_name,
                                            split_corner_elements);
 }

 // Find elements next to boundaries
 Fluid_mesh_pt->setup_boundary_element_info();

 // The following corresponds to the boundaries as specified by
 // facets in the tetgen input:

 // Fluid mesh has one inflow boundary: Boundary 0
 Inflow_boundary_id.resize(1);
 Inflow_boundary_id[0]=0;
 
 // Fluid mesh has two outflow boundaries: Boundaries 1 and 2
 Outflow_boundary_id.resize(2);
 Outflow_boundary_id[0]=1;
 Outflow_boundary_id[1]=2;
 
 // Apply BCs
 //----------
 
 // Map to indicate which boundary has been done
 std::map<unsigned,bool> done; 
  
 // Loop over inflow/outflow boundaries to impose parallel flow
 for (unsigned in_out=0;in_out<2;in_out++)
  {
   // Loop over in/outflow boundaries
   unsigned n=nfluid_inflow_traction_boundary();
   if (in_out==1) n=nfluid_outflow_traction_boundary();
   for (unsigned i=0;i<n;i++)
    {
     // Get boundary ID
     unsigned b=0;
     if (in_out==0)
      {
       b=Inflow_boundary_id[i];
      }
     else
      {
       b=Outflow_boundary_id[i];
      }

     // Number of nodes on that boundary
     unsigned num_nod=Fluid_mesh_pt->nboundary_node(b);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get the node
       Node* nod_pt=Fluid_mesh_pt->boundary_node_pt(b,inod);
       
       // Pin transverse velocities
       nod_pt->pin(0);
       nod_pt->pin(1);
      }
     
     // Done!
     done[b]=true;
    }

  } // done in and outflow
 
 
 
 // Loop over all fluid mesh boundaries and pin velocities
 // of nodes that haven't been dealt with yet
 unsigned nbound=Fluid_mesh_pt->nboundary();
 for(unsigned b=0;b<nbound;b++)
  {

   // Has the boundary been done yet?
   if (!done[b])
    {
     unsigned num_nod=Fluid_mesh_pt->nboundary_node(b);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       // Get node
       Node* nod_pt= Fluid_mesh_pt->boundary_node_pt(b,inod);
       
       // Pin all velocities
       nod_pt->pin(0); 
       nod_pt->pin(1); 
       nod_pt->pin(2); 
      }
    }

  } // done no slip elsewhere 
 
 
 // Complete the build of the fluid elements so they are fully functional
 //----------------------------------------------------------------------
 unsigned n_element = Fluid_mesh_pt->nelement();
 for(unsigned e=0;e<n_element;e++)
  {

   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Fluid_mesh_pt->element_pt(e));
   
   //Set the Reynolds number
   el_pt->re_pt() = &Global_Parameters::Re;   

  } 
 
 
 // Create meshes of fluid traction elements at inflow/outflow
 //-----------------------------------------------------------
 
 // Create the meshes
 unsigned n=nfluid_traction_boundary();
 Fluid_traction_mesh_pt.resize(n);
 for (unsigned i=0;i<n;i++)
  {
   Fluid_traction_mesh_pt[i]=new Mesh;
  } 
 
 // Populate them with elements
 create_fluid_traction_elements();
 
 
 // Combine the lot
 //----------------
 
 // Add sub meshes:

 // Fluid bulk mesh
 add_sub_mesh(Fluid_mesh_pt);
 
 // The fluid traction meshes
 n=nfluid_traction_boundary();
 for (unsigned i=0;i<n;i++)
  { 
   add_sub_mesh(Fluid_traction_mesh_pt[i]);
  }
 
 // Build global mesh
 build_global_mesh();

 // Setup equation numbering scheme
 std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 

 ///////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////

 // Now we create and setup the solvers.
 
 // Are we using an iterative linear solver?
 if(Global_Parameters::Use_iterative_lin_solver)
 {
  // We choose either trilinos or oomph-lib's GMRES as our linear solver
  if(Global_Parameters::Use_trilinos)
  {
#ifdef OOMPH_HAS_TRILINOS
    // Create the trilinos solver.
    TrilinosAztecOOSolver* trilinos_solver_pt = new TrilinosAztecOOSolver;
    trilinos_solver_pt->solver_type() = TrilinosAztecOOSolver::GMRES;

    // Store the solver pointer.
    Solver_pt = trilinos_solver_pt;
#endif
  }
  else
  {
    // Create oomph-lib iterative linear solver.
    IterativeLinearSolver* solver_pt = new GMRES<CRDoubleMatrix>;
    
    // We use RHS preconditioning. Note that by default,
    // left hand preconditioning is used.
    static_cast<GMRES<CRDoubleMatrix>*>(solver_pt)
      ->set_preconditioner_RHS();

    // Store the solver pointer.
    Solver_pt = solver_pt;
  }

  // Set tolerance
  Solver_pt->tolerance() = 1.0e-6;

 
 // Create the preconditioner
 NavierStokesSchurComplementPreconditioner* lsc_prec_pt = 
   new NavierStokesSchurComplementPreconditioner(this);
 lsc_prec_pt->set_navier_stokes_mesh(Fluid_mesh_pt);
 lsc_prec_pt->use_lsc();
 
  // Use AMG for f block?
  if(Global_Parameters::Use_amg_for_f)
   {
     F_preconditioner_pt
       = Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper::
         boomer_amg_for_3D_momentum();
       lsc_prec_pt->set_f_preconditioner(F_preconditioner_pt);
   }
 
   // Use AMG for p block?
   if(Global_Parameters::Use_amg_for_p)
   {
     P_preconditioner_pt
       = Lagrange_Enforced_Flow_Preconditioner_Subsidiary_Operator_Helper::
         boomer_amg_for_3D_poisson_problem();
     lsc_prec_pt->set_p_preconditioner(P_preconditioner_pt);
   }

  // Store the preconditioner pointers.
  Prec_pt = lsc_prec_pt;

  // Pass the preconditioner to the solver.
  Solver_pt->preconditioner_pt() = lsc_prec_pt;

  // Max linear solver iterations.
  Solver_pt->max_iter() = 200;

  // Pass the solver to the problem.
  this->linear_solver_pt() = Solver_pt;

  // Set the Newton solver tolerance.
  this->newton_solver_tolerance() = 1.0e-6;

 }
 
} // end constructor



//============start_of_fluid_traction_elements==============================
/// Create fluid traction elements 
//=======================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::create_fluid_traction_elements()
{

 // Counter for number of fluid traction meshes
 unsigned count=0;

 // Loop over inflow/outflow boundaries
 for (unsigned in_out=0;in_out<2;in_out++)
  {
   // Loop over boundaries with fluid traction elements
   unsigned n=nfluid_inflow_traction_boundary();
   if (in_out==1) n=nfluid_outflow_traction_boundary();
   for (unsigned i=0;i<n;i++)
    {
     // Get boundary ID
     unsigned b=0;
     if (in_out==0)
      {
       b=Inflow_boundary_id[i];
      }
     else
      {
       b=Outflow_boundary_id[i];
      }
     
     // How many bulk elements are adjacent to boundary b?
     unsigned n_element = Fluid_mesh_pt->nboundary_element(b);
     
     // Loop over the bulk elements adjacent to boundary b
     for(unsigned e=0;e<n_element;e++)
      {
       // Get pointer to the bulk element that is adjacent to boundary b
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
        Fluid_mesh_pt->boundary_element_pt(b,e));
       
       //What is the index of the face of the element e along boundary b
       int face_index = Fluid_mesh_pt->face_index_at_boundary(b,e);
       
       // Create new element 
       NavierStokesTractionElement<ELEMENT>* el_pt=
        new NavierStokesTractionElement<ELEMENT>(bulk_elem_pt,
                                                       face_index);
       
       // Add it to the mesh
       Fluid_traction_mesh_pt[count]->add_element_pt(el_pt);
       
       // Set the pointer to the prescribed traction function
       if (in_out==0)
        {
         el_pt->traction_fct_pt() = 
          &Global_Parameters::prescribed_inflow_traction;
        }
       else
        {
         el_pt->traction_fct_pt() = 
          &Global_Parameters::prescribed_outflow_traction;
        }
      }
     // Bump up counter
     count++;
    }
  }
 
 } // end of create_traction_elements



//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void UnstructuredFluidProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{ 

 // Create the file name string.
 std::stringstream filename;
 filename << doc_info.directory() 
           << "/"
          <<doc_info.label()<<doc_info.number()<<".dat";

 // Number of plot points
 unsigned npts=5;

 // Output solution
 std::ofstream some_file;
 some_file.open(filename.str().c_str());
 Fluid_mesh_pt->output(some_file,npts);
 some_file.close();
 
}


void print_avg_iter(const Vector<unsigned> * iters_pt,
                    std::ostringstream* results_stream_pt)
{
  // 
  const unsigned nnewtonstep = iters_pt->size();
  unsigned total_its = 0;

  for(unsigned istep = 0; istep < nnewtonstep; istep++)
  {
    total_its += (*iters_pt)[istep];
  }

  double average_its = ((double)total_its)
    / ((double)nnewtonstep);

  (*results_stream_pt) << "RAYAVGITS:\t";
  // Print to one decimal place if the average is not an exact
  // integer. Otherwise we print normally.
  std::streamsize tmp_precision = results_stream_pt->precision();
  (*results_stream_pt) << "\t" << std::fixed << std::setprecision(1)
    << average_its << "(" << nnewtonstep << ")" << "\n";

  // reset the precision
  (*results_stream_pt) << std::setprecision(tmp_precision);
}



//=============start_main=================================================
/// Demonstrate how to solve an unstructured 3D fluids problem
//========================================================================
int main(int argc, char **argv)
{
#ifdef OOMPH_HAS_MPI
  MPI_Helpers::init(argc,argv);
#endif

  // Store command line arguments
  CommandLineArgs::setup(argc,argv);

  DriverCodeHelpers::specify_command_line_flag_helper();

  // Parse the above flags.
  CommandLineArgs::parse_and_assign();
  CommandLineArgs::doc_specified_flags();

  // Label for output
  DocInfo doc_info;

  // Set up flags
  DriverCodeHelpers::setup_command_line_flags(doc_info);


  if(Global_Parameters::Use_brick)
  {
    //Set up the problem
    UnstructuredFluidProblem<QTaylorHoodElement<3> > problem;
    // Solve the problem
    problem.newton_solve();

    std::ostringstream results_stream;
    print_avg_iter(&Global_Parameters::Iterations,
        &results_stream);

    // Create an out file.
    // The output file.
    std::ofstream outfile;

    // If we want to output to a file, we create the outfile.
    std::ostringstream filename_stream;
    filename_stream <<"res_iterations/iter"
      << Global_Parameters::Doc_num;
    outfile.open(filename_stream.str().c_str());

    outfile << "\n" << results_stream.str();
    outfile.close();


  }
  else
  {
    //Set up the problem
    UnstructuredFluidProblem<TTaylorHoodElement<3> > problem;

    //Output initial guess
    // problem.doc_solution(doc_info);
    // doc_info.number()++;

    // Solve the problem
    problem.newton_solve();

    std::ostringstream results_stream;
    print_avg_iter(&Global_Parameters::Iterations,
        &results_stream);

    // Create an out file.
    // The output file.
    std::ofstream outfile;

    // If we want to output to a file, we create the outfile.
    std::ostringstream filename_stream;
    filename_stream <<"res_iterations/iter"
      << Global_Parameters::Doc_num;
    outfile.open(filename_stream.str().c_str());

    outfile << "\n" << results_stream.str();
    outfile.close();
  }

#ifdef OOMPH_HAS_MPI
  MPI_Helpers::finalize();
#endif  


} // end_of_main




