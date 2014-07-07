// Working parallel outflow.
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

#include <sstream>
#include <iomanip>
#include <ios>

// Oomph-lib includes
#include "generic.h"
#include "navier_stokes.h"
//#include "raymon.h"

// The 2D mesh
#include "meshes/quarter_circle_sector_mesh.h"

//#include "meshes/simple_rectangular_quadmesh.h"
//#include "meshes/rectangular_quadmesh.h"

//using namespace std;
using namespace oomph;

namespace rayvars
{
  /// Functional body force
//  void body_force(const double& time, const Vector<double>& x, 
//                  Vector<double>& result)
//  {
//    result[0]=0.0;
//    result[1]=-Re_invFr;
//  }

  /// Zero functional body force
  void zero_body_force(const double& time, const Vector<double>& x, 
                       Vector<double>& result)
  {
   result[0]=0.0;
   result[1]=0.0;
  }
}


struct SquareLagrangeVariables{

  unsigned Vis; // CL
  double Rey; // CL
  double Re_invFr; // preset.
  Vector<double> Gravity;
  unsigned Noel; // CL


  std::string Prob_str; // To set, no CL

  std::string Vis_str; // CL
  std::string Rey_str; // To set from CL
  std::string Noel_str; // To set from CL

  bool Loop_reynolds;

  // Setting the defaults:
  SquareLagrangeVariables() :
    Vis(0), // 0 - Simple, 1 - Stress divergence.
    Rey(100.0), // Reynolds number
    Re_invFr(100.0), // Reynolds/Froude number
    Gravity(2),
    Noel(4), // Number of elements.
    Prob_str("Qcirc2D"), // This is set inside the code, not from commandline.
    Vis_str("Sim"), // Sim - Simple, Str = Stress Divergence
    Rey_str("R100"), // Reynolds number.
    Noel_str("N4"), // Number of elements.
    Loop_reynolds(false)
  {}


 
};


namespace oomph
{

//===start_of_problem_class=============================================
//======================================================================

template<class ELEMENT>
class QuarterCircleDrivenCavityProblem : public Problem
{
public:

 /// \short Constructor: Pass number of elements in x and y directions and
 /// lengths
 QuarterCircleDrivenCavityProblem(
   NavierStokesEquations<2>::NavierStokesBodyForceFctPt,
   SquareLagrangeVariables&,
   DocInfo&);

 /// Destructor: Empty
 ~QuarterCircleDrivenCavityProblem() {}


 /// Update before solve is empty
 void actions_before_newton_solve()
 {
   // Inflow: left boundary
   unsigned ibound=2; 
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
   {
     // Parabolic inflow along the y axis
     Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
     
     double y=nod_pt->x(1);
     
     double u = y*(2.0-y);
     double v = 0.0;
     
     nod_pt->set_value(0,u);
     nod_pt->set_value(1,v);
   }

   // Outflow: bottom boundary
   ibound = 0;
   num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
   {
    
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
    
    double x=nod_pt->x(0);
    
    double u = 0.0;
    double v = -x*(2.0-x);
    
    nod_pt->set_value(0,u);
    nod_pt->set_value(1,v);
   }

   /*
   // Set by lagrange multiplier
   ibound = 1;
   num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
   {
    
    Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
    
    double u = 0.0;
    double v = 0.0;
    
    nod_pt->set_value(0,u);
    nod_pt->set_value(1,v);
   }
   */
 } // actions_before_newton_solve

 /// \short Update after solve is empty
 void actions_after_newton_solve() {}

 void actions_after_newton_step()
 {
 }

 /// Actions before adapt: Wipe the mesh of prescribed flux elements
 void actions_before_adapt();

 /// After adaptation: Unpin pressure and pin redudant pressure dofs.
 void actions_after_adapt()
  {
   create_impenetrable_lagrange_elements(1,Bulk_mesh_pt,Surface_mesh1_pt);
   
   rebuild_global_mesh();
   
   // Unpin all pressure dofs
   RefineableNavierStokesEquations<2>::
    unpin_all_pressure_dofs(Bulk_mesh_pt->element_pt());

   // Pin redundant pressure dofs
   RefineableNavierStokesEquations<2>::
    pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
   
   // Now pin the first pressure dof in the first element and set it to 0.0
   //fix_pressure(0,0,0.0);
  } // end_of_actions_after_adapt

 /// Doc the solution
 void doc_solution(DocInfo&);


 void create_impenetrable_lagrange_elements(const unsigned &b,
                                            Mesh* const &bulk_mesh_pt,
                                            Mesh* const &surface_mesh_pt);
 
 /// \short Delete Poisson flux elements and wipe the surface mesh
 void delete_impenetrable_lagrange_elements(Mesh* const &surface_mesh_pt);

private:

 /// Pointer to body force function
 NavierStokesEquations<2>::NavierStokesBodyForceFctPt Body_force_fct_pt;

 /// Fix pressure in element e at pressure dof pdof and set to pvalue
 void fix_pressure(const unsigned &e, const unsigned &pdof, 
                   const double &pvalue)
  {
   //Cast to proper element and fix pressure
   dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e))->
                          fix_pressure(pdof,pvalue);
  } // end_of_fix_pressure


 /// Pointer to the "bulk" mesh
 RefineableQuarterCircleSectorMesh<ELEMENT>* Bulk_mesh_pt;
 
 /// Pointer to the "surface" mesh 1
 Mesh* Surface_mesh1_pt;
 

 /// Pointer to the "surface" mesh
 //Mesh* Surface_mesh_T_pt;
 //Mesh* Surface_mesh_P_pt;

 // Preconditioner
 Preconditioner* Prec_pt;
 // Solver
 IterativeLinearSolver* Solver_pt;

 DocInfo* Doc_info_pt;
};

} // end of namespace oomph


//===start_of_constructor=============================================
/// Problem constructor
//====================================================================
//template<class ELEMENT> // rrrback - changed here.
//TiltedCavityProblem<ELEMENT>::TiltedCavityProblem
//(SquareLagrangeVariables &myvar, DocInfo &doc_info)
//==start_of_constructor==================================================
/// Constructor for driven cavity problem in quarter circle domain
//========================================================================
template<class ELEMENT>
QuarterCircleDrivenCavityProblem<ELEMENT>::QuarterCircleDrivenCavityProblem(
 NavierStokesEquations<2>::NavierStokesBodyForceFctPt body_force_fct_pt,
 SquareLagrangeVariables &myvar, DocInfo &doc_info) :
 Body_force_fct_pt(body_force_fct_pt)
{
 Doc_info_pt = &doc_info;
 
 // Build geometric object that parametrises the curved boundary
 // of the domain

 // Half axes for ellipse
 double a_ellipse=1.0;
 double b_ellipse=1.0;

 // Setup elliptical ring 
 GeomObject* Wall_pt=new Ellipse(a_ellipse,b_ellipse); 

 // End points for wall
 double xi_lo=0.0;
 double xi_hi=2.0*atan(1.0);

 //Now create the mesh
 double fract_mid=0.5;

  Bulk_mesh_pt = new 
  RefineableQuarterCircleSectorMesh<ELEMENT>(
   Wall_pt,xi_lo,fract_mid,xi_hi);

 Surface_mesh1_pt = new Mesh;
 
 create_impenetrable_lagrange_elements(1,Bulk_mesh_pt,Surface_mesh1_pt);
 
 add_sub_mesh(Bulk_mesh_pt);
 add_sub_mesh(Surface_mesh1_pt);
 build_global_mesh();

 // Set error estimator
  Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 dynamic_cast<RefineableQuarterCircleSectorMesh<ELEMENT>*>(
  Bulk_mesh_pt)->spatial_error_estimator_pt()=error_estimator_pt;

 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here: All boundaries are Dirichlet boundaries.
 //unsigned num_bound = Bulk_mesh_pt->nboundary();
 //for(unsigned ibound=0;ibound<num_bound;ibound++)
  {
   unsigned ibound=0;
   unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
   
   ibound=2;
   num_nod= Bulk_mesh_pt->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     // Loop over values (u and v velocities)
     for (unsigned i=0;i<2;i++)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(i); 
      }
    }
  
  } // end loop over boundaries

 //Find number of elements in mesh
 unsigned n_element = Bulk_mesh_pt->nelement();

 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to Reynolds
 // number
 for(unsigned e=0;e<n_element;e++)
  {
   // Upcast from GeneralisedElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));

   //Set the Reynolds number, etc
   el_pt->re_pt() = &myvar.Rey;
   //Set the Re/Fr
   el_pt->re_invfr_pt() = &myvar.Rey;
   //Set Gravity vector
   el_pt->g_pt() = &myvar.Gravity;
   //set body force function
   el_pt->body_force_fct_pt() = Body_force_fct_pt;

  } // end loop over elements

  std::cout << "End of loop over elements" << std::endl; 
  
 // Initial refinement level
 for (unsigned i = 0; i < myvar.Noel; i++)
 {
   refine_uniformly();
 }

 // Pin redudant pressure dofs
 RefineableNavierStokesEquations<2>::
  pin_redundant_nodal_pressures(Bulk_mesh_pt->element_pt());
 
 // Now pin the first pressure dof in the first element and set it to 0.0
 //fix_pressure(0,0,0.0);
 
 // Setup equation numbering scheme
 std::cout <<"Number of equations: " << assign_eqn_numbers() << std::endl; 


 this->newton_solver_tolerance() = 1.0e-6;
 // We use RHS preconditioning. Note that by default,
 // left hand preconditioning is used.

 // Set solver and preconditioner
}

//=====================start_of_actions_before_adapt======================
/// Actions before adapt: Wipe the mesh of prescribed flux elements
//========================================================================
template<class ELEMENT>
void QuarterCircleDrivenCavityProblem<ELEMENT>::actions_before_adapt()
{
 // Kill the flux elements and wipe surface mesh
 delete_impenetrable_lagrange_elements(Surface_mesh1_pt);
 
 // Rebuild the Problem's global mesh from its various sub-meshes
 rebuild_global_mesh();
}// end of actions_before_adapt

//============start_of_delete_flux_elements==============================
/// Delete Poisson Flux Elements and wipe the surface mesh
//=======================================================================
template<class ELEMENT>
void QuarterCircleDrivenCavityProblem<ELEMENT>::
delete_impenetrable_lagrange_elements(Mesh* const &surface_mesh_pt)
{
 // How many surface elements are in the surface mesh
 unsigned n_element = surface_mesh_pt->nelement();
 
 std::cout << "No elements on surface mesh: " << n_element << std::endl; 
 // Loop over the surface elements
 for(unsigned e=0;e<n_element;e++)
  {
   // Kill surface element
   delete surface_mesh_pt->element_pt(e);
  }

 // Wipe the mesh
 surface_mesh_pt->flush_element_and_node_storage();
} // end of delete_flux_elements


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
//template<class ELEMENT>
//void TiltedCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
template<class ELEMENT>
void QuarterCircleDrivenCavityProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{
//  if(doc_info.is_doc_solution_enabled())
  {
    std::ofstream some_file;

    std::stringstream filename;

    filename << "RESLT2/soln_old.dat";


    // Number of plot points
    unsigned npts=5;


    some_file.open(filename.str().c_str());
    Bulk_mesh_pt->output(some_file,npts);
    some_file.close();
  } // if
}




//============start_of_create_parall_outflow_lagrange_elements===========
/// Create ImposeParallelOutflowElement on the b-th boundary of the
/// Mesh object pointed to by bulk_mesh_pt and add the elements to the
/// Mesh object pointeed to by surface_mesh_pt.
//=======================================================================
//template<class ELEMENT>
//void TiltedCavityProblem<ELEMENT>::
template<class ELEMENT>
void QuarterCircleDrivenCavityProblem<ELEMENT>::
create_impenetrable_lagrange_elements(const unsigned &b,
                                      Mesh* const &bulk_mesh_pt,
                                      Mesh* const &surface_mesh_pt)
{
 // How many bulk elements are adjacent to boundary b?
 unsigned n_element = bulk_mesh_pt->nboundary_element(b);

 // Loop over the bulk elements adjacent to boundary b?
 for(unsigned e=0;e<n_element;e++)
  {
   // Get pointer to the bulk element that is adjacent to boundary b
   ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
    bulk_mesh_pt->boundary_element_pt(b,e));

   //What is the index of the face of element e along boundary b
   int face_index = bulk_mesh_pt->face_index_at_boundary(b,e);

   // Build the corresponding impose_impenetrability_element
   ImposeImpenetrabilityElement<ELEMENT>* flux_element_pt = new
    ImposeImpenetrabilityElement<ELEMENT>(bulk_elem_pt,
                                          face_index);

   // Add the prescribed-flux element to the surface mesh
   surface_mesh_pt->add_element_pt(flux_element_pt);

   // Loop over the nodes
   unsigned nnod=flux_element_pt->nnode();
   for (unsigned j=0;j<nnod;j++)
    {
     Node* nod_pt = flux_element_pt->node_pt(j);

     // Is the node also on boundary 0 or 2?
     if ((nod_pt->is_on_boundary(0))||(nod_pt->is_on_boundary(2)))
      {
       // How many nodal values were used by the "bulk" element
       // that originally created this node?
       unsigned n_bulk_value=flux_element_pt->nbulk_value(j);

       // The remaining ones are Lagrange multipliers and we pin them.
       unsigned nval=nod_pt->nvalue();
       for (unsigned j=n_bulk_value;j<nval;j++)
        {
         nod_pt->pin(j);
        }
      }
    }
  }
}

int str2int(const std::string &str)
{
  std::stringstream ss(str);
  int n;
  ss >> n;
  return n;
}

unsigned str2unsigned(const std::string &str)
{
  std::stringstream ss(str);
  unsigned n;
  ss >> n;
  return n;
}

double str2double(const std::string &str)
{
  std::stringstream ss(str);
  double n;
  ss >> n;
  return n;
}

//===start_of_main======================================================
/// Driver code
//======================================================================
int main(int argc, char* argv[])
{
 #ifdef OOMPH_HAS_MPI
 // Initialise MPI
 MPI_Helpers::init(argc,argv);
 #endif
 
 /*
 SquareLagrangeVariables problemvars;
 problemvars.Dim = 2;
 std::cout << "probvar: " << problemvars.Dim << endl;

 // Store commandline arguments
 CommandLineArgs::setup(argc,argv);
 CommandLineArgs::specify_command_line_flag("--dim", &problemvars.Dim);
cout << "new dim: " << problemvars.Dim << endl;
CommandLineArgs::parse_and_assign();
CommandLineArgs::doc_specified_flags();
cout << "new dim: " << problemvars.Dim << endl;

pause("closer");

 */
//cout << "MathematicalConstants::Pi " << MathematicalConstants::Pi << endl;
//pause("damn");

 SquareLagrangeVariables myvar = SquareLagrangeVariables();

 // Set up doc info
 DocInfo doc_info;
 //doc_info.number()=0;
 doc_info.set_directory("RESLT");

 //Doc number of gmres iterations
 //ofstream out_file;

 // Store commandline arguments
 CommandLineArgs::setup(argc,argv);

 CommandLineArgs::specify_command_line_flag("--visc", &myvar.Vis);
 CommandLineArgs::specify_command_line_flag("--rey", &myvar.Rey);
 CommandLineArgs::specify_command_line_flag("--noel", &myvar.Noel);

 CommandLineArgs::parse_and_assign();
 CommandLineArgs::doc_specified_flags();
/*
 cout << "doc_soln: "  << doc_soln << ", doc_prec: "<< doc_prec << endl;

cout << "noel: " << myvar.Noel << " Sigma: " << myvar.Scaling_sigma << endl;
cout << "Visc: " << myvar.Vis_str << " Prec: " << myvar.Prec << endl;
*/

 // Set a string to identify the problem. This is unique to each problem,
 // so we hard code this. 2DStrPo = 2 dimension, straight parallel outflow.
 // straight describes the velocity flow field. Po = Parallel outflow
 // describes the boundary type.
 myvar.Prob_str = "Qcirc2D";


 // Set the viscuous term.
 if(CommandLineArgs::command_line_flag_has_been_set("--visc"))
 {
   if (myvar.Vis == 0)
   {
     myvar.Vis_str = "Sim";
     NavierStokesEquations<2>::Gamma[0]=0.0;
     NavierStokesEquations<2>::Gamma[1]=0.0;

   }
   else if (myvar.Vis == 1)
   {
     myvar.Vis_str = "Str";
     NavierStokesEquations<2>::Gamma[0]=1.0;
     NavierStokesEquations<2>::Gamma[1]=1.0;
   } // else - setting viscuous term.
   else
   {
     std::cout << "There is no such Viscous term, using 0 = simple." 
               << std::endl; 
   }
 }

 // Set Noel_str, used for book keeping.
 if(CommandLineArgs::command_line_flag_has_been_set("--noel"))
 {
   std::ostringstream strs;
   strs << "N" << myvar.Noel;
   myvar.Noel_str = strs.str();
 }


//*
 // Solve with Taylor-Hood element, set up problem
 // Set Rey_str, used for book keeping.
 if(CommandLineArgs::command_line_flag_has_been_set("--rey"))
 {
   if(myvar.Vis < 0)
   {
     myvar.Loop_reynolds = true;
   }
   else
   {
     std::ostringstream strs;
     strs << "R" << myvar.Rey;
     myvar.Rey_str = strs.str();
   }
 }


 myvar.Re_invFr=myvar.Rey;


 // Set up downwards-Gravity vector
 // Zero gravity
 myvar.Gravity[0] = 0.0;
 myvar.Gravity[1] = 0.0;


 QuarterCircleDrivenCavityProblem <RefineableQTaylorHoodElement<2> > 
   problem(&rayvars::zero_body_force, myvar,doc_info);


///////////////////////////////////////////////////////////////////////////////

 myvar.Loop_reynolds = false;
 if(myvar.Loop_reynolds)
 {
//   double Rey_start = 0.0;
//   double Rey_end = 300.0;
//   for (myvar.Rey = Rey_start; myvar.Rey < Rey_end; myvar.Rey += 100.0) 
//   {
//     std::ostringstream strs;
//     strs << "R" << myvar.Rey;
//     myvar.Rey_str = strs.str();
//
//     // Setup the label. Used for doc solution and preconditioner.
//     if(myvar.NS_solver == 0)
//     {
//       doc_info.label() = myvar.Prob_str + myvar.W_str + myvar.NS_str + myvar.Vis_str
//         + myvar.Rey_str + myvar.Noel_str
//         + myvar.W_approx_str + myvar.Sigma_str;
//     }
//     else if(myvar.NS_solver == 1)
//     {
//       doc_info.label() = myvar.Prob_str
//         + myvar.W_str + myvar.NS_str + myvar.F_str + myvar.P_str
//         + myvar.Vis_str + myvar.Rey_str
//         + myvar.Noel_str + myvar.W_approx_str + myvar.Sigma_str;
//     }
//     else
//     {
//       pause("There is no such NS preconditioner");
//     }
//
//     time_t rawtime;
//     time(&rawtime);
//
//     std::cout << "RAYDOING: "
//       << doc_info.label()
//       << " on " << ctime(&rawtime) << std::endl;
//
//
//     // Solve the problem
//     problem.newton_solve();
//
//     //Output solution
//     problem.doc_solution(doc_info);
//
//     // We now output the iteration and time.
//     Vector<Vector<std::pair<unsigned, double> > > iters_times
//       = doc_info.iterations_and_times();
//
//     // Below outputs the iteration counts and time.
//     // Output the number of iterations
//     // Since this is a steady state problem, there is only
//     // one "time step".
//     //*
//
//     // Loop over the time step:
//     unsigned ntimestep = iters_times.size();
//
//     for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
//     {
//       // New timestep:
//       std::cout << "RAYITS:\t" << intimestep << "\t";
//       // Loop through the Newtom Steps
//       unsigned nnewtonstep = iters_times[intimestep].size();
//       unsigned sum_of_newtonstep_iters = 0;
//       for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
//           innewtonstep++)
//       {
//         sum_of_newtonstep_iters += iters_times[intimestep][innewtonstep].first;
//         std::cout << iters_times[intimestep][innewtonstep].first << " ";
//       }
//       double average_its = ((double)sum_of_newtonstep_iters)
//         / ((double)nnewtonstep);
//
//       // Print to one decimal place if the average is not an exact
//       // integer. Ohterwise we print normally.
//       ((unsigned(average_its*10))%10)?
//         std::cout << "\t"<< std::fixed << std::setprecision(1)
//         << average_its << "(" << nnewtonstep << ")" << std::endl:
//         std::cout << "\t"<< average_its << "(" << nnewtonstep << ")" << std::endl;
//
//     }
//
//     // Now doing the times
//     for(unsigned intimestep = 0; intimestep < ntimestep; intimestep++)
//     {
//       // New timestep:
//       std::cout << "RAYTIME:\t" << intimestep << "\t";
//       // Loop through the Newtom Steps
//       unsigned nnewtonstep = iters_times[intimestep].size();
//       double sum_of_newtonstep_times = 0;
//       for(unsigned innewtonstep = 0; innewtonstep < nnewtonstep;
//           innewtonstep++)
//       {
//         sum_of_newtonstep_times += iters_times[intimestep][innewtonstep].second;
//         std::cout << iters_times[intimestep][innewtonstep].second << " ";
//       }
//       double average_time = ((double)sum_of_newtonstep_times)
//         / ((double)nnewtonstep);
//
//       // Print to one decimal place if the average is not an exact
//       // integer. Ohterwise we print normally.
//       std::cout << "\t"<< average_time << "(" << nnewtonstep << ")" << std::endl;
//     } // for timesteps
///*
// // Output linear solver time
// double total_time = 0;
// cout << "RAYTIME: " << setprecision(15);
// for (unsigned j = 0; j < iter; j++)
// {
//   total_time += Global_Parameters::Linear_solver_time[j];
//   cout << Global_Parameters::Linear_solver_time[j] << " ";
// }
// double average_time = total_time/iter;
//
// // Print to one decimal place if the average its is not an exact
// // integer. Otherwise we print normally.
// std::cout << "\t" << average_time << std::endl;
//
// time(&rawtime);
// std::cout << "RAYDONE: "
//           << Global_Parameters::Current_settings
//           << " on " << ctime (&rawtime) << endl;
//*/
//  
//   
//   }

 }
 else
 {

   std::ostringstream strs;
   strs << "R" << myvar.Rey;
   myvar.Rey_str = strs.str();

   // Solve the problem
   problem.newton_solve();

   problem.doc_solution(doc_info);

 }


#ifdef OOMPH_HAS_MPI
// finalize MPI
MPI_Helpers::finalize();
#endif
 return(EXIT_SUCCESS);
} // end_of_main
