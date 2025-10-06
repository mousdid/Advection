#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include <toml/toml.hpp>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

DataFile::DataFile(std::string file_name)
: _file_name(file_name)
{
   // Lecture du fichier de données
   auto config = toml::parse(file_name);

   // Other
   const auto& other = toml::find(config, "other");
   this->_mesh_name = toml::find<std::string>(other, "mesh");
   this->_mu = toml::find<double>(other, "mu");
   this->_numerical_flux_choice = toml::find<std::string>(other, "numerical_flux");
   this->_results = toml::find<std::string>(other, "results");

   // Time
   const auto& time = toml::find(config, "time");
   this->_t0 = toml::find<double>(time, "t0");
   this->_tfinal = toml::find<double>(time, "tfinal");
   this->_dt = toml::find<double>(time, "dt");
   this->_scheme = toml::find<std::string>(time, "scheme");

   // Boundary conditions
   const auto& BC = toml::find(config, "BC");
   this->_BC_ref = toml::find<std::vector<int> >(BC, "ref");
   this->_BC_type = toml::find<std::vector<std::string> >(BC, "BC");

   // Scenarii
   const auto& scenarii = toml::find(config, "scenarii");
   this->_which_scenario = toml::find<std::string>(scenarii, "which_scenario");

   if(config.contains("physics") && config.at("physics").count("P_kettle") != 0)
   {
      const auto& physics = toml::find(config, "physics");
      this->_P_kettle = toml::find<double>(physics, "P_kettle");
   }
    if(config.contains("physics") && config.at("physics").count("U") != 0)
   {
      const auto& physics = toml::find(config, "physics");
      this->_U = toml::find<double>(physics, "U");
   }


   if ((this->_numerical_flux_choice != "centered") && (this->_numerical_flux_choice != "upwind"))
   {
      cout << "Only centered and upwind numerical flows are implemented." << endl;
      exit(0);
   }

   if ((this->_scheme != "ExplicitEuler") && (this->_scheme != "ImplicitEuler"))
   {
      cout << "Only Explicit Euler and Implicit Euler schemes are implemented." << endl;
      exit(0);
   }
   if (this->_scheme == "ExplicitEuler")
   {
      cout << "Beware to the CFL condition." << endl;
   }

   if ((this->_which_scenario == "diffusion_hom_neumann") || (this->_which_scenario == "diffusion_all_BC")
   || (this->_which_scenario == "advection_hom_neumann") ||  (this->_which_scenario == "advection_all_BC")
   || (this->_which_scenario == "diffusion_advection_all_BC") )
   {
      cout << "-------------------------------------------------" << endl;
      cout << "The test case: " << this->_which_scenario << " has been chosen." << endl;
      cout << "-------------------------------------------------" << endl;
      if (this->_mesh_name == "Meshes/square_mini.mesh")
      {
         this->_print_info = true;
      }
   }
    else if (this->_which_scenario == "diffusion_temp_neumann")
   {
      cout << "-------------------------------------------------" << endl;
      cout << "The real case: " << this->_which_scenario << " has been chosen." << endl;
      cout << "-------------------------------------------------" << endl;
   }
   else if (this->_which_scenario == "diffusion_advection_temp_neumann")
   {
      cout << "-------------------------------------------------" << endl;
      cout << "The real case: " << this->_which_scenario << " has been chosen." << endl;
      cout << "-------------------------------------------------" << endl;
   }
   else if (this->_which_scenario == "none")
   {
      cout << "-------------------------------------------------" << endl;
      cout << "It is not a test case!" << endl;
      cout << "-------------------------------------------------" << endl;
   }
   else
   {
      cout << "-------------------------------------------------" << endl;
      cout << "A scenario has to be chosen. If it is not a test case, consider none." << endl;
      cout << "-------------------------------------------------" << endl;
      exit(0);
   }


   if ( (this->_which_scenario == "advection_hom_neumann") && (fabs(this->_mu) > 1e-6) )
   {
      cout << "Only advection: mu has been fixed at 0." << endl;
      this->_mu = 0;
   }

   if (this->_which_scenario == "advection_all_BC")
   {
      this->_which_scenario = "diffusion_advection_all_BC";
      if (fabs(this->_mu) > 1e-6)
      {
         cout << "Only advection: mu has been fixed at 0." << endl;
         this->_mu = 0;
      }
   }

   // Créer le dossier de résultats
   system(("mkdir -p ./" +this->_results).c_str());
   // Supprimer les anciens résultats
   system(("rm -f ./" +this->_results + "/*.vtk").c_str());
   // Copier le fichier de données dans le dossier résultats
   system(("cp -r ./" + this->_file_name + " ./"
   + this->_results + "/params.txt").c_str());
}

#define _DATA_FILE_CPP
#endif
