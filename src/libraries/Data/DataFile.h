#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>
// Définition de la classe

class DataFile {
private:
   const std::string _file_name;
   double _t0, _tfinal, _dt;

   double _mu;

   std::string _mesh_name;
   std::string _scheme;
   std::string _numerical_flux_choice;
   std::string _results;

   bool _print_info;
   bool _if_mesh_name;
   bool _if_t0;
   bool _if_tfinal;
   bool _if_dt;
   bool _if_scheme;
   bool _if_numerical_flux_choice;
   bool _if_results;
   bool _if_mu;

   int _N_BC;
   std::vector<int> _BC_ref;
   std::vector<std::string> _BC_type;
   std::string _which_scenario;

   double _P_kettle;
   double _U;

   public: // Méthodes et opérateurs de la classe
   DataFile(std::string file_name);

   void Adapt_dt(double dt){_dt = dt;};

   const bool Print_info() const {return _print_info;};

   const double Get_t0() const {return _t0;};
   const double Get_tfinal() const {return _tfinal;};
   const double Get_dt() const {return _dt;};
   const std::string Get_mesh_name() const {return _mesh_name;};
   const std::string Get_scheme() const {return _scheme;};
   const std::string Get_numerical_flux_choice() const {return _numerical_flux_choice;};

   const double Get_mu() const {return _mu;};

   const std::string Get_results() const {return _results;};
   const std::vector<int> Get_BC_ref() const {return _BC_ref;};
   const std::vector<std::string> Get_BC_type() const {return _BC_type;};

   const std::string Get_scenario() const {return _which_scenario;};

   const double Get_P_kettle() const {return _P_kettle;};
   const double Get_U() const {return _U;};


};

#define _DATA_FILE_H
#endif
