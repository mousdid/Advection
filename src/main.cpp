#include <iostream>
#include <fstream>
#include <chrono>
#include "TimeScheme.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{

   if (argc < 2)
   {
      cout << "Please, enter the name of your data file." << endl;
      exit(0);
   }
   const string data_file_name = argv[1];

   // ----------------------- Fichier de données --------------------------------
   DataFile* data_file = new DataFile(data_file_name);
   // ---------------------------------------------------------------------------

   // ------------------Définition du nombre d'itérations------------------------
   int nb_iterations = int(ceil((data_file->Get_tfinal()-data_file->Get_t0())
   /data_file->Get_dt()));
   data_file->Adapt_dt(data_file->Get_tfinal() / nb_iterations);
   // ---------------------------------------------------------------------------

   // ---------------------------- Résolution  ----------------------------------
   Mesh2D* mesh = new Mesh2D(data_file->Get_BC_ref(),data_file->Get_BC_type());
   mesh->Read_mesh(data_file->Get_mesh_name());
   Function* fct = new Function(data_file);
   FiniteVolume* fin_vol = new FiniteVolume(fct, data_file, mesh);
   TimeScheme* time_scheme = NULL;

   if (data_file->Get_scheme() == "ExplicitEuler")
   time_scheme = new EulerScheme(data_file, fin_vol);
   else
   time_scheme = new ImplicitEulerScheme(data_file, fin_vol);

   cout << "-------------------------------------------------" << endl;
   cout << "Search u such that : " << endl;
   cout << "dt u + div (v u) - div(mu grad(u)) = f" << endl;
   cout << "-------------------------------------------------" << endl;

   // Démarrage du chrono
   auto start = chrono::high_resolution_clock::now();


   if ((data_file->Get_scenario() == "diffusion_temp_neumann") ||(data_file->Get_scenario() == "diffusion_advection_temp_neumann"))// Boucle temporelle classique
   {
      cout << "Save initial condition " << endl;
      time_scheme->Save_solution(0);
      cout << "Time Loop" << endl;
      for (int n = 1; n <= nb_iterations; n++) // Boucle en temps
      {
         cout << "Iteration " << n << endl;
         time_scheme->Advance();
         if (n==nb_iterations)
         time_scheme->Save_solution(n);
      }
   }
   else //  Si on connait la solution exacte
   {
      cout << "Save initial condition " << endl;
      double t(data_file->Get_t0());
      time_scheme->Save_solution(0); // Sauvegarde condition initiale
      time_scheme->Save_solution(fin_vol->Exact_solution(t),0,"exact_solution");

      ofstream error_file;
      const string error_file_name = data_file->Get_results()+"/errorLinfini.txt";
      error_file.open(error_file_name, ios::out);
      error_file << t << " " << 0 << endl;

      cout << "Time Loop" << endl;
      for (int n = 1; n == nb_iterations; n++) // Boucle en temps
      {
         time_scheme->Advance();
         t+=data_file->Get_dt();
         if (n==nb_iterations){
         time_scheme->Save_solution(n);
         VectorXd exact_sol = fin_vol->Exact_solution(t);
         time_scheme->Save_solution(exact_sol,n,"exact_solution");
         VectorXd approx_sol = time_scheme->Get_sol();
         double error = ((approx_sol-exact_sol).array().abs()).maxCoeff();
         cout << "Error Linfini at t " <<  t << ": " << error << endl;
         error_file << t << " " << error << endl;
         }
      }
      error_file.close();
   }
   // ---------------------------------------------------------------------------

   // Fin du chrono
   auto finish = chrono::high_resolution_clock::now();
   double t = chrono::duration_cast<chrono::seconds>(finish-start).count();
   // Affichage du résultat
   cout << "Cela a pris "<< t << " seconds" << endl;

   delete time_scheme;
   delete fin_vol;
   delete data_file;
   delete fct;

   return 0;
}
