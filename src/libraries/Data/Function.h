#ifndef _FUNCTION_H

#include "DataFile.h"

class Function {
private:
   // Pointeur de la classe DataFile pour récupérer toutes les
   // valeurs de paramètres
   const DataFile* _df;
   // Diffusion coefficient
   const double _mu;

   public: // Méthodes et opérateurs de la classe
   Function(DataFile* data_file);
   double Exact_solution(const double x, const double y, const double t) const;
   double Initial_condition(const double x, const double y) const;
   double Source_term(const double x, const double y, const double t) const;
   double Velocity_x(const double x, const double y, const double t) const;
   double Velocity_y(const double x, const double y, const double t) const;
   double Neumann_Function(const double x, const double y, const double t) const;
   double Dirichlet_Function(const double x, const double y, const double t) const;
};

#define _FUNCTION_H
#endif
