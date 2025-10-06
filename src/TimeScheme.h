#ifndef _TIME_SCHEME_H

#include "FiniteVolume.h"

class TimeScheme
{
protected:
   // Pointeur vers la classe FiniteVolume
   FiniteVolume* _fin_vol;
   // Pointeur vers la classe FiniteVolume
   DataFile* _df;
   // Vecteur initial et vecteur solution
   Eigen::VectorXd _sol;
   // Time
   double _t;

public:
   // Constructeur par défaut
   TimeScheme(DataFile* data_file, FiniteVolume* adv);
   // Destructeur par défaut - Si la classe ne contient pas de destructeur par défaut
   // alors le compilateur en génère un implicitement.
   virtual ~TimeScheme();
   // Enregistre la solution un fichier
   void Save_solution(int n) {_fin_vol->Save_sol(_sol, n, "solution");};
   // Enregistre la solution un fichier
   void Save_solution(Eigen::VectorXd s, int n, std::string st) {_fin_vol->Save_sol(s, n, st);};
   // Une étape du schéma en temps
   virtual void Advance() = 0;
   // Permet de récupérer _sol
   const Eigen::VectorXd & Get_sol() const {return _sol;};;
};

class EulerScheme : public TimeScheme
{
public:
   // Constructeur
   EulerScheme(DataFile* data_file, FiniteVolume* lap);
   // Une étape du schéma en temps
   void Advance();
};

class ImplicitEulerScheme : public TimeScheme
{
private:
   Eigen::SparseLU<Eigen::SparseMatrix<double> > _solver_direct;
public:
   // Constructeur
   ImplicitEulerScheme(DataFile* data_file, FiniteVolume* lap);
   // Une étape du schéma en temps
   void Advance();
};

#define _TIME_SCHEME_H
#endif
