#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"
#include <iostream>

using namespace Eigen;
using namespace std;

// Constructeur par défaut (ne pas oublier de mettre votre pointeur à 0 !!)
TimeScheme::TimeScheme(DataFile* data_file, FiniteVolume* adv) :
_fin_vol(adv),_df(data_file), _sol(adv->Initial_condition()), _t(_df->Get_t0())
{
}

EulerScheme::EulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
}

ImplicitEulerScheme::ImplicitEulerScheme(DataFile* data_file, FiniteVolume* adv) :
TimeScheme(data_file, adv)
{
   std::cout << "Build time scheme class." << std::endl;
   std::cout << "-------------------------------------------------" << std::endl;
}

// Destructeur (car on a des fonctions virtuelles)
TimeScheme::~TimeScheme()
{
}

// Euler Explicite
void EulerScheme::Advance()
{
   // TODO
this->_fin_vol->Build_flux_mat_and_rhs(this->_t);
this->_sol+=-this->_df->Get_dt()*(this->_fin_vol->Get_flux_matrix()*this->_sol+this->_fin_vol->Get_BC_RHS()- this->_fin_vol->Source_term(this->_t));
this->_t+=this->_df->Get_dt() ;








}

// Euler Implicite
void ImplicitEulerScheme::Advance()
{
   // TODO
Eigen::SparseMatrix<double> A;
Eigen::VectorXd b;

this->_t+=this->_df->Get_dt() ;
this->_fin_vol->Build_flux_mat_and_rhs(this->_t);
A=this->_fin_vol->Get_flux_matrix();
// cout << "matrice A" <<endl << A << endl;
// cout << "vecteur b" <<endl << this->_fin_vol->Get_BC_RHS() << endl;
// cout << "vecteur S" <<endl << this->_fin_vol->Source_term(this->_t) << endl;

A.setIdentity();
SparseLU <SparseMatrix<double> > solver;
b=this->_sol+this->_df->Get_dt()*(this->_fin_vol->Source_term(this->_t)- this->_fin_vol->Get_BC_RHS());
solver.compute(A+this->_df->Get_dt()*this->_fin_vol->Get_flux_matrix());
this->_sol= solver.solve(b);




}

#define _TIME_SCHEME_CPP
#endif
