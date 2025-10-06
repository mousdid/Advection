#ifndef _FINITEVOLUME_CPP

#include "FiniteVolume.h"
#include <fstream>
#include <iostream>


using namespace std;
using namespace Eigen;

// Constructeur
FiniteVolume::FiniteVolume(Function* function, DataFile* data_file, Mesh2D* mesh) :
_fct(function), _df(data_file), _msh(mesh) 
{
	this->_indicateur=0;
	std::cout << "Build finite volume class." << std::endl;
	std::cout << "-------------------------------------------------" << std::endl;
}

// Construit la matrice des flux
void FiniteVolume::Build_flux_mat_and_rhs(const double& t)
{
	// Matrix
	this->_mat_flux.resize(this->_msh->Get_triangles().size(),this->_msh->Get_triangles().size());
	// RHS
	this->_BC_RHS.resize(this->_msh->Get_triangles().size());
	this->_BC_RHS.setZero();////remise a 0 au debut
	vector<Triplet<double>> triplets;	triplets.clear();
	for (unsigned int i = 0; i < this->_msh->Get_edges().size(); i++)
	{
		string BC = _msh->Get_edges()[i].Get_BC(); //Renvoie "Dirichlet" ou "Neumann"
		int t1 = _msh->Get_edges()[i].Get_T1();
		int t2 = _msh->Get_edges()[i].Get_T2();

		double mu,alpha,beta,delta,volumemaille,e,v_n;
		e=_msh->Get_edges_length()[i];
		mu=this->_df->Get_mu();
		v_n=this->_fct->Velocity_x(this->_msh->Get_edges_center()(i,0),this->_msh->Get_edges_center()(i,1),t)*this->_msh->Get_edges_normal()(i,0)+this->_fct->Velocity_y(this->_msh->Get_edges_center()(i,0),this->_msh->Get_edges_center()(i,1),t)*this->_msh->Get_edges_normal()(i,1);
		

		if ((t2==-1))
		{
			//delta_ik 
		delta=2.*sqrt(pow(this->_msh->Get_triangles_center()(t1,0)-this->_msh->Get_edges_center()(i,0),2)+pow(this->_msh->Get_edges_center()(i,1)-this->_msh->Get_triangles_center()(t1,1),2));
		volumemaille=_msh->Get_triangles_area()[t1];
		if (BC==("Neumann"))
		{
			//diffusif

		this->_BC_RHS(t1)+=-e*mu*this->_fct->Neumann_Function(this->_msh->Get_edges_center()(i,0),this->_msh->Get_edges_center()(i,1),t)/volumemaille;

			//advectif


//calcul de alpha et beta selon le cas
		if (this->_df->Get_numerical_flux_choice() ==("upwind"))
		{
		if(v_n>=0.)
		{
		alpha= v_n;
		beta=0.;//beta nulle 
		}
		else
		{	
		beta= v_n;
		alpha=0.;//alpha nulle
		}
		}
		else if (this->_df->Get_numerical_flux_choice() =="centered")
		{
		alpha= v_n/2.0;
		beta= alpha;
		}
		

		this->_BC_RHS(t1)+=beta*e*delta*this->_fct->Neumann_Function(this->_msh->Get_edges_center()(i,0),this->_msh->Get_edges_center()(i,1),t)/volumemaille;
		triplets.push_back({t1,t1,(alpha+beta)*e/volumemaille});
		}



		else if (BC==("Dirichlet"))
		{
			///diffusif

		alpha= 2*mu/delta;
		beta= -2*mu/delta;
		this->_BC_RHS(t1)+=e*beta*this->_fct->Dirichlet_Function(this->_msh->Get_edges_center()(i,0),this->_msh->Get_edges_center()(i,1),t)/volumemaille;
		triplets.push_back({t1,t1,alpha*e/volumemaille});

				//advectif


//calcul de alpha et beta selon le cas
		if (this->_df->Get_numerical_flux_choice() ==("upwind"))
		{
		if(v_n>=0.)
		{
		alpha= v_n;
		beta=0.;//beta nulle 
		}
		else
		{	
		beta= v_n;
		alpha=0.;//alpha nulle
		}
		}
		else if (this->_df->Get_numerical_flux_choice() =="centered")
		{
		alpha= v_n/2.0;
		beta= alpha;
		}
		

		this->_BC_RHS(t1)+=2.*beta*e*this->_fct->Dirichlet_Function(this->_msh->Get_edges_center()(i,0),this->_msh->Get_edges_center()(i,1),t)/volumemaille;
		triplets.push_back({t1,t1,(alpha-beta)*e/volumemaille});
		}


		}

		else //general
		{
		// ///////pushback le terme diffusif Centered

		delta=sqrt(pow(this->_msh->Get_triangles_center()(t1,0)-this->_msh->Get_triangles_center()(t2,0),2)+pow(this->_msh->Get_triangles_center()(t2,1)-this->_msh->Get_triangles_center()(t1,1),2));
		alpha= mu/delta;
		beta= -mu/delta;
		volumemaille=_msh->Get_triangles_area()[t1];

		triplets.push_back({t1,t1,alpha*e/volumemaille});///pour t1 diagonale
		triplets.push_back({t1,t2,beta*e/volumemaille});///pour t1 hors diagonale

        volumemaille=_msh->Get_triangles_area()[t2];
		triplets.push_back({t2,t2,alpha*e/volumemaille});///pour t2 diagonale
		triplets.push_back({t2,t1,beta*e/volumemaille});//pour t2 hors diagonale
		
		if (this->_df->Get_numerical_flux_choice() ==("upwind"))
		{
		// 	//////pushback le terme advectif Upwind
///la normale a prendre pour t2 est - la normale pour t1 
		if(v_n>=0.)
		{
		//cout << v_n <<" positif" << endl;
		
		alpha= v_n;
		beta=0.;
		volumemaille=_msh->Get_triangles_area()[t1];
		triplets.push_back({t1,t1,alpha*e/volumemaille});///pour t1 diagonale
        volumemaille=_msh->Get_triangles_area()[t2];
		triplets.push_back({t2,t1,-alpha*e/volumemaille});///pour t2 diagonale
		}
		else
		{
			//cout << v_n << "negatif" << endl;	
		beta= v_n;
		alpha=0.;
		volumemaille=_msh->Get_triangles_area()[t2];
		triplets.push_back({t2,t2,-beta*e/volumemaille});///pour t1 hors diagonale
        volumemaille=_msh->Get_triangles_area()[t1];
		triplets.push_back({t1,t2,beta*e/volumemaille});//pour t2 hors diagonale
		
		

		}
		
		
		
		}
		else if (this->_df->Get_numerical_flux_choice() =="centered")
		{
			// 	//////pushback le terme advectif Centered
		alpha= v_n/2.0;
		beta= alpha;
		volumemaille=_msh->Get_triangles_area()[t1];
		triplets.push_back({t1,t1,alpha*e/volumemaille});///pour t1 diagonale
		triplets.push_back({t1,t2,beta*e/volumemaille});///pour t1 hors diagonale
        volumemaille=_msh->Get_triangles_area()[t2];
		triplets.push_back({t2,t2,alpha*e/volumemaille});///pour t2 diagonale
		triplets.push_back({t2,t1,beta*e/volumemaille});//pour t2 hors diagonale
		
			
		
		}
		else
		{
			std::cout << "vous avez mal entré le type du schéma pour le flux";
			std::cout << "-------------------------------------------------" << std::endl;
		}
		}
		}
	this->_mat_flux.setFromTriplets(triplets.begin(), triplets.end());
	
}


// --- Déjà implémenté ---
// Construit la condition initiale au centre des triangles
VectorXd FiniteVolume::Initial_condition()
{
	VectorXd sol0(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sol0(i) = this->_fct->Initial_condition(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1));
	}

	return sol0;
}

// Terme source au centre des triangles
VectorXd FiniteVolume::Source_term(double t)
{
	VectorXd sourceterm(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		sourceterm(i) = this->_fct->Source_term(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1), t);
	}

	return sourceterm;
}

// Solution exacte au centre des triangles
VectorXd FiniteVolume::Exact_solution(const double t)
{
	VectorXd exactsol(this->_msh->Get_triangles().size());

	for (unsigned int i = 0; i < this->_msh->Get_triangles().size(); i++)
	{
		exactsol(i) = this->_fct->Exact_solution(this->_msh->Get_triangles_center()(i,0),
		this->_msh->Get_triangles_center()(i,1), t);
	}
	return exactsol;
}

// Sauvegarde la solution
void FiniteVolume::Save_sol(const Eigen::VectorXd& sol, int n, std::string st)
{
	double norm = 0,moy=0,sum=0;
	
	
	for (unsigned int i = 0; i < sol.rows(); i++)
	{
		norm += sol(i)*sol(i)*this->_msh->Get_triangles_area()[i];
		moy += sol(i)*abs(this->_msh->Get_triangles_area()[i]);
		sum += abs(this->_msh->Get_triangles_area()[i]);
		
		

	}

	

	norm = sqrt(norm);
	moy/=sum;

	if (st == "solution")
	{
		cout << "Norme de u = " << norm << endl;
		cout << "moyenne temp = " << moy-273.15 << endl;
		cout << "maximum  temp = " << sol.maxCoeff()-273.15 << endl;
		

	}

	string name_file = this->_df->Get_results() + "/" + st + "_" + std::to_string(n) + ".vtk";
	unsigned int nb_vert = this->_msh->Get_vertices().size();
	assert(((long unsigned int)sol.size() == this->_msh->Get_triangles().size())
	&& "The size of the solution vector is not the same than the number of _triangles !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution << "# vtk DataFile Version 3.0 " << endl;
	solution << "2D Unstructured Grid" << endl;
	solution << "ASCII" << endl;
	solution << "DATASET UNSTRUCTURED_GRID" << endl;

	solution << "POINTS " << nb_vert << " float " << endl;
	for (unsigned int i = 0 ; i < nb_vert ; ++i)
	{
		solution << ((this->_msh->Get_vertices()[i]).Get_coor())[0] << " "
		<< ((this->_msh->Get_vertices()[i]).Get_coor())[1] << " 0." << endl;
	}
	solution << endl;

	solution << "CELLS " << this->_msh->Get_triangles().size() << " "
	<< this->_msh->Get_triangles().size()*4 << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 3 << " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[0]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[1]
		<< " " << ((this->_msh->Get_triangles()[i]).Get_vertices())[2] << endl;
	}
	solution << endl;

	solution << "CELL_TYPES " << this->_msh->Get_triangles().size() << endl;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << 5 << endl;
	}
	solution << endl;

	solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS sol float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	double eps = 1.0e-10;
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		solution << max(eps,sol[i]) << endl;
	}
	solution << endl;

Eigen::VectorXd vitessecarc;
vitessecarc=sol-sol;
double v_n=0.;

for (unsigned int i = 0; i < this->_msh->Get_edges().size(); i++)
	{
	int t1 = _msh->Get_edges()[i].Get_T1();
	int t2 = _msh->Get_edges()[i].Get_T2();
	v_n=this->_fct->Velocity_x(this->_msh->Get_edges_center()(i,0),this->_msh->Get_edges_center()(i,1),1.0)*this->_msh->Get_edges_normal()(i,0)+this->_fct->Velocity_y(this->_msh->Get_edges_center()(i,0),this->_msh->Get_edges_center()(i,1),1.0)*this->_msh->Get_edges_normal()(i,1);

if(t2==-1)
{vitessecarc(t1)+=v_n/3.0;}
else
{vitessecarc(t1)+=v_n/3.0;
vitessecarc(t2)+=v_n/3.0;}
	}
	//solution << "CELL_DATA " << this->_msh->Get_triangles().size() << endl;
	solution << "SCALARS CFL float 1" << endl;
	solution << "LOOKUP_TABLE default" << endl;
	// To avoid strange behaviour (which appear only with Apple)
	// with Paraview when we have very small data (e-35 for example)
	for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
	{
		//solution << max(eps,this->_df->Get_dt()*fabs(sol[i])/this->_msh->Get_triangles_length()(i)) << endl;
		solution << max(eps,this->_df->Get_dt()*fabs(vitessecarc[i])/this->_msh->Get_triangles_length()(i)) << endl;
	}
	solution << endl;

	if (this->_df->Get_mu() > 1e-10)
	{
		solution << "SCALARS Pe float 1" << endl;
		solution << "LOOKUP_TABLE default" << endl;
		// To avoid strange behaviour (which appear only with Apple)
		// with Paraview when we have very small data (e-35 for example)
		for (unsigned int i = 0 ; i < this->_msh->Get_triangles().size() ; ++i)
		{
			//solution << max(eps,this->_msh->Get_triangles_length()(i)*fabs(sol[i])/this->_df->Get_mu()) << endl;
			solution << max(eps,this->_msh->Get_triangles_length()(i)*fabs(vitessecarc[i])/this->_df->Get_mu()) << endl;
		}
		solution << endl;

	}
solution << "VECTORS vel float" << endl;
for (int i = 0 ; i < _msh->Get_triangles().size() ; ++i)
{
solution << _fct->Velocity_x(_msh->Get_triangles_center()(i,0),
_msh->Get_triangles_center()(i,1),n*_df->Get_dt())
<< " " << _fct->Velocity_y(_msh->Get_triangles_center()(i,0),
_msh->Get_triangles_center()(i,1),n*_df->Get_dt()) << " 0" << endl;
}
solution << endl;
	solution.close();
}

#define _FINITEVOLUME_CPP
#endif
