// A adapter en fonction du raffinement de maillage souhaité
h = 0.02;

// Rayon de la bouilloire
rk = 0.1;
// Centre du disque
xk = 0.0; yk = 0.0;
Point(1) = {xk, yk, 0.0, h};
Point(2) = {xk-rk, y1, 0.0, h/2};
Circle(1) = {2,1,2};

// Résistances

// Rayon des résistances
rr = 0.015;
// Centre du disque resistance 1
xr1 = 0; yr1 = 0;
Point(3) = {xr1, yr1, 0.0, h/4};
Point(4) = {xr1-rr, yr1, 0.0, h/8};
Circle(2) = {4,3,4};

Line Loop(1) = {1}; // Le rond
Line Loop(2) = {2}; // Le trou

Plane Surface(1) = {1,2}; 
