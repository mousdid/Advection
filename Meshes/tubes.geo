// Rayon des disques
r = 0.125;

// A adapter en fonction du raffinement de maillage souhaité
//h = 0.05;
h = 0.2;

Point(1) = {-2.0, -1.0, 0.0, h};
Point(2) = {-2.0, 1.0, 0.0, h};
Point(3) = {2.0, -1.0, 0.0, h};
Point(4) = {2.0, 1.0, 0.0, h};

// Raffinement en h/2 autour du premier disque
// Centre du disque
x1 = -1.0; y1=0.0;
Point(5) = {x1, y1, 0.0, h};
Point(6) = {x1-r, y1, 0.0, h/2};
// Raffinement en h/2 autour du second disque
// Centre du disque
x2 = 0.0; y2=0.5;
Point(7) = {x2, y2, 0.0, h};
Point(8) = {x2-r, y2, 0.0, h/2};
// Raffinement en h/2 autour du troisième disque
// Centre du disque
x3 = +0.5; y3=-0.25;
Point(9) = {x3, y3, 0.0, h};
Point(10) = {x3-r, y3, 0.0, h/2};

// Cercle passant par 5, de centre 6, jusqu'à 5
// Ce sera la boundary 1
Circle(1) = {6,5,6};

// Ce seront les boundaries 2,3,4,5
Line(2) = {2, 1};
Line(3) = {1, 3};
Line(4) = {3, 4};
Line(5) = {4, 2};

Circle(6) = {8,7,8};
Circle(7) = {10,9,10};

Line Loop(9) = {2, 3, 4, 5};
Line Loop(10) = {1}; // Le 1° cercle
Line Loop(11) = {6}; // Le 2° cercle
Line Loop(12) = {7}; // Le 3° cercle

Plane Surface(1) = {9, 10, 11, 12}; // La plaque et les cercles
Physical Surface(1) = {1};
