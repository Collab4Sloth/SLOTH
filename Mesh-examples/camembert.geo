// Gmsh project created on Mon Dec 20 13:54:49 2021
//+
radius=0.00465;
angle = Pi/4.;
Point(1) = {0, 0, 0, 1.e-4};
Point(2) = {radius, 0, 0, 1.e-4};
Point(3) = {radius*Cos(angle), radius*Sin(angle), 0, 1.e-4};
//+
Circle(1) = {2, 1, 3};
//+
Line(2) = {1, 2};
//+
Line(3) = {1, 3};
//+
//+
Curve Loop(1) = {3, -1, -2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("lower", 1) = {2};
//+
Physical Curve("upper", 2) = {3};
//+
Physical Curve("external", 3) = {1};
//+
Physical Surface("Pellet", 4) = {1};

// Maillage
Transfinite Line {2} = 100;
Transfinite Line {3} = 100;
Transfinite Surface "*";


//+
Coherence;
