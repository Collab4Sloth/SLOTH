//+
radius=0.00465;
angle = Pi/4.;
Point(1) = {0, 0, 0, 1.e-4};
Point(2) = {radius, 0, 0, 1.e-4};
Point(3) = {radius*Cos(angle), radius*Sin(angle), 0, 1.e-4};

//+ Lines
Line(1) = {1, 2};  
Line(2) = {1, 3};
//+ Circle
Circle(11) = {2, 1, 3};
Curve Loop(1) = {2, -11, -1};
Surface(1) = {1};


Physical Curve("lower", 1) = {1};
//+
Physical Curve("upper", 2) = {2};
//+
Physical Curve("external", 3) = {11};

// Maillage
Transfinite Line {1} = 100;
Transfinite Line {2} = 100;


//+
SetFactory("OpenCASCADE");
//+
Circle(13) = {0.002, 0.0005, 0, 0.0001, 0, 2*Pi};
//+
Circle(14) = {0.003, 0.0005, 0, 0.0001, 0, 2*Pi};
//+
Circle(16) = {0.0035, 0.001, 0, 0.0001, 0, 2*Pi};
//+
Circle(17) = {0.0035, 0.002, 0, 0.0001, 0, 2*Pi};
//+
Curve Loop(2) = {13};
//+
Surface(2) = {2};
//+
//+
Curve Loop(5) = {17};
//+
Surface(3) = {5};
//+
Curve Loop(7) = {16};
//+
Surface(4) = {7};
//+
Curve Loop(9) = {14};
//+
Surface(5) = {9};
//+
Physical Surface("UnPellet", 18) = {2, 3, 4, 5};
//+
Physical Surface("Pellet", 18) = {1};
