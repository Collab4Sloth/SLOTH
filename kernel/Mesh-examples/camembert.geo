//+
radius=0.00607; 
//+ 0.00465;
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
//+
Physical Surface("Pellet", 4) = {1};

// Maillage
Transfinite Line {1} = 250;
Transfinite Line {2} = 250;


