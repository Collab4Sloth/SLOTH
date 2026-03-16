//+
radius=0.00465;
angle = Pi/4.;
height=0.01;

//+
Point(1) = {0, 0, 0, 4.e-4};
Point(2) = {radius, 0, 0, 4.e-4};
Point(3) = {radius*Cos(angle), radius*Sin(angle), 0, 1.e-4};

//+ Lines
Line(1) = {1, 2};  
Line(2) = {1, 3};
//+ Circle
Circle(11) = {2, 1, 3};
Curve Loop(1) = {2, -11, -1};
Surface(1) = {1};


surfaceVector[] = Extrude {0, 0, 5.e-3} {
Surface{1};
Layers{50};
Recombine;
};
Physical Surface("Bottom",1) = surfaceVector[0];
Physical Surface("Top",2) = {1};
Physical Surface("Front",3) = surfaceVector[2];
Physical Surface("Rear",4) = surfaceVector[4];
Physical Surface("External",5) = surfaceVector[3];
Physical Volume("Pellet",6) = surfaceVector[1];


Transfinite Line {2} = 50;
Transfinite Line {1} = 50;