L=0.002;
N=L/50;
//////////
// Points 
//////////
Point(1) = {0, 0, 0, N};
Point(2) = {L, 0, 0, N};
Point(3) = {L, L, 0, N};
Point(4) = {0, L, 0, N};
//////////
// Lines 
//////////
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//////////////////////////////
// periodicity
/////////////////////////////
Periodic Line {1} = {-3};
Periodic Line {4} = {-2};

////////////////////
// External contour
////////////////////
Curve Loop(1) = {1, 2, 3, 4};

//////////////////////////
// Surface
//////////////////////////
Plane Surface(1) = {1};
Transfinite Surface {1};

//////////////////////////////
// Physical groups for BCs
//////////////////////////////
Physical Curve("lower", 1) = {1};
Physical Curve("right", 2) = {2};
Physical Curve("upper",3) = {3};
Physical Curve("left",4) = {4};

Physical Surface(1) = {1};


