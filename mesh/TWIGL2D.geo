L1 = 24;
L2 = 56;
L3 = 80;


//+
Point(1)  = { 0,    0, 0, 1.0};
Point(2)  = { L1,   0, 0, 1.0};
Point(3)  = { L2,   0, 0, 1.0};
Point(4)  = { L3,   0, 0, 1.0};

Point(5)  = {  0,  L1, 0, 1.0};
Point(6)  = {  0,  L2, 0, 1.0};
Point(7)  = {  0,  L3, 0, 1.0};

Point(8)  = { L1,  L1, 0, 1.0};
Point(9)  = { L1,  L2, 0, 1.0};
Point(10) = { L2,  L2, 0, 1.0};
Point(11) = { L2,  L1, 0, 1.0};

Point(12) = { L3,  L3, 0, 1.0};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 12};
//+
Line(5) = {12, 7};
//+
Line(6) = {7, 6};
//+
Line(7) = {6, 5};
//+
Line(8) = {5, 1};
//+
Line(9) = {5, 8};
//+
Line(10) = {8, 11};
//+
Line(11) = {11, 10};
//+
Line(12) = {10, 9};
//+
Line(13) = {9, 6};
//+
Line(14) = {9, 8};
//+
Line(15) = {8, 2};
//+
Line(16) = {3, 11};
//+
Curve Loop(1) = {5, 6, -13, -12, -11, -16, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 15, -1, -8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {7, 9, -14, 13};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {10, -16, -2, -15};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {12, 14, 10, 11};
//+
Plane Surface(5) = {5};
//+
Physical Surface("domain1", 10) = {5};
//+
Physical Surface("domain2", 20) = {3, 4};
//+
Physical Surface("domain3", 30) = {1, 2};
//+
Physical Curve("ext_boundary", 1) = {5, 4};
//+
Physical Curve("sym_boundary", 2) = {6, 7, 8, 1, 2, 3};

Mesh.MeshSizeFactor = 1.5;
