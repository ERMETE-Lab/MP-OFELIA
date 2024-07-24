lc = 1;

scale = 1; // the geometry is made in cm

L_active = 366 / 2 / scale;
L_plug   = 1.27 / scale;

fuel_or = 0.819 / 2 / scale;
clad_ir = 0.836 / 2 / scale;
clad_or = 0.95  / 2 / scale;

R = 1.26 / 2 / scale; // pitch

// geometry: points
Point(1)  = {-L_active - L_plug, fuel_or, 0, lc};
Point(2)  = {-L_active,          fuel_or, 0, lc};
Point(3)  = {+L_active + L_plug, fuel_or, 0, lc};
Point(4)  = {+L_active,          fuel_or, 0, lc};

Point(5)  = {-L_active - L_plug, clad_ir, 0, lc};
Point(6)  = {-L_active,          clad_ir, 0, lc};
Point(7)  = {+L_active + L_plug, clad_ir, 0, lc};
Point(8)  = {+L_active,          clad_ir, 0, lc};

Point(9)  = {-L_active - L_plug, clad_or, 0, lc};
Point(10) = {-L_active,          clad_or, 0, lc};
Point(11) = {+L_active + L_plug, clad_or, 0, lc};
Point(12) = {+L_active,          clad_or, 0, lc};

Point(13) = {-L_active - L_plug, -fuel_or, 0, lc};
Point(14) = {-L_active,          -fuel_or, 0, lc};
Point(15) = {+L_active + L_plug, -fuel_or, 0, lc};
Point(16) = {+L_active,          -fuel_or, 0, lc};

Point(17) = {-L_active - L_plug, -clad_ir, 0, lc};
Point(18) = {-L_active,          -clad_ir, 0, lc};
Point(19) = {+L_active + L_plug, -clad_ir, 0, lc};
Point(20) = {+L_active,          -clad_ir, 0, lc};

Point(21) = {-L_active - L_plug, -clad_or, 0, lc};
Point(22) = {-L_active,          -clad_or, 0, lc};
Point(23) = {+L_active + L_plug, -clad_or, 0, lc};
Point(24) = {+L_active,          -clad_or, 0, lc};//+
Line(1) = {21, 22};
//+
Line(2) = {22, 24};
//+
Line(3) = {24, 23};
//+
Line(4) = {17, 18};
//+
Line(5) = {18, 20};
//+
Line(6) = {20, 19};
//+
Line(7) = {13, 14};
//+
Line(8) = {14, 16};
//+
Line(9) = {16, 15};
//+
Line(10) = {1, 2};
//+
Line(11) = {2, 4};
//+
Line(12) = {4, 3};
//+
Line(13) = {5, 6};
//+
Line(14) = {6, 8};
//+
Line(15) = {8, 7};
//+
Line(16) = {9, 10};
//+
Line(17) = {10, 12};
//+
Line(18) = {12, 11};
//+
Line(19) = {21, 17};
//+
Line(20) = {17, 13};
//+
Line(21) = {13, 1};
//+
Line(22) = {1, 5};
//+
Line(23) = {5, 9};
//+
Line(24) = {22, 18};
//+
Line(25) = {18, 14};
//+
Line(26) = {14, 2};
//+
Line(27) = {2, 6};
//+
Line(28) = {6, 10};
//+
Line(29) = {24, 20};
//+
Line(30) = {20, 16};
//+
Line(31) = {16, 4};
//+
Line(32) = {4, 8};
//+
Line(33) = {8, 12};
//+
Line(34) = {23, 19};
//+
Line(35) = {19, 15};
//+
Line(36) = {15, 3};
//+
Line(37) = {3, 7};
//+
Line(38) = {7, 11};
//+
Curve Loop(1) = {1, 2, 3, 34, 35, 36, 37, 38, -18, -17, -16, -23, -22, -21, -20, -19};
//+
Curve Loop(2) = {14, -32, -31, -30, -5, 25, 26, 27};
//+
Plane Surface(1) = {1, 2};
//+
Curve Loop(3) = {11, 32, -14, -27};
//+
Plane Surface(2) = {3};
//+
Curve Loop(4) = {8, -30, -5, 25};
//+
Plane Surface(3) = {4};
//+
Curve Loop(5) = {11, -31, -8, 26};
//+
Plane Surface(4) = {5};
//+
Physical Surface("fuel", 1) = {4};
//+
Physical Surface("gap", 2) = {2, 3};
//+
Physical Surface("cladding", 3) = {1};
//+
Physical Curve("water_contact", 10) = {2, 1, 3, 18, 17, 16};
//+
Physical Curve("inlet", 20) = {23, 22, 21, 20, 19};
//+
Physical Curve("outlet", 30) = {34, 35, 36, 37, 38};

Mesh.MeshSizeFactor = 0.05;

Transfinite Surface {2,3,4};
// Recombine Surface {1,2,3,4};
//+
Transfinite Curve {2, 5, 8, 11, 14, 17} = 1000 Using Progression 1;