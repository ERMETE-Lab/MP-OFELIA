lc = 2; // the characteristic length of the mesh

// Point definition
Point(1) = {0,   0, 0, lc};
Point(2) = {15,  0, 0, lc};
Point(3) = {75,  0, 0, lc};
Point(4) = {105, 0, 0, lc};
Point(5) = {135, 0, 0, lc};
Point(6) = {165, 0, 0, lc};
Point(7)  = {0,   15, 0, lc};
Point(8)  = {15,  15, 0, lc};
Point(9)  = {75,  15, 0, lc};
Point(10) = {105, 15, 0, lc};
Point(11) = {0,   75,  0, lc};
Point(12) = {15,  75,  0, lc};
Point(13) = {75,  75,  0, lc};
Point(14) = {105, 75,  0, lc};
Point(15) = {135, 75,  0, lc};
Point(16) = {0,   105,  0, lc};
Point(17) = {15,  105,  0, lc};
Point(18) = {75,  105,  0, lc};
Point(19) = {105, 105,  0, lc};
Point(20) = {135, 105,  0, lc};
Point(21) = {120, 105,  0, lc};
Point(22) = {120, 120,  0, lc};
Point(23) = {105, 120,  0, lc};
Point(24) = {0, 135,  0, lc};
Point(25) = {0, 105,  0, lc};
Point(26) = {105, 135,  0, lc};
Point(27) = {0, 165,  0, lc};
Point(28) = {165, 165,  0, lc};

// Line Definition
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {7, 8};
Line(7) = {9, 10};
Line(8) = {11, 12};
Line(9) = {13, 14};
Line(10) = {14, 15};
Line(11) = {16, 17};
Line(12) = {17, 18};
Line(13) = {18, 19};
Line(14) = {19, 21};
Line(15) = {21, 20};
Line(16) = {23, 22};
Line(17) = {24, 26};
Line(18) = {27, 28};
Line(19) = {7, 1};
Line(20) = {8, 2};
Line(21) = {9, 3};
Line(22) = {10, 4};
Line(23) = {15, 5};
Line(24) = {28, 6};
Line(25) = {14, 10};
Line(26) = {11, 7};
Line(27) = {11, 16};
Line(28) = {12, 17};
Line(29) = {13, 18};
Line(30) = {14, 19};
Line(31) = {20, 15};
Line(32) = {16, 24};
Line(33) = {26, 23};
Line(34) = {23, 19};
Line(35) = {22, 21};
Line(36) = {27, 24};

// Curve Loop and surface definition
Curve Loop(1) = {6, 20, -1, -19};
Plane Surface(1) = {1};
Curve Loop(2) = {7, 22, -3, -21};
Plane Surface(2) = {2};
Curve Loop(3) = {28, -11, -27, 8};
Plane Surface(3) = {3};
Curve Loop(4) = {13, -30, -9, 29};
Plane Surface(4) = {4};
Curve Loop(5) = {12, -29, 9, 25, -7, 21, -2, -20, -6, -26, 8, 28};
Plane Surface(5) = {5};
Curve Loop(6) = {10, 23, -4, -22, -25};
Plane Surface(6) = {6};
Curve Loop(7) = {31, -10, 30, 14, 15};
Plane Surface(7) = {7};
Curve Loop(8) = {32, 17, 33, 34, -13, -12, -11};
Plane Surface(8) = {8};
Curve Loop(9) = {34, 14, -35, -16};
Plane Surface(9) = {9};
Curve Loop(10) = {18, 24, -5, -23, -31, -15, -35, -16, -33, -17, -36};
Plane Surface(10) = {10};

// Mesh parameters
Mesh.Algorithm = 6; 
Mesh.Optimize  = 1;

// Physical groups
Physical Surface("reg-1", 10) = {5};
Physical Surface("reg-2", 20) = {3, 1, 2, 4};
Physical Surface("reg-3", 30) = {8, 6};
Physical Surface("CR",    35) = {7};
Physical Surface("reg-4", 40) = {9};
Physical Surface("reg-5", 50) = {10};

Physical Curve("void", 1) = {18, 24};
Physical Curve("symmetry", 2) = {36, 32, 27, 26, 19, 1, 2, 3, 4, 5};


