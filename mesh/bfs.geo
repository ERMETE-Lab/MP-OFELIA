// H = 5.2/1000;
// h = 4.9/1000;
// L = 140/1000;
// l = 40/1000;

zeta_comp_step = 1;
H = 2 * zeta_comp_step;
h = zeta_comp_step;
L = 36 * zeta_comp_step;
l = 25 * zeta_comp_step;

Point(1) = {0, 0, 0, 1.0};
Point(2) = {-l, 0, 0, 1.0};
Point(3) = {-l, H, 0, 1.0};
Point(4) = {0, H, 0, 1.0};

Point(5) = {L, H, 0, 1.0};
Point(6) = {L, 0, 0, 1.0};
Point(7) = {L, -h, 0, 1.0};
Point(8) = {0, -h, 0, 1.0};

Line(1) = {2, 1};
Line(2) = {1, 8};
Line(3) = {8, 7};
Line(4) = {7, 6};
Line(5) = {6, 1};
Line(6) = {6, 5};
Line(7) = {5, 4};
Line(8) = {4, 1};
Line(9) = {4, 3};
Line(10) = {3, 2};

// Curve Loop(1) = {7, 9, 10, 1, -5, 6};
// Curve Loop(2) = {7, 9, 10, 1, 2, 3, 4, 6};
// Plane Surface(1) = {2};

Curve Loop(1) = {9, 10, 1, -8};
Curve Loop(2) = {7, 8, -5, 6};
Curve Loop(3) = {5, 2, 3, 4};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

Physical Surface("domain", 100) = {1, 2, 3};

Physical Curve("inlet",  10) = {10};
Physical Curve("walls",  20) = {9, 7, 3, 1, 2};
Physical Curve("outlet", 30) = {6, 4};

Mesh.MeshSizeFactor = 0.05;