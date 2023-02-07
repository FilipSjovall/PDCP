// Gmsh project created on Thu Feb  2 16:47:28 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {-0.5, 0.5, 0, 1.0};
//+
Point(2) = {-0.5, -0.5, 0, 1.0};
//+
Point(3) = {0.5, -0.5, 0, 1.0};
//+
Point(4) = {0.5, 0.5, 0, 1.0};
//+
Line(1) = {1, 4};
//+
Line(2) = {3, 4};
//+
Line(3) = {2, 3};
//+
Line(4) = {1, 2};
//+
Circle(5) = {0, 0, 0, 0.25, 0, 2*Pi};
//+
Point(6) = {0, 0.25, 0, 1.0};
//+
Point(7) = {0, 0.5, 0, 1.0};
//+
Point(8) = {0.25, 0, 0, 1.0};
//+
Point(9) = {0.5, 0, 0, 1.0};
//+
Line(6) = {6, 7};
//+
Line(7) = {5, 9};
//+
Recursive Delete {
  Curve{4}; 
}
//+
Recursive Delete {
  Curve{1}; Curve{3}; 
}
//+
Recursive Delete {
  Curve{2}; 
}
//+
Point(10) = {0.5, 0.5, 0, 1.0};
//+
Line(8) = {7, 10};
//+
Line(9) = {9, 10};
//+
Recursive Delete {
  Curve{5}; 
}
//+
Point(11) = {0, 0, 0, 1.0};
//+
Circle(10) = {0.2, 0.1, 0, 0.25, 0, 2*Pi};
//+
Recursive Delete {
  Curve{10}; 
}
//+
Circle(10) = {6, 11, 5};
//+
Curve Loop(1) = {6, 8, -9, -7, -10};
//+
Plane Surface(1) = {1};
//+
Physical Surface(11) = {1};
//+
Physical Surface(12) = {1};
//+
Physical Curve(13) -= {6};
//+
Physical Curve(13) -= {8};
//+
Physical Curve(13) -= {8};
//+
Physical Curve(13) -= {8};
//+
Physical Curve(13) -= {9};
//+
Physical Curve(13) -= {7};
//+
Physical Curve(13) -= {10};
//+
Physical Curve(13) -= {10};
//+
Physical Surface(11) -= {1};
//+
Physical Surface(12) -= {1};
//+
Physical Curve(13) -= {6};
//+
Physical Curve(13) -= {8, 9, 7, 10, 6};
//+
Physical Curve(13) -= {6, 8, 10, 9, 7};
//+
Physical Surface(13) = {1};
//+
Physical Curve(14) -= {6, 8, 9, 7, 10};
//+
Physical Surface(13) -= {1};
//+
Physical Surface(14) -= {1};
//+
Physical Surface(14) -= {1};
//+
Physical Surface(14) -= {1};
//+
Physical Surface(14) -= {1};
//+
Physical Surface(14) -= {1};
//+
Physical Surface(14) -= {1};
//+
Physical Surface(14) -= {1};
//+
Physical Surface(14) -= {1};
//+
Physical Curve(14) -= {6, 8, 9, 10, 7};
//+
Physical Surface("surf", 14) = {1};
//+
Physical Surface(" surf", 14) -= {1};
//+
Physical Surface(15) -= {1};
//+
Physical Surface(15) -= {1};
//+
Physical Curve(14) -= {6, 8, 9, 7, 10};
//+
Physical Surface("surf", 14) = {1};
