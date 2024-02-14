#include <set>
#include <gmsh.h>
#include <math.h>

int main(int argc, char **argv)
{
  gmsh::initialize();
  gmsh::model::add("t4");

  double lc = 1e-2;
  double Rad = 5e-1;
  double rad1 = 3e-1;
  double rad2 = 2.9e-1;

  // Two torus, both with guide radius (Rad). Second radius of first
  // torus equals (rad1), second radius of second torus equals (rad2)
  
  // main centers
  gmsh::model::geo::addPoint(0, 0, 0, lc, 1); // main center
  gmsh::model::geo::addPoint(0, 0, rad1, lc, 2); // upper center of big torus
  gmsh::model::geo::addPoint(0, 0, -rad1, lc, 3); // lower center of big torus
  gmsh::model::geo::addPoint(0, 0, rad2, lc, 4); // upper center of small torus
  gmsh::model::geo::addPoint(0, 0, -rad2, lc, 5); // lower center of small torus
  
  // central small circle of big torus
  gmsh::model::geo::addPoint(Rad-rad1, 0, 0, lc, 6);
  gmsh::model::geo::addPoint(0, Rad-rad1, 0, lc, 7);
  gmsh::model::geo::addPoint(-Rad+rad1, 0, 0, lc, 8);
  gmsh::model::geo::addPoint(0, -Rad+rad1, 0, lc, 9);
  // upper circle of big torus
  gmsh::model::geo::addPoint(Rad, 0, rad1, lc, 10);
  gmsh::model::geo::addPoint(0, Rad, rad1, lc, 11);
  gmsh::model::geo::addPoint(-Rad, 0, rad1, lc, 12);
  gmsh::model::geo::addPoint(0, -Rad, rad1, lc, 13);
  // central big circle of big torus
  gmsh::model::geo::addPoint(Rad+rad1, 0, 0, lc, 14);
  gmsh::model::geo::addPoint(0, Rad+rad1, 0, lc, 15);
  gmsh::model::geo::addPoint(-Rad-rad1, 0, 0, lc, 16);
  gmsh::model::geo::addPoint(0, -Rad-rad1, 0, lc, 17);
  // lower circle of big torus
  gmsh::model::geo::addPoint(Rad, 0, -rad1, lc, 18);
  gmsh::model::geo::addPoint(0, Rad, -rad1, lc, 19);
  gmsh::model::geo::addPoint(-Rad, 0, -rad1, lc, 20);
  gmsh::model::geo::addPoint(0, -Rad, -rad1, lc, 21);

  int radii[4] = {1, 2, 1, 3};
  for (int i = 0; i < 4; i++){
    int k = 4*i+5;
    gmsh::model::geo::addCircleArc(k + 1, radii[i], k + 2, 4*i+1);
    gmsh::model::geo::addCircleArc(k + 2, radii[i], k + 3, 4*i+2);
    gmsh::model::geo::addCircleArc(k + 3, radii[i], k + 4, 4*i+3);
    gmsh::model::geo::addCircleArc(k + 4, radii[i], k + 1, 4*i+4);
  }

  gmsh::model::geo::addPoint(Rad, 0, 0, lc, 22);
  gmsh::model::geo::addPoint(0, Rad, 0, lc, 23);
  gmsh::model::geo::addPoint(-Rad, 0, 0, lc, 24);
  gmsh::model::geo::addPoint(0, -Rad, 0, lc, 25);

  for (int i = 0; i < 4; i++){
    gmsh::model::geo::addCircleArc(6+i, 22+i, 10+i, 17+4*i);
    gmsh::model::geo::addCircleArc(10+i, 22+i, 14+i, 18+4*i);
    gmsh::model::geo::addCircleArc(14+i, 22+i, 18+i, 19+4*i);
    gmsh::model::geo::addCircleArc(18+i, 22+i, 6+i, 20+4*i);
  }

  for (int i = 0; i < 4; i++){
    gmsh::model::geo::addCurveLoop({1+4*i, 21+i, -((5+4*i)%16), -(17+i)}, 1+4*i);
    gmsh::model::geo::addSurfaceFilling({1+4*i}, 1+4*i);
    gmsh::model::geo::addCurveLoop({2+4*i, 25+i, -((6+4*i)%16), -(21+i)}, 2+4*i);
    gmsh::model::geo::addSurfaceFilling({2+4*i}, 2+4*i);
    gmsh::model::geo::addCurveLoop({3+4*i, 29+i, -((7+4*i)%16), -(25+i)}, 3+4*i);
    gmsh::model::geo::addSurfaceFilling({3+4*i}, 3+4*i);
    gmsh::model::geo::addCurveLoop({4+4*i, 17+i, -((8+4*i)%17 + (8+4*i)/17), -(29+i)}, 4+4*i);
    gmsh::model::geo::addSurfaceFilling({4+4*i}, 4+4*i);
  }

  // central small circle of small torus
  int p6 = gmsh::model::geo::addPoint(Rad-rad2, 0, 0, lc);
  int p7 = gmsh::model::geo::addPoint(0, Rad-rad2, 0, lc);
  int p8 = gmsh::model::geo::addPoint(-Rad+rad2, 0, 0, lc);
  int p9 = gmsh::model::geo::addPoint(0, -Rad+rad2, 0, lc);
  // upper circle of small torus
  int p10 = gmsh::model::geo::addPoint(Rad, 0, rad2, lc);
  int p11 = gmsh::model::geo::addPoint(0, Rad, rad2, lc);
  int p12 = gmsh::model::geo::addPoint(-Rad, 0, rad2, lc);
  int p13 = gmsh::model::geo::addPoint(0, -Rad, rad2, lc);
  // central big circle of small torus
  int p14 = gmsh::model::geo::addPoint(Rad+rad2, 0, 0, lc);
  int p15 = gmsh::model::geo::addPoint(0, Rad+rad2, 0, lc);
  int p16 = gmsh::model::geo::addPoint(-Rad-rad2, 0, 0, lc);
  int p17 = gmsh::model::geo::addPoint(0, -Rad-rad2, 0, lc);
  // lower circle of small torus
  int p18 = gmsh::model::geo::addPoint(Rad, 0, -rad2, lc);
  int p19 = gmsh::model::geo::addPoint(0, Rad, -rad2, lc);
  int p20 = gmsh::model::geo::addPoint(-Rad, 0, -rad2, lc);
  int p21 = gmsh::model::geo::addPoint(0, -Rad, -rad2, lc);
  
  int p22 = gmsh::model::geo::addPoint(Rad, 0, 0, lc);
  int p23 = gmsh::model::geo::addPoint(0, Rad, 0, lc);
  int p24 = gmsh::model::geo::addPoint(-Rad, 0, 0, lc);
  int p25 = gmsh::model::geo::addPoint(0, -Rad, 0, lc);
  int points[26] = {0, 1, 2, 3, 4, 5,
		    p6, p7, p8, p9,
		    p10, p11, p12, p13,
		    p14, p15, p16, p17,
		    p18, p19, p20, p21,
		    p22, p23, p24, p25};

  int centres[4] = {1, 4, 1, 5};
  int arcs[33];
  arcs[0] = 0;

  for (int i = 0; i < 4; i++){
    int k = 4*i+5;
    arcs[4*i+1] = gmsh::model::geo::addCircleArc(points[k + 1], centres[i], points[k + 2]);
    arcs[4*i+2] = gmsh::model::geo::addCircleArc(points[k + 2], centres[i], points[k + 3]);
    arcs[4*i+3] = gmsh::model::geo::addCircleArc(points[k + 3], centres[i], points[k + 4]);
    arcs[4*i+4] = gmsh::model::geo::addCircleArc(points[k + 4], centres[i], points[k + 1]);
  }

  for (int i = 0; i < 4; i++){
    arcs[17+4*i] = gmsh::model::geo::addCircleArc(points[6+i], points[22+i], points[10+i]);
    arcs[18+4*i] = gmsh::model::geo::addCircleArc(points[10+i], points[22+i], points[14+i]);
    arcs[19+4*i] = gmsh::model::geo::addCircleArc(points[14+i], points[22+i], points[18+i]);
    arcs[20+4*i] = gmsh::model::geo::addCircleArc(points[18+i], points[22+i], points[6+i]);
  }

  for (int i = 0; i < 4; i++){
    gmsh::model::geo::addCurveLoop({arcs[1+4*i], arcs[21+i], -arcs[((5+4*i)%16)],
	-arcs[(17+i)]}, 16+1+4*i);
    gmsh::model::geo::addSurfaceFilling({16+1+4*i}, 16+1+4*i);
    gmsh::model::geo::addCurveLoop({arcs[2+4*i], arcs[25+i], -arcs[((6+4*i)%16)],
	-arcs[(21+i)]}, 16+2+4*i);
    gmsh::model::geo::addSurfaceFilling({16+2+4*i}, 16+2+4*i);
    gmsh::model::geo::addCurveLoop({arcs[3+4*i], arcs[29+i], -arcs[((7+4*i)%16)],
	-arcs[(25+i)]}, 16+3+4*i);
    gmsh::model::geo::addSurfaceFilling({16+3+4*i}, 16+3+4*i);
    gmsh::model::geo::addCurveLoop({arcs[4+4*i], arcs[17+i], -arcs[((8+4*i)%17 + (8+4*i)/17)],
	-arcs[(29+i)]}, 16+4+4*i);
    gmsh::model::geo::addSurfaceFilling({16+4+4*i}, 16+4+4*i);
  }

  gmsh::model::geo::addSurfaceLoop({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
      17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32}, 1);
  gmsh::model::geo::addVolume({1});
  
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(3);

  gmsh::write("t4.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}
