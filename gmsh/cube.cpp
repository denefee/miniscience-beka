#include <set>
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("t1");

  double lc = 1e-2;
  gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
  gmsh::model::geo::addPoint(.1, 0, 0, lc, 2);
  gmsh::model::geo::addPoint(0, .1, 0, lc, 3);
  gmsh::model::geo::addPoint(0, 0, .1, lc, 4);
  gmsh::model::geo::addPoint(.1, .1, 0, lc, 5);
  gmsh::model::geo::addPoint(.1, 0, .1, lc, 6);
  gmsh::model::geo::addPoint(0, .1, .1, lc, 7);
  gmsh::model::geo::addPoint(.1, .1, .1, lc, 8);

  gmsh::model::geo::addLine(1, 2, 1);
  gmsh::model::geo::addLine(1, 3, 2);
  gmsh::model::geo::addLine(1, 4, 3);
  gmsh::model::geo::addLine(3, 5, 4);
  gmsh::model::geo::addLine(3, 7, 5);
  gmsh::model::geo::addLine(4, 7, 6);
  gmsh::model::geo::addLine(6, 8, 7);
  gmsh::model::geo::addLine(2, 5, 8);
  gmsh::model::geo::addLine(7, 8, 9);
  gmsh::model::geo::addLine(5, 8, 10);
  gmsh::model::geo::addLine(4, 6, 11);
  gmsh::model::geo::addLine(2, 6, 12);

  // for(int i = 0; i < 3; i++)
  //  gmsh::model::geo::addLine(i + 1, 4, i + 4);

  gmsh::model::geo::addCurveLoop({2, 4, -8, -1}, 1);
  gmsh::model::geo::addPlaneSurface({1}, 1);

  gmsh::model::geo::addCurveLoop({2, 5, -6, -3}, 2);
  gmsh::model::geo::addPlaneSurface({2}, 2);

  gmsh::model::geo::addCurveLoop({4, 10, -9, -5}, 3);
  gmsh::model::geo::addPlaneSurface({3}, 3);

  gmsh::model::geo::addCurveLoop({6, 9, -7, -11}, 4);
  gmsh::model::geo::addPlaneSurface({4}, 4);

  gmsh::model::geo::addCurveLoop({3, 11, -12, -1}, 5);
  gmsh::model::geo::addPlaneSurface({5}, 5);

  gmsh::model::geo::addCurveLoop({8, 10, -7, -12}, 6);
  gmsh::model::geo::addPlaneSurface({6}, 6);

  gmsh::model::geo::addSurfaceLoop({1, 2, 3, 4, 5, 6}, 1);
  gmsh::model::geo::addVolume({1});

  gmsh::model::geo::synchronize();

  gmsh::model::mesh::generate(3);

  gmsh::write("t1.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();

  return 0;
}