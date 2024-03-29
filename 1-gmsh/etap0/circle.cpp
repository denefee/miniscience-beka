// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 1
//
//  Geometry basics, elementary entities, physical groups
//
// -----------------------------------------------------------------------------

#include <set>

// The Gmsh C++ API is entirely defined in the `gmsh.h' header (which contains
// the full documentation of all the functions in the API):
#include <gmsh.h>

int main(int argc, char **argv)
{
  gmsh::initialize();
  gmsh::model::add("t2");

  double lc = 1e-2;
  gmsh::model::geo::addPoint(0, 0, 0, lc, 1);
  gmsh::model::geo::addPoint(.1, 0, 0, lc, 2);
  gmsh::model::geo::addPoint(0, .1, 0, lc, 3);
  gmsh::model::geo::addPoint(-0.1, 0, 0, lc, 4);
  gmsh::model::geo::addPoint(0, -0.1, 0, lc, 5);

  gmsh::model::geo::addCircleArc(2, 1, 3, 1);
  gmsh::model::geo::addCircleArc(3, 1, 4, 2);
  gmsh::model::geo::addCircleArc(4, 1, 5, 3);
  gmsh::model::geo::addCircleArc(5, 1, 2, 4);

  gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);

  gmsh::model::geo::addPlaneSurface({1}, 1);

  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate(2);
  gmsh::write("t2.msh");

  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();
  gmsh::finalize();

  return 0;
}