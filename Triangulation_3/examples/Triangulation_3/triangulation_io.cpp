#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>

#include <CGAL/IO/Triangulation_mesh_istream.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
typedef Epick::Point_3               Point_3;
typedef CGAL::Triangulation_3<Epick> T3;



int main(int argc, char* argv[])
{
  char* filename = (argc < 2) ? "data/anchor.mesh" : argv[1];

  std::ifstream ifs(filename);
  if (!ifs) return EXIT_FAILURE;

  T3 tr;
  CGAL::import_triangulation_3_from_mesh(ifs, tr);

  CGAL_assertion( tr.is_valid() );

  return EXIT_SUCCESS;
}
