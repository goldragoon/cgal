
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_segment_traverser_3.h>

#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/foreach.hpp>

#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/intersections.h>

//#define CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE

// Define the kernel.
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3< Kernel >                DT;
typedef DT::Cell_handle                                         Cell_handle;

typedef CGAL::Triangulation_segment_cell_iterator_3< DT >       Cell_traverser;

typedef CGAL::Exact_predicates_exact_constructions_kernel       Epeck;


Epeck::Segment_3 intersection(const Epeck::Segment_3& seg,
                              const Epeck::Tetrahedron_3& tet)
{
  bool first_point_found = false;
  bool second_point_found = false;
  Epeck::Point_3 p1, p2;

  for (unsigned int i = 0; i < 4; ++i)
  {
    Epeck::Triangle_3 tri(tet[i], tet[(i + 1) % 4], tet[(i + 2) % 4]);
    if (do_intersect(tri, seg))
    {
      CGAL::cpp11::result_of<Epeck::Intersect_3(Epeck::Segment_3, Epeck::Triangle_3)>::type
        result = intersection(tri, seg);
      if (result) {
        if (const Epeck::Segment_3* s = boost::get<Epeck::Segment_3>(&*result))
          return *s;
        else if (const Epeck::Point_3* p = boost::get<Epeck::Point_3>(&*result))
        {
          if (!first_point_found)
          {
            p1 = *p;
            first_point_found = true;
          }
          else if (!second_point_found)
          {
            p2 = *p;
            second_point_found = true;
          }
        }
      }
    }
  }
  //no segment found
  assert(first_point_found && second_point_found);
  return Epeck::Segment_3(p1, p2);
}

template <typename Big_tuple>
bool test(const DT dt, const Big_tuple& tuple)
{
  bool result = true;
  using std::get;
  const auto& p1 = get<0>(tuple);
  const auto& p2 = get<1>(tuple);
//  const auto& expected_results = get<2>(tuple);

  std::cout << "\n#\n# Query segment: ( " << p1 << " , "
    << p2 << " )\n#\n";

  Cell_traverser ct(dt, p1, p2);
  unsigned int nb_cells = 0, nb_facets = 0, nb_edges = 0, nb_vertex = 0;

  CGAL::Cartesian_converter<Kernel, Epeck> to_exact;
  const Epeck::FT exact_length
    = CGAL::approximate_sqrt(CGAL::squared_distance(to_exact(p1), to_exact(p2)));
  Epeck::Segment_3 exact_query(to_exact(p1), to_exact(p2));

  // Count the number of finite cells traversed,
  // and sum exact lengths
  Epeck::FT exact_length_sum = 0;
  unsigned int inf = 0, fin = 0;

  for( ; ct != ct.end(); ++ct )
  {
    std::cerr << "Cell ( ";
    for(int i = 0; i < 4; ++i)
      std::cerr << ct->vertex(i)->point() << "  ";
    std::cerr << " )\n";

    if( dt.is_infinite(ct) )
        ++inf;
    else
    {
      ++fin;

      if (CGAL::do_intersect(exact_query, to_exact(dt.tetrahedron(ct))))
      {
        Epeck::Segment_3 s = intersection(exact_query, to_exact(dt.tetrahedron(ct)));
        exact_length_sum += CGAL::approximate_sqrt(CGAL::squared_distance(
                              s.source(), s.target()));
      }
      else
        assert(false);

      DT::Locate_type lt;
      int li, lj;
      ct.entry(lt, li, lj);

      switch (lt)
      {
      case DT::Locate_type::FACET:
       ++nb_facets;
       break;
      case DT::Locate_type::EDGE:
       ++nb_edges;
       break;
      case DT::Locate_type::VERTEX:
       ++nb_vertex;
       break;
      default:
       /*when source is in a cell*/
        ++nb_cells;
       CGAL_assertion(lt == DT::Locate_type::CELL);
      }
    }
  }

    std::cout << "Triangulation has " << dt.number_of_vertices() << " vertices."
      << std::endl;
    std::cout << "\tnb cells     : " << nb_cells << std::endl;
    std::cout << "\tnb facets    : " << nb_facets << std::endl;
    std::cout << "\tnb edges     : " << nb_edges << std::endl;
    std::cout << "\tnb vertices  : " << nb_vertex << std::endl;

    std::cout << "\tExact length      = " << exact_length << std::endl;
    std::cout << "\tExact lengths sum = " << exact_length_sum << std::endl;

     return 0;
}

int main(int argc, char* argv[])
{
  const std::vector<Point_3> points = { { -2,  0,  0 },
                                        {  2,  0,  0 },
                                        {  0,  1, -1 },
                                        {  0, -1, -1 },
                                        {  0,  0,  1 },
                                        {  -10, -10, -10  },
                                        {  -10, 10, -10   },
                                        {  10, 10, -10    },
                                        {  10, -10, -10   },
                                        {  -10, -10, 10   },
                                        {  -10, 10, 10    },
                                        {  10, 10, 10     },
                                        {  10, -10, 10    },
  };
  std::vector<DT::Vertex_handle> vertices;
  vertices.reserve(points.size());
  DT dt;
  for (auto p : points) vertices.push_back(dt.insert(p));
  DT::Cell_handle c;
  assert(dt.is_valid());
  assert(dt.is_cell(vertices[0], vertices[2], vertices[3], vertices[4], c));
  assert(dt.is_cell(vertices[1], vertices[2], vertices[3], vertices[4], c));

  std::cerr << dt.number_of_finite_cells() << '\n';

  const std::vector < std::tuple<Point_3, Point_3>> queries = {
        {{-1, 0,  0}, { 1, 0,  0}}, // CFC
        {{-1, 0,  0}, { 2, 0,  0}}, // CFCV
        {{ 2, 0,  0}, {-1, 0,  0}}, // reverse
        {{-2, 0,  0}, { 2, 0,  0}}, // VCFCV
        {{ 2, 0,  0}, {-2, 0,  0}}, // reverse case: VCFCV
        {{-3, 0,  0}, { 3, 0,  0}}, // FVCFCVF
        {{-2, 0,  0}, { 2, 2, -2}}, // VEVF
        {{ 2, 2, -2}, {-2, 0,  0}}, // reverse case: FVEV
  };
  bool ok = true;
  for (const auto& tuple : queries) {
    if (!test(dt, tuple)) ok = false;
  }
  std::cout << "Done (" << queries.size() << " queries)\n";
  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}