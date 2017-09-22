//#define CGAL_TRIANGULATION_3_TRAVERSER_CHECK_INTERSECTION
#define CGAL_EXPERIMENT_WITH_SIMPLE_CARTESIAN

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


//#define CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE

// Define the kernel.
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::Point_3                                      Point_3;

// Define the structure.
typedef CGAL::Delaunay_triangulation_3< Kernel >    DT;
typedef DT::Cell_handle                             Cell_handle;
typedef DT::Vertex_handle                           Vertex_handle;

typedef CGAL::Triangulation_segment_cell_iterator_3<DT> Cell_traverser;

template <typename PointVector>
void scale_points(PointVector& points, const double scale)
{
  //bbox
  double xmin = points[0].x();
  double xmax = points[0].x();
  double ymin = points[0].y();
  double ymax = points[0].y();
  double zmin = points[0].z();
  double zmax = points[0].z();

  BOOST_FOREACH(Point_3 p, points)
  {
    xmin = (std::min)(xmin, p.x());
    ymin = (std::min)(ymin, p.y());
    zmin = (std::min)(zmin, p.z());
    xmax = (std::max)(xmax, p.x());
    ymax = (std::max)(ymax, p.y());
    zmax = (std::max)(zmax, p.z());
  }
  std::cout << "Bbox is [" << xmin << "; " << xmax << "] "
    << "[" << ymin << "; " << ymax << "] "
    << "[" << zmin << "; " << zmax << "]" << std::endl;

  double dx = xmax - xmin;
  double dy = ymax - ymin;
  double dz = zmax - zmin;

  for (std::size_t i = 0; i < points.size(); ++i)
  {
    const Point_3 pi = points[i];
    double nx = (pi.x() - xmin) * scale / dx;
    double ny = (pi.y() - ymin) * scale / dy;
    double nz = (pi.z() - zmin) * scale / dz;

    points[i] = Point_3(nx, ny, nz);
  }
}

int main(int argc, char* argv[])
{
  const char* fname = (argc>1) ? argv[1] : "data/blobby.xyz";
  unsigned int nb_seg = (argc > 2) ? static_cast<unsigned int>(atoi(argv[2]))
                                   : 100;

  // Reads a .xyz point set file in points.
  // As the point is the second element of the tuple (that is with index 1)
  // we use a property map that accesses the 1st element of the tuple.

  std::vector<Point_3> points;
  std::ifstream stream(fname);
  if (!stream ||
      !CGAL::read_xyz_points(stream, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  scale_points(points, 100);

  //bbox
  double xmin = points[0].x();
  double xmax = points[0].x();
  double ymin = points[0].y();
  double ymax = points[0].y();
  double zmin = points[0].z();
  double zmax = points[0].z();

  BOOST_FOREACH(Point_3 p, points)
  {
    xmin = (std::min)(xmin, p.x());
    ymin = (std::min)(ymin, p.y());
    zmin = (std::min)(zmin, p.z());
    xmax = (std::max)(xmax, p.x());
    ymax = (std::max)(ymax, p.y());
    zmax = (std::max)(zmax, p.z());
  }
  std::cout << "Bbox is [" << xmin << "; " << xmax << "] "
    << "[" << ymin << "; " << ymax << "] "
    << "[" << zmin << "; " << zmax << "]" << std::endl;

    // Construct the Delaunay triangulation.
    DT dt( points.begin(), points.end() );
    assert( dt.is_valid() );

    Vertex_handle v3 = dt.insert(Point_3(6, 6, 6));
    Vertex_handle v4 = dt.insert(Point_3(10, 10, 10));

    Cell_handle c;
    int i, j;
    if (dt.is_edge(v3, v4, c, i, j))
    {
      Point_3 p1(0, 0, 0);
      Point_3 p2(24, 24, 24);

      std::cout << "Traverser through v3 and v4" << std::endl;
      Cell_traverser ct(dt, p1, p2);

      unsigned int inf = 0, fin = 0;
      for (; ct != ct.end(); ++ct)
      {
        if (dt.is_infinite(ct))
          ++inf;
        else
        {
          ++fin;
          DT::Locate_type lt;
          int li, lj;
          ct.entry(lt, li, lj);
          std::cout << "\t entry type : " << lt << std::endl;
        }
      }
      std::cout << "Done : " << inf << " " << fin << std::endl;
    }
    else
    {
      std::cout << "[v3 ; v4] IS NOT AN EDGE" << std::endl;
    }
    CGAL::default_random = CGAL::Random(0);
    CGAL::Random rng(0);
    CGAL::Timer time;
    time.start();

    unsigned int nb_facets = 0, nb_edges = 0, nb_vertex = 0;
    for (unsigned int i = 0; i < nb_seg; ++i)
    {
      // Construct a traverser.
      Point_3 p1(rng.get_double(xmin, xmax),
                 rng.get_double(ymin, ymax),
                 rng.get_double(zmin, zmax));
      Point_3 p2(rng.get_double(xmin, xmax),
                 rng.get_double(ymin, ymax),
                 rng.get_double(zmin, zmax));

#ifdef CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE
      std::cout << "Traverser " << (i + 1)
        << "\n\t(" << p1
        << ")\n\t(" << p2 << ")" << std::endl;
#endif
      Cell_traverser ct(dt, p1, p2);

      // Count the number of finite cells traversed.
      unsigned int inf = 0, fin = 0;
      for( ; ct != ct.end(); ++ct )
      {
        if( dt.is_infinite(ct) )
            ++inf;
        else
        {
          ++fin;

          //DT::Locate_type lt;
          //int li, lj;
          //ct.entry(lt, li, lj);

          //switch (lt)
          //{
          //case DT::Locate_type::FACET:
          //  ++nb_facets;
          //  break;
          //case DT::Locate_type::EDGE:
          //  ++nb_edges;
          //  break;
          //case DT::Locate_type::VERTEX:
          //  ++nb_vertex;
          //  break;
          //default:
          //  /*when source is in a cell*/
          //  CGAL_assertion(lt == DT::Locate_type::CELL);
          //}
        }
      }

#ifdef CGAL_TRIANGULATION_3_VERBOSE_TRAVERSER_EXAMPLE
      std::cout << "While traversing from " << ct.source()
                << " to " << ct.target() << std::endl;
      std::cout << inf << " infinite and "
                << fin << " finite cells were visited." << std::endl;
      std::cout << std::endl << std::endl;
#endif
    }

    time.stop();
    std::cout << "Triangulation has " << dt.number_of_vertices() << " vertices."
      << std::endl;
    std::cout << "Traversing cells of triangulation with "
      << nb_seg << " segments took " << time.time() << " seconds."
      << std::endl;
    std::cout << "\tnb facets    : " << nb_facets << std::endl;
    std::cout << "\tnb edges     : " << nb_edges << std::endl;
    std::cout << "\tnb vertices  : " << nb_vertex << std::endl;

     return 0;
}