// (LicenseStuffHere)
//
// $URL$
// $Id$
// 
//
// Author(s)     : Stephen Kiazyk

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Random.h>

#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path_traits.h>
#include <CGAL/Polyhedron_shortest_path/Polyhedron_shortest_path.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <iterator>
#include <fstream>

#define UNUSED(X) (void)sizeof(X)

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Polyhedron_shortest_path_default_traits<Kernel, Polyhedron_3> Traits;
typedef CGAL::Polyhedron_shortest_path<Traits> Polyhedron_shortest_path;
typedef boost::graph_traits<Polyhedron_3> GraphTraits;
typedef GraphTraits::vertex_descriptor vertex_descriptor;
typedef GraphTraits::vertex_iterator vertex_iterator;
typedef GraphTraits::halfedge_descriptor halfedge_descriptor;
typedef GraphTraits::halfedge_iterator halfedge_iterator;
typedef GraphTraits::face_descriptor face_descriptor;
typedef GraphTraits::face_iterator face_iterator;

Traits::Barycentric_coordinate random_coordinate(CGAL::Random& rand)
{
  Traits::FT u = rand.uniform_01<Traits::FT>();
  Traits::FT v = rand.uniform_real(Traits::FT(0.0), Traits::FT(1.0) - u);
  return Traits::Barycentric_coordinate(u, v, Traits::FT(1.0) - u - v);
}

int main(int argc, char** argv)
{
  UNUSED(argc);
  UNUSED(argv);

  Polyhedron_3 polyhedron;
  
  std::ifstream inStream("data/elephant.off");
  
  inStream >> polyhedron;
  
  inStream.close();

  face_iterator facesStart, facesEnd;
  boost::tie(facesStart, facesEnd) = CGAL::faces(polyhedron);
  
  std::vector<face_descriptor> faceList;
  
  for (face_iterator facesCurrent = facesStart; facesCurrent != facesEnd; ++facesCurrent)
  {
    faceList.push_back(*facesCurrent);
  }
  
  CGAL::Random rand(2379912);
  const size_t numSamplePoints = 30;
  
  std::vector<Polyhedron_shortest_path::Face_location> faceLocations;
  
  for (size_t i = 0; i < numSamplePoints; ++i)
  {
    faceLocations.push_back(Polyhedron_shortest_path::Face_location(faceList[rand.get_int(0, CGAL::num_faces(polyhedron))], random_coordinate(rand)));
  }
  
  Traits traits;
  Polyhedron_shortest_path shortestPaths(polyhedron, traits);

  shortestPaths.construct_sequence_tree(faceLocations.begin(), faceLocations.end());
  
  vertex_iterator verticesCurrent, verticesEnd;
  
  std::ofstream outPaths("polylines.cgal");
  
  for (boost::tie(verticesCurrent, verticesEnd) = boost::vertices(polyhedron); verticesCurrent != verticesEnd; ++verticesCurrent)
  {
    std::vector<Traits::Point_3> points;
    
    shortestPaths.shortest_path_points_to_source_points(*verticesCurrent, std::back_inserter(points));
    
    outPaths << points.size();
    
    for (size_t i = 0; i < points.size(); ++i)
    {
      outPaths << " " << points[i];
    }
    
    outPaths << std::endl;
  }
  
  outPaths.close();
  
  return 0;
}
