
#define CGAL_MESH_3_VERBOSE 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Image_3.h>
#include <functional>
#include <CGAL/Random.h>


typedef short Image_word_type;

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;


class More {
  double iso;
public:
  More(double iso) : iso(iso) {}

  template <typename T>
  int operator()(T v) const {
    return int(v > iso);
  }
};

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main(int argc, char*argv[])
{
  const char* fname = (argc>1)?argv[1]:"data/skull_2.9.inr";

  std::cout << "Mesh " << fname << std::endl;

  Image_word_type iso = (argc > 2) ? boost::lexical_cast<Image_word_type>(argv[2]) : 1;
  double fs = (argc > 3) ? boost::lexical_cast<double>(argv[3]) : 1;
  double fd = (argc > 4) ? boost::lexical_cast<double>(argv[4]) : 0.1;
  double cs = (argc > 5) ? boost::lexical_cast<double>(argv[5]) : 1;

  // Load image
  CGAL::Image_3 image;
  if(!image.read(fname)){
    std::cerr << "Error: Cannot read file " <<  fname << std::endl;
    return EXIT_FAILURE;
  }

  // Mesh criteria
  Mesh_criteria criteria(facet_angle = 30, facet_size = fs, facet_distance = fd,
                         cell_radius_edge_ratio = 3, cell_size = cs);

  CGAL::Random rng;
//  while(true)
  {
    //int r = rng.get_int(10, 1000000);
    //std::cout << "Run # " << r << std::endl;
    //CGAL::get_default_random() = CGAL::Random(r);

    std::cout << "\tSeed is\t"
      << CGAL::get_default_random().get_seed() << std::endl;

    /// [Domain creation]
    Mesh_domain domain = Mesh_domain::create_gray_image_mesh_domain
    (image,
      image_values_to_subdomain_indices = More(iso),
      value_outside = -1000,
      p_rng = &CGAL::get_default_random());
      /// [Domain creation]

    // Meshing
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria, no_perturb(), no_exude());

//    std::cout << "#" << r << " done" << std::endl;
  }


//  // Output
//  std::ofstream medit_file("out.mesh");
//  c3t3.output_to_medit(medit_file);
//
  return 0;
}
