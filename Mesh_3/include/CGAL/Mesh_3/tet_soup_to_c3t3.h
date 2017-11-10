// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Mael Rouxel-Labbé, Maxime Gimeno
//
//******************************************************************************
//
//******************************************************************************

#ifndef CGAL_MESH_3_TET_SOUP_TO_C3T3_H
#define CGAL_MESH_3_TET_SOUP_TO_C3T3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/IO/File_medit.h>

#include <CGAL/assertions.h>

#include <boost/array.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <boost/foreach.hpp>


namespace CGAL
{
template<class Tr>
void build_vertices(Tr& tr,
                    const std::vector<typename Tr::Point>& points,
                    std::vector<typename Tr::Vertex_handle>& vertex_handle_vector)
{
  std::cout << "build_vertices...(" << points.size() << " points)" << std::endl;
  typedef typename Tr::Vertex_handle            Vertex_handle;

  vertex_handle_vector[0] = tr.tds().create_vertex(); // creates the infinite vertex
  tr.infinite_vertex() = vertex_handle_vector[0];
  CGAL_assertion(vertex_handle_vector[0] == tr.infinite_vertex());

  // build vertices
  for(std::size_t i=0; i<points.size(); ++i)
  {
    Vertex_handle vh = tr.tds().create_vertex();
    vertex_handle_vector[i+1] = vh;
    vh->set_point(points[i]);
  }

  CGAL_postcondition(vertex_handle_vector[0] == tr.infinite_vertex());
  for (std::size_t i = 0; i < points.size(); ++i)
  {
    CGAL_postcondition(vertex_handle_vector[i + 1] != tr.infinite_vertex());
  }
}

template<class Tr>
void add_facet_to_incident_cells_map(const typename Tr::Cell_handle c,
                                     int i,
                                     boost::unordered_map<std::set<typename Tr::Vertex_handle>,
                                                          std::vector<std::pair<typename Tr::Cell_handle,
                                                                                int> > >& incident_cells_map)
{
  typedef typename Tr::Vertex_handle                                Vertex_handle;
  typedef typename Tr::Cell_handle                                  Cell_handle;
  typedef std::set<Vertex_handle>                                   Facet;
  typedef std::pair<Cell_handle, int>                               Incident_cell;
  typedef boost::unordered_map<Facet, std::vector<Incident_cell> >  Incident_cells_map;

  // the opposite vertex of f in c is i
  Facet f;
  f.insert(c->vertex((i + 1) % 4));
  f.insert(c->vertex((i + 2) % 4));
  f.insert(c->vertex((i + 3) % 4));
  std::cout << "FACET is " << std::endl;
  std::cout << "[" << i << "\t" << &*(c->vertex(i)) << "]" << std::endl;
  std::cout << ((i + 1) % 4) << "\t" << &*(c->vertex((i + 1) % 4)) << std::endl;
  std::cout << ((i + 2) % 4) << "\t" << &*(c->vertex((i + 2) % 4)) << std::endl;
  std::cout << ((i + 3) % 4) << "\t" << &*(c->vertex((i + 3) % 4)) << std::endl;

  std::set<Vertex_handle> debug_v;
  debug_v.insert(c->vertex(0));
  debug_v.insert(c->vertex(1));
  debug_v.insert(c->vertex(2));
  debug_v.insert(c->vertex(3));

  CGAL_assertion(debug_v.size() == 4);
  CGAL_precondition(f.size() == 3);

  Incident_cell e = std::make_pair(c, i);

  typename Incident_cells_map::iterator fit = incident_cells_map.find(f);
  if(fit != incident_cells_map.end())  // the entry already exists in the map
  {
    // a facet must have exactly two incident cells
    CGAL_assertion(fit->second.size() == 1);
    fit->second.push_back(e);

    std::vector<Incident_cell> vec2 = fit->second;
    std::cout << &*vec2[0].first << "\t" << vec2[0].second << std::endl;
    std::cout << &*vec2[1].first << "\t" << vec2[1].second << std::endl;
    std::cout << &*e.first << "\t" << e.second << std::endl;
  }
  else
  {
    std::vector<Incident_cell> vec(1);
    vec[0] = e;
    incident_cells_map.insert(std::make_pair(f, vec));
  }

}

template<class Tr>
void build_finite_cells(Tr& tr,
    const std::vector<boost::array<int,5> >& finite_cells,
    const std::vector<typename Tr::Vertex_handle>& vertex_handle_vector,
    boost::unordered_map<std::set<typename Tr::Vertex_handle>,
                         std::vector<std::pair<typename Tr::Cell_handle, int> > >& incident_cells_map,
    const std::map<boost::array<int,3>, int>& border_facets)
{
  typedef boost::array<int, 5>              Tet_with_ref; // 4 ids + 1 reference

  typedef typename Tr::Vertex_handle                            Vertex_handle;
  typedef typename Tr::Cell_handle                              Cell_handle;

  CGAL_assertion_code(
    typename Tr::Geom_traits::Construct_point_3 wp2p =
      tr.geom_traits().construct_point_3_object();
  )

  // build the finite cells
  for(std::size_t i=0; i<finite_cells.size(); ++i)
  {
    const Tet_with_ref& tet = finite_cells[i];
    boost::array<Vertex_handle, 4> vs;

    for(int j=0; j<4; ++j)
    {
      CGAL_precondition(static_cast<std::size_t>(tet[j]) < tr.number_of_vertices() &&
                        tet[j] >= 0);
      //vertex_handle_vector[tet[j] + 1]->set_dimension(3);
      vs[j] = vertex_handle_vector[tet[j] + 1];
      CGAL_postcondition(vs[j] != Vertex_handle());
      CGAL_postcondition(vs[j] != tr.infinite_vertex());
    }

    // this assertion also tests for degeneracy
    CGAL_assertion(CGAL::orientation(wp2p(vs[0]->point()), wp2p(vs[1]->point()),
                                     wp2p(vs[2]->point()), wp2p(vs[3]->point()))
                     == POSITIVE);

    Cell_handle c = tr.tds().create_cell(vs[0], vs[1], vs[2], vs[3]);
    //c->info() = tet[4]; // the cell's info keeps the reference of the tetrahedron

    CGAL_precondition(tet[4] > 0);
    // assign cells to vertices
    for(int j=0; j<4; ++j)
    {
      if(vs[j]->cell() == Cell_handle())
        vs[j]->set_cell(c);
    }

    // build the map used for adjacency later
    for(int j=0; j<4; ++j)
    {
      add_facet_to_incident_cells_map<Tr>(c, j, incident_cells_map);
      //if(border_facets.size() !=0)
      //{
      //  boost::array<int,3> facet;
      //  facet[0]=tet[(j+ 1) % 4];
      //  facet[1]=tet[(j+ 2) % 4];
      //  facet[2]=tet[(j+ 3) % 4];
      //  //find the circular permutation that puts the smallest index in the first place.
      //  int n0 = (std::min)((std::min)(facet[0],facet[1]), facet[2]);
      //  int k=0;
      //  boost::array<int,3> f;
      //  do
      //  {
      //    f[0]=facet[(0+k)%3];
      //    f[1]=facet[(1+k)%3];
      //    f[2]=facet[(2+k)%3];
      //    ++k;
      //  }while(f[0] != n0);
      //  if(border_facets.find(f) != border_facets.end())
      //    c->set_surface_patch_index(j, border_facets.at(f));
      //  else
      //  {
      //    int temp = f[2];
      //    f[2]=f[1];
      //    f[1]=temp;
      //    if(border_facets.find(f) != border_facets.end())
      //      c->set_surface_patch_index(j, border_facets.at(f));
      //    else
      //      c->set_surface_patch_index(j, 0);
      //  }
      //}
    }
  }
}

template<class Tr>
void add_infinite_facets_to_incident_cells_map(typename Tr::Cell_handle c,
           int inf_vert_pos,
           boost::unordered_map<std::set<typename Tr::Vertex_handle>,
                                std::vector<std::pair<typename Tr::Cell_handle,
                                                      int> > >& incident_cells_map)
{
  std::set<typename Tr::Vertex_handle> debug_v;
  debug_v.insert(c->vertex(0));
  debug_v.insert(c->vertex(1));
  debug_v.insert(c->vertex(2));
  debug_v.insert(c->vertex(3));
  CGAL_assertion(debug_v.size() == 4);

  int l = (inf_vert_pos + 1) % 4;
  add_facet_to_incident_cells_map<Tr>(c, l, incident_cells_map);
  l = (inf_vert_pos + 2) % 4;
  add_facet_to_incident_cells_map<Tr>(c, l, incident_cells_map);
  l = (inf_vert_pos + 3) % 4;
  add_facet_to_incident_cells_map<Tr>(c, l, incident_cells_map);
}

template<class Tr>
void build_infinite_cells(Tr& tr,
                          boost::unordered_map<std::set<typename Tr::Vertex_handle>,
                                               std::vector<std::pair<typename Tr::Cell_handle,
                                                                     int> > >& incident_cells_map)
{
  std::cout << "** build_infinite_cells **" << std::endl;
  typedef typename Tr::Vertex_handle                               Vertex_handle;
  typedef typename Tr::Cell_handle                                 Cell_handle;
  typedef std::set<Vertex_handle>                                  Facet;
  typedef std::pair<Cell_handle, int>                              Incident_cell;
  typedef boost::unordered_map<Facet, std::vector<Incident_cell> > Incident_cells_map;

  std::vector<Cell_handle> infinite_cells;

  // check the incident cells map for facets who only have one incident cell
  // and build the infinite cell on the opposite side
  typename Incident_cells_map::iterator it = incident_cells_map.begin();
  const typename Incident_cells_map::iterator end = incident_cells_map.end();
  for(; it != end; ++it)
  {
    std::vector<Incident_cell>& inc_to_f = it->second;
    if(inc_to_f.size() == 2) // facet already has both its incident cells
      continue;
    CGAL_assertion(inc_to_f.size() == 1);

    Cell_handle c = inc_to_f[0].first;
    int i = inc_to_f[0].second;

    std::set<Vertex_handle> debug_v;
    bool f_infinite = false;
    for (std::size_t i = 0; i < 4; ++i)
    {
      if (c->vertex(i) == tr.infinite_vertex())
      {
        std::cout << "vertex " << i << " is INFINITE" << std::endl;
        f_infinite = true;
        break;
      }
      debug_v.insert(c->vertex(i));
    }
    if (f_infinite)
      continue;
    if (debug_v.size() != 4)
    {
      std::cout << "C is degenerate!" << std::endl;
      std::cout << "i = " << i << std::endl;
    }

    Cell_handle opp_c;
    // the infinite cell that we are creating needs to be well oriented...1
    int inf_vert_position_in_opp_c = 0;
    if(i == 0 || i == 2)
      opp_c = tr.tds().create_cell(tr.infinite_vertex(),
                                   c->vertex((i+1)%4),
                                   c->vertex((i+2)%4),
                                   c->vertex((i+3)%4));
    else
      opp_c = tr.tds().create_cell(tr.infinite_vertex(),
                                   c->vertex((i+1)%4),
                                   c->vertex((i+3)%4),
                                   c->vertex((i+2)%4));

    debug_v.clear();
    debug_v.insert(opp_c->vertex(0));
    debug_v.insert(opp_c->vertex(1));
    debug_v.insert(opp_c->vertex(2));
    debug_v.insert(opp_c->vertex(3));
    if (debug_v.size() != 4)
    {
      std::cout << "Infinite vertex is : " << std::endl;
      std::cout << &*(tr.infinite_vertex()) << std::endl;
      std::cout << "v0 = " << &*(opp_c->vertex(0)) << std::endl;
      std::cout << "v1 = " << &*(opp_c->vertex(1)) << std::endl;
      std::cout << "v2 = " << &*(opp_c->vertex(2)) << std::endl;
      std::cout << "v3 = " << &*(opp_c->vertex(3)) << std::endl;
    }
    CGAL_assertion(debug_v.size() == 4);

    // set the infinite_vertex's incident cell
    if(tr.infinite_vertex()->cell() == Cell_handle())
    {
      tr.infinite_vertex()->set_cell(opp_c);
      std::cout << "SET_CELL" << std::endl;
    }

    // add the facets to the incident cells map

    // the only finite facet
    inc_to_f.push_back(std::make_pair(opp_c, inf_vert_position_in_opp_c));
    CGAL_assertion(inc_to_f.size() == 2);

    CGAL_assertion(opp_c->vertex(inf_vert_position_in_opp_c) == tr.infinite_vertex());
  }

  for (std::size_t i = 0; i < infinite_cells.size(); ++i)
  {
    Cell_handle cinf = infinite_cells[i];
    CGAL_assertion(cinf->vertex(0) == tr.infinite_vertex());
    add_infinite_facets_to_incident_cells_map<Tr>(cinf, 0, incident_cells_map);
  }
}

template<class Tr>
bool assign_neighbors(Tr& tr,
    const boost::unordered_map<
        std::set<typename Tr::Vertex_handle>,
        std::vector<std::pair<typename Tr::Cell_handle, int> > >& incident_cells_map)
{
  typedef typename Tr::Vertex_handle                                 Vertex_handle;
  typedef typename Tr::Cell_handle                                   Cell_handle;
  typedef std::set<Vertex_handle>                                    Facet;
  typedef std::pair<Cell_handle, int>                                Incident_cell;
  typedef boost::unordered_map<Facet, std::vector<Incident_cell> >  Incident_cells_map;

  // 4 facets per cell, each facet shared by 2 cells
//  if (incident_cells_map.size() != tr.number_of_cells() * 2)
//    return false;

  typename Incident_cells_map::const_iterator icit = incident_cells_map.begin();
  for(; icit!=incident_cells_map.end(); ++icit)
  {
    const std::vector<Incident_cell>& adjacent_cells = icit->second;
    if(adjacent_cells.size() != 2)
      return false;

    Cell_handle c0 = adjacent_cells[0].first;
    int i0 = adjacent_cells[0].second;
    Cell_handle c1 = adjacent_cells[1].first;
    int i1 = adjacent_cells[1].second;

    tr.tds().set_adjacency(c0, i0, c1, i1);
  }
  return true;
}

template<class Tr, bool c3t3_loader_failed>
bool build_triangulation(Tr& tr,
                         const std::vector<typename Tr::Point>& points,
                         const std::vector<boost::array<int,5> >& finite_cells,
                         const std::map<boost::array<int,3>, int>& border_facets)
{
  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Cell_handle              Cell_handle;
  typedef std::set<Vertex_handle>               Facet;

  // associate to a face the two (at most) incident tets and the id of the face in the cell
  typedef std::pair<Cell_handle, int>                   Incident_cell;
  typedef boost::unordered_map<Facet, std::vector<Incident_cell> >  Incident_cells_map;

  Incident_cells_map incident_cells_map;
  std::vector<Vertex_handle> vertex_handle_vector(points.size() + 1); // id to vertex_handle

  CGAL_precondition(!points.empty());

  if(finite_cells.empty())
  {
    std::cout << "WARNING: No finite cells were provided. Only the points will be loaded."<<std::endl;
  }

  tr.tds().clear(); // not tr.clear() since it calls tr.init() which we don't want

  build_vertices<Tr>(tr, points, vertex_handle_vector);
  //BOOST_FOREACH(Vertex_handle vh, vertex_handle_vector)
  //{
  //  vh->set_dimension(-1);
  //}
  if(!finite_cells.empty())
  {
    build_finite_cells<Tr>(tr, finite_cells, vertex_handle_vector, incident_cells_map, border_facets);
    build_infinite_cells<Tr>(tr, incident_cells_map);
    tr.tds().set_dimension(3);
    if(!assign_neighbors<Tr>(tr, incident_cells_map))
      return false;
    std::cout << "built triangulation : " << std::endl;
    std::cout << tr.number_of_cells() << " cells" << std::endl;
  }
  std::cout << tr.number_of_vertices() << " vertices" << std::endl;
  if(c3t3_loader_failed)
  {
    return true;
  }
  else
    return tr.is_valid(true);
}

template<class Tr, bool c3t3_loader_failed>
bool build_triangulation_from_file(std::istream& is,
                                   Tr& tr)
{
  typedef typename Tr::Point                                  Point_3;

  typedef boost::array<int, 3> Facet; // 3 = id
  typedef boost::array<int, 5> Tet_with_ref; // first 4 = id, fifth = reference

  std::vector<Tet_with_ref> finite_cells;
  std::vector<Point_3> points;
  std::map<Facet, int> border_facets;

  // grab the vertices
  int dim;
  int nv, nf, ntet, ref;
  std::string word;

  is >> word >> dim; // MeshVersionFormatted 1
  is >> word >> dim; // Dimension 3

  CGAL_assertion(dim == 3);

  while(is >> word && word != "End")
  {
    if(word == "Vertices")
    {
      is >> nv;
      for(int i=0; i<nv; ++i)
      {
        double x,y,z;
        is >> x >> y >> z>>ref;
        points.push_back(Point_3(x,y,z));
      }
    }

    if(word == "Triangles")
    {
      is >> nf;
      for(int i=0; i<nf; ++i)
      {
        int n1, n2, n3, surface_patch_id;
        is >> n1 >> n2 >> n3 >> surface_patch_id;
        Facet facet;
        facet[0] =n1 -1;
        facet[1] =n2 -1;
        facet[2] =n3 -1;
        //find the circular permutation that puts the smallest index in the first place.
        int n0 = (std::min)((std::min)(facet[0],facet[1]), facet[2]);
        int k=0;
        Facet f;
        do
        {
          f[0]=facet[(0+k)%3];
          f[1]=facet[(1+k)%3];
          f[2]=facet[(2+k)%3];
          ++k;
        }while(f[0] != n0);
        border_facets.insert(std::make_pair(f, surface_patch_id));
      }
    }

    if(word == "Tetrahedra")
    {
      is >> ntet;
      for(int i=0; i<ntet; ++i)
      {
        int n0, n1, n2, n3, reference;
        is >> n0 >> n1 >> n2 >> n3 >> reference;
        Tet_with_ref t;
        t[0] = n0 - 1;
        t[1] = n1 - 1;
        t[2] = n2 - 1;
        t[3] = n3 - 1;
        t[4] = reference;
        finite_cells.push_back(t);
      }
    }
  }
  if(finite_cells.empty())
  {
    return false;
  }

  bool is_well_built = build_triangulation<Tr, c3t3_loader_failed>(tr, points, finite_cells, border_facets);
  return is_well_built;
}

}  // namespace CGAL

#endif // CGAL_MESH_3_TET_SOUP_TO_C3T3_H
