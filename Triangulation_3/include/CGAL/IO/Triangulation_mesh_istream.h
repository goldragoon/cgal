// Copyright (c) 2017  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:  $
// $Id:  $
//
// Author(s)     : Jane Tournois


#ifndef CGAL_TRIANGULATION_MESH_ISTREAM_3_H
#define CGAL_TRIANGULATION_OFF_ISTREAM_3_H

#include <CGAL/license/Triangulation.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Mesh_3/tet_soup_to_c3t3.h>


#include <sstream>
#include <iostream>

namespace CGAL {

  template < class GT, class TDS >
  std::istream &
  import_triangulation_3_from_mesh(std::istream & is,
                                   Triangulation_3<GT, TDS> & tr)
  {
    typedef Triangulation_3<GT, TDS> Tr;
    CGAL::build_triangulation_from_file<Tr, true>(is, tr);

    assert(tr.is_valid());

    return is;
  }

}

#endif //CGAL_TRIANGULATION_OFF_ISTREAM_3_H
