/* This file is part of Vistle.

   You can use it under the terms of the GNU Lesser General Public License
   version 2.1 or later, see lgpl-2.1.txt.

 * License: LGPL 2+ */

#ifndef VISTLE_READANSYS_CONSTANTS_H
#define VISTLE_READANSYS_CONSTANTS_H

// Constants used by both ReadAnsys and ReadRST to avoid circular dependencies
namespace vistle {
namespace ansys {
static const int V_OFFSET = 500;
static const int EX_OFFSET = 1000;
} // namespace ansys
} // namespace vistle

#endif
