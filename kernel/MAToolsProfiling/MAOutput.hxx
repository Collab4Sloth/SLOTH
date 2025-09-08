/**
 * @file MAOutput.hxx
 * @author Raphaël Prat (raphael.prat@cea.fr)
 * @brief
 * @version 0.1
 * @date 2025-09-08
 *
 * Copyright CEA (C) 2025
 *
 * This file is part of SLOTH.
 *
 * SLOTH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SLOTH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <MAToolsProfiling/MAToolsMPI.hxx>
#include <iostream>
#include <mfem.hpp>

namespace MATools {
namespace MAOutput {
using namespace MATools::MPI;

/**
 * @brief Displays one message if the current mpi rank is 0
 */
template <typename Arg>
void printMessage(Arg a_msg) {
  mfem::out << a_msg << std::endl;
  /*
          if(is_master())
          {
            std::cout <<" part2" << a_msg << std::endl;
          }
  */
}

/**
 * @brief Displays some messages if the current mpi rank is the master.
 */
template <typename Arg, typename... Args>
void printMessage(Arg a_msg, Args... a_msgs) {
  mfem::out << a_msg << " ";
  printMessage(a_msgs...);
  // if(is_master())
  //{
  //   std::cout <<" part2:"<< a_msg << " ";
  //   printMessage(a_msgs...);
  // }
}
}  // namespace MAOutput
}  // namespace MATools
