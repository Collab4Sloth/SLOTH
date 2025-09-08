/**
 * @file MATimerMPI.hxx
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
#pragma once
#include <string>
#include <vector>

#include "MAToolsProfiling/EnumTimer.hxx"
#include "MAToolsProfiling/MATimerNode.hxx"
#include "MAToolsProfiling/MATimers.hxx"
#include "MAToolsProfiling/MATimersFullTreeMode.hxx"

namespace MATools {
namespace MATimer {
namespace FullTreeMode {
using namespace MATools::MATimer;

void rec_call(MATimerNode* a_node, minimal_info* a_ptr) {
  int elem = a_ptr->m_nb_daughter;
  a_node->inc_mpi();

  for (int it = 0; it < elem; it++) {
    auto ptr = a_ptr + it + 1;
    std::string name = ptr->m_name;
    auto node = a_node->find(name);
    rec_call(node, ptr);
  }
}

void transform_to_MATimerMPI(std::vector<minimal_info>& a_in, std::vector<int> a_sizes,
                             int a_mpi_size) {
  MATimerNode*& root = get_MATimer_node<enumTimer::ROOT>();
  int acc = 0;
  int rank = MPI::get_rank();
  for (int mpi = 0; mpi < a_mpi_size; mpi++) {
    if (a_sizes[mpi] == 0) continue;
    auto local_root = a_in.data() + acc;
    rec_call(root, local_root);
    acc += a_sizes[mpi];
  }
}
};  // namespace FullTreeMode
};  // namespace MATimer
};  // namespace MATools
