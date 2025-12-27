/**
 * @file MATrace.hxx
 * @author Raphaël Prat (raphael.prat@cea.fr)
 * @brief
 * @version 0.1
 * @date 2025-12-27
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

#include <random>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

#include <MAToolsProfiling/MATrace.hxx>
#include <MAToolsProfiling/MAToolsMPI.hxx>
#include <MAToolsProfiling/MATraceColor.hxx>
#include <MAToolsProfiling/MATraceOptional.hxx>
#include <MAToolsProfiling/MATraceInstance.hxx>

namespace MATools {
namespace MATrace {
void start() {
  if (!Optional::is_MATrace_mode()) {
    return;
  }
  auto& start = get_MATrace_point();
  start = MATrace_point();
}

void stop(std::string a_name) {
  if (!Optional::is_MATrace_mode()) {
    return;
  }
  auto end = MATrace_point();
  auto& start = get_MATrace_point();
  auto& ref = get_ref_MATrace_point();
  char name[64];
  strncpy(name, a_name.c_str(), 63);
  name[63]='\0';  // warning issue without this line
  MATrace_section section(name, ref, start, end);
  auto& local_MATrace = get_local_MATrace();
  local_MATrace.push_back(section);
}

void initialize() {
  auto start = get_ref_MATrace_point();
  start.set_proc();
  auto& local_MATrace = get_local_MATrace();
  local_MATrace.reserve(12000000);  // Not really useful
}

void header(std::ofstream& out, vite_event& event) {
  out << "%EventDef PajeDefineContainerType 1" << std::endl;
  out << "% Alias string " << std::endl;
  out << "% ContainerType string " << std::endl;
  out << "% Name string " << std::endl;
  out << "%EndEventDef " << std::endl;
  out << "%EventDef PajeDefineStateType 3" << std::endl;
  out << "% Alias string " << std::endl;
  out << "% ContainerType string " << std::endl;
  out << "% Name string " << std::endl;
  out << "%EndEventDef " << std::endl;
  out << "%EventDef PajeDefineEntityValue 6" << std::endl;
  out << "% Alias string  " << std::endl;
  out << "% EntityType string  " << std::endl;
  out << "% Name string  " << std::endl;
  out << "% Color color " << std::endl;
  out << "%EndEventDef  " << std::endl;
  out << "%EventDef PajeCreateContainer 7" << std::endl;
  out << "% Time date  " << std::endl;
  out << "% Alias string  " << std::endl;
  out << "% Type string  " << std::endl;
  out << "% Container string  " << std::endl;
  out << "% Name string  " << std::endl;
  out << "%EndEventDef  " << std::endl;
  out << "%EventDef PajeDestroyContainer 8" << std::endl;
  out << "% Time date  " << std::endl;
  out << "% Name string  " << std::endl;
  out << "% Type string  " << std::endl;
  out << "%EndEventDef  " << std::endl;
  out << "%EventDef PajeSetState 10 " << std::endl;
  out << "% Time date  " << std::endl;
  out << "% Type string  " << std::endl;
  out << "% Container string  " << std::endl;
  out << "% Value string  " << std::endl;
  out << "%EndEventDef" << std::endl;
  out << "%EventDef PajeDefineVariableType 50 " << std::endl;
  out << "% Alias string" << std::endl;
  out << "% Name  string" << std::endl;
  out << "% ContainerType string " << std::endl;
  out << "%EndEventDef " << std::endl;
  out << "%EventDef PajeSetVariable 51" << std::endl;
  out << "% Time date " << std::endl;
  out << "% Type string " << std::endl;
  out << "% Container string " << std::endl;
  out << "% Value double " << std::endl;
  out << "%EndEventDef  " << std::endl;
  out << "%EventDef PajeAddVariable 52" << std::endl;
  out << "% Time date " << std::endl;
  out << "% Type string " << std::endl;
  out << "% Container string " << std::endl;
  out << "% Value double " << std::endl;
  out << "%EndEventDef  " << std::endl;
  out << "%EventDef PajeSubVariable 53" << std::endl;
  out << "% Time date " << std::endl;
  out << "% Type string " << std::endl;
  out << "% Container string " << std::endl;
  out << "% Value double " << std::endl;
  out << "%EndEventDef" << std::endl;
  event.write_items(out);
  out << "1 CT_Prog   0       'Program'" << std::endl;
  out << "1 CT_Thread CT_Prog 'Thread'" << std::endl;
  out << "3 ST_ThreadState CT_Thread 'Thread State'" << std::endl;
  event.write_colors(out);
  out << "50 V_Sem Semaphore CT_Thread" << std::endl;
  out << "7 0.000000 C_Prog CT_Prog 0 'Programme'" << std::endl;
  int mpi_size = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  for (int i = 0 ; i < mpi_size ; i++) {
    out << "7 0.000000 C_Thread" << i <<" CT_Thread C_Prog 'Thread " << i <<"'" << std::endl;
    out << "51 0.000000 V_Sem C_Thread" << i << " 0.0" << std::endl;
  }
}

void ending(std::ofstream& out, double last) {
  int mpi_size = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  for (int i = 0 ; i < mpi_size ; i++) {
    out << "8 " << last << " C_Thread" << i << " CT_Thread" << std::endl;
  }
  out <<"8 " << last << " C_Prog CT_Prog" << std::endl;
}

void write_trace() {
  //  No MATrace
  if (!Optional::is_MATrace_mode()) {
    return;
  }

  auto& local_MATrace = get_local_MATrace();
  const int local_size = local_MATrace.size();
  int local_byte_size = sizeof(MATrace_section) * local_size;

  // update proc id
  const auto my_rank = get_rank();
  for (auto& it : local_MATrace) {
    it.set_proc(my_rank);
  }
  // all data are sent to the process 0
  // we gather the sizes for each process
  int mpi_size = -1;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  assert(mpi_size >= 0);
  std::vector<int> sizes(mpi_size, 0);
  std::vector<int> dists(mpi_size, 0);
  MPI_Gather(
      &local_byte_size, 1, MPI_INT,
      sizes.data(), 1, MPI_INT,
      0, MPI_COMM_WORLD);

  int total_byte_size = 0;
  int acc = 0;
  for (size_t i = 0 ; i < sizes.size() ; i++) {
    total_byte_size += sizes[i];
    dists[i] = acc;
    acc += sizes[i];
  }
  int total_size = total_byte_size / sizeof(MATrace_section);
  std::vector<char> recv(total_byte_size, 0);
  MPI_Gatherv(
      local_MATrace.data(), local_byte_size, MPI_CHAR,
      recv.data(), sizes.data(), dists.data(), MPI_CHAR,
      0, MPI_COMM_WORLD);


  MATrace_section * ptr = reinterpret_cast<MATrace_section*>(recv.data());

  // write trace header -> trace core -> end
  if (is_master()) {
    std::string filename = "MATrace." + std::to_string(mpi_size) + ".paje";
    std::ofstream out(filename, std::ofstream::out);
    vite_event event;
    for (int it = 0 ; it < total_size ; it++) {
      event.add(ptr[it].m_name);
    }

    header(out, event);

    for (int it = 0 ; it < total_size ; it++) {
      ptr[it].write(out, event);
    }

    // last time
    double last = 0;
    for (int it = 0 ; it < total_size ; it++) {
      last = std::max(last, ptr[it].m_end);
    }

    ending(out, last);
  }
}
}  // namespace MATrace
}  // namespace MATools
