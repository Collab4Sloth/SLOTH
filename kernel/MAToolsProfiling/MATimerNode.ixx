/**
 * @file MATimerNode.ixx
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
#include <MAToolsProfiling/MATimerNode.hxx>
#include <algorithm>
#include <iomanip>

namespace MATools {
/**
 * MATimerNode is the storage class corresponding to a node of the MATimer tree.
 */
namespace MATimer {
using duration = std::chrono::duration<double>;
template <>
MATimerNode*& get_MATimer_node<enumTimer::CURRENT>() {
  static MATimerNode* __current;
  return __current;
}

template <>
MATimerNode*& get_MATimer_node<enumTimer::ROOT>() {
  static MATimerNode* __current;
  return __current;
}

/**
 * @brief default constructor.
 */
MATimerNode::MATimerNode()
    : m_daughter()  // only used for root
{
  m_name = "root";
  m_iteration = 1;
  m_level = 0;
  m_mother = nullptr;
  m_nb_mpi = 1;
}

/**
 * @brief MATimerNode constructor used to initialize a node with a node name and a mother node.
 */
MATimerNode::MATimerNode(std::string name, MATimerNode* mother) : m_daughter(), m_duration(0) {
  m_name = name;
  m_iteration = 0;
  m_level = mother->m_level + 1;
  m_mother = mother;
  m_nb_mpi = 1;
}

/**
 * @brief Updates the iteration count.
 * This function updates the iteration count of the MATimerNode by adding the provided count value.
 * @param a_count The count value to add. Default value is 1.
 */
void MATimerNode::update_count() { m_iteration++; }

/**
 * @brief This function is used to find if a daughter node is already defined with this node name.
 * If this node does not exist, a new daughter MATimerNode is added
 * @param[in] name name of the desired node
 * @return the MATimerNode desired
 */
MATimerNode* MATimerNode::find(const std::string name) {
  struct comp {
    comp(const std::string a_value) : m_value(a_value) {}

    inline bool operator()(MATimerNode* a_item) {
      if (a_item->get_name() == m_value)
        return true;
      else
        return false;
    }
    std::string m_value;
  };

  assert(this != nullptr);
  auto it = std::find_if(m_daughter.begin(), m_daughter.end(), comp(name));
  if (it != m_daughter.end()) {
    return (*it);
  }

  /*				for(auto it = m_daughter.begin() ; it < m_daughter.end() ; it++)
                                  {
                                          if((*it)->m_name == name)
                                          {
                                                  (*it)->m_iteration++;
                                                  return (*it);
                                          }
                                  }
  */
  MATimerNode* myTmp = new MATimerNode(name, this);
  m_daughter.push_back(myTmp);

  return myTmp;
}

/**
 * @brief Displays a motif several times.
 * @param[in] begin column number where the motif starts.
 * @param[in] end column number where the motif finishs.
 * @param[in] motif the replicated motif, this motif should have a length equal to 1.
 */
void MATimerNode::print_replicate(int a_begin, int a_end, std::string a_motif) {
#ifndef NDEBUG
  if (a_end < a_begin) {
    // MATools::MAOutput::printMessage("MATools_LOG: column width is probalby too small and it's not
    // possible to replicated this motif");
  }
#endif

  for (int i = a_begin; i < a_end; i++) {
    std::cout << a_motif;
  }
}

/**
 * @brief Displays a blank character.
 */
void MATimerNode::space() { std::cout << " "; }

/**
 * @brief Displays a "|".
 */
void MATimerNode::column() { std::cout << "|"; }

/**
 * @brief Displays a return line.
 */
void MATimerNode::end_line() { std::cout << std::endl; }

/**
 * @brief Displays the banner/header.
 * @param[in] shift number of blank character displayed
 */
void MATimerNode::print_banner(int shift) {
  if (m_name == "root") {
    if (MATools::MPI::is_master()) {
      MATools::MAOutput::printMessage(" MPI feature activated, rank 0:");
      std::string start_name = " |-- start timetable ";
      std::cout << start_name;
      int end = shift + nColumns * (cWidth + 1) + 1;
      print_replicate(start_name.size(), end, "-");
      column();
      end_line();
      std::string name = " |    name";
      std::cout << name;
      print_replicate(name.size(), shift + 1, " ");

      // columns name are displayed, columns change if the MPI mode is ued
      for (int i = 0; i < nColumns; i++) {
        column();
        int size = cName[i].size();
        print_replicate(0, (int(cWidth) - size - 1), " ");
        std::cout << cName[i];
        space();
      }

      column();
      end_line();
      space();
      column();
      print_replicate(2, end, "-");
      column();
      end_line();
    }
  }
}

/**
 * @brief Displays the header.
 * @param[in] shift number of blank character displayed
 */
void MATimerNode::print_ending(int shift) {
  if (m_name == "root") {
    if (MATools::MPI::is_master()) {
      shift += nColumns * (cWidth + 1) + 1;  // +1 for "|";
      std::string end_name = " |-- end timetable ";
      std::cout << end_name;
      print_replicate(end_name.size(), shift, "-");
      column();
      end_line();
    }
  }
}

/**
 * @brief Gets of the duration member
 * @return pointer of the duration member of a MATimerNode
 */
duration* MATimerNode::get_ptr_duration() { return &m_duration; }

/**
 * @brief Displays the local runtime.
 * @param[in] shift number of blank character displayed
 * @param[in] runtime local duration value
 */
void MATimerNode::print_local(int shift, double total_time) {
  assert(total_time >= 0);
  int nC = 3;
  std::string cValue[3];
  int realShift = shift;
  space();
  column();
  space();
  int currentShift = 3;
  for (int i = 0; i < int(m_level) - 1; i++) {
    int spaceSize = 3;
    for (int j = 0; j < spaceSize; j++) space();
    currentShift += spaceSize;
  }
  if (m_level > 0) {
    std::cout << "|--";
    currentShift += 3;
  }
  std::cout << "> " << m_name;
  currentShift += m_name.size() + 1;
  print_replicate(currentShift, realShift, " ");

  cValue[0] = std::to_string(m_iteration);
  cValue[1] = std::to_string(m_duration.count());
  cValue[2] = std::to_string((m_duration.count() / total_time) * 100);

  for (int i = 0; i < nC; i++) {
    column();
    int size = cValue[i].size();
    print_replicate(0, (int(cWidth) - size - 1), " ");
    std::cout << cValue[i];
    space();
  }
  column();
  end_line();
}

/**
 * @brief Displays the runtime.
 * @param[in] shift number of blank character displayed
 * @param[in] runtime duration value
 */
void MATimerNode::print(int shift, double total_time) {
  assert(total_time >= 0);
  std::string cValue[nColumns];
  if (MATools::MPI::is_master()) {
    int realShift = shift;
    space();
    column();
    space();
    int currentShift = 3;
    for (int i = 0; i < int(m_level) - 1; i++) {
      int spaceSize = 3;
      for (int j = 0; j < spaceSize; j++) space();
      currentShift += spaceSize;
    }
    if (m_level > 0) {
      std::cout << "|--";
      currentShift += 3;
    }
    std::cout << "> " << m_name;
    currentShift += m_name.size() + 1;
    print_replicate(currentShift, realShift, " ");

    cValue[0] = std::to_string(m_iteration);
  }
  double local = m_duration.count();
  int size = MATools::MPI::get_mpi_size();

  assert(size > 0);
  std::vector<double> list;

  if (MATools::MPI::is_master()) list.resize(size);

  MPI_Gather(&local, 1, MPI_DOUBLE, list.data(), 1, MPI_DOUBLE, 0,
             MPI_COMM_WORLD);  // master rank is 0

  if (MATools::MPI::is_master()) {
#ifdef UNDEFINED  // __cplusplus > 201103L
    const auto [min, max] = std::minmax_element(list.begin(), list.end());
    auto global_max = *max;
    auto global_min = *min;
    auto sum = std::accumulate(list.begin(), list.end(), double(0.));
    auto global_mean = sum / double(size);
    auto part_time = (global_max / total_time) * 100;
#else
    double min = list[0];
    double max = list[0];
    double sum = list[0];
    int size = list.size();
    for (int id = 1; id < size; id++) {
      min = std::min(min, list[id]);
      max = std::max(max, list[id]);
      sum += list[id];
    }
    auto global_max = max;
    auto global_min = min;
    auto global_mean = sum / double(size);
    auto part_time = (global_max / total_time) * 100;
#endif
    assert(global_mean >= 0);
    assert(global_min >= 0);
    assert(global_max >= 0);

    assert(global_max >= global_mean);
    assert(global_mean >= global_min);

    std::cout << std::setprecision(25);
    const int precisionVal = 25;

    cValue[1] = std::to_string(global_min)
                    .substr(0, std::to_string(global_min).find(".") + precisionVal + 1);
    cValue[2] = std::to_string(global_mean)
                    .substr(0, std::to_string(global_mean).find(".") + precisionVal + 1);
    cValue[3] = std::to_string(global_max)
                    .substr(0, std::to_string(global_max).find(".") + precisionVal + 1);
    cValue[4] = std::to_string(part_time).substr(
                    0, std::to_string(part_time).find(".") + precisionVal + 1) +
                "%";
    cValue[5] = std::to_string((global_max / global_mean) - 1)
                    .substr(0, std::to_string((global_max / global_mean) - 1).find(".") +
                                   precisionVal + 1) +
                "%";
  }
  if (MATools::MPI::is_master()) {
    for (int i = 0; i < nColumns; i++) {
      column();
      int size = cValue[i].size();
      print_replicate(0, (int(cWidth) - size - 1), " ");
      std::cout << cValue[i];
      space();
    }
    column();
    end_line();
  }
}

/**
 * @brief Retruns the MATimerNode name
 * @return name
 */
std::string MATimerNode::get_name() {
  std::string ret = m_name;
  return ret;
}

/**
 * @brief Retruns the duration
 * @return duration value
 */
double MATimerNode::get_duration() {
  double ret = m_duration.count();
  return ret;
}

/**
 * @brief Retruns the MATimerNode iteration number
 * @return the iteration number
 */
int MATimerNode::get_iteration() {
  int res = m_iteration;
  return res;
}

/**
 * @brief Retruns the MATimerNode level
 * @return level
 */
int MATimerNode::get_level() {
  int res = m_level;
  return res;
}

/**
 * @brief Retruns a vector of daughter MATimerNode pointers
 * @return daughter nodes
 */
std::vector<MATimerNode*>& MATimerNode::get_daughter() { return m_daughter; }

/**
 * @brief Retruns the mother MATimerNode pointer
 * @return mother pointer
 */
MATimerNode* MATimerNode::get_mother() {
  MATimerNode* ret = m_mother;
  return ret;
}

/** @brief This function displays information about one timer node */
void MATimerNode::debug_info() {
  std::cout << " node name : " << get_name() << std::endl;
  std::cout << " mother node address : " << get_mother() << std::endl;
  std::cout << " number of daughter nodes : " << m_daughter.size() << std::endl;
  std::cout << " level : " << get_level() << std::endl;
  std::cout << " number of iterations : " << get_iteration() << std::endl;
  std::cout << " current duration : " << get_duration() << std::endl;
}
}  // namespace MATimer
};  // namespace MATools
