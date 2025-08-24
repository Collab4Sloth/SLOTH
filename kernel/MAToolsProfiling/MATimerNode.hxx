
#pragma once

#include <chrono>
#include <vector>
#include <cassert>
#include <MAToolsProfiling/EnumTimer.hxx>
#include <MAToolsProfiling/Column.hxx>
#include <MAToolsProfiling/MAOutput.hxx>
#include <MAToolsProfiling/MAToolsMPI.hxx>

namespace MATools
{
  namespace MATimer
  {
    /**
     * MATimerNode is the storage class corresponding to a node of the MATimer tree.
     */
    class MATimerNode
    {
      using duration = std::chrono::duration<double>;

      public:

      /**
       * @brief default constructor.
       */
      MATimerNode();
      /**
       * @brief MATimerNode constructor used to initialize a node with a node name and a mother node.
       */
      MATimerNode(std::string name, MATimerNode* mother);

      /**
       * @brief Updates the iteration count.
       * This function updates the iteration count of the MATimerNode by adding the provided count value.
       * @param a_count The count value to add. Default value is 1.
       */
      void update_count();

      /**
       * @brief This function is used to find if a daughter node is already defined with this node name. If this node does not exist, a new daughter MATimerNode is added.
       * @param[in] name name of the desired node
       * @return the MATimerNode desired
       */
      MATimerNode* find(const std::string name);

      // printer functions

      /**
       * @brief Displays a motif several times.
       * @param[in] begin column number where the motif starts.
       * @param[in] end column number where the motif finishs.
       * @param[in] motif the replicated motif, this motif should have a length equal to 1.
       */
      void print_replicate(int a_begin, int a_end, std::string motif);

      /**
       * @brief Displays a blank character.
       */
      void space();

      /**
       * @brief Displays a "|".
       */
      void column();

      /**
       * @brief Displays a return line.
       */
      void end_line();

      /**
       * @brief Displays the banner/header.
       * @param[in] shift number of blank character displayed
       */
      void print_banner(int shift);

      /**
       * @brief Displays the header.
       * @param[in] shift number of blank character displayed
       */
      void print_ending(int shift);

      /**
       * @brief Gets of the duration member
       * @return pointer of the duration member of a MATimerNode
       */
      duration* get_ptr_duration();

      /**
       * @brief Displays the runtime.
       * @param[in] shift number of blank character displayed
       * @param[in] runtime duration value
       */
      void print(int shift, double runtime);

      /**
       * @brief Displays the local runtime.
       * @param[in] shift number of blank character displayed
       * @param[in] runtime local duration value
       */
      void print_local(int shift, double runtime);

      // accessors

      /**
       * @brief Retruns the MATimerNode name
       * @return name
       */
      std::string get_name();

      /**
       * @brief Retruns the MATimerNode iteration number
       * @return the iteration number
       */
      int get_iteration();

      /**
       * @brief Retruns the MATimerNode level
       * @return level
       */
      int get_level();

      /**
       * @brief Retruns a vector of daughter MATimerNode pointers
       * @return daughter nodes
       */
      std::vector<MATimerNode*>& get_daughter();

      /**
       * @brief Retruns the mother MATimerNode pointer
       * @return mother pointer
       */
      MATimerNode* get_mother();

      /**
       * @brief Retruns the duration
       * @return duration value
       */
      double get_duration();

      void inc_mpi() {m_nb_mpi++;} 

      /** @brief This function displays information about one timer node */
      void debug_info();

      private:

      /** @brief name of the measured section */
      std::string m_name;
      /** @brief number of time the measured section is called */
      int m_iteration;
      /** @brief depth of this MATimerNode */
      int m_level;
      /** @brief bunch of daughter MATimerNode */
      std::vector<MATimerNode*> m_daughter;
      /** @brief pointer on the mother MATimerNode, this value is nullptr if it's the root MATimerNode */
      MATimerNode* m_mother = nullptr;
      /** @brief duration time */
      duration m_duration;
      /** @brief this member is used when MATimerNode trees are unbalanced */
      int m_nb_mpi;
    };

    /*
     * @brief Gets the MATimerNode corresponding of an Enum value
     * @return MATimerNode
     * @see enumTimer
     */
    template<enumTimer T>
      MATimerNode*& get_MATimer_node();

    /**
     * @brief Template specialization for getting the MATimerNode pointer for enumTimer::CURRENT.
     * @return The MATimerNode pointer for enumTimer::CURRENT.
     */
    template<>
      MATimerNode*& get_MATimer_node<enumTimer::CURRENT>();

    /**
     * @brief Template specialization for getting the MATimerNode pointer for enumTimer::ROOT.
     * @return The MATimerNode pointer for enumTimer::ROOT.
     */
    template<>
      MATimerNode*& get_MATimer_node<enumTimer::ROOT>();

    /**
     * @brief Debugs the MATimerNode based on the provided enumTimer value.
     * @tparam T The enumTimer value (should be enumTimer::ROOT or enumTimer::CURRENT).
     */
    template<enumTimer T>
      void debug_MATimer_node()
      {
	static_assert(T == enumTimer::ROOT || T == enumTimer::CURRENT, "T should be enumTimer::ROOT or enumTimer::CURRENT");
	auto node = get_MATimer_node<T>();
	std::cout << node << std::endl;
	node->debug_info();
      }
  }
}

#include <MAToolsProfiling/MATimerNode.ixx>
