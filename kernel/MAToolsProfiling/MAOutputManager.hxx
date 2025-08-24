#pragma once

#include <algorithm>
#include <fstream>
#include <MAToolsProfiling/Column.hxx>
#include <MAToolsProfiling/MATimerNode.hxx>
#include <MAToolsProfiling/MATimerOptional.hxx>
#include <MAToolsProfiling/MATimerInfo.hxx>

namespace MATools
{
	namespace MAOutputManager
	{
		using MATools::MATimer::MATimerNode;

		/**
		 * @brief Build the name for output file or debug file.
		 * @return The constructed name as a string.
		 */
		std::string build_name();

		/**
		 * @brief Build the current MPI name for output file or debug file.
		 * @return The constructed name as a string.
		 */
		std::string build_current_mpi_name();

		/**
		 * @brief Print the timetable to the console.
		 */
		void print_timetable();

		/**
		 * @brief Write the timetable to a file with the default name.
		 */
		void write_file();

		/**
		 * @brief Write the timetable to a file with a custom name.
		 * @param a_name The name of the output file.
		 */
		void write_file(std::string a_name);

		/**
		 * @brief Write the debug information to a file with the default name.
		 */
		void write_debug_file();

		/**
		 * @brief Write the debug information to a file with a custom name.
		 * @param a_name The name of the debug file.
		 */
		void write_debug_file(std::string a_name);

		/**
		 * @brief Get filtered timers based on a specific name.
		 * @param a_name The name to filter the timers.
		 * @return A vector of MATimerInfo containing the filtered timers.
		 */
		std::vector<MATools::MATimer::MATimerInfo> get_filtered_timers(std::string a_name);

		/**
		 * @brief Print filtered timers based on a specific name to the console.
		 * @param a_name The name to filter the timers.
		 */
		void print_filtered_timers(std::string a_name);

		/**
		 * @brief Print the timetable for a specific enumTimer.
		 * @tparam T The enumTimer value (enumTimer::CURRENT or enumTimer::ROOT).
		 */
		template<enumTimer T>
		void print_timetable()
		{
			MATimerNode* root_timer = MATools::MATimer::get_MATimer_node<T>();
			assert(root_timer != nullptr);
			double runtime = root_timer->get_duration();
			runtime = MATools::MPI::reduce_max(runtime); // if MPI, else return runtime

			auto my_print = [](MATimerNode* a_ptr, size_t a_shift, double a_runtime)
			{
				a_ptr->print(a_shift, a_runtime);
			};

			using namespace MATools::MATimer::Optional;
			const bool sorted_by_name = is_full_tree_mode();

			// the sorting is finally disabled
			[[maybe_unused]] auto sort_comp = [sorted_by_name](MATimerNode* a_ptr, MATimerNode* b_ptr)
			{
				if (sorted_by_name)
					return a_ptr->get_name() > b_ptr->get_name();
				else
					return a_ptr->get_duration() > b_ptr->get_duration();
			};

			auto max_length = [](MATimerNode* a_ptr, size_t& a_count, size_t& a_nbElem)
			{
				size_t length = a_ptr->get_level() * 3 + a_ptr->get_name().size();
				a_count = std::max(a_count, length);
				a_nbElem++;
			};
			size_t count(0), nbElem(0);

			recursive_call(max_length, root_timer, count, nbElem);
			count += 6;
			root_timer->print_banner(count);
			recursive_call(my_print, root_timer, count, runtime);
			//recursive_sorted_call(my_print, sort_comp, root_timer, count, runtime);
			root_timer->print_ending(count);
		}

		/**
		 * @brief Recursive helper function for applying a function to each timer node.
		 * @tparam Func The type of the function.
		 * @tparam Args The type of additional arguments.
		 * @param func The function to apply.
		 * @param ptr The current timer node.
		 * @param arg Additional arguments to pass to the function func.
		 */
		template<typename Func, typename... Args>
		void recursive_call(const Func& func, MATimerNode* const ptr, Args&... arg)
		{
			func(ptr, arg...);
			auto& daughters = ptr->get_daughter();
			for (auto& it : daughters)
			{
				assert(it != nullptr);
				recursive_call(func, it, arg...);
			}
		}

		/**
		 * @brief Recursive helper function for applying a function to each timer node in a sorted order.
		 * @tparam Func The type of the function.
		 * @tparam Sort The type of the sorting predicate.
		 * @tparam Args The type of additional arguments.
		 * @param func The function to apply.
		 * @param mySort The sorting predicate.
		 * @param ptr The current timer node.
		 * @param arg Additional arguments to pass to the function.
		 */
		template<typename Func, typename Sort, typename... Args>
		void recursive_sorted_call(const Func& func, const Sort& mySort, MATimerNode* const ptr, Args&... arg)
		{
			func(ptr, arg...);
			auto& daughters = ptr->get_daughter();
			std::sort(daughters.begin(), daughters.end(), mySort);
			for (auto& it : daughters)
			{
				assert(it != nullptr);
				recursive_sorted_call(func, mySort, it, arg...);
			}
		}
	}
}

#include <MAToolsProfiling/MAOutputManager.ixx>
