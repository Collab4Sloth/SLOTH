#pragma once

#include <MAToolsProfiling/MATimerNode.hxx>

namespace MATools
{

	namespace MADebug
	{

		using namespace MATools::MAOutputManager;
		using namespace MATools::MATimer;
		using namespace MATools::MPI;

		/**
		 * return print data in a from a MATimerNode
		 * @see class MATimerNode
		 */
		template<enumTimer T>
			void debug_print()
			{
				MATimerNode*& ptr = get_MATimer_node<T>();

				auto my_print = [](MATimerNode* a_ptr, int& shift, double a_runtime)
				{
					a_ptr->print(shift, a_runtime);
				};

				auto sort_comp = [] (MATimerNode* a_ptr, MATimerNode* b_ptr)
				{
					return a_ptr->get_name() > b_ptr->get_name() ;
				};

				int shift=0;
				double d = ptr->get_duration();
				//recursive_sorted_call(my_print, sort_comp, ptr, shift, d);
				recursive_call(my_print, ptr, shift, d);
			}
		
		template<enumTimer T>
			void debug_write()
			{
				MATimerNode*& ptr = get_MATimer_node<T>();

				auto my_print = [](MATimerNode* a_ptr, int& shift, double a_runtime)
				{
					a_ptr->print(shift, a_runtime);
				};

				int shift=0;
				double d = ptr->get_duration();
				//recursive_sorted_call(my_print, sort_comp, ptr, shift, d);
				recursive_call(my_print, ptr, shift, d);
			}
		
		/**
		 * return print data in a from a MATimerNode of the master process
		 * @see class MATimerNode
		 */
		template<enumTimer T>
			void master_debug_print()
			{
				MATimerNode*& ptr = get_MATimer_node<T>();

				auto my_print = [](MATimerNode* a_ptr, int& shift, double a_runtime)
				{
					a_ptr->print_local(shift, a_runtime);
				};
				int shift=0;
				double d = ptr->get_duration();
				if(is_master()) recursive_call(my_print, ptr, shift, d);
			}
	};
};
