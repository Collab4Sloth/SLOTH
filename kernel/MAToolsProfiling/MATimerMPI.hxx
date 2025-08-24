#pragma once

#include <MAToolsProfiling/MATimers.hxx>
#include <MAToolsProfiling/MATimerNode.hxx>
#include <MAToolsProfiling/MATimersFullTreeMode.hxx>
#include <MAToolsProfiling/EnumTimer.hxx>

namespace MATools
{
	namespace MATimer
	{
		namespace FullTreeMode
		{
			using namespace MATools::MATimer;


			void rec_call(MATimerNode* a_node, minimal_info* a_ptr)
			{
				int elem = a_ptr->m_nb_daughter;
				a_node->inc_mpi();
				
				for(int it = 0 ; it < elem ; it++)
				{
					auto ptr = a_ptr + it + 1; 
					std::string name = ptr->m_name;
					auto node = a_node->find(name);
					rec_call(node, ptr);
				}
			}

			void transform_to_MATimerMPI(std::vector<minimal_info>& a_in, std::vector<int> a_sizes, int a_mpi_size)
			{
				MATimerNode*& root = get_MATimer_node<enumTimer::ROOT>();
				int acc = 0;
				int rank = MPI::get_rank();
				for(int mpi = 0 ; mpi < a_mpi_size; mpi++)
				{
					if(a_sizes[mpi] == 0) continue;
					auto local_root = a_in.data() + acc;
					rec_call(root,local_root);
					acc += a_sizes[mpi];
				}
			}
		};
	};
};
