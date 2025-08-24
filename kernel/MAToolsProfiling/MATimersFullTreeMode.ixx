#include <cstring>
#include <vector>
#include <iostream>
#include <string.h>

#include <MAToolsProfiling/MATimersFullTreeMode.hxx>
#include <MAToolsProfiling/MAToolsMPI.hxx>
#include <MAToolsProfiling/MATimerNode.hxx>
#include <MAToolsProfiling/MAOutputManager.hxx>
#include <MAToolsProfiling/MATimerMPI.hxx>
#include <MAToolsProfiling/MADebug.hxx>

namespace MATools 
{
	namespace MATimer
	{
		namespace FullTreeMode 
		{
			using namespace MATools::MPI;
			
			/**
			 * return the number of bytes in a minimal_info_size. This size doesn't change.
			 * @see class minimal_info
			 * @return int : number of bytes.
			 */
			int minimal_info_size()
			{
				int ret = sizeof(minimal_info);
				return ret;
			}

			/**
			 * default minimal_info constructor*
			 * @see class minimal_info
			 */
			minimal_info::minimal_info() : minimal_info("undefined", 666) {}

			/**
			 * minimal_info constructor*
			 * @see class minimal_info
			 */
			minimal_info::minimal_info(std::string a_name, std::size_t a_s)
			{
				_Pragma("GCC diagnostic push");
				_Pragma("GCC diagnostic ignored \"-Wstringop-truncation\"");
				strncpy(m_name, a_name.c_str(), minimal_info_name_size);
				m_name[minimal_info_name_size-1] = '\0';
				_Pragma("GCC diagnostic pop");
				m_nb_daughter = a_s;
			}

			/**
			 * print information of a minimal_info structure on the master node 
			 * @see class minimal_info
			 * @return void
			 */
			void minimal_info::print()
			{
				using namespace MATools::MAOutput;
				std::string name = m_name;
				printMessage(" (" , name , "," , m_nb_daughter, ")");
			}

			/**
			 * build a vector of minimal info per MPI process to send them (after this routine) on the master process. 
			 * @see class minimal_info
			 * @return a vector of minimal_info
			 */
			std::vector<minimal_info> build_my_tree()
			{
				typedef std::vector<minimal_info> VMI;
				VMI ret;
				MATimerNode* root_timer = MATools::MATimer::get_MATimer_node<ROOT>();

				// this function is called on every MATimerNode.
				auto build_tree = [](MATimerNode* a_ptr, VMI& a_vec)
				{
					auto& daughter = a_ptr->get_daughter();
					std::size_t dsize = daughter.size();
					a_vec.push_back(minimal_info(a_ptr->get_name(), dsize));
				};   

				MATools::MAOutputManager::recursive_call(build_tree, root_timer, ret);

				return ret;
			}


			/**
			 * build a vector of minimal info on the master process. Every trees are packed in a vector of minimal_info and they are gathered on the master node. 
			 * @see class minimal_info
			 * @return void
			 */
			void build_full_tree()
			{
				auto my_info = build_my_tree();

				int info_size = sizeof(minimal_info) * my_info.size();
				int mpi_size = mfem::Mpi::WorldSize();

				std::vector<int> sizes(mpi_size,0);
				std::vector<int> displ(mpi_size, 0);

				MPI_Allgather(&info_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

				int size = 0;
				displ[0] = 0;
				int itD = 1;
				for(auto it : sizes) {
					size += it;
					if(itD <mpi_size);
					displ[itD++]=size;
				}

				std::vector<minimal_info> recv (size/sizeof(minimal_info));  

				MPI_Allgatherv(my_info.data(), info_size, MPI_CHAR, recv.data() , sizes.data() , displ.data(), MPI_CHAR, MPI_COMM_WORLD);

				// convert to int
				for(auto& it : sizes) it /= sizeof(minimal_info);
				transform_to_MATimerMPI(recv, sizes, mpi_size);

				//MATools::MADebug::debug_print<ROOT>();
			}
		};
	};
};
