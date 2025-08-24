/**
 * @file
 * @brief This file contains the FullTreeMode namespace and related classes and functions.
 */

#pragma once

#include <cstring>
#include <string>
#include <vector>
#include <cstdint>

namespace MATools 
{
	namespace MATimer
	{
		/**
		 * @namespace FullTreeMode
		 * @biref The FullTreeMode namespace provides classes and functions related to Full Tree Mode.
		 */
		namespace FullTreeMode
		{

			constexpr size_t minimal_info_name_size = 64; /**< The size of the minimal_info name. */

			/**
			 * @class minimal_info
			 * @brief Represents minimal information about a node with the Full Tree Mode.
			 * The minimal_info class represents minimal information about a node with the Full Tree Mode,
			 * including its name and the number of daughter nodes.
			 */
			class minimal_info
			{
				public:
					minimal_info(); /**< Default constructor. */
					
					/**
					 * @brief Constructs a minimal_info object with specified values.
					 * @param a_name The name of the node.
					 * @param a_nb_daughter The number of daughter nodes.
					 */
					minimal_info(std::string a_name, std::size_t a_nb_daughter);
					
					/**
					 * @brief Prints the minimal_info object.
					 */
					void print();
					
					char m_name[minimal_info_name_size]; /**< The name of the node. */
					std::uint64_t m_nb_daughter; /**< The number of daughter nodes. */
			};

			/**
			 * @brief Gets the size of the minimal_info object.
			 * @return The size of the minimal_info object.
			 */
			int minimal_info_size();
			
			/**
			 * @brief Builds a vector of minimal_info objects representing the tree structure.
			 * @return The vector of minimal_info objects representing the tree structure.
			 */
			std::vector<minimal_info> build_my_tree();
			
			/**
			 * @brief Builds the full tree structure.
			 */
			void build_full_tree();
		};
	};
};

#include <MAToolsProfiling/MATimersFullTreeMode.ixx>
