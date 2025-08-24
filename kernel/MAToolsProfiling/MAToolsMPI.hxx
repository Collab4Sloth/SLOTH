#pragma once

#include <mpi.h>

/**
 * @namespace MATools
 * @brief Namespace containing utility tools for various purposes.
 */
namespace MATools
{
	/**
	 * @namespace MPI
	 * @brief Namespace containing MPI-related utilities.
	 */
	namespace MPI
	{
		/**
		 * @brief Initializes the MPI environment.
		 * 
		 * This function initializes the MPI (Message Passing Interface) environment. The parameters `argc` and `argv` are unused.
		 * @param argc The number of command line arguments (unused).
		 * @param argv The command line arguments (unused).
		 */
		void mpi_initialize([[maybe_unused]] int *argc, [[maybe_unused]] char ***argv);

		/**
		 * @brief Finalizes the MPI environment.
		 * This function finalizes the MPI (Message Passing Interface) environment.
		 */
		void mpi_finalize();

		/**
		 * @brief Checks if MPI is initialized using MPI_Initialized routine.
		 * This function checks if MPI (Message Passing Interface) is initialized.
		 * @return True if MPI is initialized, false otherwise.
		 */
		bool check_mpi_initialized();

		/**
		 * @brief Checks if MPI is finalized using MPI_Finalized routine.
		 * This function checks if MPI (Message Passing Interface) is finalized.
		 * @return True if MPI is finalized, false otherwise.
		 */
		bool check_mpi_finalized();

		/**
		 * @brief Checks if the current process is the master process.
		 * This function checks if the current process is the master process in the MPI environment.
		 * @return True if the current process is the master process, false otherwise.
		 */
		bool is_master();

		/**
		 * @brief Gets the rank of the current process.
		 * This function returns the rank of the current process in the MPI environment.
		 * @return The rank of the current process.
		 */
		int get_rank();

		/**
		 * @brief Gets the size of the MPI environment.
		 * This function returns the size of the MPI environment, which represents the total number of processes.
		 * @return The size of the MPI environment.
		 */
		int get_mpi_size();

		/**
		 * @brief Reduces a value to the maximum value across all processes.
		 * This function reduces the given value to the maximum value across all processes in the MPI environment.
		 * @param value The value to be reduced.
		 * @return The maximum value across all processes.
		 */
		double reduce_max(double);

		/**
		 * @brief Reduces a value to the mean value across all processes.
		 * This function reduces the given value to the mean value across all processes in the MPI environment.
		 * @param value The value to be reduced.
		 * @return The mean value across all processes.
		 */
		double reduce_mean(double);
	}
}

#include <MAToolsProfiling/MAToolsMPI.ixx>
