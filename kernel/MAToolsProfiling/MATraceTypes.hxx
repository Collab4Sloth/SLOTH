#pragma once

#include <chrono>
#include <fstream>
#include <vector>
#include <cassert>
#include <string.h>
#include <thread>
#include <map>

#include <MAToolsProfiling/MATraceColor.hxx>

namespace MATools
{
	namespace MATrace
	{
		using high_resolution_clock = std::chrono::high_resolution_clock;
		using time_point = std::chrono::time_point<high_resolution_clock>;
		struct MATrace_point
		{
			time_point m_time;
			int m_proc;
			MATrace_point();
			time_point& data();
			const time_point data() const;
			void set_proc();
			void set_proc(int a_rank);
			int get_proc() const;
		};

		constexpr int _sstart = 54;
		class vite_event
		{
			public:
				void add(char a_name[64]);
				int operator[] (char a_name[64]);
				void write_items(std::ofstream& a_out);
				void write_colors(std::ofstream& a_out);

				int m_acc = _sstart;
				std::map<std::string,int> m_data;
		};


		/**
		 * @brief The MATrace_section class represents a time section section captured (MPI or OpenMP).
		 * This class provides information such as name, start and end points, and process ID.
		 */
		class MATrace_section
		{
			public:
				/**
				 * @brief Constructor for MATrace_section.
				 * @param a_name The name of the time section captured.
				 * @param a_ref The reference MATrace_point.
				 * @param a_start The start MATrace_point.
				 * @param a_end The end MATrace_point.
				 */
				MATrace_section(char a_name[64], const MATrace_point& a_ref, const MATrace_point& a_start, const MATrace_point& a_end);

				/**
				 * @brief Writes the MATrace_section to the output file stream.
				 * This function writes the MATrace_section to the specified output file stream.
				 * @param a_out The output file stream to write to.
				 */
				void write(std::ofstream& a_out);

				/**
				 * @brief Writes the MATrace_section to the output file stream with VITE event information.
				 * This function writes the MATrace_section, along with the VITE event information, to the specified output file stream.
				 * @param a_out The output file stream to write to.
				 * @param a_ve The VITE event information.
				 */
				void write(std::ofstream& a_out, vite_event& a_ve);

				/**
				 * @brief Sets the process ID for the MATrace_section.
				 * This function sets the process ID for the MATrace_section.
				 * @param a_rank The process ID to set.
				 */
				void set_proc(int const a_rank);

				int m_proc_id; /**< The process ID for the MATrace_section. */
				double m_start; /**< The start time of the MATrace_section. */
				double m_end; /**< The end time of the MATrace_section. */
				char m_name[64]; /**< The name of the MATrace_section. */
		};	


		typedef std::vector<MATrace_section> Trace;
	}
}

#include <MAToolsProfiling/MATraceTypes.ixx>
