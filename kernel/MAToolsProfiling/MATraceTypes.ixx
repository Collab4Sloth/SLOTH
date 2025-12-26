#include <MAToolsProfiling/MATrace.hxx>
#include <MAToolsProfiling/MAToolsMPI.hxx>
#include <MAToolsProfiling/MATraceColor.hxx>
#include <random>

namespace MATools
{
	namespace MATrace
	{
		using high_resolution_clock = std::chrono::high_resolution_clock;
		using time_point = std::chrono::time_point<high_resolution_clock>;

		MATrace_point::MATrace_point()
		{   
			m_time = std::chrono::system_clock::now();
		}

		time_point& MATrace_point::data()
		{
			return m_time;
		}

		const time_point MATrace_point::data() const
		{
			return m_time;
		}

		void MATrace_point::set_proc()
		{
			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			m_proc = rank;
		}

		void MATrace_point::set_proc(int a_rank)
		{
			m_proc = a_rank;
		}

		int MATrace_point::get_proc() const
		{
			return m_proc;
		}

		void vite_event::add(char a_name[64])
		{
			std::string name = a_name;
			auto it = m_data.find(name);
			if(it == m_data.end())
				m_data[name] = m_acc++;
		}

		int vite_event::operator[] (char a_name[64])
		{
			std::string name = a_name;
			return m_data[name];
		}

		void vite_event::write_items(std::ofstream& a_out)
		{
			for(auto& it : m_data)
			{
				a_out 	<< "%EventDef PajeSetState " << it.second << std::endl;
				a_out	<< "% Time date"<< std::endl;  
				a_out	<< "% Type string"<< std::endl;
				a_out	<< "% Container string" << std::endl;
				a_out	<< "% Value string"  << std::endl;
				a_out	<< "%EndEventDef"<< std::endl;
			}
		}

		void vite_event::write_colors(std::ofstream& a_out)
		{
			MATraceRGB idle_color = get_idle_color();
			a_out << "6 idle ST_ThreadState 'idle'  '" << idle_color.r << " " << idle_color.g << " " << idle_color.b << "'" << std::endl;
			for(auto& it : m_data)
			{
				auto idx = it.second;
				auto name = it.first;
				MATraceRGB color = get_default_color(idx - _sstart); // _start is defined in MATrace.hxx and idx >= _start
				a_out << "6 " <<name<<" ST_ThreadState '"<< name <<"'  '" << color.r << " " << color.g << " " << color.b << "'" << std::endl;
			}
		}

		/**
		 * @brief Constructor for MATrace_section.
		 * @param a_name The name of the time section captured.
		 * @param a_ref The reference MATrace_point.
		 * @param a_start The start MATrace_point.
		 * @param a_end The end MATrace_point.
		 */
		MATrace_section::MATrace_section(char a_name[64], const MATrace_point& a_ref, const MATrace_point& a_start, const MATrace_point& a_end)
		{
			std::chrono::duration<double> duration ;

			// compute start
			duration = a_start.data() - a_ref.data();
			m_start = duration.count();
			assert(m_start >= 0);

			// compute end
			duration = a_end.data() - a_ref.data();
			m_end = duration.count();
			assert(m_end > m_start);

			// set proc id and name
			m_proc_id = a_ref.get_proc();
			strncpy(m_name, a_name, 63);
			m_name[63]='\0'; // warning issue without this line
		}

		/**
		 * @brief Writes the MATrace_section to the output file stream.
		 * This function writes the MATrace_section to the specified output file stream.
		 * @param a_out The output file stream to write to.
		 */
		void MATrace_section::write(std::ofstream& a_out)
		{
			a_out << m_name   << " "
				<< m_proc_id    << " "
				<< m_start      << " "
				<< m_end        << " "
				<< std::endl;
		}

		/**
		 * @brief Writes the MATrace_section to the output file stream with VITE event information.
		 * This function writes the MATrace_section, along with the VITE event information, to the specified output file stream.
		 * @param a_out The output file stream to write to.
		 * @param a_ve The VITE event information.
		 */
		void MATrace_section::write(std::ofstream& a_out, vite_event& a_ve)
		{
			a_out   << a_ve[m_name] << " "
				<< m_start    << " "
				<< "ST_ThreadState "
				<< "C_Thread" << m_proc_id << " "<< m_name
				<< std::endl;
			a_out   << "10 "
				<< m_end    << " "
				<< "ST_ThreadState "
				<< "C_Thread" << m_proc_id << " idle"
				<< std::endl;
		}

		/**
		 * @brief Sets the process ID for the MATrace_section.
		 * This function sets the process ID for the MATrace_section.
		 * @param a_rank The process ID to set.
		 */
		void MATrace_section::set_proc(const int a_rank)
		{
			m_proc_id = a_rank;
		}
	}
};
