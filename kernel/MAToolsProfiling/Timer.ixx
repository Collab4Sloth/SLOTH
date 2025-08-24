#include <MAToolsProfiling/Timer.hxx>

namespace MATools
{
	namespace MATimer
	{
		/**
		 * @brief Constructor used with MATimerNode.
		 * @param [in] pointor on a duration that should be stored in a MATimerNode.
		 */
		Timer::Timer(duration * acc) 
		{
			m_duration = acc;
			start();
		}

		/**
		 * @brief This function sets m_start to the current time point.
		 */
		void Timer::start()
		{
			m_start = high_resolution_clock::now();
		}

		/**
		 * @brief This function sets m_stop to the current time point.
		 */
		void Timer::end()
		{
			assert(m_duration != nullptr && "duration has to be initialised");
			m_stop = high_resolution_clock::now();
			*m_duration += m_stop - m_start;
			assert(m_duration->count() >= 0);
		}

		/**
		 * @brief Default destructor. The duration is incremented with the time section : m_stop - m_start
		 */
		Timer::~Timer() 
		{
			end();
			auto& current_timer = MATools::MATimer::get_MATimer_node<CURRENT>();
			current_timer = current_timer->get_mother();
		}

		/**
		 * @brief This function sets m_start to the current time point.
		 */
		void BasicTimer::start()
		{
			m_start = high_resolution_clock::now();
		}

		/**
		 * @brief This function sets m_stop to the current time point.
		 */
		void BasicTimer::end()
		{
			m_stop = high_resolution_clock::now();
		}

		/**
		 * @brief return the duration time defined between the time point m_start and m_stop.
		 * @return a duration time in seconds
		 */
		double BasicTimer::get_duration()
		{
			using duration = std::chrono::duration<double>;
			duration measure = m_stop - m_start;
			double ret = measure.count();
			return ret;
		}

		HybridTimer::HybridTimer(const std::string a_name)
		{
			using namespace MATools::MATimer;
			auto* ptr_matimer_current_node = MATools::MATimer::get_MATimer_node<CURRENT>();
			assert(ptr_matimer_current_node != nullptr && "do not use an undefined MATimerNode");
			auto* ptr_matimer_node = ptr_matimer_current_node->find(a_name);
			assert(ptr_matimer_node != nullptr && "do not use an undefined MATimerNode");
			this->m_node = ptr_matimer_node;
		}

		void HybridTimer::start_time_section()
		{
			this->start();
		}

		void HybridTimer::end_time_section()
		{
			this->end();
			this->m_node->update_count();
			this->m_node->get_ptr_duration()[0] += this->m_stop - this->m_start;
		}
	};
};
