#pragma once

#include <chrono>
#include <MAToolsProfiling/MATimerNode.hxx>

namespace MATools
{
  namespace MATimer
  {
    using duration = std::chrono::duration<double>;
    using high_resolution_clock = std::chrono::high_resolution_clock;
    using MATime_point = std::chrono::time_point<high_resolution_clock>;

    /**
     * @brief Mother class for each Timer type : CPU and GPU. This class is pure virtual.
     */
    template<typename T>  
			  class BaseTimer
    {
      public:

	virtual ~BaseTimer() = default;

	/**
	 * @brief This function sets m_start to the current time point.
	 */
	virtual void start() = 0;

	/**
	 * @brief This function sets m_stop to the current time point.
	 */
	virtual void end() = 0;

      protected:

	/** @brief timer point on the beginning of a time section */
	T m_start;

	/** @brief timer point on the end of a time section */
	T m_stop;
    };

    /**
     * @brief Timer is the BasTimer version with std::chrono time point. A storage is added to accumulate time. 
     */
    class Timer : public BaseTimer<MATime_point>
    {
      public:

	/**
	 * @brief Constructor used with MATimerNode.
	 * @param [in] pointor on a duration that should be stored in a MATimerNode.
	 */ 
	Timer(duration * acc);

	/**
	 * @brief This function sets m_start to the current time point.
	 */
	void start() override;

	/**
	 * @brief This function sets m_stop to the current time point.
	 */
	void end() override;

	/**
	 * @brief Default destructor. The duration is incremented with the time section : m_stop - m_start
	 */
	~Timer(); 

      private:

	/** @brief pointer on a duration item stored in a MATimerNode */
	duration * m_duration; 
    };

    /**
     * @brief BasicTimer is a temporary version of Timer i.e. with no memory.
     */
    class BasicTimer : public BaseTimer<MATime_point>
    {
      public:
	/**
	 * @brief This function sets m_start to the current time point.
	 */
	void start() override;

	/**
	 * @brief This function sets m_stop to the current time point.
	 */
	void end() override;

	/**
	 * @brief return the duration time defined between the time point m_start and m_stop.
	 * @return a duration time in seconds
	 */
	double get_duration();
    };  

    class HybridTimer : public BasicTimer
    {
      public:
	/**
	 * @brief This function sets m_start to the current time point.
	 */
	void start_time_section();

	/**
	 * @brief This function sets m_stop to the current time point.
	 */
	void end_time_section();

	/**
	 * @brief Constructor used with MATimerNode.
	 * @param [in] pointor on a duration that should be stored in a MATimerNode.
	 */ 
	HybridTimer(const std::string a_name);

	/** @brief pointer on a duration item stored in a MATimerNode */
	MATimerNode * m_node; 
    };

    /**
     * @brief Gets an unique Timer for an enumTimer (ROOT)
     * @return a Timer
     */
    template<enumTimer T>
      Timer*& get_timer()
      {
	static Timer* __timer;
	return __timer;
      }

    /**
     * @brief Creates the root MATimerNode tree. This function should be called one time.
     */
    template<enumTimer T>
      void start_global_timer()
      {
	assert(T == enumTimer::ROOT);
	auto& timer = get_timer<T>(); 
	auto& root_ptr = MATools::MATimer::get_MATimer_node<T>(); 
	assert(root_ptr != nullptr);
	timer = new Timer(root_ptr->get_ptr_duration());
	timer->start(); // reset start
	assert(timer != nullptr);
      }

    /**
     * @brief stop the timer corresponding to the root MATimerNode tree. This function should be called one time.
     */
    template<enumTimer T>
      void end_global_timer()
      {
	assert(T == enumTimer::ROOT);
	auto timer = get_timer<T>();
	assert(timer != nullptr);
	timer->end();
      }
  }
}

#include <MAToolsProfiling/Timer.ixx>
