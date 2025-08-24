#include <vector>
#include <cassert>
#include <fstream>
#include <array>
#include <tuple>
#include <cstring>
#include <MAToolsProfiling/MAOutput.hxx>
#include <MAToolsProfiling/MAMemory.hxx>
#include <MAToolsProfiling/MAToolsMPI.hxx>

namespace MATools
{
  namespace MAMemory
  {

    /**
     * This function creates a rusage variable
     * @return rusage data type that get memory informations
     */
    rusage make_memory_checkpoint()
    {
      rusage obj;
      int who = 0;
      [[maybe_unused]] auto res = getrusage(who, &obj);
      assert((res = -1) && "error: getrusage has failed");
      return obj;
    };

    /**
     * Default constructor that calls default std::vector constructor 
     */
    MAFootprint::MAFootprint() {}

    /**
     * This function create a memory checkpoint, i.e. a rusage, and insert it in the data storage 
     */
    void MAFootprint::add_memory_checkpoint()
    {
      auto obj = make_memory_checkpoint();
      this->push_back(obj);
    }


    /**
     * Reduce is called to get the total memory footprint among MPI processes for every memory checkpoints
     * Data are gathered/reduced on the master mpi process 0.
     * @return a vector with the total memory footprint by memory checkpoints  
     */
    std::vector<long> MAFootprint::reduce()
    {
      // one value per memory checkpoint
      std::vector<long> ret;
      int nb_points = this->size();
      auto m_data = this->data();
      ret.resize(nb_points);

      // get total memory footprint for every memory checkpoints
      for(int id = 0 ; id < nb_points ; id++)
      {
	ret[id] = 0;
	MPI_Reduce(&(m_data[id].ru_maxrss), &(ret[id]), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
      }
      return ret;
    }

    std::vector<std::tuple<long, long, long, double>> MAFootprint::get_sum_min_max_mean()
    {
      using namespace MATools::MPI;
      std::vector<std::tuple<long, long, long, double>> res;
      int nb_points = this->size();
      auto m_data = this->data();
      res.resize(nb_points);

      const int mpi_size = get_mpi_size();
      std::vector<long> buffer;
      //if(is_master()) buffer.resize(mpi_size);
      buffer.resize(mpi_size);

      // get total memory footprint for every memory checkpoints
      for(int id = 0 ; id < nb_points ; id++)
      {
	std::tuple<long, long, long, double> item;
	MPI_Gather(&(m_data[id].ru_maxrss), 1, MPI_LONG, buffer.data(), 1, MPI_LONG, 0, MPI_COMM_WORLD);
	long _sum(0);
	long _min(buffer[0]);
	long _max(buffer[0]);
	double _mean = 0.0;
	for(int mpi = 0; mpi < mpi_size ; mpi++)
	{
	  _sum += buffer[mpi];
	  _min = std::min(_min, buffer[mpi]);
	  _max = std::max(_max, buffer[mpi]);
	}
	_mean = ((double)(_sum)) / ((double)(mpi_size));
	item = std::make_tuple(_sum, _min, _max, _mean);
	res[id] = std::move(item);
      }

      return res;
    }


    /*
     * @brief Get the maximum memory usage value (rusage::.ru_maxrss) for the point number i
     * @param a_idx a_idx is the index of the memory point
     * @return return the ru_maxrss of the rusage structure
     */
    long MAFootprint::get_usage(const int a_idx)
    {
      const long ret = this->data()[a_idx].ru_maxrss;
      return ret;
    }

    /*
     * The memory footprint is printed for every memory checkpoints
     * @param f is a mafootprint object that contains memory checkpoints
     * @see mafootprint
     */
    void print_checkpoints(MAFootprint& a_f)
    {
      using namespace MATools::MPI;
      // This function extracts the maxmimal data size (footprint) and do a reduction if the MPI feature is activated
      // For every memory checkpoints
      auto obj = a_f.reduce();
      if(is_master())
      {
	std::cout << " List (maximum resident size): ";
	for(auto it : obj)
	{
	  std::cout << " " << it;
	}
	std::cout << std::endl;
      }
    };


    /*
     * The memory footprint is written for every memory checkpoints
     * @param f is a mafootprint object that contains memory checkpoints
     * @see mafootprint
     */
    void write_memory_checkpoints(MAFootprint& a_f, std::string file_name)
    {
      using namespace MATools::MPI;
      // This function extracts the maxmimal data size (footprint) and do a reduction if the MPI feature is activated
      // For every memory checkpoints
      auto obj = a_f.reduce();
      if(is_master())
      {
	std::ofstream out (file_name, std::ios::out);
	std::cout << " List (maximum resident size): ";
	for(size_t i = 0 ; i < obj.size() ; i++)
	{
	  out << i << " " << obj[i] << std::endl;
	}
      }
    };

    /*
     * The memory footprint is written for every memory checkpoints
     * @param f is a mafootprint object that contains memory checkpoints
     * @param a_labels is a list of labels
     * @see mafootprint
     */
    void write_memory_checkpoints(MAFootprint& a_f, std::vector<std::string>& a_labels,  std::string file_name)
    {
      using namespace MATools::MPI;

      // rename file
      // replace * by mpi_size
      std::string pattern =  "*";
      auto ppos = file_name.find(pattern);
      if( ppos == std::string::npos )
      {
	MATools::MAOutput::printMessage("MATools_Error: output filename should contains the pattern *");
	std::abort();
      }

      const int mpi_size = get_mpi_size();
      file_name = file_name.replace(ppos, 1, std::to_string(mpi_size));

      MATools::MAOutput::printMessage("Writting" , file_name );

      // This function extracts the maxmimal data size (footprint) and do a reduction if the MPI feature is activated
      // For every memory checkpoints
      auto obj = a_f.get_sum_min_max_mean(); 
      if(is_master())
      {
	std::ofstream out (file_name, std::ios::out);
	std::cout << " List (maximum resident size): ";
	for(size_t i = 0 ; i < obj.size() ; i++)
	{
	  auto item = obj[i];
	  long _sum = std::get<0>(item);
	  long _min = std::get<1>(item);
	  long _max = std::get<2>(item);
	  double _mean = std::get<3>(item);
	  out << a_labels[i] << " " << _sum << " " << _min << " " << _max << " " << _mean << std::endl;
	}
      }
    };

    /*
     * The memory footprint is printed where this function is called.
     * @see MAFootprint
     */
    void print_memory_footprint()
    {
      // Create an MAFootprint object that will be destroyed at the end
      MAFootprint f;
      // We fill this function with only one checkpoint
      f.add_memory_checkpoint();
      // This function extracts the maxmimal data size (footprint) and do a reduction if the MPI feature is activated
      auto obj = f.reduce();
      // As we have only one memory checkpoint, we get this element
      auto last = obj.back() * 1e-6; // conversion kb to Gb
      MATools::MAOutput::printMessage(" memory footprint: ", last, " GB");
    };


    template<class IO>
      void io_trace_memory_points_per_mpi(IO& a_stream, MAFootprint& a_mem_points)
      {
	using namespace MATools;
	//			int rank = MPI::get_rank();
	int size = MPI::get_mpi_size();
	int number_of_mem_points = a_mem_points.size();
	std::vector<std::vector<long>> mpi_trace(number_of_mem_points, std::vector<long>(size));
	// gather all data per rank
	for(int i = 0 ; i < number_of_mem_points ; i++)
	{
	  long local = a_mem_points.get_usage(i);
	  MPI_Gather(&local, 1, MPI_LONG, mpi_trace[i].data(), 1, MPI_LONG, 0, MPI_COMM_WORLD);
	}

	int point_id = 0; 
	for(auto& mem_points : mpi_trace)
	{
	  std::string buf = std::to_string(point_id++) + " "; 
	  for(auto& mpi : mem_points)
	  {
	    buf += std::to_string(mpi) + " "; 
	  }
	  a_stream << buf << std::endl;
	}
      }

    /*
     * The memory footprint is printed for every memory checkpoints according to their mpi rank.
     * @param f is a mafootprint object that contains memory checkpoints
     * @see mafootprint
     */
    void print_trace_memory_points_per_mpi(MAFootprint&& a_mem_points)
    {
      io_trace_memory_points_per_mpi(std::cout, a_mem_points);
    }

    /*
     * The memory footprint is written for every memory checkpoints according to their mpi rank.
     * @param f is a mafootprint object that contains memory checkpoints
     * @see mafootprint
     */
    void write_trace_memory_points_per_mpi(MAFootprint&& a_mem_points, std::string a_name)
    {
      std::ofstream ofs;
      ofs.open(a_name, std::ofstream::out);
      io_trace_memory_points_per_mpi(ofs, a_mem_points);
      ofs.close();
    }

    /**
     * @brief This function gets the "common" MAFootprint or create if necessary.
     */
    MAFootprint& get_MAFootprint()
    {
      static MAFootprint ret;
      return ret;
    }
  };
};
