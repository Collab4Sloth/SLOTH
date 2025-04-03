
/**
 * @file HDF54Sloth.cpp
 * @author cp273896  (clement.plumecocq@cea.fr)
 * @brief Class to define interaction between hdf5 files and sloth
 * @version 0.1
 * @date 2025-03-19
 *
 * Copyright CEA (c) 2025
 *
 */
#include <H5Cpp.h>

#include <boost/multi_array.hpp>
#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>
#pragma once

template <int T>
class HDF54Sloth {
 private:
  /* data */
 public:
  HDF54Sloth();
  void get_data_from_HDF5(const H5std_string& file_name, const std::string& dataset_name,
                          boost::multi_array<double, T>& output_multi_array);
  void get_data_from_HDF5(const H5std_string& file_name, const std::string& dataset_name,
                          std::vector<double>& output_multi_array);
  void write_data_to_HDF5(const H5std_string& file_name, const std::string& dataset_name,
                          const boost::multi_array<double, T>& input_multi_array);
  void write_data_to_HDF5(const H5std_string& file_name, const std::string& dataset_name,
                          const std::vector<double>& input_multi_array);
  ~HDF54Sloth();
};
template <int T>
HDF54Sloth<T>::HDF54Sloth() {}

/**
 * @brief Get the multi-array store in HDF5 dataset
 *
 * @tparam T
 * @param file_name
 * @param dataset_name
 * @param output_multi_array
 */
template <int T>
void HDF54Sloth<T>::get_data_from_HDF5(const H5std_string& file_name,
                                       const std::string& dataset_name,
                                       boost::multi_array<double, T>& output_multi_array) {
  bool is_hdf5 = H5::H5File::isHdf5(file_name);
  H5::H5File file(file_name, H5F_ACC_RDONLY);
  H5::DataSet dataset = file.openDataSet(dataset_name);
  H5::DataSpace dataspace = dataset.getSpace();
  int rank = dataspace.getSimpleExtentNdims();
  std::vector<hsize_t> dims(rank);
  dataspace.getSimpleExtentDims(dims.data(), NULL);
  std::vector<size_t> extents(dims.begin(), dims.end());
  output_multi_array.resize(extents);
  dataset.read(output_multi_array.data(), H5::PredType::NATIVE_DOUBLE);
}

/**
 * @brief Get the vector store in HDF5 dataset
 *
 * @tparam T
 * @param file_name
 * @param dataset_name
 * @param output_multi_array
 */
template <int T>
void HDF54Sloth<T>::get_data_from_HDF5(const H5std_string& file_name,
                                       const std::string& dataset_name,
                                       std::vector<double>& output_multi_array) {
  bool is_hdf5 = H5::H5File::isHdf5(file_name);
  H5::H5File file(file_name, H5F_ACC_RDONLY);
  H5::DataSet dataset = file.openDataSet(dataset_name);
  H5::DataSpace dataspace = dataset.getSpace();
  int rank = dataspace.getSimpleExtentNdims();
  std::vector<hsize_t> dims(rank);
  dataspace.getSimpleExtentDims(dims.data(), NULL);
  size_t total_size = 1;
  for (size_t dim : dims) {
    total_size *= dim;
  }
  output_multi_array.resize(total_size);
  dataset.read(output_multi_array.data(), H5::PredType::NATIVE_DOUBLE);
}

template <int T>
void HDF54Sloth<T>::write_data_to_HDF5(const H5std_string& file_name,
                                       const std::string& dataset_name,
                                       const boost::multi_array<double, T>& input_multi_array) {
  std::ifstream file_check(file_name.c_str());
  bool file_exists = std::filesystem::exists(file_name);

  H5::H5File file(file_name, file_exists ? H5F_ACC_RDWR : H5F_ACC_TRUNC);
  std::vector<hsize_t> dims(T);
  for (int i = 0; i < T; ++i) {
    dims[i] = input_multi_array.shape()[i];
  }
  H5::DataSpace dataspace(T, dims.data());
  H5::DataSet dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);
  dataset.write(input_multi_array.data(), H5::PredType::NATIVE_DOUBLE);
}

template <int T>
void HDF54Sloth<T>::write_data_to_HDF5(const H5std_string& file_name,
                                       const std::string& dataset_name,
                                       const std::vector<double>& input_multi_array) {
  std::ifstream file_check(file_name.c_str());
  bool file_exists = std::filesystem::exists(file_name);

  H5::H5File file(file_name, file_exists ? H5F_ACC_RDWR : H5F_ACC_TRUNC);
  std::vector<hsize_t> dims(1);
  dims[0] = input_multi_array.size();
  H5::DataSpace dataspace(1, dims.data());
  H5::DataSet dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, dataspace);
  dataset.write(input_multi_array.data(), H5::PredType::NATIVE_DOUBLE);
}

template <int T>
HDF54Sloth<T>::~HDF54Sloth() {}
