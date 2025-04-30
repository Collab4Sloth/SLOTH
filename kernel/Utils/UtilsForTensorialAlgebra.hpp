
/**
 * @file UtilsForTensorialAlgebra.hpp
 * @author cp273896 (clement.plumecocq@cea.fr)
 * @brief Class for linear algebra
 * @version 0.1
 * @date 2025-04-17
 *
 * @copyright Copyright CEA (c) 2025
 *
 */

#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#pragma once
/**
 * @brief
 *
 * @tparam T
 */
template <typename T>
class FlattenedTensor {
 private:
  std::size_t dim;
  std::vector<std::size_t> shape;
  bool is_Cstyle = true;
  std::vector<T> vec;

 public:
  FlattenedTensor();
  FlattenedTensor(std::initializer_list<T> _data);
  FlattenedTensor(const FlattenedTensor<T> &original);


  void set_dim(const std::size_t &v) { this->dim = v; }
  std::size_t get_dim() const { return this->dim; }

  void set_shape(const std::vector<std::size_t> &v) { this->shape = v; }
  std::vector<std::size_t> get_shape() const { return this->shape; }

  void set_is_Cstyle(const bool &v) { this->is_Cstyle = v; }
  bool get_is_Cstyle() const { return this->is_Cstyle; }

  std::vector<T> &get() { return vec; }
  const std::vector<T> &get() const { return vec; }
  FlattenedTensor &operator=(const FlattenedTensor &) = default;


  void emplace_back(const T &value) { vec.emplace_back(value); }
  void emplace_back_and_move(const T &value) { vec.emplace_back(std::move(value)); }
  std::size_t size() const { return vec.size(); }
  void resize(const std::size_t &newsize) { return vec.resize(newsize); }
  T &operator[](std::size_t i) { return vec[i]; }
  const T &operator[](std::size_t i) const { return vec[i]; }
  const T *data() const { return vec.data(); }
  T *data() { return vec.data(); }


  std::size_t flattened_index(std::vector<std::size_t> indices) const;
  T evaluate(std::vector<std::size_t> indices) const;
  void apply_scalling(const std::function<double(double)> &f) {
  std::transform(this->vec.begin(), this->vec.end(), this->vec.begin(), f);
  };

  ~FlattenedTensor();
};

/**
 * @brief Construct a new Flattened Tensor< T>:: Flattened Tensor object
 *
 * @tparam T
 */
template <typename T>
FlattenedTensor<T>::FlattenedTensor() {}

/**
 * @brief Construct a new Flattened Tensor< T>:: Flattened Tensor object
 *
 * @tparam T
 * @param _data
 */
template <typename T>
FlattenedTensor<T>::FlattenedTensor(std::initializer_list<T> _data) : vec(_data) {}

/**
 * @brief Construct a new Flattened Tensor< T>:: Flattened Tensor object
 *
 * @tparam T
 * @param original
 */
template <typename T>
FlattenedTensor<T>::FlattenedTensor(const FlattenedTensor<T> &original) {
  *this = original;
}

/**
 * @brief
 *
 * @tparam T
 * @param indices
 * @return std::size_t
 */
template <typename T>
std::size_t FlattenedTensor<T>::flattened_index(std::vector<std::size_t> indices) const {
  std::size_t flat_index = 0;
  std::size_t stride = 1;

  for (std::size_t i = this->dim; i-- > 0;) {
    flat_index += indices[i] * stride;
    stride *= this->shape[i];
  }

  return flat_index;
}

/**
 * @brief
 *
 * @tparam T
 * @param indices
 * @return T
 */
template <typename T>
T FlattenedTensor<T>::evaluate(std::vector<std::size_t> indices) const {
  std::size_t flat_ind = this->flattened_index(indices);
  T val = this->vec[flat_ind];
  return val;
}

/**
 * @brief Destroy the Flattened Tensor< T>:: Flattened Tensor object
 *
 * @tparam T
 */
template <typename T>
FlattenedTensor<T>::~FlattenedTensor() {}
