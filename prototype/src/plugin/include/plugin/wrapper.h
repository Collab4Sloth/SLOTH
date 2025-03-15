#pragma once

#include <memory>

namespace SlothProto
{
  template<typename T>
  struct Wrapper
  {
    std::shared_ptr<T> data;
    Wrapper() { data = nullptr; }
    bool empty() { return data == nullptr; } 
    void default_definition() { data = std::shared_ptr<T>(new T()); }

    void initialize() {
      if(empty()) { default_definition(); }
      else { std::cout << "warning, this wrapper is already initialized" << std::endl;}
    }
    void alias(std::shared_ptr<T>& object) { data = object; }
    void alias(T* object) { data = std::shared_ptr<T>(object); }
    T& get() { return *data; };
    T* get_ptr() { return data.get(); };
  };
}
