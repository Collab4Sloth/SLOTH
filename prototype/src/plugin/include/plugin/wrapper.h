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

    void wrap(std::shared_ptr<T>& object) { 
      if( !empty() ) {
        onika::lout << "warning, this wrapper is already defined" << std::endl;
      }
      data = object; 
    }

    void wrap(T* object) { 
      if( !empty() ) {
        onika::lout << "warning, this wrapper is already defined" << std::endl;
      }
      data = std::shared_ptr<T>(object); 
    }

    T& get() { return *data; };
    T* get_ptr() { return data.get(); };
  };
}
