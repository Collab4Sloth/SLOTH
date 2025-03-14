#pragma once


namespace SlothProto
{
  template<typename P>
  struct Wrapper
  {
    Wrapper() { data = new P();}
    P* data;
    P& get() {return *data;};
  };
}
