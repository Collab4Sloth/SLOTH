#pragma once


namespace MATools
{
	namespace MATrace
	{
    struct MATraceRGB{double r; double g; double b; };

		constexpr int default_color_size()
		{
			constexpr int ret = 12;
			return ret;
		}

		MATraceRGB get_idle_color()
		{
			MATraceRGB ret = {0.5,0.5,0.5};
			return ret;
		}

		MATraceRGB get_default_color(int color_id)
		{
			auto size = default_color_size();
			const int i = color_id % size;
			switch(i)
			{
				case 0 : return {1, 0, 0};
				case 1 : return {0, 0, 1};
				case 2 : return {0, 1, 0};
				case 3 : return {0.5, 0, 0};
				case 4 : return {1, 1, 0};
				case 5 : return {0.5, 0.5, 0};
				case 6 : return {0, 0.5, 0};
				case 7 : return {0, 1, 1};
				case 8 : return {0, 0.5, 0.5};
				case 9 : return {0, 0, 0.5};
				case 10 : return {1, 0, 1};
				case 11 : return {0.5, 0, 0.5};
				default : return {0,0,0};
			}
		}
	};
};
