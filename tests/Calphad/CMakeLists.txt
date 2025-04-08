add_subdirectory(AnalyticalIdealSolution)
# LIBTORCH API ??
if(SLOTH_USE_LIBTORCH)
    add_subdirectory(NeuralNetwork)
endif(SLOTH_USE_LIBTORCH)