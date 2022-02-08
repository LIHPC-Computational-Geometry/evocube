
project(gco)

include_directories(include) # TODO improve
add_library(gco include/block.h include/energy.h include/GCoptimization.h include/GCoptimization.cpp include/graph.cpp include/graph.h include/LinkedBlockList.cpp include/LinkedBlockList.h include/maxflow.cpp)