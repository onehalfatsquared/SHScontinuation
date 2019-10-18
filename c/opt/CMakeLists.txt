set(SOURCES
	continuation.cpp)

add_library(continuation ${SOURCES})
target_link_libraries(continuation potentials)