set(SOURCES
	euDist.cpp
	stickyFn.cpp
	morse.cpp
	lennardJones.cpp)

add_library(potentials ${SOURCES})
target_link_libraries(potentials dlib::dlib)

