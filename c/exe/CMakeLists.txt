add_executable(testOpt testOpt.cpp)
add_executable(runDescent runDescent.cpp)

target_link_libraries(testOpt continuation)
target_link_libraries(runDescent continuation)
