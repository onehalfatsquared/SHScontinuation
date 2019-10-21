add_executable(testOpt testOpt.cpp)
add_executable(runDescent runDescent.cpp)
add_executable(findMerge findMerge.cpp)

target_link_libraries(testOpt continuation)
target_link_libraries(runDescent continuation)
target_link_libraries(findMerge continuation)
