file(GLOB array_sources array/*.c)
file(GLOB stack_sources stack/*.c)
file(GLOB queue_sources queue/*.c)
file(GLOB linkList_sources linkList/*.c)
file(GLOB hash_sources hash/*.c)

add_library(dataStructure STATIC
        ${array_sources}
        ${stack_sources}
        ${queue_sources}
        ${linkList_sources}
        ${hash_sources}
        )

target_include_directories(dataStructure PUBLIC array stack queue linkList hash)

########### stack test ############
add_executable(stack_test stack/stack_test.c)
target_link_libraries(stack_test dataStructure)

########### queue test ############
add_executable(queue_test queue/queue_test.c)
target_link_libraries(queue_test dataStructure)

########### linkList test ############
add_executable(linkList_test linkList/linkList_test.c)
target_link_libraries(linkList_test dataStructure)

########### hash test ############
add_executable(hash_test hash/hash_test.c)
target_link_libraries(hash_test dataStructure)