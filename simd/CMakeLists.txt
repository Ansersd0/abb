file(GLOB_RECURSE fft_srcs fft/*.c)

add_library(simd STATIC
        ${fft_srcs}
        )

target_include_directories(simd PUBLIC
        fft
        )

target_link_libraries(simd PRIVATE m)

target_compile_options(simd PRIVATE -mavx2)

######################## fft test ############################

add_executable(fft_test
        fft/test.c
        )

target_link_libraries(fft_test simd)
target_compile_options(fft_test PRIVATE -mavx2)