file(GLOB_RECURSE transform_src CONFIGURE_DEPENDS
        "fft/*.c"
        "ntt/*.c"
        "fftw/*.c"
        )

add_library(transform STATIC ${transform_src}
        fft/fft1.c
        fft/fft1.h
        fft/common.h
        pft/pft.c
        )

target_include_directories(transform INTERFACE fft ntt fftw)

# set the path to the FFTW library
if (WIN32)
    set(FFTW_PATH "D:/FFTW/fftw-3.3.5-dll64/libfftw3-3.lib")
elseif (UNIX) # Linux
    set(FFTW_PATH "/usr/local/lib/libfftw3.a")
endif ()


# add the FFTW library
include_directories("D:/FFTW/fftw-3.3.5-dll64")

target_link_libraries(transform ${FFTW_PATH})

if (WIN32)
    #    include_directories("C:/Program Files (x86)/Windows Kits/10/Include/10.0.22000.0/um)
elseif (UNIX) # Linux
    target_link_options(transform PRIVATE -lm)
endif ()

############### ntt_test ################
add_executable(ntt_test ntt/ntt_test.c)

target_link_libraries(ntt_test transform)

############### fftw_test ###############

add_executable(fftw_test fftw/fftw_test.c)

target_link_libraries(fftw_test transform)


############### fft_test ###############

add_executable(fft_test fft/fft_test.c)

target_link_libraries(fft_test transform)

############### pft_test ###############

add_executable(pft_test pft/pft_test.c)



target_link_libraries(pft_test transform)
