cmake_minimum_required(VERSION 3.29)
project(Schnorr_Side_Channel C)

set(CMAKE_C_STANDARD 11)

link_directories(/usr/local/lib)
link_libraries(papi gmp)

add_executable(Schnorr_Side_Channel
    main.c
        src/barret.c
        src/barret.h
        src/size.h)
