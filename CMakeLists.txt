cmake_minimum_required(VERSION 3.29)
project(Schnorr_Side_Channel C)

set(CMAKE_C_STANDARD 11)

link_directories(/usr/local/lib)
link_libraries(papi gmp)


add_compile_options(-w -O3)

add_executable(barretM_ref
        main_barret_mult.c
        src/barret.c
        src/barret.h
        src/size.h
)



add_executable(main_our
    main_barret_reduction.c
        src/barret.c
        src/barret.h
        src/size.h
        )

add_executable(main_ref
        main_ref.c
        src/barret.c
        src/barret.h
        src/size.h)