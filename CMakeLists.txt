cmake_minimum_required(VERSION 3.29)
project(Schnorr_Side_Channel C)

set(CMAKE_C_STANDARD 11)

include_directories(src src/barret_reduction src/utils)

link_directories(/usr/local/lib)
link_libraries(papi gmp)



add_compile_options(-w -O3 -mavx2)

#add_executable(barretM_ref
#        main_barret_mult.c
#        src/barret.c
#        src/barret.h
#        src/size.h
#        src/misc/utils.h
#        src/misc/utils.c
#        src/utils/debug.c
#        src/utils/debug.h
#        src/barret_reduction/barret_reduction.c
#        src/barret_reduction/barret_reduction.h
#)


add_executable(barret_reduction_ref
        src/barret_lib.h
        src/utils/mod.c
        src/utils/debug.c
        src/barret_reduction/barret_reduction.c
        test/main_barret_reduction.c
)


add_executable(barret_multiplication_ref
        src/barret_lib.h
        src/utils/mod.c
        src/utils/debug.c
        src/barret_multiplication/barret_multiplication.c
        test/main_barret_multiplication.c
)


#
#add_executable(main_our
#    main_barret_reduction.c
#        src/barret.c
#        src/barret.h
#        src/size.h
#        )
#
#add_executable(main_ref
#        main_ref.c
#        src/barret.c
#        src/barret.h
#        src/size.h)