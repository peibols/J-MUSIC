set (SOURCES
    u_derivative.cpp
    dissipative.cpp
    main.cpp
    cell.cpp
    init.cpp
    reconst.cpp
    minmod.cpp
    music.cpp
    cornelius.cpp
    hydro_source.cpp
    pretty_ostream.cpp
    freeze.cpp
    grid_info.cpp
    grid.cpp
    util.cpp
    read_in_parameters.cpp
    freeze_pseudo.cpp
    reso_decay.cpp
    advance.cpp
    eos.cpp
    evolve.cpp
    emoji.cpp
    music_logo.cpp
    Parton.cc
    jets.cc
    jets_embedding.cc
    lund.cc
    tree.cc
    )

if (test)
    add_executable (unittest_grid.e grid.cpp)
    install(TARGETS unittest_grid.e DESTINATION ${CMAKE_HOME_DIRECTORY})
    add_executable (unittest_minmod.e minmod.cpp)
    install(TARGETS unittest_minmod.e DESTINATION ${CMAKE_HOME_DIRECTORY})
else (test)
    add_executable (mpihydro ${SOURCES})
    target_link_libraries (mpihydro ${GSL_LIBRARIES} ${PYTHIA8_LIBRARIES})
    if (APPLE)
        set_target_properties (mpihydro PROPERTIES COMPILE_FLAGS "-DAPPLE")
    endif (APPLE)
    install(TARGETS mpihydro DESTINATION ${CMAKE_HOME_DIRECTORY})
endif (test)
