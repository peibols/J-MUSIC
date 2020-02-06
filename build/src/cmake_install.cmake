# Install script for directory: /home/peibols/projects/rrg-jeon-ac/peibols/J-MUSIC/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/peibols/projects/rrg-jeon-ac/peibols/J-MUSIC/mpihydro" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/peibols/projects/rrg-jeon-ac/peibols/J-MUSIC/mpihydro")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/peibols/projects/rrg-jeon-ac/peibols/J-MUSIC/mpihydro"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/peibols/projects/rrg-jeon-ac/peibols/J-MUSIC/mpihydro")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/peibols/projects/rrg-jeon-ac/peibols/J-MUSIC" TYPE EXECUTABLE FILES "/home/peibols/projects/rrg-jeon-ac/peibols/J-MUSIC/build/src/mpihydro")
  if(EXISTS "$ENV{DESTDIR}/home/peibols/projects/rrg-jeon-ac/peibols/J-MUSIC/mpihydro" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/peibols/projects/rrg-jeon-ac/peibols/J-MUSIC/mpihydro")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/cvmfs/soft.computecanada.ca/nix/var/nix/profiles/16.09/bin/strip" "$ENV{DESTDIR}/home/peibols/projects/rrg-jeon-ac/peibols/J-MUSIC/mpihydro")
    endif()
  endif()
endif()

