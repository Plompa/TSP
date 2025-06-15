cmake_minimum_required(VERSION 3.29)
project(TSP)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS -Wall)
set(CMAKE_CXX_FLAGS -Wextra)
set(CMAKE_CXX_FLAGS -Wpedantic)
set(CMAKE_CXX_FLAGS -Werror)


add_executable(TSP main.cpp
        TSP_reader.cpp
        TSP_reader.h
        TSP_enumeration.cpp
        TSP_enumeration.h
        Graph.cpp
        Graph.h
        TSP_writer.cpp
        TSP_writer.h
        Einstiegsaufgabe.cpp
        Einstiegsaufgabe.h
)
