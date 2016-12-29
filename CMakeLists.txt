cmake_minimum_required( VERSION 2.6 )

project(ljmc)

aux_source_directory(
	"src" 
	ljmc-src
)

add_executable (
	ljmc 
	${ljmc-src}
)

add_definitions( --std=c++11 )
target_link_libraries (
	ljmc
)
