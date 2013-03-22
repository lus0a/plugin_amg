################################################################################
# AMG
#
# The algebraic multigrid plugin.
# Ruge-Stüben AMG and Filtering AMG
################################################################################

cmake_minimum_required(VERSION 2.6)

project(AMG)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


# sources
set(srcAMGBridge 	bridge/amg_bridge.cpp
					bridge/famg_bridge.cpp
					bridge/rsamg_bridge.cpp)

if(PARALLEL)
	set(srcAMG src/rsamg/rsamg_parallel_coarsening.cpp
	            src/amg_debug.cpp)
endif(PARALLEL)
	
	set(srcAMG ${srcAMG}
		src/rsamg/rsamg_coarsening.cpp	
		src/agglomeration.cpp
		)								



if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/amg ${srcAMGBridge} ${srcAMG})
else(buildEmbeddedPlugins)
	add_library(amg SHARED ${srcAMGBridge} ${srcAMG})
	target_link_libraries (amg ug4)
endif(buildEmbeddedPlugins)
