# Copyright (c) 2012-2013:  G-CSC, Goethe University Frankfurt
# Author: Martin Rupp
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

################################################################################
# AMG
#
# The algebraic multigrid plugin.
# Ruge-Stueben AMG and Filtering AMG
################################################################################

cmake_minimum_required(VERSION 2.6)

project(AMG)

# include the definitions and dependencies for ug-plugins.
include("../../../cmake/ug_plugin_includes.cmake")


# sources
set(srcAMGBridge 	bridge/amg_bridge.cpp
					bridge/famg_bridge.cpp
					bridge/semiamg_bridge.cpp
					bridge/rsamg_bridge.cpp)

if(PARALLEL)
	set(srcAMG src/rsamg/rsamg_parallel_coarsening.cpp
	            )
endif(PARALLEL)
	
	set(srcAMG ${srcAMG}
		src/rsamg/rsamg_coarsening.cpp	
		src/agglomeration.cpp
		src/famg/debug_ids.cpp	
		src/amg_debug.cpp
		src/semiamg/tensor_positions.cc	
		)								



if(buildEmbeddedPlugins)
	# add the sources to ug4's sources
	EXPORTSOURCES(plugins/experimental/amg ${srcAMGBridge} ${srcAMG})
else(buildEmbeddedPlugins)
	add_library(amg SHARED ${srcAMGBridge} ${srcAMG})
	target_link_libraries (amg ug4)
endif(buildEmbeddedPlugins)
