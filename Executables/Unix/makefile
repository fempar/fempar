# Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
#
# This file is part of FEMPAR (Finite Element Multiphysics PARallel library)
#
# FEMPAR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FEMPAR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FEMPAR. If not, see <http://www.gnu.org/licenses/>.
#
# Additional permission under GNU GPL version 3 section 7
#
# If you modify this Program, or any covered work, by linking or combining it 
# with the Intel Math Kernel Library and/or the Watson Sparse Matrix Package 
# and/or the HSL Mathematical Software Library (or a modified version of them), 
# containing parts covered by the terms of their respective licenses, the
# licensors of this Program grant you additional permission to convey the 
# resulting work. 
#
###############################################################################
#=======================================================================
# Rules to build libraries and programs
.PHONY: default all debug release config
.PHONY: clean clean_all clean_release clean_debug

default: debug TAGS

all: debug release

TAGS: ../../Sources/Drivers/Fem/$(FEMPAR_USER)/* ../../Sources/Drivers/Par/$(FEMPAR_USER)/* ../../Sources/Fem/* ../../Sources/Fem/Geometry/* ../../Sources/Fem/Integration/* ../../Sources/Fem/LinearAlgebra/* ../../Sources/Fem/Tools/* ../../Sources/Par/* ../../Sources/Par/Geometry/* ../../Sources/Par/LinearAlgebra/* ../../Sources/Par/Tools/*  ../../Sources/Generic/* 
	@rm -f TAGS
	@etags ../../Sources/Drivers/Fem/$(FEMPAR_USER)/* ../../Sources/Drivers/Par/$(FEMPAR_USER)/* ../../Sources/Fem/* ../../Sources/Fem/Geometry/* ../../Sources/Fem/Integration/* ../../Sources/Fem/LinearAlgebra/* ../../Sources/Fem/Tools/* ../../Sources/Par/* ../../Sources/Par/Geometry/* ../../Sources/Par/LinearAlgebra/* ../../Sources/Par/Tools/* ../../Sources/Generic/*

config:
	@make --no-print-directory -f make.libs config
	@make --no-print-directory -f make.drivers config

newfile:
	@make --no-print-directory -f make.libs newfile
	@make --no-print-directory -f make.drivers newfile

debug:
	@make --no-print-directory -f make.libs debug
	@make --no-print-directory -f make.drivers debug

release:
	@make --no-print-directory -f make.libs release
	@make --no-print-directory -f make.drivers release


clean: clean_debug

clean_debug:
	@make --no-print-directory  -f make.libs clean_debug
	@make --no-print-directory  -f make.drivers clean_debug

clean_release:
	@make --no-print-directory  -f make.libs clean_release
	@make --no-print-directory  -f make.drivers clean_release

clean_all:
	@make --no-print-directory -f make.libs clean_all
	@make --no-print-directory -f make.drivers clean_all

clean_drivers:
	@make --no-print-directory -f make.drivers clean_all

clean_libs:
	@make --no-print-directory -f make.libs clean_all


#=======================================================================
# Serial only entries
#
serial_config:
	@make --no-print-directory -f make.libs config_serial
	@make --no-print-directory -f make.drivers config_serial
serial:
	@make --no-print-directory -f make.libs debug_serial
	@make --no-print-directory -f make.drivers debug_serial

serial_all:
	@make --no-print-directory -f make.libs debug_serial
	@make --no-print-directory -f make.libs release_serial
	@make --no-print-directory -f make.drivers debug_serial
	@make --no-print-directory -f make.drivers release_serial

serial_clean:
	@make --no-print-directory -f make.libs clean_debug_serial
	@make --no-print-directory -f make.drivers clean_debug_serial

serial_clean_all:
	@make --no-print-directory -f make.libs clean_debug_serial
	@make --no-print-directory -f make.libs clean_release_serial
	@make --no-print-directory -f make.drivers clean_debug_serial
	@make --no-print-directory -f make.drivers clean_release_serial
