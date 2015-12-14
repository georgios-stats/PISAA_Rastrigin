#!/bin/bash

# --------------------------------------------------------------------------------
# 
# Copyrigtht 2014 Georgios Karagiannis
# 
# This file is part of PISAA_Rastrigin.
# 
# PISAA_Rastrigin is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 2 of the License.
# 
# PISAA_Rastrigin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PISAA_Rastrigin.  If not, see <http://www.gnu.org/licenses/>.
# 
# --------------------------------------------------------------------------------

# Georgios Karagiannis
#
# Postdoctoral research associate
# Department of Mathematics, Purdue University
# 150 N. University Street
# West Lafayette, IN 47907-2067, USA
#
# Telephone: +1 (765) 496-1007
#
# Email: gkaragia@purdue.edu
#
# Contact email: georgios.stats@gmail.com

CC=icc
CFLAGS=-O2
LDFLAGS=
CPPFLAGS=-D__ROTATE__=0

FUN=cost_rastrigin.c

SOURCES=pisaa.c \
	Crossover_operations.c \
	Mutation_operations.c \
	HitAndRun_update.c \
	Self_adjastment_prosedure.c \
	$(FUN) \
	rotationmatrix_salomon.c \
	uniformdirectionrng.c \
	normalrng.c \
	uniformrng.c \
	nrutil.c 
	
OBJECTS=$(SOURCES:.c=.o)

EXECUTABLE=exe

# ACTIONS

build: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# CLEAR

clean:
	rm -rf *o exe

# DETAILS

details:
	@echo CC       : $(CC)
	@echo CFLAGS   : $(CFLAGS)
	@echo LDFLAGS  : $(LDFLAGS)
	@echo FUN      : $(FUN)
	@echo CPPFLAGS : $(CPPFLAGS)

# RUN

run:
	./exe

test_run:
	./exe -ID 1 \
			-Ndim 10 \
			-Niter 10000 \
			-Npop 30 \
			-Nsam 0 \
			-Gwarm 100 \
			-Ghigh 1.0 \
			-Gpow 0.55 \
			-Hlow -0.01 \
			-Hhigh 800.0 \
			-Hsize $(shell echo "scale=0; ( 800.0 - (-0.01))/0.05+1;" | bc) \
			-Hzeta 0.0 \
			-Hconst 1.0 \
			-Twarm 1 \
			-Thigh 5.0 \
			-Tlow 0.01 \
			-Tpow 0.5 \
			-Tini 100.0 \
			-Tref 0.00001 \
			-Sini 1.0 \
			-SMO0 2.0 \
			-SMO1 2.0 \
			-SMO2 2.0 \
			-SMO3 0.05 \
			-SCO1 0.1 \
			-SCO2 0.5 \
			-SCO3 0.5 \
			-Sref 0.005


	
