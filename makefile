##
## @file		makefile
## @brief	PPDyn makefile.
## @author	Dr. Sayan Adhikari <sayan.adhikari@fys.uio.no>
##		Dr. Rupak Mukherjee <rupakm@princeton.edu>
##

# # OS DETECTION
# BASHFILE	:=
UNAME_S := $(shell uname -s)
# ifeq ($(UNAME_S),Linux)
# 	BASHFILE = bashrc
# endif
# ifeq ($(UNAME_S),Darwin)
# 	BASHFILE = zshrc
# endif



# the virtual environment directory
VENV := venv

# Source directory
SRC := src

# Current directory
# PATHS := $(shell pwd)
# mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
CURRENT_DIR := $(shell pwd)

# Figure directory
FIG := figures
DATA := data

ifeq ($(UNAME_S),Linux)
all: venv
	@echo "Creating virtual environment for running the code"
$(VENV)/bin/activate: requirements.txt
	python -m venv $(VENV)
	./$(VENV)/bin/pip install -r requirements.txt
	mkdir $(FIG) $(DATA) 2> /dev/null
	@echo "alias ppdyn='$(CURRENT_DIR)/./$(VENV)/bin/python3 $(CURRENT_DIR)/$(SRC)/main.py'" >> $${HOME}/.bashrc
	. $${HOME}/.bashrc
else ifeq ($(UNAME_S),Darwin)
all: venv
	@echo "Creating virtual environment for running the code"
$(VENV)/bin/activate: requirements.txt
	python3 -m venv $(VENV)
	./$(VENV)/bin/pip install -r requirements.txt
	mkdir $(FIG) $(DATA) 2> /dev/null
	@echo "alias ppdyn='$(CURRENT_DIR)/./$(VENV)/bin/python3 $(CURRENT_DIR)/$(SRC)/main.py'" >> $${HOME}/.zshrc
	. $${HOME}/.zshrc
else
# default target, when make executed without arguments
all: venv
	@echo "Creating virtual environment for running the code"
$(VENV)/bin/activate: requirements.txt
	python3 -m venv $(VENV)
	./$(VENV)/bin/pip install -r requirements.txt
	mkdir $(FIG) $(DATA) 2> /dev/null
endif

# venv is a shortcut target
venv: $(VENV)/bin/activate

run: venv
	@echo "==================================================================="
	@echo "Running PPDyn (Plasma Particle Dynamics)"
	@echo "Author: Dr. Sayan Adhikari, PostDoc @ UiO, Norway"
	@echo "::::::: Dr. Rupak Mukherjee, Associate Research Physicist @ PPPL, NJ"
	@echo "Input: Edit input.ini file to change the parameters for simulation"
	@echo "==================================================================="
	./$(VENV)/bin/python3 $(SRC)/main.py

clean:
	@echo "Cleaning compiled files..."
	rm -rf $(VENV)
	rm -rf $(FIG) $(DATA)
	find . -type f -name '*.pyc' -delete

cleandata:
	@echo "Cleaning data files..."
	find . -type f -name '*.hdf5' -delete

.PHONY: all venv run clean
