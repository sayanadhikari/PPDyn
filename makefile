##
## @file		makefile
## @brief	PPDyn makefile.
## @author	Dr. Sayan Adhikari <sayan.adhikari@fys.uio.no>
##		Dr. Rupak Mukherjee <rupakm@princeton.edu>
##

# the virtual environment directory
VENV := venv

# Source directory
SRC := src

# Figure directory
FIG := figures
DATA := data
# default target, when make executed without arguments
all: venv
	@echo "Creating virtual environment for running the code"
$(VENV)/bin/activate: requirements.txt
	python3 -m venv $(VENV)
	./$(VENV)/bin/pip install -r requirements.txt
#	mkdir $(FIG) $(DATA) 2> /dev/null

# venv is a shortcut target
venv: $(VENV)/bin/activate

run: venv
	@echo "==================================================================="
	@echo "Running PPDyn (Plasma Particle Dynamics)"
	@echo "Author: Dr. Sayan Adhikari, PostDoc @ UiO, Norway"
	@echo "::::::: Dr. Rupak Mukherjee, Associate Research Physicist @ PPPL, NJ"
	@echo "Input: Edit input.ini file to change the parameters for simulation"
	@echo "==================================================================="
#	find . -type f -name '*.dat' -delete
	./$(VENV)/bin/python3 $(SRC)/main.py

clean:
	@echo "Cleaning compiled files..."
	rm -rf $(VENV)
	rm -rf $(FIG) $(DATA)
	find . -type f -name '*.pyc' -delete

cleandata:
	@echo "Cleaning data files..."
	find . -type f -name '*.dat' -delete

.PHONY: all venv run clean
