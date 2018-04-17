CODE = $(shell find -name '*.f90')
MTMOD = $(shell find -name '*mtmod.f90')
METROPOLIS = $(shell find -name '*Metropolis.f90')
TRANSLATE = $(shell find -name '*Translate.f90')
PBC = $(shell find -name '*PBC.f90')
GENGEOM = $(shell find -name '*gengeom.f90')
ENERGIESMOD = $(shell find -name '*energiesmod.f90')
THETOWER = $(shell find -name '*TheTower.f90')
MAIN = $(shell find -name '*main.f90')
PARAM = $(shell find -name '*param.dat')
DATA = ener.dat traj.xyz
DATA_SCRIPT_PYTHON = ener.dat
SCRIPT_PYTHON = $(shell find -name '*data_visualisation.py')
FIGURE_SCRIPT_PYTHON = pol.png energy.png
TARGET = kmc


$(TARGET) : $(TARGET).o
	gfortran -o $(TARGET) *.o -O3

$(TARGET).o : $(CODE)
	gfortran -c $(MTMOD)
	gfortran -c $(TRANSLATE)
	gfortran -c $(ENERGIESMOD)
	gfortran -c $(PBC)
	gfortran -c $(METROPOLIS)
	gfortran -c $(GENGEOM)
	gfortran -c $(THETOWER)
	gfortran -c $(MAIN)

$(DATA) : $(TARGET) $(PARAM)
	./$(TARGET)

$(FIGURE_SCRIPT_PYTHON) : $(DATA_SCRIPT_PYTHON) $(SCRIPT_PYTHON)
	python2 $(SCRIPT_PYTHON)

## datum : generate data files about MD simulation
.PHONY : datum
datum : $(DATA)

## plot : generate figures about various magnitudes of interest in MD
.PHONY : plot
plot : $(FIGURE_SCRIPT_PYTHON)

## compilation : compilation of the program
.PHONY : compilation
compilation : $(TARGET)

## hardclean : remove ALL auto-generated files including output
.PHONY : hardclean
hardclean :
	@rm -f *.o
	@rm -f *.mod
	@rm -f $(DATA)
	@rm -f $(FIGURE_SCRIPT_PYTHON)
	@rm -f $(TARGET)
	@echo 'All files have been removed'

## clean : remove auto-generated compilation files
.PHONY : clean
clean :
	@rm -f *.o
	@rm -f *.mod
	@echo 'Compilation files removed'

## help : provide some instructions useful to use the makefile
.PHONY : help
help :
	@sed -n 's/^##//p' makefile
