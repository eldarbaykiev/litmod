CC=
CFLAGS=-Ofast -fno-automatic -fd-lines-as-code -ffixed-line-length-none -std=legacy
CFLAGS_INTF=

UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CC=gfortran
#	CFLAGS +=
  CFLAGS_INTF += -lpgplot -L/usr/lib -lX11 -fd-lines-as-comments -lpng 

endif
ifeq ($(UNAME), Darwin)
	CC=gfortran
#	CFLAGS +=
  CFLAGS_INTF += -lpgplot -L/opt/X11/lib -I/opt/X11/include -lX11 -L/usr/local/opt/pgplot/lib -I/usr/local/opt/pgplot/include  -fd-lines-as-comments -lpng 
endif

all: litmod litmod_intf litmod2

interface: litmod_intf

litmod:
	$(CC) -o litmod  src/modules.for src/LITMOD3D.for src/SUB* $(CFLAGS)
	cp src/conductionNd_serial.py conductionNd_serial.py
	cp src/temperature_solver.py temperature_solver.py

	$(CC) -o gravcalc_parallel src/modules.for src/gravcalc_parallel.for src/SUB_Geo_Grad3D.for src/SUB_GeoGrav_Grad3D.for src/SUB_Grav_Grad3D.for src/SUB_SumTan.for src/SUB_U_SECOND_DER.for
	cp src/gravity_calculator.py gravity_calculator.py
	
	cp src/periods.py periods.py
	
litmod2:
	$(CC) -o litmod2 LITMOD3D_FOR_STATOIL_last/LITMOD3D.for LITMOD3D_FOR_STATOIL_last/SUB* $(CFLAGS)

litmod_intf:
	$(CC) -o litmod_intf src_intf/LITMOD3D_INTF.f src_intf/SUB* $(CFLAGS_INTF)
	cp src_intf/minmax.py minmax.py

clean:
	rm litmod gravcalc_parallel conductionNd_serial.py temperature_solver.py periods.py gravity_calculator.py litmod_intf minmax.py litmod2
