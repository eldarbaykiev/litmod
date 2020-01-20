CC=
CFLAGS=-Ofast -fno-automatic -fd-lines-as-code -ffixed-line-length-none
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
  CFLAGS_INTF += -lpgplot -L/opt/X11/lib -I/opt/X11/include -lX11 -fd-lines-as-comments -lpng
endif

all: litmod litmod_intf

interface: litmod_intf

litmod:
	$(CC) -o litmod  src/modules.for src/LITMOD3D.for src/SUB* $(CFLAGS)
	cp src/conductionNd_serial.py conductionNd_serial.py
	cp src/temperature_solver.py temperature_solver.py

	$(CC) -o gravcalc_parallel src/modules.for src/gravcalc_parallel.for src/SUB_Geo_Grad3D.for src/SUB_GeoGrav_Grad3D.for src/SUB_Grav_Grad3D.for src/SUB_SumTan.for src/SUB_U_SECOND_DER.for
	cp src/gravity_calculator.py gravity_calculator.py

litmod_intf:
	$(CC) -o litmod_intf src_intf/LITMOD3D_INTF.f src_intf/SUB* $(CFLAGS_INTF)

clean:
	rm litmod gravcalc_parallel conductionNd_serial.py temperature_solver.py gravity_calculator.py litmod_intf
