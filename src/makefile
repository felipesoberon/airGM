CC      = g++
CFLAGS  = -O3 -std=c++11
EXE     = airGM2.1
OBJECTS = globalmodel.o species.o reaction.o commandline.o
INSTDIR = /usr/bin


$(EXE):		main.cpp $(OBJECTS)
		$(CC) $(CFLAGS) main.cpp $(OBJECTS) -o $(EXE)

commandline.o:	commandline.h commandline.cpp
		$(CC) $(CFLAGS) -c commandline.cpp

reaction.o: 	reaction.h reaction.cpp
		$(CC) $(CFLAGS) -c reaction.cpp

species.o:	species.h species.cpp
		$(CC) $(CFLAGS) -c species.cpp

globalmodel.o:	globalmodel.h globalmodel.cpp
		$(CC) $(CFLAGS) -c globalmodel.cpp

clean:
		rm *~ $(EXE) $(OBJECTS)

install:
	cp $(EXE) $(INSTDIR)

uninstall:
	rm $(INSTDIR)/$(EXE)
