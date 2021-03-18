CC         = g++ 
CFLAGS     = -g -Wall
CFLAGSROOT = `root-config --cflags`
LIBSROOT   = `root-config --glibs`

all: main

IElica.o: IElica.cpp
	$(CC) $(CFLAGS) -c IElica.cpp           $(CFLAGSROOT)

IRilevatore.o: IRilevatore.cpp
	$(CC) $(CFLAGS) -c IRilevatore.cpp           $(CFLAGSROOT)

IRiemannFitter.o: IRiemannFitter.cpp
	$(CC) $(CFLAGS) -c IRiemannFitter.cpp           $(CFLAGSROOT)

IPatternFinder.o: IPatternFinder.cpp
	$(CC) $(CFLAGS) -c IPatternFinder.cpp           $(CFLAGSROOT)

Main: main.cpp IElica.o IRilevatore.o IRiemannFitter.o IPatternFinder.o
	$(CC) $(CFLAGS) -o main main.cpp IElica.o IRilevatore.o IRiemannFitter.o IPatternFinder.o $(CFLAGSROOT) $(LIBSROOT)

clean:
	rm *.o
