# this makefile is intended for g++ on Linux

CC = g++
# CFLAGS = -c -Wall -Wpedantic
CFLAGS = -c -O3
LDFLAGS = -lm -pthread
SOURCES = iofunc.cpp filter.cpp  fourier.cpp  genfunc.cpp  logfunc.cpp rf_module.cpp fm_radio.cpp helper.cpp
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = fm_radio

all: $(EXECUTABLE) clean

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	-rm $(OBJECTS)
