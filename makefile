CC = gcc
OBJS = Main.o allvars.o funzioni.o
CFLAGS = -Wall -Wextra -O3 -ggdb3 -fno-omit-frame-pointer
all: distrGen
distrGen: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o distrGen -lm
Main.o: Main.c
	$(CC) $(CFLAGS) -c Main.c
allvars.o: build/allvars.c build/allvars.h
	$(CC) $(CFLAGS) -c build/allvars.c
funzioni.o: build/funzioni.c build/funzioni.h
	$(CC) $(CFLAGS) -c build/funzioni.c
clean:
	rm -f $(OBJS) $(FOBJS)
