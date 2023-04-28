EXECUTABLE = proj2 
CFLAGS = -std=gnu99 -Wall -Wextra -lpthread -pedantic
CC = gcc

.PHONY: zip

all: $(EXECUTABLE) 

$(EXECUTABLE): proj2.o
	$(CC) $(CFLAGS) -o $@ $^

proj2.o: proj2.c
	$(CC) $(CFLAGS) -c $^

zip: proj2.c
	zip proj2.zip proj2.c Makefile