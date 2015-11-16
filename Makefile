SHELL := /bin/bash

OPTIMISE = -O0 -g
WARN = -Wall -Wextra -pedantic -Wformat=2 -Winit-self -Wmissing-include-dirs   \
	   -Wfloat-equal -Wshadow -Wmissing-prototypes -Wmissing-declarations      \
	   -Wredundant-decls -Wnested-externs -Winline -Wno-variadic-macros        \
	   -Wcast-align -Wstrict-overflow -Wundef -Wbad-function-cast -Wcast-qual  \
	   -Wwrite-strings -Waggregate-return -Wstrict-prototypes                  \
	   -Wold-style-definition -Wpacked -Wconversion -fstrict-aliasing          \
	   -fstrict-overflow -Wdeclaration-after-statement


all: hungarian

hungarian: hungarian.c
	gcc -std=gnu99 $(OPTIMISE) $(WARN) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $<

test: hungarian
	./$<

valgrind: hungarian
	valgrind --tool=memcheck --leak-check=full ./$<

clean:
	rm -f hungarian


.PHONY: all test valgrind clean

