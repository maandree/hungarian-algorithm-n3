SHELL := /bin/bash

OPTIMISE = -O0 -g
WARN = -Wall -Wextra -pedantic -Wdouble-promotion -Wformat=2 -Winit-self -Wmissing-include-dirs  \
       -Wtrampolines -Wfloat-equal -Wshadow -Wmissing-prototypes -Wmissing-declarations          \
       -Wredundant-decls -Wnested-externs -Winline -Wno-variadic-macros -Wsync-nand              \
       -Wunsafe-loop-optimizations -Wcast-align -Wstrict-overflow                                \
       -Wundef -Wbad-function-cast -Wcast-qual -Wwrite-strings -Wlogical-op -Waggregate-return   \
       -Wstrict-prototypes -Wold-style-definition -Wpacked -Wvector-operation-performance        \
       -Wunsuffixed-float-constants -Wsuggest-attribute=const -Wsuggest-attribute=noreturn       \
       -Wsuggest-attribute=pure -Wsuggest-attribute=format -Wnormalized=nfkc -Wconversion        \
       -fstrict-aliasing -fstrict-overflow -ftree-vrp                                            \
        -funsafe-loop-optimizations
# excluded: -Wdeclaration-after-statement


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

