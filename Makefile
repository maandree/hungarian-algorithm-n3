OPTIMISE = -Og -g
WARN = -Wall -Wextra -pedantic -Wdouble-promotion -Wformat=2 -Winit-self -Wmissing-include-dirs  \
       -Wtrampolines -Wfloat-equal -Wshadow -Wmissing-prototypes -Wmissing-declarations          \
       -Wredundant-decls -Wnested-externs -Winline -Wno-variadic-macros -Wsync-nand              \
       -Wunsafe-loop-optimizations -Wcast-align -Wstrict-overflow                                \
       -Wundef -Wbad-function-cast -Wcast-qual -Wwrite-strings -Wlogical-op -Waggregate-return   \
       -Wstrict-prototypes -Wold-style-definition -Wpacked -Wvector-operation-performance        \
       -Wunsuffixed-float-constants -Wsuggest-attribute=const -Wsuggest-attribute=noreturn       \
       -Wsuggest-attribute=pure -Wsuggest-attribute=format -Wnormalized=nfkc -Wconversion        \
       -fstrict-aliasing -fstrict-overflow -fipa-pure-const -ftree-vrp -fstack-usage             \
        -funsafe-loop-optimizations
# excluded: -Wdeclaration-after-statement


all: hungarian

hungarian: hungarian.c
	gcc -std=gnu99 $(OPTIMISE) $(WARN) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $<

test:
	./"hungarian"

valgrind:
	valgrind --tool=memcheck --leak-check=full ./"hungarian"

clean:
	-rm hungarian


.PHONY: all test valgrind clean

