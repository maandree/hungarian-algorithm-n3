all:
	gcc -o "hungarian"{,.c}

nodebug:
	gcc -o "hungarian"{,.c}

test:
	./"hungarian"

valgrind:
	valgrind --tool=memcheck --leak-check=full ./"hungarian"

clean:
	if [ -f "hungarian" ]; then  unlink "hungarian";  fi


.PHONY: test valgrind clean

