all:
	gcc -g -o "hungarian"{,.c}

test:
	./"hungarian"

valgrind:
	valgrind --leak-check=full ./"hungarian"

clean:
	if [ -f "hungarian" ]; then  unlink "hungarian";  fi


.PHONY: test valgrind clean

