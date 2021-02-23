.POSIX:

CPPFLAGS = -D_DEFAULT_SOURCE -D_BSD_SOURCE -D_XOPEN_SOURCE=700
CFLAGS   = -std=c99 -g
LDFLAGS  =


all: hungarian

hungarian: hungarian.c
	$(CC) -o $@ hungarian.c $(CFLAGS) $(CPPFLAGS) $(LDFLAGS)

clean:
	-rm -f -- hungarian


.PHONY: all clean
