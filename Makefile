CFLAGS+=-Wall -Werror -pedantic -std=gnu99
LDLIBS+=-lm

all: main
main: rtweekend.o

clean:
	rm main *.o
