CFLAGS+=-Wall -Werror -pedantic -std=gnu99 -march=native -O2 -flto
LDLIBS+=-lm

all: main
main: rtweekend.o

clean:
	rm main *.o
