CFLAGS = -g -std=c99

all:	starcalib

starcalib:	starcalib.o amoeba.o hipp8.o
	cc starcalib.o amoeba.o hipp8.o -o starcalib -lm -O3

test:	starcalib
	pngtopnm DSC_0429.png > /tmp/e && \
	./starcalib -f 18 -p 7.9 /tmp/e -o /tmp/f && \
	pnmtopng /tmp/f > /tmp/f.png && \
	open /tmp/f.png

install:
	install starcalib /usr/local/bin
