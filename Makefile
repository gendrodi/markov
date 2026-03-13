#required packages on ubuntu:
#libc6 freeglut3 libgl1 libglx0 libglvnd0 libx11-6 libxi6 libxxf86vm1 libxext6 libxcb1 libxau6 libxdmcp6 libbsd0 libmd0

default: markov

CC = gcc
CFLAGS = -lm -lglut -lGLU -lGL -g

.c.o:
	$(CC) $< -c -o $(patsubst %.c,%.o,$<) $(CFLAGS)

HEADER = 			\
	includes.h		\
	ranlxd.h

SOURCE:=			\
	markov.c		\
	graphics_utils.c	\
	rng_gauss.c		\
	ranlxd.c

OBJECT:=		\
	$(patsubst %.c,%.o,$(filter %.c,$(SOURCE)))	\

$(OBJECT): $(HEADER) Makefile

markov: $(OBJECT)
	${CC} -o $@ $(OBJECT) $(CFLAGS)

clean:
	rm -rf *.o *~ *.backup

