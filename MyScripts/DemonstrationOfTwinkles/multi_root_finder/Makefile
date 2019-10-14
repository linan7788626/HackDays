CC	= gcc -Wall -L/usr/local/lib -I/usr/local/include
CFLAGS	= -O3
#
image	= ./nie_images_finder.o ran1.o

all:	a.out

a.out:	$(image)
	$(CC) $(CFLAGS) $(image) -o a.out -lm -lgsl -lgslcblas

# clean out all object code and junk files
#
clean:
	- rm -f a.out \
	*.o *~
