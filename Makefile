cpp=c++ -g -O0 -ansi -pedantic -Wall -W -Wno-long-long
#cpp=c++ -O3 -fomit-frame-pointer -funroll-loops -Wno-long-long
cc=cc -O3 -fomit-frame-pointer -funroll-loops
INCLUDE=-Iinclude
LIBRARY=
OUTEXEC=RandomAdd

# Mac OS X
# MAC_UNIVERSAL=-arch i386 -arch ppc -mmacosx-version-min=10.0
# cpp=c++ -g -O3 -fomit-frame-pointer -funroll-loops ${MAC_UNIVERSAL}
# cc=cc -O3 -fomit-frame-pointer -funroll-loops ${MAC_UNIVERSAL}
# OUTEXEC=reroot_trees.macosx

all: ${OUTEXEC}

${OUTEXEC}: main.o rmq.o
	echo 'const char *builddate = "'`date`'";' > buildversion.cpp
	${cpp} -c buildversion.cpp
	${cpp} buildversion.o main.o rmq.o ${INCLUDE} ${LIBRARY} -o ${OUTEXEC}

main.o: main.cpp common.h input.h Makefile
	${cpp} ${INCLUDE} -c $<

rmq.o: rmq.c rmq.h Makefile
	${cc} -c $<

clean:
	rm -f *.o
	rm -f buildversion.*
	rm -f ${OUTEXEC}
	rm -f *~
	rm -f core
	rm -f *.orig

format:
	for i in *.cpp *.h; do astyle --style=kr "$${i}"; done
