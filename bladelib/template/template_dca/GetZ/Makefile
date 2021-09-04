
# $(CC) $(CPPFLAGS) $(CFLAGS) -c
# $(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)

HFILES = GenerateMoments.h

CFILES = GenerateMoments.c

SOURCES = ${HFILES} ${CFILES}

OBJECTS = GenerateMoments.o

SCRIPTS = Compile.sh

default : GenerateMoments

GenerateMoments.o : GenerateMoments.c GenerateMoments.h

GenerateMoments : ${OBJECTS}
	${CC} ${LDFLAGS} ${OBJECTS} $(LOADLIBES) $(LDLIBS) -o GenerateMoments
	
.PHONY: clean
clean:
	rm -f ${OBJECTS} GenerateMoments

.PHONY: rebuild
rebuild: clean GenerateMoments

