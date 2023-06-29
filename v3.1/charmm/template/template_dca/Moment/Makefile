
# $(CC) $(CPPFLAGS) $(CFLAGS) -c
# $(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)

HFILES = ObservedMoments.h

CFILES = ObservedMoments.c

SOURCES = ${HFILES} ${CFILES}

OBJECTS = ObservedMoments.o

SCRIPTS = Compile.sh

default : ObservedMoments

ObservedMoments.o : ObservedMoments.c ObservedMoments.h

ObservedMoments : ${OBJECTS}
	${CC} ${LDFLAGS} ${OBJECTS} $(LOADLIBES) $(LDLIBS) -o ObservedMoments
	
.PHONY: clean
clean:
	rm -f ${OBJECTS} ObservedMoments

.PHONY: rebuild
rebuild: clean ObservedMoments

