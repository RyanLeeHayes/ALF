
# $(CC) $(CPPFLAGS) $(CFLAGS) -c
# $(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)

HFILES = Filter.h

CFILES = Filter.c

SOURCES = ${HFILES} ${CFILES}

OBJECTS = Filter.o

SCRIPTS = Compile.sh

default : Filter

Filter.o : Filter.c Filter.h

Filter : ${OBJECTS}
	${CC} ${LDFLAGS} ${OBJECTS} $(LOADLIBES) $(LDLIBS) -o Filter
	
.PHONY: clean
clean:
	rm -f ${OBJECTS} Filter

.PHONY: rebuild
rebuild: clean Filter

