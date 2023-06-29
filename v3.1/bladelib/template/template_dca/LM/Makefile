
# $(CC) $(CPPFLAGS) $(CFLAGS) -c
# $(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)

HFILES = LM.h

CFILES = LM.c

SOURCES = ${HFILES} ${CFILES}

OBJECTS = LM.o

SCRIPTS = Compile.sh

default : LM

LM.o : LM.c LM.h

LM : ${OBJECTS}
	${CC} ${LDFLAGS} ${OBJECTS} $(LOADLIBES) $(LDLIBS) -o LM
	
.PHONY: clean
clean:
	rm -f ${OBJECTS} LM

.PHONY: rebuild
rebuild: clean LM

