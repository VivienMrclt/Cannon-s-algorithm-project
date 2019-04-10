CC=mpicc
CFLAGS=-W -Wall -std=c99 -g
LDFLAGS= -lm
EXEC=serial cannon

all: $(EXEC)

%: obj/%.o obj/strassen.o
	$(CC) -o $@ $^ $(LDFLAGS)

obj/%.o: src/%.c
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf obj/*.o

mrproper: clean
	rm -rf $(EXEC)
