CC=mpicc
CFLAGS=-W -Wall -std=c99 -g
LDFLAGS= -lm -g
EXEC=serial cannon cannon_init_transfer

all: $(EXEC)

%: obj/%.o obj/strassen.o obj/load_matrix.o
	$(CC) -o $@ $^ $(LDFLAGS)

obj/%.o: src/%.c
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf obj/*.o

mrproper: clean
	rm -rf $(EXEC)
