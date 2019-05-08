CC=mpicc #mpiicc

CFLAGS=-W -Wall -std=c99 -g
LDFLAGS= -lm -g
OPT=-O3 #-mkl

PDC_PARAM=-o GSSAPIDelegateCredentials=yes -o GSSAPIKeyExchange=yes -o GSSAPIAuthentication=yes
EXEC=serial cannon
EXEC_MKL=cannon_mkl
PDC_ID=marcault

all: $(EXEC)

mkl: obj/load_matrix.o $(EXEC_MKL)

%: obj/%.o obj/strassen.o obj/load_matrix.o obj/utils.o
	$(CC) -o $@ $^ $(LDFLAGS) $(OPT)

obj/%.o: src/%.c
	$(CC) -o $@ -c $< $(CFLAGS) $(OPT)

.PHONY: clean mrproper

clean:
	rm -rf obj/*.o

mrproper: clean
	rm -rf $(EXEC) $(EXEC_MKL)

sendPDC:
	scp $(PDC_PARAM) src/* $(PDC_ID)@tegner.pdc.kth.se:~/Private/cannon/src
	scp $(PDC_PARAM) job.sh $(PDC_ID)@tegner.pdc.kth.se:~/Private/cannon/

connectPDC:
	ssh $(PDC_PARAM) $(PDC_ID)@tegner.pdc.kth.se
