T1 = MOEA_D
T2 = Build
CF = g++
OP = -O2 
LIBC=$(PWD)/Libs
PD=$(PWD)/Main
all: dir bld comp clean

dir: 
	mkdir -p $(T1) $(T2)

bld: dir 
	$(CF) -I $(LIBC) -c $(LIBC)/aux_structs.cpp $(OP)
	$(CF) -I $(LIBC) -c $(LIBC)/moea_d.cpp $(OP)
	mv *.o $(T2)
	cp $(LIBC)/Toolkit/*.o $(T2)
comp: bld
	$(CF) -I $(LIBC) -o $(T1)/MOEA_D $(PD)/MOEA_D.cpp $(T2)/*.o $(OP)


.PHONY: clean
clean:
	rm -r $(T2)


clean_all:
	rm -r $(T1)
	rm -r $(T2)
	clear