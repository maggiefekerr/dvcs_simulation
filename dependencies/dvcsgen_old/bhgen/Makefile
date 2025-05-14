OBJ=  bhgen.o dvcsycol.o
OBJC=  dvcs_grid_pack.o 

# -lstdc++ for Gagik's stuff
#
FOR   =  -lstdc++ -ldl
CERNLIBS =  -L/Users/sangbaek/CLAS12/cernlib  -lmathlib  -lpacklib 
 
bhgen : $(OBJ) $(OBJC)
	 gfortran -std=legacy  -o	bhgen  $(OBJ) $(OBJC) $(FOR) 
# bhgen : $(OBJ)
# 	 gfortran -std=legacy  -o	bhgen  $(OBJ) $(FOR) 
$(OBJ) : %.o: %.F
	gfortran -std=legacy   -DLinux -fno-automatic  -ffixed-line-length-none -fno-second-underscore  -c $< -o $@  
$(OBJC) : %.o: %.cc
	g++  -DLinux   -c $<  -o $@ 
clean:
	rm -f bhgen   $(OBJ)