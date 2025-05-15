OBJ=  dvcsgen.o aac.o daccepte16emcgen.o  daccepte16prmcgen.o  accep11.o dvcsycol.o jetset.o vpkutil.o accepteg1.o  elacc_newnew.o  mycernlib.o bhradgen.o
OBJC=  dvcs_grid_pack.o 

# -lstdc++ for Gagik's stuff
#
FOR   =  -lstdc++ 
CERNLIBS =  -L/apps/cernlib/x86_64_rhel6_4.7.2/2005/lib  -lmathlib  -lpacklib 
 
dvcsgen : $(OBJ) $(OBJC)
	 gfortran  -o	dvcsgen  $(OBJ) $(OBJC) $(FOR) 
$(OBJ) : %.o: %.F
	gfortran  -std=legacy  -DLinux -fno-automatic  -ffixed-line-length-none -fno-second-underscore  -c $< -o $@  
$(OBJC) : %.o: %.cc
	g++  -DLinux   -c $<  -o $@ 
clean:
	rm -f dvcsgen   $(OBJ)



##g77  -O2 -fno-automatic -finit-local-zero -ffixed-line-length-none -fno-second-underscore \
##        -DLinux \
##        -I. -I./ -I/group/clas/builds/release-4-14/packages/include -I/group/clas/builds/release-4-14/packages/inc_derived -I/apps/tcl/include -I/usr/X11R6/include -c \
##        aac.F  -o /home/avakian/w6/tmp/obj/LinuxRHEL3/generator/aac.o


















