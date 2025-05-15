# OBJ=  dvcsgen.o aac.o daccepte16emcgen.o  daccepte16prmcgen.o  accep11.o dvcsycol.o jetset.o vpkutil.o accepteg1.o  elacc_newnew.o  mycernlib.o bhradgen.o
# OBJC=  dvcs_grid_pack.o 

# # -lstdc++ for Gagik's stuff
# #
# FOR   =  -lstdc++ 
# CERNLIBS =  -L/apps/cernlib/x86_64_rhel6_4.7.2/2005/lib  -lmathlib  -lpacklib 
 
# dvcsgen : $(OBJ) $(OBJC)
# 	 gfortran  -o	dvcsgen  $(OBJ) $(OBJC) $(FOR) 
# $(OBJ) : %.o: %.F
# 	gfortran  -std=legacy  -DLinux -fno-automatic  -ffixed-line-length-none -fno-second-underscore  -c $< -o $@  
# $(OBJC) : %.o: %.cc
# 	g++  -DLinux   -c $<  -o $@ 
# clean:
# 	rm -f dvcsgen   $(OBJ)




# Compiler settings for JLab ifarm
FC = gfortran
CXX = g++
FFLAGS = -std=legacy -DLinux -fno-automatic -ffixed-line-length-none -fno-second-underscore -fno-align-commons
CXXFLAGS = -O2 -DLinux

# Object files
FOBJS = dvcsgen.o aac.o daccepte16emcgen.o daccepte16prmcgen.o accep11.o \
         dvcsycol.o jetset.o vpkutil.o accepteg1.o elacc_newnew.o mycernlib.o bhradgen.o
COBJS = dvcs_grid_pack.o

all: dvcsgen

dvcsgen: $(FOBJS) $(COBJS)
	$(FC) -o $@ $^ -lstdc++

$(COBJS): %.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(FOBJS): %.o: %.F
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -f dvcsgen *.o








