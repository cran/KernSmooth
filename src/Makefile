LIBNAME=KernSmooth
LD=ld
# if you built R with -lblas you may need
#LIBS=-lblas

OBJS=   blkest.o cp.o dgedi.o linbin.o linbin2D.o locpoly.o rlbin.o \
        sdiag.o sstdiag.o dgefa.o dgesl.o

$(LIBNAME).so:  $(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS) $(LIBS)

clean:
	@rm -f *.o *.so

realclean: clean
