.SUFFIXES: .o .c .f90
.PHONY: all clean

all: ${TARGET}

run: ${TARGET}
	./${TARGET}

help:
	@echo 'Makefile for the tracking code                                 '
	@echo '                                                               '
	@echo 'Usage:                                                         '
	@echo '    make all                              Compile tracking code'
	@echo '    make run                         Run the tracking algorithm'
	@echo '    make clean           Remove the executable and object files'
	@echo '    make purge                                Remove the output'
	@echo '                                                               '

${TARGET}: ${OBJ}
	@mkdir -p ${OUTDIR}
	${FORTRAN} ${FFLAGS} -o $@ ${OBJ} ${INCS} ${LIBS}

${OBJDIR}/%.o: ${SRCDIR}/%.f90
	@mkdir -p ${OBJDIR}
	@mkdir -p ${MODDIR}
	${FORTRAN} ${FFLAGS} ${INCS} -c $< -o $@ 

clean:
	-rm -f ${OBJDIR}/*.o
	-rm -f ${MODDIR}/*.mod
	-rm -f ${TARGET}
purge:
	-rm -f ${OUTDIR}/vor*[.txt,.dat]
