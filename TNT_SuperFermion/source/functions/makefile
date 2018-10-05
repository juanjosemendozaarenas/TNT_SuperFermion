INCROOT := ../

include ${INCROOT}makefile.inc

# Directory for obj files
OBJDIR := ../$(OBJDIR)

OBJstar := $(patsubst %,$(OBJDIR)/%,$(OBJstar))

.PHONY: functions

functions: ${OBJstar}

${OBJDIR}/functions/%.o: %.c ${DEPSstar}
	$(CC) -c -o $@ $< $(CFLAGS)
