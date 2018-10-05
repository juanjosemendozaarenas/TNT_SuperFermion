include makefile.inc

OBJMAIN := $(patsubst %,$(OBJDIR)/%,$(OBJMAIN))

.PHONY: all clean $(SUBDIRS) $(MKFILE)

$(MKFILE): all

all: ${OBJMAIN} $(SUBDIRS) $(EXDIR)/main_tebd_thermal_machine

$(EXDIR)/main_tebd_thermal_machine: $(OBJ)
	${FF} -o $@ $(OBJ) $(NAGFLAGS) $(CLIBS) $(FFLAGS)

$(SUBDIRS):
	mkdir --parents $(OBJDIR)/$@
	cd $@; $(MAKE) $(MKFLAGS) $(MKFILE)

$(OBJDIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	-rm $(OBJ) $(EXDIR)/main_tebd_thermal_machine
