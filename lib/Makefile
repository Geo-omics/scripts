all:

lib_files = liba.sh

EXTRA_DIST = Makefile

install: installdir = $(DESTDIR)$(datadir)/$(package_name)
install:
	$(info Installing lib files ...)
	mkdir -p $(installdir)
	$(INSTALL_DATA) -t $(installdir) $(lib_files)

distdir:
	$(info Copying lib files ...)
	mkdir -p ../$(dist_dir)/lib
	cp -a $(lib_files) $(EXTRA_DIST) ../$(dist_dir)/lib/
