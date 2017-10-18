all:

# get all python files
py_files = $(shell find -name "*.py")

lib_files = liba.sh

EXTRA_DIST = Makefile

install: installdir = $(DESTDIR)$(datadir)/$(package_name)
install: install-py
	$(info Installing lib files ...)
	mkdir -p -- "$(installdir)"
	$(INSTALL_DATA) -t $(installdir) $(lib_files)

install-py: installdir = $(DESTDIR)$(prefix)/lib/python3.5/site-packages
install-py:
	$(info Installing python packages ...)
	for i in $(py_files); do \
	    $(INSTALL_DATA) -D $$i $(installdir)/$$i; \
	done

distdir:
	$(info Copying lib files ...)
	mkdir -p -- "../$(dist_dir)/lib"
	cp -a $(lib_files) $(EXTRA_DIST) ../$(dist_dir)/lib/
	# copy each python package individually
	for i in $(shell find -name __init__.py -printf "%h "); do \
	    mkdir -p ../$(dist_dir)/lib/$$i && \
	    cp -p $$i/*.py ../$(dist_dir)/lib/$$i/ ; \
	done
