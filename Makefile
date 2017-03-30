export
package_name = geo-omics-scripts
version = $(strip $(shell cat VERSION))

.SILENT:

prefix = 
datadir = $(prefix)/share
docdir = $(datadir)/doc
bindir = $(prefix)/bin
mandir = $(datadir)/man
man1dir = $(mandir)/man1
man5dir = $(mandir)/man5
man7dir = $(mandir)/man7

EXTRA_DIST = \
	COPYRIGHT \
	docs \
	Makefile \
	modulefiles \
	README.md \
	VERSION

doc_files = COPYRIGHT README.md

dist_dir = $(package_name)-$(version)
tarball = $(dist_dir).tar.gz
html_dirs = \
		  html \
		  html/_sources \
		  html/_static \
		  html/_static/fonts \
		  html/_static/css \
		  html/_static/js \


INSTALL = /usr/bin/install
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA = $(INSTALL) -m 644

all: sphinx-docs scripts-man

sphinx-docs:
ifeq ($(shell hostname),csheikMP)
	cd docs && $(MAKE) html man
else
	@echo "-- skipping sphinx-based documentation --"
endif

scripts-man:
	cd scripts && $(MAKE) man

#increment_version = $(file > VERSION,)

distdir:
	mkdir -p $(dist_dir)
	cd lib && $(MAKE) $@
	cd scripts && $(MAKE) $@
	$(info Copying extra files ...)
	cp -a $(EXTRA_DIST) $(dist_dir)

dist: distdir
	$(info Creating $(dist_dir).tar.gz)
	tar czhf $(dist_dir).tar.gz $(dist_dir)
	$(info Cleaning up temporary dist directory ...)
	$(RM) -r $(dist_dir)
	#$(call increment_version)

vondamm-test:
	$(info Copying tarball to $(host) ...)
	scp -p $(tarball) vondamm.geo.lsa.umich.edu:/tmp/heinro/
	$(info Installing remotely ...)
	ssh vondamm.geo.lsa.umich.edu "\
	    cd /tmp/heinro/ && \
	    tar xfz $(tarball) && \
	    cd $(dist_dir) && \
	    make && \
	    DESTDIR=/tmp/heinro/dest make install"

# install sphinx-generated docs and file in doc_files
install-docs: dest = $(DESTDIR)$(docdir)/$(package_name)
install-docs: html_dirs = $(shell cd docs/_build && find html -type d)
install-docs: html_files = $(shell cd docs/_build && find html -type f)
install-docs:
	$(info Creating directories ...)
	mkdir -p $(dest)
	for i in $(html_dirs); do mkdir -p "$(dest)/$$i"; done
	$(info Installing html files ...)
	for i in $(html_files); do $(INSTALL_DATA) "docs/_build/$$i" $(dest)/$$(dirname $$i); done
	$(info Installing other documentation ...)
	$(INSTALL_DATA) $(doc_files) $(dest)
	$(info Installing manual page ...)
	mkdir -p $(DESTDIR)$(man7dir)
	$(INSTALL_DATA) docs/_build/man/geomics.7 $(DESTDIR)$(man7dir)

install: install-docs
	cd lib && $(MAKE) install
	cd scripts && $(MAKE) install

clean:
	$(info Cleaning sphinx-generated documentation ...)
	cd docs && $(MAKE) clean
	cd scripts && $(MAKE) clean

distclean: clean
	$(info Removing tarballs ...)
	$(RM) $(package_name)-*.tar.gz

debug:
	@echo "share: $(datadir) bin: $(bindir)"
	@echo "install: $(INSTALL)"
	@echo $(MAKE)
	@echo $(LDCONFIG)
	@echo $(TEX)

