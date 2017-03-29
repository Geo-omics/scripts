export
package_name = geo-omics-scripts
version = $(strip $(shell cat VERSION))

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

export dist_dir = $(package_name)-$(version)
tarball = $(dist_dir).tar.gz

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
	cp -a $(EXTRA_DIST) $(dist_dir)
	cd lib && $(MAKE) $@
	cd scripts && $(MAKE) $@



dist: distdir
	tar czhf $(dist_dir).tar.gz $(dist_dir)
	$(RM) -r $(dist_dir)
	#$(call increment_version)

vondamm-test:
	scp -p $(tarball) vondamm.geo.lsa.umich.edu:/tmp/heinro/
	ssh vondamm.geo.lsa.umich.edu "\
	    cd /tmp/heinro/ && \
	    tar xfz $(tarball) && \
	    cd $(dist_dir) && \
	    make && \
	    DESTDIR=/tmp/heinro/dest make install"

install-docs:
	@echo TODO: install sphinx generated docs

install: install-docs
	cd lib && $(MAKE) install
	cd scripts && $(MAKE) install

clean:
	cd docs && $(MAKE) clean
	cd scripts && $(MAKE) clean

distclean: clean
	$(RM) $(package_name)-*.tar.gz

debug:
	@echo "share: $(datadir) bin: $(bindir)"
	@echo "install: $(INSTALL)"
	@echo $(MAKE)
	@echo $(LDCONFIG)
	@echo $(TEX)

