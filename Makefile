# Copyright 2019 Regents of The University of Michigan.

# This file is part of geo-omics-scripts.

# Geo-omics-scripts is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.

# Geo-omics-scripts is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with Geo-omics-scripts.  If not, see <https://www.gnu.org/licenses/>.

package_name = geo-omics-scripts

# shell fragment to give status 0 if all changes are commited
all_committed := ! git status --porcelain 2>/dev/null | grep -q -v '^??'

git_version := $(shell git describe 2>/dev/null | sed s,^release/,,)

# date part of version if there are uncommitted changes
date := $(shell date +%Y%m%d-%H%M)
date_version_appendix := $(shell $(all_committed) || echo -n -$(date))

# Get version, try in order of priority (highest to lowest)
# 0. VERSION environment variablke
# 1. content of VERSION file
# 2. output of `git describe` (plus possibly date)
# 3. version in parent directory
version0 := $(VERSION)
version1 := $(shell cat VERSION 2>/dev/null)
version2 := $(git_version)$(date_version_appendix)
version3 := $(shell basename $$(pwd) | grep -o -P "(?<=$(package_name)-).*")

ifneq ($(MAKECMDGOALS),install-comics-local)
version := $(strip $(if $(version0),\
        $(version0),\
        $(if $(version1),\
                $(version1),\
                $(if $(version2),\
                        $(version2),\
                        $(if $(version3),\
                                $(version3),\
                                $(error "Failed to get version information")\
                        )\
                )\
        )\
))
endif

comics_temp=/tmp/comics_temp

export
.SILENT:

.NOTPARALLEL:

prefix = 
datadir := $(prefix)/share
docdir := $(datadir)/doc
bindir := $(prefix)/bin
mandir := $(datadir)/man
man1dir := $(mandir)/man1
man5dir := $(mandir)/man5
man7dir := $(mandir)/man7

EXTRA_DIST = \
	COPYING \
	Makefile \
	README.md \
	localenv \
	bash-completion \

data_files = \
	phylosiftrc \
	TruSeq3-PE-2+omics.fa \

doc_files = \
	COPYING \
	README.md \
	scripts/COPYRIGHT.tetramer_freqs_esom \

dist_dir = ../$(package_name)-$(version)
html_dirs = \
	html \
	html/_sources \
	html/_static \
	html/_static/fonts \
	html/_static/css \
	html/_static/js

version_module = lib/omics/_version.py


INSTALL = /usr/bin/install
INSTALL_PROGRAM := $(INSTALL)
INSTALL_DATA := $(INSTALL) -m 644

all: sphinx-docs hard-code-version
	cd scripts && $(MAKE) $@

sphinx-docs: custom_style = _static/css/custom.css
sphinx-docs: stylesheet = _build/html/_static/css/theme.css
sphinx-docs:
	if which sphinx-build >/dev/null; then \
	    cd docs && $(MAKE) html man latexpdf; \
	    if ! grep -q noredcode ${stylesheet}; then \
	        echo "Fixing style..."; \
	        cat ${custom_style} >> ${stylesheet}; \
	    fi; \
	else \
	    echo "[WARNING] sphinx-build not available, skipping sphinx-based documentation"; \
	fi

hard-code-version:
	echo "hard-coding version $(version) ..."
	sed -i -r "s/^VERSION.*/VERSION = \'$(version)\'/" $(version_module)

ifneq ($(MAKECMDGOALS),install-comics-local)
# version arithmetic:
# 0. check $version is compatible with `git describe` output
#    with major.minor semantic versioning tags
#    or appended date-time
git_tag_pat = "^\d+\.\d+(-\d+-g[a-f0-9]+)?(-[0-9]{8}-[0-9]{4})?$$"
version_pat = "^\d+\.\d+\.\d+$$"
_good_version := $(if $(shell echo $(version) | grep -P $(git_tag_pat)), $(version), $(error Failed to parse version i.e. output of git describe or content of file VERSION: "$(version)"))
# 1. extract semantic version numbers
_sem_versions := $(subst ., ,$(subst -, ,$(_good_version)))
major_version := $(word 1,$(_sem_versions))
minor_version := $(word 2,$(_sem_versions))
# 2. increment minor version
inc_minor_version := $(shell echo $(minor_version)+1 | bc)
inc_version := $(major_version).$(inc_minor_version)
endif

distmkdir:
	mkdir -p -- "$(dist_dir)"

distdocs: sphinx-docs
	# copy sources
	for p in docs docs/_static/css; do \
	    mkdir -p -- "$(dist_dir)/$$p"; \
	    for i in $$(find $$p -maxdepth 1 -type f -not -name '*~'); do \
	        cp -p $$i "$(dist_dir)/$$p"; \
	    done; \
	done
	# copy sphinx-generated docs to not require full sphinx/latex support on deployment sites
	cp -rp docs/_build "$(dist_dir)/docs/_build"

distdir: distmkdir distdocs
	$(all_committed) || echo "Warning: git reports uncommitted changes; will be included in distribution!"
	cd lib && $(MAKE) $@
	cd scripts && $(MAKE) $@
	$(info Making VERSION file with content: $(version))
	echo "$(version)" > $(dist_dir)/VERSION
	$(info Copying files ...)
	cp -rp $(EXTRA_DIST) $(data_files) $(dist_dir)

dist: distdir
	$(info Creating $(dist_dir).tar.gz)
	tar czhf $(dist_dir).tar.gz $(dist_dir)
	$(info Cleaning up temporary dist directory ...)
	$(RM) -r -- "$(dist_dir)"

allcommitted:
	$(all_committed) || (echo "All changes must be committed to the git repo before making a release!"; exit 1)

inc-version-tag:
	# increment patch version unless $(VERSION) was provided
	# This way we can use
	#    VERSION=1.0.0 make release
	# To make a new major release and
	#    make release
	# to make patch releases.
	$(if $(VERSION),,$(eval version := $(inc_version)))
	git tag -a "release/$(version)" -m "Release version $(version)"
	$(info Version incremented to $(version))

# make a new release:
# 1. increment the version and set as git tag
# 2. build stuff as needed for a release, i.e. sphinx docs
# 3. build the tarball
release: allcommitted inc-version-tag sphinx-docs hard-code-version dist
	# reset hard-coded version file
	git checkout -- lib/omics/_version.py

install-data: installdir = $(DESTDIR)$(datadir)/$(package_name)
install-data:
	$(info Installing data files ...)
	mkdir -p -- "$(installdir)"
	$(INSTALL_DATA) -t $(installdir) $(data_files)

# install sphinx-generated docs and file in doc_files
# html files remain in their directory structure
# all others go flat into their directory $(dest) or man dir
install-docs: dest := $(DESTDIR)$(docdir)/$(package_name)
install-docs: html_dirs = $(shell cd docs/_build && find . -type d -path "./html*")
install-docs: html_files = $(shell cd docs/_build && find . -type f -path "./html*")
install-docs: man_1_pages = $(shell find docs/_build -type f -path "docs/_build/man/*.1")
install-docs: man_7_pages = $(shell find docs/_build -type f -path "docs/_build/man/*.7")
install-docs: pdfs = $(shell find docs/_build -type f -path "docs/_build/latex/*.pdf")
install-docs:
	$(info Creating directories ...)
	mkdir -p -- "$(dest)"
	for i in $(html_dirs); do mkdir -p -- "$(dest)/$$i"; done
	$(info Installing html files ...)
	for i in $(html_files); do $(INSTALL_DATA) "docs/_build/$$i" $(dest)/$$(dirname $$i); done
	$(info Installing manual page ...)
	mkdir -p -- "$(DESTDIR)$(man1dir)"
	$(INSTALL_DATA) ${man_1_pages} $(DESTDIR)$(man1dir)
	mkdir -p -- "$(DESTDIR)$(man7dir)"
	$(INSTALL_DATA) ${man_7_pages} $(DESTDIR)$(man7dir)
	$(info Installing other documentation ...)
	$(INSTALL_DATA) $(doc_files) $(dest)
	$(if $(pdfs), $(INSTALL_DATA) $(pdfs) $(dest), $(warning PDFs missing from documentation))

install: install-docs install-data
	cd lib && $(MAKE) install
	cd scripts && $(MAKE) install
	$(info Installing bash completion ...)
	mkdir -p -- $(DESTDIR)$(datadir)/bash-completion/completions
	$(INSTALL) -t $(DESTDIR)$(datadir)/bash-completion/completions/ bash-completion/omics


standalone-comics: prog = scripts/comics
standalone-comics: insert_tmp := $(shell mktemp)
standalone-comics:
	# local installation of the comics script:
	# 1. get liba.sh, from trap command to end
	# 2. replace the "source liba.sh" call, including shellcheck comments with a mark (BORK42)
	# 3. write insert after mark
	# 4. remove mark
	$(info Building standalone $(prog) ...)
	sed -n '/^trap/,/not accessible/p' lib/liba.sh > $(insert_tmp)
	rm -f -- $(comics_temp)
	cat $(prog) \
	    | sed "/liba.sh/,/liba.sh/c BORK42" \
	    | sed "/BORK42/r $(insert_tmp)" \
	    | sed "/BORK42/d" \
	    > $(comics_temp)
	chmod +x -- $(comics_temp)
	rm -f -- $(insert_tmp)
	$(info done)

install-comics-local: dest = /usr/local/bin/comics
install-comics-local: destdir = $(shell dirname $(dest))
install-comics-local: standalone-comics
	$(info Installing $(comics_temp) to $(dest) ...)
	mkdir -p -- $(destdir)
	$(INSTALL_PROGRAM) $(comics_temp) $(dest)
	chmod +x -- $(dest)
	rm -f -- $(comics_temp)
	$(info done)

install-omics-module: src = modulefiles/flux.omics
install-omics-module: dest_base = /geomicro/data9/flux
install-omics-module: comics_base = $(dest_base)/apps
install-omics-module: module_base = $(dest_base)/modulefiles/geomicro
install-omics-module: standalone-comics
	$(info Installing omics module to $(dest_base) ...)
	[ -d $(comics_base) ] || { echo "missing $(comics_base)"; false; }
	[ -d $(module_base) ] || { echo "missing $(module_base)"; false; }
	mkdir -p -- $(comics_base)/comics
	mkdir -p -- $(comics_base)/comics/bin
	mkdir -p -- $(comics_base)/comics/man
	mkdir -p -- $(module_base)/omics
	$(INSTALL_DATA) -t $(module_base)/omics/ $(src)/1 $(src)/2
	ln -s -f $(module_base)/omics/2 $(module_base)/omics/default
	$(INSTALL_PROGRAM) -T $(comics_temp) $(comics_base)/comics/bin/comics
	rm -f -- $(comics_temp)
	$(info done)

uninstall:
	$(info Removing documentation...)
	$(RM) -r -- $(DESTDIR)/$(docdir)/$(package_name)
	$(info Removing package data...)
	$(RM) -r -- $(DESTDIR)/$(datadir)/$(package_name)
	cd scripts && $(MAKE) $@

clean:
	$(info Cleaning sphinx-generated documentation ...)
	if which sphinx-build >/dev/null; then \
	    cd docs && $(MAKE) clean; \
	fi
	cd scripts && $(MAKE) clean

distclean: clean
	$(info Removing tarballs ...)
	$(RM) -r -- $(package_name)-*

debug:
	$(info share: $(datadir) bin: $(bindir))
	$(info Version 0: $(version0))
	$(info Version 1: $(version1))
	$(info Version 2: $(version2))
	$(info Version 3: $(version3))
	$(info Version: $(version))
	$(info major version: $(major_version))
	$(info minor versions: $(minor_version) --> $(inc_minor_version))
	$(info incremented version: $(inc_version))
	$(info dist_dir: $(dist_dir))
	if $(all_committed); then echo "git: all changed committed"; else echo "git: there are uncommitted changes"; fi

.PHONY: all sphinx-docs scripts distdir dist inc_version_tag release install-docs install uninstall clean distclean debug
