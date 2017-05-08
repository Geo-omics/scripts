package_name = geo-omics-scripts
# get version from  VERSION file is available
# else use `git describe`
# but override with $VERSION in environment
version := $(strip $(if \
	$(VERSION), \
	$(VERSION), \
	$(shell cat VERSION 2>/dev/null \
		|| git describe 2>/dev/null \
		|| echo 1>&2 "Failed to get version information, no VERSION file and no git tags!" \
	) \
))

# give status 0 if there are uncommited changes
all_committed := ! git status --porcelain | grep -q -v '^??'

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
	COPYRIGHT \
	docs \
	Makefile \
	modulefiles \
	README.md

doc_files = COPYRIGHT README.md

dist_dir = $(package_name)-$(version)
html_dirs = \
	html \
	html/_sources \
	html/_static \
	html/_static/fonts \
	html/_static/css \
	html/_static/js


INSTALL = /usr/bin/install
INSTALL_PROGRAM := $(INSTALL)
INSTALL_DATA := $(INSTALL) -m 644

all: sphinx-docs scripts-man

sphinx-docs: custom_style = _static/css/custom.css
sphinx-docs: stylesheet = _build/html/_static/css/theme.css
sphinx-docs:
	if which sphinx-build >/dev/null; then \
	    cd docs && $(MAKE) html man; \
	    if ! grep -q noredcode ${stylesheet}; then \
	        echo "Fixing style..."; \
	        cat ${custom_style} >> ${stylesheet}; \
	    fi; \
	else \
	    echo "[WARNING] sphinx-build not available, skipping sphinx-based documentation"; \
	fi

scripts-man:
	cd scripts && $(MAKE) man

# version arithmetic:
# 0. check $version is compatible with `git describe` output
#    with 1.2.3 sematic versioning tags
git_tag_pat = "^\d+\.\d+\.\d+(-\d+-g[a-f0-9]+)?$$"
version_pat = "^\d+\.\d+\.\d+$$"
_good_version := $(if $(shell echo $(version) | grep -P $(git_tag_pat)), $(version), $(error Failed to parse version i.e. output of git describe or content of file VERSION: $(version)))
# 1. extract semantic version numbers
_sem_versions := $(subst ., ,$(subst -, ,$(_good_version)))
major_version := $(word 1,$(_sem_versions))
minor_version := $(word 2,$(_sem_versions))
patch_version := $(word 3,$(_sem_versions))
# 2. increment patch level
inc_patch_version := $(shell echo $(patch_version)+1 | bc)
inc_version := $(major_version).$(minor_version).$(inc_patch_version)

distdir:
	$(all_committed) || echo "Warning: git reports uncommitted changes; will be included in distribution!"
	mkdir -p -- "$(dist_dir)"
	cd lib && $(MAKE) $@
	cd scripts && $(MAKE) $@
	$(info Making VERSION file with content: $(version))
	echo "$(version)" > $(dist_dir)/VERSION
	$(info Copying extra files ...)
	cp -a $(EXTRA_DIST) $(dist_dir)

dist: distdir
	$(info Creating $(dist_dir).tar.gz)
	tar czhf $(dist_dir).tar.gz $(dist_dir)
	$(info Cleaning up temporary dist directory ...)
	$(RM) -r -- "$(dist_dir)"

inc-version-tag:
	$(eval version := $(inc_version))
	! git status --porcelain | grep -q '^A' && \
	git tag -a "$(version)" -m "Release version $(version)"
	$(info Version incremented to $(version))

release: inc-version-tag dist

# install sphinx-generated docs and file in doc_files
# html files remain in their directory structure
# all others go flat into their directory $(dest) or man dir
install-docs: dest := $(DESTDIR)$(docdir)/$(package_name)
install-docs: html_dirs := $(shell cd docs/_build && find . -type d -path "./html*")
install-docs: html_files := $(shell cd docs/_build && find . -type f -path "./html*")
install-docs: man_1_pages := $(shell find docs/_build -type f -path "docs/_build/man/*.1")
install-docs: man_7_pages := $(shell find docs/_build -type f -path "docs/_build/man/*.7")
install-docs:
	$(info Creating directories ...)
	mkdir -p -- "$(dest)"
	for i in $(html_dirs); do mkdir -p -- "$(dest)/$$i"; done
	$(info Installing html files ...)
	for i in $(html_files); do $(INSTALL_DATA) "docs/_build/$$i" $(dest)/$$(dirname $$i); done
	$(info Installing manual page ...)
	# mkdir -p -- "$(DESTDIR)$(man1dir)"
	# $(INSTALL_DATA) ${man_1_pages} $(DESTDIR)$(man1dir)
	mkdir -p -- "$(DESTDIR)$(man7dir)"
	$(INSTALL_DATA) ${man_7_pages} $(DESTDIR)$(man7dir)
	$(info Installing other documentation ...)
	$(INSTALL_DATA) $(doc_files) $(dest)

install: install-docs
	cd lib && $(MAKE) install
	cd scripts && $(MAKE) install

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
	$(RM) -- $(package_name)-*.tar.gz

debug:
	$(info share: $(datadir) bin: $(bindir))
	$(info Version: $(version))
	$(info major version: $(major_version))
	$(info minor version: $(minor_version))
	$(info patch versions: $(patch_version) --> $(inc_patch_version))
	$(info incremented version: $(inc_version))
	$(info dist_dir: $(dist_dir))
	if $(all_committed); then echo "git: all changed committed"; else echo "git: there are uncommitted changes"; fi

.PHONY: all sphinx-docs scripts-man distdir dist inc_version_tag release install-docs install uninstall clean distclean debug
