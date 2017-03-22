SUBDIRS = scripts
dist_pkgdata_DATA = COPYRIGHT README.md 

EXTRA_DIST = docs/Makefile docs/conf.py docs/*.txt Makefile.local

dist_man_MANS = docs/_build/man/geomics.7

# add extra Makefile fragments
fixme:
	echo "include $(srcdir)/Makefile.local" >> $(srcdir)Makefile.in
	echo "include $(srcdir)/Makefile.local" >> $(srcdir)/scripts/Makefile.in
	@echo "You need to re-run ./configure !!!"

.PHONY: html-local man all-local install-html-local fixme
