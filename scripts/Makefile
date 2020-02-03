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

PERL_SCRIPTS = \
    addFileName2header\
    addInfo2lrn\
    ani\
    assemble.pl\
    bamTools\
    basicHF\
    batchBlast\
    binTablesForIMG\
    blastDensityPlot\
    calcN50\
    changeClasses\
    changePGDBattribs\
    chopper\
    clusterDensity\
    consolidateJGIdata\
    contigMetadata\
    countInstances\
    coveragePerBin\
    coveragePerScaffold\
    createFastq\
    createNodes\
    CRISPR_spacer_extractor\
    curateDB\
    derep+alias\
    derep_ClusterMap\
    derep_getReadAbundance\
    dereplicate\
    embl2picture\
    esomCodonMod\
    esomTrain\
    esomWrapper\
    expandGFF\
    extract_Blast_Hits_Of_Interest\
    extractContigReads\
    extractGenbankMetadata\
    extractGenomes\
    extractSeqs\
    extractSubSeq\
    extractTranslationsFromGbk\
    fileChopper\
    findStretchesOfNs\
    folderLevelSize\
    fragRec\
    gbk2fna\
    gcSkew\
    genomeCheck\
    genomicFluidity\
    getBwaMappedReadList\
    getClassFasta\
    getCol\
    getFamilyP\
    getFastaFromAccNos\
    getGeneClusters\
    getGFF\
    getGiInfo\
    getGI\
    getLineage\
    getMasterList\
    getMyContigs\
    getRandomData\
    getSciNames\
    gff2fasta\
    gff2neo\
    gff2tbl\
    GI_info_XMLParser\
    iClust\
    img_Bin_Classifier\
    inflate\
    interleave\
    itemize\
    kmerFreq\
    legacy_consolidateJGIdata\
    length+GC\
    limit2Length\
    makeAnnotationFile\
    mapper_getQueryList\
    mapper\
    map_project_names\
    matchQueryNames\
    measureCompleteness\
    nameClassFiles\
    nsmpReport\
    parallel_antiSmash\
    parallel_getGenomesFromTaxa\
    parseBlastXML\
    parseFastq\
    patchBlastLineage\
    plotCoverage\
    postBlast\
    ppt_getGI\
    ppt_getXML\
    removeBlastSubj\
    removeCommentLines\
    renameHeaders\
    reverse_complement\
    sangerSeqParser\
    silvaTaxonAppend\
    summarize_antiSmash\
    tally\
    tally-weave\
    tallyWrap\
    taxonDist\
    test_fragRec\
    tetramer_freqs_esom\
    toMultiGBK\
    top5\
    toPhylipAndBack\
    triage\
    U2T\
    uClustHomology\
    usageStats\

OTHER_EXECUTABLES = \
    antiSmash_summary \
    assemblyModules \
    bins2fasta \
    createPhgDB \
    comics \
    do2folder \
    do2list \
    Ebot.Output.Extract.Gi.Title.Rev3 \
    extractEuks \
    firefox_already_running \
    getGIAnnotation \
    getGISummary \
    merge-covs\
    mVelvetPipe_paired \
    mVelvetPipe_singles \
    oasesPaired_pipe \
    omics \
    omics-assemble \
    omics-binning \
    omics-container \
    omics-mapping \
    omics-prep \
    omics-run \
    omics-init \
    omics-qc \
    omics-qc-check \
    omics-qc-sample \
    refseq-rna \
    remove_space_from_filenames \
    rgi-coverage \
    rgi-setup \
    setup_metapathways \
    silva-db \
    unchop-contigs \
    ESOM_binning_results_parser \
    Metabat_to_anvio_parser \
    VizBin_parser \

# scripts for which help2man can make man pages
HELP2MAN_EXECUTABLES = \
    aview\
    chop-contigs\
    dada2shared\
    downsample\
    match-dada2-mothur\
    merge-coverage\
    merge-covs\
    plot-blast-frag-cov\
    plot-megahit-log\
    plot_alignment\
    plot-shared-otu-counts\
    separate-interleaved\
    silva-align\

NON_EXECUTABLES = \
    parseTinySeqXML.xslt\
    slideshow.xml\
    twitterscript.xml\
    tinySeq2fasta.xslt\
    tinySeq2table.xslt

EXTRA_DIST = \
    Makefile \
    COPYRIGHT.tetramer_freqs_esom \

# perl scripts for which pod2man fails
BAD_PODS = \
    calcN50\
    createFastq\
    curateDB\
    derep_getReadAbundance\
    embl2picture\
    extractSeqs\
    folderLevelSize\
    genomeCheck\
    getGISummary\
    getMasterList\
    getMyContigs\
    GI_info_XMLParser\
    iClust\
    kmerFreq\
    patchBlastLineage\
    plotCoverage\
    removeBlastSubj\
    removeCommentLines\
    tally\
    tally-weave\

# compute names for manpages
bad_re = $(foreach i,$(BAD_PODS),"-F $(i)")
perl_scripts_good = $(filter-out $(BAD_PODS),$(PERL_SCRIPTS))
man_pages_pod = $(foreach i,$(perl_scripts_good),$(i).1)
man_pages_h2m = $(foreach i,$(HELP2MAN_EXECUTABLES),$(i).1)

$(man_pages_pod): name = $(subst .1,,$@)
$(man_pages_pod): version ?= unknown_version
$(man_pages_pod):
	$(info pod manpage for $(name))
	@pod2man --utf8 --release=$(version) --center $(name) $(name) $@

# help2man: for python scripts with argparse
# $(short) should be the first line of the file's __doc__ string
# The --version option can display the version as VERSION_PLACEHOLDER that can
# be replaces with the real version when distributing (see also distdir target)
$(man_pages_h2m): name = $(subst .1,,$@)
$(man_pages_h2m): short = $(shell ./$(name) --help | grep -m1 -A1 "^$$" | tail -n1)
$(man_pages_h2m): version ?= unknown_version
$(man_pages_h2m):
	$(info help2man manpage for $(name))
	@help2man --locale en_US.UTF-8 --no-info --name "$(short)" --source "geo-omics-scripts $(version)" --output $(name).1 ./$(name)

man: $(man_pages_pod) $(man_pages_h2m)

all: man

distdir: version ?= unknown_version
distdir: dst = ../$(dist_dir)/scripts/
distdir:
	$(info Copying scripts ...)
	mkdir -p -- "$(dst)"
	cp -a $(PERL_SCRIPTS) $(OTHER_EXECUTABLES) $(HELP2MAN_EXECUTABLES) $(NON_EXECUTABLES) $(EXTRA_DIST) $(dst)
	$(info Replacing version placeholder in some scripts ...)
	for i in $(HELP2MAN_EXECUTABLES); do sed "s/VERSION_PLACEHOLDER/$(version)/" $$i > $(dst)/$$i; chmod +x $(dst)/$$i; done

install: install-scripts install-man

install-scripts:
	mkdir -p -- "$(DESTDIR)$(bindir)"
	$(info Installing perl scripts ...)
	$(INSTALL_PROGRAM) $(PERL_SCRIPTS) $(DESTDIR)$(bindir)
	$(info Installing other executable scripts ...)
	$(INSTALL_PROGRAM) $(OTHER_EXECUTABLES) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(HELP2MAN_EXECUTABLES) $(DESTDIR)$(bindir)
	mkdir -p -- "$(DESTDIR)$(datadir)/$(package_name)"
	$(info Installing non-executable scripts ...)
	$(INSTALL_DATA) $(NON_EXECUTABLES) $(DESTDIR)$(datadir)/$(package_name)

install-man: man
	$(info Installing generated POD and help2man man pages ...)
	mkdir -p -- "$(DESTDIR)$(man1dir)"
	$(INSTALL_DATA) $(man_pages_pod) $(DESTDIR)$(man1dir)
	$(INSTALL_DATA) $(man_pages_h2m) $(DESTDIR)$(man1dir)

uninstall:
	$(info Removing executable scripts...)
	cd $(DESTDIR)$(bindir) && $(RM) -f -- $(PERL_SCRIPTS) $(OTHER_EXECUTABLES)
	$(info Removing executable perl script man pages...)
	cd $(DESTDIR)$(man1dir) && $(RM) -f -- $(man_pages_pod)

clean:
	$(info Removing generated man pages ...)
	$(RM) -- *.1

debug:
	@echo good Perl scripts found:
	@echo $(perl_scripts_good)
	@echo
	@echo Corresponding man pages:
	@echo $(man_pages_pod)
	@echo
	@echo Bad pods grep args:
	@echo $(bad_re)

.PHONY: all man distdir install install-scripts install-man clean debug
