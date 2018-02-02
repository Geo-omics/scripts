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
    chop-contigs \
    createPhgDB \
    do2folder \
    do2list \
    downsample \
    Ebot.Output.Extract.Gi.Title.Rev3 \
    extractEuks \
    firefox_already_running \
    getGIAnnotation \
    getGISummary \
    merge-coverage \
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
    separate-interleaved\
    setup_metapathways \
    silva-db \

NON_EXECUTABLES = \
    parseTinySeqXML.xslt\
    slideshow.xml\
    twitterscript.xml\
    tinySeq2fasta.xslt\
    tinySeq2table.xslt

EXTRA_DIST = \
    Makefile

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
man_pages = $(foreach i,$(perl_scripts_good),$(i).1)

$(man_pages): name = $(subst .1,,$@)
$(man_pages): version ?= unknown_version
$(man_pages):
	@pod2man --utf8 --release=$(version) --center $(name) $(name) $@
	$(info manpage for $(name))

man: $(man_pages)

all: man

distdir:
	$(info Copying scripts ...)
	mkdir -p -- "../$(dist_dir)/scripts"
	cp -a $(PERL_SCRIPTS) $(OTHER_EXECUTABLES) $(NON_EXECUTABLES) $(EXTRA_DIST) ../$(dist_dir)/scripts/

install: install-scripts install-man

install-scripts:
	mkdir -p -- "$(DESTDIR)$(bindir)"
	$(info Installing perl scripts ...)
	$(INSTALL_PROGRAM) $(PERL_SCRIPTS) $(DESTDIR)$(bindir)
	$(info Installing other executable scripts ...)
	$(INSTALL_PROGRAM) $(OTHER_EXECUTABLES) $(DESTDIR)$(bindir)
	mkdir -p -- "$(DESTDIR)$(datadir)/$(package_name)"
	$(info Installing non-executable scripts ...)
	$(INSTALL_DATA) $(NON_EXECUTABLES) $(DESTDIR)$(datadir)/$(package_name)

install-man: man
	$(info Installing man pages for perl scripts ...)
	mkdir -p -- "$(DESTDIR)$(man1dir)"
	$(INSTALL_DATA) $(man_pages) $(DESTDIR)$(man1dir)

uninstall:
	$(info Removing executable scripts...)
	cd $(DESTDIR)$(bindir) && $(RM) -f -- $(PERL_SCRIPTS) $(OTHER_EXECUTABLES)
	$(info Removing executable perl script man pages...)
	cd $(DESTDIR)$(man1dir) && $(RM) -f -- $(man_pages)

clean:
	$(info Removing generated man pages ...)
	$(RM) -- *.1

debug:
	@echo good Perl scripts found:
	@echo $(perl_scripts_good)
	@echo
	@echo Corresponding man pages:
	@echo $(man_pages)
	@echo
	@echo Bad pods grep args:
	@echo $(bad_re)

.PHONY: all man distdir install install-scripts install-man clean debug
