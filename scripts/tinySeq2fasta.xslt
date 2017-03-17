<?xml version='1.0' ?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<!--
This stylesheet transforms a tinySeq XML file from NCBI to a Fasta file

USAGE: xsltproc tinySeq2fasta.xslt file.xml > output.txt

AUTHOR: Sunit Jain, sunitj [AT] umich [DOT] edu
CREATED: July 2014
-->
<xsl:output method="text"/>

<!-- Match the root node -->
<xsl:template match="/">
<!-- In the next line, Change the value to "nucleotide" if your tinySeq xml contains nucleotides -->
<xsl:for-each select="//TSeq[TSeq_seqtype/@value='protein']">
	<xsl:value-of select="concat('&gt;gi|',TSeq_gi,'|',TSeq_accver,'|',TSeq_defline)"/>
	<xsl:text>&#x9;</xsl:text>
	<xsl:value-of select="TSeq_taxid"/>
	<xsl:text>&#x9;</xsl:text>
	<xsl:value-of select="TSeq_orgname"/>
	<xsl:text>&#x9;</xsl:text>
	<xsl:value-of select="TSeq_length"/>
	<xsl:text>&#10;</xsl:text>
	<xsl:value-of select="TSeq_sequence"/>
	<xsl:text>&#10;</xsl:text>
	</xsl:for-each>
</xsl:template>
</xsl:stylesheet>
