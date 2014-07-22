<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text"/>
<xsl:template match="/">
<xsl:for-each select="//TSeq[TSeq_seqtype/@value='protein']">
	<xsl:value-of select="concat('&gt;gi|',TSeq_gi,'|',TSeq_accver,'|Taxa:',TSeq_taxid,'|',TSeq_defline)"/>
	<xsl:text>&#10;</xsl:text>
	<xsl:value-of select="TSeq_sequence"/>
	<xsl:text>&#10;</xsl:text>
	</xsl:for-each>
</xsl:template>
</xsl:stylesheet>
