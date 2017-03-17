<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
<xsl:output method="text"/>
<xsl:template match="/">
<xsl:for-each select="//TSeq[TSeq_seqtype/@value='protein']">
	<xsl:value-of select="TSeq_gi"/>
        <xsl:text>&#x9;</xsl:text>
        <xsl:value-of select="TSeq_taxid"/>
        <xsl:text>&#x9;</xsl:text>
        <xsl:value-of select="TSeq_sequence"/>
        <xsl:text>&#10;</xsl:text>
</xsl:for-each>
</xsl:template>
</xsl:stylesheet>
