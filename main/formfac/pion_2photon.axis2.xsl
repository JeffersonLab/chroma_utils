<?xml version="1.0"?>
<xsl:stylesheet version = '1.0'
     xmlns:xsl='http://www.w3.org/1999/XSL/Transform'>

<!-- Text output -->
<xsl:output method="text"/>

<!-- Filter rule -->
<xsl:template match="/">

<xsl:text>
#c \cr cs 0.5 e
</xsl:text>

<!-- Loop over ff -->
<xsl:for-each select="/Pion2Photon/VaryUnknowns/elem/FormFactors/elem">
  <h2>! Qsq = </h2><xsl:value-of select="Qsq"/> 
  <xsl:text>
  </xsl:text>

  <xsl:variable name="n" select="n"/>
  <h2>#m </h2><xsl:value-of select="$n+1"/>
  <xsl:text>
  </xsl:text>

  <xsl:for-each select="Fit/elem">
    <xsl:variable name="t_i" select="t_i"/>
    <xsl:for-each select="elem[number(n)=1]">
          <xsl:value-of select="$t_i"/> <xsl:text> </xsl:text>
          <xsl:value-of select="F"/> <xsl:text> </xsl:text>
          <xsl:value-of select="F_err"/>
          <xsl:text>
          </xsl:text>
    </xsl:for-each> <!-- End foreach -->
  </xsl:for-each> <!-- End foreach -->
</xsl:for-each> <!-- End foreach -->
</xsl:template>

</xsl:stylesheet> 
