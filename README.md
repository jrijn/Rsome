---
title: "Untitled"
author: "Jorik van Rijn"
date: "7/3/2019"
output:
  pdf_document: default
  word_document: default
  html_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## EXP-19-BT0308

The analysis of this qPCR experiment was performed using the Rsome package. This package was developed in house, more information is available in the vignette from http://www.github.com/jrijn/Rsome.

# QPCR analysis

First, the samples were labeled according to the plate setup in CFX maestro. The output of CFX was saved as .txt files, and imported in R for analysis.

The working directory were the .txt files are stored was defined:

```{r import}
#setwd('//user.uu.se/bmci/IMB-Users/jorva668/Documents/Lab_data_JR/2019/EXP-19-BT0308-differentiation human organoids')

filename <- "20190703_094313_CT029447_QPCR_JORIK -  Quantification Cq Results"
meltderivative <- "20190703_094313_CT029447_QPCR_JORIK -  Melt Curve Derivative Results_SYBR"
```

Next run the Rsome qPCR pipeline:

```{r message=FALSE}
df <- Rsome::cqimport(filename)
p1 <- Rsome::cq.plot(df)

mc <- Rsome::mcimport(cqimport = df, meltderivative = meltderivative)
p2 <- Rsome::mc.plot(mc)

re <- Rsome::relex(df, household = 'HP1BP3')
p3 <- Rsome::re.plot(re)

p1
p2
p3

ggsave(plot=p1, "output/cqplot.pdf", width=20, height=10)
ggsave(plot=p2, "output/mcplot.pdf", width=20, height=10)
ggsave(plot=p3, "output/relexplot.pdf", width=20, height=10)
```

