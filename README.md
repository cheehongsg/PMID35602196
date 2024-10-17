# Code used in publication
"Reduced subgenomic RNA expression is a molecular indicator of asymptomatic SARS-CoV-2 infection" ([PMID: 35602196](https://github.com/cheehongsg/PMID35602196))


NOTE: Code package provided as-is without any warranty.


# Usage examples

Download this repository. Perl script uses samtools and bedtools and assume that they are accessible in your Unix environment.

Let's assume that you have place the STAR aligned bam file for sample1 (```sample1.bam```) and sample2 (```sample2.bam```) in the subfolder ```mydata```.

We will process sample1.

**IMPORTANT:** it is ASSUME that your **current working directory** is the subfolder ```mydata```

```bash
export BAMPREFIX=sample1

# step 1
# input = STAR alignment file = <bamprefix>.bam
# output = <bamprefix>.amplicontag.bam
perl ../scripts/artictools.pl splitreads ${BAMPREFIX}.bam
# output:
# Processing sample1.sarcov2.bam.. 0.01M.. 0.02M.. 0.03M.. 0.04M.. 0.05M.. 0.06M.. 0.07M.. 0.08M.. 0.09M.. 0.10M.. 0.11M.. 0.12M.. 0.13M.. 0.14M.. 0.15M.. 0.16M.. 0.17M.. 0.18M.. 0.19M.. 0.20M.. 0.21M.. 0.22M.. 0.23M.. 0.24M.. 0.25M.. 0.26M.. 0.27M.. 0.28M.. 0.29M.. 0.30M.. 0.31M.. 0.32M.. 0.33M.. 0.34M.. 0.35M.. 0.36M.. 0.37M.. 0.38M.. 0.39M.. 0.40M.. 0.41M.. 0.42M.. 0.43M.. 0.44M.. 0.45M.. 0.46M.. 0.47M.. 0.48M.. 0.49M.. 0.50M.. 0.51M.. 0.52M.. 0.53M.. 0.54M.. 0.55M.. 0.56M.. 0.57M.. 0.58M.. 0.59M.. 0.60M.. 0.61M.. 0.62M.. 0.63M.. 0.64M.. 0.65M.. 0.66M.. 0.67M.. 0.68M.. 0.69M.. 0.70M.. 0.71M.. 0.72M.. 0.73M.. 0.74M.. 0.75M.. 758797.. done!
#

# step 2
# input = <bamprefix>.amplicontag.bam from step 1
# output = <bamprefix>.amplicontag.xls.gz
# output = <bamprefix>.amplicontag.oddpool.bam & .bai
# output = <bamprefix>.amplicontag.evenpool.bam & .bai
# output = <bamprefix>.amplicontag.mixedpool.bam & .bai
# output = <bamprefix>.amplicontag.primerchimeric.bam & .bai
perl ../scripts/artictools.pl splitpool ${BAMPREFIX}.amplicontag.bam
# output:
# Processing sample1.sarcov2.amplicontag.bam.. 0.01M.. 0.02M.. 0.03M.. 0.04M.. 0.05M.. 0.06M.. 0.07M.. 0.08M.. 0.09M.. 0.10M.. 0.11M.. 0.12M.. 0.13M.. 0.14M.. 0.15M.. 0.16M.. 0.17M.. 0.18M.. 0.19M.. 0.20M.. 0.21M.. 0.22M.. 0.23M.. 0.24M.. 0.25M.. 0.26M.. 0.27M.. 0.28M.. 0.29M.. 0.30M.. 0.31M.. 0.32M.. 0.33M.. 0.34M.. 0.35M.. 0.36M.. 0.37M.. 0.38M.. 0.39M.. 0.40M.. 0.41M.. 0.42M.. 0.43M.. 0.44M.. 0.45M.. 0.46M.. 0.47M.. 0.48M.. 0.49M.. 0.50M.. 0.51M.. 0.52M.. 0.53M.. 0.54M.. 0.55M.. 0.56M.. 0.57M.. 0.58M.. 0.59M.. 0.60M.. 0.61M.. 0.62M.. 0.63M.. 0.64M.. 0.65M.. 0.66M.. 0.67M.. 0.68M.. 0.69M.. 0.70M.. 0.71M.. 715340.. (2022 mapq<255) done!
# Writing sample1.sarcov2.amplicontag.xls.gz.. done!
# Processing sample1.sarcov2.amplicontag.bam.. 0.01M.. 0.02M.. 0.03M.. 0.04M.. 0.05M.. 0.06M.. 0.07M.. 0.08M.. 0.09M.. 0.10M.. 0.11M.. 0.12M.. 0.13M.. 0.14M.. 0.15M.. 0.16M.. 0.17M.. 0.18M.. 0.19M.. 0.20M.. 0.21M.. 0.22M.. 0.23M.. 0.24M.. 0.25M.. 0.26M.. 0.27M.. 0.28M.. 0.29M.. 0.30M.. 0.31M.. 0.32M.. 0.33M.. 0.34M.. 0.35M.. 0.36M.. 0.37M.. 0.38M.. 0.39M.. 0.40M.. 0.41M.. 0.42M.. 0.43M.. 0.44M.. 0.45M.. 0.46M.. 0.47M.. 0.48M.. 0.49M.. 0.50M.. 0.51M.. 0.52M.. 0.53M.. 0.54M.. 0.55M.. 0.56M.. 0.57M.. 0.58M.. 0.59M.. 0.60M.. 0.61M.. 0.62M.. 0.63M.. 0.64M.. 0.65M.. 0.66M.. 0.67M.. 0.68M.. 0.69M.. 0.70M.. 0.71M.. 715340.. odd(346074,1000 mapq<255), even(270942,864 mapq<255), mixed(46460,48 mapq<255), pcd(1574,2 mapq<255), ignore(50290,108 mapq<255) done!
# Indexing sample1.sarcov2.amplicontag.oddpool.bam.. done!
# Indexing sample1.sarcov2.amplicontag.evenpool.bam.. done!
# Indexing sample1.sarcov2.amplicontag.mixedpool.bam.. done!
# Indexing sample1.sarcov2.amplicontag.primerchimeric.bam..[W::bam_hdr_read] EOF marker is absent. The input is probably truncated
#  done!
#

# step 3
# input = <bamprefix>.amplicontag.bam
# internal input = <bamprefix>.amplicontag.oddpool.bam & .bai
# internal input = <bamprefix>.amplicontag.evenpool.bam & .bai
# internal input = <bamprefix>.amplicontag.mixedpool.bam & .bai
# output = <bamprefix>.amplicontag.oddpool.bdg
# output = <bamprefix>.amplicontag.evenpool.bdg
# output = <bamprefix>.amplicontag.mixedpool.bdg
perl ../scripts/artictools.pl splitpoolcov ${BAMPREFIX}.amplicontag.bam
# output:
# Processing sample1.sarcov2.amplicontag.oddpool.bam.. 0.01M.. 0.02M.. 0.03M.. 0.04M.. 0.05M.. 0.06M.. 0.07M.. 0.08M.. 0.09M.. 0.10M.. 0.11M.. 0.12M.. 0.13M.. 0.14M.. 0.15M.. 0.16M.. 0.17M.. 0.18M.. 0.19M.. 0.20M.. 0.21M.. 0.22M.. 0.23M.. 0.24M.. 0.25M.. 0.26M.. 0.27M.. 0.28M.. 0.29M.. 0.30M.. 0.31M.. 0.32M.. 0.33M.. 0.34M.. 346074.. done!
# Accumulating reads.. 0.01M.. 0.02M.. 0.03M.. 0.04M.. 0.05M.. 0.06M.. 0.07M.. 0.08M.. 0.09M.. 0.10M.. 0.11M.. 0.12M.. 0.13M.. 0.14M.. 0.15M.. 0.16M.. 0.17M.. 173037.. done!
# Processing sample1.sarcov2.amplicontag.oddpool.bdg.. done!
# Processing sample1.sarcov2.amplicontag.evenpool.bam.. 0.01M.. 0.02M.. 0.03M.. 0.04M.. 0.05M.. 0.06M.. 0.07M.. 0.08M.. 0.09M.. 0.10M.. 0.11M.. 0.12M.. 0.13M.. 0.14M.. 0.15M.. 0.16M.. 0.17M.. 0.18M.. 0.19M.. 0.20M.. 0.21M.. 0.22M.. 0.23M.. 0.24M.. 0.25M.. 0.26M.. 0.27M.. 270942.. done!
# Accumulating reads.. 0.01M.. 0.02M.. 0.03M.. 0.04M.. 0.05M.. 0.06M.. 0.07M.. 0.08M.. 0.09M.. 0.10M.. 0.11M.. 0.12M.. 0.13M.. 135471.. done!
# Processing sample1.sarcov2.amplicontag.evenpool.bdg.. done!
# Processing sample1.sarcov2.amplicontag.mixedpool.bam.. 0.01M.. 0.02M.. 0.03M.. 0.04M.. 46460.. done!
# Accumulating reads.. 0.01M.. 0.02M.. 23230.. done!
# Processing sample1.sarcov2.amplicontag.mixedpool.bdg.. done!
#

# step 4
# input = <bamprefix>.amplicontag.bam
# internal input = <bamprefix>.amplicontag.oddpool.bdg
# internal input  = <bamprefix>.amplicontag.evenpool.bdg
# output = <bamprefix>.amplicontag.coverage.statistics
# otuput = <bamprefix>.amplicontag.coverage.junctions
perl ../scripts/artictools.pl tabulatecov ${BAMPREFIX}.amplicontag.bam
# output:
# Processing sample1.sarcov2.amplicontag.oddpool.bam.. 0.01M.. 0.02M.. 0.03M.. 0.04M.. 0.05M.. 0.06M.. 0.07M.. 0.08M.. 0.09M.. 0.10M.. 0.11M.. 0.12M.. 0.13M.. 0.14M.. 0.15M.. 0.16M.. 0.17M.. 0.18M.. 0.19M.. 0.20M.. 0.21M.. 0.22M.. 0.23M.. 0.24M.. 0.25M.. 0.26M.. 0.27M.. 0.28M.. 0.29M.. 0.30M.. 0.31M.. 0.32M.. 0.33M.. 0.34M.. 346074.. done!
# Processing sample1.sarcov2.amplicontag.evenpool.bam.. 0.01M.. 0.02M.. 0.03M.. 0.04M.. 0.05M.. 0.06M.. 0.07M.. 0.08M.. 0.09M.. 0.10M.. 0.11M.. 0.12M.. 0.13M.. 0.14M.. 0.15M.. 0.16M.. 0.17M.. 0.18M.. 0.19M.. 0.20M.. 0.21M.. 0.22M.. 0.23M.. 0.24M.. 0.25M.. 0.26M.. 0.27M.. 270942.. done!
#

# step 5
# <minReadSupport> <minSampleSupport> <reportAuxResult=0> <reportPrefix> [<.amplicontag.bam1> [<.amplicontag.bam2> ...]]
# input = <bamprefix>.amplicontag.bam
# output = <reportPrefix>.grid.absolute.expression.ge<minReadSupport>readsupport.ge<minSampleSupport>sample.xls
perl ../scripts/artictools.pl gridexpress 1 1 0 oneSample ${BAMPREFIX}.amplicontag.bam
# output:
# 


```

We may process sample2 by repeating the same process after the following setting:
```bash
export BAMPREFIX=sample2
```
and repeat the bash code above (step 1-5).

To create a grid from both sample1 and sample2:
```bash
# input = sample1.amplicontag.bam and sample2.amplicontag.bam
# output = twoSamples.grid.absolute.expression.ge<minReadSupport>readsupport.ge<minSampleSupport>sample.xls
perl ../scripts/artictools.pl gridexpress 1 1 0 twoSamples sample1.amplicontag.bam sample2.amplicontag.bam
```

To get the translation products:
```bash
# translation
tail -n +7 twoSamples.grid.absolute.expression.ge1readsupport.ge1sample.xls \
| cut -f 1-2 \
| perl ../scripts/translate.pl \
1> twoSamples.grid.absolute.expression.ge1readsupport.ge1sample.canonical.product.xls \
2> twoSamples.grid.absolute.expression.ge1readsupport.ge1sample.canonical.product.log 
# output:
#
```

# Others
[Long-Read Workshop 2021 Slides](https://cheehongsg.github.io/web/LRW2021/SAR-CoV-2.html)

