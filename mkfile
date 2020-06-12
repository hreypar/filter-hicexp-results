# DESCRIPTION:
# mk module to filter hicexp objects and
# obtain significant region pairs
#
# USAGE:
# Single target execution: `mk <TARGET>` where TARGET is
# any line printed by the script `bin/mk-targets`
#
# Multiple target execution in tandem `bin/mk-targets | xargs mk`
#
# AUTHOR: HRG
#
< config.mk
# Run R script to filter significant pairs.
#
results/%.significantpairs.Rds:	data/%.qlf.cycnorm.hicexp.Rds
	mkdir -p `dirname $target`
	bin/filter-hicexp-list.R \
		--input $prereq \
		--pvalue $PVAL \
		--output $target

