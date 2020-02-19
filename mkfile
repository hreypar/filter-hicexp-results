# DESCRIPTION:
# mk module to filter significant regions from
# a normalized and compared hicexp object 
#
# USAGE:
# Single target execution: `mk <TARGET>` where TARGET is
# any line printed by the script `bin/mk-targets`
#
# Multiple target execution in tandem `bin/mk-targets | xargs mk`
#
# AUTHOR: HRG
#
# Run R script to produce hicexp objects with significant regions.
#
results/%.cycnorm.glm.significant.Rds:	data/%.cycnorm.glm.Rds
	mkdir -p `dirname $target`
	bin/filter-hicexp.R \
		--vanilla \
		$prereq \
		$target

