############################################################################
#  misc/iostat-plot.mk
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2008-2011 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

#
# to record the data, include this file from your Makefile
# and use a target like
#
# my_results.out:
#	$(IOSTAT_PLOT_RECORD_DATA) -p $(@:.out=) \
#	my_program arg1 ... argn > $@
#
# and then run
#
# make my_results.out
# make my_results.io.plot
# make my_results.cpu.plot
#

empty	?=#
space	?= $(empty) $(empty)
comma	?= ,

ifndef IOSTAT_PLOT_BINDIR
IOSTAT_PLOT_BINDIR		:= $(dir $(lastword $(MAKEFILE_LIST)))
endif
IOSTAT_PLOT_RECORD_DATA		?= $(IOSTAT_PLOT_BINDIR)/record-load-iostat
IOSTAT_PLOT_CONCAT_LINES	?= $(IOSTAT_PLOT_BINDIR)/concat-lines
IOSTAT_PLOT_FLOATING_AVERAGE	?= $(IOSTAT_PLOT_BINDIR)/floating-average

IOSTAT_PLOT_LINE_IDENTIFIER	?= ^sd
IOSTAT_PLOT_LINE_IGNORE_PATTERN	?=
IOSTAT_PLOT_CPU_LINE_IDENTIFIER	?= ^$(space)$(space)$(space)
IOSTAT_PLOT_DISKS		?= 8
IOSTAT_PLOT_CPUS		?= 8
IOSTAT_PLOT_AVERAGE		?= 1
IOSTAT_PLOT_AVERAGE.io		?= $(IOSTAT_PLOT_AVERAGE)
IOSTAT_PLOT_AVERAGE.cpu		?= $(IOSTAT_PLOT_AVERAGE)
IOSTAT_PLOT_LOAD_REDUCE		?= 0
IOSTAT_PLOT_DISK_LIST		?= sda sdb sdc sdd sde sdf sdg sdh sdi sdj
IOSTAT_PLOT_IO_WITH_UTILIZATION	?= no

IOSTAT_PLOT_Y_LABEL.io		?= Bandwidth [MiB/s]
IOSTAT_PLOT_Y_LABEL.cpu		?= CPU Usage [%]

GNUPLOT_PS_COLOR		?= color solid
GNUPLOT_PS_STRING_ESCAPE	?= $(subst _,\\_,$(subst @,\\@,$(subst ~,\\~,$1)))
GNUPLOTFILEINFO			?= set label "$(call GNUPLOT_PS_STRING_ESCAPE,$(HOST):$(patsubst $(HOME)/%,~/%,$(CURDIR))/)" at character 0,-1 font ",6"

ECHO				?= echo

# $1 = iostat file
# (12)     Device:         rrqm/s   wrqm/s     r/s     w/s    rMB/s    wMB/s avgrq-sz avgqu-sz   await  svctm  %util
# with -p: Device:         rrqm/s   wrqm/s     r/s     w/s    rMB/s    wMB/s avgrq-sz avgqu-sz   await r_await w_await  svctm  %util
define get-deviceline
__deviceline	:= $(shell grep ^Device $1 | head -n 1)
endef
offset_read	 = $(if $(filter rMB/s,$(word 6,$(__deviceline))),6,8)
offset_write	 = $(if $(filter wMB/s,$(word 7,$(__deviceline))),7,9)
offset_util	 = $(if $(filter %util,$(word 12,$(__deviceline))),12,14)

# $1 = first, $2 = increment, $3 = last
define gnuplot-column-sequence-sum
$(subst $(space),+,$(patsubst %,$$%,$(shell seq $1 $2 $3)))
endef

# $1 = io | cpu, $2 = numdisks (io) | numcpus (cpu)
define template-iostat-gnuplot
	$(RM) $@
	$(ECHO) 'set title "$(subst _, ,$*) (avg=$(IOSTAT_PLOT_AVERAGE.$(strip $1)))"' >> $@
	$(ECHO) 'set xlabel "Time [s]"' >> $@
	$(ECHO) 'set ylabel "$(IOSTAT_PLOT_Y_LABEL.$(strip $1))"' >> $@
$(if $(and $(filter io,$1),$(filter yes,$(IOSTAT_PLOT_IO_WITH_UTILIZATION))),
	$(ECHO) 'set y2label "Utilization [%]"' >> $@
	$(ECHO) 'set ytics nomirror' >> $@
	$(ECHO) 'set y2tics' >> $@
,$(if $(wildcard $*.waitlog),
	$(ECHO) 'set y2label "Wait Time [s]"' >> $@
	$(ECHO) 'set ytics nomirror' >> $@
	$(ECHO) 'set y2tics' >> $@
))
#	$(ECHO) 'set style data linespoints' >> $@
	$(ECHO) 'set style data lines' >> $@
	$(ECHO) 'set macros' >> $@
	$(ECHO) 'set pointsize 0.4' >> $@
#	$(ECHO) 'set samples 1000' >> $@
	$(ECHO) 'set key top left' >> $@
	$(ECHO) 'set yrange [0:]' >> $@
	$(ECHO) '' >> $@
$(if $(filter io,$1),
	$(ECHO) 'read = "$(call gnuplot-column-sequence-sum,2,4,$(shell expr $2 '*' 4))"' >> $@
	$(ECHO) 'write = "$(call gnuplot-column-sequence-sum,3,4,$(shell expr $2 '*' 4))"' >> $@
	$(ECHO) 'utilize = "($(call gnuplot-column-sequence-sum,4,4,$(shell expr $2 '*' 4))) / $2"' >> $@
	$(ECHO) '' >> $@
	$(ECHO) 'plot \' >> $@
	$(ECHO) '	"$<" using 0:(@read + @write) title "Read + Write" ls 3$(comma) \' >> $@
	$(ECHO) '	"$<" using 0:(@read) title "Read" ls 2$(comma) \' >> $@
	$(ECHO) '	"$<" using 0:(@write) title "Write" ls 1$(comma) \' >> $@
$(if $(filter yes,$(IOSTAT_PLOT_IO_WITH_UTILIZATION)),
	$(ECHO) '	"$<" using 0:(@utilize) title "Utilization" ls 5 axes x1y2$(comma) \' >> $@
,$(if $(wildcard $*.waitlog),
	$(ECHO) '	"$*.waitlog" using 1:4 title "Wait Read" ls 5 axes x1y2$(comma) \' >> $@
	$(ECHO) '	"$*.waitlog" using 1:5 title "Wait Write" ls 4 axes x1y2$(comma) \' >> $@
))
	$(ECHO) '	"not.existing.dummy" using 8:15 notitle' >> $@
)
$(if $(filter cpu,$1),
	$(ECHO) 'plot \' >> $@
	$(ECHO) '	"$(word 2,$^)" using 0:(100*($$1-$(IOSTAT_PLOT_LOAD_REDUCE))/$(strip $2)) title "Load (100% = $2)" ls 5$(comma) \' >> $@
	$(ECHO) '	"$<" using 0:($$1+$$2+$$3+$$4) title "Total" ls 4$(comma) \' >> $@
	$(ECHO) '	"$<" using 0:2 title "Nice" ls 6$(comma) \' >> $@
	$(ECHO) '	"$<" using 0:4 title "Wait" ls 2$(comma) \' >> $@
	$(ECHO) '	"$<" using 0:1 title "User" ls 1$(comma) \' >> $@
	$(ECHO) '	"$<" using 0:3 title "System" ls 3$(comma) \' >> $@
	$(ECHO) '	"$(word 2,$^)" using 0:(100*($$1-$(IOSTAT_PLOT_LOAD_REDUCE))/$(strip $2)) notitle ls 5$(comma) \' >> $@
	$(ECHO) '	"not.existing.dummy" using 8:15 notitle' >> $@ 
)
	$(ECHO) '' >> $@
	$(ECHO) 'pause -1' >> $@
	$(ECHO) '' >> $@
	$(ECHO) 'set terminal postscript enhanced $(GNUPLOT_PS_COLOR) 10' >> $@
	$(ECHO) 'set output "$*.$(strip $1).eps"' >> $@
	$(ECHO) '$(GNUPLOTFILEINFO)' >> $@
	$(ECHO) 'replot' >> $@
	$(ECHO) '' >> $@
endef

define iostat-to-dat
	$(eval $(call get-deviceline,$<))
	$(pipefail) \
	cat $< | \
	$(if $(IOSTAT_PLOT_LINE_IGNORE_PATTERN),grep -v '$(IOSTAT_PLOT_LINE_IGNORE_PATTERN)' |) \
	awk '{ print $$1 " " $$$(offset_read) " " $$$(offset_write) " " $$$(offset_util) }' | \
	$(IOSTAT_PLOT_CONCAT_LINES) '$(IOSTAT_PLOT_LINE_IDENTIFIER)' | \
	grep '$(IOSTAT_PLOT_LINE_IDENTIFIER)' | \
	tail -n +2 | \
	$(IOSTAT_PLOT_FLOATING_AVERAGE) $(IOSTAT_PLOT_AVERAGE.io) > $@

endef

%.io-$(IOSTAT_PLOT_AVERAGE.io).dat: %.iostat $(MAKEFILE_LIST)
	$(iostat-to-dat)

define per-disk-plot-template
%.$(disk).io-$$(IOSTAT_PLOT_AVERAGE.io).dat: IOSTAT_PLOT_LINE_IDENTIFIER=^$(disk)
%.$(disk).io-$$(IOSTAT_PLOT_AVERAGE.io).plot: IOSTAT_PLOT_DISKS=1
%.$(disk).io-$$(IOSTAT_PLOT_AVERAGE.io).plot: IOSTAT_PLOT_IO_WITH_UTILIZATION=yes

%.$(disk).io-$$(IOSTAT_PLOT_AVERAGE.io).dat: %.iostat $$(MAKEFILE_LIST)
	$$(iostat-to-dat)
endef
$(foreach disk, $(IOSTAT_PLOT_DISK_LIST),$(eval $(per-disk-plot-template)))

%.cpu-$(IOSTAT_PLOT_AVERAGE.cpu).dat: %.iostat $(MAKEFILE_LIST)
	$(pipefail) \
	grep "$(IOSTAT_PLOT_CPU_LINE_IDENTIFIER)" $< | \
	tail -n +2 | \
	$(IOSTAT_PLOT_FLOATING_AVERAGE) $(IOSTAT_PLOT_AVERAGE.cpu) > $@

%.io-$(IOSTAT_PLOT_AVERAGE.io).plot: %.io-$(IOSTAT_PLOT_AVERAGE.io).dat $(MAKEFILE_LIST)
	$(call template-iostat-gnuplot,io,$(IOSTAT_PLOT_DISKS))

%.io.plot: %.io-$(IOSTAT_PLOT_AVERAGE.io).plot
	@$(ECHO) Your plot file is: $<

%.cpu-$(IOSTAT_PLOT_AVERAGE.cpu).plot: %.cpu-$(IOSTAT_PLOT_AVERAGE.cpu).dat %.loadavg $(MAKEFILE_LIST)
	$(call template-iostat-gnuplot,cpu,$(IOSTAT_PLOT_CPUS))

%.cpu.plot: %.cpu-$(IOSTAT_PLOT_AVERAGE.cpu).plot
	@$(ECHO) Your plot file is: $<

%.io.xplot: %.io-$(IOSTAT_PLOT_AVERAGE.io).plot
	gnuplot $<

%.cpu.xplot: %.cpu-$(IOSTAT_PLOT_AVERAGE.cpu).plot
	gnuplot $<

.SECONDARY:
.DELETE_ON_ERROR:
