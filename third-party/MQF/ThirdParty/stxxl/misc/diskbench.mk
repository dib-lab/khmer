############################################################################
#  misc/diskbench.mk
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2008-2011 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

HOST		?= unknown
FILE_SIZE	?= $(or $(SIZE),100)	# GiB
BLOCK_SIZE	?= $(or $(STEP),256)	# MiB
BATCH_SIZE	?= 1	# blocks
DIRECT_IO	?= yes	# unset to disable O_DIRECT
SYNC_IO		?= no   # set to 'yes' to enable O_SYNC|O_DSYNC|O_RSYNC
FILE_TYPE	?= syscall
FILE_RESIZE	?= no

disk2file	?= /stxxl/sd$1/stxxl

DISKS_1by1	?= a b c d
DISKS_a		?= a
DISKS_ab	?= a b
DISKS_abcd	?= a b c d

DISKS_a_ro	?= a
FLAGS_a_ro	?= R

DISKS		?= abcd $(DISKS_1by1) ab a_ro

SIP_NUM_BLOCKS	?= 0
SIP_CHUNK_BLOCKS?= 1
SIP_BLOCK_SIZE	?= 0

DISKBENCH_TITLE	?= STXXL Disk Benchmark $(DISKNAME) B=$(strip $(BATCH_SIZE))x$(call format_block_size,$(BLOCK_SIZE)) @ $(HOST)
DISKAVG_TITLE	?= STXXL Disk Benchmark $(DISKNAME) @ $(HOST)$(if $(DISKNAME2),\n$(DISKNAME2))

ifndef MISC_BINDIR
MISC_BINDIR	:= $(dir $(lastword $(MAKEFILE_LIST)))
endif
DISKBENCH_BINDIR?= $(MISC_BINDIR)/../io
DISKBENCH	?= benchmark_disks.stxxl.bin
SCATTERINPLACE	?= iobench_scatter_in_place.stxxl.bin

ECHO		?= echo
pipefail	?= set -o pipefail;

$(foreach d,$(DISKS_1by1),$(eval DISKS_$d ?= $d))

ifeq ($(SHELL),/bin/sh)
SHELL		 = bash
endif

DISKBENCH_FLAGS	+= $(if $(filter y yes Y YES,$(DIRECT_IO)),,--no-direct)
DISKBENCH_FLAGS	+= $(if $(filter y yes Y YES,$(SYNC_IO)),--sync)
DISKBENCH_FLAGS	+= --file-type=$(strip $(FILE_TYPE))
DISKBENCH_FLAGS	+= $(if $(filter y yes Y YES,$(FILE_RESIZE)),--resize)

define do-some-disks
	$(if $(filter ???,$(strip $(BLOCK_SIZE))),$(error ERROR: BLOCK_SIZE=$(strip $(BLOCK_SIZE))))
	-$(pipefail) \
	$(if $(IOSTAT_PLOT_RECORD_DATA),$(IOSTAT_PLOT_RECORD_DATA) -p $(@:.log=)) \
	$(DISKBENCH_BINDIR)/$(DISKBENCH) $(DISKBENCH_FLAGS) 0 $(strip $(FILE_SIZE)) $(strip $(BLOCK_SIZE)) $(strip $(BATCH_SIZE)) $(FLAGS_$*) $(FLAGS_EX) $(foreach d,$(DISKS_$*),$(call disk2file,$d)) | tee $@

endef

define do-some-disks-sip
	-$(pipefail) \
	$(if $(IOSTAT_PLOT_RECORD_DATA),$(IOSTAT_PLOT_RECORD_DATA) -p $(@:.log=)) \
	$(DISKBENCH_BINDIR)/$(SCATTERINPLACE) $(strip $(SIP_NUM_BLOCKS)) $(strip $(SIP_CHUNK_BLOCKS)) $(strip $(SIP_BLOCK_SIZE)) $(foreach d,$(DISKS_$*),$(call disk2file,$d)) | tee $@

endef

$(HOST)-%.cr.log:
	$(if $(keep-old-file),,$(RM) $(foreach d,$(DISKS_$*),$(call disk2file,$d)))
	$(do-some-disks)

$(HOST)-%.crx.log: FLAGS_EX = W
$(HOST)-%.crx.log:
	$(if $(keep-old-file),,$(RM) $(foreach d,$(DISKS_$*),$(call disk2file,$d)))
	$(do-some-disks)

# interleaved write-read-test
$(HOST)-%.wr.log:
	$(do-some-disks)

# scanning write
$(HOST)-%.wrx.log: FLAGS_EX = W
$(HOST)-%.wrx.log:
	$(do-some-disks)

# scanning read
$(HOST)-%.rdx.log: FLAGS_EX = R
$(HOST)-%.rdx.log:
	$(do-some-disks)

# scatter-in-place
$(HOST)-%.sip.log:
	$(do-some-disks-sip)

all: crx wr ex

crx: $(foreach d,$(DISKS_1by1),$(HOST)-$d.crx.log)
cr: $(foreach d,$(DISKS_1by1),$(HOST)-$d.cr.log)
cr+: $(foreach d,$(DISKS),$(HOST)-$d.crx.log)
wr: $(foreach d,$(DISKS),$(HOST)-$d.wr.log)
wrx: $(foreach d,$(DISKS_1by1),$(HOST)-$d.wrx.log)
wr+: $(foreach d,$(DISKS),$(HOST)-$d.wrx.log)
rdx: $(foreach d,$(DISKS_1by1),$(HOST)-$d.rdx.log)
rd+: $(foreach d,$(DISKS),$(HOST)-$d.rdx.log)
ex: $(foreach d,$(DISKS_1by1),$(HOST)-$d.wrx.log $(HOST)-$d.rdx.log)
ex+: $(foreach d,$(DISKS),$(HOST)-$d.wrx.log $(HOST)-$d.rdx.log)
sip: $(foreach d,$(DISKS_1by1),$(HOST)-$d.sip.log)

all-sizes-targets	?= cr+ wr ex+
all-sizes:
	for d in $(wildcard 0016iMB 0064MB 0004MB 0256MB 0001MB 1024MB *MB) ; do make -C $$d $(all-sizes-targets) ; done

all-sizes-raw: keep-old-file=1
all-sizes-raw: all-sizes

plot: $(HOST).gnuplot
	gnuplot $<

dotplot: $(HOST).d.gnuplot
	gnuplot $<

avgplot: $(HOST)-avg.gnuplot
	gnuplot $<

avg3plot: $(HOST)-avg3.gnuplot
	gnuplot $<

# $1 = logfile, $2 = column
extract_average	= $(if $(wildcard $1),$(shell grep ' Average over ' $1 | awk '{ print $$($2+1) }'),......)

# $1 = logfile, $2 = disk, $3 = column, $4 = label
# (does not plot if avg = nan)
define plotline
	$(if $(wildcard $1),$(if $(filter nan,$(call extract_average,$1,$3)),,$(ECHO) '        "$1" using ($$3/1024):($$$3) w l title "$2 $4 ($(call extract_average,$1,$3))", \' >> $@))

endef

# $1 = logfile, $2 = disk
define plotline-cr1
	$(call plotline,$1,$2,7,cr1)
endef
define plotline-cr
	$(call plotline,$1,$2,7,cr)
endef
define plotline-crx
	$(call plotline,$1,$2,7,crx)
endef
define plotline-wr
	$(call plotline,$1,$2,7,wr)
endef
define plotline-rd
	$(call plotline,$1,$2,14,rd)
endef
define plotline-wrx
	$(call plotline,$1,$2,7,wrx)
endef
define plotline-rdx
	$(call plotline,$1,$2,14,rdx)
endef

# $1 = disk letter
disk2label	?= sd$1
disks2label	?= sd[$1]

DISKNAME	?= unknown disk
DISKNAME2	?=# optional second line describing the test environment
PLOTXMAX	?= 475
PLOTYMAX	?= 120
AVGPLOTYMAX	?= $(PLOTYMAX)

fmt_block_size_2560000B		?= 2.5
fmt_block_size_12800000B	?= 12.5
fmt_block_size_51200000B	?= 50
format_block_size = $(or $(fmt_block_size_$(strip $1)),$(strip $1))MiB

$(HOST).gnuplot: $(MAKEFILE_LIST) $(wildcard *.log)
	$(RM) $@
	$(ECHO) 'set title "$(DISKBENCH_TITLE)"' >> $@
	$(ECHO) 'set xlabel "Disk offset [GiB]"' >> $@
	$(ECHO) 'set ylabel "Bandwidth per disk [MiB/s]"' >> $@
	$(ECHO) '' >> $@

	$(ECHO) 'plot [0:$(PLOTXMAX)] [0:$(PLOTYMAX)] \' >> $@
	$(foreach d,$(DISKS_1by1),\
		$(call plotline-cr1,$(HOST)-$d.cr1.log,$(call disk2label,$d)) \
		$(call plotline-crx,$(HOST)-$d.crx.log,$(call disk2label,$d)) \
		$(call plotline-cr,$(HOST)-$d.cr.log,$(call disk2label,$d)) \
		$(call plotline-wr,$(HOST)-$d.wr1.log,$(call disk2label,$d)) \
		$(call plotline-rd,$(HOST)-$d.wr1.log,$(call disk2label,$d)) \
		$(call plotline-wr,$(HOST)-$d.wr.log,$(call disk2label,$d)) \
		$(call plotline-rd,$(HOST)-$d.wr.log,$(call disk2label,$d)) \
		$(call plotline-wrx,$(HOST)-$d.wrx.log,$(call disk2label,$d)) \
		$(call plotline-rdx,$(HOST)-$d.rdx.log,$(call disk2label,$d)) \
	)
	$(foreach d,$(filter-out $(DISKS_1by1),$(DISKS)),\
		$(call plotline-crx,$(HOST)-$d.crx.log,$(call disks2label,$d)) \
		$(call plotline-wr,$(HOST)-$d.wr.log,$(call disks2label,$d)) \
		$(call plotline-rd,$(HOST)-$d.wr.log,$(call disks2label,$d)) \
		$(call plotline-wrx,$(HOST)-$d.wrx.log,$(call disks2label,$d)) \
		$(call plotline-rdx,$(HOST)-$d.rdx.log,$(call disks2label,$d)) \
	)
	$(ECHO) '        "nothing" notitle' >> $@

	$(ECHO) '' >> $@
	$(ECHO) 'pause -1' >> $@
	$(ECHO) '' >> $@
	$(ECHO) 'set title "$(call GNUPLOT_PS_STRING_ESCAPE,$(DISKBENCH_TITLE))"' >> $@
	$(ECHO) 'set term postscript enhanced $(GNUPLOT_PS_COLOR) 10' >> $@
	$(ECHO) 'set output "$(HOST).ps"' >> $@
	$(ECHO) '$(GNUPLOTFILEINFO)' >> $@
	$(ECHO) 'replot' >> $@

$(HOST).d.gnuplot: $(HOST).gnuplot
	sed -e 's/ w l / w d lw 2 /' $< > $@

$(HOST)-avg.dat: $(MISC_BINDIR)/diskbench-avgdat.sh $(wildcard *KB/*.log *MB/*.log)
	$(MISC_BINDIR)/diskbench-avgdat.sh $(wildcard *KB *MB) > $@

$(HOST)-avg.gnuplot: $(HOST)-avg.dat $(MAKEFILE_LIST)
	$(RM) $@
	$(ECHO) 'set title "$(DISKAVG_TITLE)"' >> $@
	$(ECHO) 'set xlabel "Block Size [MiB]"' >> $@
	$(ECHO) 'set ylabel "Average Sequential Bandwidth [MiB/s]"' >> $@
	$(ECHO) 'set key bottom' >> $@
	$(ECHO) '' >> $@

	$(ECHO) 'plot [] [0:$(AVGPLOTYMAX)] \' >> $@
	$(ECHO) '        "$(HOST)-avg.dat" using 0:2:xtic(1) w lp lt 1 pt 1 title "crx", \' >> $@
	$(ECHO) '        "$(HOST)-avg.dat" using 0:3:xtic(1) w lp lt 2 pt 2 title "wr", \' >> $@
	$(ECHO) '        "$(HOST)-avg.dat" using 0:4:xtic(1) w lp lt 3 pt 3 title "rd", \' >> $@
	$(ECHO) '        "$(HOST)-avg.dat" using 0:5:xtic(1) w lp lt 4 pt 4 title "wrx", \' >> $@
	$(ECHO) '        "$(HOST)-avg.dat" using 0:6:xtic(1) w lp lt 5 pt 5 title "rdx", \' >> $@
	$(ECHO) '        "nothing" notitle' >> $@

	$(ECHO) '' >> $@
	$(ECHO) 'pause -1' >> $@
	$(ECHO) '' >> $@
	$(ECHO) 'set term png size 800,600' >> $@
	$(ECHO) 'set output "$(HOST)-avg.png"' >> $@
	$(ECHO) 'replot' >> $@
	$(ECHO) '' >> $@
	$(ECHO) 'set title "$(call GNUPLOT_PS_STRING_ESCAPE,$(DISKAVG_TITLE))"' >> $@
	$(ECHO) 'set term postscript enhanced $(GNUPLOT_PS_COLOR)' >> $@
	$(ECHO) 'set output "$(HOST)-avg.ps"' >> $@
	$(ECHO) '$(GNUPLOTFILEINFO)' >> $@
	$(ECHO) 'replot' >> $@

$(HOST)-avg3.gnuplot: $(HOST)-avg.gnuplot
	grep -v -E '0:[34]:xtic' $< | sed -e 's/"crx"/"create"/g;s/"wrx"/"write"/g;s/"rdx"/"read"/g;/set output/s/-avg\./-avg3./g' > $@

-include iostat-plot.mk

