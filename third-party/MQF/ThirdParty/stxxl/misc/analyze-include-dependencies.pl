#!/usr/bin/perl -w
############################################################################
#  misc/analyze-include-dependencies.pl
#
#  A tool to analyze #include dependencies and find unwanted cycles.
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2007-2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

use File::Basename;

$debug = 0;
$fakeheaders = 'fakeinclude';
# fake headers currently needed to prevent canonical name errors:
# backward_warning.h -> backward/backward_warning.h
# io.h -> #nothing
# platform.hpp -> boost/thread/detail/platform.hpp
# thread_data.hpp -> boost/thread/pthread/thread_data.hpp
# thread_heap_alloc.hpp -> boost/thread/detail/thread_heap_alloc.hpp
# vstring.tcc -> ext/vstring.tcc

$mcstl = 0;
$mcstlpath = 'c++';

$stxxl = 1;
$stxxlpath = 'include';

#$CXX = 'g++-4.2';
$CXX = 'g++-4.4 -std=c++0x';
$cxxtarget = 'foo-bar-gnu';
$cxxheaderpath = undef;
$gcctargetheaderpath = undef;

@includepath = ();
push @includepath, $fakeheaders if $fakeheaders;
push @includepath, $mcstlpath if $mcstl;
push @includepath, $stxxlpath if $stxxl;

%breakloops = qw(
	libio.h			_G_config.h
	sys/cdefs.h		features.h
	sys/ucontext.h		signal.h
	debug/formatter.h	debug/debug.h
	debug/bitset		bitset
	debug/deque		deque
	debug/list		list
	debug/map		map
	debug/set		set
	debug/vector		vector
	bits/sstream.tcc	sstream
	debug/hash_map		ext/hash_map
	debug/hash_set		ext/hash_set
	debug/unordered_map	unordered_map
	debug/unordered_set	unordered_set
	xmmintrin.h		emmintrin.h
	ext/vstring.h		vstring.tcc
	wctype.h		wchar.h
	stdatomic.h		cstdatomic
	algorithm		parallel/algorithm
	numeric			parallel/numeric
	bits/stl_algobase.h	parallel/algobase.h
	parallel/partition.h	parallel/sort.h
	boost/type_traits/is_class.hpp	boost/type_traits/is_scalar.hpp
	boost/type_traits/msvc/remove_cv.hpp	boost/type_traits/is_pointer.hpp
	boost/preprocessor/list/fold_left.hpp	boost/preprocessor/control/while.hpp
	boost/preprocessor/list/fold_right.hpp	boost/preprocessor/control/while.hpp
	);

%seen = ();
@todo = ();
%in = ();
%out = ();
$out{'MISSING'} = [];
%canonical = ();

sub get_file_list($;@)
{
	my $path = shift;
	my @patterns = @_;
	@patterns = ('*') unless scalar @patterns;

	my @l;
	foreach my $p (@patterns){
		foreach (glob "$path/$p") {
			if (-f $_) {
				my $x = $path;
				$x =~ s/([+])/\\$1/g;
				s|^$x/||;
				push @l, $_;
			}
		}
	}
	print "GLOB $path @patterns: @l\n" if $debug;
	return @l;
}

sub get_cxx_include_paths($)
{
	my $cxx = shift;
	open CXXOUT, "$cxx -E -v -xc++ /dev/null 2>&1 >/dev/null |" or die $!;
	my $ok = 0;
	while (<CXXOUT>) {
		chomp;
		$cxxtarget = $1 if /^Target: (.*)/;
		$ok = 1 if /^#include .\.\.\.. search starts here:$/;
		$ok = 0 if /^End of search list\.$/;
		if ($ok && s/^ //) {
			next if /backward/;
			push @includepath, $_;
			unless ($cxxheaderpath) {
				$cxxheaderpath = $_ if m|/c\+\+|;
			}
			unless ($gcctargetheaderpath) {
				$gcctargetheaderpath = $_ if (m|/$cxxtarget/| && ! m|/c\+\+|);
			}
		}
	}
	close CXXOUT;
	print "TARGET: \t$cxxtarget\n";
	print "HEADERS c++: \t$cxxheaderpath\n";
	print "HEADERS gcc: \t$gcctargetheaderpath\n";
	print "INCLUDES: \t@includepath\n";
}

sub find_header($;$)
{
	my $header = shift;
	my $relpath = dirname(shift || '.');
	foreach $_ (@includepath, $relpath) {
		my $file = "$_/$header";
		if (-f $file) {
			if (exists $canonical{$file} && $canonical{$file} ne $header) {
				print "CANONICAL MISMATCH: $file $header $canonical{$file}\n";
			}
			$canonical{$file} = $header;
			print "FOUND: $header as $file\n" if $debug;
			return [ $header, "$file" ];
		}
	}
	print "NOT FOUND: $header\n";
	return [$header, undef];
}

sub parse_header($)
{
	my $arg = shift;
	my $header = $$arg[0];
	my $file = $$arg[1];
	return if $seen{$file || $header};
	$seen{$file || $header} = 1;
	if ($file) {
		print "PROCESSING: $header \t($file)\n";
		print "OUCH: $header\n" if exists $$out{$header};
		$out{$header} = [];
		open HEADER,"<$file" or die $!;
		while ($_ = <HEADER>) {
			if (/^\s*#\s*include\s*["<]([^">]*)[">]/) {
				my $dep_header = $1;
				print "DEP: $header \t$dep_header\n";
				push @{$out{$header}}, $dep_header;
				push @{$in{$dep_header}}, $header;
				push @todo, find_header($dep_header, $file);
			}
		}
		close HEADER;
	} else {
		print "NOT FOUND: $header\n";
		push @{$out{$header}}, 'MISSING';
		push @{$in{'MISSING'}}, $header;
	}
}


get_cxx_include_paths($CXX);

@cxxheaders = get_file_list($cxxheaderpath);
#@cxxbwheaders = get_file_list("$cxxheaderpath/backward");
push @cxxbwheaders, get_file_list($cxxheaderpath, "backward/*");
@cxxextheaders = get_file_list($cxxheaderpath, 'ext/*');
@cxxbitsheaders = get_file_list($cxxheaderpath, 'bits/*');
@cxxtargetheaders = get_file_list("$cxxheaderpath/$cxxtarget", '*', 'bits/*');
#@cxxtr1headers = get_file_list($cxxheaderpath, 'tr1/*', 'tr1_impl/*');
@gcctargetheaders = get_file_list($gcctargetheaderpath, '*.h', 'bits/*.h');

if ($mcstl) {
@mcstlheaders = get_file_list($mcstlpath);
@mcstlmetaheaders = get_file_list($mcstlpath, 'meta/*');
@mcstlbitsheaders = get_file_list($mcstlpath, 'bits/*');

push @todo, find_header($_) foreach (
	sort(@mcstlheaders),
	sort(@mcstlmetaheaders),
	sort(@mcstlbitsheaders),
	);
}

if ($stxxl) {
@stxxlheaders = get_file_list($stxxlpath, '*', 'stxxl/*');
@stxxlbitsheaders = get_file_list($stxxlpath, 'bits/*', 'stxxl/bits/*', 'stxxl/bits/*/*');

push @todo, find_header($_) foreach (
	sort(@stxxlheaders),
	sort(@stxxlbitsheaders),
	);
}

push @todo, find_header($_) foreach (
	sort(@cxxheaders),
	sort(@cxxbwheaders),
	sort(@cxxextheaders),
	sort(@cxxbitsheaders),
	sort(@cxxtargetheaders),
	#sort(@cxxtr1headers),
	sort(@gcctargetheaders),
	);

while (@todo) {
	parse_header(shift @todo);
}

%odgr = ();
@zodgr = ();

if ($debug) {
foreach (sort keys %out) {
	print "\t$_\n\t=> (".(scalar @{$out{$_}}).")\t";
	foreach (@{$out{$_}}) {
		print " $_";
	}
	print "\n\t<= (".(scalar @{$in{$_}}).")\t";
	foreach (@{$in{$_}}) {
		print " $_";
	}
	print "\n";
}
}

foreach (sort keys %out) {
	$odgr{$_} = scalar @{$out{$_}};
	push @zodgr, $_ if $odgr{$_} == 0;
}

while (my ($h, $l) = each %breakloops) {
	next unless exists $out{$h};
	foreach (@{$out{$h}}) {
		if ($l eq $_) {
			print "BREAK LOOP: $h --> $l\n";
			--$odgr{$h};
			push @zodgr, $h if $odgr{$h} == 0;
			$_ .= ".UNLOOP";
		}
	}
	foreach (@{$in{$l}}) {
		if ($h eq $_) {
			$_ .= ".UNLOOP";
		}
	}
}

print "TOPSORT:\n";
while (@zodgr) {
	$curr = shift @zodgr;
	next unless exists $out{$curr};
	print "\t$curr";
	foreach (@{$in{$curr}}) {
		--$odgr{$_};
		push @zodgr, $_ if $odgr{$_} == 0;
		print " $_(".$odgr{$_}.")" if $debug;
	}
	print "\n";
	delete $out{$curr};
}

print "CYCLIC(".(scalar keys %out)."):\n";
foreach (sort keys %out) {
	print "\t$_($odgr{$_}):";
	foreach (@{$out{$_}}) {
		print " $_" if exists $out{$_};
	}
	print "\n";
}

