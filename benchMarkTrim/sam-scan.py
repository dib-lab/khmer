#! /usr/bin/env python2
import sys
import argparse
import screed
import math

def ignore_at(iter):
    for item in iter:
        if item.startswith('@'):
            continue
        yield item

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('genome')
    parser.add_argument('samfile')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)

    args = parser.parse_args()

    genome_dict = dict([ (record.name, record.sequence) for record in \
                        screed.open(args.genome) ])

    n = 0
    n_skipped = 0
    n_rev = n_fwd = 0

    for samline in ignore_at(open(args.samfile)):
        n += 1
        if n % 100000 == 0:
            print >>sys.stderr, '...', n

        readname, flags, refname, refpos, _, _, _, _, _, seq = \
                  samline.split('\t')[:10]
        if refname == '*' or refpos == '*':
            # (don't count these as skipped)
            continue
        
        refpos = int(refpos)
        try:
            ref = genome_dict[refname][refpos-1:refpos+len(seq) - 1]
        except KeyError:
            print >>sys.stderr, "unknown refname: %s; ignoring (read %s)" % (refname, readname)
            n_skipped += 1
            continue

        errors = []
        for pos, (a, b) in enumerate(zip(ref, seq)):
            if a.upper() != b.upper():
                # see http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output - '16' is the revcomp flag.
                if int(flags) & 16:
                    pos = len(seq) - pos - 1
                    n_rev += 1
                else:
                    n_fwd += 1
                errors.append(pos)

        print >>args.outfile, readname, ",".join(map(str, errors))

    # avoid log errors via pseudocount
    n_fwd += 1
    n_rev += 1
    
    print >>sys.stderr, 'logratio of fwd to rev: %.2f' % (math.log(n_fwd / float(n_rev), 2))
    if n_skipped / float(n) > .01:
        raise Exception, "Error: too many reads ignored! %d of %d" % \
              (n_skipped, n)

if __name__ == '__main__':
    main()
