# find-nullomers

Find k-mers that have an abundance of 0 in the input
files. Canonicalizes k-mers so that forward and reverse complement
k-mers are correctly accounted for; this means that 'TTT' will count
as 'AAA', and only 'AAA' will be output if neither 'AAA' nor 'TTT'
exist.

To test, run:

```
./find-nullomers.py -k 3 tst.fa
```

and you should see the output

```
AAA
CCC
```
