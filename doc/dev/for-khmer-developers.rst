A guide for khmer committers
============================

This document is for people with commit rights to github.com/dib-lab/khmer.

----

If you have commit privileges to the dib-lab/khmer repository, here are a
few useful tips.

First, never merge something unless it's been through a review!  This
rule can be broken under specific conditions when doing a release; see
:doc:`release`.

Second, need to force another continuous integration run? Put "test
this please" in a comment.  This can be used to ask our continuous
integration system to run on someone else's pull request -- by
default, it only runs on commits from people who have write privileges
to khmer, so you may need to do this if you're reviewing someone else's
pull request.

Third, we ask that all contributors set up standing Pull Requests
while they are working something.  (This is a **requirement** if
you're in the GED lab.)  This lets us track what's going on. On the
flip side, please do not review pull requests until they are indicated
as "ready for review".
