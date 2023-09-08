Store source code here that is not being used in EPIC,
but is too clever to throw away.

-----
Spring 2013

cdisort.nonscat.c
cdisort.nonscat.h

These were an attempt by T. Dowling to develop a faster
two-stream algorithm than c_twostr() in the case of no
scattering.  Alas, this code only beat the original for
the non-scattering, non-thermal emission case
(ds->flag.planck == FALSE), which became the EPIC
function beer_law_only(). Nevertheless, there is some
code here worth keeping.
-----
