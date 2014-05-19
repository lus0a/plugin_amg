/*
 * documentation.h
 *
 *  Created on: 03.04.2014
 *      Author: mrupp
 */

/*


A is additive
The Master/Slave Layout is with minmal overlap, that means only one vertex overlaps.

x  x  x  x  x  x  x  x  x  x  x

x  x  x  x  x  M
               S  x  x  x  x  x


How to construct P:

Construct the prolongation for all master or inner nodes.
If that one is finished, P is additive unique then.

OVERLAPs

Note that you will need mostly a Master Overlap of 1 for this:

This is a master overlap 1 (this is what RSAMG needs)

x  x  x  x  x  M0 S1
               S0 M1  x  x  x  x


This is a slave overlap 1

x  x  x  x  M1 M0
            S1 S0 x  x  x  x  x

This is a slave AND master overlap 1

x  x  x  x  M1 M0 S1
            S1 S0 M0 x  x  x  x


This is a slave 1 master overlap 2
(needed by FAMG, since it needs 2-neighborhood of every inner/master node)

x  x  x  x  M1 M0 S1 S2
            S1 S0 M0 M1 x  x  x

General rule: if you need a N-Neighborhood of every inner/master node, you need a Master N, Slave N-1 overlap.

ADDITIVE to CONSISTENT

After that, we need to make P consistent, since we need to calculate A_H = R A P, and since A ist additive, P and R need to be consistent.

How to make P consistent













 */
