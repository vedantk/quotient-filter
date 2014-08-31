What are quotient filters?
==========================

A quotient filter is a compact structure which efficiently supports the
following operations [1]:

- Insert(qf, key)
- May-Contain(qf, key)
- Remove(qf, key)
- Merge(qf1, qf2) -> qfout

Like Bloom filters, quotient filters support approximate membership tests with
no false negatives. Advantages of quotient filters over Bloom filters include:

- Superior data locality
- Lessened reliance on hash functions
- Deterministic key removal
- The ability to merge filters without rehashing data

These properties make quotient filters useful as on-disk data structures. See
[2] for a specific example with the log-structured merge tree.

Advantages of Bloom filters over quotient filters include:

- Significantly lower implementation complexity
- Constant query cost w.r.t the load factor of the filter

Semantics of the QF metadata bits
=================================

- is_occupied

  The purpose of the `is_occupied' bit is to quickly let the search algorithm
  determine whether a fingerprint s.t Fq = x exists in the QF, where x is an
  index in the QF table. I.e, the `is_occupied' bit is attached to a slot in
  the table and not to any ephemeral fingerprint.

- is_continuation

  Set to false for the first fingerprint in a run.
  Set to true for subsequent fingerprints in a run.

- is_shifted

  Set to true when a fingerprint is not in its canonical slot (i.e when
  Fq != x, where x is the slot index). Set to false otherwise. Note that this
  metadata bit is affixed to the fingerprint, and not to a particular slot.

References
==========

[1] Don’t Thrash: How to Cache Your Hash on Flash (Bender et. al.)

[2] Building Workload-Independent Storage with VT-Trees (Shetty et. al.)
