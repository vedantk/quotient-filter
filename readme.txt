What is this?
=============

This is a working in-memory quotient filter written in C.

qf.c: Implementation
qf.h: API and documentation
test.cc: Comprehensive randomized tester

What are quotient filters?
==========================

A quotient filter is a compact structure which efficiently supports the
following operations [1]:

- Insert(qf, key)
- May-Contain(qf, key)
- Remove(qf, key) (with a caveat, see the documentation in qf.h)
- Merge(qf1, qf2) -> qfout
- Iterate(qf)

Like Bloom filters, quotient filters support approximate membership tests with
no false negatives. Advantages of quotient filters over Bloom filters include:

- Improved locality of accesses
- Lessened reliance on hash functions
- Deterministic (or `correct') key removal
- The ability to merge filters without rehashing data
- The ability to iterate through hashes in the filter (modulo 2^(q+r))

These properties make quotient filters useful on-disk data structures. See [2]
for a specific example involving the log-structured merge tree.

Advantages of Bloom filters over quotient filters include:

- Significantly lower implementation complexity
- Constant query cost w.r.t the load factor of the filter

Example usage
=============

  struct quotient_filter qf;
  qf_init(&qf, 16, 8); // 2^16 slots, 8-bit remainders, ~0.7 MB.

  uint64_t hash = hash24();
  if (qf_insert(&qf, hash)) { // Fails only if the QF is already full.
    assert(qf_may_contain(&qf, hash));
  }

  qf_remove(&qf, hash);

Semantics of the QF metadata bits
=================================

This is a point of confusion which [3] attempts to clarify. To summarize:

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
  metadata bit is attached to the fingerprint, and not to a particular slot.

References
==========

[1] Don’t Thrash: How to Cache Your Hash on Flash (Bender et. al.)

[2] Building Workload-Independent Storage with VT-Trees (Shetty et. al.)

[3] Wikipedia, https://en.wikipedia.org/wiki/Quotient_filter
