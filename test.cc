/*
 * test.cc
 *
 * Copyright (c) 2014 Vedant Kumar <vsk@berkeley.edu>
 */

extern "C" {
  #include "qf.c"
}

#include <set>
#include <vector>
#include <cassert>
#include <cstdio>

using namespace std;

static void fail(struct quotient_filter *qf, const char *s)
{
  fprintf(stderr, "qf(q=%u, r=%u): %s\n", qf->qf_qbits, qf->qf_rbits, s);
  abort();
}

static uint64_t rand64()
{
  return (((uint64_t) rand()) << 32) | ((uint64_t) rand());
}

static void qf_print(struct quotient_filter *qf)
{
  char buf[32];
  uint32_t pad = qf->qf_qbits + 2;

  for (uint32_t i = 0; i < pad; ++i) {
    printf(i == qf->qf_qbits ? "|" : " ");
  }
  puts("is_shifted | is_continuation | is_occupied | remainder");

  for (uint64_t idx = 0; idx < qf->qf_max_size; ++idx) {
    snprintf(buf, sizeof(buf), "%llu", idx);
    printf("%s", buf);

    int fillspace = qf->qf_qbits - strlen(buf);
    for (int i = 0; i < fillspace; ++i) {
      printf(" ");
    }
    printf("| ");

    uint64_t elt = get_elem(qf, idx);
    printf("%d          | ", !!is_shifted(elt));
    printf("%d               | ", !!is_continuation(elt));
    printf("%d           | ", !!is_occupied(elt));
    printf("%llu\n", get_remainder(elt));
  }
}

/* Check QF structural invariants. */
static void qf_consistent(struct quotient_filter *qf)
{
  assert(qf->qf_qbits);
  assert(qf->qf_rbits);
  assert(qf->qf_qbits + qf->qf_rbits <= 64);
  assert(qf->qf_elem_bits == (qf->qf_rbits + 3));

  uint64_t idx;
  uint64_t start;
  uint64_t size = qf->qf_max_size;
  assert(qf->qf_entries <= size);
  uint64_t last_run_elt;

  for (start = 0; start < size; ++start) {
    if (is_cluster_start(get_elem(qf, start))) {
      break;
    }
  }

  idx = start;
  do {
    uint64_t elt = get_elem(qf, idx);

    /* Check for invalid metadata bits. */
    if (is_continuation(elt)) {
      assert(is_shifted(elt));

      /* Check that this is actually a continuation. */
      uint64_t prev = get_elem(qf, decr(qf, idx));
      assert(!is_empty(prev));
    }

    /* Check that remainders within runs are sorted. */
    if (!is_empty(elt)) {
      uint64_t rem = get_remainder(elt);
      if (is_continuation(elt)) {
        assert(rem > last_run_elt);
      }
      last_run_elt = rem;
    }

    idx = incr(qf, idx);
  } while (idx != start);
}

/* Generate a random 64-bit hash. If @clrhigh, clear the high (64-p) bits. */
static uint64_t genhash(struct quotient_filter *qf, bool clrhigh,
    set<uint64_t> &keys)
{
  uint64_t hash;
  uint64_t mask = clrhigh ? LOW_MASK(qf->qf_qbits + qf->qf_rbits) : ~0UL;
  uint64_t size = qf->qf_max_size;

  /* If the QF is overloaded, use a linear scan to find an unused hash. */
  if (keys.size() > (3 * (size / 4))) {
    uint64_t probe;
    uint64_t start = rand64() & qf->qf_index_mask;
    for (probe = incr(qf, start); probe != start; probe = incr(qf, probe)) {
      uint64_t hi = clrhigh ? 0 : (rand64() & mask);
      hash = hi | probe;
      if (!keys.count(hash)) {
        return hash;
      }
    }
  }

  /* Find a random unused hash. */
  do {
    hash = rand64() & mask;
  } while (keys.count(hash));
  return hash;
}

/* Insert a random p-bit hash into the QF. */
static void ht_put(struct quotient_filter *qf, set<uint64_t> &keys)
{
  uint64_t hash = genhash(qf, true, keys);
  assert(qf_insert(qf, hash));
  keys.insert(hash);
}

/* Decode the p-bit hash at QF[idx] and remove it from the table. */
static uint64_t ht_erase(struct quotient_filter *qf, uint64_t idx)
{
  uint64_t c;
  uint64_t run_delta = 0;
  uint64_t elt = get_elem(qf, idx);
  for (c = idx; !is_cluster_start(get_elem(qf, c)); c = decr(qf, c)) {
    if (is_run_start(get_elem(qf, c))) {
      ++run_delta;
    }
  }

  uint64_t hi = c + run_delta;
  uint64_t lo = get_remainder(elt);
  uint64_t hash = (hi << qf->qf_rbits) | lo;
  qf_remove(qf, hash);
  assert(!qf_may_contain(qf, hash));
  return hash;
}

/* Remove a hash from the filter. */
static void ht_del(struct quotient_filter *qf, set<uint64_t> &keys)
{
  uint64_t idx;
  uint64_t start = ((uint64_t) rand()) & qf->qf_index_mask;
  for (idx = incr(qf, start); idx != start; idx = incr(qf, idx)) {
    uint64_t elt = get_elem(qf, idx);
    if (!elt) {
      continue;
    }

    uint64_t hash = ht_erase(qf, idx);
    keys.erase(hash);
    return;
  }
}

/* Check that a set of keys are in the QF. */
static void ht_check(struct quotient_filter *qf, set<uint64_t> &keys)
{
  qf_consistent(qf);
  set<uint64_t>::iterator it;
  for (it = keys.begin(); it != keys.end(); ++it) {
    uint64_t hash = *it;
    assert(qf_may_contain(qf, hash));
  }
}

static void qf_test(struct quotient_filter *qf)
{
  /* Basic get/set tests. */
  uint64_t idx;
  uint64_t size = qf->qf_max_size;
  for (idx = 0; idx < size; ++idx) {
    assert(get_elem(qf, idx) == 0);
    set_elem(qf, idx, idx);
  }
  for (idx = 0; idx < size; ++idx) {
    assert(get_elem(qf, idx) == idx);
  }
  qf_clear(qf);

  /* Random get/set tests. */
  vector<uint64_t> elements(size, 0);
  for (idx = 0; idx < size; ++idx) {
    uint64_t slot = rand64() % size;
    uint64_t hash = rand64();
    set_elem(qf, slot, hash);
    elements[slot] = hash & qf->qf_elem_mask;
  }
  for (idx = 0; idx < elements.size(); ++idx) {
    assert(get_elem(qf, idx) == elements[idx]);
  }
  qf_clear(qf);

  /* Check: forall x, insert(x) => may-contain(x). */
  set<uint64_t> keys;
  for (idx = 0; idx < size; ++idx) {
    uint64_t elt = genhash(qf, false, keys);
    assert(qf_insert(qf, elt));

    // XXX
    qf_print(qf);

    keys.insert(elt);
  }
  ht_check(qf, keys);
  keys.clear();
  qf_clear(qf);

  // Enable once qf_remove() has been implemented.
#if 0
  /* Check that the QF works like a hash set when all keys are p-bit values. */
  for (idx = 0; idx < 1000000; ++idx) {
    while (qf->qf_entries < size) {
      ht_put(qf, keys);
    }
    while (qf->qf_entries > (size / 2)) {
      ht_del(qf, keys);
    }
    ht_check(qf, keys);

    struct qf_iterator qfi;
    qfi_start(qf, &qfi);
    while (!qfi_done(qf, &qfi)) {
      uint64_t hash = qfi_next(qf, &qfi);
      assert(keys.count(hash));
    }
  }
#endif
}

/* Fill up the QF (at least partially). */
static void random_fill(struct quotient_filter *qf)
{
  set<uint64_t> keys;
  uint64_t elts = ((uint64_t) rand()) % qf->qf_max_size;
  while (elts) {
    ht_put(qf, keys);
    --elts;
  }
  qf_consistent(qf);
}

/* Check if @lhs is a subset of @rhs. */
static void subsetof(struct quotient_filter *lhs, struct quotient_filter *rhs)
{
  struct qf_iterator qfi;
  qfi_start(lhs, &qfi);
  while (!qfi_done(lhs, &qfi)) {
    uint64_t hash = qfi_next(lhs, &qfi);
    assert(qf_may_contain(rhs, hash));
  }
}

int main()
{
  srand(0);
  for (uint32_t q = 1; q < 24; ++q) {
    for (uint32_t r = 1; r < 32; ++r) {
      struct quotient_filter qf;
      if (!qf_init(&qf, q, r)) {
        fail(&qf, "init-1");
      }
      qf_test(&qf);
      qf_destroy(&qf);
    }
  }

  for (uint32_t q1 = 1; q1 < 24; ++q1) {
    for (uint32_t r1 = 1; r1 < 32; ++r1) {
      for (uint32_t q2 = 1; q2 < 24; ++q2) {
        for (uint32_t r2 = 1; r2 < 32; ++r2) {
          if ((q1 + r1) != (q2 + r2)) {
            continue;
          }

          struct quotient_filter qf;
          struct quotient_filter qf1, qf2;
          if (!qf_init(&qf1, q1, r1) || !qf_init(&qf2, q2, r2)) {
            fail(&qf1, "init-2");
          }

          random_fill(&qf1);
          random_fill(&qf2);
          assert(qf_merge(&qf1, &qf2, &qf));
          qf_consistent(&qf);
          subsetof(&qf1, &qf);
          subsetof(&qf2, &qf);
          qf_destroy(&qf1);
          qf_destroy(&qf2);
          qf_destroy(&qf);
        }
      }
    }
  }

  puts("[PASSED] basic qf tests");
  return 0;
}
