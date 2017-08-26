/*
 * qf.c
 *
 * Copyright (c) 2014 Vedant Kumar <vsk@berkeley.edu>
 */

#include <stdlib.h>
#include <string.h>

#include "qf.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define LOW_MASK(n) ((1ULL << (n)) - 1ULL)

bool qf_init(struct quotient_filter *qf, uint32_t q, uint32_t r)
{
	if (q == 0 || r == 0 || q + r > 64) {
		return false;
	}

	qf->qf_qbits = q;
	qf->qf_rbits = r;
	qf->qf_elem_bits = qf->qf_rbits + 3;
	qf->qf_index_mask = LOW_MASK(q);
	qf->qf_rmask = LOW_MASK(r);
	qf->qf_elem_mask = LOW_MASK(qf->qf_elem_bits);
	qf->qf_entries = 0; 
	qf->qf_max_size = 1 << q;
	qf->qf_table = (uint64_t *) calloc(qf_table_size(q, r), 1);
	return qf->qf_table != NULL;
}

/* Return QF[idx] in the lower bits. */
static uint64_t get_elem(struct quotient_filter *qf, uint64_t idx)
{
	uint64_t elt = 0;
	size_t bitpos = qf->qf_elem_bits * idx;
	size_t tabpos = bitpos / 64;
	size_t slotpos = bitpos % 64;
	int spillbits = (slotpos + qf->qf_elem_bits) - 64;
	elt = (qf->qf_table[tabpos] >> slotpos) & qf->qf_elem_mask;
	if (spillbits > 0) {
		++tabpos;
		uint64_t x = qf->qf_table[tabpos] & LOW_MASK(spillbits);
		elt |= x << (qf->qf_elem_bits - spillbits);
	}
	return elt;
}

/* Store the lower bits of elt into QF[idx]. */
static void set_elem(struct quotient_filter *qf, uint64_t idx, uint64_t elt)
{
	size_t bitpos = qf->qf_elem_bits * idx;
	size_t tabpos = bitpos / 64;
	size_t slotpos = bitpos % 64;
	int spillbits = (slotpos + qf->qf_elem_bits) - 64;
	elt &= qf->qf_elem_mask;
	qf->qf_table[tabpos] &= ~(qf->qf_elem_mask << slotpos);
	qf->qf_table[tabpos] |= elt << slotpos;
	if (spillbits > 0) {
		++tabpos;
		qf->qf_table[tabpos] &= ~LOW_MASK(spillbits);
		qf->qf_table[tabpos] |= elt >> (qf->qf_elem_bits - spillbits);
	}
}

static inline uint64_t incr(struct quotient_filter *qf, uint64_t idx)
{
	return (idx + 1) & qf->qf_index_mask;
}

static inline uint64_t decr(struct quotient_filter *qf, uint64_t idx)
{
	return (idx - 1) & qf->qf_index_mask;
}

static inline int is_occupied(uint64_t elt)
{
	return elt & 1;
}

static inline uint64_t set_occupied(uint64_t elt)
{
	return elt | 1;
}

static inline uint64_t clr_occupied(uint64_t elt)
{
	return elt & ~1;
}

static inline int is_continuation(uint64_t elt)
{
	return elt & 2;
}

static inline uint64_t set_continuation(uint64_t elt)
{
	return elt | 2;
}

static inline uint64_t clr_continuation(uint64_t elt)
{
	return elt & ~2;
}

static inline int is_shifted(uint64_t elt)
{
	return elt & 4;
}

static inline uint64_t set_shifted(uint64_t elt)
{
	return elt | 4;
}

static inline uint64_t clr_shifted(uint64_t elt)
{
	return elt & ~4;
}

static inline uint64_t get_remainder(uint64_t elt)
{
	return elt >> 3;
}

static inline bool is_empty_element(uint64_t elt)
{
	return (elt & 7) == 0;
}

static inline bool is_cluster_start(uint64_t elt)
{
	return is_occupied(elt) && !is_continuation(elt) && !is_shifted(elt);
}

static inline bool is_run_start(uint64_t elt)
{
	return !is_continuation(elt) && (is_occupied(elt) || is_shifted(elt));
}

static inline uint64_t hash_to_quotient(struct quotient_filter *qf,
		uint64_t hash)
{
	return (hash >> qf->qf_rbits) & qf->qf_index_mask;
}

static inline uint64_t hash_to_remainder(struct quotient_filter *qf,
		uint64_t hash)
{
	return hash & qf->qf_rmask;
}

/* Find the start index of the run for fq (given that the run exists). */
static uint64_t find_run_index(struct quotient_filter *qf, uint64_t fq)
{
	/* Find the start of the cluster. */
	uint64_t b = fq;
	while (is_shifted(get_elem(qf, b))) {
		b = decr(qf, b);
	}

	/* Find the start of the run for fq. */
	uint64_t s = b;
	while (b != fq) {
		do {
			s = incr(qf, s);
		} while (is_continuation(get_elem(qf, s)));

		do {
			b = incr(qf, b);
		} while (!is_occupied(get_elem(qf, b)));
	}
	return s;
}

/* Insert elt into QF[s], shifting over elements as necessary. */
static void insert_into(struct quotient_filter *qf, uint64_t s, uint64_t elt)
{
	uint64_t prev;
	uint64_t curr = elt;
	bool empty;

	do {
		prev = get_elem(qf, s);
		empty = is_empty_element(prev);
		if (!empty) {
			/* Fix up `is_shifted' and `is_occupied'. */
			prev = set_shifted(prev);
			if (is_occupied(prev)) {
				curr = set_occupied(curr);
				prev = clr_occupied(prev);
			}
		}
		set_elem(qf, s, curr);
		curr = prev;
		s = incr(qf, s);
	} while (!empty);
}

bool qf_insert(struct quotient_filter *qf, uint64_t hash)
{
	if (qf->qf_entries >= qf->qf_max_size) {
		return false;
	}

	uint64_t fq = hash_to_quotient(qf, hash);
	uint64_t fr = hash_to_remainder(qf, hash);
	uint64_t T_fq = get_elem(qf, fq);
	uint64_t entry = (fr << 3) & ~7;

	/* Special-case filling canonical slots to simplify insert_into(). */
	if (is_empty_element(T_fq)) {
		set_elem(qf, fq, set_occupied(entry));
		++qf->qf_entries;
		return true;
	}

	if (!is_occupied(T_fq)) {
		set_elem(qf, fq, set_occupied(T_fq));
	}

	uint64_t start = find_run_index(qf, fq);
	uint64_t s = start;

	if (is_occupied(T_fq)) {
		/* Move the cursor to the insert position in the fq run. */
		do {
			uint64_t rem = get_remainder(get_elem(qf, s));
			if (rem == fr) {
				return true;
			} else if (rem > fr) {
				break;
			}
			s = incr(qf, s);
		} while (is_continuation(get_elem(qf, s)));

		if (s == start) {
			/* The old start-of-run becomes a continuation. */
			uint64_t old_head = get_elem(qf, start);
			set_elem(qf, start, set_continuation(old_head));
		} else {
			/* The new element becomes a continuation. */
			entry = set_continuation(entry);
		}
	}

	/* Set the shifted bit if we can't use the canonical slot. */
	if (s != fq) {
		entry = set_shifted(entry);
	}

	insert_into(qf, s, entry);
	++qf->qf_entries;
	return true;
}

bool qf_may_contain(struct quotient_filter *qf, uint64_t hash)
{
	uint64_t fq = hash_to_quotient(qf, hash);
	uint64_t fr = hash_to_remainder(qf, hash);
	uint64_t T_fq = get_elem(qf, fq);

	/* If this quotient has no run, give up. */
	if (!is_occupied(T_fq)) {
		return false;
	}

	/* Scan the sorted run for the target remainder. */
	uint64_t s = find_run_index(qf, fq);
	do {
		uint64_t rem = get_remainder(get_elem(qf, s));
		if (rem == fr) {
			return true;
		} else if (rem > fr) {
			return false;
		}
		s = incr(qf, s);
	} while (is_continuation(get_elem(qf, s)));
	return false;
}

/* Remove the entry in QF[s] and slide the rest of the cluster forward. */
static void delete_entry(struct quotient_filter *qf, uint64_t s, uint64_t quot)
{
	uint64_t next;
	uint64_t curr = get_elem(qf, s);
	uint64_t sp = incr(qf, s);
	uint64_t orig = s;

	/*
	 * FIXME(vsk): This loop looks ugly. Rewrite.
	 */
	while (true) {
		next = get_elem(qf, sp);
		bool curr_occupied = is_occupied(curr);

		if (is_empty_element(next) || is_cluster_start(next) || sp == orig) {
			set_elem(qf, s, 0);
			return;
		} else {
			/* Fix entries which slide into canonical slots. */
			uint64_t updated_next = next;
			if (is_run_start(next)) {
				do {
					quot = incr(qf, quot);
				} while (!is_occupied(get_elem(qf, quot)));

				if (curr_occupied && quot == s) {
					updated_next = clr_shifted(next);
				}
			}

			set_elem(qf, s, curr_occupied ?
					set_occupied(updated_next) :
					clr_occupied(updated_next));
			s = sp;
			sp = incr(qf, sp);
			curr = next;
		}
	}
}

bool qf_remove(struct quotient_filter *qf, uint64_t hash)
{
	uint64_t highbits = hash >> (qf->qf_qbits + qf->qf_rbits);
	if (highbits) {
		return false;
	}

	uint64_t fq = hash_to_quotient(qf, hash);
	uint64_t fr = hash_to_remainder(qf, hash);
	uint64_t T_fq = get_elem(qf, fq);

	if (!is_occupied(T_fq) || !qf->qf_entries) {
		return true;
	}

	uint64_t start = find_run_index(qf, fq);
	uint64_t s = start;
	uint64_t rem;

	/* Find the offending table index (or give up). */
	do {
		rem = get_remainder(get_elem(qf, s));
		if (rem == fr) {
			break;
		} else if (rem > fr) {
			return true;
		}
		s = incr(qf, s);
	} while (is_continuation(get_elem(qf, s)));
	if (rem != fr) {
		return true;
	}

	uint64_t kill = (s == fq) ? T_fq : get_elem(qf, s);
	bool replace_run_start = is_run_start(kill);

	/* If we're deleting the last entry in a run, clear `is_occupied'. */
	if (is_run_start(kill)) {
		uint64_t next = get_elem(qf, incr(qf, s));
		if (!is_continuation(next)) {
			T_fq = clr_occupied(T_fq);
			set_elem(qf, fq, T_fq);
		}
	}

	delete_entry(qf, s, fq);

	if (replace_run_start) {
		uint64_t next = get_elem(qf, s);
		uint64_t updated_next = next;
		if (is_continuation(next)) {
			/* The new start-of-run is no longer a continuation. */
			updated_next = clr_continuation(next);
		}
		if (s == fq && is_run_start(updated_next)) {
			/* The new start-of-run is in the canonical slot. */
			updated_next = clr_shifted(updated_next);
		}
		if (updated_next != next) {
			set_elem(qf, s, updated_next);
		}
	}

	--qf->qf_entries;
	return true;
}

bool qf_merge(struct quotient_filter *qf1, struct quotient_filter *qf2,
		struct quotient_filter *qfout)
{
	uint32_t q = 1 + MAX(qf1->qf_qbits, qf2->qf_qbits);
	uint32_t r = MAX(qf1->qf_rbits, qf2->qf_rbits);
	if (!qf_init(qfout, q, r)) {
		return false;
	}

	struct qf_iterator qfi;
	qfi_start(qf1, &qfi);
	while (!qfi_done(qf1, &qfi)) {
		qf_insert(qfout, qfi_next(qf1, &qfi));
	}
	qfi_start(qf2, &qfi);
	while (!qfi_done(qf2, &qfi)) {
		qf_insert(qfout, qfi_next(qf2, &qfi));
	}
	return true;
}

void qf_clear(struct quotient_filter *qf)
{
	qf->qf_entries = 0;
	memset(qf->qf_table, 0, qf_table_size(qf->qf_qbits, qf->qf_rbits));
}

size_t qf_table_size(uint32_t q, uint32_t r)
{
	size_t bits = (1 << q) * (r + 3);
	size_t bytes = bits / 8;
	return (bits % 8) ? (bytes + 1) : bytes;
}

void qf_destroy(struct quotient_filter *qf)
{
	free(qf->qf_table);
}

void qfi_start(struct quotient_filter *qf, struct qf_iterator *i)
{
	/* Mark the iterator as done. */
	i->qfi_visited = qf->qf_entries;

	if (qf->qf_entries == 0) {
		return;
	}

	/* Find the start of a cluster. */
	uint64_t start;
	for (start = 0; start < qf->qf_max_size; ++start) {
		if (is_cluster_start(get_elem(qf, start))) {
			break;
		}
	}

	i->qfi_visited = 0;
	i->qfi_index = start;
}

bool qfi_done(struct quotient_filter *qf, struct qf_iterator *i)
{
	return qf->qf_entries == i->qfi_visited;
}

uint64_t qfi_next(struct quotient_filter *qf, struct qf_iterator *i)
{
	while (!qfi_done(qf, i)) {
		uint64_t elt = get_elem(qf, i->qfi_index);

		/* Keep track of the current run. */
		if (is_cluster_start(elt)) {
			i->qfi_quotient = i->qfi_index;
		} else {
			if (is_run_start(elt)) {
				uint64_t quot = i->qfi_quotient;
				do {
					quot = incr(qf, quot);
				} while (!is_occupied(get_elem(qf, quot)));
				i->qfi_quotient = quot;
			}
		}

		i->qfi_index = incr(qf, i->qfi_index);

		if (!is_empty_element(elt)) {
			uint64_t quot = i->qfi_quotient;
			uint64_t rem = get_remainder(elt);
			uint64_t hash = (quot << qf->qf_rbits) | rem;
			++i->qfi_visited;
			return hash;
		}
	}

	abort();
}
