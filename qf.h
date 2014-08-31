/*
 * qf.h
 *
 * Copyright (c) 2014 Vedant Kumar <vsk@berkeley.edu>
 */

#pragma once

#include <stdint.h>
#include <stdbool.h>

struct quotient_filter {
	uint32_t qf_qbits;
	uint32_t qf_rbits;
	uint32_t qf_elem_bits;
	uint32_t qf_entries;
	uint64_t qf_index_mask;
	uint64_t qf_rmask;
	uint64_t qf_elem_mask;
	uint64_t qf_max_size;
	uint64_t *qf_table;
};

struct qf_iterator {
	uint64_t qfi_index;
	uint64_t qfi_start_cluster;
	uint64_t qfi_quotient;
	uint64_t qfi_visited;
};

/*
 * Initializes a quotient filter with capacity 2^q.
 * Increasing r improves the filter's accuracy but uses more space.
 * 
 * Returns false if q == 0, r == 0, q+r > 64, or on ENOMEM.
 */
bool qf_init(struct quotient_filter *qf, uint32_t q, uint32_t r);

/*
 * Inserts a hash into the QF.
 *
 * Returns false if the QF is full.
 */
bool qf_insert(struct quotient_filter *qf, uint64_t hash);

/*
 * Returns true if the QF may contain the hash.
 * Returns false otherwise.
 */
bool qf_may_contain(struct quotient_filter *qf, uint64_t hash);

/*
 * Removes the hash from the QF.
 *
 * Caution: If your hash function emits more than q+r bits, some queries may
 * report false negatives.
 */
void qf_remove(struct quotient_filter *qf, uint64_t hash);

/*
 * Initializes qfout and copies over all elements from qf1 and qf2.
 * Caution: qfout is allocated twice as much space as either qf1 or qf2.
 *
 * Returns false if q1+r1 != q2+r2, or on ENOMEM.
 */
bool qf_merge(struct quotient_filter *qf1, struct quotient_filter *qf2,
	struct quotient_filter *qfout);

/*
 * Resets the QF table.
 *
 * This function does not deallocate any memory.
 */
void qf_clear(struct quotient_filter *qf);

/*
 * Finds the size (in bytes) of a QF table.
 *
 * Caution: sizeof(struct quotient_filter) is not included.
 */
size_t qf_table_size(uint32_t q, uint32_t r);

/*
 * Finds the expected false positive rate for the QF.
 *
 * The rate is in the interval [0, 1). See qf_init for ways to lower it.
 */
float qf_false_positive_rate(struct quotient_filter *qf);

/*
 * Deallocates the QF table.
 */
void qf_destroy(struct quotient_filter *qf);

/*
 * Initialize an iterator for the QF.
 */
void qfi_start(struct quotient_filter *qf, struct qf_iterator *i);

/*
 * Returns true if there are no elements left to visit.
 */
bool qfi_done(struct quotient_filter *qf, struct qf_iterator *i);

/*
 * Returns the next (q+r)-bit fingerprint in the QF.
 *
 * Caution: Do not call this routine is qfi_done() == true.
 */
uint64_t qfi_next(struct quotient_filter *qf, struct qf_iterator *i);
