/*
 * qf.h
 *
 * Copyright (c) 2014 Vedant Kumar <vsk@berkeley.edu>
 */

#pragma once

#include <stdint.h>
#include <stdbool.h>

struct quotient_filter {
	uint8_t qf_qbits;
	uint8_t qf_rbits;
	uint8_t qf_elem_bits;
	uint32_t qf_entries;
	uint64_t qf_index_mask;
	uint64_t qf_rmask;
	uint64_t qf_elem_mask;
	uint64_t qf_max_size;
	uint64_t *qf_table;
};

struct qf_iterator {
	uint64_t qfi_index;
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
 * Only the lowest q+r bits are actually inserted into the QF table.
 *
 * Returns false if the QF is full.
 */
bool qf_insert(struct quotient_filter *qf, uint64_t hash);

/*
 * Returns true if the QF may contain the hash. Returns false otherwise.
 */
bool qf_may_contain(struct quotient_filter *qf, uint64_t hash);

/*
 * Removes a hash from the QF.
 *
 * Caution: If you plan on using this function, make sure that your hash
 * function emits no more than q+r bits. Consider the following scenario;
 *
 *	insert(qf, A:X)   # X is in the lowest q+r bits.
 *	insert(qf, B:X)   # This is a no-op, since X is already in the table.
 *	remove(qf, A:X)   # X is removed from the table.
 *
 * Now, may-contain(qf, B:X) == false, which is a ruinous false negative.
 *
 * Returns false if the hash uses more than q+r bits.
 */
bool qf_remove(struct quotient_filter *qf, uint64_t hash);

/*
 * Initializes qfout and copies over all elements from qf1 and qf2.
 * Caution: qfout holds twice as many entries as either qf1 or qf2.
 *
 * Returns false on ENOMEM.
 */
bool qf_merge(struct quotient_filter *qf1, struct quotient_filter *qf2,
	struct quotient_filter *qfout);

/*
 * Resets the QF table. This function does not deallocate any memory.
 */
void qf_clear(struct quotient_filter *qf);

/*
 * Finds the size (in bytes) of a QF table.
 *
 * Caution: sizeof(struct quotient_filter) is not included.
 */
size_t qf_table_size(uint32_t q, uint32_t r);

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
 * Caution: Do not call this routine if qfi_done() == true.
 */
uint64_t qfi_next(struct quotient_filter *qf, struct qf_iterator *i);
