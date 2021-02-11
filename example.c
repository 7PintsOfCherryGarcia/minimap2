// To compile:
//   gcc -g -O2 example.c libminimap2.a -lz

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include "minimap.h"
#include "kseq.h"
#include "sqzlib.h"

//Changing zlib for sqzlib
//old zlib
//KSEQ_INIT(gzFile, gzread)
//new sqzlib
KSEQ_INIT(sqz_File*, sqz_sqzread)

int main(int argc, char *argv[])
{
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	int n_threads = 3;

	mm_verbose = 2; // disable message output to stderr
	mm_set_opt(0, &iopt, &mopt);
	mopt.flag |= MM_F_CIGAR; // perform alignment

	if (argc < 3) {
		fprintf(stderr, "Usage: minimap2-lite <target.fa> <query.fa>\n");
		return 1;
	}

	// open query file for reading; you may use your favorite FASTA/Q parser
  //Changing zlib for sqzlib
  //old zlib
  //gzFile f = gzopen(argv[2], "r");
  //new sqzlib
  sqz_File f = sqz_sqzopen(argv[2]);
	assert(f.sqz);
  assert(f.blk);
	kseq_t *ks = kseq_init(&f);
  fprintf(stderr, "PENE\n");
	// open index reader
	mm_idx_reader_t *r = mm_idx_reader_open(argv[1], &iopt, 0);
	mm_idx_t *mi;
  // traverse each part of the index
  while ((mi = mm_idx_reader_read(r, n_threads)) != 0) {
    fprintf(stderr, "PENE2\n");
    // this sets the maximum minimizer occurrence;
    //TODO: set a better default in mm_mapopt_init()!
    mm_mapopt_update(&mopt, mi);
    // thread buffer; for multi-threading, allocate one tbuf for each thread
		mm_tbuf_t *tbuf = mm_tbuf_init();
    // open query file for reading; you may use your favorite FASTA/Q parser
    //Changing zlib for sqzlib
    //old zlib
    //gzrewind(f);
    //new sqzlib
    sqzrewind(&f);
    kseq_rewind(ks);
    fprintf(stderr, "PENE3\n");
    // each kseq_read() call reads one query sequence
    while (kseq_read(ks) >= 0) {
      fprintf(stderr, "^^^%s\n", ks->name.s);
      mm_reg1_t *reg;
			int j, i, n_reg;
      // get all hits for the query
			reg = mm_map(mi, ks->seq.l, ks->seq.s, &n_reg, tbuf, &mopt, 0);
      // traverse hits and print them out
      for (j = 0; j < n_reg; ++j) {
        // with MM_F_CIGAR, this should not be NULL
        mm_reg1_t *r = &reg[j];
				assert(r->p);
				printf("%s\t%ld\t%d\t%d\t%c\t",
               ks->name.s, ks->seq.l, r->qs, r->qe, "+-"[r->rev]);
				printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tcg:Z:",
               mi->seq[r->rid].name, mi->seq[r->rid].len,
               r->rs, r->re, r->mlen, r->blen, r->mapq);
        // IMPORTANT: this gives the CIGAR in the aligned regions.
        // NO soft/hard clippings!
        for (i = 0; i < r->p->n_cigar; ++i)
					printf("%d%c", r->p->cigar[i]>>4, "MIDNSH"[r->p->cigar[i]&0xf]);
				putchar('\n');
				free(r->p);
			}
			free(reg);
		}
		mm_tbuf_destroy(tbuf);
		mm_idx_destroy(mi);
	}
  // close the index reader
	mm_idx_reader_close(r);
  // close the query file
	kseq_destroy(ks);
	sqz_sqzclose(f);
	return 0;
}
