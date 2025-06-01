#include <zlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <unistd.h>
#include "cgranges.h"
#include "ketopt.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

#define BEDTK_VERSION "1.0-r32-dirty"

/*****************
 * Faster printf *
 *****************/

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

static void mm_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[16]; // for integer to string conversion
	const char *p, *q;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) str_copy(s, q, p);
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				str_enlarge(s, l);
				for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				str_copy(s, r, r + strlen(r));
			} else if (*p == 'c') {
				str_enlarge(s, 1);
				s->s[s->l++] = va_arg(ap, int);
			} else abort();
			q = p + 1;
		}
	}
	if (p > q) str_copy(s, q, p);
	va_end(ap);
	s->s[s->l] = 0;
}

/***************
 * BED3 parser *
 ***************/

typedef struct {
	int32_t l;
	char *s;
} bed_rest1_t;

typedef struct {
	int64_t n, m;
	bed_rest1_t *a;
} bed_rest_t;

static char *parse_vcf(char *s, int32_t *st_, int32_t *en_, char **r) // r points to the end of contig name
{
	char *p, *q, *ctg = 0;
	int32_t i, st = -1, en = -1;
	*r = 0;
	if (s[0] == '#') {
		*r = s;
		return 0;
	}
	for (i = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
			if (i == 0) {
				ctg = p;
				*r = q + 1;
			} else if (i == 1) st = atol(p) - 1;
			else if (i == 3) en = st + (q - p);
			else if (i == 7) {
				char *s = 0;
				if (strncmp(p, "END=", 4) == 0) s = p + 1;
				else {
					s = strstr(p, ";END=");
					if (s) s += 5;
				}
				en = strtol(s, &s, 10);
			}
			if (i != 0) *q = c;
			++i, p = q + 1;
			if (i == 7 || c == '\0') break;
		}
	}
	*st_ = st, *en_ = en;
	return i >= 7? ctg : 0;
}

static char *parse_paf(char *s, int32_t *st_, int32_t *en_, char **r)
{
	char *p, *q, *ctg = 0;
	int32_t i, st = -1, en = -1;
	*r = 0;
	for (i = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
			if (i == 5) {
				ctg = p;
				*r = q + 1;
			} else if (i == 7) {
				st = atol(p);
			} else if (i == 8) {
				en = atol(p);
			}
			if (i != 5) *q = c;
			++i, p = q + 1;
			if (i == 9 || c == '\0') break;
		}
	}
	*st_ = st, *en_ = en;
	return i >= 9? ctg : 0;
}

static char *parse_bed3b(char *s, int32_t *st_, int32_t *en_, char **r)
{
	char *p, *q, *ctg = 0;
	int32_t i, st = -1, en = -1;
	if (r) *r = 0;
	for (i = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
			if (i == 0) ctg = p;
			else if (i == 1) st = atol(p);
			else if (i == 2) {
				en = atol(p);
				if (r && c != 0) *r = q, *q = c;
			}
			++i, p = q + 1;
			if (i == 3 || c == '\0') break;
		}
	}
	*st_ = st, *en_ = en;
	return i >= 3? ctg : 0;
}

static char *parse_bed3(char *s, int32_t *st_, int32_t *en_)
{
	return parse_bed3b(s, st_, en_, 0);
}

static cgranges_t *read_bed3b(const char *fn, bed_rest_t *r, const char *fn_order)
{
	gzFile fp;
	cgranges_t *cr;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int64_t k = 0;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) {
		fprintf(stderr, "ERROR: failed to open the input file\n");
		return 0;
	}
	cr = cr_init();
	if (fn_order) {
		gzFile fp;
		if ((fp = gzopen(fn_order, "r")) == 0) {
			fprintf(stderr, "ERROR: failed to open the list file\n");
			return 0;
		}
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
			char *p;
			for (p = str.s; *p && !isspace(*p); ++p);
			*p = 0;
			cr_add_ctg(cr, str.s, 0);
		}
		ks_destroy(ks);
		gzclose(fp);
	}
	ks = ks_init(fp);
	if (r) r->m = r->n = 0, r->a = 0;
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		char *ctg, *rest;
		int32_t st, en;
		ctg = parse_bed3b(str.s, &st, &en, &rest);
		if (ctg) {
			cr_add(cr, ctg, st, en, k);
			if (r) {
				bed_rest1_t *p;
				if (r->n == r->m) {
					int64_t old_m = r->m;
					r->m = r->m? r->m<<1 : 16;
					r->a = (bed_rest1_t*)realloc(r->a, r->m * sizeof(bed_rest1_t));
					memset(&r->a[old_m], 0, (r->m - old_m) * sizeof(bed_rest1_t));
				}
				p = &r->a[r->n++];
				p->l = rest? str.l - (rest - str.s) : 0;
				if (rest) {
					p->s = (char*)malloc(p->l + 1);
					memcpy(p->s, rest, p->l);
					p->s[p->l] = 0;
				}
			}
			++k;
		}
	}
	if (k > INT32_MAX)
		fprintf(stderr, "WARNING: more than %d records; some functionality may not work!\n", INT32_MAX);
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);
	return cr;
}

static cgranges_t *read_bed3(const char *fn)
{
	return read_bed3b(fn, 0, 0);
}

/*********
 * Tools *
 *********/

int main_isec(int argc, char *argv[])
{
	cgranges_t *cr, *qr;
	ketopt_t o = KETOPT_INIT;
	int64_t i, m_b = 0, *b = 0, n_b;
	int c;
	char *fn_order = 0;

	while ((c = ketopt(&o, argc, argv, 1, "s:", 0)) >= 0) {
		if (c == 's') fn_order = o.arg;
	}

	if (argc - o.ind < 1 || (argc - o.ind < 2 && isatty(0))) {
		printf("Usage: bedtk isec [options] <A.bed> <B.bed>\n");
		printf("Options:\n");
		printf("  -s FILE   list of contig IDs to specify the output order []\n");
		return 1;
	}

	cr = read_bed3b(argv[o.ind], 0, fn_order);
	assert(cr);
	cr_index2(cr, 1);

	qr = read_bed3b(o.ind + 1 < argc? argv[o.ind + 1] : 0, 0, fn_order);
	assert(qr);
	if (!cr_is_sorted(qr)) cr_sort(qr);
	cr_merge_pre_index(qr);
	for (i = 0; i < qr->n_r; ++i) {
		cr_intv_t *q = &qr->r[i];
		int32_t st1 = (int32_t)q->x, en1 = q->y, cov_st = 0, cov_en = 0;
		int64_t j;
		char *ctg = qr->ctg[q->x>>32].name;
		n_b = cr_overlap(cr, ctg, st1, en1, &b, &m_b);
		for (j = 0; j < n_b; ++j) {
			cr_intv_t *r = &cr->r[b[j]];
			int32_t st0 = cr_st(r), en0 = cr_en(r);
			if (st0 < st1) st0 = st1;
			if (en0 > en1) en0 = en1;
			if (st0 > cov_en) {
				if (cov_en > cov_st) printf("%s\t%d\t%d\n", ctg, cov_st, cov_en);
				cov_st = st0, cov_en = en0;
			} else cov_en = cov_en > en0? cov_en : en0;
		}
		if (cov_en > cov_st) printf("%s\t%d\t%d\n", ctg, cov_st, cov_en);
	}
	free(b);
	cr_destroy(qr);
	cr_destroy(cr);
	return 0;
}

int main_flt(int argc, char *argv[])
{
	cgranges_t *cr;
	ketopt_t o = KETOPT_INIT;
	int64_t m_b = 0, *b = 0, n_b;
	int c, win = 0, vcf_in = 0, paf_in = 0, test_con = 0, non_sat = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0}, out = {0,0,0};

	while ((c = ketopt(&o, argc, argv, 1, "cw:Cvp", 0)) >= 0) {
		if (c == 'c') vcf_in = 1;
		else if (c == 'p') paf_in = 1;
		else if (c == 'C') test_con = 1;
		else if (c == 'v') non_sat = 1;
		else if (c == 'w') win = atol(o.arg);
	}

	if (argc - o.ind < 1 || (argc - o.ind < 2 && isatty(0))) {
		printf("Usage: bedtk flt [options] <loaded.bed> <streamed.bed>\n");
		printf("Options:\n");
		printf("  -c      the second input is VCF\n");
		printf("  -p      the second input is PAF\n");
		printf("  -C      print records contained in the union of <loaded.bed>\n");
		printf("  -v      print non-satisfying records\n");
		printf("  -w INT  window size [0]\n");
		return 1;
	}

	if (vcf_in && paf_in) {
		fprintf(stderr, "ERROR: -c and -p can't be applied at the same time\n");
		return 1;
	}

	cr = read_bed3(argv[o.ind]);
	assert(cr);
	cr_index2(cr, 1);

	fp = o.ind + 1 < argc && strcmp(argv[o.ind + 1], "-")? gzopen(argv[o.ind + 1], "r") : gzdopen(0, "r");
	assert(fp);
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		int32_t st1, en1, st2, en2, sat = 0;
		int64_t i;
		char *ctg, *rest;
		if (paf_in) {
			ctg = parse_paf(str.s, &st1, &en1, &rest);
		} else if (vcf_in) {
			ctg = parse_vcf(str.s, &st1, &en1, &rest);
			if (ctg == 0 && rest && rest[0] == '#')
				puts(rest);
		} else ctg = parse_bed3b(str.s, &st1, &en1, &rest);
		if (ctg == 0) continue;
		st2 = st1 - win, en2 = en1 + win;
		if (st2 < 0) st2 = 0;
		n_b = cr_overlap(cr, ctg, st2, en2, &b, &m_b);
		if (test_con) {
			for (i = 0, sat = 0; i < n_b; ++i) {
				cr_intv_t *r = &cr->r[b[i]];
				if (cr_st(r) <= st2 && en2 <= cr_en(r)) {
					sat = 1;
					break;
				}
			}
		} else sat = (n_b > 0);
		out.l = 0;
		if ((sat && !non_sat) || (!sat && non_sat)) {
			if (paf_in) {
				*(rest - 1) = '\t';
				fwrite(str.s, 1, str.l, stdout);
				putchar('\n');
			} else if (vcf_in) {
				mm_sprintf_lite(&out, "%s\t", ctg);
				fwrite(out.s, 1, out.l, stdout);
				puts(rest);
			} else {
				mm_sprintf_lite(&out, "%s\t%d\t%d", ctg, st1, en1);
				fwrite(out.s, 1, out.l, stdout);
				if (rest) puts(rest);
				else putchar('\n');
			}
		}
	}
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);

	free(b);
	cr_destroy(cr);
	return 0;
}

int main_cov(int argc, char *argv[])
{
	cgranges_t *cr;
	gzFile fp;
	ketopt_t o = KETOPT_INIT;
	kstream_t *ks;
	kstring_t str = {0,0,0}, out = {0,0,0};
	int64_t m_b = 0, *b = 0, n_b;
	int c, cnt_only = 0, contained = 0, print_depth = 0;

	while ((c = ketopt(&o, argc, argv, 1, "cCd", 0)) >= 0) {
		if (c == 'c') cnt_only = 1;
		else if (c == 'C') contained = 1;
		else if (c == 'd') print_depth = 1;
	}

	if (argc - o.ind < 1 || (argc - o.ind < 2 && isatty(0))) {
		printf("Usage: bedtk cov [options] <loaded.bed> <streamed.bed>\n");
		printf("Options:\n");
		printf("  -c       only count; no breadth of depth\n");
		printf("  -C       containment only\n");
		return 1;
	}

	cr = read_bed3(argv[o.ind]);
	assert(cr);
	cr_index(cr);

	fp = o.ind+1 < argc && strcmp(argv[o.ind + 1], "-")? gzopen(argv[o.ind + 1], "r") : gzdopen(0, "r");
	assert(fp);
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		int32_t st1, en1;
		char *ctg;
		int64_t j, cnt = 0, cov = 0, cov_st = 0, cov_en = 0, depth = 0;
		ctg = parse_bed3(str.s, &st1, &en1);
		if (ctg == 0) continue;
		if (contained)
			n_b = cr_contain(cr, ctg, st1, en1, &b, &m_b);
		else
			n_b = cr_overlap(cr, ctg, st1, en1, &b, &m_b);
		out.l = 0;
		if (!cnt_only) {
			for (j = 0; j < n_b; ++j) {
				cr_intv_t *r = &cr->r[b[j]];
				int32_t st0 = cr_st(r), en0 = cr_en(r);
				if (st0 < st1) st0 = st1;
				if (en0 > en1) en0 = en1;
				if (st0 > cov_en) {
					cov += cov_en - cov_st;
					cov_st = st0, cov_en = en0;
				} else cov_en = cov_en > en0? cov_en : en0;
				++cnt, depth += en0 - st0;
			}
			cov += cov_en - cov_st;
			if (print_depth)
				mm_sprintf_lite(&out, "%s\t%d\t%d\t%d\t%d\t%d\n", ctg, st1, en1, (int)cnt, (int)cov, (int)depth);
			else
				mm_sprintf_lite(&out, "%s\t%d\t%d\t%d\t%d\n", ctg, st1, en1, (int)cnt, (int)cov);
		} else {
			mm_sprintf_lite(&out, "%s\t%d\t%d\t%d\n", ctg, st1, en1, (int)n_b);
		}
		fputs(out.s, stdout);
	}
	free(b);
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);

	cr_destroy(cr);
	return 0;
}

int main_sub(int argc, char *argv[])
{
	cgranges_t *cr;
	gzFile fp;
	ketopt_t o = KETOPT_INIT;
	kstream_t *ks;
	kstring_t str = {0,0,0}, out = {0,0,0};
	int64_t m_b = 0, *b = 0, n_b;
	int32_t c, min_len = 0;

	while ((c = ketopt(&o, argc, argv, 1, "l:", 0)) >= 0) {
		if (c == 'l') min_len = atoi(o.arg);
	}

	if (argc - o.ind < 2) {
		printf("Usage: bedtk sub <minuend.bed> <subtrahend.bed>\n");
		printf("Note: <subtrahend.bed> is loaded into memory.\n");
		return 1;
	}

	cr = read_bed3(argv[o.ind + 1]);
	assert(cr);
	cr_index2(cr, 1);

	fp = strcmp(argv[o.ind], "-")? gzopen(argv[o.ind], "r") : gzdopen(0, "r");
	assert(fp);
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
		int32_t st1, en1, x;
		char *ctg, *rest;
		int64_t j;
		ctg = parse_bed3b(str.s, &st1, &en1, &rest);
		if (ctg == 0) continue;
		n_b = cr_overlap(cr, ctg, st1, en1, &b, &m_b);
		for (j = 0, x = st1; j < n_b; ++j) {
			cr_intv_t *r = &cr->r[b[j]];
			int32_t st0 = cr_st(r), en0 = cr_en(r);
			if (st0 < st1) st0 = st1;
			if (en0 > en1) en0 = en1;
			if (st0 > x && st0 - x >= min_len) {
				out.l = 0;
				mm_sprintf_lite(&out, "%s\t%d\t%d", ctg, x, st0);
				if (rest) mm_sprintf_lite(&out, "%s", rest);
				puts(out.s);
			}
			x = en0;
		}
		if (x < en1 && en1 - x >= min_len) {
			out.l = 0;
			mm_sprintf_lite(&out, "%s\t%d\t%d", ctg, x, en1);
			if (rest) mm_sprintf_lite(&out, "%s", rest);
			puts(out.s);
		}
	}
	free(b);
	free(out.s);
	free(str.s);
	ks_destroy(ks);
	gzclose(fp);

	cr_destroy(cr);
	return 0;
}

int main_merge(int argc, char *argv[])
{
	cgranges_t *cr;
	ketopt_t o = KETOPT_INIT;
	int c, assume_srt = 0;

	while ((c = ketopt(&o, argc, argv, 1, "s", 0)) >= 0) {
		if (c == 's') assume_srt = 1;
	}

	if (argc - o.ind < 1 && isatty(0)) {
		printf("Usage: bedtk merge [options] <in.bed>\n");
		printf("Options:\n");
		printf("  -s       assume the input is sorted (NOT implemented yet)\n");
		return 1;
	}

	if (assume_srt) {
		fprintf(stderr, "ERROR: NOT implemented yet\n");
	} else {
		int64_t i;
		cr = read_bed3(o.ind < argc? argv[o.ind] : 0);
		assert(cr);
		if (!cr_is_sorted(cr)) cr_sort(cr);
		cr_merge_pre_index(cr);
		for (i = 0; i < cr->n_r; ++i) {
			const cr_intv_t *r = &cr->r[i];
			printf("%s\t%d\t%d\n", cr->ctg[r->x>>32].name, (int32_t)r->x, r->y);
		}
		cr_destroy(cr);
	}

	return 0;
}

int main_sum(int argc, char *argv[])
{
	cgranges_t *cr;
	ketopt_t o = KETOPT_INIT;
	int c, merge = 0;
	int64_t sum = 0;

	while ((c = ketopt(&o, argc, argv, 1, "m", 0)) >= 0) {
		if (c == 'm') merge = 1;
	}

	if (argc - o.ind < 1 && isatty(0)) {
		printf("Usage: bedtk sum [options] <in.bed>\n");
		printf("Options:\n");
		printf("  -m       merge overlapping regions\n");
		return 1;
	}

	if (!merge) {
		gzFile fp;
		kstream_t *ks;
		kstring_t str = {0,0,0};
		fp = o.ind < argc && strcmp(argv[o.ind], "-")? gzopen(argv[o.ind], "r") : gzdopen(0, "r");
		assert(fp);
		ks = ks_init(fp);
		while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
			int32_t st1, en1;
			char *ctg;
			ctg = parse_bed3(str.s, &st1, &en1);
			if (ctg == 0) continue;
			sum += en1 - st1;
		}
		free(str.s);
		ks_destroy(ks);
		gzclose(fp);
	} else {
		int64_t i;
		cr = read_bed3(o.ind < argc? argv[o.ind] : 0);
		assert(cr);
		if (!cr_is_sorted(cr)) cr_sort(cr);
		cr_merge_pre_index(cr);
		for (i = 0; i < cr->n_r; ++i)
			sum += cr->r[i].y - (int32_t)cr->r[i].x;
		cr_destroy(cr);
	}
	printf("%lld\n", (long long)sum);
	return 0;
}

int main_sort(int argc, char *argv[])
{
	cgranges_t *cr;
	ketopt_t o = KETOPT_INIT;
	int c;
	int64_t i;
	bed_rest_t rest;
	char *fn_order = 0;
	kstring_t str = {0,0,0};

	while ((c = ketopt(&o, argc, argv, 1, "s:", 0)) >= 0)
		if (c == 's') fn_order = o.arg;

	if (argc - o.ind < 1 && isatty(0)) {
		printf("Usage: bedtk sort [options] <in.bed>\n");
		printf("Options:\n");
		printf("  -s FILE   list of contig IDs to specify the order []\n");
		return 1;
	}

	cr = read_bed3b(o.ind < argc? argv[o.ind] : 0, &rest, fn_order);
	assert(cr);
	if (!cr_is_sorted(cr)) cr_sort(cr);
	for (i = 0; i < cr->n_r; ++i) {
		const cr_intv_t *r = &cr->r[i];
		str.l = 0;
		mm_sprintf_lite(&str, "%s\t%d\t%d", cr->ctg[r->x>>32].name, (int32_t)r->x, r->y);
		if (rest.a[r->label].l) mm_sprintf_lite(&str, "%s", rest.a[r->label].s);
		puts(str.s);
	}
	free(str.s);
	cr_destroy(cr);

	for (i = 0; i < rest.n; ++i) free(rest.a[i].s);
	free(rest.a);
	return 0;
}

/*****************
 * Main function *
 *****************/

static int usage(FILE *fp)
{
	fprintf(fp, "Usage: bedtk <command> <arguments>\n");
	fprintf(fp, "Command:\n");
	fprintf(fp, "  isec      intersection (bedtools intersect)\n");
	fprintf(fp, "  flt       filter BED/VCF file (bedtools intersect/window)\n");
	fprintf(fp, "  cov       breadth of coverage (bedtools coverage)\n");
	fprintf(fp, "  sub       subtraction (bedtools subtract)\n");
	fprintf(fp, "  merge     merge overlapping regions (bedtools merge)\n");
	fprintf(fp, "  sort      sort regions (bedtools sort)\n");
	fprintf(fp, "  sum       total region length\n");
	fprintf(fp, "  version   version number\n");
	return fp == stdout? 0 : 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage(stdout);
	if (strcmp(argv[1], "isec") == 0) return main_isec(argc-1, argv+1);
	else if (strcmp(argv[1], "flt") == 0) return main_flt(argc-1, argv+1);
	else if (strcmp(argv[1], "cov") == 0) return main_cov(argc-1, argv+1);
	else if (strcmp(argv[1], "sub") == 0) return main_sub(argc-1, argv+1);
	else if (strcmp(argv[1], "merge") == 0) return main_merge(argc-1, argv+1);
	else if (strcmp(argv[1], "sum") == 0) return main_sum(argc-1, argv+1);
	else if (strcmp(argv[1], "sort") == 0) return main_sort(argc-1, argv+1);
	else if (strcmp(argv[1], "version") == 0) {
		puts(BEDTK_VERSION);
		return 0;
	} else {
		fprintf(stderr, "ERROR: unrecognized command '%s'. Abort!\n", argv[1]);
		return 1;
	}
}
