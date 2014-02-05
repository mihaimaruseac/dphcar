/**
 * Differentially-private high-confidence association rule extractor.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gmp.h>
#include <mpfr.h>

struct {
	/** total number of items */
	int n;
	/** number of items passing min threshold */
	int N;
	/** total number of rules */
	double R;
	/** number of transactions */
	int t;
	/** \alpha, \beta parameters */
	int alpha, beta;
	/** a and b parameters */
	int a, b;
	/** how many exponential-mechanism calls */
	int k;
	/** how many rules */
	int K;
	/** min-conf */
	float c;
	/** \epsilon */
	float epsilon;
	/** summarized itemsets filename */
	char *fname;
} args;

static void usage(char *prg)
{
	fprintf(stderr, "Usage: %s OPTIONS\n", prg);
	fprintf(stderr, "\n");
	fprintf(stderr, "\tOPTIONS:\n");
	fprintf(stderr, "\t\t-n n      \t(REQ, int)\tnumber of items in the original dataset\n");
	fprintf(stderr, "\t\t-t t      \t(opt, int)\tnumber of transactions in the original dataset\n");
	fprintf(stderr, "\t\t-a a      \t(opt, int)\n");
	fprintf(stderr, "\t\t-b b      \t(opt, int)\tinitial selection center\n");
	fprintf(stderr, "\t\t-c c      \t(opt, float)\tmin confidence value\n");
	fprintf(stderr, "\t\t-A alpha  \t(REQ, int)\tx-cut-off\n");
	fprintf(stderr, "\t\t-B beta   \t(opt, int)\ty-cut-off\n");
	fprintf(stderr, "\t\t-f file   \t(REQ, filename)\tthe file with summarized itemsets above a threshold\n");
	fprintf(stderr, "\t\t-e epsilon\t(REQ, float)\tthe privacy budget\n");
	fprintf(stderr, "\t\t-k k      \t(opt, int)\tthe number of exponential-mechanism invocations\n");
	fprintf(stderr, "\t\t-K K      \t(opt, int)\tnumber of required rules\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "\tCONSTRAINTS and DEFAULTS:\n");
	fprintf(stderr, "\t\t-no given value should be 0 or negative (ignored otherwise)\n");
	fprintf(stderr, "\t\t-if t is not given it is read from the first line of file\n");
	fprintf(stderr, "\t\t-either (a,b) or c must be given\n");
	fprintf(stderr, "\t\t\t-if c is given (a,b) is chosen in the middle of the min-conf line\n");
	fprintf(stderr, "\t\t\t-if (a,b) is given then c=b/a\n");
	fprintf(stderr, "\t\t\t-(a,b) overwrites any c definition\n");
	fprintf(stderr, "\t\t-alpha and beta should be at least the threshold used to generate file (invalid results otherwise)\n");
	fprintf(stderr, "\t\t-if beta is not given it is equal to alpha\n");
	fprintf(stderr, "\t\t-one of k and K must be given\n");
	fprintf(stderr, "\t\t\t-k is computed as half of K if not given\n");
	fprintf(stderr, "\t\t\t-k overwrites any K definition\n");
	exit(EXIT_FAILURE);
}

static void parse_arguments(int argc, char **argv)
{
	int opt;
	char extra;

	while ((opt = getopt(argc, argv, "n:t:a:b:c:A:B:f:e:k:K:")) != -1)
		switch (opt) {
		case 'n':
			if (sscanf(optarg, "%d%c", &args.n, &extra) != 1)
				goto bad_call;
			break;
		case 't':
			if (sscanf(optarg, "%d%c", &args.t, &extra) != 1)
				goto bad_call;
			break;
		case 'a':
			if (sscanf(optarg, "%d%c", &args.a, &extra) != 1)
				goto bad_call;
			break;
		case 'b':
			if (sscanf(optarg, "%d%c", &args.b, &extra) != 1)
				goto bad_call;
			break;
		case 'c':
			if (sscanf(optarg, "%f%c", &args.c, &extra) != 1)
				goto bad_call;
			break;
		case 'A':
			if (sscanf(optarg, "%d%c", &args.alpha, &extra) != 1)
				goto bad_call;
			break;
		case 'B':
			if (sscanf(optarg, "%d%c", &args.beta, &extra) != 1)
				goto bad_call;
			break;
		case 'f':
			args.fname = strdup(optarg);
			break;
		case 'e':
			if (sscanf(optarg, "%f%c", &args.epsilon, &extra) != 1)
				goto bad_call;
			break;
		case 'k':
			if (sscanf(optarg, "%d%c", &args.k, &extra) != 1)
				goto bad_call;
			break;
		case 'K':
			if (sscanf(optarg, "%d%c", &args.K, &extra) != 1)
				goto bad_call;
			break;
		default: goto bad_call;
		}

	/* consistency checks */
	if (optind != argc)
		goto bad_call; /* extra arguments */
	if (args.n <= 0) {
		fprintf(stderr, "Missing or bad number of items in the original dataset\n");
		goto bad_call;
	}
	if ((args.c <= 0) && ((args.a <= 0) || (args.b <= 0))) {
		fprintf(stderr, "Missing or bad (a, b) or c for selection center\n");
		goto bad_call;
	}
	if (((args.a <= 0) || (args.b <= 0)) && (args.a + args.b != 0)) {
		fprintf(stderr, "Only one of a,b given. Ignored in favor of c\n");
		args.a = args.b = 0;
	}
	if (args.alpha <= 0) {
		fprintf(stderr, "Missing or bad number for alpha\n");
		goto bad_call;
	}
	if (args.beta <= 0)
		args.beta = args.alpha;
	if (!args.fname) {
		fprintf(stderr, "Missing itemsets filename\n");
		goto bad_call;
	}
	if (args.epsilon <= 0) {
		fprintf(stderr, "Missing or bad privacy budget\n");
		goto bad_call;
	}
	if ((args.k <= 0) && (args.K <= 0)) {
		fprintf(stderr, "Missing or bad k or K values\n");
		goto bad_call;
	}
	if (args.k <= 0)
		args.k = args.K / 2;

	/** all good */
	return;

bad_call:
	usage(argv[0]);
}

#define PREC 2048

int main(int argc, char **argv)
{
	mpfr_rnd_t rnd_mode = MPFR_RNDN;
	mpfr_t d, d1, tmp, tmp1;
	mpz_t r, tmpz;
	int i;

	parse_arguments(argc, argv);

	mpfr_init2(d, PREC);
	mpfr_init2(d1, PREC);
	mpfr_init2(tmp, PREC);
	mpfr_init2(tmp1, PREC);
	mpz_init_set_ui(r, 1);
	mpz_init_set_ui(tmpz, 3);

	/* compute total number of rules */
	mpz_pow_ui(tmpz, tmpz, args.n);
	mpz_add(r, r, tmpz);
	mpz_set_ui(tmpz, 2);
	mpz_pow_ui(tmpz, tmpz, args.n + 1);
	mpz_sub(r, r, tmpz);

	/* compute distance over sensibility (power = 1) */
	mpfr_set_ui(d, 2, rnd_mode);
	mpfr_rec_sqrt(d, d, rnd_mode);
	mpfr_set_ui(tmp, ((long)args.a) * ((long)args.a) + ((long)args.b * (long)args.b), rnd_mode);
	mpfr_sqrt(tmp, tmp, rnd_mode);
	mpfr_mul(d, d, tmp, rnd_mode);

	for (i = 1; i < 10; i++) {
		/* tmp = -epsilon/2 */
		mpfr_set_d(tmp, -args.epsilon / 2, rnd_mode);
		mpfr_set_d(tmp1, 1 + (i + 0.0) / 10, rnd_mode);
		mpfr_pow(d1, d, tmp1, rnd_mode);

		/* compute final distance */
		mpfr_mul(tmp, d1, tmp, rnd_mode);
		mpfr_exp(tmp, tmp, rnd_mode);

		/* multiply with number of rules */
		mpfr_mul_z(tmp, tmp, r, rnd_mode);

		printf("Power %d: ", i); mpfr_out_str(stdout, 10, 0, tmp, rnd_mode); printf("\n");
		printf("Power %d: ", i); mpfr_out_str(stdout, 10, 0, d1, rnd_mode); printf("\n");
	}

	mpfr_clear(d);
	mpfr_clear(d1);
	mpfr_clear(tmp);
	mpfr_clear(tmp1);
	mpz_clear(r);
	mpz_clear(tmpz);

	return 0;
}
