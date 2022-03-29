/ based on Rabiner (1989)
/ implemented with scaling

\d .hmm

/ a      transition matrix
/ b      observation probability matrix (flipped)
/ pi     initial state probability
/ S      state space
/ V      observation space
/ O      observations
/ alpha  forward matrix
/ beta   backward matrix
/ g      global state matrix symbol for dynamic algorithms

forward:{[pi;a;b;O]
	/ returns a tuple of alpha and scale factors
	.hmm.g:();
	.hmm.f:();
	i:(pi;pi*b first O);
	.hmm.g,:i%sum each i;
	.hmm.f,:reciprocal each sum each i;
	forwardInduction[`.hmm.g;`.hmm.f;b;a;]each 1_O;
	r:1_.hmm.g;
	f:1_.hmm.f;
	delete g, f from `.hmm;
	(r;f)}

forwardInduction:{[g;f;b;a;o]r:b[o]*last[value g]mmu a;m:reciprocal sum r;f upsert m;g upsert r*m}

backward:{[f;a;b;O]
	/ expects scale [f]actors
	.hmm.g:();
	i:enlist count[first b]#1f;
	.hmm.g,:i%f 0; 
	backwardInduction[`.hmm.g;a;b;] each reverse O,'f;
	r:1_reverse .hmm.g;
	delete g from `.hmm;
	r}

backwardInduction:{[g;a;b;o]r:sum each a*\:b[o 0]*last[value g];g upsert r*o 1}

/ NOTE this is not the standard Viterbi definition
viterbi:{[pi;a;b;O;S]S first each idesc each forward[pi;a;b;O] 0}

gamma:{[alpha;beta]a:alpha*beta;a%sum each a}

/ flip each utility for flipping submatrices
fe:flip each

xi:{[alpha;beta;a;b;O]
	r:fe(1_beta)*'fe(-1_alpha)*'fe flip[a]*/:b@/:-1_next O;
	r%sum each raze each r}

PI:{[alpha;beta]gamma[alpha;beta]0}

A:{[alpha;beta;a;b;O]sum[xi[alpha;beta;a;b;O]]%sum -1_gamma[alpha;beta]}

B:{[alpha;beta;S;V;O]
	m:gamma[alpha;beta];
	l:{sum x[y 0]where x[0]=y 1};
	r:(0N;j)#l[flip[O,'m];]each(1+til count S)cross til j:count V;
	flip r%sum m}

reestimate:{[pi;a;b;S;V;O]
	alpha:forward[pi;a;b;O];
	beta:backward[alpha 1;a;b;O];
	`alpha`beta`pi`a`b!(
		alpha;
		beta;
		PI[alpha 0;beta];
		A[alpha 0;beta;a;b;O];
		B[alpha 0;beta;S;V;O])}

baumWelch:{[pi;a;b;S;V;O;t;i]
	/ iterate until 
	/ difference in observation probability greater than [t]hreshold
	/ iterations less than max [i]terations
	r:`alpha`beta`pi`a`b!(forward[pi;a;b;O];backward[a;b;O];pi;a;b);
	m:-0Wf;
	n:neg sum log r[`alpha] 1;
	while[(t<n-m)&-1<i-:1;
		m:n;
		r:reestimate[r`pi;r`a;r`b;S;V;O];
		n:neg sum log forward[r`pi;r`a;r`b;O] 1];
	r}

\d .
