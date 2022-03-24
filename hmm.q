/ based on Rabiner (1989)

/ a      transition matrix
/ b      observation probability matrix
/ pi     initial state probability
/ S      state space
/ V      observation space
/ O      observations
/ alpha  forward matrix
/ beta   backward matrix
/ g      global state matrix symbol for dynamic algorithms

\d .hmm

forward:{[pi;a;b;O]
	.hmm.g:();
	b:flip b;
	.hmm.g,:(pi;pi*b first O);
	forwardInduction[`.hmm.g;b;a;]each 1_O;
	r:1_.hmm.g;
	delete g from `.hmm;
	r}

forwardInduction:{[g;b;a;o]g upsert b[o]*last[value g]mmu a}

backward:{[a;b;O]
	.hmm.g:();
	b:flip b;
	.hmm.g,:enlist count[first b]#1f;
	backwardInduction[`.hmm.g;a;b;] each reverse O;
	r:1_reverse .hmm.g;
	delete g from `.hmm;
	r}

backwardInduction:{[g;a;b;o]g upsert sum each a*\:b[o]*last[value g]}

viterbi:{[pi;a;b;O;S]S first each idesc each forward[pi;a;b;O]}

gamma:{[alpha;beta]a:alpha*beta;a%sum each a}

xi:{[alpha;beta;a;b;O]
	r:flip each(1_beta)*'flip each(-1_alpha)*'flip each flip[a]*/:flip[b]@/:-1_next O;
	r%sum each raze each r}

PI:{[alpha;beta]gamma[alpha;beta]0}

A:{[alpha;beta;a;b;O]sum[xi[alpha;beta;a;b;O]]%sum -1_gamma[alpha;beta]}

B:{[alpha;beta;S;V;O]
	m:gamma[alpha;beta];
	l:{sum x[y 0]where x[0]=y 1};
	r:(0N;j)#l[flip[O,'m];]each(1+til count S)cross til j:count V;
	r%sum m}

reestimate:{[pi;a;b;S;V;O]
	alpha:forward[pi;a;b;O];
	beta:backward[a;b;O];
	`alpha`beta`pi`a`b!(alpha;beta;PI[alpha;beta];A[alpha;beta;a;b;O];B[alpha;beta;S;V;O])}

baumWelch:{[pi;a;b;S;V;O;t;i]
	/ iteration criteria
	/ difference in observation probability greater than [t]hreshold
	/ iterations less than max [i]terations

	r:`alpha`beta`pi`a`b!(forward[pi;a;b;O];backward[a;b;O];pi;a;b);
	m:-0Wf;
	n:sum last r`alpha;
	while[(t<n-m)&-1<i-:1;
		m:n;
		r:reestimate[r`pi;r`a;r`b;S;V;O];
		n:sum last forward[r`pi;r`a;r`b;O]];
	r}

\d .

/ / example from wikipedia entry for Baum-Welch 
/ S:0 1;
/ V:`E`N;
/ b:(.3 .7;.8 .2);
/ a:(.5 .5;.3 .7);
/ pi:.2 .8;
/ O:V?seq:`N`N`N`N`N`E`E`N`N`N; 
/ show .hmm.forward[pi;a;b;O]
/ show .hmm.backward[a;b;O]
/ show .hmm.viterbi[pi;a;b;O;S]
/ show each .hmm.baumWelch[pi;a;b;S;V;O;.001;100]
