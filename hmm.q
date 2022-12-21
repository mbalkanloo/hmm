/ based on Rabiner (1989)
/ for discrete observation distribution

\d .hmm

/ a      transition matrix
/ b      observation probability matrix
/ pi     initial state probability
/ S      state space
/ V      observation space
/ O      observations
/ alpha  forward matrix
/ beta   backward matrix

/ scale factor utility functions
sf:{reciprocal sum x}
sfr:{sf raze x}

/ forward 
fwd:{[pi;a;b;O]
	i*:c:sf i:pi*b[;first O];
	`.hmm.c set enlist c;
	l:{[a;b;x;y].hmm.c,:c:sf r:b[;y]*x mmu a;r*c};
	r:enlist[pi*c:sf pi],enlist[i],l[a;b;;]\[i;1_O];
	c:enlist[c],.hmm.c;
	delete c from `.hmm;
	(r;c)}

/ backward
bwd:{[a;b;O]
	i*:sf i:count[first flip b]#1f;
	l:{[a;b;x;y]r*sf r:sum each a*\:b[;y]*x};
	reverse enlist[i],l[a;b;;]\[i;reverse O]}

viterbi:{[pi;a;b;O;S]S first each idesc each 1_fwd[pi;a;b;O] 0}

gamma:{[alpha;beta]1_a*sf each a:alpha*beta}

xi:{[alpha;beta;a;b;O]m:(-1_1_alpha)*\:a;n:2#/:enlist each(flip[b]@/:1_O)*2_beta;r*sfr each r:m*n}

PI:{[alpha;beta]gamma[alpha;beta]0}
A:{[alpha;beta;a;b;O]An[alpha;beta;a;b;O]%Ad[alpha;beta]}
An:{[alpha;beta;a;b;O]sum xi[alpha;beta;a;b;O]}
Ad:{[alpha;beta]sum -1_gamma[alpha;beta]}
B:{[alpha;beta;S;V;O]g:gamma[alpha;beta];Bn[g;S;V;O]%Bd g}
Bn:{[g;S;V;O](0N;j)#Bnl[g;O;]each til[count S]cross til j:count V}
Bnl:{[g;O;x]sum g[;x 0] where O=x 1}
Bd:{[g]sum g}

reestimate:{[pi;a;b;S;V;O]
	alpha:fwd[pi;a;b;O] 0;
	beta:bwd[a;b;O];
	`pi`a`b!(PI[alpha;beta];A[alpha;beta;a;b;O];B[alpha;beta;S;V;O])}

baumWelch:{[pi;a;b;S;V;O;t;i]
	/ iterate until 
	/ consecutive observation probability delta less than [t]hreshold
	/ iterations equal to max [i]terations
	r:`pi`a`b!(pi;a;b);
	m:-0Wf;
	n:neg sum log fwd[pi;a;b;O] 1;
	I:i-1;
	while[(t<n-m)&(not n=m)&-1<i-:1;
		m:n;
		r:reestimate[r`pi;r`a;r`b;S;V;O];
		n:neg sum log fwd[r`pi;r`a;r`b;O] 1];
	r[`iterations]:I-i;
	r}
