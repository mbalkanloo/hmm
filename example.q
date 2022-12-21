\c 500 500
\l hmm.q

/ from wikipedia entry for forward-backward algorithm
S:`R`N
V:`U`N
a:(.7 .3;.3 .7)
b:(.9 .1;.2 .8)
pi:.5 .5
O:V?seq:`U`U`N`U`U

show "forward"
show flip alpha:.hmm.fwd[pi;a;b;O]
show "backward"
show beta:.hmm.bwd[a;b;O]
show "forward-backward gamma"
show .hmm.gamma[alpha 0;.hmm.bwd[a;b;O]]
show "viterbi"
show `observed`state!flip seq,'.hmm.viterbi[pi;a;b;O;S]
show "baum-welch"
show r:.hmm.baumWelch[pi;a;b;S;V;O;.03;100]
show `observed`state!flip seq,'.hmm.viterbi[r`pi;r`a;r`b;O;S]

/ example from mit slides from Emilio Frazzoli 

S:`LA`NY
V:`LA`NY`
a:(.5 .5;.5 .5)
b:(.4 .1 .5;.1 .5 .4)
pi:1. .0
O:V?seq:``LA`LA``NY``NY`NY`NY``NY`NY`NY`NY`NY```LA`LA`NY

show "forward"
show flip alpha:.hmm.fwd[pi;a;b;O]
show "backward"
show beta:.hmm.bwd[a;b;O]
show "forward-backward gamma"
show .hmm.gamma[alpha 0;.hmm.bwd[a;b;O]]
show "viterbi"
show `observed`state!flip seq,'.hmm.viterbi[pi;a;b;O;S]
show "baum-welch"
show r:.hmm.baumWelch[pi;a;b;S;V;O;.03;100]
show `observed`state!flip seq,'.hmm.viterbi[r`pi;r`a;r`b;O;S]

/ example from wikipedia entry for Baum-Welch 

S:0 1
V:`E`N
a:(.5 .5;.3 .7) 
b:(.3 .7;.8 .2)
pi:.1 .9
O:V?seq:`N`N`N`N`N`E`E`N`N`N 

show "forward"
show flip alpha:.hmm.fwd[pi;a;b;O]
show "backward"
show beta:.hmm.bwd[a;b;O]
show "forward-backward gamma"
show .hmm.gamma[alpha 0;.hmm.bwd[a;b;O]]
show "viterbi"
show `observed`state!flip seq,'.hmm.viterbi[pi;a;b;O;S]
show "baum-welch"
show r:.hmm.baumWelch[pi;a;b;S;V;O;.03;100]
show `observed`state!flip seq,'.hmm.viterbi[r`pi;r`a;r`b;O;S]


