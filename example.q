\l hmm.q

/ example from wikipedia entry for Baum-Welch 
S:0 1;
V:`E`N;
b:(.3 .7;.8 .2);
a:(.5 .5;.3 .7);
pi:.2 .8;
O:V?seq:`N`N`N`N`N`E`E`N`N`N; 
show .hmm.forward[pi;a;b;O]
show .hmm.backward[a;b;O]
show .hmm.viterbi[pi;a;b;O;S]
show r:.hmm.baumWelch[pi;a;b;S;V;O;.001;100]
show .hmm.viterbi[r`pi;r`a;r`b;O;S]
