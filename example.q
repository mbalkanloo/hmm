\l hmm.q

/ example from wikipedia entry for Baum-Welch 

/ state space
S:0 1
/ observation space
V:`E`N
/ observation probability matrix
b:flip (.3 .7;.8 .2)
/ state transition probability matrix
a:(.5 .5;.3 .7)
/ initial state probability
pi:.2 .8
/ observations (indices of observation space)
O:V?sequence:`N`N`N`N`N`E`E`N`N`N

show .hmm.forward[pi;a;b;O]
show .hmm.backward[a;b;O]
show .hmm.viterbi[pi;a;b;O;S]
show r:.hmm.baumWelch[pi;a;b;S;V;O;.001;100]
show .hmm.viterbi[r`pi;r`a;r`b;O;S]
