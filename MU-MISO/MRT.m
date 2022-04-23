function P = MRT(H,pow)
pre = H';
P = sqrt(pow/trace(pre*pre'))*pre;