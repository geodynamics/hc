# generate a synthetic density anomaly model in the short SH format
# (see README.TXT): three layers, random coefficients with power
# decreasing with degree. used by propagator_stability.sh
#
# usage: gawk -v L=63 -f gen_synth_dens.awk > dens.ab
#
BEGIN{
  srand(1);
  nlayer = 3;
  split("2000 1000 300",z," ");	# depths in km
  print nlayer;
  for(i=1;i <= nlayer;i++){
    print z[i],L;
    for(l=0;l <= L;l++){
      for(m=0;m <= l;m++){
	s = (l > 0) ? 1.0/(l*sqrt(2*l+1)) : 0;
	a = (2*rand()-1)*s;
	b = (m > 0) ? (2*rand()-1)*s : 0.0;
	printf("%.8e %.8e\n",a,b);
      }
    }
  }
}
