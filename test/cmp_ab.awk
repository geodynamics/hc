# compare two ASCII files with numeric content (e.g. spherical
# harmonic expansions written by hc) token by token
#
# usage: gawk -f cmp_ab.awk -v tol=1e-6 file1 file2
#
# the error measure is max |a-b| / max |a|, taken over all numeric
# tokens of file1, i.e. relative to the largest value in the file
#
# exit status 0 if error <= tol, 1 otherwise
#
BEGIN{
  if(tol == "")
    tol = 1e-6;
  n = 0;
  amax = 0;
}
FNR == NR{
  for(i=1;i <= NF;i++){
    n++;
    a[n] = $i;
    if($i+0 == $i){		# numeric
      v = ($i < 0) ? -$i : $i;
      if(v > amax)
	amax = v;
    }
  }
  next;
}
{
  for(i=1;i <= NF;i++){
    m++;
    b[m] = $i;
  }
}
END{
  if(n != m){
    printf("cmp_ab: token count differs: %i vs %i\n",n,m);
    exit(1);
  }
  dmax = 0;
  for(i=1;i <= n;i++){
    if((a[i]+0 == a[i]) && (b[i]+0 == b[i])){
      d = a[i] - b[i];
      if(d < 0)
	d = -d;
      if(d > dmax)
	dmax = d;
    }else if(a[i] != b[i]){
      printf("cmp_ab: non-numeric token differs: %s vs %s\n",a[i],b[i]);
      exit(1);
    }
  }
  if(amax == 0)
    amax = 1;
  err = dmax/amax;
  printf("cmp_ab: max rel. diff %.3e (tol %.1e) ",err,tol);
  if(err > tol){
    printf("FAILED\n");
    exit(1);
  }else{
    printf("ok\n");
    exit(0);
  }
}
