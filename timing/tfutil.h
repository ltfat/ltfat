static inline int positiverem(int a,int b)
{
  int c;
  c=a%b;
  return(c<0 ? c+b : c);
}
