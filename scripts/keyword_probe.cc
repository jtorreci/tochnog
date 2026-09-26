#include <cstdio>
#include <cstring>
void db_initialize(long int[],long int[]);
long int db_number(const char[]);
static long int a[200000], b[200000];
int main(){ db_initialize(a,b); char buf[512];
  while(fgets(buf,sizeof buf,stdin)){ buf[strcspn(buf,"\n")]=0; if(!*buf) continue;
    printf("%s %s\n", db_number(buf)>=0?"OK":"NO", buf);} }
