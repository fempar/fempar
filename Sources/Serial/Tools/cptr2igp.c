//#include <stdio.h>
#include <stdlib.h>
void cptr2igp(void * cptr, long int * igp)
{
  /* printf("%p\n", cptr);  */
  /* printf("%d\n", *((int *)cptr));  */
  *igp = (long int) cptr;
  /* printf("%lu\n", *igp); */
}

/* int main ( int argc, char ** argv ) */
/* { */
/*   int i = 3; */
/*   long int  igp; */
   
/*   printf("%p\n",&i); */
/*   cptr2igp((void * ) &i, &igp); */
/*   exit(0); */
/* } */
