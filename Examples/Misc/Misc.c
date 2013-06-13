

#include <stdio.h>
#include <Lgm_Misc.h>

int main(void) {
  char *str1 = "I am a String";
  char str2[50];

  printf("Orig: %s\n", str1);
  
  Lgm_ReplaceSubString( str2, str1, "am", "was" );
  printf("New: %s\n", str2);




  return(0);
}





