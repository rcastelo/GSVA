// 2012, Georg Sauthoff <mail@georg.so>
// GPLv3+

#include <stdio.h>
#include <stdlib.h>


#include <sys/types.h>
#include <unistd.h>


static void fill_some(size_t n, char ch)
{
  size_t bytes = n * (size_t) 1024 * (size_t) 1024;
  char *b = malloc(bytes);
  if (!b) {
    fprintf(stderr, "Failed to allocate %zu MiB\n", n);
    exit(1);
  }
  printf("Allocating %zu MiBs\n", n);
  for (size_t i = 0; i<bytes; i+=1024) {
    b[i] = ch;
  }
  sleep(1);
}

int main(int argc, char **argv)
{
  char ch = argv[1][0];

  printf("Parent PID is %zd\n", (ssize_t)getpid());
  fill_some(atoi(argv[2]), ch);

  for (int r = 3; r<argc; ++r) {
    pid_t ret = fork();
    if (ret < 0) {
      perror("Fork failed");
      exit(2);
    }
    if (!ret) {
      printf("New Child: %zd\n", (ssize_t)getpid());
      fill_some(atoi(argv[r]), ch);
      exit(0);
    }
  }
  sleep(argc);
  return 0;
}
