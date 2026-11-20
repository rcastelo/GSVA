// cgmemtime - high-water memory usage of a command including descendant processes
//
// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: © 2022 Georg Sauthoff <mail@gms.tf>
//
// MODIFIED VERSION. This file was modified in November 2026 by the authors of
// the GSVA Bioconductor package, to calibrate its memory model, from the
// version of https://github.com/gsauthof/cgmemtime at commit 5bb2029. It is
// distributed under the same license, the GNU General Public License version 3
// or any later version, see the file COPYING. The changes are:
//   - the option -M LIMIT sets a hard memory limit on the cgroup (memory.max),
//     with swap disabled (memory.swap.max=0), as workload managers such as
//     SLURM do with their jobs, and adds to the output the number of
//     processes killed by the OOM killer;
//   - the options -s FILE and -i MSEC write samples of the memory of the
//     cgroup (memory.current) and of its anonymous memory (anon in
//     memory.stat), which excludes the page cache, to FILE every MSEC
//     milliseconds;
//   - the processes left running in the cgroup when the command ends, such
//     as the forked workers of a killed R process, are killed;
//   - the check of the memory controller accepts a list of controllers that
//     ends with a newline;
//   - the result of getopt() is held in an int instead of a char, which can't
//     hold its -1 where char is unsigned, e.g., on ARM;
//   - a short usage hint is shown when no command is given.

// for asprintf(), clone3(), ...
#define _GNU_SOURCE

#include <errno.h>
#include <fcntl.h>        // open(), ...
#include <sys/resource.h> // struct rusage
#include <linux/sched.h>  // clone3(), ...
#include <sched.h>        // clone3(), ...
#include <signal.h>       // kill(), ...
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>        // fprintf(), ...
#include <stdlib.h>       // exit(), ...
#include <string.h>
#include <sys/stat.h>     // mkdirat(), ...
#include <sys/syscall.h>  // SYS_*, ...
#include <sys/wait.h>     // waitid macros, ...
#if __GLIBC__ == 2 && __GLIBC_MINOR__ < 36
    #include <linux/wait.h>   // waitid macros, ...
#endif
#include <time.h>         // clock_gettime(), ...
#include <unistd.h>       // getopt(), ...


static int clone3(struct clone_args *args)
{
    long r = syscall(SYS_clone3, args, sizeof *args);
    return r;
}

static int raw_waitid(idtype_t idtype, id_t id, siginfo_t *infop, int options,
        struct rusage *rusage)
{
    long r = syscall(SYS_waitid, idtype, id, infop, options, rusage);
    return r;
}


struct Args {
    const char *cg_fs_dir;
    char *cg_dir;
    char **child_argv;
    bool machine_readable;
    char delim;
    bool no_systemd_run;
    const char *mem_limit;
    size_t mem_limit_bytes;
    const char *samples;
    long interval_ms;
};
typedef struct Args Args;
#define SYNOPSIS "[OPTION..] COMMAND [COMMAND_OPTION..]"

// short usage, pointing to the help of -h, when the arguments have no command
// or are invalid
static void usage_hint(FILE *o, const char *argv0)
{
    fprintf(o, "Usage: %s %s\n"
            "Try '%s -h' for more information.\n", argv0, SYNOPSIS, argv0);
}

static void help(FILE *o, const char *argv0)
{
    fprintf(o, "%s - high-water memory usage of a command including descendants\n"
            "Usage: %s " SYNOPSIS "\n"
            "\n"
            "Options:\n"
            "    -h            Display this help text\n"
            "    -m CGFS       Cgroup v2 base (default: /sys/fs/cgroup)\n"
            "    -c BASE_DIR   Cgroup directory, must end in XXXXXX\n"
            "                  (default: $CGFS/user.slice/user-$UID.slice/user@$UID.service/cgmt-$RANDOM)\n"
            "    -t            machine readable output (delimited columns)\n"
            "    -d CHAR       column delimiter (default: ';')\n"
            "    -Z            disable falling back to systemd-run\n"
            "    -M LIMIT      hard memory limit of the cgroup (memory.max), in bytes\n"
            "                  or with suffix K, M, G or T (powers of 1024), with swap\n"
            "                  disabled (memory.swap.max=0), which adds to the output\n"
            "                  the number of processes killed by the OOM killer\n"
            "    -s FILE       sample the memory of the cgroup (memory.current) every\n"
            "                  -i milliseconds, writing each sample to FILE as a line\n"
            "                  with the time in seconds since the epoch, the memory\n"
            "                  and its anonymous memory (anon in memory.stat), which\n"
            "                  excludes the page cache, in KiB\n"
            "    -i MSEC       sampling interval of -s in milliseconds (default: 100)\n"
            "\n"
            "2022, Georg Sauthoff, GPLv3+\n"
            "\n"
            , argv0, argv0);
}

// parses a memory size given in bytes or with suffix K, M, G or T
static int parse_mem(const char *s, size_t *bytes)
{
    char *end = 0;
    double x = strtod(s, &end);
    if (end == s || x <= 0)
        return -1;
    double fac = 1;
    switch (*end) {
        case 0:                         break;
        case 'K': case 'k': fac = 1024.0; break;
        case 'M': case 'm': fac = 1024.0 * 1024; break;
        case 'G': case 'g': fac = 1024.0 * 1024 * 1024; break;
        case 'T': case 't': fac = 1024.0 * 1024 * 1024 * 1024; break;
        default: return -1;
    }
    if (*end && end[1])
        return -1;
    *bytes = (size_t) (x * fac);
    return 0;
}

static void parse_args(int argc, char **argv, Args *args)
{
    *args = (const Args){
        .cg_fs_dir = "/sys/fs/cgroup",
        .delim = ';',
        .interval_ms = 100
    };
    // int, because getopt() gives -1 at the end, which a char can't hold
    // where it is unsigned, e.g., on ARM
    int c = 0;
    bool found_cmd = false;
    // '-' prefix: no reordering of arguments, non-option (positional) arguments are
    // returned as argument to the 1 option
    // ':': preceding option takes a mandatory argument
    while (!found_cmd && (c = getopt(argc, argv, "-c:d:i:m:M:s:htZ")) != -1) {
        switch (c) {
            case 1:
                found_cmd = true;
                break;
            case '?':
                fprintf(stderr, "unexpected option character: %c\n", optopt);
                usage_hint(stderr, argv[0]);
                exit(125);
                break;
            case 'm':
                args->cg_fs_dir = optarg;
                break;
            case 'c':
                args->cg_dir = optarg;
                break;
            case 'd':
                args->delim = *optarg;
                break;
            case 'h':
                help(stdout, argv[0]);
                exit(0);
                break;
            case 't':
                args->machine_readable = true;
                break;
            case 's':
                args->samples = optarg;
                break;
            case 'i':
                args->interval_ms = atol(optarg);
                if (args->interval_ms < 1) {
                    fprintf(stderr, "invalid sampling interval: %s\n", optarg);
                    exit(125);
                }
                break;
            case 'M':
                args->mem_limit = optarg;
                if (parse_mem(optarg, &args->mem_limit_bytes)) {
                    fprintf(stderr, "invalid memory limit: %s\n", optarg);
                    exit(125);
                }
                break;
            case 'Z':
                args->no_systemd_run = true;
        }
    }
    if (!found_cmd) {
        if (argc > 1)
            fprintf(stderr, "No command found in the arguments!\n");
        usage_hint(stderr, argv[0]);
        exit(125);
    }
    int i = optind - 1;
    args->child_argv = argv + i;
}

static int current_cgroup(char *s, size_t n, const char **dest, size_t *dest_n)
{
    int fd = open("/proc/self/cgroup", O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Can't open /proc/self/cgroup (errno: %d - %s)\n",
                errno, strerrorname_np(errno));
        return -1;
    }
    ssize_t l = read(fd, s, n - 1);
    if (l == -1) {
        fprintf(stderr, "Can't read /proc/self/cgroup (errno: %d - %s)\n",
                errno, strerrorname_np(errno));
        close(fd);
        return -2;
    }
    if (!l) {
        fprintf(stderr, "Read nothing from /proc/self/cgroup\n");
        return -3;
    }
    --l;
    s[l] = 0;
    if (close(fd) == -1) {
        fprintf(stderr, "Can't close /proc/self/cgroup (errno: %d - %s)\n",
                errno, strerrorname_np(errno));
        return -4;
    }
    char *p = memchr(s, '/', l);
    if (!p) {
        fprintf(stderr, "Cgroup %s doesn't contain a slash\n", s);
        return -5;
    }
    *dest = p;
    *dest_n = l - (p - s);
    return 0;
}

static int search_cgroup_dir(const char *cg_fs_dir, char **dir)
{
    size_t a = strlen(cg_fs_dir);
    char buf[1024] = {0};
    const char *d = 0;
    size_t n = 0;
    if (current_cgroup(buf, sizeof buf, &d, &n))
        return -1;
    char *p = memmem(d, n, ".service", 8);
    if (!p) {
        return 1;
    }
    p[8] = 0;
    size_t b = (p + 8 - d);
    size_t m = a + b + 12 + 1;
    char *t = calloc(m, 1);
    if (!t) {
        fprintf(stderr, "Couldn't alloc cgroup directory string\n");
        return -2;
    }
    memcpy(mempcpy(mempcpy(t, cg_fs_dir, a), d, b), "/cgmt-XXXXXX", 12);
    *dir = t;
    return 0;
}

static int add_memory_controller(int cg_fd, const Args *args)
{
    const char fn[] = "cgroup.subtree_control";
    int fd = openat(cg_fd, fn, O_WRONLY);
    if (fd == -1) {
        fprintf(stderr, "Can't open %s/%s (errno: %d - %s)\n", args->cg_dir, fn,
                errno, strerrorname_np(errno));
        return -1;
    }
    const char s[] = "+memory";
    ssize_t l = write(fd, s, sizeof s - 1);
    if (l == -1) {
        fprintf(stderr, "Write to %s/%s failed (errno: %d - %s)\n", args->cg_dir, fn,
                errno, strerrorname_np(errno));
        return -2;
    }
    if (l != sizeof s - 1) {
        fprintf(stderr, "Partial write to %s/%s\n", args->cg_dir, fn);
        return -3;
    }
    if (close(fd) == -1) {
        fprintf(stderr, "Can't close %s/%s (errno: %d - %s)\n", args->cg_dir, fn,
                errno, strerrorname_np(errno));
        return -4;
    }
    return 0;
}

static int check_cgroupfs_memory(const Args *args, const char *filename, const char *verb, const char *hint)
{
    char *fn = 0;
    if (asprintf(&fn, "%s/%s", args->cg_fs_dir, filename) == -1) {
        fprintf(stderr, "can't allocate %s string\n", filename);
        return -1;
    }
    int r = 0;
    int fd = open(fn, O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Can't open %s (errno: %d - %s)\n",
                fn, errno, strerrorname_np(errno));
        r = -2;
        goto error;
    }
    char buf[1024] = {0};
    ssize_t l = read(fd, buf, sizeof buf - 1);
    if (l == -1) {
        fprintf(stderr, "Can't read %s (errno: %d - %s)\n",
                fn, errno, strerrorname_np(errno));
        r = -3;
        goto error;
    }
    const char *p = strstr(buf, "memory");
    // the list of controllers ends with a newline
    if (!p || !(p[6] == 0 || p[6] == ' ' || p[6] == '\n')) {
        fprintf(stderr, "Cgroup memory controller isn't %s in %s"
                " (%s)\n", verb, fn, hint);
        r = -23;
        goto error;
    }
error:
    if (fd != -1)
        close(fd);
    free(fn);
    return r;
}

static int check_cgroupfs(const Args *args)
{
    int r = check_cgroupfs_memory(args, "cgroup.controllers", "available", "perhaps disable via kernel parameter?");
    if (r)
        return r;
    r = check_cgroupfs_memory(args, "cgroup.subtree_control", "enabled", "systemd should enable it, by default");
    return r;
}

static int write_cg_file(int cg_fd, const char *fn, const char *s, const Args *args)
{
    int fd = openat(cg_fd, fn, O_WRONLY);
    if (fd == -1) {
        fprintf(stderr, "Can't open %s/leaf/%s (errno: %d - %s)\n", args->cg_dir, fn,
                errno, strerrorname_np(errno));
        return -1;
    }
    size_t n = strlen(s);
    ssize_t l = write(fd, s, n);
    if (l == -1 || (size_t) l != n) {
        fprintf(stderr, "Write to %s/leaf/%s failed (errno: %d - %s)\n", args->cg_dir, fn,
                errno, strerrorname_np(errno));
        close(fd);
        return -2;
    }
    if (close(fd) == -1) {
        fprintf(stderr, "Can't close %s/leaf/%s (errno: %d - %s)\n", args->cg_dir, fn,
                errno, strerrorname_np(errno));
        return -3;
    }
    return 0;
}

// sets the hard memory limit of the leaf cgroup and disables swap, as workload
// managers such as SLURM do with their jobs, e.g., with ConstrainRAMSpace=yes
// and AllowedSwapSpace=0
static int set_mem_limit(int cg_fd, const Args *args)
{
    char buf[32] = {0};
    snprintf(buf, sizeof buf, "%zu", args->mem_limit_bytes);
    if (write_cg_file(cg_fd, "memory.max", buf, args))
        return -1;
    // memory.swap.max is missing without swap accounting, e.g., swapaccount=0
    if (faccessat(cg_fd, "memory.swap.max", F_OK, 0) == 0) {
        if (write_cg_file(cg_fd, "memory.swap.max", "0", args))
            return -2;
    } else
        fprintf(stderr, "Warning: no memory.swap.max in %s/leaf, swap can't be disabled\n",
                args->cg_dir);
    return 0;
}

static int setup_cgroup(const Args *args, int *ofd)
{
    if (!mkdtemp(args->cg_dir)) {
        fprintf(stderr, "Can't create %s (errno: %d - %s)\n", args->cg_dir,
                errno, strerrorname_np(errno));
        return -1;
    }
    int fd = open(args->cg_dir, O_PATH);
    if (fd == -1) {
        fprintf(stderr, "Can't open %s (errno: %d - %s)\n", args->cg_dir,
                errno, strerrorname_np(errno));
        return -2;
    }

    // otherwise, without the nested setup we can't add a process to the parent cgroup
    // because we also need to write its cgroup.subtree_control file Cgroup v2
    // disallows doing both (yields EBUSY) - cf. https://unix.stackexchange.com/a/713343/1131
    const char fn[] = "leaf";
    if (mkdirat(fd, fn, 0700) == -1) {
        fprintf(stderr, "Can't mdir %s/%s (errno: %d - %s)\n", args->cg_dir, fn,
                errno, strerrorname_np(errno));
        return -3;
    }

    if (add_memory_controller(fd, args))
         return -3;

    int cg_fd = openat(fd, fn, O_PATH);
    if (fd == -1) {
        fprintf(stderr, "Can't open %s/%s (errno: %d - %s)\n", args->cg_dir, fn,
                errno, strerrorname_np(errno));
        return -4;
    }

    if (args->mem_limit && set_mem_limit(cg_fd, args))
        return -5;

    *ofd = fd;
    return cg_fd;
}

struct Result {
  struct timeval  child_user;
  struct timeval  child_sys;
  struct timespec child_wall;

  size_t child_rss_highwater;
  size_t cg_rss_highwater;
  size_t oom_kills;
};
typedef struct Result Result;


static void print_timeval(FILE *o, const struct timeval *t)
{
    fprintf(o, "%7.3f s",
            (double) t->tv_sec + (double)t->tv_usec / (1000.0 * 1000.0));
}

static void print_timeval_m(FILE *o, const struct timeval *t)
{
    fprintf(o, "%g",
            (double) t->tv_sec + (double)t->tv_usec / (1000.0 * 1000.0));
}

static void print_timespec(FILE *o, const struct timespec *t)
{
    fprintf(o, "%7.3f s",
            (double) t->tv_sec + (double)t->tv_nsec / 1000000000.0);
}

static void print_timespec_m(FILE *o, const struct timespec *t)
{
    fprintf(o, "%g",
            (double) t->tv_sec + (double)t->tv_nsec / 1000000000.0);
}

static int pretty_print(FILE *o, const Args *args, const Result *res)
{
    fputc('\n', o);
    fputs("user: ", o);
    print_timeval(o, &res->child_user);
    fputc('\n', o);
    fputs("sys:  ", o);
    print_timeval(o, &res->child_sys);
    fputc('\n', o);
    fputs("wall: ", o);
    print_timespec(o, &res->child_wall);
    fputc('\n', o);
    fprintf(o, "child_RSS_high: %10zu KiB\n",
            res->child_rss_highwater / 1024
           );
    fprintf(o, "group_mem_high: %10zu KiB\n",
            res->cg_rss_highwater / 1024
           );
    if (args->mem_limit)
        fprintf(o, "oom_kill:       %10zu\n", res->oom_kills);
    return 0;
}

static int machine_print(FILE *o, const Args *args, const Result *res)
{
    print_timeval_m(o, &res->child_user);
    fputc(args->delim, o);
    print_timeval_m(o, &res->child_sys);
    fputc(args->delim, o);
    print_timespec_m(o, &res->child_wall);
    fputc(args->delim, o);
    fprintf(o, "%zu%c%zu",
            res->child_rss_highwater / 1024,
            args->delim,
            res->cg_rss_highwater    / 1024
           );
    if (args->mem_limit)
        fprintf(o, "%c%zu", args->delim, res->oom_kills);
    fputc('\n', o);
    return 0;
}

static void print_result(FILE *o, const Args *args, const Result *res)
{
    if (args->machine_readable)
        machine_print(o, args, res);
    else
        pretty_print(o, args, res);
}


static struct timespec sub_ts(const struct timespec *a, const struct timespec *b)
{
    struct timespec r = {
        .tv_sec  = a->tv_sec  - b->tv_sec,
        .tv_nsec = a->tv_nsec - b->tv_nsec
    };
    if (r.tv_nsec < 0) {
        --r.tv_sec;
        r.tv_nsec += 1000000000l;
    }
    return r;
}

static int read_cg_rss_high(int cg_fd, size_t *rss)
{
    int fd = openat(cg_fd, "memory.peak", O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Can't open memory.peak (requires Kernel 5.19 or later): %d - %s\n",
                errno, strerrorname_np(errno));
        return -1;
    }
    char b[21] = {0};
    ssize_t l = read(fd, b, sizeof b - 1);
    if (l == -1) {
        fprintf(stderr, "Can't read memory.peak: %d - %s\n", errno, strerrorname_np(errno));
        close(fd);
        return -2;
    }
    if (close(fd) == -1) {
        fprintf(stderr, "Can't close memory.peak: %d - %s\n", errno, strerrorname_np(errno));
        return -3;
    }
    *rss = atol(b);
    return 0;
}

// current memory of the cgroup, in bytes
static int read_cg_current(int cg_fd, size_t *mem)
{
    int fd = openat(cg_fd, "memory.current", O_RDONLY);
    if (fd == -1)
        return -1;
    char b[21] = {0};
    ssize_t l = read(fd, b, sizeof b - 1);
    close(fd);
    if (l == -1)
        return -2;
    *mem = strtoul(b, 0, 10);
    return 0;
}

// current anonymous memory of the cgroup, in bytes, i.e., the memory of its
// processes not backed by files, which excludes the page cache, from the
// line 'anon' of memory.stat
static int read_cg_anon(int cg_fd, size_t *anon)
{
    int fd = openat(cg_fd, "memory.stat", O_RDONLY);
    if (fd == -1)
        return -1;
    char b[8192] = {0};
    ssize_t l = read(fd, b, sizeof b - 1);
    close(fd);
    if (l <= 0)
        return -2;
    // the line 'anon' is the first one of memory.stat, but it is searched
    // for at the start of any line
    for (const char *p = b; p && *p; p = strchr(p, '\n'), p = p ? p + 1 : 0) {
        if (!strncmp(p, "anon ", 5)) {
            *anon = strtoul(p + 5, 0, 10);
            return 0;
        }
    }
    return -3;
}

// writes a sample of the memory of the cgroup to 'samples', together with
// its anonymous memory
static void write_sample(FILE *samples, int cg_fd)
{
    size_t mem = 0, anon = 0;
    if (read_cg_current(cg_fd, &mem) || read_cg_anon(cg_fd, &anon))
        return;
    struct timespec now = {0};
    clock_gettime(CLOCK_REALTIME, &now);
    fprintf(samples, "%.3f,%zu,%zu\n",
            (double) now.tv_sec + (double) now.tv_nsec / 1000000000.0,
            mem / 1024, anon / 1024);
}

// number of processes of the cgroup killed by the OOM killer
static int read_cg_oom_kills(int cg_fd, size_t *kills)
{
    int fd = openat(cg_fd, "memory.events", O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Can't open memory.events: %d - %s\n", errno, strerrorname_np(errno));
        return -1;
    }
    char b[1024] = {0};
    ssize_t l = read(fd, b, sizeof b - 1);
    close(fd);
    if (l == -1) {
        fprintf(stderr, "Can't read memory.events: %d - %s\n", errno, strerrorname_np(errno));
        return -2;
    }
    const char *p = strstr(b, "oom_kill ");
    *kills = p ? strtoul(p + 9, 0, 10) : 0;
    return 0;
}

// number of processes in the cgroup, from cgroup.procs, killing each of them
// when 'kill_them' is true, or -1 when it can't be read
static int cg_procs(int cg_fd, bool kill_them)
{
    int fd = openat(cg_fd, "cgroup.procs", O_RDONLY);
    if (fd == -1)
        return -1;
    int n = 0;
    char b[4096];
    char pid[21] = {0};
    size_t k = 0;
    ssize_t l;
    while ((l = read(fd, b, sizeof b)) > 0) {
        for (ssize_t i = 0; i < l; ++i) {
            if (b[i] != '\n') {
                if (k < sizeof pid - 1)
                    pid[k++] = b[i];
                continue;
            }
            pid[k] = 0;
            k = 0;
            ++n;
            if (kill_them)
                kill((pid_t) strtol(pid, 0, 10), SIGKILL);
        }
    }
    close(fd);
    if (l == -1)
        return -1;
    return n;
}

// kills the processes left in the cgroup after the command ended, e.g.,
// forked workers of a command killed for exceeding the memory limit, which
// would otherwise keep running, and holding the output of cgmemtime open,
// and the cgroup could not be removed. it writes to cgroup.kill (Linux 5.14
// or later), which kills them all at once, also those that they fork, or
// otherwise kills them one by one from cgroup.procs, and waits up to 10 s
// for the cgroup to be empty. it gives the number of processes left, or -1
// when they can't be known
static int kill_cg_leftovers(const Args *args, int cg_fd)
{
    int n = cg_procs(cg_fd, false);
    if (n <= 0)
        return n;
    fprintf(stderr, "Killing %d process%s left in %s/leaf after the command ended\n",
            n, n == 1 ? "" : "es", args->cg_dir);
    bool cgkill = false;
    int fd = openat(cg_fd, "cgroup.kill", O_WRONLY);
    if (fd != -1) {
        cgkill = write(fd, "1", 1) == 1;
        close(fd);
    }
    struct timespec interval = { .tv_sec = 0, .tv_nsec = 100000000l };
    for (int i = 0; i < 100; ++i) {
        int m = cg_procs(cg_fd, !cgkill);
        if (m == 0)
            return n;
        if (m < 0)
            break;
        nanosleep(&interval, 0);
    }
    fprintf(stderr, "Processes still left in %s/leaf\n", args->cg_dir);
    return n;
}

static int execute(const Args *args, int cg_fd, Result *res)
{
    int pid_fd = -1;
    struct clone_args ca = {
        .flags = CLONE_INTO_CGROUP | CLONE_PIDFD | CLONE_VFORK,
        .pidfd = (uintptr_t) &pid_fd,
        // otherwise, there is nothing to wait on
        .exit_signal = SIGCHLD,
        .cgroup = cg_fd
    };
    struct timespec start = {0};
    if (clock_gettime(CLOCK_MONOTONIC, &start) == -1) {
        fprintf(stderr, "gettime start failed: %d - %s\n", errno, strerrorname_np(errno));
        return -1;
    }
    int pid = clone3(&ca);
    if (pid == -1) {
        fprintf(stderr, "clone3 failed: %d - %s\n", errno, strerrorname_np(errno));
        return -2;
    }
    if (!pid) { // child
        execvp(*args->child_argv, args->child_argv);
        fprintf(stderr, "execvp failed: %d - %s\n", errno, strerrorname_np(errno));
        if (errno == ENOENT)
            _exit(127);
        _exit(126);
    } else { // parent
        struct sigaction sa = {
            .sa_handler = SIG_IGN
        };
        // otherwise, Ctrl+C/+] also kill cgmemtime before it has a chance printing its summary
        if (sigaction(SIGINT, &sa, 0) == -1)
            fprintf(stderr, "failed to ignore SIGINT: %d - %s\n", errno, strerrorname_np(errno));
        if (sigaction(SIGQUIT, &sa, 0) == -1)
            fprintf(stderr, "failed to ignore SIGQUIT: %d - %s\n", errno, strerrorname_np(errno));

        FILE *samples = 0;
        if (args->samples && !(samples = fopen(args->samples, "w")))
            fprintf(stderr, "Can't open %s (errno: %d - %s)\n", args->samples,
                    errno, strerrorname_np(errno));
        struct timespec interval = {
            .tv_sec  = args->interval_ms / 1000,
            .tv_nsec = (args->interval_ms % 1000) * 1000000l
        };
        siginfo_t info = {0};
        struct rusage usg = {0};
        // with samples, poll the child between samples instead of blocking
        for (;;) {
            if (samples)
                write_sample(samples, cg_fd);
            info = (siginfo_t){0};
            usg = (struct rusage){0};
            if (raw_waitid(P_PIDFD, pid_fd, &info,
                           WEXITED | (samples ? WNOHANG : 0), &usg) == -1) {
                if (errno == EINTR)
                    continue;
                fprintf(stderr, "waitid failed: %d - %s\n", errno, strerrorname_np(errno));
                return -3;
            }
            if (info.si_pid != 0)
                break;
            nanosleep(&interval, 0);
        }
        if (samples)
            fclose(samples);
        struct timespec stop = {0};
        if (clock_gettime(CLOCK_MONOTONIC, &stop) == -1) {
            fprintf(stderr, "gettime stop failed: %d - %s\n", errno, strerrorname_np(errno));
            return -4;
        }
        res->child_user = usg.ru_utime;
        res->child_sys  = usg.ru_stime;
        res->child_wall = sub_ts(&stop, &start);
        res->child_rss_highwater = usg.ru_maxrss * 1024;
        if (read_cg_rss_high(cg_fd, &res->cg_rss_highwater))
            return -5;
        if (args->mem_limit && read_cg_oom_kills(cg_fd, &res->oom_kills))
            return -6;
        if (info.si_code == CLD_KILLED)
            return 128 + info.si_status;
        return info.si_status;
    }
    return 0;
}

static int teardown_cgroup(const Args *args, int cgp_fd, int cg_fd)
{
    const char fn[] = "leaf";
    if (close(cg_fd) == -1) {
        fprintf(stderr, "Can't close %s/%s (errno: %d - %s)\n", args->cg_dir, fn,
                errno, strerrorname_np(errno));
        return -1;
    }
    if (unlinkat(cgp_fd, fn, AT_REMOVEDIR) == -1) {
        fprintf(stderr, "Can't remove %s/%s (errno: %d - %s)\n", args->cg_dir, fn,
                errno, strerrorname_np(errno));
        return -2;
    }
    if (close(cgp_fd) == -1) {
        fprintf(stderr, "Can't close %s (errno: %d - %s)\n", args->cg_dir,
                errno, strerrorname_np(errno));
        return -3;
    }
    if (rmdir(args->cg_dir) == -1) {
        fprintf(stderr, "Can't remove %s (errno: %d - %s)\n", args->cg_dir,
                errno, strerrorname_np(errno));
        return -4;
    }
    return 0;
}

// The reexec is done in cases where the user's session doesn't run
// under something like /sys/fs/cgroup/.../user@$UID.service/...,
// e.g. when being logged in via ssh.
// In that case the process in a cgroup such as /user.slice/user-1000.slice/session-46.scope
// for which user 1000 doesn't have write permissions on cgroup.procs.
// Since /user.slice/user-$UID.slice/user@$UID.service/cgroup.procs usually
// is set up with write permissions we need to move the process to it because
// of the Cgroup v2 'Cgroup delegation containment rules' (cf. cgroups(7)), especially:
//
// > The writer has write permission on the cgroup.procs file in the nearest
// > common ancestor of the source  and  destination cgroups.
//
// Because systemd-run runs the process under /user.slice/user-1000.slice/user@1000.service/app.slice/run-$RANDOM.scope (with privilleged help) we can follow that rule again (without extra privilleges).
//
// NB: A Wayland/X11 terminal already runs under .../user@UID.service/... such that
//     there is no need to reexec like that.
//
// See also: https://github.com/systemd/systemd/issues/3388
static int reexec_with_systemd_run(int argc, char **argv)
{
    char **new_argv = calloc(argc + 5 + 1, sizeof new_argv[0]);
    if (!new_argv) {
        fprintf(stderr, "Can't allocate new argv\n");
        return -1;
    }
    new_argv[0] = strdup("systemd-run");
    if (!new_argv[0]) {
        fprintf(stderr, "Can't allocate new argv[0]\n");
        return -2;
    }
    new_argv[1] = strdup("--user");
    if (!new_argv[1]) {
        fprintf(stderr, "Can't allocate new argv[1]\n");
        return -3;
    }
    new_argv[2] = strdup("--scope");
    if (!new_argv[2]) {
        fprintf(stderr, "Can't allocate new argv[2]\n");
        return -4;
    }
    new_argv[3] = strdup("--quiet");
    if (!new_argv[3]) {
        fprintf(stderr, "Can't allocate new argv[3]\n");
        return -5;
    }
    new_argv[4] = argv[0];
    new_argv[5] = strdup("-Z");
    if (!new_argv[5]) {
        fprintf(stderr, "Can't allocate new argv[5]\n");
        return -6;
    }
    for (int i = 1; i < argc; ++i)
        new_argv[5 + i] = argv[i];
    execvp(new_argv[0], new_argv);
    fprintf(stderr, "execvp failed: %d - %s\n", errno, strerrorname_np(errno));
    if (errno == ENOENT) {
        fprintf(stderr, "Couldn't find systemd-run\n");
        return -7;
    }
    fprintf(stderr, "Execing with stemd-run failed (errno: %d - %s)\n",
            errno, strerrorname_np(errno));
    return -8;
}

int main(int argc, char **argv)
{
    Args args = {0};
    parse_args(argc, argv, &args);

    if (check_cgroupfs(&args)) {
        return 124;
    }

    if (!args.cg_dir) {
        int r = search_cgroup_dir(args.cg_fs_dir, &args.cg_dir);
        if (r < 0)
            return 120;
        if (r == 1) {
            if (args.no_systemd_run) {
                fprintf(stderr, "Couldn't find user@$UID.service cgroup - cf. -c option\n");
                return 119;
            } else {
                if (reexec_with_systemd_run(argc, argv))
                    return 118;
            }
        }
    }

    int cgp_fd = -1;
    int cg_fd = setup_cgroup(&args, &cgp_fd);
    if (cg_fd < 0)
        return 123;

    Result res = {0};
    int r = execute(&args, cg_fd, &res);
    kill_cg_leftovers(&args, cg_fd);
    print_result(stdout, &args, &res);

    if (teardown_cgroup(&args, cgp_fd, cg_fd))
        return 122;

    if (r < 0)
        return 121;
    return r;
}
