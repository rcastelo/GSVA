// macmemtime - high-water memory usage of a command including descendants
//              on macOS, mimicking the interface and output of cgmemtime
//              (https://github.com/gsauthof/cgmemtime), which is Linux-only.
//
// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: © 2026 the authors of the GSVA Bioconductor package
//
// It was written for the GSVA package to calibrate its memory model, after
// the interface and output of cgmemtime, which is distributed under the GNU
// General Public License version 3 or any later version, and is distributed
// under the same license, see the file COPYING.
//
// macOS has no cgroups, so the group high-water mark is obtained by sampling:
// every -i milliseconds the process tree rooted at the command is traversed
// (proc_listchildpids()) and the physical footprint (proc_pid_rusage(),
// ri_phys_footprint) of all its processes is summed up. The physical footprint
// is the metric shown by Activity Monitor and includes compressed and swapped
// memory. It excludes the page cache, but the anonymous pages that a forked
// process shares copy-on-write with its parent, and has not written to, are
// charged to both once the forked process reads them, e.g., the footprint of
// each of four workers forked by parallel::mclapply() that only read a vector
// of 2 GB includes the 2 GB. To count those pages once, as cgroups do on
// Linux, every process but the command counts at most its anonymous pages not
// shared with other processes, resident or compressed, from its memory regions
// (proc_pidinfo(), PROC_PIDREGIONPATHINFO), which include the pages that it
// has copied by writing to them. Pages that a worker still shares with its
// parent after the parent released them are then not counted.
// The child RSS high-water mark is the ru_maxrss reported by wait4(), i.e.,
// the maximum RSS among the command and its waited-for descendants.
//
// macOS has no hard memory limit for a process tree either, so the option -M
// emulates the one of a workload manager, such as SLURM, by killing (SIGKILL)
// the whole process tree as soon as a sample of its footprint exceeds LIMIT.
//
// A process whose parent ends is adopted by launchd and leaves the tree, e.g.,
// a forked worker of a command that was killed, so the processes seen in the
// tree are remembered, by their pid and start time, and those still running
// are counted, and killed by -M, also after leaving it, as cgroups do on
// Linux. Those still running when the command ends are killed, because they
// would otherwise keep running, and holding the output of macmemtime open.

#include <errno.h>
#include <libproc.h>
#include <sys/proc_info.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#define MAX_PIDS 4096

struct Args {
    char delim;
    int machine;
    long interval_ms;
    const char *mem_limit;
    uint64_t mem_limit_bytes;
    const char *samples;
    char **argv;
};
typedef struct Args Args;

#define SYNOPSIS "[OPTIONS] COMMAND [ARG...]"

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
            "\n"
            "Usage: %s " SYNOPSIS "\n"
            "\n"
            "Options:\n"
            "    -t            machine readable output (delimited columns)\n"
            "    -d CHAR       column delimiter (default: ';')\n"
            "    -i MSEC       sampling interval in milliseconds (default: 100)\n"
            "    -M LIMIT      memory limit, in bytes or with suffix K, M, G or T\n"
            "                  (powers of 1024), above which the process tree is\n"
            "                  killed, which adds to the output the number of kills\n"
            "    -s FILE       write each sample to FILE as a line with the time in\n"
            "                  seconds since the epoch and the memory in KiB, twice,\n"
            "                  as the memory and the anonymous memory of the output\n"
            "                  of cgmemtime, because the physical footprint already\n"
            "                  excludes the page cache\n"
            "    -h            show this help\n"
            "\n"
            "Columns (-t): user,sys,wall,child_RSS_high,group_mem_high\n"
            "(times in seconds, memory in KiB)\n",
            argv0, argv0);
}

// parses a memory size given in bytes or with suffix K, M, G or T
static int parse_mem(const char *s, uint64_t *bytes)
{
    char *end = 0;
    double x = strtod(s, &end);
    if (end == s || x <= 0)
        return -1;
    double fac = 1;
    switch (*end) {
        case 0:                           break;
        case 'K': case 'k': fac = 1024.0; break;
        case 'M': case 'm': fac = 1024.0 * 1024; break;
        case 'G': case 'g': fac = 1024.0 * 1024 * 1024; break;
        case 'T': case 't': fac = 1024.0 * 1024 * 1024 * 1024; break;
        default: return -1;
    }
    if (*end && end[1])
        return -1;
    *bytes = (uint64_t) (x * fac);
    return 0;
}

static int parse_args(int argc, char **argv, Args *args)
{
    *args = (Args) { .delim = ';', .interval_ms = 100 };
    int c;
    while ((c = getopt(argc, argv, "+d:hi:M:s:t")) != -1) {
        switch (c) {
            case 'd': args->delim = *optarg; break;
            case 'h': help(stdout, argv[0]); exit(0);
            case 'i': args->interval_ms = atol(optarg); break;
            case 'M':
                args->mem_limit = optarg;
                if (parse_mem(optarg, &args->mem_limit_bytes)) {
                    fprintf(stderr, "invalid memory limit: %s\n", optarg);
                    return -1;
                }
                break;
            case 's': args->samples = optarg; break;
            case 't': args->machine = 1; break;
            default:  help(stderr, argv[0]); return -1;
        }
    }
    if (optind >= argc) {
        if (argc > 1)
            fprintf(stderr, "No command found in the arguments!\n");
        usage_hint(stderr, argv[0]);
        return -1;
    }
    if (args->interval_ms < 1) {
        fprintf(stderr, "Sampling interval must be a positive integer\n");
        return -1;
    }
    args->argv = argv + optind;
    return 0;
}

// stores in 'tree' pid and all its descendants, returning their number
static int process_tree(pid_t pid, pid_t *tree)
{
    pid_t stack[MAX_PIDS];
    pid_t kids[MAX_PIDS];
    int n = 0, nt = 0;

    stack[n++] = pid;
    while (n > 0 && nt < MAX_PIDS) {
        pid_t p = stack[--n];
        tree[nt++] = p;
        int nb = proc_listchildpids(p, kids, sizeof kids);
        if (nb <= 0)
            continue;
        // depending on the macOS version, the return value is either the
        // number of pids or the number of bytes
        int nk = nb;
        if (nk > MAX_PIDS)
            nk = nb / (int) sizeof(pid_t);
        for (int i = 0; i < nk && n < MAX_PIDS; i++)
            if (kids[i] > 0)
                stack[n++] = kids[i];
    }
    return nt;
}

// a process seen in the tree, identified by its pid and start time, because
// pids are reused
struct Seen {
    pid_t pid;
    uint64_t start_sec, start_usec;
};
typedef struct Seen Seen;

// start time of the process 'pid', returning -1 when it is not running
static int start_time(pid_t pid, uint64_t *sec, uint64_t *usec)
{
    struct proc_bsdinfo bi;
    if (proc_pidinfo(pid, PROC_PIDTBSDINFO, 0, &bi, sizeof bi) != sizeof bi)
        return -1;
    *sec = bi.pbi_start_tvsec;
    *usec = bi.pbi_start_tvusec;
    return 0;
}

// whether the process seen 's' is still running
static int still_running(const Seen *s)
{
    uint64_t sec, usec;
    return start_time(s->pid, &sec, &usec) == 0 && sec == s->start_sec &&
           usec == s->start_usec;
}

// adds to 'seen', of 'ns' processes, those of 'tree' not in it, after
// removing those no longer running, returning the new number of processes
// seen
static int add_seen(Seen *seen, int ns, const pid_t *tree, int nt)
{
    int k = 0;
    for (int j = 0; j < ns; j++)
        if (still_running(&seen[j]))
            seen[k++] = seen[j];
    ns = k;
    for (int i = 0; i < nt && ns < MAX_PIDS; i++) {
        Seen s = { .pid = tree[i] };
        if (start_time(s.pid, &s.start_sec, &s.start_usec))
            continue;
        int known = 0;
        for (int j = 0; j < ns && !known; j++)
            known = seen[j].pid == s.pid && seen[j].start_sec == s.start_sec &&
                    seen[j].start_usec == s.start_usec;
        if (!known)
            seen[ns++] = s;
    }
    return ns;
}

// anonymous pages (bytes) of the process 'pid', i.e., not backed by a file,
// that it does not share with other processes, either resident or compressed
static uint64_t private_anon(pid_t pid)
{
    static long pagesize = 0;
    if (!pagesize)
        pagesize = sysconf(_SC_PAGESIZE);
    uint64_t pages = 0, addr = 0;
    struct proc_regionwithpathinfo r;
    while (proc_pidinfo(pid, PROC_PIDREGIONPATHINFO, addr, &r, sizeof r) ==
           sizeof r) {
        if (!r.prp_vip.vip_path[0])
            pages += r.prp_prinfo.pri_private_pages_resident +
                     r.prp_prinfo.pri_pages_swapped_out;
        addr = r.prp_prinfo.pri_address + r.prp_prinfo.pri_size;
    }
    return pages * (uint64_t) pagesize;
}

// sum of the physical footprint (bytes) of the processes seen that are still
// running, counting the pages that they share copy-on-write only in the
// command 'root'
static uint64_t seen_footprint(const Seen *seen, int ns, pid_t root)
{
    uint64_t total = 0;
    for (int i = 0; i < ns; i++) {
        struct rusage_info_v4 ri;
        if (still_running(&seen[i]) &&
            proc_pid_rusage(seen[i].pid, RUSAGE_INFO_V4,
                            (rusage_info_t *) &ri) == 0) {
            uint64_t m = ri.ri_phys_footprint;
            if (seen[i].pid != root) {
                uint64_t priv = private_anon(seen[i].pid);
                if (priv < m)
                    m = priv;
            }
            total += m;
        }
    }
    return total;
}

// kills (SIGKILL) the processes seen that are still running, other than
// 'except', returning their number
static int kill_seen(const Seen *seen, int ns, pid_t except)
{
    int n = 0;
    for (int i = 0; i < ns; i++)
        if (seen[i].pid != except && still_running(&seen[i]) &&
            kill(seen[i].pid, SIGKILL) == 0)
            n++;
    return n;
}

static double timespec_diff(const struct timespec *a, const struct timespec *b)
{
    return (double) (b->tv_sec - a->tv_sec) +
           (double) (b->tv_nsec - a->tv_nsec) / 1e9;
}

static double timeval_s(const struct timeval *t)
{
    return (double) t->tv_sec + (double) t->tv_usec / 1e6;
}

int main(int argc, char **argv)
{
    Args args;
    if (parse_args(argc, argv, &args))
        return 1;

    struct timespec start, stop;
    if (clock_gettime(CLOCK_MONOTONIC, &start)) {
        fprintf(stderr, "gettime start failed: %s\n", strerror(errno));
        return 1;
    }
    pid_t pid = fork();
    if (pid == -1) {
        fprintf(stderr, "fork failed: %s\n", strerror(errno));
        return 1;
    }
    if (pid == 0) {
        execvp(args.argv[0], args.argv);
        fprintf(stderr, "execvp failed: %s\n", strerror(errno));
        _exit(127);
    }
    signal(SIGINT, SIG_IGN);
    signal(SIGQUIT, SIG_IGN);

    FILE *samples = NULL;
    if (args.samples && !(samples = fopen(args.samples, "w")))
        fprintf(stderr, "Can't open %s: %s\n", args.samples, strerror(errno));

    uint64_t grp_high = 0;
    size_t kills = 0;
    pid_t tree[MAX_PIDS];
    static Seen seen[MAX_PIDS];
    int ns = 0;
    struct timespec ts = { .tv_sec = args.interval_ms / 1000,
                           .tv_nsec = (args.interval_ms % 1000) * 1000000L };
    int status;
    struct rusage ru;
    for (;;) {
        int nt = process_tree(pid, tree);
        ns = add_seen(seen, ns, tree, nt);
        uint64_t m = seen_footprint(seen, ns, pid);
        if (m > grp_high)
            grp_high = m;
        if (samples) {
            struct timespec now;
            clock_gettime(CLOCK_REALTIME, &now);
            fprintf(samples, "%.3f,%llu,%llu\n",
                    now.tv_sec + now.tv_nsec / 1e9,
                    (unsigned long long) (m / 1024),
                    (unsigned long long) (m / 1024));
        }
        if (args.mem_limit && m > args.mem_limit_bytes && !kills) {
            kill(pid, SIGKILL);
            kill_seen(seen, ns, pid);
            kills = 1;
            fprintf(stderr, "macmemtime: memory footprint %llu KiB above the "
                    "limit %s, process tree killed\n",
                    (unsigned long long) (m / 1024), args.mem_limit);
        }
        pid_t r = wait4(pid, &status, WNOHANG, &ru);
        if (r == pid)
            break;
        if (r == -1 && errno != EINTR) {
            fprintf(stderr, "wait4 failed: %s\n", strerror(errno));
            return 1;
        }
        nanosleep(&ts, NULL);
    }
    clock_gettime(CLOCK_MONOTONIC, &stop);
    if (samples)
        fclose(samples);
    // the processes left after the command ended, which left its tree when
    // it ended, except those forked since the last sample, which are missed
    int nleft = kill_seen(seen, ns, pid);
    if (nleft > 0)
        fprintf(stderr, "macmemtime: killed %d process%s left after the "
                "command ended\n", nleft, nleft == 1 ? "" : "es");

    double user = timeval_s(&ru.ru_utime);
    double sys  = timeval_s(&ru.ru_stime);
    double wall = timespec_diff(&start, &stop);
    size_t child_rss_kib = (size_t) ru.ru_maxrss / 1024; // bytes on macOS
    size_t grp_mem_kib   = (size_t) (grp_high / 1024);

    if (args.machine) {
        fprintf(stdout, "%g%c%g%c%g%c%zu%c%zu", user, args.delim, sys,
                args.delim, wall, args.delim, child_rss_kib, args.delim,
                grp_mem_kib);
        if (args.mem_limit)
            fprintf(stdout, "%c%zu", args.delim, kills);
        fputc('\n', stdout);
    } else {
        fprintf(stdout, "\nuser: %7.3f s\nsys:  %7.3f s\nwall: %7.3f s\n"
                "child_RSS_high: %10zu KiB\ngroup_mem_high: %10zu KiB\n",
                user, sys, wall, child_rss_kib, grp_mem_kib);
        if (args.mem_limit)
            fprintf(stdout, "oom_kill:       %10zu\n", kills);
    }

    if (WIFEXITED(status))
        return WEXITSTATUS(status);
    return 128 + WTERMSIG(status);
}
