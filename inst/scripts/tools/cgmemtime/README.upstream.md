cgmemtime measures the high-water RSS+CACHE memory usage of a process
and its descendant processes, using Linux Control Group v2.

To be able to do so it puts the process into its own
[cgroup][cg].

For example, process A allocates 10 [MiB][mib] and forks a child B that
allocates 20 MiB and that forks a child C that allocates 30 MiB.
All three processes share a time window where their allocations
result in a certain aggregate [RSS (resident set size)][rss] memory usage.

The question now is: How much memory is actually used as a result
of running A?

Answer: 60 MiB (assuming that A and B are still running when C
allocates its memory)

cgmemtime is the tool to answer such questions.

(It also measures the runtime.)

Last significant update: 2022-11-20

## Usage

Now you can use cgmemtime like this:

	$ ./cgmemtime ./testa x 10 20 30
    [..]
	child_RSS_high:      11808 KiB
	group_mem_high:      62164 KiB

Or to produce machine readable output:

    $ ./cgmemtime -t ./testa x 10 20 30

It also has some options (cf. `-h`).


## Dependencies

cgmemtime requires a Linux kernel with Control Group v2 support, including the
`memory.peak` feature, that means Linux 5.19 or newer.
For example, Fedora 36 and 37 work fine.

Other than that you need a C compiler, GNU make and the usual
development headers.

Enterprise Linux distributions might backport `memory.peak` to their
nominally 'frozen' and lower versioned kernels.
However, as of December, 2022, RHEL 9.1 (with
5.14.0-162.6.1.el9_1.x86_64) and Ubuntu 22.04.1 (with their 5.15
kernel) don't support it. FWIW, cgmemtime works on Ubuntu 22.04.1
when running their 6.0 'oem' Kernel.

By default, cgmemtime creates a temporary cgroup under the default systemd user
service cgroup, which doesn't require any special setup or root privileges.  If
you don't use systemd you can come up with a similar scheme and use the `-m`
and `-c` options.

See also older cgmemtime versions if you need one that supports
Linux Control Group v1.


## Compile

Just:

    $ make

Which creates `cgmemtime` and `testa`. `testa` is a small forking
allocation test program.


## Testing

You can run the test suite:

    $ bash test.sh


## FAQ

### Why is the accumulated RSS sometimes lower than the child RSS?

The thing is that the child number and the accumulated number
come from different subsystems in the kernel - which
have slightly different trade-offs/approximations of
the RSS of a process.

A simple test case:

    $ ./cgmemtime python -c 'import time; import os; print(os.getpid()); time.sleep(300)'
	35595
    [..]
	child_RSS_high:       9060 KiB
	group_mem_high:       3860 KiB

The first number is consistent to what GNU time (`/usr/bin/time`)
reports.  With both GNU time/cgmemtime, the number doesn't come
from the cgroups subsystem.

You can also approximate it with [something
like](http://unix.stackexchange.com/a/33388):

    $ awk '/Rss:/{ sum += $2 } END { print sum }' /proc/24131/smaps
    6388

The 2nd number comes from the cgroup subsystem. You can
approximate it via excluding some shared library mappings, e.g.:

    $ grep '^[0-9a-f]\|Rss:' /proc/24131/smaps | tr -d '\n' \
      | sed 's/ kB/ kB\n/g' | grep -v '.so' | sed 's/^.*Rss://' \
      | awk '{a+=$1} END {print a}'
    2760

Hypothesis: Linux cgroup doesn't account for the shared library
mappings and the effect is easy to demonstrate with Python
because it loads such a large amount of shared libraries.


## Contact

Don't hesitate to mail feedback (comments, questions, ...) to:

    Georg Sauthoff <mail@gms.tf>

## Licence

[GPLv3+][gpl]

## Accuracy

The reported high-water RSS+CACHE usage values are as accurate as the
[`memory.peak` value][peak] exported by the cgroup memory resource
controller.

The Control Group v2 documentation doesn't say much about its accuracy,
but probably similar caveats apply as to the similar cgroup v1 measure,
as detailed in kernel documentation:

> For efficiency, as other kernel components, memory cgroup uses
some optimization to avoid unnecessary cacheline false sharing.
usage_in_bytes is affected by the method and doesn't show 'exact'
value of memory(and swap) usage, it's an fuzz value for efficient
access. (Of course, when necessary, it's synchronized.) If you
want to know more exact memory usage, you should use
RSS+CACHE(+SWAP) value in memory.stat(see 5.2).

([Section 5.5](https://elixir.bootlin.com/linux/v5.19/source/Documentation/admin-guide/cgroup-v1/memory.rst#L617))

We can't use memory.stat because it does not include high-water
memory usage information and we don't want to poll it.

Doing some tests with e.g. `./testa` the reported values seem to
be exact enough, though.

## Limitations

The `memory.peak` measure reports the sum of RSS and CACHE
usage.  Thus, you can't measure the high-water RSS-without-CACHE
usage. In a program that does a lot of IO the CACHE part then
dominates the high-water RSS+CACHE value.

For example:

    $ cgmemtime dd if=test.img | dd of=out
    # vs.
    $ cgmemtime dd if=test.img of=out
    $ cgmemtime dd if=test.img of=out

(for a large test.img the 2nd command has a large RSS+CACHE
high-water value, i.e. 2 times the `test.img` size or so - while the 3rd command
yields a high-water usage of pretty much the test.img size, iff the input is
still part of the buffer cache ...)

Currently, I am not aware of a cgroup way to just derive the
RSS-only high-water mark.

FWIW, for some IO access patterns it makes sense to advise the kernel on how
it should cache file data (cf. `madvise()` and  `posix_fadvise()`).


## Implementation Details

Cgmemtime uses modern Linux specific syscalls, including ones for
which glibc lacks wrappers. Notably, is uses `clone3()` in order
to directly spawn the child process into the fresh measurement
cgroup and obtain its [PIDFD][pidfd]. While at it, it also
specifies the vfork flag to avoid superfluous COW setup, since
the child immediately execs a command.

For waiting on the child and obtaining some usage attributes the
extended Linux `waitid()` syscall is used. Besides obtaining the
resource usage, the parent waits on the child through a
[PIDFD][pidfd], because it's possible. Note that waiting on
the child's PID is as good here, since a terminating child stays
around as zombie after it terminates such that its PID can't be
recycled and a process only can wait on its child, anyways.


## Related Work

There are also other tools available which measure memory usage
of processes. One way to categorize them is a two-fold
classification: tools that use polling and tools that don't.

In that context - when you are only interested in the high-water
usage - polling is the inferior approach. As described in
previous sections, cgmemtine does not use polling. At the time of
writing, I am not aware of any other tool that uses Linux Control
Groups for memory measurements.

### Non-Polling based tools

#### Highwater measurements

- [GNU time][gtime] - uses something like `wait4()` or `waitpid()` and
  `getrusage()`, thus on systems where available it is able to
  display the high-water RSS usage of a single child process,
  when using the verbose mode.
- [tstime][tstime] - uses the taskstructs API of the Linux kernel to get
  the high-water RSS _and_ the highwater VMEM usage of a child.
  Does not follow descendant processes. Provides also a process
  monitor mode that displays stats for all exiting processes. But
  the taskstats API is kind of cumbersome to use and on current
  kernels only accessible as root.
- [dtmemtime][dtmt] - Dtrace Memtime, i.e. for Solaris built using Dtrace.
  One could probably implement something similar, on Linux, using
  bpftrace or even BPF directly, however, it would require root
  privileges

#### Snapshot measurements

- [smem][smem] - Tool written in Python that analyses proc files
  like `/proc/$$/smaps` and generates a memory usage report of
  one ore multiple processes for one point in time. It is designed
  to provide a system-wide view, but one can also filter processes
  (or even loaded libraries) by various criteria. Smem distributes
  shared memory between all dependent processes (the result is
  called proportional set size - PSS - of a process). It does not
  take swapped-out memory into account.

### Polling based tools

- [memtime][memtime] ([mirror][memtime2]) - Uses polling of `/proc/$PID/stat` to
  measure high-water RSS/VMEM usage of a child. It supports Linux
  and Solaris styles of `/proc`.  Polling is in general a
  sub-optimal solution (e.g. short-running processes are not
  accurately measured, it wastes resources etc.). memtime is not
  maintained and has 64 Bit issues (last release 2002).
- [tmem][tmem] - Polls `/proc/$PID/status`, thus has access to
  more detailed memory measures, e.g.  VmPeak, VmSize, VmLck,
  VmPin, VmHWM, VmRSS, VmData, VmStk, VmExe, VmLib, VMPTE and
  VMSwap.
- [memusg][memusg] - Python script that polls the VmSize values
  of a group of processes via the command `ps` and displays its
  high-water mark.  That means that it forks/execs `ps` and parses
  its output 10 times a second. For a given command line it creates
  a new session (via setsid()) and executes it in that session.
  Thus, children of the watched process are likely part of that
  session, too. Memusg then sums the VMSize value of each process
  of that session up and returns the maximum when the session
  leader exits. Note, that this method is not reliable, because
  child processes may still be alive after the session leader has
  exited and they may also create new sessions during their runtime,
  thus escaping the measurement via memusg.


[mib]: http://en.wikipedia.org/wiki/Mebibyte
[cg]: https://docs.kernel.org/admin-guide/cgroup-v2.html
[cgmem]: https://docs.kernel.org/admin-guide/cgroup-v2.html#memory
[systemd]: http://en.wikipedia.org/wiki/Systemd
[rss]: http://en.wikipedia.org/wiki/Resident_set_size
[gpl]: http://www.gnu.org/licenses/gpl.html
[gtime]: http://www.gnu.org/software/time/
[tstime]: https://bitbucket.org/gsauthof/tstime
[memtime]: http://www.update.uu.se/~johanb/memtime/
[memtime2]: https://github.com/phuseman/memtime
[tmem]: http://locklessinc.com/articles/memory_usage/
[smem]: http://www.selenic.com/smem/
[memusg]: https://github.com/jhclark/memusg
[peak]: https://elixir.bootlin.com/linux/v5.19/source/Documentation/admin-guide/cgroup-v2.rst#L1232
[dtmt]: https://github.com/gsauthof/utility/#dtmemtime
[pidfd]: https://lwn.net/Articles/801319/

