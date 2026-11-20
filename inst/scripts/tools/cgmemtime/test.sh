# 2012, Georg Sauthoff <mail@georg.so>
# GPLv3+

set -e
set -u

TIME=/usr/bin/time
CGMEMTIME=./cgmemtime
TESTA=./testa

test_linux()
{
  if uname -o | grep -q Linux ; then
    true
  else
    echo -n Not a Linux system
    return 1
  fi
  return 0
}

test_time()
{
  if $TIME --verbose true 2>&1 | grep -q 'Maximum resident set size' ; then
    true
  else
    echo -n '(GNU?)' time does not report maximum RSS
    return 1
  fi
  return 0
}

test_simple()
{
  for a in 4 5; do
  x=`$CGMEMTIME -t $TESTA x 10 2>&1 | grep ';' | cut -d';' -f $a`
  if [ $((x/1024)) -lt 10 -o $((x/1024)) -gt 11 ]; then
    echo -n "Not 10 MiB (x=$x, a=$a)"
    return 1
  fi
  done
  return 0
}

test_accum()
{
  o=`$CGMEMTIME -t $TESTA x 10 20 30 40 2>&1 | grep ';'`
 
  x=`echo $o | cut -d';' -f 4`
  if [ $((x/1024)) -lt 10 -o $((x/1024)) -gt 11 ]; then
    echo -n "Not 10 MiB (x=$x)"
    return 1
  fi
  x=`echo $o | cut -d';' -f 5`
  if [ $((x/1024)) -lt 100 -o $((x/1024)) -gt 102 ]; then
    echo -n "Not 100 MiB (x=$x)"
    return 1
  fi
  return 0
}

test_wall()
{
  x=`$CGMEMTIME -t sleep 3 2>&1 | grep ';' | cut -d';' -f 3 | cut -d'.' -f1`
  if [ "$x" -ne 3 ]; then
    echo -n "Not 3 seconds wall clock time (x=$x)"
    return 1
  fi
  return 0
}

test_notfound()
{
  $CGMEMTIME ./notfound 2>/dev/null
  x=$?
  if [ $x -ne 127 ]; then
    echo -n "Exit status is not 127 (x=$x)"
    return 1
  fi
  return 0
}

test_exitstatus()
{
  $CGMEMTIME false 2>/dev/null
  x=$?
  if [ $x -ne 1 ]; then
    echo -n "Exit status is not 1 (x=$x)"
    return 1
  fi
  return 0
}

tests="
test_linux
test_time
test_simple
test_accum
test_wall
test_notfound
test_exitstatus
"
failed=0

set +e
for i in $tests; do
  echo -n Executing $i ...
  $i
  if [ $? -eq 0 ]; then
    echo ' OK'
  else
    echo ' FAIL'
    failed=$((failed+1))
  fi
done
set -e

echo ================================================================================
echo ================================================================================
if [ $failed -eq 0 ]; then
  echo '=> Everything went just fine!'
  exit 0
else
  echo '=> '$failed' tests failed'
  exit 1
fi


