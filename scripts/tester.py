#!/usr/bin/env python3

import os, sys
import warnings
from dataclasses import dataclass
from subprocess import PIPE, Popen
from queue import Empty, Queue
from threading import Thread
import time
from time import sleep
import re
import argparse
import math
from sys import platform

warnings.filterwarnings("ignore", category=RuntimeWarning)

BINARY = "pace"
TMP_FILE = "tmp.res"
VERBOSE = 0
TIMEOUT = 300

# FOLDER_NAME=os.path.join(os.getcwd(), "tiny")
# FOLDER_NAME=os.path.join(os.getcwd(), "medium")
# FOLDER_NAME=os.path.join(os.getcwd(), "large")
# FOLDER_NAME=os.path.join(os.getcwd(), "exact")
# FOLDER_NAME=os.path.join(os.getcwd(), "exact-hard")
# FOLDER_NAME=os.path.join(os.getcwd(), "cutwidth")
FOLDER_NAME=os.path.join(os.getcwd(), "heuristic")
# FOLDER_NAME=os.path.join(os.getcwd(), "heuristic-hard")

# ALG_PARAMS = ""

# heuristic
# ALG_PARAMS = "-bp-iters=-1 -time-limit=290"
# exact
# ALG_PARAMS = "-bp-iters=-1 -time-limit=290 -confidence=30"
# cutwidth
# ALG_PARAMS = "-bp-iters=-1 -time-limit=30"
# lite
# ALG_PARAMS = "-bp-iters=1 -max-docs-lb=32 -post-tune-int=0 -post-tune-all=0 -part-bp-iters=0 -leaf-interval=12"

# ALG_PARAMS = "-bp-iters=1 -max-docs-lb=32 -post-tune-int=0 -post-tune-all=0 -part-bp-iters=0"
ALG_PARAMS = "-bp-iters=1 -max-docs-lb=5120 -opt-interval-size=18 -post-tune-all=1 -post-tune-int=0 -part-bp-iters=0"


@dataclass
class TestResult:
  runtime: float = None
  memory: float = None
  score: float = None
  crossings: int = None
  status: str = None # OK, FAILED, UNKNOWN_OPT
  optimal: bool = False # provably optimal
  reached_lb: bool = False # reached lb given in the solution
  has_lb: bool = False
  num_edges: int = None
  num_nodes: int = None


def LOG(msg):
  sys.stdout.write("{}".format(msg))
  sys.stdout.flush()


def LOG_IF(verbose, msg):
  if verbose > 0:
    LOG(msg)


def enqueue_output(out, queue):
    for line in iter(out.readline, b""):
        queue.put(line)
    out.close()


def enqueue_outputs(proc, queue):
    # enqueue_output(proc.stdout, queue)
    enqueue_output(proc.stderr, queue)


def average(values):
  return sum(values) / float(len(values))


def percentile(values, p):
  data = sorted(values)
  k = (len(data) - 1) * (p / 100)
  f = math.floor(k)
  c = math.ceil(k)
  if f == c:
    return data[int(k)]
  d0 = data[int(f)] * (c - k)
  d1 = data[int(c)] * (k - f)
  return d0 + d1


def count_crossings(graph_file:str, res_file:str, verify_only:bool):
  # read the graph
  nA = None
  nB = None
  m = None
  edges = []
  with open(os.path.join(FOLDER_NAME, graph_file), "r") as f:
    for line in f.readlines():
      if "c " in line:
        continue
      if "p ocr " in line:
        assert nA is None and nB is None and m is None
        tmp = line.split()
        nA = int(tmp[2])
        nB = int(tmp[3])
        m = int(tmp[4])
      tmp = line.split()
      if len(tmp) == 2:
        u = int(tmp[0])
        v = int(tmp[1])
        assert 1 <= u and u <= nA
        assert 1 + nA <= v and v <= nB + nA
        edges.append((u - 1, v - 1 - nA))
  assert nA is not None and nB is not None and m is not None
  assert len(edges) == m

  # prepare data
  b2a = [[] for _ in range(nB)]
  for (u, v) in edges:
    assert 0 <= u and u < nA
    assert 0 <= v and v < nB
    b2a[v].append(u)

  # read the permutation
  if not os.path.isfile(res_file):
    res_file = os.path.join(FOLDER_NAME, res_file)
  if not os.path.isfile(res_file):
    return None, m, nA + nB
  order = []
  with open(res_file, "r") as f:
    for line in f.readlines():
      if "c " in line:
        continue
      tmp = line.split()
      if len(tmp) == 1:
        v = int(tmp[0])
        # if not (1 + nA <= v and v <= nB + nA):
        #   print(v)
        assert 1 + nA <= v and v <= nB + nA
        order.append(v)
  assert len(order) == nB

  if verify_only:
    return None, m, nA + nB

  # count crossings
  def addST(n, stree, pos):
    pos += n
    stree[pos] += 1
    while pos > 1:
      stree[pos>>1] = stree[pos] + stree[pos^1]
      pos >>= 1

  def sumST(n, stree, l, r):
    assert l <= r
    res = 0
    l += n
    r += n
    while l < r:
      if (l&1):
        res += stree[l]
        l += 1 
      if (r&1): 
        r -= 1
        res += stree[r]
      l >>= 1 
      r >>= 1
    return res

  n = nA
  stree = [0 for _ in range(2 * n)]
  num_crossings = 0
  for vv in order:
    v = vv - nA - 1
    assert 0 <= v and v < nB
    for u in b2a[v]:
      assert 0 <= u and u < n
      num_crossings += sumST(n, stree, u + 1, n)
    for u in b2a[v]:
      addST(n, stree, u)

  return num_crossings, m, nA + nB


def read_opt(file:str):
  filepath = os.path.join(FOLDER_NAME, file)
  opt_crossings = None
  lines = []
  with open(filepath, "r") as f:
    lines = f.readlines()
    for line in lines:
      if "c opt" in line:
        tmp = re.findall(r"[-+]?(?:[\d,.]*\d)", line)
        assert len(tmp) == 1
        opt_crossings = int(tmp[0])
  return (opt_crossings, lines)


def overwrite_opt(file:str, crossings:int):
  filepath = os.path.join(FOLDER_NAME, file)
  opt_crossings, lines = read_opt(file)

  if opt_crossings is None:
    LOG("  \033[96msetting opt in '{}' to {:,d}\033[0m\n".format(file, crossings))
  else:
    if opt_crossings <= crossings:
      return
    LOG("  \033[96mreducing opt in '{}' from {:,d} to {:,d}\033[0m\n".format(file, opt_crossings, crossings))

  with open(filepath, "w") as f:
    if opt_crossings is None:
      f.write("c opt {}\n".format(crossings))
    for line in lines:
      if "c opt" in line:
        assert opt_crossings is not None
        f.write("c opt {}\n".format(crossings))
      else:
        f.write(line)


def extract_log_data(log_file:str, result: TestResult):
  assert os.path.exists(log_file)
  with open(log_file, "r") as f:
    for line in f.readlines():
      if "score = " in line and "running time is" in line:
        ll = line[line.find("running time is"):]
        tmp = re.findall(r"[-+]?(?:[\d,.]*\d)", ll)
        if "score = unknown" in ll:
          assert len(tmp) == 2
          result.score = 0.0
          result.crossings = int(tmp[1].replace(",", ""))
        else:
          assert len(tmp) == 3
          result.score = float(tmp[1])
          result.crossings = int(tmp[2].replace(",", ""))

      if "delta to lower bound is 0 " in line:
        result.optimal = True
        result.has_lb = True

      if "strong2 lower bound" in line or "original lower bound" in line or "weak lower bound" in line:
        result.has_lb = True

      if "and |E| = " in line:
        tmp = re.findall(r"[-+]?(?:[\d,.]*\d)", line)
        assert len(tmp) >= 3
        result.num_nodes = int(tmp[-3].replace(",", "")) + int(tmp[-2].replace(",", ""))
        result.num_edges = int(tmp[-1].replace(",", ""))

      if "maximum resident set size" in line:
        tmp = re.findall(r"[-+]?(?:[\d,.]*\d)", line)
        assert len(tmp) == 1
        if "(KB)" in line:
          result.memory = float(tmp[0]) / (1024.0)
        else:
          result.memory = float(tmp[0]) / (1024.0 * 1024.0)


def process_pace(graph_file, args):
  res_file = TMP_FILE
  if os.path.exists(res_file):
      os.remove(res_file)

  if VERBOSE == 0:
    CMD = "./{} -verbose={} -i={} {} {} > {}".format(BINARY, VERBOSE, os.path.join(FOLDER_NAME, graph_file), ALG_PARAMS, args.extra_params, res_file)
  else:
    CMD = "./{} -verbose={} -o={} -i={} {} {}".format(BINARY, VERBOSE, res_file, os.path.join(FOLDER_NAME, graph_file), ALG_PARAMS, args.extra_params)

  # use timeout
  if TIMEOUT > 0:
    CMD = "timeout -s SIGTERM --preserve-status {} {}".format(TIMEOUT, CMD)

  # measure memory on macos
  if platform == "darwin":
    CMD = "/usr/bin/time -l {}".format(CMD)
  elif platform == "linux":
    CMD = "/usr/bin/time -f '%M maximum resident set size (KB)' {}".format(CMD)

  log_file = "log_0"
  with open(log_file, "w") as f:
      f.write("stderr for graph " + str(graph_file) + "\n")

  p = Popen(CMD, stdout=PIPE, stderr=PIPE, bufsize=1, shell=True)
  q = Queue()
  t = Thread(target=enqueue_outputs, args=(p, q))
  t.daemon = True  # thread dies with the program
  t.start()

  opt_crossings, _ = read_opt(graph_file)

  SLEEP = 0.03 if args.verbose > 0 else 0.005
  result = TestResult()
  start_time = time.time()

  # read line without blocking
  while p.poll() is None:
      sleep(SLEEP)
      try:
          line = q.get_nowait()
      except Empty:
          line = ""
      else:  # got line
          with open(log_file, "a") as f:
            f.write(line.decode())
            f.write("\n")
          last_line = line.strip()
          while not q.empty():
            line = q.get().strip()
            if len(line) > 0:
              last_line = line
            # logging
            with open(log_file, "a") as f:
              f.write(line.decode())
              f.write("\n")

          last_line = last_line.decode().strip()
          LOG_IF(args.verbose > 1, "  " + last_line + "\n")

  end_time = time.time()
  result.runtime = end_time - start_time

  if p.returncode != 0:
    LOG("  \033[91m%s\033[0m\n" % ('failed with exit code ' + str(p.returncode)))
    result.status = "failed"
    return result

  # extract data from the logs
  extract_log_data(log_file, result)

  # recompute the crossings, if needed
  if result.crossings is None:
    result.crossings, result.num_edges, result.num_nodes = count_crossings(graph_file, res_file, False)
    if opt_crossings is None or result.crossings is None:
      result.score = 0
    else:
      result.score = float(opt_crossings + 1) / float(result.crossings + 1)
      if opt_crossings < result.crossings:
        result.score = min(result.score, 0.99999)
  else:
    if result.runtime < 5.0:
      num_crossings, num_edges, result.num_nodes = count_crossings(graph_file, res_file, False)
      assert num_crossings == result.crossings
    else:
      count_crossings(graph_file, res_file, True)

  # check the optimal solution
  if result.crossings is not None:
    if args.overwrite and (opt_crossings is None or result.crossings < opt_crossings):
      overwrite_opt(graph_file, result.crossings)
    if opt_crossings is not None and opt_crossings >= result.crossings:
      result.reached_lb = True

  assert p.returncode == 0

  result.status = "ok"
  if result.optimal:
    LOG_IF(args.verbose > 0, " \033[92m{}\033[0m ".format("opt"))
  elif result.score == 1.0:
    LOG_IF(args.verbose > 0, " \033[92m{}\033[0m  ".format("ok"))
  else:
    LOG_IF(args.verbose > 0, " ok  ")
  LOG_IF(args.verbose > 0, "[score = {:.5f}".format(result.score))
  LOG_IF(args.verbose > 0, "; runtime = {:5.1f}sec".format(result.runtime))
  LOG_IF(args.verbose > 0 and result.memory is not None, "; memory = {:3.0f}MB".format(result.memory or 0))
  LOG_IF(args.verbose > 0, "; crossings = {:10,}".format(result.crossings or 0))
  LOG_IF(args.verbose > 0, "]")
  LOG_IF(args.verbose > 0, "   \033[90m|V| = {:6,}; |E| = {:6,}\033[0m".format(result.num_nodes, result.num_edges))
  LOG_IF(args.verbose > 0, "\n")

  return result


def process_sol(graph_file, sol_filename, args):
  if args.verbose > 0:
    LOG("processing '{}' with '{}' ...".format(graph_file, sol_filename))

  opt_crossings, _ = read_opt(graph_file)

  result = TestResult()
  result.status = "ok"
  result.runtime = 0
  result.has_lb = True
  result.crossings = count_crossings(graph_file, sol_filename, False)
  assert result.crossings is not None
  if opt_crossings is None or result.crossings <= opt_crossings:
    result.optimal = True

  if args.overwrite and (opt_crossings is None or result.crossings < opt_crossings):
    overwrite_opt(graph_file, result.crossings)

  min_crossings = result.crossings
  if opt_crossings is not None:
    min_crossings = min(min_crossings, opt_crossings)
  result.score = float(min_crossings + 1) / float(result.crossings + 1)

  return result


def aggregate(infos, total_runtime, args):
  num_ok = sum(1 for i in infos if i.status == "ok")
  num_failed = sum(1 for i in infos if i.status == "failed")
  num_opt = sum(1 for i in infos if i.optimal)
  num_lb = sum(1 for i in infos if i.reached_lb)
  num_has_lb = sum(1 for i in infos if i.has_lb)
  LOG("Total {} graphs [{}ok = {}{}, {}failed = {}{}, opt = {}/{}, lb = {}/{}]. Total runtime: {:.1f} sec\n".format(
    len(info),
    "\033[92m" if num_ok == len(info) else "", 
    num_ok,
    "\033[0m" if num_ok == len(info) else "",
    "\033[91m" if num_failed > 0 else "", 
    num_failed,
    "\033[0m" if num_failed > 0 else "",
    num_opt,
    num_has_lb,
    num_lb,
    num_ok,
    total_runtime,
  ))
  if num_ok == 0:
    return
  # score
  scores = [i.score for i in infos if i.status == "ok"]
  LOG("  {:10s}: {:15.5f} / {:15.5f} / {:15.5f} / {:15.5f} / {:15.5f}  (avg / p50 / p25 / p05 / min)\n".format(
    "score",
    average(scores),
    percentile(scores, 50),
    percentile(scores, 25),
    percentile(scores, 5),
    percentile(scores, 0),
  ))
  # crossings
  crossings = [i.crossings for i in infos if i.status == "ok"]
  LOG("  {:10s}: {:15,d} / {:15,d} / {:15,d} / {:15,d} / {:15,d}  (avg / p50 / p75 / p95 / max)\n".format(
    "crossings",
    int(average(crossings)),
    int(percentile(crossings, 50)),
    int(percentile(crossings, 75)),
    int(percentile(crossings, 95)),
    int(percentile(crossings, 100)),
  ))
  # runtime p50, p95, max
  runtimes = [i.runtime for i in infos if i.status == "ok"]
  LOG("  {:10s}: {:15.2f} / {:15.2f} / {:15.2f} / {:15.2f} / {:15.2f}  (avg / p50 / p75 / p95 / max)\n".format(
    "runtime",
    average(runtimes),
    percentile(runtimes, 50),
    percentile(runtimes, 75),
    percentile(runtimes, 95),
    percentile(runtimes, 100),
  ))


def parse_args():
  parser = argparse.ArgumentParser(description="PACE'24 Tester")
  parser.add_argument("-v", "--verbose", type=int, default=0, help="Verbosity level")
  parser.add_argument("-o", "--overwrite", action='store_true', help="Overwrite best results")
  parser.add_argument("-p", "--pattern", type=str, default="", help="Pattern for file matching")
  parser.add_argument("-s", "--solution", action='store_true', help="Read solution files")
  parser.add_argument("--extra-params", type=str, default="", help="Extra pace params")
  return parser.parse_args()

################################################################################
################################ MAIN ##########################################
################################################################################

args = parse_args()

# TODO: read the option from the command-line
# Make sure the folder exists
assert os.path.isdir(FOLDER_NAME)

# read graph instances
info = []
all_gr_files = (filename for _, _, files in os.walk(FOLDER_NAME) for filename in files if filename.endswith(".gr"))
# sorted_files = sorted(all_gr_files, key = lambda f : os.path.getsize(os.path.join(FOLDER_NAME, f)))
sorted_files = sorted(all_gr_files)
sol_files = {}
if args.solution:
  sol_files = {filename for _, _, files in os.walk(FOLDER_NAME) for filename in files if filename.endswith(".sol")}

LOG("processing {} instances from {}\n".format(len(sorted_files), FOLDER_NAME));
LOG("alg_params: {} {}\n".format(ALG_PARAMS, args.extra_params));

start_time = time.time()

# process every instance
for idx, filename in enumerate(sorted_files):
  if not os.path.isfile(os.path.join(FOLDER_NAME, filename)):
    continue
  if args.pattern != "" and not re.compile(args.pattern).match(filename):
    continue

  log_msg = "processing {:6s} ({:3d}/{}) ...".format(filename, idx + 1, len(sorted_files))
  if args.verbose > 1:
    LOG(log_msg + "\n")
  elif args.verbose > 0:
    LOG(log_msg)

  if args.solution:
    assert filename.endswith(".gr")
    sol_filename = filename[:-3] + ".sol"
    assert sol_filename in sol_files
    info.append(process_sol(filename, sol_filename, args))
  else:
    info.append(process_pace(filename, args))

end_time = time.time()
total_runtime = end_time - start_time

aggregate(info, total_runtime, args)
