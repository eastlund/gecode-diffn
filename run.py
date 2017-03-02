import glob, os, subprocess, sys, pdb

output_file = "data/output.txt"
bench_path = "/home/eagiz/Tools/minizinc-benchmarks/"
data_location = {'rect_packing_mznc2014': "data_square"}
benchmarks = [("rectangle-packing/", "rect_packing_mznc2014")]
implementations = ["sweep", "original"]



timeout = 1
#ex = "mzn2fzn -G newgecode filter.mzn ar_1_3.dzn"
mzn2fzn = "mzn2fzn -G newgecode"
fzn_executable = "/home/eagiz/Tools/gecode-5.0.0/tools/flatzinc/fzn-gecode"
fzn_flags = " -s -time "+str(timeout*1000)

def reset_data():
    if not os.path.exists("./data"):
        os.makedirs("./data")
    with open(output_file, "w") as f:
        f.write("")

def write_data(implementation, model, inp, output, stderr):
    with open(output_file, "a") as f:
        f.write(implementation + "\n")
        f.write(model+ "\n")
        f.write(inp+ "\n")
        f.write("stdout: " + output + "\n")
        f.write("stderr: " + stderr + "\n")
        f.write("\n")

def scenario_generator(benchmark):
    l = []
    data_files = glob.glob(bench_path+benchmark[0]+data_location[benchmark[1]]+"/*.dzn")
    for data_file in data_files:
        command = mzn2fzn+" "+bench_path+benchmark[0]+benchmark[1]+".mzn"+" "+data_file
        l.append(command)
    return l

# probe a scenario for execution (generate new .fzn file)
def prepare(scenario):
    p = subprocess.Popen(scenario, shell=True, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()

def run_scenario(benchmark, scenario):
    data_file = scenario.rpartition('/')[2]
    benchmark_full = bench_path+benchmark[0]+benchmark[1]+".fzn"
    command = fzn_executable+fzn_flags+" "+benchmark_full
    p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    out, err = p.communicate()
    #print out
    #print err
    return (str(data_file), out, err)

def parse_output(sat_problem, implementation, outputs, d):
    for data_file, output, err in outputs:
        if not err:
            if len(output) == 0:
                optimal = False
                time = "t/o"
                final_sol = "none"
            else:
                try:
                    time_split = output.split("%%  runtime:")
                    #print output
                    final_time = time_split[-1]

                    solution_split = output.rpartition('----------')[0] #output.split("Sol:")
                    if sat_problem:
                        solution_split = solution_split[1]
                    else:
                        solution_split = solution_split[-1]
                    # final_sol = solution_split[:solution_split.find("\n")]
                    # optimal = False
                    # if "===" in final_time:
                    #     optimal = True
                    final_time
                    time = final_time.split('ms')[0][1:][:-1]
                    time = time.split('(')[1]
                    if sat_problem:
                        final_sol = "sat"
                        optimal = False
                except:
                    time = "t/o"
                    final_sol = "none"
                    optimal = False
        if not d[implementation].get(data_file):
            d[implementation][data_file] = {}
        d[implementation][data_file] = (optimal, time, final_sol)

def sort_alg_rect_packing(x):
    temp = x.rpartition('_')[0]
    n = temp.rpartition('rpp')[2]
    try:
        k = int(n)
        if "false" in x:
            return k
        else:
            return k*50
    except:
        return 0

sort_algs = {'rect_packing_mznc2014': sort_alg_rect_packing}

def latex_print(benchmark, d):
    print benchmark
    for implementation, results in sorted(d.iteritems(), key=lambda x: x[0]):
        print("--- {0} ---".format(implementation))
        for data_file in sorted(results.iteritems(), key=lambda x: sort_algs[benchmark](x[0])):
            print data_file[0]
            print("\\\\ \\texttt{{{0}}}".format(data_file[0]))
            #print("\\\\ \\texttt{{n={0}}}".format(n))
            (t, (optimal, time, sol)) = data_file
            if optimal:
                print("&\t\\textbf{%s} & \\textbf{%s}" % (sol, time))
            else:
                print("&\t{0} &\t {1}".format(sol, time))
            # for backend, (optimal, time, sol) in sorted(data_file.iteritems(), key=lambda x: x[0]):
            #     if optimal:
            #         print("&\t\\textbf{%s} & \\textbf{%s}" % (sol, time))
            #     else:
            #         print("&\t{0} &\t {1}".format(sol, time))
        print("\\\\")

def main():
    d = {}
    reset_data()
    print "hej"
    
        # command = "stdbuf -oL timeout {0} {1}{2}{3} {5}{6} --output-time".format(timeout, ex, bench_path, benchmark[0], model_path, model, data_path, data_file)
        # Pre-processing
    
        #command = "{0} {1}{2}{3} {5}{6} --output-time".format(ex, bench_path, benchmark[0], benchmark[1], data_file)
    for benchmark in benchmarks:
        scenarios = scenario_generator(benchmark)
        l = []
        for scenario in scenarios:
            #print scenario
            prepare(scenario)
            l.append(run_scenario(benchmark, scenario))
        d["sweep"] = {}
        parse_output(True, "sweep", l, d)
        print d
        latex_print(benchmark[1], d)

main()

