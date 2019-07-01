#! /usr/bin/env python3
import argparse
import random
from math import e
from shutil import copyfile

def main():
    parser = argparse.ArgumentParser(description="Run a simulation given the parameters in the par file.")
    parser.add_argument("-p", help=".par file to attach", dest="par", type=str, required=True)
    parser.set_defaults(func=parse_par)
    args = parser.parse_args()
    args.func(args)

def parse_par(args):
    filename = args.par
    params = {}
    with open(filename) as f:
        for line in f:
            parts = line.split(":")
            parts = [part.strip() for part in parts]
            print(parts)
            if True not in [bool(part) for part in parts]:
                continue
            assert len(parts) == 2, line + "Format must be 'parameter_name : parameter_value'"
            params[parts[0].strip()] = parts[1].strip()

    for param in ['n', 'ancestorAgeno', 'ancestorBgeno', 'ancestorAind', 'ancestorBind', 'ancestorAsnp', 
            'ancestorBsnp', 'trackancestry', 'theta', 'outputfilename', 'lambda', 'haploidoutput']:
        assert param in params.keys(), "Must give required parameter: " + param

    params["lambda"] = float(params["lambda"])
    params["theta"] = float(params["theta"])
    params["n"] = int(params["n"])
    params["trackancestry"] = False if params["trackancestry"] == 'False' else True
    params["haploidoutput"] = False if params["haploidoutput"] == 'False' else True

    print("Loaded par file, parameters are:")
    for param in ['n', 'ancestorAgeno', 'ancestorBgeno', 'ancestorAind', 'ancestorBind', 'ancestorAsnp', 
    'ancestorBsnp', 'trackancestry', 'theta', 'outputfilename', 'lambda', 'haploidoutput']:
        print(param, ":", params[param])


    run_simulation(params)

def readfile(filename):
    lines = []
    with open(filename) as f:
        for line in f:
            if line != None:
                lines.append(line.strip())
    return lines

def run_simulation(params):
    geno = {'A' : readfile(params['ancestorAgeno']), 'B' : readfile(params['ancestorBgeno'])}
    ind  = {'A' : readfile(params['ancestorAind']) , 'B' : readfile(params['ancestorBind']) }
    snp  = {'A' : readfile(params['ancestorAsnp']) , 'B' : readfile(params['ancestorBsnp']) }

    assert len(geno['A']) == len(snp['A']), "Pool A geno data length and number of snps don't line up!"
    assert len(geno['B']) == len(snp['B']), "Pool B geno data length and number of snps don't line up!"
    for item in geno['A']:
        assert len(item) == len(ind['A']), "Pool A geno number of individuals and ind file length don't line up!"
    for item in geno['B']:
        assert len(item) == len(ind['B']), "Pool B geno number of individuals and ind file length don't line up!"

    assert len(geno['A']) == len(geno['B']), "Pool A and B must have the same number of snps"

    haploid = params['haploidoutput']
    N = params['n'] if haploid else params['n'] * 2
    theta = params['theta']
    lamb = params['lambda']
    index, free = {}, {}
    index['A'] = list(range(len(ind['A'])))
    index['B'] = list(range(len(ind['B'])))
    active, free['A'], free['B'] = initialize_states(N, theta, index)

    data = [data_row(geno, active, 0, haploid)]
    ancestors = [ancestry(geno, active, 0, haploid)]
    interval = int(len(snp['A']) / 100)
    distance = lambda cur, prev : float(cur.split()[2]) - float(prev.split()[2])
    frac_done = 0
    for i in range(1, len(snp['A'])):
        d = distance(snp['A'][i], snp['A'][i-1])
        if snp['A'][i].split()[1] == snp['B'][i-1].split()[1]: 
            crossovers = [random.uniform(0, 1) > e**(-lamb * d) for _ in active]
            active, free = shuffle_states(crossovers, active, free, theta)
        else:
            index, free = {}, {}
            index['A'] = list(range(len(ind['A'])))
            index['B'] = list(range(len(ind['B'])))
            active, free['A'], free['B'] = initialize_states(N, theta, index)


        data.append(data_row(geno, active, i, haploid))
        if params["trackancestry"]:
            ancestors.append(ancestry(geno, active, i, haploid))
        if (i // interval != frac_done):
            frac_done = i // interval
            #print("Simulation is " + str(frac_done) + "% done", end='\r')
    
    print()


    output = params["outputfilename"]
    suffix = ".phgeno" if haploid else ".geno"
    print("Dumping results to " + output + suffix)
    with open(output + suffix, 'w') as f:
        for entry in data:
            f.write(entry + "\n")

    suffix = ".phind" if haploid else ".ind"
    print("Dumping results to " + output + suffix)
    with open(output + suffix, 'w') as f:
        indices = ["{0:05}".format(i) for i in range(params['n'])]
        i = 0
        k = 0
        while i < params['n']:
            phind = ":B" if i % 2 else ":A"
            if haploid:
                f.write("NA" + indices[k] + phind + " U        Simulation\n")
            else:
                f.write("NA" + indices[k] + " U        Simulation\n")

            k = k + i % 2 if haploid else k + 1
            i += 1

    suffix = ".phsnp" if haploid else ".snp"
    print("Dumping results to " + output + suffix)
    copyfile(params['ancestorAsnp'], output + suffix)

    if params["trackancestry"]:
        print("Dumping results to " + output + ".ancestry")
        with open(output + ".ancestry", 'w') as f:
            for entry in ancestors:
                f.write(entry + "\n")

    print("Job done!")


def initialize_states(N, theta, index):
    """
    Returns a length N list of dictionary elements in the form:
    [ { "A" : i_1 }, { "A" : i_2 }, {"B" : i_3 } ... ]
    Where each key is "A" with probability theta, and B otherwise,
    and each corresponding index is drawn from the corresponding index pool
    """
    active = []
    for _ in range(N):
        key = "A" if random.uniform(0, 1) < theta else "B"
        random.shuffle(index[key])
        val = index[key].pop(0)
        active.append({key : val})
    
    return active, set(index['A']), set(index['B'])

   
def data_row(geno, active, i, haploid=True):
    row = ''
    if haploid:
        for item in active:
            key, val = list(item.items())[0]
            row = row + geno[key][i][val]
    else:
        for j in range(0, len(active), 2):
            key1, val1 = list(active[j].items())[0]
            key2, val2 = list(active[j+1].items())[0]
            item1 = int(geno[key1][i][val1]) 
            item2 = int(geno[key2][i][val2]) 
            row = row + str(item1 + item2)
    return row

    
def ancestry(geno, active, i, haploid=True):
    anc_row = ''
    if haploid:
        for item in active:
            key, val = list(item.items())[0]
            anc_row = anc_row + '-' + key
    else:
        for j in range(0, len(active), 2):
            key1, val1 = list(active[j].items())[0]
            key2, val2 = list(active[j+1].items())[0]
            anc_row = anc_row + '-' + key1 + key2
    return anc_row


def shuffle_states(crossovers, active, free, theta):
    # print('A: ', len(free['A']), 'B:', len(free['B']))
    # trues = [item for item in crossovers if item ]
    # print(len(free['A']) + len(free['B']), len(trues))
    for i in range(len(crossovers)):
        if not crossovers[i]:
            continue
        key, val = list(active[i].items())[0]
        free[key].add(val)
    # print('A: ', len(free['A']), 'B:', len(free['B']))
    # print(len(free['A']) + len(free['B']))
    newActive = []
    for i in range(len(crossovers)):
        if not crossovers[i]:
            newActive.append(active[i])
            continue
        key, curVal = list(active[i].items())[0]
        # print(key, curVal)
        pool = "A" if random.uniform(0, 1) < theta else "B"
        newVal = random.sample(free[pool] - {curVal}, 1)[0]
        # print(pool, newVal)
        newActive.append({pool : newVal})
        free[pool].remove(newVal)
    return newActive, free



if __name__=="__main__":
    main()
