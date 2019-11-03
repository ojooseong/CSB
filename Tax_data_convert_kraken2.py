########################################################################################################################
##
##  This tool convert taxonomic profiling data to table form for visualizing using R packages.
##  Contact : Jooseong Oh (ojooseong@gmail.com)
##
########################################################################################################################

import sys
import os
import argparse as ap
import operator

## Clade name return
def get_cladeName(clade_char):
	if "k" == clade_char.lower(): return "Kingdom"
	elif "p" == clade_char.lower(): return "Phylum"
	elif "c" == clade_char.lower(): return "Class"
	elif "o" == clade_char.lower(): return "Order"
	elif "f" == clade_char.lower(): return "Family"
	elif "g" == clade_char.lower(): return "Genus"
	elif "s" == clade_char.lower(): return "Species"
	else: return "Uknown"

## list make
def list_initiation(list_size):
	return [0] * list_size

## Relative abundance merge.
## Clade that average of relative abundance smaller than relative_abundance_cutoff is integrated to Others.
def integrate_relative_abundance(outputTable, cladeName, relativeAbundance):
	if cladeName in outputTable.keys():
		outputTable[cladeName] = list(map(operator.add, outputTable[cladeName], relativeAbundance))
	else:
		outputTable[cladeName] = relativeAbundance

	return outputTable

## Average of relative abundance calculate.
def average(relAbun):
	return sum(relAbun) / len(relAbun)

## Input parameter parse.
def parse_parameter(args):
	p = ap.ArgumentParser(description=" DESCRIPTION ", formatter_class=ap.RawTextHelpFormatter, add_help=False)
	arg = p.add_argument
	arg('-i', '--input', dest='inputfile', type=ap.FileType('r'), nargs="+", required=True)
	arg('-o', '--output', dest='outputfile', type=str, default=sys.stdout)
	arg('-c', '--relative_abundance_cutoff', dest='abundanceCutoff', type=float, nargs='?', default=3.0)

	return p.parse_args()

def main(args):
	param = parse_parameter(args)

	inputfile_list = param.inputfile
	outputfile = param.outputfile
	abundCutoff = param.abundanceCutoff

	numberOfSamples = len(inputfile_list)
	taxonomyLevel = set()
	inputTable = dict()
	header = ["Clade_name"]
	numberOfSpecies_inOthers = 0
	for idx in range(0, len(inputfile_list)):
		inputfile = inputfile_list[idx]
		header.append(inputfile.name)
		inputfile.readline()

		for line in inputfile:
			columns = line.strip().split("\t")
			cladeName = columns[0]
			taxonomyLevel.add(columns[2])
			relativeAbundance = float(columns[6]) * 100

			if cladeName in inputTable.keys():
				inputTable[cladeName][idx] = relativeAbundance
			else:
				inputTable[cladeName] = list_initiation(numberOfSamples)
				inputTable[cladeName][idx] = relativeAbundance
	
	if len(taxonomyLevel) != 1:
		sys.stderr("Input data is not consistancy for taxonomy level.")
	else:
		taxonomyLevel = get_cladeName(taxonomyLevel.pop())

	outputTable = dict()
	for key in inputTable.keys():
		average_relativeAbundance = average(inputTable[key])
		if average_relativeAbundance <= float(abundCutoff):
			outputTable = integrate_relative_abundance(outputTable, "Others", inputTable[key])
			numberOfSpecies_inOthers += 1
		else:
			outputTable[key] = inputTable[key]

	with open(outputfile, "w") as of:
		of.write("\t".join(header) + "\n")
		for key in sorted(outputTable.keys()):
			if "Others" != key: of.write(key + "\t" + "\t".join(list(map(str, outputTable[key][:]))) + "\n")
			else: of.write(key + "(" + str(numberOfSpecies_inOthers) + ")" + "\t" + "\t".join(list(map(str, outputTable[key][:]))) + "\n")

if __name__ == '__main__':
	main(sys.argv[1:])
