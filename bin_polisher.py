#!/usr/bin/python
import argparse
import numpy #using numpy for mean and stdev calculation, as this script requires biopython anyway, and numpy is already a dependancy of biopython
import itertools
import sys
from Bio import SeqIO


version="v0.01 (July 2016)"
parser=argparse.ArgumentParser(description = "bin_polisher {} : Purify pre-generated bins (e.g. using Maxbin) based on z-score differences in contig coverage using mutliple sequencing datasets".format(version))
input_args = parser.add_argument_group("input")
input_args.add_argument("-if", "--input_fasta", action = "store", nargs = "+", dest = "input_fasta_list", required = True, help = "Input fasta file(s) of a (single) bin")
input_args.add_argument("-ic", "--input_coverage", action = "store", nargs = "+", dest = "input_coverage_list", required = True, help = "Abundance files (in maxbin \".abund\"-format)")
#TODO: add option to generate abundance info in this script using SAM/BAM-files (or even from fastqs). Remove the above "required"-option and change above nargs to "*" in that case
parser.add_argument("-pr", "--pre-remove", action = "store", choices = ["high", "low", "both", "none"], dest = "pre_remove", default = "none", help = "remove extreme values based on upper (99%%) and lower (1%%) percentiles, or both. default = none")
parser.add_argument("-mi", "--max_iterations", action = "store", type = int, dest = "max_iterations", default = 10, help = "maximum number of iterations for recalculating and comparing coverage average, standard deviation and z-score values and removing difference-outliers. Default = 10")
parser.add_argument("--out_bad", action = "store_true", dest = "out_bad", default = False, help = "create seperate output files for rejected reads also. default = False")
parser.add_argument("-op", "--out_prefix", action = "store", dest = "out_prefix", default = "bin_polisher", help = "prefix for output file(s). default = \"bin_polisher\"")
parser.add_argument("-v", "--version", action = "store_true", dest = "show_version", help = "show version info and quit")
args = parser.parse_args()


def get_compression_type(infasta):
	if fastafile.endswith(".gz"):
		return "gz"
	elif fastafile.endswith(".bz2"):
		return "bzip2"
	elif fastafile.endswith(".zip"):
		return "zip"
	else:
		return "none"

def read_fasta(infasta):
	compression_type = get_compression_type(infasta)
	try:
		if compression_type in ["gzip", "gz"]:
			readfile = gzip.open(fastafile, 'r')
		elif compression_type in ["bzip2", "bz2"]:
			readfile = bz2.BZ2File(fastafile, 'r')
		elif compression_type == "zip":
			myzipfile = zipfile.ZipFile(fastafile, 'r')
			if len(myzipfile.namelist()) > 1:
				raise IOError("TOO MANY FILES IN ZIPFILE")
			else:
				readfile = myzipfile.open(myzipfile.namelist()[0])
		else:
			readfile = open(fastafile, 'r')
	except Exception, ex:
		sys.stderr.write(ex.__class__.__name__ + " : " + str(ex))
		return None
	else:
		in_fasta_iterator = SeqIO.parse(readfile, "fasta")
		bindict = { record.id : {"record": record} for record in in_fasta_iterator}#coverage values will be added for each record in the below "read_coverage" function
		return bin_dict

class bin_object(object):
	def __init__(self, input_fasta_list, input_coverage_list):
		self.bindict = {}
		for fasta in input_fasta_list:
			self.bindict.update(read_fasta(fasta))
		self.cov_datasets = []
		for incov in input_coverage_list:
			self.cov_datasets.append(coverage_dataset(incov, bindict))
		
	def pre_remove_upper_extremes(self):
		del_list = []
		for c_d in self.cov_datasets:
			del_list.extend(c_d.pre_remove_upper_extremes())
		del_list = list(set(del_list)) #remove duplicates from list
		
		for c_d in self.cov_datasets: #bring all datasets to the same level
			c_d.remove_and_update_covstats(del_list)
				
		del_records = []
		for key in del_list:
			del_records.append(self.bindict.pop(key))
		return del_records #return deleted records for optional exporting to fasta-file (to double check WHAT is being removed)
		
	def pre_remove_lower_extremes(self):
		del_list = []
		for c_d in self.cov_datasets:
			del_list.extend(c_d.pre_remove_lower_extremes())
		
		del_list = list(set(del_list)) #remove duplicates from list
		
		for c_d in self.cov_datasets: #bring all datasets to the same level
			c_d.remove_and_update_covstats(del_list)
		
		del_records = []
		for key in del_list:
			del_records.append(self.bindict.pop(key))
		return del_records #return deleted records for optional exporting to fasta-file (to double check WHAT is being removed)
		
	def pre_remove_both_extremes(self):
		del_records = []
		del_records.extend(self.pre_remove_upper_extremes())
		del_records.extend(self.pre_remove_lower_extremes())
		return del_records #return deleted records for optional exporting to fasta-file (to double check WHAT is being removed)
		
	def filter_zscore_differences(self):
		dataset_combinations = itertools.combinations(self.cov_datasets, 2)
		del_list = []
		for combination in dataset_combinations:
			del_list.extend(combination[0].compare_zscores(combination[1]))
		
		del_list = list(set(del_list)) #remove duplicates from list
		
		for c_d in self.cov_datasets: #bring all datasets to the same level
			c_d.remove_and_update_covstats(del_list)
		
		del_records = []
		for key in del_list:
			del_records.append(self.bindict.pop(key))
		return del_records #return deleted records for optional exporting to fasta-file (to double check WHAT is being removed)
		
class coverage_dataset(object):
	def __init__(self, incov, bindict):
		self.cov_value = self.read_coverage(incov, bindict)
		self.average, self.stdev, self.zscore, self.p99, self.p1 = self.calc_zscore(self.cov_value)
		
	def read_coverage(self, incov, bindict):
		#bla
		incov_file = open(incov, "r")
		covdict = {line.split()[0] : line.split()[1] for line in incov_file if line.split()[0] in bindict} #add only coverage info of contigs that are actually in the bin
		incov_file.close()
		return covdict
		
	def calc_zscore(self, covdict):
		cov_values = [covdict[key] for key in covdict]
		average = numpy.mean(cov_values) 
		stdev = numpy.std(cov_values)
		zscoredict = {key : ((covdict[key] - average)/stdev) for key in covdict}
		p99 = numpy.percentile(cov_values, 99) # upper extremes
		p1 = numpy.percentile(cov_values, 1) # lower extremes
		return average, stdev, zscoredict, p99, p1
		
	def remove_and_update_covstats(self, del_list):
		for key in del_list:
			if key in self.cov_value:
				self.cov_value.pop(key)
		self.average, self.stdev, self.zscore, self.p99, self.p1 = self.calc_zscore(self.cov_value)
		
	def pre_remove_upper_extremes(self):
		del_list = []
		for key in cov_value:
			if cov_values[key] > p99:
				del_list.append(key)
		self.remove_and_update_covstats(del_list)
		return del_list
		
	def pre_remove_lower_extremes(self):
		del_list = []
		for key in cov_value:
			if cov_values[key] < p1:
				del_list.append(key)
		self.remove_and_update_covstats(del_list)
		return del_list
		
	def pre_remove_both_extremes(self):
		del_list = []
		del_list.extend(self.pre_remove_upper_extremes())
		del_list.extend(self.pre_remove_lower_extremes())
		return del_list
		
	def compare_zscores(self, other_coverage_dataset):
		del_list = []
		for key in self.zscore:
			this_zscore = self.zscore[key]
			other_zscore = other_coverage_dataset.zscore[key]
			if abs(this_zscore - other_zscore) > 4:
				del_list.append(key)
		#self.remove_and_update_covstats(del_list) #the update is best done from the calling function AFTER all combinations have been compared. otherwise the order of the input influences the outcome
		#other_coverage_dataset.remove_and_update_covstats(del_list) #see comment above
		return del_list

def pre_remove_extremes(mybin):
	assert args.pre_remove in ["high", "low", "both"], "Error: Illegal calling of function \"preremove_extremes\" with \"-pr {}".format(args.pre_remove)
	preremove_record_list = []
	if args.pre_remove == "high":
		preremove_record_list = mybin.pre_remove_upper_extremes()
	if args.pre_remove == "low":
		preremove_record_list = mybin.pre_remove_lower_extremes()
	if args.pre_remove == "both":
		preremove_record_list = mybin.pre_remove_both_extremes()
	if args.out_bad:
		SeqIO.write(preremove_record_list, "{}_preremoved.fasta".format(args.out_prefix), "fasta")
	return len(preremove_record_list)

def main():
	my_bin = bin_object(args.input_fasta_list, args.input_coverage_list)
	counter = 0
	if args.pre_remove != "none":
		counter += pre_remove_extremes(mybin)
		sys.stderr.write("pre-remove extremes set to \"{}\" : removed {} records\n".format(args.pre_remove, counter))
	
	for i in range(1, args.max_iterations + 1):
		remove_record_list, out_good_list = [], []
		outname_bad = "{}_REMOVED_FilterIteration{}.fasta".format(args.out_prefix, i)
		outname_good = "{}_KEPT_FilterIteration{}.fasta".format(args.out_prefix, i)
		remove_record_list = my_bin.filter_zscore_differences()
		if args.out_bad:
			SeqIO.write(remove_record_list, outname, "fasta")
		sys.stderr.write("Iteration {} : removed {} records --> {} removed in total\n".format(i, len(remove_record_list), counter))
		out_good_list = my_bin.bindict.values()
		sys.stderr.write("\twriting {} remaining filtered contigs to {}\n".format(len(out_good_list), outname_good))
	sys.stderr.write("===FINISHED!===")

main()
