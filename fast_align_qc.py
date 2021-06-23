#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Ziwei Xue
# Modified by YU Zhejian

import argparse
import gzip

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import pysam
from matplotlib import gridspec

BLACK = "\033[0;30m"
DGRAY = "\033[1;30m"
RED = "\033[0;31m"
LRED = "\033[1;31m"
GREEN = "\033[0;32m"
LGREEN = "\033[1;32m"
ORANGE = "\033[0;33m"
YELLOW = "\033[1;33m"
BLUE = "\033[0;34m"
LBLUE = "\033[1;34m"
PURPLE = "\033[0;35m"
LPURPLE = "\033[1;35m"
CYAN = "\033[0;36m"
LCYAN = "\033[1;36m"
LGRAY = "\033[0;37m"
WHITE = "\033[1;37m"
BRED = "\033[0;37;41m"
BGREEN = "\033[0;37;42m"
BYELLOW = "\033[0;37;43m"
BBLUE = "\033[0;37;44m"
NC = "\033[0m"

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fastq', type=str, help='input read file in (GZipped) fastq format')
parser.add_argument('-b', '--bam', type=str, help='input alignment file in sam/bam. Index needed for bam')
parser.add_argument('-o', '--out', type=str, help='output prefix')

args = parser.parse_args()

all_reads = dict()  # Store the fasta file in readId-read pair
perc_alignment = 0

print(LGREEN + "Reading\t" + NC + args.fasta, end="")

flag = False
if args.fasta.endswith(".gz"):
	read_id = ""
	with gzip.open(args.fasta, 'r') as fp:
		while s := fp.readline():
			s = s.decode("utf-8")
			if flag:
				all_reads[read_id] = s.strip()
				flag = False  # Read one line only
			if s.startswith("@"):
				# @84443138-6ebe-4314-9eb0-372e5d1220c6 runid=5c36a9dc99329515d5cdea36d6fbb6fd327cd511 read=19358 ch=1108 start_time=2020-10-31T15:01:48Z flow_cell_id=PAF26057 protocol_group_id=FONR20H101138-1A sample_id=FONR20H101138-1A barcode=barcode23
				read_id = s.split(" ")[0][1:]
				flag = True

elif args.fasta.endswith(".fastq") or args.fasta.endswith(".fq") or args.fasta.endswith(
		".fasta") or args.fasta.endswith(".fa"):
	with open(args.fasta, 'r') as fp:
		while s := fp.readline():
			if flag:
				all_reads[read_id] = s.strip()
			elif s.startswith("@"):
				read_id = s.split(" ")[0][1:]
			else:
				pass
else:
	raise TypeError("The read sequence should be either suffixed by fastq, fq, fasta or fa.")

if len(all_reads) == 0:
	raise TypeError("No readable reads identified. Do you passed in a reference genome?")
print(str(len(all_reads)) + " was read")

print(LGREEN + "Reading\t" + NC + args.bam)

if (args.bam.endswith("bam")):
	mapped_reads = pysam.AlignmentFile(args.bam, 'rb')
elif (args.bam.endswith("sam")):
	mapped_reads = pysam.AlignmentFile(args.bam, 'r')
else:
	raise TypeError("The alignment file should be either suffixed by bam or sam.")
all_reads_tmp = set(all_reads.keys())  # All readIds already identified in fastq
mapped_reads_tmp = set()  # readIds already identified in BAM
multi_mapped = set()  # All multi-mapped readIds
alignment_stats = []
perc_bases = 0
mapped_length = []  # Mapped reads length distribution
unmapped_length = []  # Unmapped reads length distribution
mapped_num = 0

# Step through all of the reads
# FIXME: Need to be implemented
for i in mapped_reads.fetch():  # Get all reads from BAM
	if i.qname in all_reads_tmp and not i.is_unmapped:
		if i.qname in mapped_reads_tmp:  # Reads already identified, which means it is multi-mapped (secondary)
			multi_mapped.add(i.qname)
			mapped_length.append(["Multiply Aligned", len(all_reads[i.qname])])
		else:  # Reads not identified
			all_reads_tmp.remove(i.qname)
			mapped_reads_tmp.add(i.qname)
			print(i.qname)
			cigar = i.cigartuples
			aligned_bases = sum(map(lambda x: x[1], filter(lambda y: y[0] == 0, cigar)))  # ?
			perc_bases = (perc_bases * (mapped_num - 1) + aligned_bases / i.query_length) / mapped_num  # ?
		mapped_num += 1

for i in mapped_reads_tmp:
	unmapped_length.append(["Uniquely Aligned", len(all_reads[i])])

for i in all_reads_tmp:
	unmapped_length.append(["Not Aligned", len(all_reads[i])])

perc_alignment = mapped_num / len(all_reads)

print(YELLOW + str(round(perc_alignment * 100, 2))+ "%" + NC + "\tof reads aligned")

df = pd.DataFrame(mapped_length + unmapped_length)
df.columns = ["Align", "Length"]
df_agg = df.groupby("Align")
vals = {i: df["Length"].values.tolist() for i, df in df_agg}

print(YELLOW + str(len(vals["Uniquely Aligned"]) / mapped_num * 100) + "%" + NC + "\tof uniquely aligned reads")
if "Multiply Aligned" in vals.keys():
	print(YELLOW + str(len(vals["Multiply Aligned"]) / mapped_num * 100) + "%" + NC + "\tof multiply aligned reads")
else:
	print(YELLOW + "0%" + NC + "\tof multiply aligned reads")

if len(vals.keys()) == 2:
	colors = ["#777777", "#FFFFFF"]
	vals = [vals["Uniquely Aligned"], vals["Not Aligned"]]
else:
	colors = ["#777777", "#057DFF", "#FFFFFF"]
	vals = [vals["Uniquely Aligned"], vals["Multiply Aligned"], vals["Not Aligned"]]

matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['font.size'] = '16'
matplotlib.rcParams['font.weight'] = 100
matplotlib.rcParams['axes.linewidth'] = 2.5
matplotlib.rcParams['axes.edgecolor'] = "#000000"


def createFig(ax):
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	ax.spines['bottom'].set_color('none')
	# ax.spines['bottom'].set_color('none')
	# ax.spines['left'].set_color('none')
	for line in ax.yaxis.get_ticklines():
		line.set_markersize(5)
		line.set_color("#585958")
		line.set_markeredgewidth(2.5)
	for line in ax.xaxis.get_ticklines():
		line.set_markersize(5)
		line.set_markeredgewidth(2.5)
		line.set_color("#585958")
	ax.set_xbound(0, 10)
	ax.set_ybound(0, 10)


gs = gridspec.GridSpec(1, 3, width_ratios=[1, 6, 1])
fig = plt.Figure(figsize=(12, 8))
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
ax3 = fig.add_subplot(gs[2])
createFig(ax1)
createFig(ax2)
createFig(ax3)
ax1.hist(vals, bins=1, range=[0, max(max(vals[0]), max(vals[0]))], stacked=True, density=False, color=colors,
		 edgecolor='black', linewidth=2.5)
ax2.hist(vals, bins=[0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000], range=[-500, 4000], stacked=True, density=False,
		 color=colors, edgecolor='black', linewidth=2.5)
ax3.hist([list(filter(lambda x: x > 4000, vals[0])), list(filter(lambda x: x > 4000, vals[1]))], bins=1,
		 range=[0, max(max(vals[0]), max(vals[0]))], stacked=True, density=False, color=colors, edgecolor='black',
		 linewidth=2.5)
# ax1.ticklabel_format(useOffset=False, style='plain')
# ax2.ticklabel_format(useOffset=False, style='plain')
# ax3.ticklabel_format(useOffset=False, style='plain')
ax1.tick_params(axis='y', labelrotation=90)
ax2.tick_params(axis='y', labelrotation=90)
ax3.tick_params(axis='y', labelrotation=90)
ax1.set_ylabel("Number of Reads")
ax1.set_xticks([])
ax3.set_xticks([])
ax1.set_xlabel("All")
ax3.set_xlabel(">4000")

if args.out:
	fig.savefig(args.out + "pdf")
else:
	fig.savefig("alignqc.pdf")

print(LGREEN + "Finished \t" + NC + "AlignQC")
