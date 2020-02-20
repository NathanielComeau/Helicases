#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

def asciiToScore(score):
	temp = list(score)
	return [ord(x) - 33 for x in score]

def writeAverageQualityScores(qaFile, outputFile):
	# Writes a list of quality scores to a text file
	with open(qaFile, 'r') as f1:
		with open(outputFile, 'w') as f2:
			line = f1.readline()
			while line:
				line = line.rstrip('\n')
				line = asciiToScore(line)
				line = sum(line)/len(line)
				f2.write(str(line)+'\n')
				line = f1.readline()

writeAverageQualityScores(sys.argv[1], sys.argv[2])
