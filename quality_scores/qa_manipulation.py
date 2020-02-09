#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

def asciiToScore(score):
	temp = list(score)
	return [ord(x) - 33 for x in score]

def writeQualityScores(qaFile, outputFile):
	# Writes a list of quality scores to a text file
	with open(qaFile, 'r') as f:
		scores = f.readlines()
	scores = [x.rstrip('\n') for x in scores]
	# Convert scores to ints
	int_scores = [asciiToScore(x) for x in scores]
	#print(int_scores)

	with open(outputFile, 'w') as f:
		for scores in int_scores:
			str_scores = [str(x) for x in scores]
			f.write(",".join(str_scores)+'\n')

def plotQualityScores(qaFile, i, j):
	with open(qaFile, 'r') as f:
		scores = f.readlines()
	scores = [x.rstrip('\n') for x in scores]
	# Convert scores to ints
	int_scores = [asciiToScore(x) for x in scores]

	#print(scores[i])
	#print(int_scores[i])

	plt.rcdefaults()
	fig, ax = plt.subplots()
	for k in range(i,j):
		x = np.arange(100) # Hard-coded read length of 100
		y = np.array(int_scores[k])
		ax.plot(x, y)
	plt.show()


def plotFromConverted(qaFile, i, j):
	with open(qaFile, 'r') as f:
		scores = f.readlines()
	scores = [x.rstrip('\n').split(',') for x in scores]
	#print(scores)
	plt.rcdefaults()
	fig, ax = plt.subplots()
	for k in range(i,j):
		x = np.arange(100) # Hard-coded read length of 100
		int_scores = [int(x) for x in scores[k]]
		y = np.array(int_scores)
		ax.plot(x, y)
	plt.show()


# Plot quality score 0 to 1000 from "quality_scores.txt"
#plotQualityScores("quality_scores.txt", 0, 100)
#writeQualityScores("quality_scores.txt", 'converted_quality_scores.txt')
plotFromConverted("converted_quality_scores.txt", 0, 1)
plotFromConverted("converted_quality_scores.txt", 0, 10)
plotFromConverted("converted_quality_scores.txt", 0, 100)
