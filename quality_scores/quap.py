#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt

def asciiToScore(score):
	temp = list(score)
	return [ord(x) - 33 for x in score]


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


# Plot quality score 0 to 1000 from "quality_scores.txt"
plotQualityScores("quality_scores.txt", 0, 10)


