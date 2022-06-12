# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 13:29:34 2021

@author: jaspe
"""

import re

def tree2Tuple(treeStr):
    if treeStr[0] != '(':
        leaf = treeStr.split(":")
        leaf[1] = float(leaf[1])
        return tuple(leaf)
    else:
        bracketStatus = 0
        split_idx = None
        end_idx = None
        for i, char in enumerate(treeStr):
            if char == '(':
                bracketStatus += 1
            elif char == ')':
                bracketStatus -= 1
            if bracketStatus == 1 and char == ',':
                split_idx = i
            if bracketStatus == 0:
                end_idx = i
                break
        nodeLength = treeStr[end_idx+1:]
        if nodeLength[0] == ';':
            nodeLength = 0
        else:
            nodeLength = float(nodeLength.split(":")[1])
        return ((tree2Tuple(treeStr[1:split_idx]),tree2Tuple(treeStr[split_idx+1:end_idx])),nodeLength)    

def tuple2Tree(tupleTree, fullTree=True):
    if type(tupleTree[0]) is str: 
        return tupleTree[0]+":"+str(tupleTree[1])
    else:
        if fullTree:
            return "("+tuple2Tree(tupleTree[0][0],fullTree=False)+","+tuple2Tree(tupleTree[0][1],fullTree=False)+");"
        else:
            return "("+tuple2Tree(tupleTree[0][0],fullTree=False)+","+tuple2Tree(tupleTree[0][1],fullTree=False)+"):"+str(tupleTree[1])

def getBLS(tree, species, finished=True):
    BLS = 0
    foundSp = set()
    if type(tree[0]) is str:
        if tree[0] in species:
            foundSp.add(tree[0])
            BLS = tree[1]
    else:
        BLS1, foundSp1 = getBLS(tree[0][0], species, finished=False)
        BLS2, foundSp2 = getBLS(tree[0][1], species, finished=False)
        BLS, foundSp = BLS1+BLS2, foundSp1.union(foundSp2)
        nbOfFoundSp = len(species.intersection(foundSp))
        if nbOfFoundSp != 0 and nbOfFoundSp != len(species):
            BLS += tree[1]
    if finished:
        return BLS
    else:
        return BLS, foundSp

def getTotalScore(tree):
    score = 0
    if type(tree[0]) is str:
        score += tree[1]
    else:
        score += (getTotalScore(tree[0][0]) + getTotalScore(tree[0][1]) + tree[1])
    return score

def rescaleTupleTree(tree, totalScore=None):
    if totalScore is None:
        totalScore = getTotalScore(tree)
    rescaledScore = tree[1]/totalScore
    if type(tree[0]) is str:
        subtree = tree[0]
    else:
        subtree = (rescaleTupleTree(tree[0][0],totalScore=totalScore),rescaleTupleTree(tree[0][1],totalScore=totalScore)) 
    return (subtree,rescaledScore)

def rescaleTree(treeStr):
    scaledTupleTreeStr = str(rescaleTupleTree(tree2Tuple(treeStr)))
    scaledScores = []
    for scoreMatch in re.finditer(', [0-9.]*e?-?[0-9]*\)',scaledTupleTreeStr):
        scaledScores.append(scoreMatch.group()[2:-1])
    scaledTreeStr = []
    scoreChars = False # Whether we are looking at chars representing a score
    ScoreCounter = 0
    for char in treeStr:
        if not scoreChars:
            scaledTreeStr.append(char)
        elif char == ',' or char == ')': # while scoreChars is True
            scaledTreeStr.append(char)
            scoreChars = False
        if char == ':':
            scaledTreeStr.append(scaledScores[ScoreCounter])
            scoreChars = True
            ScoreCounter += 1
    return "".join(scaledTreeStr)

def cleanTree(treeStr, rmChildLabels=True):
    if rmChildLabels:
        return tuple2Tree(rescaleTupleTree(tree2Tuple(treeStr)))
    else:
        return rescaleTree
