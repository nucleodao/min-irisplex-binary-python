#!/usr/bin/env python3.10
'''
irisplex.py: irisplex simply implementation
'''
import math
import os
import argparse
from numpy import append
import requests

VERSION = "0.1.0"
PROJECT_NAME = "irisplex"
CURDIR = os.getcwd()

eye_model_data = """
Constant	X	2.575575265	0.306416526
rs12913832	T	-5.385897557	-2.302942706
rs1800407	A	1.330922472	0.9785432025
rs12896399	T	0.7636299425	0.2541024725
rs16891982	C	-1.529280117	-0.9342328224
rs1393350	T	0.4340394174	0.2087410752
rs12203592	T	0.653488393	0.6457702022
"""

eye_colors = ["Blue", "Intermediate", "Brown"]

class EyeModel:
  def __init__(self, rsid, allele, weights):
    self.rsid = rsid
    self.allele = allele
    self.weights = weights

class Genome:
  def __init__(self, rsid, chromosome, position, allele):
    self.rsid = rsid
    self.chromosome = chromosome
    self.position = position
    self.allele = allele

def parse_args():
    parser = argparse.ArgumentParser(prog='irisplex', description='irisplex v{}'.format(VERSION))
    parser.add_argument("-V", '--version', action='version', version='%(prog)s {}'.format(VERSION))
    parser.add_argument("-u", '--url', action='store', help="url to parse")
    parser.add_argument("-p", '--path', action='store', help="file path to parse")
    parser.add_argument("-c", '--content', action='store', help="raw genome file string content")
    
    return parser.parse_args()

def parse_eye_model():
  eye_models_dic = {}
  for line in eye_model_data.splitlines():
    segments = line.split('\t', -1)
    if len(segments) < 3: continue
    rsid = segments[0]
    allele = segments[1]
    weights = list(map(lambda s:float(s),segments[2:]))
    eye_models_dic[rsid] = EyeModel(rsid, allele, weights)

  return eye_models_dic

def parse_genome(url, path, content):
  genome_dic = {}
  lines = []
  if url is not None:
    response = requests.get(url)
    path = 'genome.txt'
    open(path, 'wb').write(response.content)
  if path is not None:
    file = open(path, 'r')
    content = file.read()
  if content is not None:
    lines = list(filter(lambda l:l.startswith('#') == False, content.splitlines()))

  for line in lines:
    segments = line.split('\t', -1)
    if len(segments) != 4: continue
    rsid = segments[0]
    chromosome = segments[1]
    position = segments[2]
    allele = segments[3]
    genome_dic[rsid] = Genome(rsid, chromosome, position, allele)
  
  return genome_dic

def analyze():
    args = parse_args()
    eye_model_dic = parse_eye_model()
    genome_dic = parse_genome(args.url, args.path, args.content)

    sum_weights = eye_model_dic['Constant'].weights
    for rsid, genome in genome_dic.items():
      eye_model = eye_model_dic.get(rsid)
      
      if eye_model is None: continue

      allele = genome.allele
      allele_number = 0.0
      for _allele in allele:
        if _allele == eye_model.allele: allele_number += 1.0

      # calculate the adjusted weights and add to the weights vector
      count = 0
      for weight in eye_model.weights:
        sum_weights[count] += weight * allele_number
        count += 1

      # calculate the exponents
      exp_weights = []
      total_exp = 1.0
      for weight in sum_weights:
        exp = math.exp(weight)
        exp_weights.append(exp)
        total_exp += exp
      
      # calculate the probs
      probs = []
      total_prob = 0.0
      for exp in exp_weights:
        prob = exp / total_exp
        probs.append(prob)
        total_prob += prob
      probs.append(1.0 - total_prob)

      # display answers
      answers = {}
      count = 0
      for eye_color in eye_colors:
        answers[eye_color] = probs[count]
        count += 1
      
      return answers

if __name__ == '__main__':
    print(analyze())