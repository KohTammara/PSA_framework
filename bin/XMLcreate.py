#!/usr/bin/env python
import sys
import argparse
import lxml.etree as ET

parser = argparse.ArgumentParser(description='Creates an XML file making use of a template and parameters provided by user for RosettaSimpleThreading mover')

parser.add_argument("-name", type=str, help='Name you would like to provide to mover', required=True)
parser.add_argument("-pack_neighbors", type=str, help='', required=True)
parser.add_argument("-neighbor_dis", type=int, help='')
parser.add_argument("-start_position", help='')##is this required
parser.add_argument("-thread_sequence", help='')
parser.add_argument("-scorefxn", help='')
parser.add_argument("-skip_unknown_mutant", help='')
parser.add_argument("-pack_rounds", help='')
parser.add_argument("-sequence_mode", help='')

args = parser.parse_args()

#Check if all variables are populated and assemble the input as necessary
mover = "<SimpleThreadingMover "
if args.name != "":
  name = 'name="'+ args.name +'" ' 

if args.pack_neighbors != "":
  pack_neighbors = 'name="'+ args.pack_neighbors +'" ' 

if args.neighbor_dis != "":
  neighbor_dis = 'name="'+ args.neighbor_dis +'" ' 

if args.scorefxn != "":
  scorefxn = 'name="'+args.scorefxn +'" ' 

if args.skip_unknown_mutant != "":
  skip_unknown_mutant = 'name="'+args.skip_unknown_mutant +'" ' 

if args.pack_rounds != "":
  pack_rounds = 'name="'+args.pack_rounds +'" ' 

if args.sequence_mode != "":
  sequence_mode = 'name="'+args.sequence_mode +'" ' 

mover = mover+name + pack_neighbors + neighbor_dis + scorefxn + skip_unknown_mutant + pack_rounds + sequence_mode + '/>'

#using lxml instead of xml preserved the comments
#adding the encoding when the file is opened and written is needed to avoid a charmap error
with open('template.xml', encoding="utf8") as f:
  tree = ET.parse(f)
  root = tree.getroot()
  for elem in root.getiterator():
    try:
      elem.text = elem.text.replace('input', mover)
    #   elem.text = elem.text.replace('FEATURE NUMBER', '123456')
    except AttributeError:
      pass
#tree.write('output.xml', encoding="utf8")
# Adding the xml_declaration and method helped keep the header info at the top of the file.
tree.write('output.xml', xml_declaration=True, method='xml', encoding="utf8")