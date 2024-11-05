#!/usr/bin/env python
import sys
import argparse
import lxml.etree as ET

parser = argparse.ArgumentParser(description='Creates an XML file making use of a template and parameters provided by user for RosettaSimpleThreading mover')

parser.add_argument("-name", type=str, default="", help='Name you would like to provide to mover')
parser.add_argument("-sequence", type=str, required=True, help='Path to Sequence to be threaded')
parser.add_argument("-start_pos", type=str, required=True, help='Startong positon of sequence to be threaded')
parser.add_argument("-pack_neighbors", type=str, default="", help='Boolean option to pack neighbors while threading. By default false.')
parser.add_argument("-neighbor_dis", type=str, default="", help='Distance to repack neighbor side chains. Default is 6.0 Angstroms.')
parser.add_argument("-scorefxn", type=str, default="", help='Optional score function name passed')
parser.add_argument("-skip_unknown_mutant", default="", type=str, help='Boolean option to skip unknown amin acids in the sequence string instead of throwing exception')
parser.add_argument("-pack_rounds", type=str, default="", help='Number of packing rounds for threading.')
parser.add_argument("-sequence_mode", type=str, default="", help='Format for input sequence. (oneletter/threeletter/basename/fullname)')
parser.add_argument("-template", type=str, default="", help='Path to template xml file')
parser.add_argument("-weight", type=str, required=True, help='The desired weights for scoring')

args = parser.parse_args()


#Check if all variables are populated and assemble the input as necessary
mover = '<SimpleThreadingMover '
if args.name != "NONE":
    name = 'name="'+ args.name +'" ' 
    mover = mover+name
if args.sequence != "NONE":
    sequence = 'thread_sequence="'+ args.sequence +'" ' 
    mover = mover+sequence
if args.start_pos != "NONE":
    name = 'start_position="'+ args.name +'" ' 
    mover = mover+name
if args.pack_neighbors != "NONE":
    pack_neighbors = 'pack_neighbors="'+ args.pack_neighbors +'" ' 
    mover = mover + pack_neighbors
if args.neighbor_dis != "NONE":
    neighbor_dis = 'neighbor_dis="'+ args.neighbor_dis +'" ' 
    mover = mover + neighbor_dis
if args.scorefxn != "NONE":
    scorefxn = 'scorefxn="'+args.scorefxn +'" ' 
    mover = mover + scorefxn
if args.skip_unknown_mutant != "NONE":
    skip_unknown_mutant = 'skip_unknown_mutant="'+args.skip_unknown_mutant +'" ' 
    mover = mover + skip_unknown_mutant
if args.pack_rounds != "NONE":
    pack_rounds = 'pack_rounds="'+args.pack_rounds +'" ' 
    mover = mover + pack_rounds
if args.sequence_mode != "NONE":
    sequence_mode = 'sequence_mode="'+args.sequence_mode +'" ' 
    mover = mover + sequence_mode
mover = mover + '/>'

score = '<ScoreFunction '
if args.weight != "NONE":
    weight = 'weights="'+args.weight +'" ' 
    score_name = 'name="'+args.scorefxn +'" ' 
    score = score + score_name + weight + '/>'

#using lxml instead of xml preserved the comments
#adding the encoding when the file is opened and written is needed to avoid a charmap error
with open(args.template, encoding="utf8") as f:
  tree = ET.parse(f)
  root = tree.getroot()
  for elem in root.getiterator():
    try:
      elem.text = elem.text.replace('threader', mover)
      elem.text = elem.text.replace('score_input', score)
    except AttributeError:
      pass
# Serialize the XML tree to a string without escaping special characters
xml_string = ET.tostring(root, encoding="unicode")

# Replace &lt; and &gt; with < and > in the XML string
xml_string = xml_string.replace('&lt;', '<').replace('&gt;', '>')

# Write the modified XML string to the output file
with open('thread.xml', 'w', encoding="utf8") as output_file:
    # Include the XML declaration and write the modified XML string
    output_file.write(f'<?xml version="1.0" encoding="utf-8"?>\n{xml_string}')
