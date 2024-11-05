#!/usr/bin/env python
from string import Template
import argparse

parser = argparse.ArgumentParser(description='Creates an XML file making use of a template and parameters provided by user for Maestro')
parser.add_argument("-pdb_path", type=str, required=True)
parser.add_argument("-prefix", type=str, default="pdb")
parser.add_argument("-postfix", type=str, default=".ent.gz")
parser.add_argument("-tolower", type=str, default="true")
parser.add_argument("-toupper", type=str, default="false")
parser.add_argument("-bu", type=str, default="false")

args = parser.parse_args()

pdb_path = args.pdb_path
prefix = args.prefix
postfix = args.postfix
tolower = args.tolower
toupper = args.toupper
bu = args.bu

xml_template = '''<?xml version="1.0"?>
<sef>
    <evaluate>
        <pdbhome dir="${pdb_path}" prefix="${prefix}" postfix="${postfix}" tolower="${tolower}" toupper="${toupper}" bu="${bu}"/>
        <!-- <pdbhome dir="/home/bill/biounits/" postfix=".pdb1.gz" tolower="true" bu="true"/> -->
    </evaluate>
    <!--  "Do not modify anything beyond this line unless you really know what you do"  -->
    <energy_functions>
        <energy_function class="PairedSEF" file="[CONFIG_ROOT]/effiles/128944.ef.bin" format="bin" use="stab"/>
        <energy_function class="PairedSEF" file="[CONFIG_ROOT]/effiles/137642.ef.bin" format="bin" use="stab"/>
        <energy_function class="ContactSEF" file="[CONFIG_ROOT]/effiles/92318.ef.bin" format="bin" use="stab">
            <evaluate>
                <score_values adjust="zscore" scale="1pm"/>
            </evaluate>
        </energy_function>
        <energy_function class="PairedSEF" file="[CONFIG_ROOT]/effiles/55787.ef.bin" format="bin" use="struct"/>
        <energy_function class="PairedSEF" file="[CONFIG_ROOT]/effiles/58894.ef.bin" format="bin" use="struct"/>
        <energy_function class="ContactSEF" file="[CONFIG_ROOT]/effiles/237.ef.bin" format="bin" use="struct"/>
    </energy_functions>
    <ddg_council mode="default">
        <pca filename="[CONFIG_ROOT]/council/pca.xml"/>
        <additional_input>nresidue,asa,mass,hydrophilicity,isoelectric_point,secstr,</additional_input>
        <scale file="[CONFIG_ROOT]/council/scale.xml">mean_sdev</scale>
        <adjust file="[CONFIG_ROOT]/council/member_adjustment.xml" method="none"/>
        <add_constant n="1" value="1"/>
        <ph>true</ph>
        <cmode>1</cmode>
        <members>
            <neuralnet_master file="[CONFIG_ROOT]/council/neural_nets.xml" netpath="[CONFIG_ROOT]/council/neuralnets/" pca="false" specialist="no"/>
            <neuralnet_master file="[CONFIG_ROOT]/council/stab_neural_nets.xml" netpath="[CONFIG_ROOT]/council/neuralnets/stab/" pca="false" specialist="stab"/>
            <neuralnet_master file="[CONFIG_ROOT]/council/destab_neural_nets.xml" netpath="[CONFIG_ROOT]/council/neuralnets/destab/" pca="false" specialist="destab"/>
            <svm file="[CONFIG_ROOT]/council/svm.model" pca="false" specialist="no"/>
            <svm file="[CONFIG_ROOT]/council/stab_svm.model" pca="false" specialist="stab"/>
            <svm file="[CONFIG_ROOT]/council/destab_svm.model" pca="false" specialist="destab"/>
            <regression comb_score="false" file="[CONFIG_ROOT]/council/regression.xml" pca="false" specialist="no"/>
        </members>
    </ddg_council>
        <ddg_council mode="ssbond">
        <pca filename="[CONFIG_ROOT]/council/ss/pca.xml"/>
        <additional_input>nresidue,asa,mass,hydrophilicity,isoelectric_point,secstr,</additional_input>
        <scale file="[CONFIG_ROOT]/council/ss/scale.xml">mean_sdev</scale>
        <adjust file="[CONFIG_ROOT]/council/ss/member_adjustment.xml" method="none"/>
        <add_constant n="1" value="1"/>
        <ph>true</ph>
        <cmode>1</cmode>
        <members>
            <neuralnet_master file="[CONFIG_ROOT]/council/ss/neural_nets.xml" pca="false" specialist="no"/>
            <neuralnet_master file="[CONFIG_ROOT]/council/ss/stab_neural_nets.xml" pca="false" specialist="stab"/>
            <neuralnet_master file="[CONFIG_ROOT]/council/ss/destab_neural_nets.xml" pca="false" specialist="destab"/>
            <svm file="[CONFIG_ROOT]/council/ss/svm.model" pca="false" specialist="no"/>
            <svm file="[CONFIG_ROOT]/council/ss/stab_svm.model" pca="false" specialist="stab"/>
            <svm file="[CONFIG_ROOT]/council/ss/destab_svm.model" pca="false" specialist="destab"/>
            <regression comb_score="false" file="[CONFIG_ROOT]/council/ss/regression.xml" pca="false" specialist="no"/>
        </members>
    </ddg_council>
</sef>'''
template = Template(xml_template)
xml = template.substitute(pdb_path = pdb_path, prefix = prefix, postfix = postfix, tolower = tolower, bu = bu)
print(xml)