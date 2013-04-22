#!/usr/bin/env python
"""
                               xml_icd

a collection of classes and functions for parsing the XML produced by SALT's
TCS ICD.

Author                     Version             Date
--------------------------------------------------------
TE Pickering                 0.2             20121012

TODO
--------------------------------------------------------
will need to expand as more of the ICD is made available

Updates
--------------------------------------------------------
20121112 - switched from xml.dom.minidom to lxml.etree for XML parsing.
           etree is MUCH cleaner and easier to use.

"""
import xml
import urllib2
from lxml import etree


class ICD_EW:
    """
    ICD_EW(choices, index)

    Create an ICD_EW object from an array of choices and an index. The choices
    and index come from an <EW> Enum tag in SALT's XML ICD.  An example:
    <EW>
      <Name>tcs mode</Name>
      <Choice>OFF</Choice>
      <Choice>INITIALISE</Choice>
      <Choice>READY</Choice>
      <Choice>SLEW</Choice>
      <Choice>TRACK</Choice>
      <Choice>GUIDE</Choice>
      <Choice>MAINTENANCE</Choice>
      <Choice>CALIBRATION</Choice>
      <Choice>MAJOR FAULT</Choice>
      <Choice>SHUTDOWN</Choice>
      <Val>6</Val>
    </EW>

    Parameters
    ----------
    choices : list
        A list of choices pulled from an <EW> tag in the XML ICD.
    index : int
        The <Val> of the <EW> denoting the index of the current choice

    Provides
    -------
    ICD_EW.choices - list of choices (as input)
    ICD_EW.index - index of current choice (as input)
    ICD_EW.val - choice denoted by provided index

    """
    def __init__(self, choices, index):
        self.choices = choices
        self.index = index
        if type(choices) == list:
            if index < len(choices):
                self.val = choices[index]
            else:
                self.val = None
        else:
            self.val = None


def safeType(x, typ):
    """
    need some way to take a type as an argument and safely transform
    first argument into it
    """
    try:
        return typ(x)
    except:
        return None


def parseElement(e):
    """
    almost every element has a Name and a Val so pull these out and
    handle null Vals
    """
    k = e.find("Name").text
    if e.find("Val") is not None:
        v = e.find("Val").text
    else:
        v = None
    return (k, v)


def parseICD(url="http://sgs.salt/xml/salt-tcs-icd.xml"):
    """
    parser to take the XML ICD and turn it into a dict of clusters.
    each cluster in turn is a dict of values within the cluster.
    the values can be in one of six different data types:
       U32 - mapped to a python int
       DBL - mapped to a python float
       String - left as python unicode
       Boolean - mapped to python bool
       EW - mapped to a ICD_EW object defined here
       Array - mapped to a python list
    """
    doc = etree.parse(urllib2.urlopen(url, timeout=3))

    # the root element will be a Cluster containing all other clusters
    root = doc.getroot()

    # main dict we'll populate and return
    tcs = {}

    # easy types where we just turn a string into the type we want, safely
    types = ["U32", "DBL", "String", "Boolean"]
    lambdas = [lambda x: safeType(x, int),
               lambda x: safeType(x, float),
               lambda x: safeType(x, str),
               lambda x: safeType(safeType(x, int), bool),
               ]
    simples = zip(types, lambdas)

    # get clusters. first is the root cluster so leave it and take the rest.
    clusters = [c for c in root.iter("Cluster")]

    # loop through each cluster
    for cluster in clusters:
        cls_name = cluster.find("Name").text
        tcs[cls_name] = {}

        # pull out the complex datatypes we handle differently.
        lists = cluster.findall("EW")
        arrays = cluster.findall("Array")

        # go through the simple data types
        for s in simples:
            (typ, func) = s
            tags = cluster.findall(typ)
            for t in tags:
                (key, val) = parseElement(t)
                tcs[cls_name][key] = func(val)

        # loop through the EW elements and make them into ICD_EW objects
        for l in lists:
            (key, val) = parseElement(l)
            choices = []
            for c in l.findall("Choice"):
                choices.append(c.text)
            tcs[cls_name][key] = ICD_EW(choices, safeType(val, int))

        # loop through the arrays
        for a in arrays:
            (key, tag) = parseElement(a)
            vals = []
            for s in simples:
                (typ, func) = s
                tags = a.findall(typ)
                for t in tags:
                    (k, v) = parseElement(t)
                    vals.append(func(v))

            tcs[cls_name][key] = vals

    # the BMS returns temperatures as an array.  let's make this a dict
    # so we know what's what.
    temp_map = ["2m", "5m", "10m", "15m", "20m", "25m", "30m"]
    try:
        temps = tcs['bms external conditions']['Temperatures']
        tcs['bms external conditions']['Temperatures'] = dict(zip(temp_map,
                                                                  temps))
    except:
        pass
    return tcs
