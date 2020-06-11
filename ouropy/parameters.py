# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 14:32:33 2017

@author: DanielM
"""


class Parameter(object):

    def __init__(self, mech_name, sec_name, value):
        self.mech_name = str(mech_name)
        self.sec_name = str(sec_name)
        self.value = float(value)
        self._i = 0

    def get_mech(self):
        split = self.mech_name.split('_')
        if len(split) < 2:
            return None
        else:
            return split[1]

    def __iter__(self):
        return self

    def next(self):
        if self._i == 0:
            return self
        else:
            self._i = 0
            raise StopIteration()

    def __repr__(self):
        return ("Parameter('%s','%s',%f)"
                % (self.mech_name, self.sec_name, self.value))

    def __str__(self):
        return ("mech_name: %s, sec_name: %s, value: %d"
                % (self.mech_name, self.sec_name, self.value))


class ParameterSet(list):

    def __init__(self, param_list):
        self.param_list = param_list
        self._i = 0

    def get_mechs(self):
        mechs = {}
        for x in self.param_list:
            if x.sec_name not in mechs.keys():
                mechs[x.sec_name] = set()
            """if x.get_mech() and (not (x.get_mech() in mechs[x.sec_name])):
                mechs[x.sec_name].append(x.get_mech())"""
            if x.get_mech():
                mechs[x.sec_name].add(x.get_mech())

        return mechs

    def __str__(self):
        mylist = [x.__str__() for x in self.param_list]
        return ';'.join(map(str, mylist))

    def __iter__(self):
        return self

    def __next__(self):
        if self._i < (len(self.param_list)):
            i = self._i
            self._i += 1
            return self.param_list[i]
        else:
            self._i = 0
            raise StopIteration()

    def next(self):
        return self.__next__()

def read_parameters(path):

    reader = open(path, 'r')

    data = reader.read()

    first_split = data.split('\n')
    second_split = [x.split('\t') for x in first_split]
    second_split = [x for x in second_split if bool(x[0])]

    parameter_list = [Parameter(x[0], x[1], x[2]) for x in second_split]
    return ParameterSet(parameter_list)
