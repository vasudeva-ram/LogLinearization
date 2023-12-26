#!/usr/bin/env python
# coding: utf-8

# Some preliminaries
import numpy as np
from collections import namedtuple
import sympy as smp
import re
from sympy.abc import _clash, _clash2
from sympy.printing import pprint
from IPython.display import display, Math
from IPython.core.debugger import set_trace
from typing import Any, List, Optional, NamedTuple
smp.init_session(use_unicode=True)
set(_clash)


# # General Convenience Classes

class VarSet(NamedTuple):
    """ Defines all the auxiliary vars created
    for each endogenous variable in the economy
    """    
    varx : str
    backVar : str
    forwVar : str
    ssVar : str
    hatVar : str
    hatbackVar : str
    hatforVar : str


class Variable(smp.Symbol):
    def __init__(self, 
                 strval: str
                ) -> None:
        self.modname: str = strval
        super()


class VarFunction(smp.Function):
    def __init__(self, 
                 strval: str
                ) -> None:
        self.modname: str = strval
        super()


## ModEconomy Class


class ModEconomy(object):
    def __init__(self, 
                 varDict: dict, #contains names of variables and their tex representations
                 paramDict: dict, #contains names of parameters and their tex representations
                 ssvals: dict, # contains steady state values of variables
                 funcDict: dict = None #contains names of all externally defined functions and their tex reps.
                ) -> None:
        
        self._allvarDict = self._make_symbol_vars(varDict)
        self._varList = self._make_varList(self.allvarDict)
        self._ssDict = self._make_ssDict(self.allvarDict)
        self._hatDict = self._make_hatDict(self.allvarDict)
        self._funcList = self._make_funcs(funcDict) if funcDict is not None else None
        self._params = self._make_paramList(paramDict)
        self._simpvals = self._make_simpvals(ssvals)
        self._eqsymbols = self._make_eqsymbols()
        #self._make_properties()
        
        
    def _make_symbol_vars(self,
                          varDict: dict
                         ) -> dict:
        """
        Making symbols of all relevant variables, their forward and backward states, 
        steady states and log deviations
        """
        allvarDict = {}
        for item in varDict.items():
            varz, names = self._make_allvars(item)
            singlevardict = {}
            for key in varz.keys():
                sym = Variable(varz[key])
                sym.name = names[key]
                setattr(self, sym.modname, sym) # makes each variable callable as a property
                singlevardict.update({key:sym})
            allvarDict.update({item[0]:singlevardict})

        return allvarDict
    
    def _make_allvars(self, 
                      item):
        varx = item[0]
        vardict = {'varx' : varx,
                   'backVar' : varx+'bb',
                    'forwVar' : varx+'ff',
                    'ssVar' : 'ss'+varx,
                    'hatVar' : 'hat'+varx,
                    'hatbackVar' : 'hat'+varx+'bb',
                    'hatforVar' : 'hat'+varx+'ff'
                  }

        varname = item[1]
        namedict = {'varx' : varname+'_t',
                   'backVar' : varname + r'_{t-1}',
                    'forwVar' : varname + r'_{t+1}',
                    'ssVar' : r'\bar{' + varname + r'}',
                    'hatVar' : r'\hat{' + varname + r'}_t',
                    'hatbackVar' : r'\hat{' + varname + r'}_{t-1}',
                    'hatforVar' : r'\hat{' + varname + r'}_{t+1}'
                   }

        return vardict, namedict
    
    
    def _make_varList(self,
                     allvarDict: dict
                    ) -> list:
        """Makes a list of all variables, defined as {t}, {t-1}
        and {t+1} states of each variable 
        """
        varList = []
        for key in allvarDict.keys():
            varx = allvarDict[key]['varx']
            backVar = allvarDict[key]['backVar']
            forwVar = allvarDict[key]['forwVar']
            varList.extend([varx, backVar, forwVar])

        return varList
        
    def _make_ssDict(self,
                    allvarDict: dict) -> dict:
        """ Makes a dictionary mapping each variable (as defined in varList)
        to its steady state representation
        """

        ssDict = {}
        for key in allvarDict.keys():
            varx = allvarDict[key]['varx']
            backVar = allvarDict[key]['backVar']
            forwVar = allvarDict[key]['forwVar']
            ssvar = allvarDict[key]['ssVar']
            ssDict.update({varx:ssvar, backVar:ssvar, forwVar: ssvar})

        return ssDict
    
    def _make_hatDict(self,
                      allvarDict: dict
                     ) -> dict:
        """ Makes a dictionary mapping all variables (as defined in varList)
        with their log deviation representations
        """
        hatDict = {}
        for key in allvarDict.keys():
            varx = allvarDict[key]['varx']
            backVar = allvarDict[key]['backVar']
            forwVar = allvarDict[key]['forwVar']
            hatVar = allvarDict[key]['hatVar']
            hatbackVar = allvarDict[key]['hatbackVar']
            hatforVar = allvarDict[key]['hatforVar'] 
            hatDict.update({varx:hatVar, backVar:hatbackVar, forwVar:hatforVar})

        return hatDict
    
    def _make_paramList(self,
                        paramDict: dict
                       ) -> list:
        """Makes a list of all symbolized parameters
        """
        paramlist = []
        for key in paramDict.keys():
            prm = Variable(key)
            prm.name = paramDict[key]
            setattr(self, prm.modname, prm)
            paramlist.append(prm)
            
        return paramlist
    
    def _make_simpvals(self,
                       ssvals: dict
                      ) -> dict:
        """
        Takes a dictionary of steady state values and prepares it
        for simplification purposes later.
        """
        binaries = {key:val for (key, val) in ssvals.items() if (val == 0 or val == 1)}
        avd = self.allvarDict
        simp_vars = {key:val for (key,val) in binaries.items() if (key in avd.keys())}
        simpvals = {avd[key]['ssVar']:val for (key, val) in simp_vars.items()}
        
        return simpvals
    
    def _make_funcs(self,
                    funcDict: dict
                   ) -> list:
        funcList = []
        for fn in funcDict.keys():
            fnvar = VarFunction(fn)
            fnvar.name = funcDict[fn]
            setattr(self, fn, fnvar) # make function accessible by property
            funcList.append(fnvar)
            
        return funcList
    
    def _make_eqsymbols(self) -> dict:
        """Namespace dictionary of all relevant variables, parameters and functions
        likely to be encountered by an Equation object
        """
        avd = self.allvarDict
        eqsyms = {}
        
        for key in avd.keys():
            for subkey in avd[key].keys():
                sym = avd[key][subkey]
                eqsyms.update({sym.modname: sym})
                
        for prm in self.params:
            eqsyms.update({prm.modname: prm})

        if self.funcList is not None:
            for fn in self.funcList:
                eqsyms.update({fn.modname: fn})
                
        return eqsyms
                        

    @property
    def allvarDict(self):
        return self._allvarDict
    
    @property
    def varList(self):
        return self._varList
    
    @property
    def params(self):
        return self._params
    
    @property
    def ssDict(self):
        return self._ssDict
    
    @property
    def hatDict(self):
        return self._hatDict
    
    @property
    def simpvals(self):
        return self._simpvals
    
    @property
    def funcList(self):
        return self._funcList
    
    @property
    def eqsymbols(self):
        return self._eqsymbols
    
    


## Equation Class

class Equation(object):
    def __init__(self,
                 expr: str, #expression from dynare file with LHS and RHS (remove the ";" at the end)
                 mod: ModEconomy # ModEconomy object with variables, parameters and ss-values
                ) -> None:
        self.mod = mod
        self._spliteqn: dict = self._make_splits(expr)
        self._expr: smp.Equality = smp.Eq(self.spliteqn['lhs'], self.spliteqn['rhs'], evaluate=False)
        self._ss = self._make_ss()
    
        
    def _make_splits(self,
                     expr: str
                    ) -> dict:
        timefix = expr.replace("(+1)", "ff").replace("(-1)", "bb").replace("^", "**")
        l, r = timefix.split("=")
        mod = self.mod
        lhs = smp.sympify(l, mod.eqsymbols)
        rhs = smp.sympify(r, mod.eqsymbols)
        
        return {'lhs': lhs, 'rhs': rhs}
    
    
    def loglin(self) -> smp.Expr:
        """ Returns the loglinearized version of an equation
        """
        ssv_lhs, ll_lhs = self._single_loglin(self.spliteqn['lhs'])
        ssv_rhs, ll_rhs = self._single_loglin(self.spliteqn['rhs'])
        
        simp_l = self._simplifier(ssv_lhs, ll_lhs)
        simp_r = self._simplifier(ssv_rhs, ll_rhs)
        
        return smp.Eq(simp_l, simp_r, evaluate=False)
        
    def _simplifier(self,
                    ssval: Any,
                    ll: Any
                   ) -> smp.Expr:
        """ Simplifies loglinearized equations using steady state
        expressions and other Sympy solvers
        """
        simp_vals = self.mod.simpvals
        hats = self.mod.hatDict.values()
        
        stage1 = ll - ssval
        stage2 = stage1.subs(simp_vals)
        stage3 = smp.expand(stage2, mul=True)
        stage4 = smp.powsimp(stage3, force=True)
        stage5 = smp.simplify(stage4)
        final = smp.collect(stage5, hats)
        
        return final
               
    
    def _single_loglin(self,
                       expr: Any
                      ) -> Any:
        
        mod = self.mod
        varList = mod.varList
        ssDict = mod.ssDict
        hatDict = mod.hatDict
        avd = mod.allvarDict
        
        # getting first term ready
        ssval = smp.log(expr).subs(ssDict)
        # exp_ss = smp.exp(all_ss, evaluate=False)

        terms = [ssval]
        for endog in varList:
            if endog in expr.free_symbols:
                tempdict = {key:value for (key,value) in ssDict.items() if key!=endog}
                single_var = expr.subs(tempdict)
                all_ss = smp.expand_log(smp.log(single_var), force=True)
                taylor = all_ss.series(endog, ssDict[endog], n=2).removeO()
                hatit = taylor.subs(endog - ssDict[endog], hatDict[endog] * ssDict[endog])
                clean = smp.logcombine(hatit, force=True)
                term = clean.coeff(hatDict[endog])
                term = smp.simplify(term)
                terms.append(term * hatDict[endog])
            else:
                continue

        return ssval, smp.Add(*tuple(terms))  
    
    def _make_ss(self):
        """ Returns the steady state version of the equation
        """
        expr = self._expr.subs(self.mod.ssDict)
        ss1 = smp.simplify(expr)
        ss2 = smp.expand(ss1, mul=True)
        ss3 = smp.powsimp(ss2, force=True)
        ss = smp.simplify(ss3)
        
        return ss
        
    
    @property
    def spliteqn(self):
        return self._spliteqn
    
    @property
    def ss(self):
        return self._ss
    
    @property
    def expr(self):
        return self._expr
    
        
        
        


## ModFile Parser

def clean_varblock(varblock: list) -> dict:
    clean = {}
    for i in varblock:
        if i.startswith(('var','varexo','parameters')):
            i = i.split(' ', 1)[1]
        tst = i.strip().split('(')[0]
        fin = re.split(' +', tst)
        pyname = fin[0]
        sname = fin[1].split('$')[1]
        clean.update({pyname:sname})

    return clean


def get_blocks(filename: str,
               block_type: str # can be "var ", "varexo " or "parameters "
              ) -> list:
    
    with open(filename) as f:
        lines = [line for line in f.readlines() if line.strip()]

        block_lines = []
        l_iter = iter(lines)
        for line in l_iter:
            if line.startswith(("%", "//")):
                continue
            
            if line.startswith(block_type):
                while ';' not in line:
                    if line.startswith(("%", "//")):
                        pass
                    else:
                        block_lines.append(line)
                    line = next(l_iter)
                break
        
    block_dict = clean_varblock(block_lines)

    return block_dict


def get_eqns(filename: str,
             block_type: str # can be "var ", "varexo ", "model" or "parameters "
             ) -> list:
    
    with open(filename) as f:
        lines = [line for line in f.readlines() if line.strip()]

        eqns = []
        l_iter = iter(lines)
        for line in l_iter:
            if line.startswith("%") or line.startswith("//"):
                continue
            
            if line.startswith(block_type):
                while 'end;' not in line:
                    if line.startswith(("%", "//", "[", "#")):
                        pass
                    else:
                        eqns.append(line)
                    line = next(l_iter)
                break

        return eqns

def cleaned_eqns(filename: str,
                 block_type: str # can be "var ", "varexo " or "parameters "
                 ) -> list:
    
    eqns = get_eqns(filename, block_type)
    clean_eqns = []
    for line in eqns:
        eq = line.split("//")[0]
        eq = eq.replace(";", "").replace('\n', "")
        clean_eqns.append(eq)
        
    return clean_eqns[1:]


def createModEconomy(filename: str) -> ModEconomy:
    varonlyDict = get_blocks(filename, "var ")
    varexoDict = get_blocks(filename, "varexo ")
    varDict = {**varonlyDict, **varexoDict}
    paramDict = get_blocks(filename, "parameters ")
    ssvals = get_blocks(filename, "steady_state_model ")
    # funcDict = get_blocks(filename, "functions ")
    mod = ModEconomy(varDict, paramDict, ssvals)
    mod_eqns = cleaned_eqns(filename, "model")
    eqns = [Equation(eq, mod) for eq in mod_eqns]
    
    return mod, eqns

