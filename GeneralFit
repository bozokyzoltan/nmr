#!/usr/bin/env pthon
# -*- coding: utf-8 -*-
"""
Created by Zoltan Bozoky on 2013.11.08.
Under GPL licence.

A recoursive and general gradient descent algorithm using user provided error function

Purpose:
========
* Multiparamter fitting algorithm
* User defined cost function
* Plot data generation
* Fitting error calculation
* Capable for paralell computing


Usage:
======
import GeneralFit as GF
myfit = GF.Fit_general(x,y,[p1,p2],my_error_function)
print myfit.Value, myfit.Chi2

"""

import math
import threading


class Fit_general():
    """
    This is a general fitting algoritm for parameters up to 5.
    Arguments:
        x          = x data
        y          = y data
        parameters = array [1..5] of the starting values
        function   = any function calculates a y value depending on the parameters and an x value
        z          = additional data passed to the function, if needed

    Note1: It is better to overestimate the starting value, than under, because of speed.
    Note2: The order of the parameters must be the same order as the function accepts it.
    Note3: If you use more than 3 parameters, it is worth decreasing the IterationMaximum and increaseing the FittingErrorLimit

    Usage:

        def R2decay(x,parameters):
            return parameters[0]*math.exp(-x/float(parameters[1]))

        x = [0.0001,0.006,0.012]
        y = [1.11448e+07,7.50362e+06,6.64364e+06]
        I0 = 10792583.4557
        T2 = 1E-2

        FitR2 = Fit_general(x,y,[I0,T2],R2decay)
        print FitR2.Value[1]
        print FitR2.Error

        a,b = FitR2.GenerateCurve(0,0.1)
    """
    #############################
    def __init__(self,x,y,parameters,function, **kw):
        """
        Parameters:
        ===========
            * NoErrors
            * FittingErrorLimit
            * IterationMaximum
            * Sparse
            * DecreamentStep
            * Log
            * ExtendedLog
            * z
        """
        ##########################
        ### Keyword arguments ####
        ##########################
        #
        # Calculation errors need time
        fastprocess = 'NoErrors' in kw
        #
        # The FittingErrorLimit defines the valueable digits 1E-8 means 8 digits
        if 'FittingErrorLimit' in kw:
            self.Errorlimit = kw['FittingErrorLimit']
        else:
            self.Errorlimit = 1E-8
        # The maximum deepness of each iteration
        if 'IterationMaximum' in kw:
            self.IterationMax = kw['IterationMaximum']
        else:
            self.IterationMax = 1000
        # Sparse = 2 means every second element will be used during the fitting
        if 'Sparse' in kw:
            Sparse = kw['Sparse']
        else:
            Sparse = 0
        # the fit is getting closer by this step every turn
        if 'DecreamentStep' in kw:
            self.StepDecreament = kw['DecreamentStep']
        else:
            self.StepDecreament = 1E+1
        if 'Log' in kw:
            self.log = kw['Log']
        else:
            self.log = True
        if 'ExtendedLog' in kw:
            self.extendedlog = kw['ExtendedLog']
        else:
            self.extendedlog = False
        if 'z' in kw:
            z = kw['z']
        else:
            z = []
        #
        ###########################
        ###  ####  ####  ####  ####
        ###########################
        if self.log:
            print '### Fitting has started ###'
            print 'Note: Number of parameters: '+str(len(parameters))
            print 'Initial parameters:',parameters
        #
        if Sparse > 0:
            self._x = []
            self._y = []
            for i in range(len(x)):
                if i % Sparse == 0:
                    self._x.append(x[i])
                    self._y.append(x[i])
            if self.log:
                print 'Note: Data length reduction from '+str(len(x))+' to '+str(len(self._x))
        else:
            self._x = x
            self._y = y
        #
        if z != []:
            self.AdditionalParameter = True
            self._z = z
            if Sparse > 0:
                self._z = []
                for i in range(len(x)):
                    if i % (Sparse) == 0:
                        self._z.append(z[i])
        else:
            self.AdditionalParameter = False
        #
        #
        self.externalfunction    = function
        self.NumberOfParameters  = len(parameters)
        #
        self.initparameters = parameters
        #
        #########################
        ### Run the fitting  ####
        #########################
        self.results = {}
        self.errors = []
        if self.NumberOfParameters <= 7:
            if fastprocess:
                self.results[0] = self.fit(self._x,self._y,parameters,self.NumberOfParameters,self.Errorlimit,self.log,self.extendedlog)
            else:
                self.fitparameter(self._x,self._y,parameters,self.NumberOfParameters,self.Errorlimit)
            
        else:
            print 'Too many parameters...'
            exit()
        #
        if self.log:
            print '### The fitting is done... ###'
        #
        return None
    #############################
    def ErrorFunction(self,x,y,parameters):
        """
        Calculates the squared error (rmsd) based on an external function
        """
        #
        error = 0.0
        for i in range(len(x)):
            #
            if self.AdditionalParameter:
                parameters.append(self._z[i])
            try:
                previouserror = error
                error += (y[i]-self.externalfunction(x[i],parameters))**2/float(y[i])
                # if the devision causes a problem
                if math.isnan(error) or math.isinf(error):
                    error = previouserror
            except OverflowError:
                error += 9.99E+99
            except ZeroDivisionError:
                #error += 9.99E+99
                pass
            #
            if self.AdditionalParameter:
                parameters.pop()
        return [error]
    #############################
    def skipone(self,array,position):
        """
        Returns an array, length is len(array) - 1, skipping the value of the
        position position
        """
        return array[0 : position] + array[position + 1 : len(array)]
    #############################
    def mean_std(self, values):
        """
        returns the mean and the standard deviation of the values
        """
        mean = 0.0
        diff = 0.0
        for element in values:
            mean += element
            diff += element**2
        mean /= float(len(values))
        diff /= float(len(values))
        try:
            diff = math.sqrt(diff - mean**2)
        except ValueError:
            diff = math.sqrt(abs(diff - mean**2))
        return [mean, diff]
    #############################
    def fitparameter(self, data_x, data_y, parameters, number, finalerrorlimit):
        """
        """
        maxnumberofnodes = 50
        #
        number_of_nodes = len(data_x)+1
        if number_of_nodes > maxnumberofnodes:
            number_of_nodes = maxnumberofnodes
        #
        nodes = {}
        nodes[0] = threading.Thread(target = self.submit_a_fit,args = (0,data_x,data_y,parameters,number,finalerrorlimit,self.log,self.extendedlog))
        for i in range(number_of_nodes):
            nodes[i+1] = threading.Thread(target = self.submit_a_fit,args = (i+1,self.skipone(data_x,i),self.skipone(data_y,i),parameters,number,finalerrorlimit,False,False))
        # Run the fitting
        for i in range(number_of_nodes):
            nodes[i].start()
        # Wait for all the thread to be finish
        for i in range(number_of_nodes):
            nodes[i].join()
        for parameternumber in range(self.NumberOfParameters):
            parametervalues = []
            for result in self.results:
                parametervalues.append(self.results[result][parameternumber])
            self.errors.append(self.mean_std(parametervalues)[-1])

        return None
    #############################
    def submit_a_fit(self,index,data_x,data_y,parameters,fitting_level,finalerrorlimit,log,exendedlog):
        self.results[index] = self.fit(data_x,data_y,parameters,fitting_level,finalerrorlimit,log,exendedlog)
        return None
    #############################
    def fit(self,data_x,data_y,initial_parameters,fitting_level,finalerrorlimit,log,extendedlog):
        """
        The fitting process runs till the error difference between a two tested
        values are equal or less than the final error limit. The fitting_level
        specifies the number of parameters to fit in each run - This is a
        recoursive procedure calling itself with a decreased fitting_level.

        Parameters:
        ===========
            * x = x data
            * y = y data
            * initial_parameters = fitting parameters for the external function
            * fitting_level = number of fitting level, maximum is the number of
                       parameters
            * finalerrorlimit = error limit for stopping the calculations
            * log = print out parameters during the fitting
            * extendedlog = print out the fit parameters of the first child
                      function
        Returns:
        ========
            * Fitted values in the same order of initial_parameters plus
              the squared fitting error

        """
        ######################
        ### Initialization ###
        ######################
        # test_value is defined by the fitting level, order from the back
        test_value = initial_parameters[-1*fitting_level]
        # delta specifies the step going closer and closer to the fitted value
        delta = test_value / self.StepDecreament
        # error limit for the child functions
        errorlimit = finalerrorlimit
        # Storage for all the tested values
        all_tested_values = {}
        #
        if fitting_level == 1:
            # Level ground call the errorfunction
            all_tested_values[test_value] = self.ErrorFunction(data_x,data_y,initial_parameters)
        else:
            # call a child
            all_tested_values[test_value] = self.fit(data_x,data_y,initial_parameters,fitting_level-1,errorlimit,log,extendedlog)
        # iteration limit
        close_to_fit = False
        # remember the last change and if it was long time ago increase the steps
        last_change = 0
        # remember the best scored value
        best = test_value
        #
        ######################
        ###      LOOP      ###
        ######################
        iteration_number = 0
        # Iterate till either it reached the iteration limit, or we are closer to the value than the Errorlimit
        while not close_to_fit and (iteration_number < self.IterationMax):
            ###########################
            ###  TEST test_value 1  ###
            ###########################
            if not test_value in all_tested_values:
                initial_parameters[-1*fitting_level] = test_value
                if fitting_level == 1:
                    all_tested_values[test_value] = self.ErrorFunction(data_x,data_y,initial_parameters)
                else:
                    all_tested_values[test_value] = self.fit(data_x,data_y,initial_parameters,fitting_level-1,errorlimit,log,extendedlog)
                if all_tested_values[test_value][-1] < all_tested_values[best][-1]:
                    best = test_value
            ###########################
            ###  TEST test_value 2  ###
            ###########################
            if test_value + delta in all_tested_values:
                test_plusz_delta_error = all_tested_values[test_value + delta][-1]
            else:
                initial_parameters[-1*fitting_level] = test_value + delta
                if fitting_level == 1:
                    test_plusz_delta_error = self.ErrorFunction(data_x,data_y,initial_parameters)[-1]
                else:
                    test_plusz_delta_error = self.fit(data_x,data_y,initial_parameters,fitting_level-1,errorlimit,log,extendedlog)[-1]
            ####################################
            ###  test and setup the limits  ###
            ####################################
            try:
                if abs(all_tested_values[test_value][-1] - test_plusz_delta_error) / float(all_tested_values[test_value][-1]) < finalerrorlimit:
                    close_to_fit = True
                else:
                    if errorlimit < finalerrorlimit:
                        errorlimit = finalerrorlimit
                    else:
                        errorlimit = (abs(all_tested_values[test_value][-1] - test_plusz_delta_error) / 100.0) / float(all_tested_values[test_value][-1])
            except ZeroDivisionError:
                close_to_fit = True
            #######################
            ###  log if needed  ###
            #######################
            if log:
                if self.NumberOfParameters == fitting_level:
                    print '{0:4d})'.format(iteration_number+1),test_value,
                    for value in all_tested_values[test_value]:
                        print value,
                    print ''
                if extendedlog and self.NumberOfParameters == fitting_level + 1:
                    print '    {0:2d}:{1:4d})'.format(fitting_level,iteration_number+1),test_value,
                    for value in all_tested_values[test_value]:
                        print value,
                    print ''
            #########################
            ###  Step up or down  ###
            #########################
            if test_plusz_delta_error < all_tested_values[test_value][-1]:
                change = delta
            else:
                change = -1*delta
            ###########################################
            ###  Change the step or the test value  ###
            ###########################################
            if test_value + change in all_tested_values:
                delta /= self.StepDecreament
                test_value = best
                # Reset the counter of the last change
                last_change = 0
            else:
                test_value += change
                # count the last change
                last_change += 1
                # Increase the step if the the last change was ong time ago
                if last_change == 15:
                    delta *= self.StepDecreament
                    last_change = 0
            #
            iteration_number += 1
        #########################
        ### Return the values ###
        #########################
        all_tested_values[best].insert(0,best)
        return all_tested_values[best]
    #############################
    def CalculateZ(self, x):
        i = 0
        while (i<len(self._x)) and (self._x[i]<= x):
            i += 1
        if (0 <= i) and (i < len(self._x)-1):
            zz = self._z[i] + (((x-self._x[i])*(self._z[i+1]-self._z[i]))/(self._x[i+1]-self._x[i]))
        elif x < self._x[0]:
            zz = self._z[0]
        else:
            zz = self._z[-1]
        return zz
    #############################
    def GenerateCurve(self, xmin, xmax, NumberOfSteps = 1000):
        x = []
        y = []
        for i in range(NumberOfSteps):
            x.append(xmin + i*(xmax-xmin)/float(NumberOfSteps-1))
            parameters = []
            for v in self.Value:
                parameters.append(v)
            if self.AdditionalParameter:
                parameters.append(self.CalculateZ(x[-1]))
                
            y_value = self.externalfunction(x[-1],parameters)
            if abs(y_value) > 1E+150:
                y_value = 1E+150
            y.append(y_value)
        return x,y
    #############################
    def GenerateCurveEx(self, xmin, xmax, externalparameters, NumberOfSteps = 1000):
        x = []
        y = []
        for i in range(NumberOfSteps):
            x.append(xmin + i*(xmax-xmin)/float(NumberOfSteps-1))
            #
            y.append(self.externalfunction(x[-1],externalparameters))
            #
        return x,y
    #############################
    def GetRMSD(self):
        try:
            error = math.sqrt(self.results[0][-1])
        except ValueError:
            error = 9.99E+99
        return error
    #############################
    def GetValue(self):
        return self.results[0][0:-1]
    #############################
    def GetError(self):
        return self.errors
    #############################
    def GetChi2(self):
        """
        chi2 = sum((Observed - Expected)^2/ Expected)
        """
        error = 0.0
        for i in range(len(self._x)):
            #
            try:
                observed = self._y[i]
                # 
                parameters = self.Value
                if self.AdditionalParameter:
                    parameters.append(self._z[i])
                    
                expected = self.externalfunction(self._x[i],parameters)

                if self.AdditionalParameter:
                    parameters.pop()
                #
                error += (observed-expected)**2 / expected
                #
            except OverflowError:
                error += 9.99E+99
            except ZeroDivisionError:
                pass
            #
        return error
    #############################
    Error = property(GetError)
    Value = property(GetValue)
    Chi2  = property(GetChi2)
#######################################################################################

