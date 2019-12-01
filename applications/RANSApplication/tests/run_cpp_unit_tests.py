from KratosMultiphysics import *
from KratosMultiphysics.RANSApplication import *


def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS 
    Tester.RunTestSuite("KratosRansFastSuite")
    #Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS)
    Tester.RunTestSuite("KratosRansFastSuite_1")
if __name__ == '__main__':
    run()
