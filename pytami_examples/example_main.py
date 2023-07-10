import torch
import pytami
import examples as ex
import time

def main():
    example_2()
    example_4()
    example_6()

def example_2():

    print("\n\n-_-_-_ Example - Second Order _-_-_-")

    # start example
    print("\n-----Constructing Fermi Tree format -----\n")
    
    # torch device to perform calculations on
    device = torch.device("cuda")

    # class instance
    ami = pytami.TamiBase(device)

    # second order self energy setup
    R0 = ex.construct_example2()

    # external variables
    avars = ex.construct_ext_example2(ami)

    # timing info for setup
    t1 = time.time()

    # make ft_terms object
    ftout = pytami.TamiBase.ft_terms()

    # Integration/Evaluation parameters
    E_REG = 0 # numberical regulator for small energies.  If inf/nan results try E_REG=1e-8 
    N_INT = int(2) # number of matsubara sums to perform 
    test_amiparms = pytami.TamiBase.ami_parms(N_INT, E_REG) # SHOULD BE (0, 0) ?

    # now construct!
    ami.construct(N_INT, R0, ftout)

    # time construct leg and start evaluate time
    t2 = time.time()

    # Evaluate the integrand for ext parms in avars
    calc_result = ami.evaluate(test_amiparms, ftout, avars)
    
    # time the end of eval
    t_end = time.time()

    # microsecond result
    diff1 = (t2 - t1) * 1000000
    diff2 = (t_end - t2) * 1000000

    # print results
    torch.set_printoptions(precision=10)
    print(f"Result was {calc_result}")
    print(f"Construction took {diff1} microseconds")
    print(f"Evaluation took {diff2} microseconds")

    # end example

def example_4():

    print("\n\n-_-_-_ Example - Fourth Order _-_-_-")

    # start example
    print("\n-----Constructing Fermi Tree format -----\n")
    
    # torch device to perform calculations on
    device = torch.device("cuda")

    # class instance
    ami = pytami.TamiBase(device)

    # second order self energy setup
    R0 = ex.construct_example4()

    # external variables
    avars = ex.construct_ext_example4(ami)

    # timing info for setup
    t1 = time.time()

    # make ft_terms object
    ftout = pytami.TamiBase.ft_terms()

    # Integration/Evaluation parameters
    E_REG = 0 # numberical regulator for small energies.  If inf/nan results try E_REG=1e-8 
    N_INT = 4 # number of matsubara sums to perform 
    test_amiparms = pytami.TamiBase.ami_parms(N_INT, E_REG) # SHOULD BE (0, 0) ?

    # now construct!
    ami.construct(N_INT, R0, ftout)

    # time construct leg and start evaluate time
    t2 = time.time()

    # Evaluate the integrand for ext parms in avars
    calc_result = ami.evaluate(test_amiparms, ftout, avars)
    
    # time the end of eval
    t_end = time.time()

    # microsecond result
    diff1 = (t2 - t1) * 1000000
    diff2 = (t_end - t2) * 1000000

    # print results
    torch.set_printoptions(precision=10)
    print(f"Result was {calc_result}")
    print(f"Construction took {diff1} microseconds")
    print(f"Evaluation took {diff2} microseconds")

def example_6():

    print("\n\n-_-_-_ Example - Sixth Order _-_-_-")

    # start example
    print("\n-----Constructing Fermi Tree format -----\n")
    
    # torch device to perform calculations on
    device = torch.device("cuda")

    # class instance
    ami = pytami.TamiBase(device)

    # second order self energy setup
    R0 = ex.construct_example6()

    # external variables
    avars = ex.construct_ext_example6(ami)

    # timing info for setup
    t1 = time.time()

    # make ft_terms object
    ftout = pytami.TamiBase.ft_terms()

    # Integration/Evaluation parameters
    E_REG = 0 # numberical regulator for small energies.  If inf/nan results try E_REG=1e-8 
    N_INT = 6 # number of matsubara sums to perform 
    test_amiparms = pytami.TamiBase.ami_parms(N_INT, E_REG) # SHOULD BE (0, 0) ?

    # now construct!
    ami.construct(N_INT, R0, ftout)

    # time construct leg and start evaluate time
    t2 = time.time()

    # Evaluate the integrand for ext parms in avars
    calc_result = ami.evaluate(test_amiparms, ftout, avars)
    
    # time the end of eval
    t_end = time.time()

    # microsecond result
    diff1 = (t2 - t1) * 1000000
    diff2 = (t_end - t2) * 1000000

    # print results
    torch.set_printoptions(precision=10)
    print(f"Result was {calc_result}")
    print(f"Construction took {diff1} microseconds")
    print(f"Evaluation took {diff2} microseconds")


if __name__ == '__main__':
    main()
