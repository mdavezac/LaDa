""" NLEP package """


def fit(system_params=None):
    import mpopt
    print "fit called"
    mpopt.main(system_params)

    
