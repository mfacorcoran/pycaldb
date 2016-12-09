__author__ = 'corcoran'


def read_model_xcm(xcmo_file):
    """
    This function reads a model file created wiht the xspec "save model <filename>" command
    and returns a pyxspec model object
    @param xcmo_file:
    @return:
    """
    import xspec
    par = []
    modelfound = False
    with open(xcmo) as file:
        for line in file:
            # print line
            if 'model' in line:
                xspecmodel = xspec.Model(line[5:].strip())
                modelfound = True
            if modelfound:
                par.append(line.strip('\n').strip())
    print model
    par = par[1:]
    for i in arange(len(par)) + 1:
        xspecmodel(i).values = par[i - 1]
    return xspecmodel