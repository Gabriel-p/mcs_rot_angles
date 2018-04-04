
from astropy.io import ascii


def main(rpath):
    '''
    Read files with data values for each cluster for each galaxy.
    '''
    # smc_data = ascii.read(rpath + '/input/valid_sol/smc_input_' + zvar + '.dat', delimiter=' ')
    # lmc_data = ascii.read(rpath + '/input/valid_sol/lmc_input_' + zvar + '.dat', delimiter=' ')

    smc_data = ascii.read(rpath + '/input/smc_input.dat', delimiter=' ')
    lmc_data = ascii.read(rpath + '/input/lmc_input.dat', delimiter=' ')

    return smc_data, lmc_data
