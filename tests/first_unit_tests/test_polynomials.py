import numpy as np

import openmc


def test_zernike_radial():
    coeff = np.asarray([1.3, -3.0, 9e-1, -6e-1, 0.11])
    zn_rad = openmc.ZernikeRadial(coeff)
    assert zn_rad.order == 8
    assert zn_rad.radius == 1

    coeff = np.asarray([1.3, -3.0, 9e-1, -6e-1, 0.11, 0.222])
    zn_rad = openmc.ZernikeRadial(coeff, 0.392)
    assert zn_rad.order == 10
    assert zn_rad.radius == 0.392
    norm_vec = (2 * np.arange(6) + 1) / (np.pi * 0.392 ** 2)
    norm_coeff = norm_vec * coeff

    rho = 0.5
    # Reference solution from running the Fortran implementation
    raw_zn = np.array([
        1.00000000e+00, -5.00000000e-01, -1.25000000e-01,
        4.37500000e-01, -2.89062500e-01, -8.98437500e-02])

    ref_vals = np.sum(norm_coeff * raw_zn)

    test_vals = zn_rad(rho * zn_rad.radius)

    assert ref_vals == test_vals

    rho = [0.2, 0.5]
    # Reference solution from running the Fortran implementation
    raw_zn1 = np.array([
        1.00000000e+00, -9.20000000e-01, 7.69600000e-01,
        -5.66720000e-01, 3.35219200e-01, -1.01747000e-01])
    raw_zn2 = np.array([
        1.00000000e+00, -5.00000000e-01, -1.25000000e-01,
        4.37500000e-01, -2.89062500e-01, -8.98437500e-02])

    ref_vals = [np.sum(norm_coeff * raw_zn1), np.sum(norm_coeff * raw_zn2)]

    test_vals = zn_rad([i * zn_rad.radius for i in rho])

    assert np.allclose(ref_vals, test_vals)


def test_zernike():
    import openmc.lib as lib
     
    coeff = np.asarray([1.1e-1, -3.2e2, 5.3, 7.4, -9.5, 0.005])
    zn_azimuthal = openmc.Zernike(coeff)
    assert zn_azimuthal.order == 2
    assert zn_azimuthal.radius == 1

    coeff = np.asarray([1.5, -3.6, 9.7e-1, -6.8e-1, 0.11, 0.33e2, 0.002, 13.75, 
                        3.1, -7.3, 7.8e-1, -1.1e-1, 2.56, 5.25e3, 0.123])
    zn_azimuthal = openmc.Zernike(coeff, 0.392)
    assert zn_azimuthal.order == 4
    assert zn_azimuthal.radius == 0.392
    norm_vec = np.array([1, 4, 4, 6, 3, 6, 8, 8, 8, 8, 
                        10, 10, 5, 10, 10]) / (np.pi * 0.392 ** 2)
    norm_coeff = norm_vec * coeff 
    
    rho = 0.5
    
    theta = np.radians(45) 
    # Reference solution from running the C API for calc_zn
    raw_zn = lib.calc_zn(zn_azimuthal.order, rho, theta)
 
    ref_vals = np.sum(norm_coeff * raw_zn)

    test_vals = zn_azimuthal(rho * zn_azimuthal.radius, theta)

    assert ref_vals == test_vals

    rho = [0.2, 0.5]
    
    theta = np.radians(30) 
    #Reference solution from running the C API for calc_zn
    raw_zn1 = lib.calc_zn(zn_azimuthal.order, rho[0], theta)
    
    raw_zn2 = lib.calc_zn(zn_azimuthal.order, rho[1], theta)
 
    ref_vals = [np.sum(norm_coeff * raw_zn1), np.sum(norm_coeff * raw_zn2)]

    test_vals = zn_azimuthal([i * zn_azimuthal.radius for i in rho], theta)
    
    assert np.allclose(ref_vals, test_vals)    
    
    rho = 0.2
    
    theta = np.radians([30, 60]) 
    #Reference solution from running the C API for calc_zn
    raw_zn1 = lib.calc_zn(zn_azimuthal.order, rho, theta[0])
    
    raw_zn2 = lib.calc_zn(zn_azimuthal.order, rho, theta[1])

    ref_vals = [np.sum(norm_coeff * raw_zn1), np.sum(norm_coeff * raw_zn2)]

    test_vals = zn_azimuthal(rho * zn_azimuthal.radius, [j for j in theta])
    
    assert np.allclose(ref_vals, test_vals)    
    
    rho = [0.2, 0.5]
    
    theta = np.radians([30, 60]) 
    #Reference solution from running the C API for calc_zn
    raw_zn1 = lib.calc_zn(zn_azimuthal.order, rho[0], theta[0])
    
    raw_zn2 = lib.calc_zn(zn_azimuthal.order, rho[1], theta[0])
    
    raw_zn3 = lib.calc_zn(zn_azimuthal.order, rho[0], theta[1])
    
    raw_zn4 = lib.calc_zn(zn_azimuthal.order, rho[1], theta[1])

    ref_vals = [np.sum(norm_coeff * raw_zn1), np.sum(norm_coeff * raw_zn2),
                np.sum(norm_coeff * raw_zn3), np.sum(norm_coeff * raw_zn4)]

    test_vals = zn_azimuthal([i * zn_azimuthal.radius for i in rho], [j for j in theta])
    
    test_vals = np.ravel(test_vals)    
    
    assert np.allclose(ref_vals, test_vals)    
    

