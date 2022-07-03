import unittest
import numpy as np

def car2sph(R):
    """ Cartesian to spherical coordinates """
    R = np.asarray(R)

    if np.size(R) > 3:
        R = np.transpose(R)

    r = np.linalg.norm(R, axis=0)
    rho = np.linalg.norm(R[0:2], axis=0)
    X, Y, Z = [R[0], R[1], R[2]]
    R_polar = [r,
               np.arctan2(rho, Z),
               np.arctan2(Y, X)]

    return np.asarray(R_polar).T

class TestStringMethods(unittest.TestCase):
    def test_field(self):
        import kmag
        print('KAMG imported successfully.')
        # print(kmag.__doc__)
        # import matplotlib.pyplot as plt

        # Inputs
        R = [10, 0, 0]      # Test location in RS units
        ETIME = 1483228800  # (01 Jan 2017)
        EPOCH = 'j2000'     # ('j2000' | 'ctime')
        BY_IMF = -0.2
        BZ_IMF = 0.1
        Dp = 0.017
        IN_COORD =  'DIS'   # S3C/DIS/DIP/KSM/KSO
        OUT_COORD = 'DIS'   # S3C/DIS/DIP/KSM/KSO

        # Rotate position vector to S3C/DIS spherical coordinates
        R_S3C_CAR = kmag.krot(IN_COORD, 'S3C', R, ETIME, EPOCH)
        R_S3C = car2sph(R_S3C_CAR)

        # Return the magnetic field in S3C/DIS spherical coordinates
        LT, BR, BTH, BPHI = kmag.kmag(ETIME, EPOCH, R_S3C[0], R_S3C[1], R_S3C[2], BY_IMF, BZ_IMF, Dp)

        # Rotate the magnetic field vector to the requested coordinates
        BX, BY, BZ = kmag.sph2car_mag(BR, BTH, BPHI, R_S3C[1], R_S3C[2])
        B = kmag.krot('S3C', OUT_COORD, [BX, BY, BZ], ETIME, EPOCH)

        # Reference Magnetic Field value
        ref_B = [ -2.98050199, -0.09383225, -15.72154196]
        ref_B = np.round(np.asarray(ref_B), 6)
        B = np.round(np.asarray(B), 6)
        print(f'Input Position, R = {R} RS')
        print(f'Magnetic Field, B = {B} nT')
        print(f'Reference Field, B = {ref_B} nT')

        self.assertTrue(np.allclose(B, ref_B))

if __name__ == '__main__':
    unittest.main()
