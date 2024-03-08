import sys
sys.path.append('/public/home/jlyang/quantum/program/MBE/src')
from mbe import MBE
import numpy as np
def geometry_carbon_ring(shift=20.0):
    nat = 18
    shift =shift*2*np.pi/360
    R = 7.31/2
    geometry = []
    angle = 0.0
    for i in range( nat // 2 ):
        geometry.append(('C', (R * np.cos(angle        ), R * np.sin(angle        ), 0.0)))
        geometry.append(('C', (R * np.cos(angle + shift), R * np.sin(angle + shift), 0.0)))
        angle += 4.0 * np.pi / nat
    return geometry

def test():
    mbe_energy_1=[]
    mbe_energy_2=[]
    mbe_energy_3=[]
    scf_energy = []
    mbe_energy_4=[]
    natom = 18
    for shift in np.arange(18.0,18.31,0.1):
        geometry = geometry_carbon_ring(shift)
        natom_per_fragment=2
        n_fragment = int(np.ceil(natom/natom_per_fragment))
        fragment = []
        for i in range(n_fragment):
            if (i+1)*natom_per_fragment <=natom:
                fragment.append(np.arange(i*natom_per_fragment,(i+1)*natom_per_fragment))
            else:
                fragment.append(np.arange(i*natom_per_fragment,natom))
        mbe = MBE(geometry,fragment,solver = 'ccsd',basis = 'cc-pvdz',periodic = True, link_atom =True)
        mbe_energy_1.append(mbe.get_mbe_energy(1))
        mbe_energy_2.append(mbe.get_mbe_energy(2))
        mbe_energy_3.append(mbe.get_mbe_energy(3))
        mbe_energy_4.append(mbe.get_mbe_energy(4))
        scf_energy.append(mbe.mf.e_tot)
    print('many body expansion calculation finished')
    print(scf_energy)
    print('mbe_energy_1:')
    print(mbe_energy_1)
    print('mbe_energy_2:')
    print(mbe_energy_2)
    print('mbe_energy_3:')
    print(mbe_energy_3)
    print('mbe_energy_4:')
    print(mbe_energy_4)
if __name__=="__main__":
    test()
