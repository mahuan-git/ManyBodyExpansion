import sys 
sys.path.append('/public/home/jlyang/quantum/program/MBE/src')
from mbe import geometry_h_chain, MBE
import numpy as np
def test():
    natom = 10
    mbe_energy_1=[]
    mbe_energy_2=[]
    mbe_energy_3=[]
    scf_energy = []
    for bondlength in np.arange(0.6,3.0,0.1):
        geometry = geometry_h_chain(natom,bondlength)
        natom_per_fragment=2
        n_fragment = int(np.ceil(natom/natom_per_fragment))
        fragment = []
        for i in range(n_fragment):
            if (i+1)*natom_per_fragment <=natom:
                fragment.append(np.arange(i*natom_per_fragment,(i+1)*natom_per_fragment))
            else:
                fragment.append(np.arange(i*natom_per_fragment,natom))
        mbe = MBE(geometry,fragment,solver = 'vqe',basis = 'sto-3g',periodic = False,link_atom = False)
        mbe_energy_1.append(mbe.get_mbe_energy(1))
        mbe_energy_2.append(mbe.get_mbe_energy(2))
        mbe_energy_3.append(mbe.get_mbe_energy(3))
        scf_energy.append(mbe.mf.e_tot)
    print('many body expansion calculation finished')
    print(scf_energy)
    print('mbe_energy_1:')
    print(mbe_energy_1)
    print('mbe_energy_2:')
    print(mbe_energy_2)
    print('mbe_energy_3:')
    print(mbe_energy_3)
if __name__=="__main__":
    test()
