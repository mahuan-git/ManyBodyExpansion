'''
Many Body Expansion 
'''
import numpy as np
from itertools import combinations
from pyscf import gto, scf

#from solver import pyscf_uhf
#class combinations_new(combinations):
    
class MBE():
    def __init__(   self,
                    geometry: list,
                    fragment: list,   # information about atom indices every fragment holds
                    solver  : str = 'uhf',
                    basis : str = 'sto_3g',
                    periodic: bool = False,
                    link_atom :bool = False   #used when bond break and need to saturate bonds
                    ):
        self.geometry = geometry
        self.fragment = fragment
        self.solver = solver.lower()
        self.mbe_series = []
        self.mbe_energy=None
        self.mbe_1=None
        self.mbe_2=None
        self.mbe_3=None
        self.mbe_4=None
        self.mbe_5=None
        self.basis = basis
        self.periodic = periodic
        self.link_atom=link_atom
        mol=gto.Mole()
        mol.atom=self.geometry
        if mol.nelectron%2==0:
            mol.spin=0
        else:
            mol.spin=1
        mol.basis = self.basis
        mol.build()
        mf = scf.UHF(mol)
        mf.verbose = 3
        mf.max_cycle = 1000
        mf.scf(dm0=None)
        self.mf = mf

    def fractorial(self,n : int):
        if n==1 or n==0:
            return 1
        else:
            return n*self.fractorial(n-1)

    def get_mbe_energy_approx( self,
                        n:int, # MBE order
                        ):
        calculated_mbe_order = len(self.mbe_series)
        for i in np.arange(calculated_mbe_order,n):
            print('Calculating the order %d mbe energy'%(i))
            self.get_n_th_order_mbe(n=i)

    def get_mbe_energy( self,
                        order :int =2
                        ):
        self.mbe_series=[]
        if order >5:
            print('err: order larger than 5 not supported')
        elif order ==1:
            if self.mbe_1==None:
                self.get_mbe_1()
            #self.mbe_series.append(sum(self.mbe_1))
            mbe_energy=sum(self.mbe_1)
        elif order ==2:
            if self.mbe_1==None:
                self.get_mbe_1()
            if self.mbe_2==None:
                self.get_mbe_2()
            #self.mbe_series.append(sum(self.mbe_2)-(len(self.fragment)-2)*sum(self.mbe_1))
            mbe_energy=sum(self.mbe_2)-(len(self.fragment)-2)*sum(self.mbe_1)
        elif order ==3:
            if self.mbe_1==None:
                self.get_mbe_1()
            if self.mbe_2==None:
                self.get_mbe_2()
            if self.mbe_3==None:
                self.get_mbe_3()
            #self.mbe_series.append(sum(self.mbe_2)-(len(self.fragment)-2)*sum(self.mbe_1))
            mbe_energy=sum(self.mbe_3)-(len(self.fragment)-3)*sum(self.mbe_2)\
                        +0.5*(len(self.fragment)-2)*(len(self.fragment)-3)*sum(self.mbe_1)
        elif order ==4:
            if self.mbe_1==None:
                self.get_mbe_1()
            if self.mbe_2==None:
                self.get_mbe_2()
            if self.mbe_3==None:
                self.get_mbe_3()
            if self.mbe_4==None:
                self.get_mbe_4()
            #self.mbe_series.append(sum(self.mbe_2)-(len(self.fragment)-2)*sum(self.mbe_1))
            mbe_energy=sum(self.mbe_4)-(len(self.fragment)-4)*sum(self.mbe_3)\
                        +0.5*(len(self.fragment)-3)*(len(self.fragment)-4)*sum(self.mbe_2)\
                            -(1/6)*(len(self.fragment)-2)*(len(self.fragment)-3)*(len(self.fragment)-4)*sum(self.mbe_1)
        elif order ==5:
            if self.mbe_1==None:
                self.get_mbe_1()
            if self.mbe_2==None:
                self.get_mbe_2()
            if self.mbe_3==None:
                self.get_mbe_3()
            if self.mbe_4==None:
                self.get_mbe_4()
            if self.mbe_5==None:
                self.get_mbe_5()
            #self.mbe_series.append(sum(self.mbe_2)-(len(self.fragment)-2)*sum(self.mbe_1))
            mbe_energy=sum(self.mbe_5)-(len(self.fragment)-5)*sum(self.mbe_4)\
                        +0.5*(len(self.fragment)-4)*(len(self.fragment)-5)*sum(self.mbe_3)\
                            -(1/6)*(len(self.fragment)-3)*(len(self.fragment)-4)*(len(self.fragment)-5)*sum(self.mbe_2)\
                                +(1/24)*(len(self.fragment)-2)*(len(self.fragment)-3)*(len(self.fragment)-4)*(len(self.fragment)-5)*sum(self.mbe_1)
        self.mbe_energy = mbe_energy
        return(mbe_energy)


    def get_n_th_order_mbe( self,
                            n:int, #MBE order
                            ):
        if n==0:
            mbe_energy_0=0
            for i in range(len(self.fragment)):
                mbe_energy_0+=self.get_energy(self.fragment[i])
            self.mbe_series.append(mbe_energy_0)
        elif n>len(self.fragment):
            print('exceeds highest order')
            exit()
        elif n>len(self.mbe_series):
            print('lower order energy needed')
            exit()
        else:
            num_fragment = len(self.fragment)
            mbe_energy=0
            for atom_list in combinations(self.fragment,n+1):
                print(atom_list)
                atom_list_re = []
                for frag in atom_list:
                    for atom_idx in frag:
                        atom_list_re.append(atom_idx)
                print(atom_list_re)
                mbe_energy += self.get_energy(atom_list_re)
            for i in range(n):
                print(self.fractorial(num_fragment-i-1))
                mbe_energy -= self.mbe_series[i]*self.fractorial(num_fragment-i-1)/(self.fractorial(num_fragment-n-1)*self.fractorial(n-i))
                print(mbe_energy)
            self.mbe_series.append(mbe_energy)
    
    def get_mbe_1(self):
        self.mbe_1=[]
        if self.periodic == True:
            energy_1 = self.get_energy(self.fragment[0])
            for i in range(len(self.fragment)):
                self.mbe_1.append(energy_1)
        else:
            for i in range(len(self.fragment)):
                self.mbe_1.append(self.get_energy(self.fragment[i]))
    
    def get_mbe_2(self):
        self.mbe_2=[]
        if self.periodic == True:
            mbe_2_tmp=[]
            for i in np.arange(1,len(self.fragment)):
                atom_list = []
                for atom_idx in self.fragment[0]:
                    atom_list.append(atom_idx)
                for atom_idx in self.fragment[i]:
                    atom_list.append(atom_idx)
                #print(self.fragment[0])
                #print(self.fragment[i])
                #print(atom_list)
                mbe_2_tmp.append(self.get_energy(atom_list))
            for i in range(len(self.fragment)-1):
                for j in np.arange(i,len(self.fragment)-1):
                    self.mbe_2.append(mbe_2_tmp[j])
            #print(self.mbe_2)
            #exit()
        else:
            for atom_list in combinations(self.fragment,2):
                #print(atom_list)
                atom_list_re = []
                for frag in atom_list:
                    for atom_idx in frag:
                        atom_list_re.append(atom_idx)
                #print(atom_list_re)
                self.mbe_2.append(self.get_energy(atom_list_re))
    
    def get_mbe_3(self):
        self.mbe_3=[]
        for atom_list in combinations(self.fragment,3):
            #print(atom_list)
            atom_list_re = []
            for frag in atom_list:
                for atom_idx in frag:
                    atom_list_re.append(atom_idx)
            #print(atom_list_re)
            self.mbe_3.append(self.get_energy(atom_list_re))

    def get_mbe_4(self):
        self.mbe_4=[]
        for atom_list in combinations(self.fragment,4):
            #print(atom_list)
            atom_list_re = []
            for frag in atom_list:
                for atom_idx in frag:
                    atom_list_re.append(atom_idx)
            #print(atom_list_re)
            self.mbe_4.append(self.get_energy(atom_list_re))
    def get_mbe_5(self):
        self.mbe_5=[]
        for atom_list in combinations(self.fragment,5):
            #print(atom_list)
            atom_list_re = []
            for frag in atom_list:
                for atom_idx in frag:
                    atom_list_re.append(atom_idx)
            #print(atom_list_re)
            self.mbe_5.append(self.get_energy(atom_list_re))
    
    def get_energy(self,atom_list):
        if self.solver=='uhf':
            from solver import pyscf_uhf
            energy = pyscf_uhf( geometry=self.geometry,
                                atom_list=atom_list,
                                basis = self.basis,
                                link_atom=self.link_atom)
        elif self.solver=='rhf':
            from solver import pyscf_rhf
            energy = pyscf_rhf( geometry=self.geometry,
                                atom_list=atom_list,
                                basis = self.basis,
                                link_atom=self.link_atom)
        elif self.solver=='vqe':
            from solver import run_vqechem
            energy=run_vqechem( geometry=self.geometry,
                                atom_list=atom_list,
                                basis = self.basis,
                                link_atom=self.link_atom)
        elif self.solver == 'ccsd':
            from solver import pyscf_ccsd
            energy = pyscf_ccsd( geometry=self.geometry,
                                atom_list=atom_list,
                                basis = self.basis,
                                link_atom=self.link_atom)
        elif self.solver == 'dmrg' or self.solver == 'chemps2' :
            from solver import chemps2
            energy = chemps2(geometry=self.geometry,
                            atom_list=atom_list,
                            basis = self.basis,
                            link_atom=self.link_atom)
        else:
            print('can not find solver, please check solver option')
        
        return energy

def geometry_h_ring(nat:int = 10,  #number of hydrogen atoms in the ring
                    bondlength : float=1.0
                    ):
    geometry = []
    r = 0.5 * bondlength / np.sin(np.pi/nat)
    for i in range(nat):
        theta = i * (2*np.pi/nat)
        geometry.append(('H', (r*np.cos(theta), r*np.sin(theta), 0)))
    return geometry

def geometry_h_chain(nat:int = 10,  #number of hydrogen atoms in the ring
                    bondlength : float=1.0
                    ):
    geometry = []
    for i in range(nat):
        geometry.append(('H', (0, 0, i*bondlength)))
    return geometry

def geometry_Be_ring(nat:int = 30,  #number of atoms in the ring
                    bondlength : float=2.0
                    ):
    geometry = []
    r = 0.5 * bondlength / np.sin(np.pi/nat)
    for i in range(nat):
        theta = i * (2*np.pi/nat)
        geometry.append(('Be', (r*np.cos(theta), r*np.sin(theta), 0)))
    return geometry

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
    natom = 10
    geometry = geometry_h_ring(natom)
    natom_per_fragment=2
    n_fragment = int(np.ceil(natom/natom_per_fragment))
    fragment = []
    for i in range(n_fragment):
        if (i+1)*natom_per_fragment <=natom:
            fragment.append(np.arange(i*natom_per_fragment,(i+1)*natom_per_fragment))
        else:
            fragment.append(np.arange(i*natom_per_fragment,natom))
    mbe = MBE(geometry,fragment,solver = 'uhf')
    return mbe
    mbe.get_n_th_order_mbe(0)
    print(mbe.mbe_series)
    
#if __name__=="__main__":
#    test()