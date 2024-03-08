'''
solvers for mbe
'''
from pyscf import gto, scf
from pyscf.cc import ccsd
from itertools import combinations
import numpy as np
import os

def get_distance(coordinate_1,coordinate_2):
    assert(len(coordinate_1)==len(coordinate_2)==3)
    distance = 0
    for i in range(len(coordinate_1)):
        distance +=(coordinate_1[i]-coordinate_2[i])**2
    distance = np.sqrt(distance)
    return distance

def get_link_atom_coordinate(   geometry : list,  
                                end_atoms : list,
                                add_atom_idx : int,
                                mode : str = 'extend',
                                ):
    '''
    get the coordinate of the added link atom
    geometry : the geometry of the original molecule
    end_atoms: the index of the two atom at the end of a fragment
    add_atom_idx : the atom index of the added atom.
    mode : two mode to add a link atom
    '''   
    natoms = len(geometry)    
    end_atoms = list(end_atoms)
    end_atoms.sort()     
    bondlength = 1.0
    link_atom_coordinate = [0,0,0]
    if mode == 'extend':
        coordinate_1 = geometry[end_atoms[0]%natoms][1]
        coordinate_2 = geometry[end_atoms[1]%natoms][1]
        distance = get_distance(coordinate_1,coordinate_2)
        delta_x = (coordinate_2[0]-coordinate_1[0])*bondlength/distance
        delta_y = (coordinate_2[1]-coordinate_1[1])*bondlength/distance
        delta_z = (coordinate_2[2]-coordinate_1[2])*bondlength/distance
        #print([delta_x,delta_y,delta_z])
        if add_atom_idx == (end_atoms[1]+1):
            link_atom_coordinate[0] = coordinate_2[0] + delta_x
            link_atom_coordinate[1] = coordinate_2[1] + delta_y
            link_atom_coordinate[2] = coordinate_2[2] + delta_z
        elif add_atom_idx == (end_atoms[0]-1):
            link_atom_coordinate[0] = coordinate_1[0] - delta_x
            link_atom_coordinate[1] = coordinate_1[1] - delta_y
            link_atom_coordinate[2] = coordinate_1[2] - delta_z
        else:
            print('add atom index does not match, check please')
            exit()

    elif mode =='origin':
        '''to be finished later
        '''
        pass
    else:
        print('mode to add the link atom not recognized, please check')
        exit()
    
    if (False):
        print('coordinates:')
        print(coordinate_1)
        print(coordinate_2)
        print(link_atom_coordinate)
        exit()
    return link_atom_coordinate


def add_link_atoms( geometry,
                    atom_list,
                    mode='extend' 
                    ):
    '''
    Do not use link_atom if there is only one atom in the fragment
    mode can either be extend or origin
    extend 为附加原子位于两原子的延长线上
    origin 为附加原子位于分子原有原子的位置附近
    '''
    natoms = len(geometry)
    link_atom_coordinate=[]
    for atoms in combinations(atom_list,2):
        if abs(atoms[1]-atoms[0])==1:   ## two selected atoms are neighbours
            atom_idx_1 = max(atoms)+1
            atom_idx_2 = min(atoms)-1
            if atom_idx_1%natoms not in atom_list:
                coordinate = get_link_atom_coordinate(geometry,atoms,atom_idx_1,mode = mode)
                link_atom_coordinate.append(coordinate)
            if atom_idx_2%natoms not in atom_list:
                coordinate = get_link_atom_coordinate(geometry,atoms,atom_idx_2,mode = mode)
                link_atom_coordinate.append(coordinate)

        elif abs(atoms[1]-atoms[0])==natoms-1:   ##the two selected atoms are at the head and the tail of the molecule respectively
            distance = get_distance(geometry[atoms[1]][1],geometry[atoms[0]][1])
            if distance >3.0:  ## indicate that the head and the tail are not connected
                pass
            else:          ## the head and the tail are connected. 
                atom_idx_1 = 1
                atom_idx_2 = -2
                atoms = [-1,0]
                if (atom_idx_1%natoms) not in atom_list:
                    coordinate = get_link_atom_coordinate(geometry,atoms,atom_idx_1,mode = mode)
                    link_atom_coordinate.append(coordinate)
                if (atom_idx_2%natoms) not in atom_list:
                    coordinate = get_link_atom_coordinate(geometry,atoms,atom_idx_2,mode = mode)
                    link_atom_coordinate.append(coordinate)
        else: ## the two atoms chosen are not connected
            pass
    
    return link_atom_coordinate

def pyscf_uhf(  geometry:list,
                atom_list:list,
                basis:str='sto-3g',
                link_atom : bool = False):
    mol=gto.Mole()
    mol.atom=[]
    for i in atom_list:
        mol.atom.append(geometry[i])
    if link_atom == True:
        H_atom_coordinates = add_link_atoms(geometry,atom_list,mode = 'extend')
        for coordinate in H_atom_coordinates:
            mol.atom.append(('H',coordinate))
    #print('molecule after link H atoms are added')
    #print(mol.atom)
    if mol.nelectron%2==0:
        mol.spin=0
    else:
        mol.spin=1
    mol.basis = basis
    mol.build()
    mf = scf.UHF(mol)
    mf.verbose = 3
    mf.max_cycle = 1000
    mf.scf(dm0=None)
    return mf.e_tot

def pyscf_rhf(  geometry:list,
                atom_list:list,
                basis:str='sto-3g',
                link_atom :bool= False):
    mol=gto.Mole()
    mol.atom=[]
    for i in atom_list:
        mol.atom.append(geometry[i])
    if link_atom == True:
        H_atom_coordinates = add_link_atoms(geometry,atom_list,mode = 'extend')
        for coordinate in H_atom_coordinates:
            mol.atom.append(('H',coordinate))
    if mol.nelectron%2==0:
        mol.spin=0
    else:
        mol.spin=1
    mol.basis = basis
    mol.build()
    mf = scf.RHF(mol)
    mf.verbose = 3
    mf.max_cycle = 1000
    mf.scf(dm0=None)
    return mf.e_tot


def run_vqechem(geometry:list,
                atom_list:list,
                basis: str ='sto-3g',
                link_atom : bool = False):
    mol=gto.Mole()
    mol.atom=[]
    for i in atom_list:
        mol.atom.append(geometry[i])
    if link_atom == True:
        H_atom_coordinates = add_link_atoms(geometry,atom_list,mode = 'extend')
        for coordinate in H_atom_coordinates:
            mol.atom.append(('H',coordinate))

    if mol.nelectron%2==0:
        mol.spin=0
    else:
        mol.spin=1
    mol.basis = basis
    mol.build()
    import sys
    sys.path.append('/public/home/jlyang/quantum/program/vqechem/src')
    from algorithms import run_vqe
    #from scf_from_pyscf import 
    import scipy
    import math 

    options = {
               'vqe' : {'algorithm':'adapt-vqe'},
               'scf' : {'ncas':None,'ncore':None,'shift':0.5},
               'ops' : {'class':'fermionic','spin_sym':'sa'},
               'ansatz' : {'method':'adapt','form':'taylor','Nt':10},
               'opt' : {'maxiter':300}
              }
    ansatz = run_vqe(mol,options)
    return ansatz._energy

def pyscf_ccsd( geometry:list,
                atom_list:list,
                basis:str='sto-3g',
                link_atom : bool = False):
    mol=gto.Mole()
    mol.atom=[]
    for i in atom_list:
        mol.atom.append(geometry[i])
    if link_atom == True:
        H_atom_coordinates = add_link_atoms(geometry,atom_list,mode = 'extend')
        for coordinate in H_atom_coordinates:
            mol.atom.append(('H',coordinate))
    #print('molecule after link H atoms are added')
    #print(mol.atom)
    if mol.nelectron%2==0:
        mol.spin=0
    else:
        mol.spin=1
    mol.basis = basis
    mol.build()
    mf = scf.RHF(mol)
    mf.verbose = 3
    mf.max_cycle = 1000
    mf.scf(dm0=None)
    ccsolver = ccsd.CCSD( mf )
    ccsolver.verbose = 5
    ECORR, t1, t2 = ccsolver.ccsd()
    ERHF = mf.hf_energy
    ECCSD = ERHF + ECORR
    return ECCSD


def chemps2(  geometry:list,
                atom_list:list,
                basis:str='sto-3g',
                link_atom : bool = False):
    import sys
    sys.path.append('/public/home/jlyang/quantum/program/vqechem/src')
    from scf_from_pyscf import PySCF, pyscf_interface , get_1e_integral, get_2e_integral
    import ctypes
    import PyCheMPS2
    mol=gto.Mole()
    mol.atom=[]
    for i in atom_list:
        mol.atom.append(geometry[i])
    if link_atom == True:
        H_atom_coordinates = add_link_atoms(geometry,atom_list,mode = 'extend')
        for coordinate in H_atom_coordinates:
            mol.atom.append(('H',coordinate))
    #print('molecule after link H atoms are added')
    #print(mol.atom)
    if mol.nelectron%2==0:
        mol.spin=0
    else:
        mol.spin=1
    mol.basis = basis
    mol.build()
    
    options =   {
                         'mapping'        : 'JW',
                         'ncas'           : None,
                         'ncore'          : None,
                         'mo_list'        : None,
                         'shift'          : 0
                }

    hf = pyscf_interface(mol,options)
    FOCK = hf._h1
    CONST = hf._Enuc
    TEI = get_2e_integral(hf)
    Norb = hf._nmo
    Nel = sum(hf._mol.nelec)

    Initializer = PyCheMPS2.PyInitialize()
    Initializer.Init()

    # Setting up the Hamiltonian
    Group = 0
    orbirreps = np.zeros([ Norb ], dtype=ctypes.c_int)
    HamCheMPS2 = PyCheMPS2.PyHamiltonian(Norb, Group, orbirreps)
    HamCheMPS2.setEconst( CONST )
    for cnt1 in range(Norb):
        for cnt2 in range(Norb):
            HamCheMPS2.setTmat(cnt1, cnt2, FOCK[cnt1, cnt2])
            for cnt3 in range(Norb):
                for cnt4 in range(Norb):
                    HamCheMPS2.setVmat(cnt1, cnt2, cnt3, cnt4, TEI[cnt1, cnt3, cnt2, cnt4]) #From chemist to physics notation
    '''HamCheMPS2.save()
    exit(123)'''
    # Killing output if necessary
    # to be fixed later
    if ( False ):
        sys.stdout.flush()
        old_stdout = sys.stdout.fileno()
        new_stdout = os.dup(old_stdout)
        devnull = os.open('/dev/null', os.O_WRONLY)
        os.dup2(devnull, old_stdout)
        os.close(devnull)
    #if ( Norb <= 10 ):
    if (Norb <10):
    #if (False):
        # FCI ground state calculation
        #this part will raise Segmentation fault, the fault is in: /public/home/jlyang/quantum/anaconda3/envs/vqechem/lib/python3.7/site-packages/PyCheMPS2.cpython-37m-x86_64-linux-gnu.so ,this part is going to be fixed later, or skip this part
        assert( Nel % 2 == 0 )
        Nel_up       = Nel / 2
        Nel_down     = Nel / 2
        Irrep        = 0
        maxMemWorkMB = 100.0
        FCIverbose   = 2
        
        theFCI = PyCheMPS2.PyFCI( HamCheMPS2, Nel_up, Nel_down, Irrep, maxMemWorkMB, FCIverbose )
        GSvector = np.zeros( [ theFCI.getVecLength() ], dtype=ctypes.c_double )
        theFCI.FillRandom( theFCI.getVecLength() , GSvector ) # Random numbers in [-1,1[
        GSvector[ theFCI.LowestEnergyDeterminant() ] = 12.345 # Large component for quantum chemistry
        print('start FCI.GSDavidson')
        
        EnergyCheMPS2 = theFCI.GSDavidson( GSvector )
        
        print('end FCI.GCDavidson')
        #SpinSquared = theFCI.CalcSpinSquared( GSvector )
        TwoRDM = np.zeros( [ Norb**4 ], dtype=ctypes.c_double )
        
        ############################
        theFCI.Fill2RDM( GSvector, TwoRDM )# the source of segmentation fault
        ############################
        
        TwoRDM = TwoRDM.reshape( [Norb, Norb, Norb, Norb], order='F' )
        TwoRDM = np.swapaxes( TwoRDM, 1, 2 ) #From physics to chemistry notation
        del theFCI
    else:
    
        # DMRG ground state calculation
        assert( Nel % 2 == 0 )
        TwoS  = 0
        Irrep = 0
        Prob  = PyCheMPS2.PyProblem( HamCheMPS2, TwoS, Nel, Irrep )

        OptScheme = PyCheMPS2.PyConvergenceScheme(3) # 3 instructions
        #OptScheme.setInstruction(instruction, D, Econst, maxSweeps, noisePrefactor)
        OptScheme.setInstruction(0,  500, 1e-10,  3, 0.05)
        OptScheme.setInstruction(1, 1000, 1e-10,  3, 0.05)
        OptScheme.setInstruction(2, 1000, 1e-10, 10, 0.00) # Last instruction a few iterations without noise
        theDMRG = PyCheMPS2.PyDMRG( Prob, OptScheme )
        EnergyCheMPS2 = theDMRG.Solve()

    return EnergyCheMPS2