
import sys
sys.path.insert(0, '../' )

from utils import randombytes as rd


from Crypto.Cipher import AES
#from Crypto.Util   import Counter

def preon_keygen():
  aes128_key = rd.randombytes( 16 )
  aes128_plaintext = rd.randombytes( 16 )
  aes128_ciphertext = AES.new( aes128_key , AES.MODE_ECB ).encrypt( aes128_plaintext )
  return (aes128_plaintext , aes128_ciphertext ), (aes128_plaintext , aes128_key)


from aesR1CS import aes128R1CS_z as R1CSz
from aesR1CS import r1cs
import numpy as np
from utils import hash as H
from aurora import m_aurora as aurora


def preon_sign( sk , mesg ):
    z = R1CSz.get_vec_z( *sk ,64)
    MATs = r1cs.get_r1cs( R1CSz.aes128R1CS_num_constrains , len(z) , 64, 128)  # AES128
    mat_a , mat_b , mat_c = np.array(MATs['A']) , np.array(MATs['B']) , np.array(MATs['C'])
    h_state = H.gen( bytes([1]) , mesg )
    sig = aurora.generate_proof( (mat_a,mat_b,mat_c,z,64) , h_state )
    return sig

def preon_verify( pk , sig , mesg ):
    z = R1CSz.get_vec_z( *pk , 64 )
    z = R1CSz.get_vec_1v( *pk , len(z) )
    MATs = r1cs.get_r1cs( R1CSz.aes128R1CS_num_constrains , len(z) , 64, 128)  # AES128
    mat_a , mat_b , mat_c = np.array(MATs['A']) , np.array(MATs['B']) , np.array(MATs['C'])
    h_state = H.gen( bytes([1]) , mesg )
    return aurora.verify_proof( sig , (mat_a,mat_b,mat_c,z[:64],64) , h_state )
