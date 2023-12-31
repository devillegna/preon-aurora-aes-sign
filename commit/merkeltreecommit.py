

import sys
sys.path.insert(0, '../' )

from utils import randombytes as utrd


LAMBDA = 256

#from Crypto.Hash import SHA3_256 as G
import hashlib as G

def commit( msgList ):
    num = len(msgList)
    assert( 2 <= num )
    assert( 0==(num&(num-1)) )   # len(msgList) is a power of 2
    r = [ utrd.randombytes( LAMBDA//8 ) for i in range(num) ]
    mktree = [ [ G.sha3_256( msgList[i]+r[i] ).digest() for i in range(num) ]  ]
    while( num > 2 ):
        last_layer = mktree[-1]
        mktree.append([ G.sha3_256( last_layer[i*2]+last_layer[i*2+1] ).digest() for i in range(num>>1) ])
        num = num//2
    last_layer = mktree[-1]
    rt = G.sha3_256( b''.join([last_layer[0],last_layer[1]]) ).digest()
    return rt, r , mktree

def open( msg , idx , r , mktree ):
    _idx = idx
    auth_path = [ msg , r[idx] ]
    for layer in mktree :
        auth_path.append( layer[idx^1] )
        idx = idx//2
    return auth_path

def batchopen( idxes , mesgs , r , mktree ):
    return [ open(mesgs[idx], idx , r , mktree) for idx in idxes ]

def verify( rt , idx , auth_path ):
    state = G.sha3_256( auth_path[0]+auth_path[1] ).digest()
    for i in range(2,len(auth_path)):
        if (idx&1) : state = G.sha3_256( auth_path[i]+state ).digest()
        else       : state = G.sha3_256( state+auth_path[i] ).digest()
        idx = idx//2
    return state == rt

def batchverify( idxes , rt , auth_paths ):
    vv = [ verify(rt,idx,auth_paths[i]) for i,idx in enumerate(idxes) ]
    return all(vv)


if '__main__' == __name__ :
   rt, r , mktree = commit( [ b'123' , b'456' , b'789' , b'012' ] )
   auth_path = open( b'012' , 3 , r , mktree )
   print( "PASS?" , verify( rt , 3 , auth_path ) )
