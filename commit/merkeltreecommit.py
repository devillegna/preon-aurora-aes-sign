

from Crypto.Random import get_random_bytes

def _randombytes( bsize ) :
    return get_random_bytes( bsize )


LAMBDA = 256

from Crypto.Hash import SHA3_256 as G

def commit( msgList ):
    num = len(msgList)
    assert( 2 <= num )
    assert( 0==(num&(num-1)) )   # len(msgList) is a power of 2
    r = [ _randombytes( LAMBDA//8 ) for i in range(num) ]
    mktree = [ [ G.new(b''.join([msgList[i] , r[i]])).digest() for i in range(num) ]  ]
    while( num > 2 ):
        last_layer = mktree[-1]
        mktree.append([ G.new( b''.join([last_layer[i*2],last_layer[i*2+1]]) ).digest() for i in range(num>>1) ])
        num = num//2
    last_layer = mktree[-1]
    rt = G.new( b''.join([last_layer[0],last_layer[1]]) ).digest()
    return rt, r , mktree

def open( msg , idx , r , mktree ):
    _idx = idx
    auth_path = [ msg , r[idx] ]
    for layer in mktree :
        auth_path.append( layer[idx^1] )
        idx = idx//2
    return auth_path

def verify( rt , idx , auth_path ):
    state = G.new( b''.join([auth_path[0],auth_path[1]]) ).digest()
    for i in range(2,len(auth_path)):
        if (idx&1) : state = G.new( b''.join([auth_path[i],state]) ).digest()
        else       : state = G.new( b''.join([state,auth_path[i]]) ).digest()
        idx = idx//2
    return state == rt