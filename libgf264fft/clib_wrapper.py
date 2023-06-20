import ctypes
import pathlib


_clib = ctypes.CDLL( pathlib.Path().resolve()/ "../polyeval/btfy/gf264/libgf264btfy.so" )

#uint64_t index_to_gf264( uint64_t fft_position_index );
_clib.index_to_gf264.argtypes = [ ctypes.c_uint64 ]
_clib.index_to_gf264.restype = ctypes.c_uint64

#uint64_t gf264_mul( uint64_t a , uint64_t b );
_clib.gf264_mul.argtypes = [ ctypes.c_uint64 , ctypes.c_uint64 ]
_clib.gf264_mul.restype = ctypes.c_uint64

#void btfy_64( uint64_t * poly , unsigned log_polylen , uint64_t scalar_a );
_clib.btfy_64.argtypes = [ ctypes.c_void_p , ctypes.c_uint , ctypes.c_uint64 ]
_clib.btfy_64.restype = None

#void ibtfy_64( uint64_t * poly , unsigned log_polylen , uint64_t scalar_a );
_clib.ibtfy_64.argtypes = [ ctypes.c_void_p , ctypes.c_uint , ctypes.c_uint64 ]
_clib.ibtfy_64.restype = None

def _log2(n):
    r = -1
    while(n):
        r = r+1
        n >>= 1
    return r

def index_to_gf264( idx ):
    r : int = _clib.index_to_gf264( idx )
    return r


def gf264_mul( a , b ):
    r : int = _clib.gf264_mul( a , b )
    return r


def btfy( list_u64 , offset_idx ):
    ss = len(list_u64)
    if 0 == ss : return
    assert( 0 == (ss&(ss-1)) )
    buffer = ctypes.create_string_buffer( b''.join( map( lambda x: x.to_bytes(8,'little') , list_u64 ) ) , 8*ss )
    _clib.btfy_64( buffer , _log2(ss) , offset_idx )
    return [ int.from_bytes( buffer[i*8:i*8+8] , 'little' ) for i in range(ss) ]


def ibtfy( list_u64 , offset_idx ):
    ss = len(list_u64)
    if 0 == ss : return
    assert( 0 == (ss&(ss-1)) )
    buffer = ctypes.create_string_buffer( b''.join( map( lambda x: x.to_bytes(8,'little') , list_u64 ) ) , 8*ss )
    _clib.ibtfy_64( buffer , _log2(ss) , offset_idx )
    return [ int.from_bytes( buffer[i*8:i*8+8] , 'little' ) for i in range(ss) ]



if __name__ == '__main__':
    a = [ i for i in range(8)]
    print( "test bfty(" , a , ")" )
    b = btfy( a , 0 )
    print( "-->" , b )
    c = ibtfy( b , 0 )
    print( "<--" , c )

